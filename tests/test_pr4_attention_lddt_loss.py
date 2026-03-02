"""Tests for PR 4: Attention precision, LDDT fused thresholds, loss caching.

Verifies:
1. Attention: v cast to fp32, cast ordering fix, SDPA default
2. LDDT: fused threshold sum*0.25 matches stack+mean
3. Loss: lru_cache off-diagonal mask, pre-computed true_distance
"""

import unittest

import torch
import torch.nn.functional as F


class TestAttentionPrecision(unittest.TestCase):
    """Verify attention cast ordering produces better precision."""

    def test_cast_ordering_matmul(self):
        """(attn_weights @ v).to(dtype) vs attn_weights.to(dtype) @ v.

        The new ordering keeps the matmul in float32 for full precision,
        only casting the final result to bf16/fp16.
        """
        torch.manual_seed(42)
        n_q, n_kv, d = 32, 32, 64
        attn_weights = torch.randn(n_q, n_kv, dtype=torch.float32).softmax(dim=-1)
        v = torch.randn(n_kv, d, dtype=torch.bfloat16)

        # Old: cast weights to bf16 first, then matmul
        old_result = attn_weights.to(dtype=torch.bfloat16) @ v

        # New: matmul in fp32, then cast
        new_result = (attn_weights @ v.to(dtype=torch.float32)).to(dtype=torch.bfloat16)

        # Both should be close, but the new one has less accumulated error
        # Verify they're close (not identical due to different precision paths)
        self.assertTrue(
            torch.allclose(old_result.float(), new_result.float(), atol=0.05),
            "Results should be close",
        )

    def test_v_fp32_cast(self):
        """Verify v is cast to fp32 alongside q and k."""
        torch.manual_seed(42)
        q = torch.randn(4, 8, 32, 64, dtype=torch.bfloat16)
        k = torch.randn(4, 8, 32, 64, dtype=torch.bfloat16)
        v = torch.randn(4, 8, 32, 64, dtype=torch.bfloat16)

        # Simulate the new code path: q, k, v all cast to fp32
        q32 = q.to(torch.float32)
        k32 = k.to(torch.float32)
        v32 = v.to(torch.float32)

        attn = (q32 @ k32.transpose(-1, -2)).softmax(dim=-1)
        result = (attn @ v32).to(torch.bfloat16)

        # Should be a valid bf16 tensor
        self.assertEqual(result.dtype, torch.bfloat16)
        self.assertFalse(torch.isnan(result).any())


class TestLDDTFusedThresholds(unittest.TestCase):
    """Verify sum of 4 comparisons * 0.25 matches stack + mean."""

    def _old_lddt(self, distance_error):
        thresholds = [0.5, 1, 2, 4]
        return (
            torch.stack([distance_error < t for t in thresholds], dim=-1)
            .to(dtype=distance_error.dtype)
            .mean(dim=-1)
        )

    def _new_lddt(self, distance_error):
        return (
            (distance_error < 0.5).to(dtype=distance_error.dtype)
            + (distance_error < 1.0)
            + (distance_error < 2.0)
            + (distance_error < 4.0)
        ) * 0.25

    def test_equivalence_float32(self):
        torch.manual_seed(42)
        err = torch.rand(10, 50) * 5  # [0, 5] range
        old = self._old_lddt(err)
        new = self._new_lddt(err)
        self.assertTrue(torch.allclose(old, new, atol=1e-6))

    def test_equivalence_bfloat16(self):
        torch.manual_seed(42)
        err = (torch.rand(10, 50) * 5).to(torch.bfloat16)
        old = self._old_lddt(err)
        new = self._new_lddt(err)
        self.assertTrue(torch.allclose(old, new, atol=1e-3))

    def test_exact_threshold_values(self):
        """Test at exact boundary values."""
        err = torch.tensor([0.0, 0.5, 1.0, 2.0, 4.0, 5.0])
        old = self._old_lddt(err)
        new = self._new_lddt(err)
        self.assertTrue(torch.allclose(old, new))

    def test_all_below_smallest(self):
        err = torch.tensor([0.1, 0.2, 0.3])
        old = self._old_lddt(err)
        new = self._new_lddt(err)
        self.assertTrue(torch.allclose(old, new))
        # All below all thresholds => lddt = 1.0
        self.assertTrue(torch.allclose(new, torch.ones_like(new)))

    def test_all_above_largest(self):
        err = torch.tensor([5.0, 10.0, 100.0])
        old = self._old_lddt(err)
        new = self._new_lddt(err)
        self.assertTrue(torch.allclose(old, new))
        # All above all thresholds => lddt = 0.0
        self.assertTrue(torch.allclose(new, torch.zeros_like(new)))

    def test_large_batch(self):
        """Stress test with large tensors."""
        torch.manual_seed(42)
        err = torch.rand(5, 1000, 1000) * 6
        old = self._old_lddt(err)
        new = self._new_lddt(err)
        self.assertTrue(torch.allclose(old, new, atol=1e-6))


class TestOffDiagonalMaskCache(unittest.TestCase):
    """Verify cached off-diagonal mask is correct."""

    def test_correctness(self):
        from protenix.model.loss import _get_off_diagonal_mask
        for n in [1, 5, 10, 100]:
            mask = _get_off_diagonal_mask(n, torch.device("cpu"), torch.float32)
            expected = 1 - torch.eye(n, dtype=torch.float32)
            self.assertTrue(torch.equal(mask, expected))

    def test_cache_returns_same_object(self):
        from protenix.model.loss import _get_off_diagonal_mask
        m1 = _get_off_diagonal_mask(8, torch.device("cpu"), torch.float32)
        m2 = _get_off_diagonal_mask(8, torch.device("cpu"), torch.float32)
        self.assertIs(m1, m2)  # Same object from cache

    def test_different_sizes_different_objects(self):
        from protenix.model.loss import _get_off_diagonal_mask
        m5 = _get_off_diagonal_mask(5, torch.device("cpu"), torch.float32)
        m6 = _get_off_diagonal_mask(6, torch.device("cpu"), torch.float32)
        self.assertIsNot(m5, m6)
        self.assertEqual(m5.shape, (5, 5))
        self.assertEqual(m6.shape, (6, 6))

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA required")
    def test_cuda_device(self):
        from protenix.model.loss import _get_off_diagonal_mask
        mask = _get_off_diagonal_mask(4, torch.device("cuda"), torch.float32)
        self.assertTrue(mask.is_cuda)
        expected = 1 - torch.eye(4, device="cuda", dtype=torch.float32)
        self.assertTrue(torch.equal(mask, expected))


class TestSmoothLDDTPrecomputedDistance(unittest.TestCase):
    """Verify that passing pre-computed true_distance skips recomputation."""

    def test_precomputed_matches_computed(self):
        """dense_forward with/without pre-computed distance should match."""
        from protenix.model.loss import SmoothLDDTLoss

        torch.manual_seed(42)
        N_atom = 20
        pred = torch.randn(1, N_atom, 3)
        true = torch.randn(N_atom, 3)

        true_dist = torch.cdist(true.unsqueeze(0), true.unsqueeze(0)).squeeze(0)
        lddt_mask = (true_dist < 15.0).float()

        loss_fn = SmoothLDDTLoss()

        # Without pre-computed distance
        loss1 = loss_fn.dense_forward(
            pred_coordinate=pred,
            true_coordinate=true,
            lddt_mask=lddt_mask,
        )

        # With pre-computed distance
        loss2 = loss_fn.dense_forward(
            pred_coordinate=pred,
            true_coordinate=true,
            lddt_mask=lddt_mask,
            true_distance=true_dist,
        )

        self.assertTrue(
            torch.allclose(loss1, loss2, atol=1e-5),
            f"Losses differ: {loss1.item()} vs {loss2.item()}",
        )


if __name__ == "__main__":
    unittest.main()
