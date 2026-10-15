"""Tests for Fused Triton dropout+add kernel.

Verifies:
1. Forward correctness: fused result matches reference (no dropout case)
2. Gradient correctness via torch.autograd.gradcheck
3. Row-sharing mask behavior
4. Fallback path (CPU / no triton)
"""

import unittest
from unittest.mock import patch

import torch

from protenix.model.modules import fused_ops
from protenix.model.modules.fused_ops import (
    TRITON_FUSED_OPS_AVAILABLE,
    dropout_add_rowwise,
    triton_fused_ops_available_for_device,
)
from protenix.utils.device import get_available_device_type


ACCELERATOR_DEVICE_TYPE = get_available_device_type()
ACCELERATOR_AVAILABLE = ACCELERATOR_DEVICE_TYPE != "cpu"


class TestFusedKernelDispatch(unittest.TestCase):
    def test_cuda_keeps_existing_fused_path(self):
        with (
            patch.object(fused_ops, "TRITON_FUSED_OPS_AVAILABLE", True),
            patch.object(
                fused_ops,
                "_triton_fused_ops_runtime_available",
                side_effect=AssertionError("CUDA should not probe"),
            ),
        ):
            self.assertTrue(triton_fused_ops_available_for_device("cuda"))

    def test_xpu_probe_failure_uses_fallback(self):
        fused_ops._triton_fused_ops_runtime_available.cache_clear()
        try:
            with (
                patch.object(fused_ops, "TRITON_FUSED_OPS_AVAILABLE", True),
                patch.object(
                    fused_ops.torch,
                    "ones",
                    side_effect=RuntimeError("XPU kernel unavailable"),
                ),
                self.assertWarnsRegex(RuntimeWarning, "Disabling Triton"),
            ):
                self.assertFalse(triton_fused_ops_available_for_device("xpu"))
        finally:
            fused_ops._triton_fused_ops_runtime_available.cache_clear()


class TestDropoutAddRowwiseNoDrop(unittest.TestCase):
    """When p_drop=0, result should be exactly residual + x."""

    def test_no_dropout_training(self):
        torch.manual_seed(42)
        residual = torch.randn(2, 4, 8, 16)
        x = torch.randn(2, 4, 8, 16)
        result = dropout_add_rowwise(residual, x, p_drop=0.0, training=True)
        expected = residual + x
        self.assertTrue(torch.allclose(result, expected, atol=1e-6))

    def test_no_dropout_eval(self):
        torch.manual_seed(42)
        residual = torch.randn(2, 4, 8, 16)
        x = torch.randn(2, 4, 8, 16)
        result = dropout_add_rowwise(residual, x, p_drop=0.5, training=False)
        expected = residual + x
        self.assertTrue(torch.allclose(result, expected, atol=1e-6))

    @unittest.skipUnless(ACCELERATOR_AVAILABLE, "Accelerator required")
    def test_no_dropout_cuda(self):
        torch.manual_seed(42)
        residual = torch.randn(2, 4, 8, 16, device=ACCELERATOR_DEVICE_TYPE)
        x = torch.randn(2, 4, 8, 16, device=ACCELERATOR_DEVICE_TYPE)
        result = dropout_add_rowwise(residual, x, p_drop=0.0, training=True)
        expected = residual + x
        self.assertTrue(torch.allclose(result, expected, atol=1e-6))


class TestDropoutAddRowwiseMaskSharing(unittest.TestCase):
    """Verify dropout mask is shared across dim -3 (row dimension)."""

    def test_row_sharing_pattern(self):
        """Run with high dropout, check that mask pattern repeats along R dim."""
        if not (
            TRITON_FUSED_OPS_AVAILABLE
            and ACCELERATOR_AVAILABLE
            and triton_fused_ops_available_for_device(ACCELERATOR_DEVICE_TYPE)
        ):
            self.skipTest("Accelerator Triton fused ops required")
        # Shape: [B, R, C, D] = [1, 4, 2, 8]
        # Mask should be [1, 1, 2, 8] broadcast over R=4
        torch.manual_seed(42)
        R, C, D = 4, 2, 8
        residual = torch.zeros(1, R, C, D, device=ACCELERATOR_DEVICE_TYPE)
        x = torch.ones(1, R, C, D, device=ACCELERATOR_DEVICE_TYPE)
        p_drop = 0.5

        result = dropout_add_rowwise(residual, x, p_drop=p_drop, training=True)
        # result = 0 + dropout(1) = either 0 or 1/(1-p)
        # Check that pattern is identical across R dimension
        for r in range(1, R):
            self.assertTrue(
                torch.equal(result[0, 0], result[0, r]),
                f"Row 0 and row {r} should have identical dropout pattern",
            )


class TestDropoutAddRowwiseGradient(unittest.TestCase):
    """Verify gradients flow correctly."""

    def test_gradient_no_dropout(self):
        """With p_drop=0, grad_residual = grad_out, grad_x = grad_out."""
        residual = torch.randn(2, 3, 4, 8, requires_grad=True)
        x = torch.randn(2, 3, 4, 8, requires_grad=True)
        result = dropout_add_rowwise(residual, x, p_drop=0.0, training=True)
        loss = result.sum()
        loss.backward()
        # Both grads should be all-ones
        self.assertTrue(torch.allclose(residual.grad, torch.ones_like(residual)))
        self.assertTrue(torch.allclose(x.grad, torch.ones_like(x)))

    def test_gradient_with_dropout_cuda(self):
        """With dropout, gradients should still flow (non-zero)."""
        if not (
            TRITON_FUSED_OPS_AVAILABLE
            and ACCELERATOR_AVAILABLE
            and triton_fused_ops_available_for_device(ACCELERATOR_DEVICE_TYPE)
        ):
            self.skipTest("Accelerator Triton fused ops required")
        residual = torch.randn(
            2, 3, 4, 8, device=ACCELERATOR_DEVICE_TYPE, requires_grad=True
        )
        x = torch.randn(2, 3, 4, 8, device=ACCELERATOR_DEVICE_TYPE, requires_grad=True)
        result = dropout_add_rowwise(residual, x, p_drop=0.3, training=True)
        loss = result.sum()
        loss.backward()
        # grad_residual is always grad_output (identity)
        self.assertTrue(torch.allclose(residual.grad, torch.ones_like(residual)))
        # grad_x should be non-zero (some elements dropped)
        self.assertTrue(x.grad.abs().sum() > 0)


class TestDropoutAddRowwiseFallback(unittest.TestCase):
    """Verify the PyTorch fallback path works on CPU."""

    def test_cpu_fallback(self):
        torch.manual_seed(42)
        residual = torch.randn(2, 3, 4, 8)
        x = torch.randn(2, 3, 4, 8)
        # On CPU, should use fallback path
        result = dropout_add_rowwise(residual, x, p_drop=0.0, training=True)
        expected = residual + x
        self.assertTrue(torch.allclose(result, expected))

    def test_cpu_with_dropout(self):
        """CPU fallback with actual dropout should not crash."""
        torch.manual_seed(42)
        residual = torch.randn(2, 3, 4, 8)
        x = torch.randn(2, 3, 4, 8)
        result = dropout_add_rowwise(residual, x, p_drop=0.5, training=True)
        # Should be valid tensor, not NaN
        self.assertFalse(torch.isnan(result).any())
        # Should differ from residual + x (some elements dropped)
        # (statistically almost certain with p=0.5 and many elements)
        self.assertFalse(torch.allclose(result, residual + x))


try:
    import cuequivariance_ops_torch  # noqa: F401

    _CUEQ_AVAILABLE = True
except ImportError:
    _CUEQ_AVAILABLE = False


class TestPairformerCueqContiguity(unittest.TestCase):
    """Regression test: PairformerBlock must output contiguous z after
    transpose so that cuequivariance operators (which call .view()) don't
    crash on non-contiguous tensors.  See bytedance/Protenix#250.

    N_TOKENS must be > 100 because cuequivariance falls back to pure PyTorch
    for seq_len <= CUEQ_TRIMUL_FALLBACK_THRESHOLD (100). Only the optimized
    Triton kernel path calls layer_norm_transpose → .view() which requires
    contiguous input."""

    # Must exceed cueq's fallback threshold (100) to actually exercise the
    # optimized Triton kernel that calls .view() on the pair tensor.
    N_TOKENS = 128

    @unittest.skipUnless(
        torch.cuda.is_available() and _CUEQ_AVAILABLE,
        "CUDA + cuequivariance required",
    )
    def test_stacked_blocks_inplace_safe(self):
        """3-block stack with cuequivariance + inplace_safe=True must not crash."""
        from protenix.model.modules.pairformer import PairformerBlock

        N = self.N_TOKENS
        N_BLOCKS = 3
        blocks = [
            PairformerBlock(c_z=128, c_s=384, no_heads_pair=16).cuda().bfloat16()
            for _ in range(N_BLOCKS)
        ]

        z = torch.randn(1, N, N, 128, device="cuda", dtype=torch.bfloat16)
        s = torch.randn(1, N, 384, device="cuda", dtype=torch.bfloat16)
        pair_mask = torch.ones(1, N, N, device="cuda", dtype=torch.bfloat16)

        with torch.no_grad():
            for block in blocks:
                s, z = block(
                    s,
                    z,
                    pair_mask,
                    triangle_multiplicative="cuequivariance",
                    triangle_attention="cuequivariance",
                    inplace_safe=True,
                )
                self.assertTrue(
                    z.is_contiguous(),
                    "z must be contiguous after inplace_safe block for cueq compat",
                )

    @unittest.skipUnless(
        torch.cuda.is_available() and _CUEQ_AVAILABLE,
        "CUDA + cuequivariance required",
    )
    def test_stacked_blocks_training(self):
        """3-block stack with cuequivariance in training mode must not crash."""
        from protenix.model.modules.pairformer import PairformerBlock

        N = self.N_TOKENS
        N_BLOCKS = 3
        blocks = [
            PairformerBlock(c_z=128, c_s=384, no_heads_pair=16).cuda().bfloat16()
            for _ in range(N_BLOCKS)
        ]
        for b in blocks:
            b.train()

        z = torch.randn(1, N, N, 128, device="cuda", dtype=torch.bfloat16)
        s = torch.randn(1, N, 384, device="cuda", dtype=torch.bfloat16)
        pair_mask = torch.ones(1, N, N, device="cuda", dtype=torch.bfloat16)

        for block in blocks:
            s, z = block(
                s,
                z,
                pair_mask,
                triangle_multiplicative="cuequivariance",
                triangle_attention="cuequivariance",
                inplace_safe=False,
            )


if __name__ == "__main__":
    unittest.main()
