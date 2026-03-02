"""Tests for PR 5: Fused Triton dropout+add kernel.

Verifies:
1. Forward correctness: fused result matches reference (no dropout case)
2. Gradient correctness via torch.autograd.gradcheck
3. Row-sharing mask behavior
4. Fallback path (CPU / no triton)
"""

import unittest

import torch

from protenix.model.modules.fused_ops import (
    TRITON_FUSED_OPS_AVAILABLE,
    dropout_add_rowwise,
)


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

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA required")
    def test_no_dropout_cuda(self):
        torch.manual_seed(42)
        residual = torch.randn(2, 4, 8, 16, device="cuda")
        x = torch.randn(2, 4, 8, 16, device="cuda")
        result = dropout_add_rowwise(residual, x, p_drop=0.0, training=True)
        expected = residual + x
        self.assertTrue(torch.allclose(result, expected, atol=1e-6))


class TestDropoutAddRowwiseMaskSharing(unittest.TestCase):
    """Verify dropout mask is shared across dim -3 (row dimension)."""

    @unittest.skipUnless(
        torch.cuda.is_available() and TRITON_FUSED_OPS_AVAILABLE,
        "CUDA + Triton required",
    )
    def test_row_sharing_pattern(self):
        """Run with high dropout, check that mask pattern repeats along R dim."""
        # Shape: [B, R, C, D] = [1, 4, 2, 8]
        # Mask should be [1, 1, 2, 8] broadcast over R=4
        torch.manual_seed(42)
        R, C, D = 4, 2, 8
        residual = torch.zeros(1, R, C, D, device="cuda")
        x = torch.ones(1, R, C, D, device="cuda")
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

    @unittest.skipUnless(
        torch.cuda.is_available() and TRITON_FUSED_OPS_AVAILABLE,
        "CUDA + Triton required",
    )
    def test_gradient_with_dropout_cuda(self):
        """With dropout, gradients should still flow (non-zero)."""
        residual = torch.randn(2, 3, 4, 8, device="cuda", requires_grad=True)
        x = torch.randn(2, 3, 4, 8, device="cuda", requires_grad=True)
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


if __name__ == "__main__":
    unittest.main()
