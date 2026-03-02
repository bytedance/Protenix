# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Fused Triton kernels for common operations.

Provides fused dropout + residual add with row-wise mask sharing,
replacing the separate Dropout → multiply → add pattern to reduce
memory traffic and intermediate allocations.
"""

from typing import Optional

import torch

try:
    import triton
    import triton.language as tl

    TRITON_FUSED_OPS_AVAILABLE = True
except (ImportError, RuntimeError):
    TRITON_FUSED_OPS_AVAILABLE = False


if TRITON_FUSED_OPS_AVAILABLE:

    @triton.jit
    def _dropout_add_rowwise_fwd_kernel(
        RESIDUAL_ptr,
        X_ptr,
        OUT_ptr,
        N_total,  # B * R * C * D
        CD,  # C * D (mask repeats every CD elements within a batch)
        RCD,  # R * C * D (stride for batch dimension)
        p_drop,
        seed,
        BLOCK: tl.constexpr,
    ):
        """Forward kernel: out = residual + dropout_rowwise(x).

        The dropout mask is shared across the R (row) dimension, i.e. for a
        tensor of shape [..., R, C, D], the mask has shape [..., 1, C, D] and
        is broadcast over R.  We achieve this by computing the mask index as
        ``batch_idx * CD + (flat_idx % CD)`` which strips the row component.
        """
        pid = tl.program_id(0)
        offsets = pid * BLOCK + tl.arange(0, BLOCK)
        valid = offsets < N_total

        # Compute mask index: skip R dimension for row-sharing
        batch_idx = offsets // RCD
        cd_idx = offsets % CD
        mask_idx = batch_idx * CD + cd_idx

        # Load
        r = tl.load(RESIDUAL_ptr + offsets, mask=valid, other=0.0)
        x = tl.load(X_ptr + offsets, mask=valid, other=0.0)

        # Dropout mask (row-shared via mask_idx)
        rand_vals = tl.rand(seed, mask_idx)
        keep = rand_vals >= p_drop
        scale = 1.0 / (1.0 - p_drop)

        # Fused dropout + residual add
        out = r + tl.where(keep, x * scale, tl.zeros_like(x))
        tl.store(OUT_ptr + offsets, out, mask=valid)

    @triton.jit
    def _dropout_add_rowwise_bwd_kernel(
        GRAD_OUT_ptr,
        GRAD_X_ptr,
        N_total,
        CD,
        RCD,
        p_drop,
        seed,
        BLOCK: tl.constexpr,
    ):
        """Backward kernel: grad_x = grad_out * dropout_mask * scale.

        Regenerates the same mask from the saved seed so we don't need to
        store the mask tensor.  grad_residual = grad_out (identity, handled
        outside this kernel).
        """
        pid = tl.program_id(0)
        offsets = pid * BLOCK + tl.arange(0, BLOCK)
        valid = offsets < N_total

        batch_idx = offsets // RCD
        cd_idx = offsets % CD
        mask_idx = batch_idx * CD + cd_idx

        grad_out = tl.load(GRAD_OUT_ptr + offsets, mask=valid, other=0.0)

        rand_vals = tl.rand(seed, mask_idx)
        keep = rand_vals >= p_drop
        scale = 1.0 / (1.0 - p_drop)

        grad_x = tl.where(keep, grad_out * scale, tl.zeros_like(grad_out))
        tl.store(GRAD_X_ptr + offsets, grad_x, mask=valid)

    class _DropoutAddRowwise(torch.autograd.Function):
        """Autograd wrapper around the fused dropout + residual add kernels."""

        @staticmethod
        def forward(
            ctx,
            residual: torch.Tensor,
            x: torch.Tensor,
            p_drop: float,
        ) -> torch.Tensor:
            assert x.ndim >= 3, "Need at least 3 dims for rowwise dropout"
            assert residual.shape == x.shape

            residual = residual.contiguous()
            x = x.contiguous()

            *batch_dims, R, C, D = x.shape
            B = 1
            for s in batch_dims:
                B *= s
            N_total = B * R * C * D
            CD = C * D
            RCD = R * C * D

            seed = torch.randint(0, 2**31 - 1, (1,), device="cpu").item()
            out = torch.empty_like(x)

            BLOCK = 1024
            grid = ((N_total + BLOCK - 1) // BLOCK,)
            _dropout_add_rowwise_fwd_kernel[grid](
                residual,
                x,
                out,
                N_total,
                CD,
                RCD,
                p_drop,
                seed,
                BLOCK=BLOCK,
            )

            ctx.save_for_backward(torch.tensor([seed], dtype=torch.int64))
            ctx.p_drop = p_drop
            ctx.shape_info = (N_total, CD, RCD)
            return out

        @staticmethod
        def backward(ctx, grad_output):
            (seed_tensor,) = ctx.saved_tensors
            seed = seed_tensor.item()
            p_drop = ctx.p_drop
            N_total, CD, RCD = ctx.shape_info

            # grad_residual is just grad_output (addition is identity)
            grad_residual = grad_output

            grad_output = grad_output.contiguous()
            grad_x = torch.empty_like(grad_output)

            BLOCK = 1024
            grid = ((N_total + BLOCK - 1) // BLOCK,)
            _dropout_add_rowwise_bwd_kernel[grid](
                grad_output,
                grad_x,
                N_total,
                CD,
                RCD,
                p_drop,
                seed,
                BLOCK=BLOCK,
            )

            return grad_residual, grad_x, None


def dropout_add_rowwise(
    residual: torch.Tensor,
    x: torch.Tensor,
    p_drop: float,
    training: bool,
) -> torch.Tensor:
    """Fused dropout (row-shared mask) + residual add.

    Replaces the pattern ``z = z + dropout_rowwise(x)`` with a single fused
    operation that avoids materializing the intermediate dropout result.

    The dropout mask is shared along dim -3 (row-wise), matching the behavior
    of ``DropoutRowwise``.

    Args:
        residual: The tensor to add to. Shape ``[..., R, C, D]``.
        x: The tensor to apply dropout to. Same shape as *residual*.
        p_drop: Dropout probability.
        training: Whether in training mode.

    Returns:
        ``residual + dropout_rowwise(x)`` with shape ``[..., R, C, D]``.
    """
    if not training or p_drop == 0.0:
        return residual + x

    if TRITON_FUSED_OPS_AVAILABLE and x.is_cuda:
        return _DropoutAddRowwise.apply(residual, x, p_drop)

    # Fallback: standard PyTorch path
    shape = list(x.shape)
    shape[-3] = 1
    mask = x.new_ones(shape)
    mask = torch.nn.functional.dropout(mask, p=p_drop, training=True)
    return residual + x * mask
