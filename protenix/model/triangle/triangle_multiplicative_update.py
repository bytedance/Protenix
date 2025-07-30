# SPDX-FileCopyrightText: Copyright (c) 2024 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: LicenseRef-NvidiaProprietary
#
# NVIDIA CORPORATION, its affiliates and licensors retain all intellectual
# property and proprietary rights in and to this material, related
# documentation and any modifications thereto. Any use, reproduction,
# disclosure or distribution of this material and related documentation
# without an express license agreement from NVIDIA CORPORATION or
# its affiliates is strictly prohibited.
import math
from typing import Optional

import numpy as np
import torch
from scipy.stats import truncnorm

from protenix.model.triangle.precision import Precision
from protenix.model.triangle.fused_layer_norm import layer_norm_transpose
from protenix.model.triangle.gated_gemm import (
    fused_sigmoid_gated_dual_gemm,
    fused_sigmoid_gated_dual_gemm_dual_x,
)


def _prod(nums):
    out = 1
    for n in nums:
        out = out * n
    return out


def ensure_dims(ten: torch.Tensor, n: int) -> torch.Tensor:
    """Ensure tensor has at least n dimensions by adding 1-sized dimensions at the beginning.

    Args:
        ten: Input tensor
        n: Target number of dimensions

    Returns:
        Tensor with at least n dimensions
    """
    while len(ten.shape) < n:
        ten = ten.unsqueeze(0)
    return ten


def _calculate_fan(linear_weight_shape, fan="fan_in"):
    fan_out, fan_in = linear_weight_shape
    if fan == "fan_in":
        f = fan_in
    elif fan == "fan_out":
        f = fan_out
    elif fan == "fan_avg":
        f = (fan_in + fan_out) / 2
    else:
        raise ValueError("Invalid fan option")
    return f


def trunc_normal_init_(weights, scale=1.0, fan="fan_in"):
    shape = weights.shape
    f = _calculate_fan(shape, fan)
    scale = scale / max(1, f)
    a = -2
    b = 2
    std = math.sqrt(scale) / truncnorm.std(a=a, b=b, loc=0, scale=1)
    size = _prod(shape)
    samples = truncnorm.rvs(a=a, b=b, loc=0, scale=std, size=size)
    samples = np.reshape(samples, shape)
    with torch.no_grad():
        weights.copy_(torch.tensor(samples, device=weights.device))


def lecun_normal_init_(weights):
    trunc_normal_init_(weights, scale=1.0)


def bias_init_zero_(bias):
    with torch.no_grad():
        bias.fill_(0.0)


def bias_init_one_(bias):
    with torch.no_grad():
        bias.fill_(1.0)


def triangle_multiplicative_update(
    x: torch.Tensor,
    direction: str = "outgoing",
    mask: Optional[torch.Tensor] = None,
    norm_in_weight: Optional[torch.Tensor] = None,
    norm_in_bias: Optional[torch.Tensor] = None,
    p_in_weight: Optional[torch.Tensor] = None,
    g_in_weight: Optional[torch.Tensor] = None,
    norm_out_weight: Optional[torch.Tensor] = None,
    norm_out_bias: Optional[torch.Tensor] = None,
    p_out_weight: Optional[torch.Tensor] = None,
    g_out_weight: Optional[torch.Tensor] = None,
    eps: float = 1e-5,
    precision: Optional[Precision] = None,
) -> torch.Tensor:
    """Apply triangle multiplicative update operation.

    This function performs a triangle multiplicative update operation, which is a key component
    in the AlphaFold2 architecture. The operation consists of:

    1. Input normalization and gating
    2. Triangular projection (either outgoing or incoming)
    3. Output normalization and gating

    The function supports both ahead-of-time (AOT) tuning and just-in-time (JIT) tuning.
    Auto-tuning behavior can be controlled through environment variables:

    - Quick testing: Default configuration where tuning configs, if existent, are looked-up. If not, then falls back to default kernel parameters. No tuning is performed.
    - On-Demand tuning: Set `CUEQ_TRITON_TUNING_MODE = "ONDEMAND"` to auto-tune for new shapes encountered on first run (may take several minutes)
    - AOT tuning: Set `CUEQ_TRITON_TUNING_MODE = "AOT"` to perform full ahead-of-time tuning for optimal performance **(may take several hours)**
    - Ignore user cache: Set CUEQ_TRITON_IGNORE_EXISTING_CACHE to ignore both the default settings that come with the package and any user-local settings previously saved with AOT/ONDEMAND tuning. May be used to regenerate optimal settings for a particular setup.
    - Cache directory: Set `CUEQ_TRITON_CACHE_DIR` to specify where tuning configurations are stored
    - Note: When using Docker with default or on-demand tuning enabled, commit the container to persist tuning changes

    Args:
        x (torch.Tensor): Input tensor of shape (B, N, N, D) where:
            B is the batch size
            N is the sequence length
            D is the hidden dimension
        direction (str): Direction of the triangular projection. Must be either "outgoing" or "incoming".
        mask (torch.Tensor): Optional Mask tensor of shape (B, N, N) for masking the output.
        norm_in_weight (torch.Tensor): Weight tensor for input normalization of shape (D,).
        norm_in_bias (torch.Tensor): Bias tensor for input normalization of shape (D,).
        p_in_weight (torch.Tensor): Weight tensor for input projection of shape (2D, D).
        g_in_weight (torch.Tensor): Weight tensor for input gating of shape (2D, D).
        norm_out_weight (torch.Tensor): Weight tensor for output normalization of shape (D,).
        norm_out_bias (torch.Tensor): Bias tensor for output normalization of shape (D,).
        p_out_weight (torch.Tensor): Weight tensor for output projection of shape (D, D).
        g_out_weight (torch.Tensor): Weight tensor for output gating of shape (D, D).
        eps (float, optional): Small constant for numerical stability in normalization. Defaults to 1e-5.
        precision (Precision, optional): Precision mode for matrix multiplications. If None, uses TF32 if enabled in PyTorch using torch.backends.cuda.matmul.allow_tf32, otherwise uses default precision.
            Available options:
            - DEFAULT: Use default precision setting of triton.language.dot
            - TF32: Use TensorFloat-32 precision
            - TF32x3: Use TensorFloat-32 precision with 3x accumulation
            - IEEE: Use IEEE 754 precision

    Returns:
        Output tensor of shape (batch_size, seq_len, seq_len, hidden_dim)

    Notes:
        (1) Context is saved for backward pass. You don't need to save it manually.
        (2) Kernel precision (fp32, bf16, fp16) is based on input dtypes. For tf32, set it from torch global scope using torch.backends.cuda.matmul.allow_tf32
        (3) **Limitation**: Currently only supports hidden_dim values that are multiples of 32.

    Example:
        >>> import torch
        >>> from cuequivariance_torch import triangle_multiplicative_update
        >>> if torch.cuda.is_available():  # doctest: +SKIP
        ...     device = torch.device("cuda")
        ...     batch_size, seq_len, hidden_dim = 1, 128, 128
        ...     # Create input tensor
        ...     x = torch.randn(batch_size, seq_len, seq_len, hidden_dim, requires_grad=True, device=device)
        ...     # Create mask (1 for valid positions, 0 for masked)
        ...     mask = torch.ones(batch_size, seq_len, seq_len, device=device)
        ...     # Perform triangular multiplication
        ...     output = triangle_multiplicative_update(
        ...         x=x,
        ...         direction="outgoing",  # or "incoming"
        ...         mask=mask,
        ...     )
        ...     print(output.shape)  # torch.Size([1, 128, 128, 128])
        ...     # Create gradient tensor and perform backward pass
        ...     grad_out = torch.randn_like(output)
        ...     output.backward(grad_out)
        ...     # Access gradients
        ...     print(x.grad.shape)  # torch.Size([1, 128, 128, 128])
        torch.Size([1, 128, 128, 128])
        torch.Size([1, 128, 128, 128])
    """
    # Set precision based on torch.backends.cuda.matmul.allow_tf32 if not provided
    if precision is None:
        precision = (
            Precision.TF32
            if torch.backends.cuda.matmul.allow_tf32
            else Precision.DEFAULT
        )
    # Input validation
    if direction not in ["outgoing", "incoming"]:
        raise ValueError("direction must be either 'outgoing' or 'incoming'")

    # Ensure x has 4 dimensions
    x = ensure_dims(x, 4)
    if x.dim() != 4:
        raise ValueError(
            "x must be 4-dimensional (batch_size, seq_len, seq_len, hidden_dim) or lower dimensional where first dimensions with size 1 are omitted"
        )

    if mask is not None:
        # Ensure mask has 3 dimensions
        mask = ensure_dims(mask, 3)
        if mask.dim() != 3:
            raise ValueError(
                "mask must be 3-dimensional (batch_size, seq_len, seq_len) or lower dimensional where first dimensions with size 1 are omitted"
            )

    # Initialize default weights if not provided
    hidden_dim = x.shape[-1]
    if norm_in_weight is None:
        norm_in_weight = torch.empty(hidden_dim, device=x.device, dtype=x.dtype)
        bias_init_one_(norm_in_weight)
    if norm_in_bias is None:
        norm_in_bias = torch.empty(hidden_dim, device=x.device, dtype=x.dtype)
        bias_init_zero_(norm_in_bias)
    if p_in_weight is None:
        p_in_weight = torch.empty(
            2 * hidden_dim, hidden_dim, device=x.device, dtype=x.dtype
        )
        lecun_normal_init_(p_in_weight)
    if g_in_weight is None:
        g_in_weight = torch.empty(
            2 * hidden_dim, hidden_dim, device=x.device, dtype=x.dtype
        )
        lecun_normal_init_(g_in_weight)
    if norm_out_weight is None:
        norm_out_weight = torch.empty(hidden_dim, device=x.device, dtype=x.dtype)
        bias_init_one_(norm_out_weight)
    if norm_out_bias is None:
        norm_out_bias = torch.empty(hidden_dim, device=x.device, dtype=x.dtype)
        bias_init_zero_(norm_out_bias)
    if p_out_weight is None:
        p_out_weight = torch.empty(
            hidden_dim, hidden_dim, device=x.device, dtype=x.dtype
        )
        lecun_normal_init_(p_out_weight)
    if g_out_weight is None:
        g_out_weight = torch.empty(
            hidden_dim, hidden_dim, device=x.device, dtype=x.dtype
        )
        lecun_normal_init_(g_out_weight)

    # Input normalization
    x = layer_norm_transpose(
        x, norm_in_weight, norm_in_bias, eps=eps, layout="bijd->bijd"
    )
    x_in = x

    # Gated dual gemm
    ab = fused_sigmoid_gated_dual_gemm(
        x, g_in_weight, p_in_weight, mask, transpose_out=True, precision=precision
    )
    a, b = torch.chunk(ab, 2, dim=0)

    # Triangular projection
    if direction == "outgoing":
        x = torch.einsum("dbik,dbjk->dbij", a, b)
    else:
        x = torch.einsum("dbki,dbkj->dbij", a, b)

    # Output normalization
    x_out = layer_norm_transpose(
        x, norm_out_weight, norm_out_bias, eps=eps, layout="dbij->bijd"
    )
    # Output gating
    x = fused_sigmoid_gated_dual_gemm_dual_x(
        x_in, x_out, g_out_weight, p_out_weight, precision=precision
    )

    return x
