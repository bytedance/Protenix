# SPDX-FileCopyrightText: Copyright (c) 2024 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: LicenseRef-NvidiaProprietary
#
# NVIDIA CORPORATION, its affiliates and licensors retain all intellectual
# property and proprietary rights in and to this material, related
# documentation and any modifications thereto. Any use, reproduction,
# disclosure or distribution of this material and related documentation
# without an express license agreement from NVIDIA CORPORATION or
# its affiliates is strictly prohibited.


import torch
import triton

from protenix.model.triangle.triton.fused_layer_norm import (
    Layout,
    layer_norm_transpose_backward_kernel,
    layer_norm_transpose_forward_kernel,
)


class LayerNorm(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, w, b, eps=1e-5, elementwise_affine=True, layout=Layout.BND_BND):
        # Handle autocast conversion
        if torch.is_autocast_enabled():
            dtype = torch.get_autocast_dtype("cuda")
            x = x.to(dtype=dtype)
            w = w.to(dtype=dtype)
            b = b.to(dtype=dtype)

        assert x.dim() == 3
        x = x.contiguous()
        w = w.contiguous()
        b = b.contiguous()

        if layout == Layout.BND_BND:
            B, N, D = x.shape
            out = torch.empty(B, N, D, dtype=x.dtype, device=x.device)
            # TODO tune
            TILE_N = 16
            TILE_D = 64
        elif layout == Layout.BDN_BND:
            B, D, N = x.shape
            out = torch.empty(B, N, D, dtype=x.dtype, device=x.device)
            # TODO tune
            TILE_N = 32
            TILE_D = 32
        elif layout == Layout.BND_BDN:
            B, N, D = x.shape
            out = torch.empty(B, D, N, dtype=x.dtype, device=x.device)
            TILE_N = 32
            TILE_D = 32
        elif layout == Layout.DBN_BND:
            D, B, N = x.shape
            out = torch.empty(B, N, D, dtype=x.dtype, device=x.device)
            TILE_N = 32
            TILE_D = 32
        elif layout == Layout.BND_DBN:
            B, N, D = x.shape
            out = torch.empty(D, B, N, dtype=x.dtype, device=x.device)
            TILE_N = 32
            TILE_D = 32
        else:
            raise ValueError

        mean = torch.empty(B, N, dtype=x.dtype, device=x.device)
        rstd = torch.empty(B, N, dtype=x.dtype, device=x.device)
        assert D % TILE_D == 0
        grid = (triton.cdiv(N, TILE_N), B, 1)
        layer_norm_transpose_forward_kernel[grid](
            x,
            out,
            w,
            b,
            mean,
            rstd,
            B,
            N,
            D=D,
            EPS=eps,
            TILE_N=TILE_N,
            TILE_D=TILE_D,
            ELEMENTWISE_AFFINE=elementwise_affine,
            LAYOUT=layout.value,
        )

        ctx.save_for_backward(x, w, mean, rstd)
        ctx.elementwise_affine = elementwise_affine
        ctx.layout = layout
        ctx.B = B
        ctx.N = N
        ctx.D = D

        return out

    @staticmethod
    def backward(ctx, grad_out):
        B = ctx.B
        N = ctx.N
        D = ctx.D

        grad_out = grad_out.contiguous()
        x, w, mean, rstd = ctx.saved_tensors
        elementwise_affine = ctx.elementwise_affine
        grad_x = torch.empty_like(x)

        if ctx.layout == Layout.BND_BND:
            # TODO tune
            TILE_N = 16
            TILE_D = 64
        elif ctx.layout == Layout.BDN_BND:
            # TODO tune
            TILE_N = 32
            TILE_D = 32
        elif ctx.layout == Layout.DBN_BND:
            # TODO tune
            TILE_N = 32
            TILE_D = 32
        elif ctx.layout == Layout.BND_BDN:
            # TODO tune
            TILE_N = 32
            TILE_D = 32
        elif ctx.layout == Layout.BND_DBN:
            # TODO tune
            TILE_N = 32
            TILE_D = 32
        else:
            raise ValueError

        assert ctx.D % TILE_D == 0
        num_tiles = triton.cdiv(N, TILE_N)
        grad_w = torch.empty((B, num_tiles, D), dtype=w.dtype, device=w.device)
        grad_b = torch.empty((B, num_tiles, D), dtype=w.dtype, device=w.device)

        grid = (num_tiles, B, 1)
        layer_norm_transpose_backward_kernel[grid](
            grad_out,
            grad_x,
            grad_w,
            grad_b,
            x,
            w,
            mean,
            rstd,
            B,
            N,
            D=D,
            TILE_N=TILE_N,
            TILE_D=TILE_D,
            ELEMENTWISE_AFFINE=elementwise_affine,
            LAYOUT=ctx.layout.value,
        )

        grad_w = torch.sum(grad_w, dim=0)
        grad_b = torch.sum(grad_b, dim=0)

        return grad_x, grad_w, grad_b, None, None, None


def layer_norm_transpose(
    x, w, b, eps=1e-5, elementwise_affine=True, layout: str = "nd->nd"
):
    """Apply fused layer normalization with support for various input layouts.

    This function performs layer normalization on the input tensor with optional
    elementwise affine transformation. It supports various input layouts and can
    transform between different tensor shapes.

    The normalization process consists of two steps:
    1. Normalize the input by subtracting mean and dividing by standard deviation
    2. Apply an affine transformation: output = weight * normalized_input + bias

    Args:
        x (torch.Tensor): Input tensor. Shape depends on the layout.
        w (torch.Tensor): Weight tensor for scaling the normalized values. Shape should be (D,).
            These weights allow the network to learn the optimal scale for each feature.
            Only used if elementwise_affine=True.
        b (torch.Tensor): Bias tensor for shifting the normalized values. Shape should be (D,).
            These biases allow the network to learn the optimal offset for each feature.
            Only used if elementwise_affine=True.
        eps (float, optional): Small constant for numerical stability. Defaults to 1e-5.
        elementwise_affine (bool, optional): Whether to apply elementwise affine transformation.
            If False, weight and bias are not used (equivalent to weight=1, bias=0).
            Defaults to True.
        layout (str, optional): Input/output layout specification. Defaults to "nd->nd".
            Supported layouts:
            - "nd->nd": (N, D) -> (N, D)
            - "nd->dn": (N, D) -> (D, N)
            - "bnd->bnd": (B, N, D) -> (B, N, D)
            - "bdn->bnd": (B, D, N) -> (B, N, D)
            - "bnd->bdn": (B, N, D) -> (B, D, N)
            - "dbn->bnd": (D, B, N) -> (B, N, D)
            - "bnd->dbn": (B, N, D) -> (D, B, N)
            - "bijd->bijd": (B, I, J, D) -> (B, I, J, D)
            - "bijd->bdij": (B, I, J, D) -> (B, D, I, J)
            - "bdij->bijd": (B, D, I, J) -> (B, I, J, D)
            - "dbij->bijd": (D, B, I, J) -> (B, I, J, D)
            - "bijd->dbij": (B, I, J, D) -> (D, B, I, J)

    Returns:
        torch.Tensor: Normalized tensor with shape determined by the output layout.

    Raises:
        ValueError: If the specified layout is not supported.

    Example:
        >>> x = torch.randn(1, 128, 128, device="cuda", dtype=torch.float, requires_grad=True)
        >>> w = torch.randn(128, device="cuda", dtype=torch.float, requires_grad=True)
        >>> b = torch.randn(128, device="cuda", dtype=torch.float, requires_grad=True)
        >>> out = layer_norm_transpose(x, w, b, layout="bnd->bnd")
    """
    supported_layouts = (
        "nd->nd",
        "nd->dn",
        "dn->nd",
        "bnd->bnd",
        "bnd->bdn",
        "bdn->bnd",
        "dbn->bnd",
        "bnd->dbn",
        "bijd->bijd",
        "bijd->bdij",
        "bdij->bijd",
        "dbij->bijd",
        "bijd->dbij",
    )

    if layout == "nd->nd":
        N, D = x.shape
        B = 1
        x = x.view(1, N, D)
        out_shape = (N, D)
        layout = Layout.BND_BND

    elif layout == "nd->dn":
        N, D = x.shape
        B = 1
        x = x.view(1, N, D)
        out_shape = (D, N)
        layout = Layout.BND_BDN

    elif layout == "bnd->bnd":
        B, N, D = x.shape
        out_shape = (B, N, D)
        layout = Layout.BND_BND

    elif layout == "bdn->bnd":
        B, D, N = x.shape
        out_shape = (B, N, D)
        layout = Layout.BDN_BND

    elif layout == "bnd->bdn":
        B, N, D = x.shape
        out_shape = (B, D, N)
        layout = Layout.BND_BDN

    elif layout == "dbn->bnd":
        D, B, N = x.shape
        out_shape = (B, N, D)
        layout = Layout.DBN_BND

    elif layout == "bnd->dbn":
        B, N, D = x.shape
        out_shape = (D, B, N)
        layout = Layout.BND_DBN

    elif layout == "bijd->bijd":
        B, II, J, D = x.shape
        out_shape = (B, II, J, D)
        x = x.view(B, II * J, D)
        layout = Layout.BND_BND

    elif layout == "bijd->bdij":
        B, II, J, D = x.shape
        out_shape = (B, D, II, J)
        x = x.view(B, II * J, D)
        layout = Layout.BND_BDN

    elif layout == "bdij->bijd":
        B, D, II, J = x.shape
        out_shape = (B, II, J, D)
        x = x.view(B, D, II * J)
        layout = Layout.BDN_BND

    elif layout == "dbij->bijd":
        D, B, II, J = x.shape
        out_shape = (B, II, J, D)
        x = x.view(D, B, II * J)
        layout = Layout.DBN_BND

    elif layout == "bijd->dbij":
        B, II, J, D = x.shape
        out_shape = (D, B, II, J)
        x = x.view(B, II * J, D)
        layout = Layout.BND_DBN

    else:
        raise ValueError(
            f"layout {layout} not supported. supported layouts are: {supported_layouts}"
        )

    out = LayerNorm.apply(x, w, b, eps, elementwise_affine, layout)
    return out.view(*out_shape)
