# SPDX-FileCopyrightText: Copyright (c) 2024 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: LicenseRef-NvidiaProprietary
#
# NVIDIA CORPORATION, its affiliates and licensors retain all intellectual
# property and proprietary rights in and to this material, related
# documentation and any modifications thereto. Any use, reproduction,
# disclosure or distribution of this material and related documentation
# without an express license agreement from NVIDIA CORPORATION or
# its affiliates is strictly prohibited.

import enum

import triton
import triton.language as tl


class Layout(enum.IntEnum):
    BND_BND = 0
    BDN_BND = 1
    BND_BDN = 2
    DBN_BND = 3
    BND_DBN = 4


@triton.jit
def layer_norm_transpose_forward_kernel(
    x_ptr,
    out_ptr,
    w_ptr,
    b_ptr,
    mean_ptr,
    rstd_ptr,
    B,
    N,
    D: tl.constexpr,
    EPS: tl.constexpr,
    TILE_N: tl.constexpr,
    TILE_D: tl.constexpr,
    ELEMENTWISE_AFFINE: tl.constexpr,
    LAYOUT: tl.constexpr,
):
    pid_n = tl.program_id(0)
    pid_b = tl.program_id(1)

    offs_n = pid_n * TILE_N + tl.arange(0, TILE_N)
    offs_d = tl.arange(0, TILE_D)

    if LAYOUT == 0:  # bnd->bnd
        x_ptrs = x_ptr + pid_b * N * D + offs_n[:, None] * D + offs_d[None, :]
    elif LAYOUT == 1:  # bdn->bnd
        x_ptrs = x_ptr + pid_b * D * N + offs_d[None, :] * N + offs_n[:, None]
    elif LAYOUT == 2:  # bnd->bdn
        x_ptrs = x_ptr + pid_b * N * D + offs_n[:, None] * D + offs_d[None, :]
    elif LAYOUT == 3:  # dbn->bnd
        x_ptrs = x_ptr + offs_d[None, :] * B * N + pid_b * N + offs_n[:, None]
    elif LAYOUT == 4:  # bnd->dbn
        x_ptrs = x_ptr + pid_b * N * D + offs_n[:, None] * D + offs_d[None, :]

    mean_ptrs = mean_ptr + pid_b * N + offs_n
    rstd_ptrs = rstd_ptr + pid_b * N + offs_n
    mask_n = offs_n < N

    num_tiles = D // TILE_D

    _mean = tl.zeros([TILE_N, TILE_D], dtype=tl.float32)
    for _ in range(0, num_tiles):
        x = tl.load(x_ptrs, mask=mask_n[:, None], other=0.0).to(tl.float32)
        _mean += x

        if LAYOUT == 0:  # bnd->bnd
            x_ptrs += TILE_D
        elif LAYOUT == 1:  # bdn->bnd
            x_ptrs += TILE_D * N
        elif LAYOUT == 2:  # bnd->bdn
            x_ptrs += TILE_D
        elif LAYOUT == 3:  # dbn->bnd
            x_ptrs += TILE_D * B * N
        elif LAYOUT == 4:  # bnd->dbn
            x_ptrs += TILE_D

    mean = tl.sum(_mean, axis=1) / D
    tl.store(mean_ptrs, mean, mask=mask_n)

    if LAYOUT == 0:  # bnd->bnd
        x_ptrs -= D
    elif LAYOUT == 1:  # bdn->bnd
        x_ptrs -= D * N
    elif LAYOUT == 2:  # bnd->bdn
        x_ptrs -= D
    elif LAYOUT == 3:  # dbn->bnd
        x_ptrs -= D * B * N
    elif LAYOUT == 4:  # bnd->dbn
        x_ptrs -= D

    _var = tl.zeros([TILE_N, TILE_D], dtype=tl.float32)
    for d in range(0, num_tiles):
        x = tl.load(x_ptrs, mask=mask_n[:, None], other=0.0).to(tl.float32)
        x = x - mean[:, None]
        _var += x * x

        if LAYOUT == 0:  # bnd->bnd
            x_ptrs += TILE_D
        elif LAYOUT == 1:  # bdn->bnd
            x_ptrs += TILE_D * N
        elif LAYOUT == 2:  # bnd->bdn
            x_ptrs += TILE_D
        elif LAYOUT == 3:  # dbn->bnd
            x_ptrs += TILE_D * B * N
        elif LAYOUT == 4:  # bnd->dbn
            x_ptrs += TILE_D

    var = tl.sum(_var, axis=1) / D
    rstd = 1.0 / tl.sqrt(var + EPS)
    tl.store(rstd_ptrs, rstd, mask=mask_n)

    if LAYOUT == 0:  # bnd->bnd
        x_ptrs -= D
        out_ptrs = out_ptr + pid_b * N * D + offs_n[:, None] * D + offs_d[None, :]
    elif LAYOUT == 1:  # bdn->bnd
        x_ptrs -= D * N
        out_ptrs = out_ptr + pid_b * N * D + offs_n[:, None] * D + offs_d[None, :]
    elif LAYOUT == 2:  # bnd->bdn
        x_ptrs -= D
        out_ptrs = out_ptr + pid_b * N * D + offs_d[None, :] * N + offs_n[:, None]
    elif LAYOUT == 3:  # dbn->bnd
        x_ptrs -= D * B * N
        out_ptrs = out_ptr + pid_b * N * D + offs_n[:, None] * D + offs_d[None, :]
    elif LAYOUT == 4:  # bnd->dbn
        x_ptrs -= D
        out_ptrs = out_ptr + offs_d[None, :] * B * N + pid_b * N + offs_n[:, None]

    if ELEMENTWISE_AFFINE:
        w_ptrs = w_ptr + offs_d
        b_ptrs = b_ptr + offs_d

    for _ in range(0, num_tiles):
        if ELEMENTWISE_AFFINE:
            w = tl.load(w_ptrs)
            b = tl.load(b_ptrs)
        else:
            w = 1.0
            b = 0.0

        x = tl.load(x_ptrs, mask=mask_n[:, None], other=0.0).to(tl.float32)
        x_hat = (x - mean[:, None]) * rstd[:, None]
        y = x_hat * w[None, :] + b[None, :]
        tl.store(out_ptrs, y, mask=mask_n[:, None])

        if LAYOUT == 0:  # bnd->bnd
            x_ptrs += TILE_D
            out_ptrs += TILE_D
        elif LAYOUT == 1:  # bdn->bnd
            x_ptrs += TILE_D * N
            out_ptrs += TILE_D
        elif LAYOUT == 2:  # bnd->bdn
            x_ptrs += TILE_D
            out_ptrs += TILE_D * N
        elif LAYOUT == 3:  # dbn->bnd
            x_ptrs += TILE_D * B * N
            out_ptrs += TILE_D
        elif LAYOUT == 4:  # bnd->dbn
            x_ptrs += TILE_D
            out_ptrs += TILE_D * B * N

        if ELEMENTWISE_AFFINE:
            w_ptrs += TILE_D
            b_ptrs += TILE_D


@triton.jit
def layer_norm_transpose_backward_kernel(
    grad_out_ptr,
    grad_x_ptr,
    grad_w_ptr,
    grad_b_ptr,
    x_ptr,
    w_ptr,
    mean_ptr,
    rstd_ptr,
    B,
    N,
    D: tl.constexpr,
    TILE_N: tl.constexpr,
    TILE_D: tl.constexpr,
    ELEMENTWISE_AFFINE: tl.constexpr,
    LAYOUT: tl.constexpr,
):
    pid_n = tl.program_id(0)
    pid_b = tl.program_id(1)

    num_tiles = D // TILE_D
    num_tiles_n = tl.cdiv(N, TILE_N)

    offs_d = tl.arange(0, TILE_D)
    offs_n = pid_n * TILE_N + tl.arange(0, TILE_N)
    mask_n = offs_n < N

    mean_ptrs = mean_ptr + pid_b * N + offs_n
    rstd_ptrs = rstd_ptr + pid_b * N + offs_n
    mean = tl.load(mean_ptrs, mask=mask_n, other=0.0).to(tl.float32)
    rstd = tl.load(rstd_ptrs, mask=mask_n, other=0.0).to(tl.float32)

    if LAYOUT == 0:  # bnd->bnd
        x_base_ptrs = x_ptr + pid_b * N * D + offs_n[:, None] * D
        grad_x_base_ptrs = grad_x_ptr + pid_b * N * D + offs_n[:, None] * D
        grad_out_base_ptrs = grad_out_ptr + pid_b * N * D + offs_n[:, None] * D
    elif LAYOUT == 1:  # bdn->bnd
        x_base_ptrs = x_ptr + pid_b * D * N + offs_n[:, None]
        grad_x_base_ptrs = grad_x_ptr + pid_b * D * N + offs_n[:, None]
        grad_out_base_ptrs = grad_out_ptr + pid_b * N * D + offs_n[:, None] * D
    elif LAYOUT == 2:  # bnd->bdn
        x_base_ptrs = x_ptr + pid_b * N * D + offs_n[:, None] * D
        grad_x_base_ptrs = grad_x_ptr + pid_b * N * D + offs_n[:, None] * D
        grad_out_base_ptrs = grad_out_ptr + pid_b * N * D + offs_n[:, None]
    elif LAYOUT == 3:  # dbn->bnd
        x_base_ptrs = x_ptr + pid_b * N + offs_n[:, None]
        grad_x_base_ptrs = grad_x_ptr + pid_b * N + offs_n[:, None]
        grad_out_base_ptrs = grad_out_ptr + pid_b * N * D + offs_n[:, None] * D
    elif LAYOUT == 4:  # bnd->dbn
        x_base_ptrs = x_ptr + pid_b * N * D + offs_n[:, None] * D
        grad_x_base_ptrs = grad_x_ptr + pid_b * N * D + offs_n[:, None] * D
        grad_out_base_ptrs = grad_out_ptr + pid_b * N + offs_n[:, None]

    grad_w_base_ptrs = grad_w_ptr + pid_b * num_tiles_n * D + pid_n * D
    grad_b_base_ptrs = grad_b_ptr + pid_b * num_tiles_n * D + pid_n * D

    c1 = tl.zeros([TILE_N, TILE_D], dtype=tl.float32)
    c2 = tl.zeros([TILE_N, TILE_D], dtype=tl.float32)

    for _ in range(num_tiles):
        if ELEMENTWISE_AFFINE:
            w_ptrs = w_ptr + offs_d
            w = tl.load(w_ptrs).to(tl.float32)
        else:
            w = 1.0

        if LAYOUT == 0:  # bnd->bnd
            x_ptrs = x_base_ptrs + offs_d[None, :]
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :]
        elif LAYOUT == 1:  # bdn->bnd
            x_ptrs = x_base_ptrs + offs_d[None, :] * N
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :]
        elif LAYOUT == 2:  # bnd->bdn
            x_ptrs = x_base_ptrs + offs_d[None, :]
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :] * N
        elif LAYOUT == 3:  # dbn->bnd
            x_ptrs = x_base_ptrs + offs_d[None, :] * B * N
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :]
        elif LAYOUT == 4:  # bnd->dbn
            x_ptrs = x_base_ptrs + offs_d[None, :]
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :] * B * N

        x = tl.load(x_ptrs, mask=mask_n[:, None], other=0.0).to(tl.float32)
        grad_out = tl.load(grad_out_ptrs, mask=mask_n[:, None], other=0.0).to(
            tl.float32
        )

        xhat = (x - mean[:, None]) * rstd[:, None]

        if ELEMENTWISE_AFFINE:
            grad_b = grad_out
            grad_w = grad_out * xhat

            grad_b = tl.sum(grad_b, axis=0)
            grad_w = tl.sum(grad_w, axis=0)

            grad_w_ptrs = grad_w_base_ptrs + offs_d
            grad_b_ptrs = grad_b_base_ptrs + offs_d

            tl.store(grad_w_ptrs, grad_w)
            tl.store(grad_b_ptrs, grad_b)

        wdo = w * grad_out

        c1 += xhat * wdo
        c2 += wdo

        offs_d += TILE_D

    c1_dot = tl.sum(c1, axis=1) / D
    c2_dot = tl.sum(c2, axis=1) / D

    offs_d -= TILE_D * num_tiles

    for _ in range(num_tiles):
        if ELEMENTWISE_AFFINE:
            w_ptrs = w_ptr + offs_d
            w = tl.load(w_ptrs).to(tl.float32)
        else:
            w = 1.0

        if LAYOUT == 0:  # bnd->bnd
            x_ptrs = x_base_ptrs + offs_d[None, :]
            grad_x_ptrs = grad_x_base_ptrs + offs_d[None, :]
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :]
        elif LAYOUT == 1:  # bdn->bnd
            x_ptrs = x_base_ptrs + offs_d[None, :] * N
            grad_x_ptrs = grad_x_base_ptrs + offs_d[None, :] * N
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :]
        elif LAYOUT == 2:  # bnd->bdn
            x_ptrs = x_base_ptrs + offs_d[None, :]
            grad_x_ptrs = grad_x_base_ptrs + offs_d[None, :]
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :] * N
        elif LAYOUT == 3:  # dbn->bnd
            x_ptrs = x_base_ptrs + offs_d[None, :] * B * N
            grad_x_ptrs = grad_x_base_ptrs + offs_d[None, :] * B * N
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :]
        elif LAYOUT == 4:  # bnd->dbn
            x_ptrs = x_base_ptrs + offs_d[None, :]
            grad_x_ptrs = grad_x_base_ptrs + offs_d[None, :]
            grad_out_ptrs = grad_out_base_ptrs + offs_d[None, :] * B * N

        x = tl.load(x_ptrs, mask=mask_n[:, None], other=0.0).to(tl.float32)
        grad_out = tl.load(grad_out_ptrs, mask=mask_n[:, None], other=0.0).to(
            tl.float32
        )

        xhat = (x - mean[:, None]) * rstd[:, None]
        wdo = w * grad_out

        dx = (wdo - (xhat * c1_dot[:, None] + c2_dot[:, None])) * rstd[:, None]
        tl.store(grad_x_ptrs, dx, mask=mask_n[:, None])

        offs_d += TILE_D
