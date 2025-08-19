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

import triton
import triton.language as tl

from protenix.model.tri_attention.autotune import autotune
from protenix.model.tri_attention.autotune_helpers import (
    _attention_bwd_dq_configs,
    _attention_bwd_dkdv_configs,
    _attention_bwd_preprocess_configs,
    _attention_bwd_dbias2_configs,
)


@autotune(
    configs=_attention_bwd_preprocess_configs,
    key=["H", "HEAD_DIM", "TAG_N"],
)
@triton.jit
def _attention_bwd_preprocess(
    O,
    DO,
    Delta,
    N,
    S,
    H,
    HEAD_DIM: tl.constexpr,
    TAG_N: tl.constexpr,
    BLOCK_M: tl.constexpr,
):
    bn = tl.program_id(2)
    h = tl.program_id(1)
    s = tl.program_id(0) * BLOCK_M
    b = bn // N
    n = bn % N

    Bs = N * S * H * HEAD_DIM
    Ns = S * H * HEAD_DIM
    Ss = H * HEAD_DIM
    Hs = HEAD_DIM

    o_offset = b * Bs + n * Ns + h * Hs
    delta_offset = b * N * H * S + n * S * H + h

    o_ptr = tl.make_block_ptr(
        base=O + o_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    do_ptr = tl.make_block_ptr(
        base=DO + o_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    delta_ptr = tl.make_block_ptr(
        base=Delta + delta_offset,
        shape=(S, 1),
        strides=(H, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, 1),
        order=(0, 1),
    )

    o = tl.load(o_ptr, boundary_check=(0, 1), padding_option="zero")
    do = tl.load(do_ptr, boundary_check=(0, 1), padding_option="zero")
    delta = tl.sum(o.to(tl.float32) * do.to(tl.float32), axis=1)
    tl.store(delta_ptr, delta[:, None], boundary_check=(0, 1))


@autotune(
    configs=_attention_bwd_dkdv_configs,
    key=["H", "HEAD_DIM", "TAG_N"],
)
@triton.jit
def _attention_bwd_dkdv(
    Q,
    K,
    V,
    Bias1,
    Bias2,
    M,
    Delta,
    DK,
    DV,
    DO,
    N,
    S,
    H,
    HEAD_DIM: tl.constexpr,
    TAG_N: tl.constexpr,
    deterministic: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
):

    input_dtype = Q.dtype.element_ty

    bn = tl.program_id(2)
    h = tl.program_id(1)
    s = tl.program_id(0) * BLOCK_M
    b = bn // N
    n = bn % N

    Bs = N * S * H * HEAD_DIM
    Ns = S * H * HEAD_DIM
    Ss = H * HEAD_DIM
    Hs = HEAD_DIM

    qkvo_offset = b * Bs + n * Ns + h * Hs
    bias1_offset = b * N * S + n * S
    bias2_offset = b * S * S * H + h * S * S
    delta_m_offset = b * N * H * S + n * S * H + h

    q_ptr = tl.make_block_ptr(
        base=Q + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(0, 0),
        block_shape=(BLOCK_N, HEAD_DIM),
        order=(1, 0),
    )

    k_ptr = tl.make_block_ptr(
        base=K + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    dk_ptr = tl.make_block_ptr(
        base=DK + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    v_ptr = tl.make_block_ptr(
        base=V + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    dv_ptr = tl.make_block_ptr(
        base=DV + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    do_ptr = tl.make_block_ptr(
        base=DO + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(0, 0),
        block_shape=(BLOCK_N, HEAD_DIM),
        order=(0, 1),
    )

    bias1_ptr = tl.make_block_ptr(
        base=Bias1 + bias1_offset,
        shape=(S, S),
        strides=(0, 1),
        offsets=(0, s),
        block_shape=(BLOCK_N, BLOCK_M),
        order=(1, 0),
    )

    bias2_ptr = tl.make_block_ptr(
        base=Bias2 + bias2_offset,
        shape=(S, S),
        strides=(S, 1),
        offsets=(0, s),
        block_shape=(BLOCK_N, BLOCK_M),
        order=(1, 0),
    )

    m_ptr = tl.make_block_ptr(
        base=M + delta_m_offset,
        shape=(S, 1),
        strides=(H, 1),
        offsets=(0, 0),
        block_shape=(BLOCK_N, 1),
        order=(0, 1),
    )

    delta_ptr = tl.make_block_ptr(
        base=Delta + delta_m_offset,
        shape=(S, 1),
        strides=(H, 1),
        offsets=(0, 0),
        block_shape=(BLOCK_N, 1),
        order=(0, 1),
    )

    k = tl.load(
        k_ptr, boundary_check=(0, 1), padding_option="zero"
    )  # (HEAD_DIM, BLOCK_M)
    v = tl.load(
        v_ptr, boundary_check=(0, 1), padding_option="zero"
    )  # (BLOCK_M, HEAD_DIM)
    dv = tl.zeros([BLOCK_M, HEAD_DIM], dtype=tl.float32)  # (HEAD_DIM, BLOCK_M)
    dk = tl.zeros([BLOCK_M, HEAD_DIM], dtype=tl.float32)  # (BLOCK_M, HEAD_DIM)
    qk_scale = 1.4426950408889634

    for _ in range(tl.cdiv(S, BLOCK_N)):
        q = tl.load(
            q_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, HEAD_DIM)
        bias1 = tl.load(
            bias1_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, BLOCK_M)
        bias2 = tl.load(
            bias2_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, BLOCK_M)
        m = tl.load(m_ptr, boundary_check=(0, 1), padding_option="zero")  # (BLOCK_N, 1)
        qk = bias1.to(tl.float32) + bias2.to(tl.float32)
        qk = (
            tl.dot(
                q, tl.trans(k), qk, input_precision="ieee" if deterministic else "tf32"
            )
        ) * qk_scale  # (BLOCK_N, BLOCK_M)
        p = tl.math.exp2(qk - m)  # (BLOCK_N, BLOCK_M)
        do = tl.load(
            do_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, HEAD_DIM)
        dv = tl.dot(
            tl.trans(p.to(input_dtype)),
            do,
            dv,
            input_precision="ieee" if deterministic else "tf32",
        )

        dp = tl.dot(
            do, tl.trans(v), input_precision="ieee" if deterministic else "tf32"
        )  # (BLOCK_N, HEAD_DIM)x(HEAD_DIM, BLOCK_M) *  -> (BLOCK_N, BLOCK_M)
        delta = tl.load(
            delta_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, 1)
        dqk = p * (dp - delta)  # (BLOCK_N, BLOCK_M)
        dk = tl.dot(
            tl.trans(dqk).to(input_dtype),
            q,
            dk,
            input_precision="ieee" if deterministic else "tf32",
        )

        q_ptr = tl.advance(q_ptr, (BLOCK_N, 0))
        do_ptr = tl.advance(do_ptr, (BLOCK_N, 0))
        bias1_ptr = tl.advance(bias1_ptr, (BLOCK_N, 0))
        bias2_ptr = tl.advance(bias2_ptr, (BLOCK_N, 0))
        m_ptr = tl.advance(m_ptr, (BLOCK_N, 0))
        delta_ptr = tl.advance(delta_ptr, (BLOCK_N, 0))
    tl.store(dv_ptr, dv.to(input_dtype), boundary_check=(0, 1))
    tl.store(dk_ptr, dk.to(input_dtype), boundary_check=(0, 1))


@autotune(
    configs=_attention_bwd_dq_configs,
    key=["H", "HEAD_DIM", "TAG_N"],
)
@triton.jit
def _attention_bwd_dq(
    Q,
    K,
    V,
    Bias1,
    Bias2,
    M,
    Delta,
    DQ,
    DO,
    N,
    S,
    H,
    HEAD_DIM: tl.constexpr,
    TAG_N: tl.constexpr,
    deterministic: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
):

    input_dtype = Q.dtype.element_ty

    bn = tl.program_id(2)
    h = tl.program_id(1)
    s = tl.program_id(0) * BLOCK_M
    b = bn // N
    n = bn % N

    Bs = N * S * H * HEAD_DIM
    Ns = S * H * HEAD_DIM
    Ss = H * HEAD_DIM
    Hs = HEAD_DIM

    qkvo_offset = b * Bs + n * Ns + h * Hs
    bias1_offset = b * N * S + n * S
    bias2_offset = b * S * S * H + h * S * S
    delta_m_offset = b * N * H * S + n * S * H + h

    q_ptr = tl.make_block_ptr(
        base=Q + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(1, 0),
    )

    dq_ptr = tl.make_block_ptr(
        base=DQ + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(1, 0),
    )

    k_ptr = tl.make_block_ptr(
        base=K + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(0, 0),
        block_shape=(BLOCK_N, HEAD_DIM),
        order=(0, 1),
    )

    v_ptr = tl.make_block_ptr(
        base=V + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(0, 0),
        block_shape=(BLOCK_N, HEAD_DIM),
        order=(0, 1),
    )

    do_ptr = tl.make_block_ptr(
        base=DO + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(0, 1),
    )

    bias1_ptr = tl.make_block_ptr(
        base=Bias1 + bias1_offset,
        shape=(S, S),
        strides=(0, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, BLOCK_N),
        order=(1, 0),
    )

    bias2_ptr = tl.make_block_ptr(
        base=Bias2 + bias2_offset,
        shape=(S, S),
        strides=(S, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, BLOCK_N),
        order=(1, 0),
    )

    m_ptr = tl.make_block_ptr(
        base=M + delta_m_offset,
        shape=(S, 1),
        strides=(H, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, 1),
        order=(0, 1),
    )

    delta_ptr = tl.make_block_ptr(
        base=Delta + delta_m_offset,
        shape=(S, 1),
        strides=(H, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, 1),
        order=(0, 1),
    )

    q = tl.load(q_ptr, boundary_check=(0, 1), padding_option="zero")
    do = tl.load(do_ptr, boundary_check=(0, 1), padding_option="zero")
    delta = tl.load(
        delta_ptr, boundary_check=(0, 1), padding_option="zero"
    )  # (BLOCK_M, 1)
    m = tl.load(m_ptr, boundary_check=(0, 1), padding_option="zero")  # (BLOCK_M, 1)
    dq = tl.zeros([BLOCK_M, HEAD_DIM], dtype=tl.float32)  # (BLOCK_M, HEAD_DIM)
    qk_scale = 1.4426950408889634
    for _ in range(tl.cdiv(S, BLOCK_N)):
        k = tl.load(k_ptr, boundary_check=(0, 1), padding_option="zero")
        bias1 = tl.load(bias1_ptr, boundary_check=(0, 1), padding_option="zero")
        bias2 = tl.load(bias2_ptr, boundary_check=(0, 1), padding_option="zero")

        qk = bias1.to(tl.float32) + bias2.to(tl.float32)
        qk = (
            tl.dot(
                q, tl.trans(k), qk, input_precision="ieee" if deterministic else "tf32"
            )
            * qk_scale
        )  # (BLOCK_M, BLOCK_N)
        p = tl.math.exp2(qk - m)  # (BLOCK_M, BLOCK_N)
        v = tl.load(
            v_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, HEAD_DIM)
        dp = tl.dot(
            do, tl.trans(v), input_precision="ieee" if deterministic else "tf32"
        )  # (BLOCK_M, BLOCK_N)

        dqk = p * (dp - delta)  # (BLOCK_M, BLOCK_N)
        dq = tl.dot(
            dqk.to(input_dtype),
            k,
            dq,
            input_precision="ieee" if deterministic else "tf32",
        )

        k_ptr = tl.advance(k_ptr, (BLOCK_N, 0))
        v_ptr = tl.advance(v_ptr, (BLOCK_N, 0))
        bias1_ptr = tl.advance(bias1_ptr, (0, BLOCK_N))
        bias2_ptr = tl.advance(bias2_ptr, (0, BLOCK_N))
    tl.store(dq_ptr, dq.to(input_dtype), boundary_check=(0, 1))


@autotune(
    configs=_attention_bwd_dbias2_configs,
    key=["H", "HEAD_DIM", "TAG_N"],
)
@triton.jit
def _attention_bwd_dbias2(
    Q,
    K,
    V,
    Bias1,
    Bias2,
    M,
    Delta,
    DBias2,
    DO,
    N,
    S,
    H,
    HEAD_DIM: tl.constexpr,
    TAG_N: tl.constexpr,
    deterministic: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
):

    input_dtype = Q.dtype.element_ty

    s1 = tl.program_id(0) * BLOCK_M
    s2 = tl.program_id(1) * BLOCK_N
    bh = tl.program_id(2)
    b = bh // H
    h = bh % H

    Bs = N * S * H * HEAD_DIM
    Ns = S * H * HEAD_DIM
    Ss = H * HEAD_DIM
    Hs = HEAD_DIM

    bias1_offset = b * N * S
    bias2_offset = b * S * S * H + h * S * S

    qkvo_offset = b * Bs + h * Hs

    delta_m_offset = b * N * H * S + h

    bias2_ptr = tl.make_block_ptr(
        base=Bias2 + bias2_offset,
        shape=(S, S),
        strides=(S, 1),
        offsets=(s1, s2),
        block_shape=(BLOCK_M, BLOCK_N),
        order=(1, 0),
    )

    dbias2_ptr = tl.make_block_ptr(
        base=DBias2 + bias2_offset,
        shape=(S, S),
        strides=(S, 1),
        offsets=(s1, s2),
        block_shape=(BLOCK_M, BLOCK_N),
        order=(1, 0),
    )

    bias2 = tl.load(bias2_ptr, boundary_check=(0, 1), padding_option="zero")
    dbias2 = tl.zeros([BLOCK_M, BLOCK_N], dtype=tl.float32)  # (BLOCK_M, BLOCK_N)
    qk_scale = 1.4426950408889634
    for n in range(0, N, 1):
        q_ptr = tl.make_block_ptr(
            base=Q + qkvo_offset + n * Ns,
            shape=(S, HEAD_DIM),
            strides=(Ss, 1),
            offsets=(s1, 0),
            block_shape=(BLOCK_M, HEAD_DIM),
            order=(1, 0),
        )

        k_ptr = tl.make_block_ptr(
            base=K + qkvo_offset + n * Ns,
            shape=(S, HEAD_DIM),
            strides=(Ss, 1),
            offsets=(s2, 0),
            block_shape=(BLOCK_N, HEAD_DIM),
            order=(1, 0),
        )

        v_ptr = tl.make_block_ptr(
            base=V + qkvo_offset + n * Ns,
            shape=(S, HEAD_DIM),
            strides=(Ss, 1),
            offsets=(s2, 0),
            block_shape=(BLOCK_N, HEAD_DIM),
            order=(1, 0),
        )

        do_ptr = tl.make_block_ptr(
            base=DO + qkvo_offset + n * Ns,
            shape=(S, HEAD_DIM),
            strides=(Ss, 1),
            offsets=(s1, 0),
            block_shape=(BLOCK_M, HEAD_DIM),
            order=(1, 0),
        )

        bias1_ptr = tl.make_block_ptr(
            base=Bias1 + bias1_offset + n * S,
            shape=(S, S),
            strides=(0, 1),
            offsets=(s1, s2),
            block_shape=(BLOCK_M, BLOCK_N),
            order=(1, 0),
        )

        m_ptr = tl.make_block_ptr(
            base=M + delta_m_offset + n * S * H,
            shape=(S, 1),
            strides=(H, 1),
            offsets=(s1, 0),
            block_shape=(BLOCK_M, 1),
            order=(0, 1),
        )

        delta_ptr = tl.make_block_ptr(
            base=Delta + delta_m_offset + n * S * H,
            shape=(S, 1),
            strides=(H, 1),
            offsets=(s1, 0),
            block_shape=(BLOCK_M, 1),
            order=(0, 1),
        )

        q = tl.load(q_ptr, boundary_check=(0, 1), padding_option="zero")
        k = tl.load(k_ptr, boundary_check=(0, 1), padding_option="zero")
        v = tl.load(v_ptr, boundary_check=(0, 1), padding_option="zero")
        do = tl.load(do_ptr, boundary_check=(0, 1), padding_option="zero")
        delta = tl.load(delta_ptr, boundary_check=(0, 1), padding_option="zero")
        m = tl.load(m_ptr, boundary_check=(0, 1), padding_option="zero")
        bias1 = tl.load(bias1_ptr, boundary_check=(0, 1), padding_option="zero")

        qk = bias1.to(tl.float32) + bias2.to(tl.float32)
        qk = (
            tl.dot(
                q, tl.trans(k), qk, input_precision="ieee" if deterministic else "tf32"
            )
            * qk_scale
        )  # (BLOCK_M, BLOCK_N)
        p = tl.math.exp2(qk - m)  # (BLOCK_M, BLOCK_N)
        v = tl.load(
            v_ptr, boundary_check=(0, 1), padding_option="zero"
        )  # (BLOCK_N, HEAD_DIM)
        dp = tl.dot(
            do, tl.trans(v), input_precision="ieee" if deterministic else "tf32"
        )  # (BLOCK_M, BLOCK_N)

        dqk = p * (dp - delta)  # (BLOCK_M, BLOCK_N)
        dbias2 += dqk

    tl.store(dbias2_ptr, dbias2.to(input_dtype), boundary_check=(0, 1))
