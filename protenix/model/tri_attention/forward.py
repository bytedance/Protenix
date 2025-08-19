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
from protenix.model.tri_attention.autotune_helpers import _attention_fwd_configs


@autotune(
    configs=_attention_fwd_configs,
    key=["H", "HEAD_DIM", "TAG_N"],
)
@triton.jit
def _attention_fwd(
    Q,
    K,
    V,
    O,
    Bias1,
    Bias2,
    M,
    N,
    S,
    H,
    HEAD_DIM: tl.constexpr,
    TAG_N: tl.constexpr,
    deterministic: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
):
    """
    B: batch
    N: number
    S: sequence, equal to N
    H: number of head
    HEAD_DIM: dimension of head

    Q: (B, N, S, H, HEAD_DIM)
    K: (B, N, S, H, HEAD_DIM)
    V: (B, N, S, H, HEAD_DIM)
    O: (B, N, S, H, HEAD_DIM)
    Bias1: (B, N, 1, 1, S), Bias1 is mask which doesn't need grad
    Bias2: (B, 1, H, S, S)

    M: (B, N, S, H)
    """
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
    m_offset = b * N * H * S + n * S * H + h

    q_ptr = tl.make_block_ptr(
        base=Q + qkvo_offset,
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
        order=(1, 0),
    )

    o_ptr = tl.make_block_ptr(
        base=O + qkvo_offset,
        shape=(S, HEAD_DIM),
        strides=(Ss, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, HEAD_DIM),
        order=(1, 0),
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
        base=M + m_offset,
        shape=(S, 1),
        strides=(H, 1),
        offsets=(s, 0),
        block_shape=(BLOCK_M, 1),
        order=(0, 1),
    )

    qk_scale = 1.4426950408889634

    m_i = tl.zeros([BLOCK_M], dtype=tl.float32) - float("inf")
    l_i = tl.zeros([BLOCK_M], dtype=tl.float32) + 1
    acc = tl.zeros([BLOCK_M, HEAD_DIM], dtype=tl.float32)

    q = tl.load(q_ptr, boundary_check=(0, 1), padding_option="zero")

    offset_q = (s + tl.arange(0, BLOCK_M))[:, None] < S
    offset_k = tl.arange(0, BLOCK_N)[None, :]
    for _ in range(tl.cdiv(S, BLOCK_N)):
        k = tl.load(k_ptr, boundary_check=(0, 1), padding_option="zero")
        bias1 = tl.load(bias1_ptr, boundary_check=(0, 1), padding_option="zero")
        bias2 = tl.load(bias2_ptr, boundary_check=(0, 1), padding_option="zero")
        mask = offset_q * (offset_k < S)
        qk = bias1.to(tl.float32) + bias2.to(tl.float32)
        qk = (
            tl.dot(
                q, tl.trans(k), qk, input_precision="ieee" if deterministic else "tf32"
            )
            * qk_scale
        )
        qk = tl.where(mask, qk, -1e6)
        m_ij = tl.maximum(tl.max(qk, 1), m_i)
        qk = qk - m_ij[:, None]
        alpha = tl.math.exp2(m_i - m_ij)
        p = tl.math.exp2(qk)
        l_ij = l_i * alpha + tl.sum(p, 1)
        acc = acc * alpha[:, None]

        v = tl.load(v_ptr, boundary_check=(0, 1), padding_option="zero")
        p = p.to(input_dtype)
        acc = tl.dot(p, v, acc, input_precision="ieee" if deterministic else "tf32")

        l_i = l_ij
        m_i = m_ij

        k_ptr = tl.advance(k_ptr, (BLOCK_N, 0))
        v_ptr = tl.advance(v_ptr, (BLOCK_N, 0))
        bias1_ptr = tl.advance(bias1_ptr, (0, BLOCK_N))
        bias2_ptr = tl.advance(bias2_ptr, (0, BLOCK_N))
        offset_k += BLOCK_N

    acc = acc / l_i[:, None]
    m_i += tl.math.log2(l_i)
    acc = acc.to(input_dtype)
    tl.store(o_ptr, acc, boundary_check=(0, 1))
    tl.store(m_ptr, m_i[:, None], boundary_check=(0, 1))
