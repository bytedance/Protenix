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

import torch
import math
import triton
from protenix.model.tri_attention.forward import _attention_fwd
from protenix.model.tri_attention.backward import (
    _attention_bwd_preprocess,
    _attention_bwd_dkdv,
    _attention_bwd_dq,
    _attention_bwd_dbias2,
)


def get_tag(x: int) -> int:
    if x > 4096:
        return 4096
    return math.ceil(x / 32) * 32


class TriAttentionFunction(torch.autograd.Function):

    @staticmethod
    def forward(ctx, Q, K, V, Bias1, Bias2, deterministic=False):
        Q = Q.contiguous()
        K = K.contiguous()
        V = V.contiguous()
        Bias1 = Bias1.contiguous()
        Bias2 = Bias2.contiguous()

        B, N, S, H, D = Q.shape
        TAG_N = get_tag(S)
        O = torch.empty_like(Q)
        M = torch.empty((B, N, S, H), device=Q.device, dtype=torch.float32)

        grid = lambda args: (triton.cdiv(S, args["BLOCK_M"]), H, B * N)
        _attention_fwd[grid](
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
            HEAD_DIM=D,
            TAG_N=TAG_N,
            deterministic=deterministic,
        )

        ctx.save_for_backward(Q, K, V, Bias1, Bias2, O, M)
        ctx.deterministic = deterministic
        return O

    def backward(ctx, DO):
        Q, K, V, Bias1, Bias2, O, M = ctx.saved_tensors
        deterministic = ctx.deterministic
        B, N, S, H, D = Q.shape
        TAG_N = get_tag(S)
        DQ = torch.empty_like(Q)
        DK = torch.empty_like(K)
        DV = torch.empty_like(V)
        DBias1 = torch.empty_like(Bias1)
        DBias2 = torch.empty_like(Bias2)
        Delta = torch.empty_like(M)
        grid = lambda args: (triton.cdiv(S, args["BLOCK_M"]), H, B * N)
        _attention_bwd_preprocess[grid](O, DO, Delta, N, S, H, HEAD_DIM=D, TAG_N=TAG_N)

        grid = lambda args: (triton.cdiv(S, args["BLOCK_M"]), H, B * N)
        _attention_bwd_dkdv[grid](
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
            HEAD_DIM=D,
            TAG_N=TAG_N,
            deterministic=deterministic,
        )

        grid = lambda args: (triton.cdiv(S, args["BLOCK_M"]), H, B * N)
        _attention_bwd_dq[grid](
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
            HEAD_DIM=D,
            TAG_N=TAG_N,
            deterministic=deterministic,
        )

        grid = lambda args: (
            triton.cdiv(S, args["BLOCK_M"]),
            triton.cdiv(S, args["BLOCK_N"]),
            B * H,
        )
        _attention_bwd_dbias2[grid](
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
            HEAD_DIM=D,
            TAG_N=TAG_N,
            deterministic=deterministic,
        )
        return DQ, DK, DV, DBias1, DBias2, None


class TriAttention(torch.nn.Module):
    def __init__(self, deterministic=False):
        super(TriAttention, self).__init__()
        self.deterministic = deterministic

    def forward(self, Q, K, V, Bias1, Bias2):
        """
        Args:
            Q: (B, S, S, H, D)
            K: (B, S, S, H, D)
            V: (B, S, S, H, D)
            Bias1: (B, S, 1, 1, S)
            Bias2: (B, 1, H, S, S)
        Returns:
            output: (B, S, S, H, D)
        """
        output = TriAttentionFunction.apply(Q, K, V, Bias1, Bias2, self.deterministic)
        return output
