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

import math
from functools import partial
from typing import Optional, Union

import torch
import torch.nn as nn
import torch.nn.functional as F

from protenix.model.utils import (
    flatten_final_dims,
    move_final_dim_to_dim,
    pad_at_dim,
    reshape_at_dim,
)
from protenix.openfold_local.model.primitives import LayerNorm, trunc_normal_init_
from protenix.openfold_local.utils.chunk_utils import chunk_layer


class Linear(nn.Linear):

    def __init__(
        self,
        in_features: int,
        out_features: int,
        bias: bool = True,
        device: torch.device = None,
        dtype: torch.dtype = None,
        precision: torch.dtype = None,
        initializer: str = "default",
    ):
        """Linear module with customized initialization.

        Args:
            precision (torch.dtype, optional): _description_. Defaults to None.
            initializer (str, optional): initializer: choose one from ['default', 'relu', 'zeros']. Defaults to "default".
        """

        self.use_bias = bias
        self.precision = precision
        self.initializer = initializer
        super().__init__(
            in_features=in_features,
            out_features=out_features,
            bias=bias,
            device=device,
            dtype=dtype,
        )

        self._init_params()

    @torch.no_grad()
    def _init_params(self):
        if self.use_bias:
            nn.init.zeros_(self.bias)  # zero-init bias

        if self.initializer == "default":
            trunc_normal_init_(self.weight, scale=1.0)
        elif self.initializer == "relu":
            trunc_normal_init_(self.weight, scale=2.0)
        elif self.initializer == "zeros":
            nn.init.zeros_(self.weight)
        else:
            raise ValueError(f"Invalid initializer: {self.initializer}.")

    def forward(self, input: torch.Tensor) -> torch.Tensor:
        if self.precision is not None:
            input_dtype = input.dtype
            with torch.cuda.amp.autocast(enabled=False):
                bias = (
                    self.bias.to(dtype=self.precision)
                    if self.bias is not None
                    else None
                )
                return F.linear(
                    input.to(dtype=self.precision),
                    self.weight.to(dtype=self.precision),
                    bias,
                ).to(dtype=input_dtype)
        else:
            return F.linear(input, self.weight, self.bias)


LinearNoBias = partial(Linear, bias=False)


class AdaptiveLayerNorm(nn.Module):
    """
    Implements Algorithm 26 in AF3
    """

    def __init__(self, c_a: int = 768, c_s: int = 384) -> None:
        """
        Args:
            c_a (int, optional): the embedding dim of a(single feature aggregated atom info). Defaults to 768.
            c_s (int, optional):  hidden dim [for single embedding]. Defaults to 384.
        """
        super(AdaptiveLayerNorm, self).__init__()
        self.layernorm_a = LayerNorm(c_a, create_scale=False, create_offset=False)
        self.layernorm_s = LayerNorm(c_s, create_offset=False)
        self.linear_s = Linear(in_features=c_s, out_features=c_a, initializer="zeros")
        self.linear_nobias_s = LinearNoBias(
            in_features=c_s, out_features=c_a, initializer="zeros"
        )

    def forward(self, a: torch.Tensor, s: torch.Tensor) -> torch.Tensor:
        """
        Args:
            a (torch.Tensor): the single feature aggregate per-atom representation
                [..., N_token, c_a]
            s (torch.Tensor): single embedding
                [..., N_token, c_s]

        Returns:
            torch.Tensor: the updated a from AdaLN
                [..., N_token, c_a]
        """
        a = self.layernorm_a(a)
        s = self.layernorm_s(s)
        a = torch.sigmoid(self.linear_s(s)) * a + self.linear_nobias_s(s)
        return a


class BiasInitLinear(Linear):
    """Support biasinit for nn.Linear Called just like torch.nn.Linear."""

    def __init__(
        self,
        in_features: int,
        out_features: int,
        bias: bool = True,
        biasinit: float = 0.0,
        **kwargs,
    ) -> None:
        """
        Args:
            in_features (int): in_features
            out_features (int): out_features
            bias (bool, optional): whether add bias. Defaults to True.
            biasinit (float, optional): the initial bias value. Defaults to 0.0.
        """
        super(BiasInitLinear, self).__init__(
            in_features=in_features, out_features=out_features, bias=bias, **kwargs
        )
        nn.init.zeros_(tensor=self.weight)
        if bias:
            nn.init.constant_(tensor=self.bias, val=biasinit)


class Transition(nn.Module):
    """
    Implements Algorithm 11 in AF3
    """

    def __init__(self, c_in: int, n: int) -> None:
        """
        Args:
            c_in (int, optional): the input dimension.
            n (int, optional): factor by which c_in is multiplied to obtain hidden dimension.
        """
        super(Transition, self).__init__()
        self.n = n
        self.c_in = c_in
        self.layernorm1 = LayerNorm(c_in)
        self.linear_no_bias_a = LinearNoBias(
            in_features=c_in, out_features=n * c_in, initializer="relu"
        )
        self.linear_no_bias_b = LinearNoBias(
            in_features=c_in, out_features=n * c_in, initializer="relu"
        )
        self.linear_no_bias = LinearNoBias(
            in_features=n * c_in, out_features=c_in, initializer="zeros"
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x (torch.Tensor): the input tensor
                [..., c]

        Returns:
            torch.Tensor: the output tensor as the same shape of x
                [..., c]
        """
        if self.training:
            x = self.layernorm1(x)
            a = self.linear_no_bias_a(x)
            b = self.linear_no_bias_b(x)
            x = self.linear_no_bias(F.silu(a) * b)
            return x
        else:
            other_dims = x.shape[:-1]
            dim_size = x.shape[-1]
            size = x.shape[-2]
            x = x.reshape(-1, dim_size)
            chunk_num = 1 if size < 3200 else 8
            chunks = torch.chunk(x, chunk_num, dim=-2)
            outputs = torch.empty(
                (x.shape[0], self.c_in), dtype=x.dtype, device=x.device
            )
            start = 0
            for chunk in chunks:
                y = self.layernorm1(chunk)
                a = self.linear_no_bias_a(y)
                a = F.silu(a, True)
                b = self.linear_no_bias_b(y)
                del y
                b *= a
                del a
                b = self.linear_no_bias(b)
                outputs[start : start + b.shape[0]] = b
                start += b.shape[0]
                del b
            outputs = outputs.reshape(*other_dims, self.c_in)
            return outputs


def _attention(
    q: torch.Tensor,
    k: torch.Tensor,
    v: torch.Tensor,
    attn_bias: Optional[torch.Tensor] = None,
    use_efficient_implementation: bool = False,
    inplace_safe: bool = False,
) -> torch.Tensor:
    """Attention.

    Args:
        q (torch.Tensor): query tensor of shape [..., n_q, d]
        k (torch.Tensor): key tensor of shape [..., n_kv, d]
        v (torch.Tensor): value tensor of shape[..., n_kv, d]
        attn_bias (torch.Tensor, optional): attention bias tensor of shape [..., n_q, n_kv]. Defaults to None.
        use_efficient_implementation (bool): whether to use the torch.nn.functional.scaled_dot_product_attention, Defaults to False.

    Returns:
        torch.Tensor: output of tensor [..., n_q, d]
    """
    assert k.shape == v.shape

    # Upcast to compute attn_weights
    input_dtype = q.dtype
    q = q.to(dtype=torch.float32)
    k = k.to(dtype=torch.float32)
    if attn_bias is not None:
        attn_bias = attn_bias.to(dtype=torch.float32)

    if use_efficient_implementation:
        attn_output = F.scaled_dot_product_attention(
            query=q,
            key=k,
            value=v,
            attn_mask=attn_bias,
            scale=1.0,
        )
        return attn_output

    with torch.cuda.amp.autocast(enabled=False):
        # [..., n_kv, d] -> [..., d, n_kv]
        k = k.transpose(-1, -2)

        # [..., n_q, d], [..., d, n_kv] -> [..., n_q, n_kv]
        attn_weights = q @ k

        if attn_bias is not None:
            if inplace_safe:
                attn_weights += attn_bias
            else:
                attn_weights = attn_weights + attn_bias

        # [..., n_q, n_kv]
        attn_weights = F.softmax(attn_weights, dim=-1)

    # [..., n_q, n_kv], [..., n_kv, d] -> [..., n_q, d]
    attn_output = attn_weights.to(dtype=input_dtype) @ v

    return attn_output


def rearrange_qk_to_dense_trunk(
    q: Union[torch.Tensor, list[torch.Tensor]],
    k: Union[torch.Tensor, list[torch.Tensor]],
    dim_q: Union[int, list[int]],
    dim_k: Union[int, list[int]],
    n_queries: int = 32,
    n_keys: int = 128,
    compute_mask: bool = True,
) -> tuple[Union[torch.Tensor, list[torch.Tensor]]]:
    """Rearrange q/k into blocked tensors for local operations.

    Args:
        q (torch.Tensor): query tensor. Could be a tensor or a list of tensors.
            [..., n_q, ...] (n_q is at dimension dim_q)
        k (torch.Tensor | List[torch.Tensor]): key tensor. Could be a tensor or a list of tensors.
            [..., n_k, ...] (n_k is at dimension dim_k)
        dim_q (int): along which dimension to build the trunks. Could be an int or a list of int.
        dim_k (int): along which dimension to build the trunks. Could be an int or a list of int.
        n_queries (int, optional): local window size of query tensor.
        n_keys (int, optional): local window size of key/value tensor.

    Returns:
        tuple[Union[torch.Tensor, list[torch.Tensor]]]:
            q_trunked: torch.Tensor or list of tensors. Same as the input type.
                [..., n_trunks, n_queries, ...]
            k_trunked: torch.Tensor or list of tensors. Same as the input type.
                [..., n_trunks, n_keys, ...]
            padding_info (dict):
                mask_trunked: torch.Tensor
                    [n_trunks, n_queries, n_keys]
                q_pad: query padded dimension
    """

    assert n_keys >= n_queries
    assert n_queries & 0x01 == 0
    assert n_keys & 0x01 == 0

    def basic_checks(x, dim_x):
        if isinstance(x, list):
            x_is_list = True
            assert isinstance(dim_x, list)
        else:
            x_is_list = False
            x = [x]
            dim_x = [dim_x]
        n_x = x[0].size(dim_x[0])
        for i in range(len(dim_x)):
            if dim_x[i] < 0:
                dim_x[i] = len(x[i].shape) + dim_x[i]
            assert x[i].size(dim_x[i]) == n_x
        return x, dim_x, x_is_list, n_x, len(x)

    q, dim_q, q_is_list, n, num_q = basic_checks(q, dim_q)
    k, dim_k, k_is_list, n_k, num_k = basic_checks(k, dim_k)

    assert n == n_k
    n_trunks = int(math.ceil(n / n_queries))
    q_pad_length = n_trunks * n_queries - n

    q_new = [
        pad_at_dim(q[i], dim=dim_q[i], pad_length=(0, q_pad_length))
        for i in range(num_q)
    ]
    q_trunked = [
        reshape_at_dim(q_new[i], dim=dim_q[i], target_shape=(n_trunks, n_queries))
        for i in range(num_q)
    ]

    pad_left = (n_keys - n_queries) // 2
    pad_right = int((n_trunks - 1 / 2) * n_queries + n_keys / 2 - n + 1 / 2)

    k_new = [
        pad_at_dim(k[i], dim=dim_k[i], pad_length=(pad_left, pad_right))
        for i in range(num_k)
    ]
    k_trunked = [
        k_new[i].unfold(dim_k[i], size=n_keys, step=n_queries) for i in range(num_k)
    ]
    k_trunked = [
        move_final_dim_to_dim(k_trunked[i], dim=dim_k[i] + 1) for i in range(num_k)
    ]

    if compute_mask:
        pad_mask = q[0].new_ones(
            *(1,) * len(q[0].shape[:-2]),
            n + q_pad_length,
            n + pad_left + pad_right,
            requires_grad=False,
        )
        pad_mask[..., :n, 0:pad_left] = 0
        pad_mask[..., :n, pad_left + n : :] = 0
        pad_mask[..., n::, :] = 0

        concat_split_data = optimized_concat_split(pad_mask, n_queries)
        pad_mask_trunked = (
            concat_split_data.unfold(
                -1, n_keys, pad_mask.size(-1) + n_queries
            ).transpose(-2, -3)
        ).bool()
    else:
        pad_mask_trunked = None

    if not q_is_list:
        q_trunked = q_trunked[0]
    if not k_is_list:
        k_trunked = k_trunked[0]

    padding_info = {
        "mask_trunked": pad_mask_trunked,
        "q_pad": q_pad_length,
        "k_pad_left": pad_left,
        "k_pad_right": pad_right,
    }

    return q_trunked, k_trunked, padding_info


def optimized_concat_split(attn_bias: torch.Tensor, n_queries: int) -> torch.Tensor:
    """Optimized concatenation and splitting of attention bias tensor.

    Args:
        attn_bias (torch.Tensor): The attention bias tensor.
            Shape: [..., D, E]
        n_queries (int): The number of queries in each split.

    Returns:
        torch.Tensor: The reshaped and permuted attention bias tensor.
            Shape: [..., n_queries, D // n_queries * E]
    """
    D = attn_bias.size(-2)
    E = attn_bias.size(-1)
    assert D % n_queries == 0
    num_splits = D // n_queries
    reshaped = attn_bias.reshape(*attn_bias.shape[:-2], num_splits, n_queries, E)
    permuted = reshaped.permute(*range(reshaped.dim() - 3), -2, -3, -1)
    output = permuted.reshape(*attn_bias.shape[:-2], n_queries, num_splits * E)
    return output


def rearrange_to_dense_trunk(
    q: torch.Tensor,
    k: torch.Tensor,
    v: torch.Tensor,
    n_queries: int,
    n_keys: int,
    attn_bias: Optional[torch.Tensor] = None,
    inf: float = 1e10,
) -> tuple[Union[torch.Tensor, int]]:
    """Rearrange q/k/v/bias into blocked tensors for local attention.

    Args:
        q (torch.Tensor): query tensor
            [..., n_q, d]
        k (torch.Tensor): key tensor
            [..., n_kv, d]
        v (torch.Tensor): value tensor
            [..., n_kv, d]
        attn_bias (torch.Tensor, optional): attention bias
            [..., n_q, n_kv] or None
        n_queries (int, optional): local window size of query tensor.
        n_keys (int, optional): local window size of key/value tensor.
        inf (float, optional): used for attention masking. Defaults to 1e10.

    Returns:
        tuple[Union[torch.Tensor, int]]:
            q_trunked
                [..., n_trunks, n_queries, d]
            k_trunked / v_trunked
                [..., n_trunks, n_keys, d]
            attn_bias_trunked:  padded position filled with -inf
                [..., n_trunks, n_queries, n_keys]
            q_pad_length: query padded dimension
    """
    assert n_keys >= n_queries
    assert n_queries & 0x01 == 0
    assert n_keys & 0x01 == 0

    n, d = q.shape[-2:]

    q_trunked, kv_trunked, padding_info = rearrange_qk_to_dense_trunk(
        q=q,
        k=[k, v],
        dim_q=-2,
        dim_k=[-2, -2],
        n_queries=n_queries,
        n_keys=n_keys,
        compute_mask=False,
    )
    q_pad_length, pad_left, pad_right = (
        padding_info["q_pad"],
        padding_info["k_pad_left"],
        padding_info["k_pad_right"],
    )

    # Padded_width = n + pad_left + pad_right
    if attn_bias is None:
        attn_bias = q.new_zeros(
            *(1,) * len(q.shape[:-2]), n + q_pad_length, n + pad_left + pad_right
        )
        attn_bias[..., :n, 0:pad_left] = -inf
        attn_bias[..., :n, pad_left + n : :] = -inf
        attn_bias[..., n::, :] = -inf
    else:
        attn_bias = F.pad(attn_bias, (pad_left, pad_right, 0, q_pad_length), value=-inf)

    concat_split_data = optimized_concat_split(attn_bias, n_queries)
    attn_bias_trunked = concat_split_data.unfold(
        -1, n_keys, attn_bias.shape[-1] + n_queries
    ).transpose(-2, -3)
    return q_trunked, kv_trunked[0], kv_trunked[1], attn_bias_trunked, q_pad_length


def _local_attention(
    q: torch.Tensor,
    k: torch.Tensor,
    v: torch.Tensor,
    n_queries: int,
    n_keys: int,
    attn_bias: Optional[torch.Tensor] = None,
    trunked_attn_bias: Optional[torch.Tensor] = None,
    inf: float = 1e10,
    use_efficient_implementation: bool = False,
    inplace_safe: bool = False,
    chunk_size: Optional[int] = None,
) -> torch.Tensor:
    """Local attention

    Args:
        q (torch.Tensor): query tensor
            [..., Q, d]
        k (torch.Tensor): key tensor
            [..., K, d]
        v (torch.Tensor): value tensor
            [..., K, d]
        n_queries (int): local window size of query.
        n_keys (int): local window size of key/value.
        attn_bias (torch.Tensor, optional): the input biases for attention. Defaults to None.
            [..., Q, K]
        trunked_attn_bias (torch.Tensor, optional): the input biases where shape has been rearranged to dense trunks. Defaults to None.
            [..., n_trunks, n_queries, n_keys]
        inf (float): inf number used for attention bias. Defaults to 1e10.
        use_efficient_implementation (bool): whether to use the torch.nn.functional.scaled_dot_product_attention, Defaults to False.
    Returns:
        torch.Tensor: standard attention output
            [..., Q, d]
    """
    assert q.shape == k.shape == v.shape  # local attention doesn't make sense if Q != K

    # Prepare for attention qkv, q: [..., n_trunks, n_queries, d], kv: [..., n_trunks, n_keys, d]

    # Rerrange to dense trunks
    # q: [*, n, d] -> [*, n_trunks, n_queries, d]
    # kv: [*, n, d] -> [*, n_trunks, n_keys, d]
    # attn_bias: [*, n, d] -> [*, n_trunks, n_queries, n_keys]
    q_trunked, k_trunked, v_trunked, attn_bias_trunked, q_pad_length = (
        rearrange_to_dense_trunk(
            q=q,
            k=k,
            v=v,
            n_queries=n_queries,
            n_keys=n_keys,
            attn_bias=attn_bias,
            inf=inf,
        )
    )

    # Apply attention
    # [..., n_trunks, n_queries, d]
    if trunked_attn_bias is not None:
        attn_bias_trunked = attn_bias_trunked + trunked_attn_bias

    if chunk_size is not None:
        attn_inputs = {
            "q": q_trunked,
            "k": k_trunked,
            "v": v_trunked,
            "attn_bias": attn_bias_trunked,
        }
        out = chunk_layer(
            partial(
                _attention,
                use_efficient_implementation=use_efficient_implementation,
                inplace_safe=inplace_safe,
            ),
            attn_inputs,
            chunk_size=chunk_size,
            no_batch_dims=len(attn_bias_trunked.shape[:-2]),
            _out=None,
        )
    else:
        out = _attention(
            q=q_trunked,
            k=k_trunked,
            v=v_trunked,
            attn_bias=attn_bias_trunked,
            use_efficient_implementation=use_efficient_implementation,
            inplace_safe=inplace_safe,
        )

    # Revert back to orignal shape and remove q_pad_length
    # [..., n_trunks, n_queries, d] ->  [..., n_trunks * n_queries, d] ->  [..., n, d]
    out = out.reshape(*out.shape[:-3], -1, out.shape[-1])
    if q_pad_length > 0:
        out = out[..., :-q_pad_length, :]
    return out


def create_local_attn_bias(
    n: int, n_queries: int, n_keys: int, inf: float = 1e10, device: torch.device = None
) -> torch.Tensor:
    """Create local attention bias based on query window n_queries and kv window n_keys.

    Args:
        n (int): the length of quiries
        n_queries (int): window size of quiries
        n_keys (int): window size of keys/values
        inf (float, optional): the inf to mask attention. Defaults to 1e10.
        device (torch.device, optional): cuda|cpu|None. Defaults to None.

    Returns:
        torch.Tensor: the diagonal-like global attention bias
    """
    n_trunks = int(math.ceil(n / n_queries))
    padded_n = n_trunks * n_queries
    attn_mask = torch.zeros(padded_n, padded_n, device=device)
    for block_index in range(0, n_trunks):
        i = block_index * n_queries
        j1 = max(0, n_queries * block_index - (n_keys - n_queries) // 2)
        j2 = n_queries * block_index + (n_queries + n_keys) // 2
        attn_mask[i : i + n_queries, j1:j2] = 1.0
    attn_bias = (1 - attn_mask) * -inf
    return attn_bias.to(device=device)[:n, :n]


class Attention(nn.Module):
    """Standard multi-head attention
    Ref to openfold:
    https://github.com/aqlaboratory/openfold/blob/feb45a521e11af1db241a33d58fb175e207f8ce0/openfold/model/primitives.py#L340
    """

    def __init__(
        self,
        c_q: int,
        c_k: int,
        c_v: int,
        c_hidden: int,
        num_heads: int,
        gating: bool = True,
        q_linear_bias: bool = True,
        local_attention_method: str = "global_attention_with_bias",
        use_efficient_implementation: bool = False,
        zero_init: bool = True,
    ) -> None:
        """

        Args:
            c_q (int): Input dimension of query data
            c_k (int): Input dimension of key data
            c_v (int): Input dimension of value data
            c_hidden (int): Per-head hidden dimension
            num_heads (int): Number of attention heads
            gating (bool, optional): Whether the output should be gated using query data. Defaults to True.
            q_linear_bias (bool, optional): whether use Linear with bias as in AF3. Defaults to False.
            local_attention_method (str, optional): local attention method, options:
              - global_attention_with_bias: use full size global attention with sparse attention bias
              - local_cross_attention: use local cross attention to minimize computation
            use_efficient_implementation (bool): whether to use the torch.nn.functional.scaled_dot_product_attention, Defaults to False.

        Notes:
            if use_efficient_implementation == True, torch.nn.functional.scaled_dot_product_attention will
            be used to compute attention efficiently
            There are currently three supported implementations of scaled dot product attention:
                1. FlashAttention: Fast and Memory-Efficient Exact Attention with IO-Awareness

                2. Memory-Efficient Attention

                3. A PyTorch implementation defined in C++ matching the above formulation

            The function may call optimized kernels for improved performance when using the CUDA backend.
            For all other backends, the PyTorch implementation will be used.All implementations are enabled by default.
            Scaled dot product attention attempts to automatically select the most optimal implementation based on the inputs.
        """
        super(Attention, self).__init__()
        self.c_q = c_q
        self.c_k = c_k
        self.c_v = c_v
        self.c_hidden = c_hidden
        self.num_heads = num_heads
        self.gating = gating
        self.local_attention_method = local_attention_method
        self.use_efficient_implementation = use_efficient_implementation

        # DISCREPANCY: c_hidden is not the per-head channel dimension, as
        # stated in the supplement, but the overall channel dimension.
        if q_linear_bias:
            # Attention in AF3
            self.linear_q = Linear(
                in_features=self.c_q, out_features=self.c_hidden * self.num_heads
            )
        else:
            # Vanilla attention
            self.linear_q = LinearNoBias(self.c_q, self.c_hidden * self.num_heads)
        self.linear_k = LinearNoBias(self.c_k, self.c_hidden * self.num_heads)
        self.linear_v = LinearNoBias(self.c_v, self.c_hidden * self.num_heads)
        self.linear_o = LinearNoBias(self.c_hidden * self.num_heads, self.c_q)
        self.linear_g = None
        if self.gating:
            self.linear_g = LinearNoBias(
                self.c_q, self.c_hidden * self.num_heads, initializer="zeros"
            )
            self.sigmoid = nn.Sigmoid()

        self.zero_init = zero_init
        if self.zero_init:
            # zero init the output layer
            nn.init.zeros_(self.linear_o.weight)

    def _prep_qkv(
        self, q_x: torch.Tensor, kv_x: torch.Tensor, apply_scale: bool = True
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Prepare qkv

        Args:
            q_x (torch.Tensor): the input x for q
                [..., c_q]
            kv_x (torch.Tensor): the input x for kv
                [..., c_k]
                [..., c_v]
            apply_scale (bool, optional): apply scale to dot product qk. Defaults to True.

        Returns:
            tuple[torch.Tensor, torch.Tensor, torch.Tensor]: the return q/k/v
                # [..., H, Q/K/V, C_hidden]
        """
        # [*, Q/K/V, H * C_hidden]
        q = self.linear_q(q_x)
        k = self.linear_k(kv_x)
        v = self.linear_v(kv_x)

        # [*, Q/K/V, H, C_hidden]
        q = q.view(q.shape[:-1] + (self.num_heads, -1))
        k = k.view(k.shape[:-1] + (self.num_heads, -1))
        v = v.view(v.shape[:-1] + (self.num_heads, -1))

        # [*, H, Q/K/V, C_hidden]
        q = q.transpose(-2, -3)
        k = k.transpose(-2, -3)
        v = v.transpose(-2, -3)

        if apply_scale:
            q = q / math.sqrt(self.c_hidden)

        return q, k, v

    def _wrap_up(self, o: torch.Tensor, q_x: torch.Tensor) -> torch.Tensor:
        """

        Args:
            o (torch.Tensor): the output of attention
                [..., G/Q, H, C_hidden]
            q_x (torch.Tensor): the input for gated g
                [..., Q, c_q]

        Returns:
            torch.Tensor: the output of attention
        """
        if self.linear_g is not None:
            g = self.sigmoid(self.linear_g(q_x))

            # [*, G/Q, H, C_hidden]
            g = g.view(g.shape[:-1] + (self.num_heads, -1))
            o = o * g

        # [*, Q, H * C_hidden]
        o = flatten_final_dims(o, num_dims=2)

        # [*, Q, C_q]
        o = self.linear_o(o)

        return o

    def forward(
        self,
        q_x: torch.Tensor,
        kv_x: torch.Tensor,
        attn_bias: Optional[torch.Tensor] = None,
        trunked_attn_bias: Optional[torch.Tensor] = None,
        n_queries: Optional[int] = None,
        n_keys: Optional[int] = None,
        inf: Optional[float] = 1e10,
        inplace_safe: bool = False,
        chunk_size: Optional[int] = None,
    ) -> torch.Tensor:
        """

        Args:
            q_x (torch.Tensor): the input x for q
                [..., Q, C_q]
            kv_x (torch.Tensor): the input x for k/v
                [..., K, C_k]
            attn_bias (torch.Tensor, optional): the input biases for attention. Defaults to None.
                [..., H, Q, K] or [..., Q, K]
            trunked_attn_bias (torch.Tensor, optional): the input biases where shape has been rearranged to dense trunks. Defaults to None.
                [..., H, n_trunks, n_queries, n_keys] or [..., n_trunks, n_queries, n_keys]
            n_queries (int, optional): local window size of query tensor. If not None, will perform local attention. Defaults to None.
            n_keys (int, optional): local window size of key tensor. Defaults to None.

        Returns:
            torch.Tensor: attention update
                [*, Q, C_q]
        """

        q, k, v = self._prep_qkv(q_x=q_x, kv_x=kv_x, apply_scale=True)

        if attn_bias is not None:
            if len(attn_bias.shape) == len(q.shape):
                assert attn_bias.shape[:-2] == q.shape[:-2]
            else:
                assert len(attn_bias.shape) == len(q.shape) - 1
                assert attn_bias.shape[:-2] == q.shape[:-3]
                # Expand at head dim, got shape [..., 1, Q, K]
                attn_bias = attn_bias.unsqueeze(dim=-3)

        if trunked_attn_bias is not None:
            # NOTE: trunked_attn_bias can only be used with "local_cross_attention" method
            assert n_queries and n_keys
            assert self.local_attention_method == "local_cross_attention"

            if len(trunked_attn_bias.shape) == len(q.shape) + 1:
                assert trunked_attn_bias.shape[:-3] == q.shape[:-2]
            else:
                assert len(trunked_attn_bias.shape) == len(q.shape)
                # Expand at head dim, got shape [..., 1, n_trunks, n_queries, n_keys]
                trunked_attn_bias = trunked_attn_bias.unsqueeze(dim=-4)

        if n_queries and n_keys:
            if self.local_attention_method == "global_attention_with_bias":
                local_attn_bias = create_local_attn_bias(
                    q.shape[-2], n_queries, n_keys, inf=inf, device=q.device
                )
                # Expand to same shape as attn_bias
                local_attn_bias = local_attn_bias.reshape(
                    (1,) * (len(q.shape[:-2])) + local_attn_bias.shape
                )
                if attn_bias is not None:
                    if inplace_safe:
                        local_attn_bias += attn_bias
                    else:
                        local_attn_bias = local_attn_bias + attn_bias
                o = _attention(
                    q=q,
                    k=k,
                    v=v,
                    attn_bias=local_attn_bias,
                    use_efficient_implementation=self.use_efficient_implementation,
                    inplace_safe=inplace_safe,
                )

            elif self.local_attention_method == "local_cross_attention":
                o = _local_attention(
                    q=q,
                    k=k,
                    v=v,
                    n_queries=n_queries,
                    n_keys=n_keys,
                    attn_bias=attn_bias,
                    trunked_attn_bias=trunked_attn_bias,
                    inf=inf,
                    use_efficient_implementation=self.use_efficient_implementation,
                    inplace_safe=inplace_safe,
                    chunk_size=chunk_size,
                )
            else:
                raise ValueError(
                    f"Invalid local attention method: {self.local_attention_method}"
                )
        else:
            o = _attention(
                q=q,
                k=k,
                v=v,
                attn_bias=attn_bias,
                use_efficient_implementation=self.use_efficient_implementation,
                inplace_safe=inplace_safe,
            )  # [*, H, Q, C_hidden]
        o = o.transpose(-2, -3)  # o: [*, Q, H, C_hidden]
        o = self._wrap_up(o, q_x)  # q_x: [*, Q, c_q]

        return o


def gather_pair_embedding_in_dense_trunk(
    x: torch.Tensor, idx_q: torch.Tensor, idx_k: torch.Tensor
):
    """
    Selectively gather elements from a tensor using two sets of indices.

        x: [..., N_token, N_token, d]
        idx_q: [N_b, N_q]
        idx_k: [N_b, N_k]

    Return:
        y: [..., N_b, N_q, N_k, d]
            where y[..., b, i, j, :] = x[..., idx_q[b, i], idx_k[b, j], :]
    """
    idx_q = idx_q.long()
    idx_k = idx_k.long()
    assert len(idx_q.shape) == len(idx_k.shape) == 2

    # Get the shape parameters
    N_b, N_q = idx_q.shape
    N_k = idx_k.shape[1]

    # Expand idx_q and idx_k to match the shape required for advanced indexing
    idx_q_expanded = idx_q.unsqueeze(-1).expand(-1, -1, N_k)
    idx_k_expanded = idx_k.unsqueeze(1).expand(-1, N_q, -1)

    # Use advanced indexing to gather the desired elements
    y = x[..., idx_q_expanded, idx_k_expanded, :]

    return y


def broadcast_token_to_local_atom_pair(
    z_token: torch.Tensor,
    atom_to_token_idx: torch.Tensor,
    n_queries: int,
    n_keys: int,
    compute_mask: bool = True,
) -> torch.Tensor:
    """Broadcast token pair embedding to atom pair embedding

    Args:
        z_token (torch.Tensor): token pair embedding
            [..., N_token, N_token, d]
        atom_to_token_idx (torch.Tensor): map atom idx to token idx
            [N_atom]

    Returns:
        z_gathered_blocked (torch.Tensor): atom pair embedding, with local blocked shape
            [..., n_trunks, n_queries, n_keys, d]
        pad_mask (torch.Tensor):
            [n_trunks, n_queries, n_keys]
        q_pad_length (int)
    """

    # [N_atom] -> [n_trunks, n_queries] and [n_trunks, n_keys]
    atom_to_token_idx_q, atom_to_token_idx_k, pad_info = rearrange_qk_to_dense_trunk(
        atom_to_token_idx,
        atom_to_token_idx,
        dim_q=-1,
        dim_k=-1,
        n_queries=n_queries,
        n_keys=n_keys,
        compute_mask=compute_mask,
    )

    z_gathered_blocked = gather_pair_embedding_in_dense_trunk(
        z_token, idx_q=atom_to_token_idx_q, idx_k=atom_to_token_idx_k
    )

    return z_gathered_blocked, pad_info


# from: https://github.com/huggingface/pytorch-image-models/blob/main/timm/layers/drop.py
def drop_path(
    x, drop_prob: float = 0.0, training: bool = False, scale_by_keep: bool = True
):
    """Drop paths (Stochastic Depth) per sample (when applied in main path of residual blocks).

    This is the same as the DropConnect impl I created for EfficientNet, etc networks, however,
    the original name is misleading as 'Drop Connect' is a different form of dropout in a separate paper...
    See discussion: https://github.com/tensorflow/tpu/issues/494#issuecomment-532968956 ... I've opted for
    changing the layer and argument names to 'drop path' rather than mix DropConnect as a layer name and use
    'survival rate' as the argument.

    """
    if drop_prob == 0.0 or not training:
        return x
    keep_prob = 1 - drop_prob
    shape = (x.shape[0],) + (1,) * (x.ndim - 1)
    random_tensor = x.new_empty(shape).bernoulli_(keep_prob)
    if keep_prob > 0.0 and scale_by_keep:
        random_tensor.div_(keep_prob)
    return x * random_tensor


class DropPath(nn.Module):
    """Drop paths (Stochastic Depth) per sample  (when applied in main path of residual blocks)."""

    def __init__(self, drop_prob: float = 0.0, scale_by_keep: bool = True):
        super(DropPath, self).__init__()
        self.drop_prob = drop_prob
        self.scale_by_keep = scale_by_keep

    def forward(self, x):
        return drop_path(x, self.drop_prob, self.training, self.scale_by_keep)

    def extra_repr(self):
        return f"drop_prob={round(self.drop_prob,3):0.3f}"
