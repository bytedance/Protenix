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
# Copyright 2021 AlQuraishi Laboratory


from abc import ABC, abstractmethod
from functools import partial, partialmethod
from typing import List, Optional

import torch
import torch.nn as nn

from protenix.model.triangular.layers import Attention, LayerNorm, OpenfoldLinear
from protenix.model.utils import chunk_layer, is_fp16_enabled, permute_final_dims


def kernel_triangular_mult(
    x: torch.Tensor,
    direction: str,
    mask: torch.Tensor,
    norm_in_weight: torch.Tensor,
    norm_in_bias: torch.Tensor,
    p_in_weight: torch.Tensor,
    g_in_weight: torch.Tensor,
    norm_out_weight: torch.Tensor,
    norm_out_bias: torch.Tensor,
    p_out_weight: torch.Tensor,
    g_out_weight: torch.Tensor,
    eps: float,
) -> torch.Tensor:
    from cuequivariance_torch.primitives.triangle import triangle_multiplicative_update

    """
        This function performs a triangle multiplicative update operation, which is a key component
    in the AlphaFold2 architecture. The operation consists of:

    1. Input normalization and gating
    2. Triangular projection (either outgoing or incoming)
    3. Output normalization and gating

    The function supports both ahead-of-time (AOT) tuning and just-in-time (JIT) tuning.
    Auto-tuning behavior can be controlled through environment variables:

    - Quick testing: Default configuration where tuning configs, if existent, are looked-up.
        If not, then falls back to default kernel parameters. No tuning is performed.
    - On-Demand tuning: Set `CUEQ_TRITON_TUNING= "ONDEMAND"` to auto-tune for new shapes encountered on first run (may take several minutes)
    - AOT tuning: Set `CUEQ_TRITON_TUNING= "AOT"` to perform full ahead-of-time tuning for optimal performance **(may take several hours)**
    - Ignore user cache: Set CUEQ_TRITON_IGNORE_EXISTING_CACHE to ignore both the default settings that come with the package
        and any user-local settings previously saved with AOT/ONDEMAND tuning. May be used to regenerate optimal settings for a particular setup.
    - Cache directory: Set `CUEQ_TRITON_CACHE_DIR` to specify where tuning configurations are stored
    - Note: When using Docker with default or on-demand tuning enabled, commit the container to persist tuning changes

    Notes:
        (1) Context is saved for backward pass. You don't need to save it manually.
        (2) Kernel precision (fp32, bf16, fp16) is based on input dtypes. For tf32, set it from torch global
            scope using torch.backends.cuda.matmul.allow_tf32
        (3) **Limitation**: Currently only supports hidden_dim values that are multiples of 32.
        (4) We have moved away from the default round-towards-zero (RZ) implementation to round-nearest (RN)
            for better tf32 accuracy in cuex.triangle_multiplicative_update. In rare circumstances,
            this may cause minor differences in results observed.
        (5) When using torch compile, use `cueuivariance_ops_torch.init_triton_cache()` to initialize
            triton cache before calling torch compiled triangular multiplicative update.
        (6) Although the example demonstrates the most common case of one batch dimension,
            the API supports variable number of leading batch dimensions.
    """
    return triangle_multiplicative_update(
        x,
        direction=direction,
        mask=mask,
        norm_in_weight=norm_in_weight,
        norm_in_bias=norm_in_bias,
        p_in_weight=p_in_weight,
        g_in_weight=g_in_weight,
        norm_out_weight=norm_out_weight,
        norm_out_bias=norm_out_bias,
        p_out_weight=p_out_weight,
        g_out_weight=g_out_weight,
        eps=eps,
    )


class BaseTriangleMultiplicativeUpdate(nn.Module, ABC):
    """
    Implements Algorithms 11 and 12.

    Args:
        c_z:
            Input channel dimension
        c_hidden:
            Hidden channel dimension
        _outgoing:
            Whether this is an outgoing update
    """

    @abstractmethod
    def __init__(self, c_z: int, c_hidden: int, _outgoing: bool) -> None:
        super(BaseTriangleMultiplicativeUpdate, self).__init__()
        self.c_z = c_z
        self.c_hidden = c_hidden
        self._outgoing = _outgoing

        self.linear_g = OpenfoldLinear(self.c_z, self.c_z, bias=False, init="gating")
        self.linear_z = OpenfoldLinear(
            self.c_hidden, self.c_z, bias=False, init="final"
        )

        self.layer_norm_in = LayerNorm(self.c_z)
        self.layer_norm_out = LayerNorm(self.c_hidden)

        self.sigmoid = nn.Sigmoid()

    def _combine_projections(
        self,
        a: torch.Tensor,
        b: torch.Tensor,
        _inplace_chunk_size: Optional[int] = None,
    ) -> torch.Tensor:
        if self._outgoing:
            a = permute_final_dims(a, (2, 0, 1))
            b = permute_final_dims(b, (2, 1, 0))
        else:
            a = permute_final_dims(a, (2, 1, 0))
            b = permute_final_dims(b, (2, 0, 1))

        if _inplace_chunk_size is not None:
            # To be replaced by torch vmap
            for i in range(0, a.shape[-3], _inplace_chunk_size):
                a_chunk = a[..., i : i + _inplace_chunk_size, :, :]
                b_chunk = b[..., i : i + _inplace_chunk_size, :, :]
                a[..., i : i + _inplace_chunk_size, :, :] = torch.matmul(
                    a_chunk,
                    b_chunk,
                )

            p = a
        else:
            p = torch.matmul(a, b)

        return permute_final_dims(p, (1, 2, 0))

    @abstractmethod
    def forward(
        self,
        z: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        inplace_safe: bool = False,
        _add_with_inplace: bool = False,
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, N_res, N_res, C_z] input tensor
            mask:
                [*, N_res, N_res] input mask
        Returns:
            [*, N_res, N_res, C_z] output tensor
        """
        pass


class TriangleMultiplicativeUpdate(BaseTriangleMultiplicativeUpdate):
    """
    Implements Algorithms 11 and 12.

    Args:
        c_z:
            Input channel dimension
        c_hidden:
            Hidden channel dimension
        _outgoing:
            Whether this is an outgoing update. Defaults to True.
    """

    def __init__(self, c_z: int, c_hidden: int, _outgoing: bool = True) -> None:
        super(TriangleMultiplicativeUpdate, self).__init__(
            c_z=c_z, c_hidden=c_hidden, _outgoing=_outgoing
        )

        self.linear_a_p = OpenfoldLinear(self.c_z, self.c_hidden, bias=False)
        self.linear_a_g = OpenfoldLinear(
            self.c_z, self.c_hidden, bias=False, init="gating"
        )
        self.linear_b_p = OpenfoldLinear(self.c_z, self.c_hidden, bias=False)
        self.linear_b_g = OpenfoldLinear(
            self.c_z, self.c_hidden, bias=False, init="gating"
        )

    def _inference_forward(
        self,
        z: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        inplace_chunk_size: Optional[int] = None,
        with_add: bool = True,
    ) -> torch.Tensor:
        """
        Args:
            z:
                A [*, N, N, C_z] pair representation
            mask:
                A [*, N, N] pair mask
            inplace_chunk_size:
                Size of chunks used in the main computation. Increase to trade
                memory for speed.
            with_add:
                If True, z is overwritten with (z + update). Otherwise, it is
                overwritten with (update).
        Returns:
            A reference to the overwritten z

        More memory-efficient, inference-only version of the forward function.
        Uses in-place operations, fusion of the addition that happens after
        this module in the Evoformer, a smidge of recomputation, and
        a cache of overwritten values to lower peak memory consumption of this
        module from 5x the size of the input tensor z to 2.5x its size. Useful
        for inference on extremely long sequences.

        It works as follows. We will make reference to variables used in the
        default forward implementation below. Naively, triangle multiplication
        attention requires the manifestation of 5 tensors the size of z:
        1) z, the "square" input tensor, 2) a, the first projection of z,
        3) b, the second projection of b, 4) g, a z-sized mask, and 5) a
        z-sized tensor for intermediate computations. For large N, this is
        prohibitively expensive; for N=4000, for example, z is more than 8GB
        alone. To avoid this problem, we compute b, g, and all intermediate
        tensors in small chunks, noting that the chunks required to compute a
        chunk of the output depend only on the tensor a and corresponding
        vertical and horizontal chunks of z. This suggests an algorithm that
        loops over pairs of chunks of z: hereafter "columns" and "rows" of
        z, even though each "column" and "row" in fact contains
        inplace_chunk_size contiguous true columns and rows of z. Writing
        output chunks to a new tensor would bring total memory consumption
        down to 3x the size of z. However, more memory can be saved by writing
        output chunks directly to z in-place. WLOG, we choose to write output
        chunks vertically, overwriting the ith "column" of z at the end of
        the ith iteration of the main loop. Despite this overwriting, the
        ith column is always one column ahead of previously overwritten columns
        and can be recovered directly from z. After the first iteration,
        however, the ith row of z is always at least partially overwritten. For
        this reason, we introduce the z-cache, a tensor one-half the size of
        z. The z-cache initially contains the left half (2nd and 3rd quadrants)
        of z. For 0 < i < N/2, the missing left part of the ith row of z is
        recovered from this cache at the beginning of the ith iteration. Once i
        exceeds n/2, the cache is "reoriented" to encompass the 3rd and 4th
        quadrants of z instead. Though the 3rd quadrant of the original z is
        entirely overwritten at this point, it can be recovered from the z-cache
        itself. Thereafter, the ith row of z can be recovered in its entirety
        from the reoriented z-cache. After the final iteration, z has been
        completely overwritten and contains the triangular multiplicative
        update. If with_add is True, it instead contains the sum of z and the
        triangular multiplicative update. In either case, peak memory
        consumption is just 2.5x the size of z, disregarding memory used for
        chunks and other small variables.
        """
        if mask is None:
            mask = z.new_ones(z.shape[:-1])

        mask = mask.unsqueeze(-1)

        def compute_projection_helper(pair, mask, a=True):
            if a:
                linear_g = self.linear_a_g
                linear_p = self.linear_a_p
            else:
                linear_g = self.linear_b_g
                linear_p = self.linear_b_p

            pair = self.layer_norm_in(pair)
            p = linear_g(pair)
            p.sigmoid_()
            p *= linear_p(pair)
            p *= mask
            p = permute_final_dims(p, (2, 0, 1))
            return p

        def compute_projection(pair, mask, a=True, chunked=True):
            need_transpose = self._outgoing ^ a
            if not chunked:
                p = compute_projection_helper(pair, mask, a)
                if need_transpose:
                    p = p.transpose(-1, -2)
            else:
                # This computation is chunked so as not to exceed our 2.5x
                # budget with a large intermediate tensor
                linear_g = self.linear_a_g if a else self.linear_b_g
                c = linear_g.weight.shape[0]
                out_shape = pair.shape[:-3] + (c,) + pair.shape[-3:-1]
                p = pair.new_zeros(out_shape)
                for i in range(0, pair.shape[-3], inplace_chunk_size):
                    pair_chunk = pair[..., i : i + inplace_chunk_size, :, :]
                    mask_chunk = mask[..., i : i + inplace_chunk_size, :, :]
                    pair_chunk = compute_projection_helper(
                        pair[..., i : i + inplace_chunk_size, :, :],
                        mask[..., i : i + inplace_chunk_size, :, :],
                        a,
                    )
                    if need_transpose:
                        pair_chunk = pair_chunk.transpose(-1, -2)
                        p[..., i : i + inplace_chunk_size] = pair_chunk
                    else:
                        p[..., i : i + inplace_chunk_size, :] = pair_chunk

                    del pair_chunk

            return p

        # We start by fully manifesting a. In addition to the input, this
        # brings total memory consumption to 2x z (disregarding size of chunks)
        # [*, N, N, c]
        a = compute_projection(z, mask, True, chunked=True)

        if inplace_chunk_size is not None:
            n = a.shape[-1]
            half_n = n // 2 + n % 2
            row_dim = -3
            col_dim = -2
            b_chunk_dim = row_dim if self._outgoing else col_dim

            def empty_slicer(t):
                return [slice(None) for _ in t.shape]

            def slice_tensor(t, start, end, dim):
                # Slices start:end from the dim dimension of t
                s = empty_slicer(t)
                s[dim] = slice(start, end)
                return t[tuple(s)]

            def flip_z_cache_(z_cache, z):
                # "Reorient" the z_cache (see below), filling it with quadrants
                # 3---recovered from the z_cache---and 4---recovered from z---
                # of the input tensor z.
                quadrant_3 = slice_tensor(z_cache, half_n, None, row_dim)
                z_cache = z_cache.transpose(row_dim, col_dim)

                # If n is odd, we need to shrink the z_cache by one row
                z_cache = z_cache[..., : (n // 2), :, :]

                # Move the 3rd quadrant of z into the
                first_half_slicer = empty_slicer(z_cache)
                first_half_slicer[col_dim] = slice(0, half_n)
                z_cache[tuple(first_half_slicer)] = quadrant_3

                # Get the fourth quadrant of z
                quadrant_4 = slice_tensor(z, half_n, None, row_dim)
                quadrant_4 = slice_tensor(quadrant_4, half_n, None, col_dim)

                # Insert said quadrant into the rotated z-cache
                quadrant_3_slicer = empty_slicer(z_cache)
                quadrant_3_slicer[col_dim] = slice(half_n, None)

                z_cache[tuple(quadrant_3_slicer)] = quadrant_4

                return z_cache

            # Initialize the z cache to the left half of z.
            z_cache_shape = list(z.shape)
            z_cache_shape[col_dim] = half_n
            z_cache = z.new_zeros(z_cache_shape)
            z_cache_slicer = empty_slicer(z_cache)
            z_cache_slicer[col_dim] = slice(0, half_n)
            z_cache.copy_(z[tuple(z_cache_slicer)])
            z_cache_rotated = False

            # We need to reorient the z-cache at the halfway point, and we
            # don't want a single chunk to straddle that point. We contract one
            # of the chunks in the middle to address that problem.
            i_range = list(range(0, half_n, inplace_chunk_size))
            initial_offsets = [
                i_2 - i_1 for i_1, i_2 in zip(i_range, i_range[1:] + [half_n])
            ]
            after_half = list(range(half_n, n, inplace_chunk_size))
            after_half_offsets = [inplace_chunk_size for _ in after_half]
            combined_range_with_offsets = zip(
                i_range + after_half, initial_offsets + after_half_offsets
            )
            for i, offset in combined_range_with_offsets:
                if not z_cache_rotated and i >= half_n:
                    z_cache = flip_z_cache_(z_cache, z)
                    z_cache_rotated = True

                z_chunk_b = slice_tensor(
                    z,
                    i,
                    i + offset,
                    b_chunk_dim,
                )
                mask_chunk = slice_tensor(
                    mask,
                    i,
                    i + offset,
                    b_chunk_dim,
                )

                z_chunk_b = z_chunk_b.clone()
                if b_chunk_dim == col_dim:
                    z_chunk_b = slice_tensor(z, i, i + offset, col_dim)
                else:  # b_chunk_dim == row_dim
                    # In this case, the b-dimension (b_chunk_dim) is partially
                    # overwritten at the end of each iteration. We need to
                    # restore the missing component from the z-cache.
                    if not z_cache_rotated:
                        z_chunk_slicer = empty_slicer(z_chunk_b)
                        z_chunk_slicer[col_dim] = slice(0, half_n)
                        z_chunk_b[tuple(z_chunk_slicer)] = slice_tensor(
                            z_cache,
                            i,
                            i + offset,
                            row_dim,
                        )
                    else:
                        z_cache_offset = i - half_n
                        z_chunk_b = slice_tensor(
                            z_cache, z_cache_offset, z_cache_offset + offset, row_dim
                        )

                b_chunk = compute_projection(
                    z_chunk_b, mask_chunk, a=False, chunked=False
                )
                del z_chunk_b

                x_chunk = torch.matmul(
                    a,
                    b_chunk,
                )
                x_chunk = permute_final_dims(x_chunk, (1, 2, 0))
                x_chunk = self.layer_norm_out(x_chunk)
                x_chunk = self.linear_z(x_chunk)

                # The g dimension (col_dim) is parallel to and ahead of the
                # overwrites in z. We can extract the g chunk normally.
                z_chunk_g = slice_tensor(z, i, i + offset, col_dim)
                g_chunk = self.linear_g(self.layer_norm_in(z_chunk_g))
                g_chunk.sigmoid_()
                del z_chunk_g

                x_chunk *= g_chunk

                # Write the columns into z in-place
                z_slicer = empty_slicer(z)
                z_slicer[col_dim] = slice(i, i + offset)
                if with_add:
                    z[tuple(z_slicer)] += x_chunk
                else:
                    z[tuple(z_slicer)] = x_chunk
        else:
            b = compute_projection(z, mask, False, False)
            x = torch.matmul(a, b)
            x = self.layer_norm_out(x)
            x = self.linear_z(x)
            g = self.linear_g(z)
            g.sigmoid_()
            x *= g
            if with_add:
                z += x
            else:
                z = x

        return z

    def forward(
        self,
        z: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        inplace_safe: bool = False,
        _add_with_inplace: bool = False,
        _inplace_chunk_size: Optional[int] = 256,
        triangle_multiplicative: str = "torch",
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, N_res, N_res, C_z] input tensor
            mask:
                [*, N_res, N_res] input mask
        Returns:
            [*, N_res, N_res, C_z] output tensor
        """
        _input_inplace_safe = inplace_safe is True

        # Note: cuequivariance requires that the hidden dimension c must equal c_z.
        # If this condition is not met, an AssertionError will be raised.
        # Therefore, we include a check here: if c != c_z, we fall back to using plain PyTorch.
        # This situation may occur in our template module.
        if triangle_multiplicative == "cuequivariance" and (self.c_z == self.c_hidden):
            if _input_inplace_safe and _add_with_inplace:
                z_in = z.clone()
            z = kernel_triangular_mult(
                z[None],
                direction="outgoing" if self._outgoing else "incoming",
                mask=z.new_ones(z.shape[:-1])[None] if mask is None else mask,
                norm_in_weight=self.layer_norm_in.weight,
                norm_in_bias=self.layer_norm_in.bias,
                p_in_weight=torch.cat(
                    [self.linear_a_p.weight, self.linear_b_p.weight], 0
                ),
                g_in_weight=torch.cat(
                    [self.linear_a_g.weight, self.linear_b_g.weight], 0
                ),
                norm_out_weight=self.layer_norm_out.weight,
                norm_out_bias=self.layer_norm_out.bias,
                p_out_weight=self.linear_z.weight,
                g_out_weight=self.linear_g.weight,
                eps=1e-5,  # In BF16, we use the default eps of 1e-5.
            )[0]
            if _input_inplace_safe and _add_with_inplace:
                return z + z_in
            else:
                return z
        elif (triangle_multiplicative == "torch") or (self.c_z != self.c_hidden):
            if inplace_safe:
                x = self._inference_forward(
                    z,
                    mask,
                    inplace_chunk_size=_inplace_chunk_size,
                    with_add=_add_with_inplace,
                )
                return x

            if mask is None:
                mask = z.new_ones(z.shape[:-1])

            mask = mask.unsqueeze(-1)

            if _input_inplace_safe and _add_with_inplace:
                z_in = z.clone()

            z = self.layer_norm_in(z)
            a = mask
            a = a * self.sigmoid(self.linear_a_g(z))
            a = a * self.linear_a_p(z)
            b = mask
            b = b * self.sigmoid(self.linear_b_g(z))
            b = b * self.linear_b_p(z)

            # Prevents overflow of torch.matmul in combine projections in
            # reduced-precision modes
            a_std = a.std()
            b_std = b.std()
            if is_fp16_enabled() and a_std != 0.0 and b_std != 0.0:
                a = a / a.std()
                b = b / b.std()

            if is_fp16_enabled():
                with torch.amp.autocast("cuda", enabled=False):
                    x = self._combine_projections(a.float(), b.float())
            else:
                x = self._combine_projections(a, b)

            del a, b
            x = self.layer_norm_out(x)
            x = self.linear_z(x)
            g = self.sigmoid(self.linear_g(z))
            x = x * g
            if _input_inplace_safe and _add_with_inplace:
                x = x + z_in
            return x
        else:
            raise ValueError(
                f"triangle_multiplicative must be 'cuequivariance' or 'torch', but got {triangle_multiplicative}"
            )


class TriangleMultiplicationOutgoing(TriangleMultiplicativeUpdate):
    """
    Implements Algorithm 11.
    """

    __init__ = partialmethod(TriangleMultiplicativeUpdate.__init__, _outgoing=True)


class TriangleMultiplicationIncoming(TriangleMultiplicativeUpdate):
    """
    Implements Algorithm 12.
    """

    __init__ = partialmethod(TriangleMultiplicativeUpdate.__init__, _outgoing=False)


class TriangleAttention(nn.Module):
    """
    Triangle attention.

    Args:
        c_in:
            Input channel dimension
        c_hidden:
            Overall hidden channel dimension (not per-head)
        no_heads:
            Number of attention heads
        starting:
            Whether this is a starting node attention. Defaults to True.
        inf:
            Value for attention masking. Defaults to 1e9.
    """

    def __init__(
        self,
        c_in: int,
        c_hidden: int,
        no_heads: int,
        starting: bool = True,
        inf: float = 1e9,
    ) -> None:
        super(TriangleAttention, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.starting = starting
        self.inf = inf

        self.layer_norm = LayerNorm(self.c_in)

        self.linear = OpenfoldLinear(c_in, self.no_heads, bias=False)

        self.mha = Attention(
            self.c_in, self.c_in, self.c_in, self.c_hidden, self.no_heads
        )

    @torch.jit.ignore
    def _chunk(
        self,
        x: torch.Tensor,
        biases: List[torch.Tensor],
        chunk_size: int,
        triangle_attention: str = "torch",
        inplace_safe: bool = False,
    ) -> torch.Tensor:
        "triangle! triangle!"
        mha_inputs = {
            "q_x": x,
            "kv_x": x,
            "biases": biases,
        }

        return chunk_layer(
            partial(
                self.mha,
                triangle_attention=triangle_attention,
            ),
            mha_inputs,
            chunk_size=chunk_size,
            no_batch_dims=len(x.shape[:-2]),
            _out=x if inplace_safe else None,
        )

    def forward(
        self,
        x: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        chunk_size: Optional[int] = None,
        triangle_attention: str = "torch",
        inplace_safe: bool = False,
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, I, J, C_in] input tensor (e.g. the pair representation)
        Returns:
            [*, I, J, C_in] output tensor
        """
        if mask is None:
            # [*, I, J]
            mask = x.new_ones(
                x.shape[:-1],
            )

        if not self.starting:
            x = x.transpose(-2, -3)
            mask = mask.transpose(-1, -2)

        # [*, I, J, C_in]
        x = self.layer_norm(x)

        # [*, I, 1, 1, J]
        mask_bias = (self.inf * (mask - 1))[..., :, None, None, :]

        # [*, H, I, J]
        triangle_bias = permute_final_dims(self.linear(x), (2, 0, 1))

        # [*, 1, H, I, J]
        triangle_bias = triangle_bias.unsqueeze(-4)

        biases = [mask_bias, triangle_bias]

        if chunk_size is not None:
            x = self._chunk(
                x,
                biases,
                chunk_size,
                triangle_attention=triangle_attention,
                inplace_safe=inplace_safe,
            )
        else:
            x = self.mha(
                q_x=x,
                kv_x=x,
                biases=biases,
                triangle_attention=triangle_attention,
            )

        if not self.starting:
            x = x.transpose(-2, -3)

        return x


# Implements Algorithm 13
TriangleAttentionStartingNode = TriangleAttention


class TriangleAttentionEndingNode(TriangleAttention):
    """
    Implements Algorithm 14.
    """

    __init__ = partialmethod(TriangleAttention.__init__, starting=False)
