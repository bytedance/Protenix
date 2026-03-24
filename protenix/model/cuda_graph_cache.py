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

"""
CUDA Graph caching for DiffusionTransformer with multi-shape support.

Inspired by vLLM's approach: pre-allocate max-size static buffers at fixed
GPU addresses, capture a separate CUDAGraph per (bucket, batch_shape) key.
Each graph records operations on the same static buffers — only the data
content changes between replays.

Enable with: PROTENIX_CUDA_GRAPHS=1
"""

import logging
import os
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

logger = logging.getLogger(__name__)

BUCKET_SIZES = (
    list(range(32, 256, 32))      # every 32: [32, 64, 96, ..., 224]
    + list(range(256, 2048 + 1, 64))  # every 64: [256, 320, 384, ..., 2048]
)
MAX_BUCKET = BUCKET_SIZES[-1]


_CUDA_GRAPHS_AUTO_MAX_DIFFUSION_TOKENS = 3000


def _cuda_graphs_mode() -> int:
    """Return CUDA graphs mode: 0=disabled, 1=auto, 2=always."""
    return int(os.environ.get("PROTENIX_CUDA_GRAPHS", "0"))


def is_cuda_graphs_enabled(
    n_token: Optional[int] = None,
    n_batch: int = 1,
) -> bool:
    """Check if CUDA graphs should be used for the given workload.

    Auto mode uses the total token-instances going through the diffusion
    transformer: n_batch * n_token, where n_batch includes both seed and
    sample batch dims. Breakeven on H200 is ~3000 total token-instances:
        - 1 seed × 5 samples × 500 tokens = 2500 (use graphs)
        - 5 seeds × 5 samples × 100 tokens = 2500 (use graphs)
        - 5 seeds × 5 samples × 200 tokens = 5000 (skip graphs)

    Args:
        n_token: Number of tokens per protein. None = unknown, assume small.
        n_batch: Total batch size (n_seeds * n_samples) from diffusion input shape.

    Modes (PROTENIX_CUDA_GRAPHS env var):
        0: disabled (default)
        1: auto — enabled when n_batch * n_token <= 3000 (recommended)
        2: always enabled regardless of workload size
    """
    mode = _cuda_graphs_mode()
    if mode == 0:
        return False
    if mode == 2:
        return True
    # mode == 1 (auto)
    if n_token is None:
        return True
    total = n_batch * n_token
    use = total <= _CUDA_GRAPHS_AUTO_MAX_DIFFUSION_TOKENS
    if not use:
        logger.info(
            f"CUDA graphs auto-disabled: n_token={n_token}, n_batch={n_batch}, "
            f"total={total} > {_CUDA_GRAPHS_AUTO_MAX_DIFFUSION_TOKENS}"
        )
    return use


def get_bucket_size(n: int) -> int:
    """Return the smallest bucket size >= n.

    Falls back to rounding up to the nearest multiple of 256 if n exceeds
    the largest pre-defined bucket.
    """
    for bucket in BUCKET_SIZES:
        if n <= bucket:
            return bucket
    return ((n + 255) // 256) * 256


def _pad(t: torch.Tensor, target: int, dim: int) -> torch.Tensor:
    """Zero-pad tensor ``t`` along ``dim`` so that dim size equals ``target``.

    Returns ``t`` unchanged if it already has the target size.
    """
    cur = t.shape[dim]
    if cur == target:
        return t
    ndim = t.ndim
    if dim < 0:
        dim = ndim + dim
    pad_spec = [0] * (2 * ndim)
    pad_spec[2 * (ndim - 1 - dim) + 1] = target - cur
    return F.pad(t, pad_spec)


def _make_pad_attn_mask(n_real: int, n_padded: int, device: torch.device,
                         dtype: torch.dtype = torch.float32) -> torch.Tensor:
    """Additive attention mask: 0 for valid pairs, -1e9 for padded.

    The large negative value ensures that after softmax, padded positions
    receive effectively zero weight (exp(-1e9) ≈ 0).

    Returns [1, 1, n_padded, n_padded]."""
    mask = torch.zeros(n_padded, device=device, dtype=dtype)
    mask[n_real:] = -1e9
    return (mask.unsqueeze(0) + mask.unsqueeze(1)).unsqueeze(0).unsqueeze(0)


def _replace_fused_layernorm(module: nn.Module) -> None:
    """Replace FusedLayerNorm with standard torch.nn.LayerNorm in-place.

    The custom fused CUDA kernel is not CUDA-graph-safe (produces different
    results on capture vs replay). Standard PyTorch LayerNorm works correctly.
    """
    from protenix.model.layer_norm.layer_norm import FusedLayerNorm

    for name, child in list(module.named_modules()):
        if isinstance(child, FusedLayerNorm):
            parts = name.split(".")
            parent = module
            for p in parts[:-1]:
                parent = getattr(parent, p)
            has_w = child.weight is not None
            has_b = child.bias is not None
            new_ln = torch.nn.LayerNorm(
                child.normalized_shape,
                elementwise_affine=has_w,
                bias=has_b if has_w else False,
            )
            if has_w:
                new_ln.weight.data.copy_(child.weight.data)
            if has_b:
                new_ln.bias.data.copy_(child.bias.data)
            new_ln = new_ln.to(
                device=child.weight.device if has_w else next(module.parameters()).device
            )
            setattr(parent, parts[-1], new_ln)


class _GraphEntry:
    """One captured CUDA graph for a specific (bucket, batch) shape."""
    __slots__ = ("graph", "static_a", "static_s", "static_z",
                 "static_mask", "static_output")

    def __init__(self):
        self.graph: Optional[torch.cuda.CUDAGraph] = None
        self.static_a: Optional[torch.Tensor] = None
        self.static_s: Optional[torch.Tensor] = None
        self.static_z: Optional[torch.Tensor] = None
        self.static_mask: Optional[torch.Tensor] = None
        self.static_output: Optional[torch.Tensor] = None


class DiffusionTransformerGraphCache:
    """Multi-shape CUDA graph cache for DiffusionTransformer.

    vLLM-style: one CUDAGraph per unique (bucket, batch_shape) key.
    Static buffers are allocated per key. The model is called directly
    (no make_graphed_callables), so it stays unmodified.
    """

    def __init__(self) -> None:
        self._entries: dict[tuple, _GraphEntry] = {}
        self._warmup_done: set[tuple] = set()
        self._ln_replaced = False

    def clear(self) -> None:
        self._entries.clear()
        self._warmup_done.clear()
        # Note: _ln_replaced is intentionally NOT reset — the LayerNorm
        # replacement is permanent on the module.

    def forward(
        self,
        transformer: nn.Module,
        a: torch.Tensor,
        s: torch.Tensor,
        z: torch.Tensor,
        **kwargs,
    ) -> torch.Tensor:
        # Replace FusedLayerNorm once (not CUDA-graph-safe)
        if not self._ln_replaced:
            _replace_fused_layernorm(transformer)
            self._ln_replaced = True
            logger.info("Replaced FusedLayerNorm with torch.nn.LayerNorm for CUDA graph safety")

        n_token = a.shape[-2]
        n_token_b = get_bucket_size(n_token)
        batch_dims = tuple(a.shape[:-2])
        key = (n_token_b, *batch_dims)

        # Pad inputs to bucket size
        a_pad = _pad(a, n_token_b, dim=-2)
        s_pad = _pad(s, n_token_b, dim=-2)
        z_pad = _pad(z, n_token_b, dim=-2)
        z_pad = _pad(z_pad, n_token_b, dim=-3)
        mask = _make_pad_attn_mask(n_token, n_token_b, a.device, a.dtype)

        # --- Replay ---
        if key in self._entries:
            entry = self._entries[key]
            entry.static_a.copy_(a_pad)
            entry.static_s.copy_(s_pad)
            entry.static_z.copy_(z_pad)
            entry.static_mask.copy_(mask)
            entry.graph.replay()
            return entry.static_output[..., :n_token, :].clone()

        # --- Warmup (first call for this key) ---
        if key not in self._warmup_done:
            logger.info(
                f"CUDA graph warmup: N_token={n_token}→{n_token_b}, "
                f"batch={batch_dims}"
            )
            self._warmup_done.add(key)
            # Set mask on module for blocks to read
            transformer._cuda_graph_pad_attn_mask = mask
            out = transformer(a_pad, s_pad, z_pad)
            transformer._cuda_graph_pad_attn_mask = None
            return out[..., :n_token, :].clone()

        # --- Capture (second call for this key) ---
        logger.info(
            f"CUDA graph capture: N_token={n_token}→{n_token_b}, "
            f"batch={batch_dims}, {len(transformer.blocks)} blocks"
        )

        entry = _GraphEntry()

        # Allocate static buffers at fixed addresses
        entry.static_a = a_pad.clone()
        entry.static_s = s_pad.clone()
        entry.static_z = z_pad.clone()
        entry.static_mask = mask.clone()

        # Set mask on module (will be read during graph capture)
        transformer._cuda_graph_pad_attn_mask = entry.static_mask

        # Warmup with static buffers (3 iterations like vLLM)
        for _ in range(3):
            entry.static_a.copy_(a_pad)
            entry.static_s.copy_(s_pad)
            entry.static_z.copy_(z_pad)
            entry.static_mask.copy_(mask)
            torch.cuda.synchronize()
            _ = transformer(entry.static_a, entry.static_s, entry.static_z)
            torch.cuda.synchronize()

        # Re-copy (warmup may have modified in-place)
        entry.static_a.copy_(a_pad)
        entry.static_s.copy_(s_pad)
        entry.static_z.copy_(z_pad)
        entry.static_mask.copy_(mask)
        torch.cuda.synchronize()

        # Capture
        entry.graph = torch.cuda.CUDAGraph()
        with torch.cuda.graph(entry.graph):
            entry.static_output = transformer(
                entry.static_a, entry.static_s, entry.static_z
            )

        self._entries[key] = entry

        # Clean up module attribute (entries have their own static mask)
        # Don't set to None — leave the static_mask reference so the
        # graph replay reads from it.

        return entry.static_output[..., :n_token, :].clone()


_cache: Optional[DiffusionTransformerGraphCache] = None


def get_diff_transformer_cache() -> DiffusionTransformerGraphCache:
    global _cache
    if _cache is None:
        _cache = DiffusionTransformerGraphCache()
    return _cache
