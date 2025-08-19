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

import os
from pathlib import Path

import torch
import triton

from protenix.utils.distributed import DIST_WRAPPER

device_capability = torch.cuda.get_device_capability()
device_capability = f"{device_capability[0]}-{device_capability[1]}"

device_name = torch.cuda.get_device_name().replace(" ", "-")


def get_config_dir() -> Path:
    current_dir = Path(__file__).parent
    tri_cache_root = current_dir / "TriAttentionCache"
    rank = f"{DIST_WRAPPER.rank}"
    target_dir = tri_cache_root / rank
    if not target_dir.exists():
        target_dir.mkdir(parents=True, exist_ok=False)
        precache_dir = tri_cache_root / "0"
        if precache_dir.exists() and any(precache_dir.iterdir()):
            ret = os.system(f"cp {precache_dir}/* {target_dir}")
    return target_dir


config_dir = get_config_dir()


def config_to_dict(config: triton.Config) -> dict:
    # This assume we are not making use of `pre_hook` in the `triton.Config`
    return {
        "kwargs": config.kwargs,
        "num_warps": config.num_warps,
        "num_stages": config.num_stages,
    }


def dict_to_config(d: dict) -> triton.Config:
    return triton.Config(
        kwargs=d["kwargs"],
        num_warps=d["num_warps"],
        num_stages=d["num_stages"],
    )


# Base configs that should be ~ok for things <= 512.
_attention_fwd_configs = [
    triton.Config(kwargs={"BLOCK_M": 16, "BLOCK_N": 32}, num_warps=1, num_stages=2),
    triton.Config(kwargs={"BLOCK_M": 32, "BLOCK_N": 16}, num_warps=1, num_stages=2),
    triton.Config(kwargs={"BLOCK_M": 64, "BLOCK_N": 32}, num_warps=4, num_stages=2),
    triton.Config(kwargs={"BLOCK_M": 64, "BLOCK_N": 16}, num_warps=2, num_stages=3),
    triton.Config(kwargs={"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=1, num_stages=1),
    triton.Config(kwargs={"BLOCK_M": 128, "BLOCK_N": 16}, num_warps=2, num_stages=2),
    triton.Config(kwargs={"BLOCK_M": 128, "BLOCK_N": 32}, num_warps=2, num_stages=2),
    triton.Config(kwargs={"BLOCK_M": 32, "BLOCK_N": 64}, num_warps=2, num_stages=2),
    triton.Config(kwargs={"BLOCK_M": 64, "BLOCK_N": 32}, num_warps=2, num_stages=3),
    triton.Config(kwargs={"BLOCK_M": 32, "BLOCK_N": 16}, num_warps=1, num_stages=4),
]

_attention_bwd_dq_configs = [
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 32}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 16}, num_warps=1, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 64}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=1, num_stages=1),
    triton.Config({"BLOCK_M": 64, "BLOCK_N": 16}, num_warps=2, num_stages=3),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 16}, num_warps=1, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 16}, num_warps=1, num_stages=3),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 64}, num_warps=1, num_stages=1),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 64}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 64, "BLOCK_N": 32}, num_warps=2, num_stages=2),
]

_attention_bwd_dkdv_configs = [
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=1, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 16}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 64, "BLOCK_N": 16}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=1, num_stages=3),
    triton.Config({"BLOCK_M": 64, "BLOCK_N": 32}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 16}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 16}, num_warps=1, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=1, num_stages=1),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 32}, num_warps=2, num_stages=1),
]

_attention_bwd_dbias2_configs = [
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 16}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 16}, num_warps=2, num_stages=3),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 32}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=1, num_stages=2),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 64}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 64}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 16}, num_warps=2, num_stages=1),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 16, "BLOCK_N": 16}, num_warps=4, num_stages=3),
    triton.Config({"BLOCK_M": 32, "BLOCK_N": 32}, num_warps=2, num_stages=3),
]

_attention_bwd_preprocess_configs = [
    triton.Config({"BLOCK_M": 16}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 16}, num_warps=2, num_stages=3),
    triton.Config({"BLOCK_M": 16}, num_warps=4, num_stages=4),
    triton.Config({"BLOCK_M": 16}, num_warps=4, num_stages=5),
    triton.Config({"BLOCK_M": 32}, num_warps=2, num_stages=2),
    triton.Config({"BLOCK_M": 32}, num_warps=2, num_stages=3),
    triton.Config({"BLOCK_M": 32}, num_warps=4, num_stages=4),
    triton.Config({"BLOCK_M": 32}, num_warps=4, num_stages=5),
]
