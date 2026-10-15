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

from __future__ import annotations

from contextlib import nullcontext
from typing import Any, Optional

import torch


ACCELERATOR_DEVICE_TYPES = ("cuda", "xpu")


def is_xpu_available() -> bool:
    return hasattr(torch, "xpu") and torch.xpu.is_available()


def get_available_device_type() -> str:
    if torch.cuda.is_available():
        return "cuda"
    if is_xpu_available():
        return "xpu"
    return "cpu"


def get_default_device(local_rank: int = 0) -> torch.device:
    device_type = get_available_device_type()
    if device_type == "cpu":
        return torch.device("cpu")
    return torch.device(f"{device_type}:{local_rank}")


def get_device_count(device_type: Optional[str] = None) -> int:
    device_type = device_type or get_available_device_type()
    if device_type == "cuda":
        return torch.cuda.device_count()
    if device_type == "xpu":
        return torch.xpu.device_count() if hasattr(torch, "xpu") else 0
    if device_type == "cpu":
        return 0
    raise ValueError(f"Unsupported device type: {device_type}")


def set_device(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.set_device(device)
    elif device.type == "xpu":
        torch.xpu.set_device(device)


def _extract_device_type(obj: Any) -> Optional[str]:
    if isinstance(obj, torch.Tensor):
        return obj.device.type
    if isinstance(obj, torch.device):
        return obj.type
    if isinstance(obj, str) and obj.split(":", 1)[0] in {
        "cpu",
        "cuda",
        "xpu",
    }:
        return torch.device(obj).type
    if isinstance(obj, dict):
        found_cpu = False
        for value in obj.values():
            device_type = _extract_device_type(value)
            if device_type in ACCELERATOR_DEVICE_TYPES:
                return device_type
            found_cpu |= device_type == "cpu"
        return "cpu" if found_cpu else None
    if isinstance(obj, (list, tuple, set)):
        found_cpu = False
        for value in obj:
            device_type = _extract_device_type(value)
            if device_type in ACCELERATOR_DEVICE_TYPES:
                return device_type
            found_cpu |= device_type == "cpu"
        return "cpu" if found_cpu else None
    return None


def infer_device_type(*objects: Any, default: Optional[str] = None) -> str:
    found_cpu = False
    for obj in objects:
        device_type = _extract_device_type(obj)
        if device_type in ACCELERATOR_DEVICE_TYPES:
            return device_type
        found_cpu |= device_type == "cpu"
    return "cpu" if found_cpu else default or get_available_device_type()


def autocast_context(
    *,
    device_type: Optional[str] = None,
    dtype: Optional[torch.dtype] = None,
    enabled: bool = True,
    cache_enabled: Optional[bool] = None,
):
    device_type = device_type or get_available_device_type()
    if device_type not in ACCELERATOR_DEVICE_TYPES:
        return nullcontext()
    kwargs: dict[str, Any] = {"device_type": device_type, "enabled": enabled}
    if dtype is not None:
        kwargs["dtype"] = dtype
    if cache_enabled is not None:
        kwargs["cache_enabled"] = cache_enabled
    return torch.autocast(**kwargs)


def disable_autocast(device: torch.device):
    if not isinstance(device, torch.device):
        raise TypeError("device must be a torch.device")
    if device.type not in ACCELERATOR_DEVICE_TYPES:
        return nullcontext()
    return torch.amp.autocast(device.type, enabled=False)


def empty_cache(device: torch.device) -> None:
    if not isinstance(device, torch.device):
        raise TypeError("device must be a torch.device")
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()


def synchronize(device: Any = None) -> None:
    device_type = infer_device_type(device, default=get_available_device_type())
    if device_type == "cuda":
        torch.cuda.synchronize()
    elif device_type == "xpu":
        torch.xpu.synchronize()


def manual_seed_all(seed: int) -> None:
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if is_xpu_available():
        torch.xpu.manual_seed_all(seed)


def is_autocast_enabled(device_type: str) -> bool:
    try:
        return torch.is_autocast_enabled(device_type)
    except TypeError:
        if device_type == "cuda":
            return torch.is_autocast_enabled()
        return False


def get_autocast_dtype(device_type: str) -> torch.dtype:
    if hasattr(torch, "get_autocast_dtype"):
        try:
            return torch.get_autocast_dtype(device_type)
        except TypeError:
            pass
    if device_type == "cuda" and hasattr(torch, "get_autocast_gpu_dtype"):
        return torch.get_autocast_gpu_dtype()
    raise RuntimeError(f"Autocast dtype is unavailable for device type {device_type}")
