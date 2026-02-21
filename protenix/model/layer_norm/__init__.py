# Copyright 2024 ByteDance and/or its affiliates.
#
# Copyright 2021- HPC-AI Technology Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http:#www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
LayerNorm implementations.

Note: importing the CUDA fused implementation eagerly can trigger JIT compilation
of a C++/CUDA extension. To avoid import-time side effects (e.g. during test
discovery), we expose `FusedLayerNorm` via a lazy attribute.
"""

from __future__ import annotations

from typing import Any

__all__ = ["FusedLayerNorm"]


def __getattr__(name: str) -> Any:
    if name == "FusedLayerNorm":
        from .layer_norm import FusedLayerNorm  # local import to keep it lazy

        return FusedLayerNorm
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
