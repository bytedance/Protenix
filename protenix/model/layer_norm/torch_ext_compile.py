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

import json
import logging
import os
import re
import subprocess
import warnings
from pathlib import Path
from typing import Any, List, Optional, Tuple

from torch.utils.cpp_extension import load

_LOGGER = logging.getLogger(__name__)

_CUDA_ROOT = Path("/usr/local/cuda")
_NVCCVersion = Tuple[int, int]
_NVCC_RELEASE_RE = re.compile(r"release\s+(\d+\.\d+)", re.IGNORECASE)
_FALLBACK_ARCH_LIST = "7.5;8.0;8.6;9.0"

# nvcc supported version range for each architecture (arch, min_ver, max_ver)
# max_ver=None indicates that there is currently no known upper limit
_ARCH_SUPPORT = [
    ("7.0", (10, 0), (13, 0)),  # Volta  V100
    ("7.5", (10, 0), None),  # Turing T4/2080
    ("8.0", (11, 1), None),  # Ampere A100
    ("8.6", (11, 1), None),  # Ampere RTX 3090
    ("8.7", (11, 4), None),  # Ampere Jetson Orin
    ("8.9", (11, 8), None),  # Ada RTX 4090
    ("9.0", (11, 8), None),  # Hopper H100
    ("10.0", (12, 8), None),  # Blackwell B100
    ("12.0", (12, 8), None),  # Blackwell Next
]


def _format_nvcc_version(version: _NVCCVersion) -> str:
    return f"{version[0]}.{version[1]}"


def _parse_major_minor(value: str) -> Optional[_NVCCVersion]:
    match = re.search(r"(\d+)\.(\d+)", value)
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def _parse_cuda_version_file(path: Path) -> Optional[_NVCCVersion]:
    try:
        if path.suffix == ".json":
            version_str = (
                json.loads(path.read_text()).get("cuda", {}).get("version", "")
            )
        else:
            version_str = path.read_text()
    except Exception:
        return None

    return _parse_major_minor(version_str)


def _probe_cuda_version_from_files() -> Optional[_NVCCVersion]:
    for version_file in (_CUDA_ROOT / "version.json", _CUDA_ROOT / "version.txt"):
        parsed = _parse_cuda_version_file(version_file)
        if parsed is not None:
            return parsed
    return None


def _get_nvcc_version() -> Optional[_NVCCVersion]:
    """
    Detection sequence:
    1. nvcc --version — Most direct, the actual compiler
    2. /usr/local/cuda/version.{json,txt} — Approximate value when nvcc is not in PATH
    """
    try:
        out = subprocess.check_output(
            ["nvcc", "--version"],
            stderr=subprocess.STDOUT,
            timeout=10,
        ).decode()
        m = _NVCC_RELEASE_RE.search(out)

        if m:
            parsed = _parse_major_minor(m.group(1))
            if parsed is not None:
                return parsed
    except FileNotFoundError:
        pass
    except Exception as e:
        warnings.warn(f"nvcc probe failed: {e}", RuntimeWarning, stacklevel=3)

    return _probe_cuda_version_from_files()


def _build_arch_list(nvcc_ver: Optional[_NVCCVersion]) -> str:
    if nvcc_ver is None:
        warnings.warn(
            "Cannot determine nvcc version, "
            f"using conservative fallback: {_FALLBACK_ARCH_LIST}\n"
            "Override by setting TORCH_CUDA_ARCH_LIST manually.",
            RuntimeWarning,
            stacklevel=2,
        )
        return _FALLBACK_ARCH_LIST

    supported = [
        arch
        for arch, min_v, max_v in _ARCH_SUPPORT
        if nvcc_ver >= min_v and (max_v is None or nvcc_ver < max_v)
    ]

    if not supported:
        nvcc_ver_str = _format_nvcc_version(nvcc_ver)
        raise RuntimeError(
            f"No supported GPU architectures for nvcc {nvcc_ver_str}. "
            "Please update _ARCH_SUPPORT or set TORCH_CUDA_ARCH_LIST manually."
        )
    return ";".join(supported)


def compile(
    name: str,
    sources: List[str],
    extra_include_paths: List[str],
    build_directory: str,
) -> Any:
    arch_list = os.environ.get("TORCH_CUDA_ARCH_LIST")
    if arch_list is not None:
        _LOGGER.info(
            f"[compile_extension] Using user-specified TORCH_CUDA_ARCH_LIST: "
            f"{arch_list}"
        )
    else:
        nvcc_ver = _get_nvcc_version()
        arch_list = _build_arch_list(nvcc_ver)
        os.environ["TORCH_CUDA_ARCH_LIST"] = arch_list
        nvcc_ver_display = (
            _format_nvcc_version(nvcc_ver) if nvcc_ver is not None else "unknown"
        )
        _LOGGER.info(
            f"[compile_extension] nvcc={nvcc_ver_display}, "
            f"TORCH_CUDA_ARCH_LIST={arch_list}"
        )

    return load(
        name=name,
        sources=sources,
        extra_include_paths=extra_include_paths,
        extra_cflags=[
            "-O3",
            "-DVERSION_GE_1_1",
            "-DVERSION_GE_1_3",
            "-DVERSION_GE_1_5",
        ],
        extra_cuda_cflags=[
            "-O3",
            "--use_fast_math",
            "-DVERSION_GE_1_1",
            "-DVERSION_GE_1_3",
            "-DVERSION_GE_1_5",
            "-std=c++17",
            "-maxrregcount=32",
            "-U__CUDA_NO_HALF_OPERATORS__",
            "-U__CUDA_NO_HALF_CONVERSIONS__",
            "--expt-relaxed-constexpr",
            "--expt-extended-lambda",
        ],
        verbose=True,
        build_directory=build_directory,
    )
