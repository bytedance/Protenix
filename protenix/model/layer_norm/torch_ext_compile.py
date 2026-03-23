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
from typing import Any, Optional

from torch.utils.cpp_extension import load


def compile(
    name: str,
    sources: list[str],
    extra_include_paths: list[str],
    build_directory: Optional[str] = None,
) -> Any:
    # Query supported architectures from nvcc (resolved via PyTorch's
    # CUDA_HOME so we use the same toolchain as cpp_extension.load).
    import re
    import shutil
    import subprocess

    from torch.utils.cpp_extension import CUDA_HOME

    _nvcc = shutil.which("nvcc")
    if CUDA_HOME:
        _candidate = os.path.join(CUDA_HOME, "bin", "nvcc")
        if os.path.isfile(_candidate):
            _nvcc = _candidate

    _supported = set()
    try:
        out = subprocess.check_output(
            [_nvcc, "--list-gpu-arch"], text=True, stderr=subprocess.STDOUT
        )
        _supported = set(re.findall(r"compute_(\d+)", out))
    except Exception:
        _supported = {"70", "80", "86", "90"}  # safe defaults

    _wanted = [
        ("70", "70"),
        ("80", "80"),
        ("86", "86"),
        ("89", "89"),
        ("90", "90"),
        ("100", "100"),
    ]
    gencode_flags = []
    for compute, sm in _wanted:
        if compute in _supported:
            gencode_flags += ["-gencode", f"arch=compute_{compute},code=sm_{sm}"]
    if not gencode_flags:
        gencode_flags = ["-gencode", "arch=compute_80,code=sm_80"]

    # Build TORCH_CUDA_ARCH_LIST dynamically from supported architectures
    _arch_list = [f"{int(c)//10}.{int(c)%10}" for c, _ in _wanted if c in _supported]
    os.environ["TORCH_CUDA_ARCH_LIST"] = ";".join(_arch_list) if _arch_list else "8.0"

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
        ]
        + gencode_flags,
        verbose=True,
        build_directory=build_directory,
    )
