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
from typing import Any

DIRECT_TEMPLATE_SUFFIXES = {".cif", ".mmcif"}
SEARCH_TEMPLATE_SUFFIXES = {".a3m", ".hhr"}
SUPPORTED_INFERENCE_TEMPLATE_SUFFIXES = (
    DIRECT_TEMPLATE_SUFFIXES | SEARCH_TEMPLATE_SUFFIXES | {".json"}
)


def template_paths_need_mmcif_dir(inputs: list[dict[str, Any]]) -> bool:
    """Return whether configured template paths need a PDB mmCIF mirror."""
    for sample in inputs:
        for sequence in sample.get("sequences", []):
            protein_chain = sequence.get("proteinChain")
            if not protein_chain:
                continue
            template_path = protein_chain.get("templatesPath")
            if not template_path:
                continue
            suffix = os.path.splitext(str(template_path))[1].lower()
            if suffix in SEARCH_TEMPLATE_SUFFIXES:
                return True
    return False
