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

# pylint: disable=C0114
import os

from protenix.config.extend_types import ListValue, RequiredValue

current_file_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file_path)
code_directory = os.path.dirname(current_directory)
# The model will be download to the following dir if not exists:
# "./release_data/checkpoint/model_v0.5.0.pt"
inference_configs = {
    "model_name": "protenix_base_default_v0.5.0",  # inference model selection
    "seeds": ListValue([101]),
    "dump_dir": "./output",
    "need_atom_confidence": False,
    "sorted_by_ranking_score": True,
    "input_json_path": RequiredValue(str),
    "load_checkpoint_dir": os.path.join(code_directory, "./release_data/checkpoint/"),
    "num_workers": 16,
    "use_msa": True,
    "enable_tf32": True,
    "enable_efficient_fusion": True,
    "enable_diffusion_shared_vars_cache": True,
}
