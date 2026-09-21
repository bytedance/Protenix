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
import unittest
from unittest.mock import patch


_previous_layer_norm_type = os.environ.get("LAYERNORM_TYPE")
os.environ["LAYERNORM_TYPE"] = "torch"

from protenix.model.triangular.layers import _use_fast_layer_norm

if _previous_layer_norm_type is None:
    del os.environ["LAYERNORM_TYPE"]
else:
    os.environ["LAYERNORM_TYPE"] = _previous_layer_norm_type


class TestLayerNormSelection(unittest.TestCase):
    def test_falls_back_to_torch_on_sm120(self) -> None:
        with (
            patch.dict(os.environ, {}, clear=True),
            patch("torch.cuda.is_available", return_value=True),
            patch("torch.cuda.get_device_capability", return_value=(12, 0)),
        ):
            self.assertFalse(_use_fast_layer_norm())

    def test_keeps_fast_layer_norm_on_other_cuda_architectures(self) -> None:
        with (
            patch.dict(os.environ, {}, clear=True),
            patch("torch.cuda.is_available", return_value=True),
            patch("torch.cuda.get_device_capability", return_value=(9, 0)),
        ):
            self.assertTrue(_use_fast_layer_norm())

    def test_explicit_fast_layer_norm_overrides_sm120_fallback(self) -> None:
        with (
            patch.dict(os.environ, {"LAYERNORM_TYPE": "fast_layernorm"}, clear=True),
            patch("torch.cuda.is_available", return_value=True),
            patch("torch.cuda.get_device_capability", return_value=(12, 0)),
        ):
            self.assertTrue(_use_fast_layer_norm())

    def test_explicit_torch_layer_norm_is_preserved(self) -> None:
        with patch.dict(os.environ, {"LAYERNORM_TYPE": "torch"}, clear=True):
            self.assertFalse(_use_fast_layer_norm())


if __name__ == "__main__":
    unittest.main()
