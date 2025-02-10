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

import time
import unittest

import torch

from protenix.model.modules.transformer import ConditionedTransitionBlock


class TestConditionedTransitionBlock(unittest.TestCase):
    def setUp(self) -> None:
        self._start_time = time.time()
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        super().setUp()

    def get_model(self, c_a: int = 768, c_s: int = 384, n: int = 2):

        model = ConditionedTransitionBlock(c_a=c_a, c_s=c_s, n=n).to(self.device)

        return model

    def test_shape(self) -> None:

        c_a = 5 * 55
        c_s = 123

        N_token = 135
        bs_dims = (2, 3, 5)

        inputs = {
            "a": torch.rand(size=(*bs_dims, N_token, c_a)).to(self.device),
            "s": torch.rand(size=(*bs_dims, N_token, c_s)).to(self.device),
        }

        model = self.get_model(c_a=c_a, c_s=c_s)

        out = model(**inputs)
        target_shape = (*bs_dims, N_token, c_a)
        self.assertEqual(out.shape, out.reshape(target_shape).shape)

    def tearDown(self):
        elapsed_time = time.time() - self._start_time
        print(f"Test {self.id()} took {elapsed_time:.6f}s")


if __name__ == "__main__":
    unittest.main()
