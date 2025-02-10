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

import math
import time
import unittest

import torch

from protenix.model.modules.transformer import DiffusionTransformer


class TestDiffusionTransformer(unittest.TestCase):
    def setUp(self) -> None:
        self._start_time = time.time()
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        super().setUp()

    def get_model(
        self,
        c_a: int = 128,
        c_s: int = 384,
        c_z: int = 64,
        n_blocks: int = 3,
        n_heads: int = 4,
    ):

        model = DiffusionTransformer(
            c_a=c_a, c_s=c_s, c_z=c_z, n_blocks=n_blocks, n_heads=n_heads
        ).to(self.device)

        return model

    def test_shape(self) -> None:

        n_heads = 2
        c_a = 13 * n_heads
        c_s = 23
        c_z = 17

        N = 45
        bs_dims = (2, 3)

        inputs = {
            "a": torch.rand(size=(*bs_dims, N, c_a)).to(self.device),
            "s": torch.rand(size=(*bs_dims, N, c_s)).to(self.device),
            "z": torch.rand(size=(*bs_dims, N, N, c_z)).to(self.device),
            "n_queries": None,
            "n_keys": None,
        }

        model = self.get_model(c_a=c_a, c_s=c_s, c_z=c_z, n_heads=n_heads)

        out = model(**inputs)
        target_shape = (*bs_dims, N, c_a)
        self.assertEqual(out.shape, out.reshape(target_shape).shape)

        N_q = 32
        N_k = 128
        N_blocks = math.ceil(N / N_q)

        inputs = {
            "a": torch.rand(size=(*bs_dims, N, c_a)).to(self.device),
            "s": torch.rand(size=(*bs_dims, N, c_s)).to(self.device),
            "z": torch.rand(size=(*bs_dims, N_blocks, N_q, N_k, c_z)).to(self.device),
            "n_queries": 32,
            "n_keys": 128,
        }

        out = model(**inputs)
        target_shape = (*bs_dims, N, c_a)
        self.assertEqual(out.shape, out.reshape(target_shape).shape)

    def tearDown(self):
        elapsed_time = time.time() - self._start_time
        print(f"Test {self.id()} took {elapsed_time:.6f}s")


if __name__ == "__main__":
    unittest.main()
