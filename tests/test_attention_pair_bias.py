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

from protenix.model.modules.transformer import AttentionPairBias


class TestAttentionPairBias(unittest.TestCase):
    def setUp(self) -> None:
        self._start_time = time.time()
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        super().setUp()

    def get_model(
        self,
        has_s=True,
        n_heads: int = 16,
        c_a: int = 768,
        c_s: int = 384,
        c_z: int = 128,
    ):

        model = AttentionPairBias(
            has_s=has_s, n_heads=n_heads, c_a=c_a, c_s=c_s, c_z=c_z
        ).to(self.device)

        return model

    def test_shape(self) -> None:
        """
        Args:
            a (torch.Tensor): the single feature aggregate per-atom representation
                [..., N_token, c_a]
            s (torch.Tensor): single embedding
                [..., N_token, c_s]
            z (torch.Tensor): pair embedding
                [..., N_token, N_token, c_z]
            n_queries (int, optional): local window size of query tensor. If not None, will perform local attention. Defaults to None.
            n_keys (int, optional): local window size of key tensor. Defaults to None.

        Returns:
            torch.Tensor: the updated a from AttentionPairBias
                [..., N_token, c_a]
        """
        n_heads = 3
        c_a = 3 * 55
        c_s = 23
        c_z = 17

        N_token = 135
        bs_dims = (2, 3)

        inputs = {
            "a": torch.rand(size=(*bs_dims, N_token, c_a)).to(self.device),
            "s": torch.rand(size=(*bs_dims, N_token, c_s)).to(self.device),
            "z": torch.rand(size=(*bs_dims, N_token, N_token, c_z)).to(self.device),
        }

        model = self.get_model(c_a=c_a, c_s=c_s, c_z=c_z, n_heads=n_heads)

        out = model(**inputs)
        target_shape = (*bs_dims, N_token, c_a)
        self.assertEqual(out.shape, out.reshape(target_shape).shape)

    def test_local_attention_shape(self) -> None:
        """Used by Algorithm 24, with beta_ij being the local mask. Used in AtomTransformer.

        Args:
            a (torch.Tensor): atom embedding
                [..., N_atom, c_a]
            s (torch.Tensor): atom embedding
                [..., N_atom, c_s]
            z (torch.Tensor): atom-atom pair embedding, in trunked dense shape. Used for computing pair bias.
                [..., n_blocks, n_queries, n_keys, c_z]
            n_queries (int, optional): local window size of query tensor. Defaults to 32.
            n_keys (int, optional): local window size of key tensor. Defaults to 128.

        Returns:
            torch.Tensor: the updated a from AttentionPairBias
                [..., N_atom, c_a]
        """
        n_heads = 3
        c_a = 3 * 27
        c_s = 23
        c_z = 17

        N_token = 128 * 2 + 45

        bs_dims = (2, 3)

        N_q = 32
        N_k = 128
        N_blocks = math.ceil(N_token / N_q)

        inputs = {
            "a": torch.rand(size=(*bs_dims, N_token, c_a)).to(self.device),
            "s": torch.rand(size=(*bs_dims, N_token, c_s)).to(self.device),
            "z": torch.rand(size=(*bs_dims, N_blocks, N_q, N_k, c_z)).to(self.device),
            "n_queries": 32,
            "n_keys": 128,
        }

        model = self.get_model(c_a=c_a, c_s=c_s, c_z=c_z, n_heads=n_heads)

        out = model(**inputs)
        target_shape = (*bs_dims, N_token, c_a)
        self.assertEqual(out.shape, out.reshape(target_shape).shape)

    def tearDown(self):
        elapsed_time = time.time() - self._start_time
        print(f"Test {self.id()} took {elapsed_time:.6f}s")


if __name__ == "__main__":
    unittest.main()
