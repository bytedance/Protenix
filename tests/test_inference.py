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

import unittest

import numpy as np
import torch
from biotite.structure import AtomArray

from runner.inference import fix_cterminal_carboxyl_oxygens


class TestFixCterminalCarboxylOxygens(unittest.TestCase):
    def test_fix_cterminal_carboxyl_oxygens(self):
        atom_array = AtomArray(13)
        atom_array.chain_id = np.array(["A"] * 9 + ["B"] * 4)
        atom_array.res_id = np.array([1] * 4 + [2] * 5 + [7] * 4)
        atom_array.atom_name = np.array(
            [
                "CA",
                "C",
                "O",
                "OXT",
                "N",
                "CA",
                "C",
                "O",
                "OXT",
                "CA",
                "C",
                "O",
                "OXT",
            ]
        )

        coordinates = torch.arange(2 * 13 * 3, dtype=torch.float64).reshape(
            2, 13, 3
        )

        # Chain A uses the x-axis through C; each sample keeps a different oxygen.
        coordinates[0, 5:9] = torch.tensor(
            [
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [10.0, 0.0, 0.0],
            ]
        )
        coordinates[1, 5:9] = torch.tensor(
            [
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [0.25, 0.0, 2.0],
            ]
        )

        # Chain B uses the y-axis through C and has a translated origin.
        coordinates[0, 9:13] = torch.tensor(
            [
                [1.0, 2.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 1.25, 1.0],
                [8.0, 8.0, 8.0],
            ]
        )
        coordinates[1, 9:13] = torch.tensor(
            [
                [-1.0, 3.0, 3.0],
                [-1.0, 2.0, 3.0],
                [9.0, 9.0, 9.0],
                [-1.0, 2.5, 4.5],
            ]
        )

        original = coordinates.clone()
        result = fix_cterminal_carboxyl_oxygens(coordinates, atom_array)

        expected = original.clone()
        expected[0, 7] = torch.tensor([0.5, 1.0, 0.0])
        expected[0, 8] = torch.tensor([0.5, -1.0, 0.0])
        expected[1, 7] = torch.tensor([0.25, 0.0, 2.0])
        expected[1, 8] = torch.tensor([0.25, 0.0, -2.0])
        expected[0, 11] = torch.tensor([2.0, 1.25, 1.0])
        expected[0, 12] = torch.tensor([0.0, 1.25, 1.0])
        expected[1, 11] = torch.tensor([-1.0, 2.5, 4.5])
        expected[1, 12] = torch.tensor([-1.0, 2.5, 1.5])

        torch.testing.assert_close(result, expected)
        torch.testing.assert_close(coordinates, original)
        self.assertEqual(result.dtype, coordinates.dtype)


if __name__ == "__main__":
    unittest.main()
