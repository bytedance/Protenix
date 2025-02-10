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

from protenix.model.modules.frames import (
    expressCoordinatesInFrame,
    gather_frame_atom_by_indices,
)


class TestFrame(unittest.TestCase):
    def setUp(self):
        self._start_time = time.time()

    def test_express_coordinates_in_frame(self):
        N_atom = 10
        N_frame = 5
        bs_dims = (2, 3)
        coordinates = torch.rand(size=(*bs_dims, N_atom, 3))
        frames = torch.rand(size=(*bs_dims, N_frame, 3, 3))
        x_transformed = expressCoordinatesInFrame(coordinate=coordinates, frames=frames)
        # shape check
        self.assertEqual(
            x_transformed.shape, torch.Size((*bs_dims, N_frame, N_atom, 3))
        )

        # invarient to batch order
        x_transformed12 = expressCoordinatesInFrame(
            coordinate=coordinates[1][2], frames=frames[1][2]
        )
        self.assertTrue(torch.allclose(x_transformed12, x_transformed[1][2]))

        # value check
        frames = torch.tensor([[[1, 0, 0], [0, 0, 0], [0, 1, 0]]]).float()
        coordinates = torch.tensor([[1, 0, 0]]).float()
        x_transformed = expressCoordinatesInFrame(coordinate=coordinates, frames=frames)
        self.assertTrue(
            torch.allclose(
                x_transformed,
                torch.tensor([[0.7071, -0.7071, 0.0000]]),
                atol=1e-3,
                rtol=1e-3,
            )
        )  # math.sqrt(2)/2

    def test_gather_frame_atom_by_indices(self):
        N_atom = 10
        N_frame = 5
        bs_dims = (2, 3)
        coordinates = torch.rand(size=(*bs_dims, N_atom, 3))
        indexes = torch.randint(size=(*bs_dims, N_frame, 3), low=0, high=10)
        out = gather_frame_atom_by_indices(
            coordinate=coordinates, frame_atom_index=indexes
        )
        # shape check
        self.assertEqual(out.shape, torch.Size((*bs_dims, N_frame, 3, 3)))
        coordinates = torch.rand(size=(*bs_dims, N_atom, 3))
        indexes = torch.randint(size=(N_frame, 3), low=0, high=10)
        out = gather_frame_atom_by_indices(
            coordinate=coordinates, frame_atom_index=indexes
        )
        # shape check [naive mode]
        self.assertEqual(out.shape, torch.Size((*bs_dims, N_frame, 3, 3)))

    def tearDown(self):
        elapsed_time = time.time() - self._start_time
        print(f"Test {self.id()} took {elapsed_time:.6f}s")


if __name__ == "__main__":
    unittest.main()
