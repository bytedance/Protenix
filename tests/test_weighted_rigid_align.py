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

from protenix.model.loss import weighted_rigid_align
from protenix.utils.seed import seed_everything


class TestWeightedRigidAlign(unittest.TestCase):

    def setUp(self) -> None:
        self._start_time = time.time()
        self.device = "cpu"
        super().setUp()

    def test_batch(self):
        batch_dims = (4,)
        num_atom = 56

        x = (
            torch.randn(size=(*batch_dims, num_atom, 3))
            + torch.rand(size=(*batch_dims, num_atom, 3))
        ).to(self.device)
        x_target = (
            torch.randn(size=(*batch_dims, num_atom, 3))
            + torch.rand(size=(*batch_dims, num_atom, 3))
        ).to(self.device)

        seed_everything(42, True)
        weight = torch.rand(size=(*batch_dims, num_atom)).to(self.device)
        mask = torch.rand(size=(*batch_dims, num_atom)).to(self.device) > 0.85
        mask2 = torch.rand(size=(*batch_dims, num_atom)).to(self.device) > 0.5

        # check self-self alignment
        x_x_aligned = weighted_rigid_align(
            x=x, x_target=x, atom_weight=torch.ones_like(x[..., 0])
        )
        self.assertTrue(torch.allclose(x, x_x_aligned, atol=1e-5))

        x_x_aligned = weighted_rigid_align(x=x, x_target=x, atom_weight=mask.float())
        self.assertTrue(torch.allclose(x, x_x_aligned, atol=1e-6))

        x_x_aligned = weighted_rigid_align(
            x=x, x_target=x, atom_weight=mask.float() + mask2.float()
        )
        self.assertTrue(
            torch.allclose(
                x * (mask + mask2).float()[..., None],
                x_x_aligned * (mask + mask2).float()[..., None],
                atol=1e-6,
            )
        )

        # check batch order
        out_batch = weighted_rigid_align(
            x=x, x_target=x_target, atom_weight=weight * mask
        )
        out1 = weighted_rigid_align(
            x=x[1][None, ...],
            x_target=x_target[1][None, ...],
            atom_weight=weight[1][None, ...] * mask[1][None, ...],
        )
        out2 = weighted_rigid_align(
            x=x[2][None, ...],
            x_target=x_target[2][None, ...],
            atom_weight=weight[2][None, ...] * mask[2][None, ...],
        )
        self.assertTrue(torch.allclose(out_batch[1], out1[0]))
        self.assertTrue(torch.allclose(out_batch[2], out2[0]))

        # more dims
        batch_dims = (2, 4)
        num_atom = 56

        x = (
            torch.randn(size=(*batch_dims, num_atom, 3))
            + torch.rand(size=(*batch_dims, num_atom, 3))
        ).to(self.device)
        x_target = (
            torch.randn(size=(*batch_dims, num_atom, 3))
            + torch.rand(size=(*batch_dims, num_atom, 3))
        ).to(self.device)
        weight = torch.rand(size=(*batch_dims, num_atom)).to(self.device)
        mask = torch.rand(size=(*batch_dims, num_atom)).to(self.device) > 0.85

        out_batch = weighted_rigid_align(
            x=x,
            x_target=x_target,
            atom_weight=weight,  # * mask
        )
        out11 = weighted_rigid_align(
            x=x[1, 1][None, ...],
            x_target=x_target[1, 1][None, ...],
            atom_weight=weight[1, 1][None, ...],  # * mask[1, 1][None, ...],
        )
        out03 = weighted_rigid_align(
            x=x[0, 3][None, ...],
            x_target=x_target[0, 3][None, ...],
            atom_weight=weight[0, 3][None, ...],  # * mask[0, 3][None, ...],
        )

        self.assertTrue(torch.allclose(out_batch[1, 1], out11[0]))
        self.assertTrue(torch.allclose(out_batch[0, 3], out03[0]))

    def tearDown(self):
        elapsed_time = time.time() - self._start_time
        print(f"Test {self.id()} took {elapsed_time:.6f}s")


if __name__ == "__main__":
    unittest.main()
