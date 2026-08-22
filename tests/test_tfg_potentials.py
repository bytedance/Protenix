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

import torch

from protenix.tfg.potentials import _flat_bottom_parabolic, PairwiseDistancePotential


def test_flat_bottom_parabolic_with_finite_bounds():
    value = torch.tensor([0.0, 1.5, 2.0, 2.5, 4.0], requires_grad=True)
    k = torch.full_like(value, 2.0)
    lower = torch.tensor([1.0])
    upper = torch.tensor([3.0])

    energy, analytic_grad = _flat_bottom_parabolic(value, k, lower, upper)
    (autograd_grad,) = torch.autograd.grad(energy.sum(), value)

    torch.testing.assert_close(energy, torch.tensor([1.0, 0.0, 0.0, 0.0, 1.0]))
    torch.testing.assert_close(analytic_grad, torch.tensor([-2.0, 0.0, 0.0, 0.0, 2.0]))
    torch.testing.assert_close(autograd_grad, analytic_grad)


def test_pairwise_clash_autograd_matches_analytic_gradient():
    coords = torch.tensor([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]], requires_grad=True)
    ref_element = torch.zeros(2, 118)
    ref_element[:, 5] = 1.0
    feats = {
        "pairwise_distance_index": torch.tensor([[0], [1]]),
        "pairwise_distance_lower_bound": torch.tensor([2.0]),
        "pairwise_distance_upper_bound": torch.tensor([3.0]),
        "pairwise_distance_is_bond": torch.tensor([0]),
        "pairwise_distance_is_angle": torch.tensor([0]),
        "ref_element": ref_element,
    }
    potential = PairwiseDistancePotential(
        default_params={
            "bond_buffer": 0.0,
            "angle_buffer": 0.0,
            "clash_buffer": 0.0,
        }
    )

    energy = potential.energy(coords, feats)
    (autograd_grad,) = torch.autograd.grad(energy.sum(), coords)
    analytic_energy, analytic_grad = potential.energy_and_grad(coords.detach(), feats)

    assert energy.item() > 0
    assert torch.isfinite(autograd_grad).all()
    torch.testing.assert_close(energy, analytic_energy)
    torch.testing.assert_close(autograd_grad, analytic_grad)
