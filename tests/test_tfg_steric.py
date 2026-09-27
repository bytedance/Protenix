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

import pytest
import torch

from protenix.tfg import parse_tfg_config, TFGEngine
from protenix.tfg.potentials import VinaStericPotential


def make_features(chain_ids):
    ref_element = torch.zeros((len(chain_ids), 128))
    ref_element[:, 6] = 1  # Carbon atoms.
    return {
        "asym_id": torch.tensor(chain_ids, dtype=torch.long),
        "atom_to_token_idx": torch.arange(len(chain_ids)),
        "ref_element": ref_element,
        "interchain_bond_index": torch.empty((2, 0), dtype=torch.long),
    }


@pytest.mark.parametrize("chain_ids", [[], [0, 0, 0], [7, 7, 7]])
@pytest.mark.parametrize("batch_shape", [(), (2, 3)])
@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_steric_energy_without_interchain_pairs(chain_ids, batch_shape, dtype):
    coords = torch.zeros((*batch_shape, len(chain_ids), 3), dtype=dtype)
    feats = make_features(chain_ids)
    expected_energy = torch.zeros(batch_shape, dtype=dtype)

    # Exercise both public APIs on a cold cache, then again on the cached result.
    energy_only = VinaStericPotential()
    with_gradient = VinaStericPotential()
    for _ in range(2):
        torch.testing.assert_close(energy_only.energy(coords, feats), expected_energy)
        energy, gradient = with_gradient.energy_and_grad(coords, feats)
        torch.testing.assert_close(energy, expected_energy)
        torch.testing.assert_close(gradient, torch.zeros_like(coords))


def test_steric_cache_can_transition_between_single_and_multiple_chains():
    potential = VinaStericPotential()
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.0, 0.0], [0.5, 1.0, 0.0]]
    )
    single_chain = make_features([0, 0, 0, 0])
    multiple_chains = make_features([0, 0, 1, 1])
    torch.testing.assert_close(
        potential.energy(coords, single_chain), torch.tensor(0.0)
    )

    energy, gradient = potential.energy_and_grad(coords, multiple_chains)
    expected_energy, expected_gradient = VinaStericPotential().energy_and_grad(
        coords, multiple_chains
    )
    assert energy > 0
    assert torch.isfinite(gradient).all()
    assert torch.count_nonzero(gradient) > 0
    torch.testing.assert_close(energy, expected_energy)
    torch.testing.assert_close(gradient, expected_gradient)
    torch.testing.assert_close(
        potential.energy(coords, single_chain), torch.tensor(0.0)
    )


def test_tfg_step_leaves_single_chain_coordinates_unchanged_for_steric_only_guidance():
    config = parse_tfg_config(
        {
            "enable": True,
            "mu": 0.1,
            "steps": {"tfg_inner": 2, "projection_outer": 0},
            "terms": {"VinaStericPotential": {"weight": 1.0}},
        }
    )
    coords = torch.arange(18, dtype=torch.float32).reshape(2, 3, 3)
    engine = TFGEngine(config, device=coords.device, dtype=coords.dtype)

    def identity_denoiser(*, x_noisy, **kwargs):
        return x_noisy

    result = engine.step(
        identity_denoiser,
        x=coords,
        t_hat=torch.ones(2),
        c_tau=torch.zeros(2),
        step_scale_eta=1.0,
        step_i=0,
        num_diffusion_steps=1,
        input_feature_dict=make_features([0, 0, 0]),
        s_inputs=None,
        s_trunk=None,
        z_trunk=None,
        pair_z=None,
        p_lm=None,
        c_l=None,
        chunk_size=None,
        inplace_safe=False,
        enable_efficient_fusion=False,
    )
    torch.testing.assert_close(result, coords)
