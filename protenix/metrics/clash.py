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

import logging
from typing import Optional

import torch
import torch.nn as nn

from protenix.data.constants import rdkit_vdws

RDKIT_VDWS = torch.tensor(rdkit_vdws)
ID2TYPE = {0: "UNK", 1: "lig", 2: "prot", 3: "dna", 4: "rna"}


def get_vdw_radii(elements_one_hot):
    """get vdw radius for each atom according to their elements"""
    element_order = elements_one_hot.argmax(dim=1)
    return RDKIT_VDWS.to(element_order.device)[element_order]


class Clash(nn.Module):
    def __init__(
        self,
        af3_clash_threshold=1.1,
        vdw_clash_threshold=0.75,
        compute_af3_clash=True,
        compute_vdw_clash=True,
    ):
        super().__init__()
        self.af3_clash_threshold = af3_clash_threshold
        self.vdw_clash_threshold = vdw_clash_threshold
        self.compute_af3_clash = compute_af3_clash
        self.compute_vdw_clash = compute_vdw_clash

    def forward(
        self,
        pred_coordinate,
        asym_id,
        atom_to_token_idx,
        is_ligand,
        is_protein,
        is_dna,
        is_rna,
        mol_id: Optional[torch.Tensor] = None,
        elements_one_hot: Optional[torch.Tensor] = None,
    ):
        chain_info = self.get_chain_info(
            asym_id=asym_id,
            atom_to_token_idx=atom_to_token_idx,
            is_ligand=is_ligand,
            is_protein=is_protein,
            is_dna=is_dna,
            is_rna=is_rna,
            mol_id=mol_id,
            elements_one_hot=elements_one_hot,
        )
        return self._check_clash_per_chain_pairs(
            pred_coordinate=pred_coordinate, **chain_info
        )

    def get_chain_info(
        self,
        asym_id,
        atom_to_token_idx,
        is_ligand,
        is_protein,
        is_dna,
        is_rna,
        mol_id: Optional[torch.Tensor] = None,
        elements_one_hot: Optional[torch.Tensor] = None,
    ):
        # Get chain info
        asym_id = asym_id.long()
        asym_id_to_asym_mask = {
            aid.item(): asym_id == aid for aid in torch.unique(asym_id)
        }
        N_chains = len(asym_id_to_asym_mask)
        # Make sure it is from 0 to N_chain-1
        assert N_chains == asym_id.max() + 1

        # Check and compute chain_types
        chain_types = []
        mol_id_to_asym_ids, asym_id_to_mol_id = {}, {}
        atom_type = (1 * is_ligand + 2 * is_protein + 3 * is_dna + 4 * is_rna).long()
        if self.compute_vdw_clash:
            assert mol_id is not None
            assert elements_one_hot is not None

        for aid in range(N_chains):
            atom_chain_mask = asym_id_to_asym_mask[aid][atom_to_token_idx]
            atom_type_i = atom_type[atom_chain_mask]
            assert len(atom_type_i.unique()) == 1
            if atom_type_i[0].item() == 0:
                logging.warning(
                    "Unknown asym_id type: not in ligand / protein / dna / rna"
                )
            chain_types.append(ID2TYPE[atom_type_i[0].item()])
            if self.compute_vdw_clash:
                # Check if all atoms in a chain are from the same molecule
                mol_id_i = mol_id[atom_chain_mask].unique().item()
                mol_id_to_asym_ids.setdefault(mol_id_i, []).append(aid)
                asym_id_to_mol_id[aid] = mol_id_i

        chain_info = {
            "N_chains": N_chains,
            "atom_to_token_idx": atom_to_token_idx,
            "asym_id_to_asym_mask": asym_id_to_asym_mask,
            "atom_type": atom_type,
            "mol_id": mol_id,
            "elements_one_hot": elements_one_hot,
            "chain_types": chain_types,
        }

        if self.compute_vdw_clash:
            chain_info.update({"asym_id_to_mol_id": asym_id_to_mol_id})

        return chain_info

    def get_chain_pair_violations(
        self,
        pred_coordinate,
        violation_type,
        chain_1_mask,
        chain_2_mask,
        elements_one_hot: Optional[torch.Tensor] = None,
    ):
        chain_1_coords = pred_coordinate[chain_1_mask, :]
        chain_2_coords = pred_coordinate[chain_2_mask, :]
        pred_dist = torch.cdist(chain_1_coords, chain_2_coords)
        if violation_type == "af3":
            clash_per_atom_pair = (
                pred_dist < self.af3_clash_threshold
            )  # [ N_atom_chain_1, N_atom_chain_2]
            clashed_col, clashed_row = torch.where(clash_per_atom_pair)
            clash_atom_pairs = torch.stack((clashed_col, clashed_row), dim=-1)
        else:
            assert elements_one_hot is not None
            vdw_radii_i, vdw_radii_j = get_vdw_radii(
                elements_one_hot[chain_1_mask, :]
            ), get_vdw_radii(elements_one_hot[chain_2_mask, :])
            vdw_sum_pair = (
                vdw_radii_i[:, None] + vdw_radii_j[None, :]
            )  # [N_atom_chain_1, N_atom_chain_2]
            relative_vdw_distance = pred_dist / vdw_sum_pair
            clash_per_atom_pair = (
                relative_vdw_distance < self.vdw_clash_threshold
            )  # [N_atom_chain_1, N_atom_chain_2]
            clashed_col, clashed_row = torch.where(clash_per_atom_pair)
            clash_rel_dist = relative_vdw_distance[clashed_col, clashed_row]
            clashed_global_col = torch.where(chain_1_mask)[0][clashed_col]
            clashed_global_row = torch.where(chain_2_mask)[0][clashed_row]
            clash_atom_pairs = torch.stack(
                (clashed_global_col, clashed_global_row, clash_rel_dist), dim=-1
            )
        return clash_atom_pairs

    def _check_clash_per_chain_pairs(
        self,
        pred_coordinate,
        atom_to_token_idx,
        N_chains,
        atom_type,
        chain_types,
        elements_one_hot,
        asym_id_to_asym_mask,
        mol_id: Optional[torch.Tensor] = None,
        asym_id_to_mol_id: Optional[torch.Tensor] = None,
    ):
        device = pred_coordinate.device
        N_sample = pred_coordinate.shape[0]

        # initialize results
        if self.compute_af3_clash:
            has_af3_clash_flag = torch.zeros(
                N_sample, N_chains, N_chains, device=device, dtype=torch.bool
            )
            af3_clash_details = torch.zeros(
                N_sample, N_chains, N_chains, 2, device=device, dtype=torch.bool
            )
        if self.compute_vdw_clash:
            has_vdw_clash_flag = torch.zeros(
                N_sample, N_chains, N_chains, device=device, dtype=torch.bool
            )
            vdw_clash_details = {}

        skipped_pairs = []
        for sample_id in range(N_sample):
            for i in range(N_chains):
                if chain_types[i] == "UNK":
                    continue
                atom_chain_mask_i = asym_id_to_asym_mask[i][atom_to_token_idx]
                N_chain_i = torch.sum(atom_chain_mask_i).item()
                for j in range(i + 1, N_chains):
                    if chain_types[j] == "UNK":
                        continue
                    chain_pair_type = set([chain_types[i], chain_types[j]])
                    # Skip potential bonded ligand to polymers
                    skip_bonded_ligand = False
                    if (
                        self.compute_vdw_clash
                        and "lig" in chain_pair_type
                        and len(chain_pair_type) > 1
                        and asym_id_to_mol_id[i] == asym_id_to_mol_id[j]
                    ):
                        common_mol_id = asym_id_to_mol_id[i]
                        logging.warning(
                            f"mol_id {common_mol_id} may contain bonded ligand to polymers"
                        )
                        skip_bonded_ligand = True
                        skipped_pairs.append((i, j))
                    atom_chain_mask_j = asym_id_to_asym_mask[j][atom_to_token_idx]
                    N_chain_j = torch.sum(atom_chain_mask_j).item()
                    if self.compute_vdw_clash and not skip_bonded_ligand:
                        vdw_clash_pairs = self.get_chain_pair_violations(
                            pred_coordinate=pred_coordinate[sample_id, :, :],
                            violation_type="vdw",
                            chain_1_mask=atom_chain_mask_i,
                            chain_2_mask=atom_chain_mask_j,
                            elements_one_hot=elements_one_hot,
                        )
                        if vdw_clash_pairs.shape[0] > 0:
                            vdw_clash_details[(sample_id, i, j)] = vdw_clash_pairs
                            has_vdw_clash_flag[sample_id, i, j] = True
                            has_vdw_clash_flag[sample_id, j, i] = True
                    if (
                        chain_types[i] == "lig" or chain_types[j] == "lig"
                    ):  # AF3 clash only consider polymer chains
                        continue
                    if self.compute_af3_clash:
                        af3_clash_pairs = self.get_chain_pair_violations(
                            pred_coordinate=pred_coordinate[sample_id, :, :],
                            violation_type="af3",
                            chain_1_mask=atom_chain_mask_i,
                            chain_2_mask=atom_chain_mask_j,
                        )
                        total_clash = af3_clash_pairs.shape[0]
                        relative_clash = total_clash / min(N_chain_i, N_chain_j)
                        af3_clash_details[sample_id, i, j, 0] = total_clash
                        af3_clash_details[sample_id, i, j, 1] = relative_clash
                        has_af3_clash_flag[sample_id, i, j] = (
                            total_clash > 100 or relative_clash > 0.5
                        )
                        af3_clash_details[sample_id, j, i, :] = af3_clash_details[
                            sample_id, i, j, :
                        ]
                        has_af3_clash_flag[sample_id, j, i] = has_af3_clash_flag[
                            sample_id, i, j
                        ]
        return {
            "summary": {
                "af3_clash": has_af3_clash_flag if self.compute_af3_clash else None,
                "vdw_clash": has_vdw_clash_flag if self.compute_vdw_clash else None,
                "chain_types": chain_types,
                "skipped_pairs": skipped_pairs,
            },
            "details": {
                "af3_clash": af3_clash_details if self.compute_af3_clash else None,
                "vdw_clash": vdw_clash_details if self.compute_vdw_clash else None,
            },
        }
