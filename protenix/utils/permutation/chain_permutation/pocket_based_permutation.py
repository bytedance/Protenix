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

from protenix.metrics.rmsd import rmsd
from protenix.model.utils import expand_at_dim
from protenix.utils.distributed import traverse_and_aggregate
from protenix.utils.logger import get_logger
from protenix.utils.permutation.chain_permutation.utils import (
    apply_transform,
    get_optimal_transform,
    num_unique_matches,
)

logger = get_logger(__name__)


def permute_pred_to_optimize_pocket_aligned_rmsd(
    pred_coord: torch.Tensor,  # [N_sample, N_atom, 3]
    true_coord: torch.Tensor,  # [N_atom, 3]
    true_coord_mask: torch.Tensor,
    true_pocket_mask: torch.Tensor,
    true_ligand_mask: torch.Tensor,
    atom_entity_id: torch.Tensor,  # [N_atom]
    atom_asym_id: torch.Tensor,  # [N_atom]
    mol_atom_index: torch.Tensor,  # [N_atom]
    use_center_rmsd: bool = False,
):
    """

    Returns:
        permute_pred_indices (list[torch.Tensor]): A list of LongTensor.
            The list contains N_sample elements.
            Each elements is a LongTensor of shape = [N_atom].
        permuted_aligned_pred_coord (torch.Tensor): permuted and aligned coordinates of pred_coord.
            [N_sample, N_atom, 3]
    """

    log_dict = {}
    atom_entity_id = atom_entity_id.long()
    atom_asym_id = atom_asym_id.long()
    mol_atom_index = mol_atom_index.long()
    true_coord_mask = true_coord_mask.bool()
    true_pocket_mask = true_pocket_mask.bool()
    true_ligand_mask = true_ligand_mask.bool()
    assert pred_coord.size(-2) == true_coord.size(-2), "Atom numbers are difference."
    assert pred_coord.dim() == 3

    # find entity_id/asym_id of pocket and ligand chains
    def _get_entity_and_asym_id(atom_mask):

        masked_asym_id = atom_asym_id[atom_mask]
        masked_entity_id = atom_entity_id[atom_mask]
        assert (masked_asym_id[0] == masked_asym_id).all()
        assert (masked_entity_id[0] == masked_entity_id).all()
        return masked_asym_id[0].item(), masked_entity_id[0].item()

    pocket_asym_id, pocket_entity_id = _get_entity_and_asym_id(true_pocket_mask)
    ligand_asym_id, ligand_entity_id = _get_entity_and_asym_id(true_ligand_mask)

    candidate_pockets = {}
    for i in torch.unique(atom_asym_id[atom_entity_id == pocket_entity_id]):
        i = i.item()
        pocket_mask = atom_asym_id == i
        pocket_mask = pocket_mask * torch.isin(
            mol_atom_index, mol_atom_index[true_pocket_mask]
        )
        assert pocket_mask.sum() == true_pocket_mask.sum()
        candidate_pockets[i] = pocket_mask.clone()
    candidate_ligands = {}
    for j in torch.unique(atom_asym_id[atom_entity_id == ligand_entity_id]):
        j = j.item()
        lig_mask_j = atom_asym_id == j
        if lig_mask_j.sum() != true_ligand_mask.sum():
            logger.warning(
                f"The ligand selected by 'mol_id' has {lig_mask_j.sum().item()} atoms."
                + f"The true ligand selected by 'asym_id' has {true_ligand_mask.sum().item()} atoms."
            )
            lig_mask_j = lig_mask_j * torch.isin(
                mol_atom_index, mol_atom_index[true_ligand_mask]
            )
        assert lig_mask_j.sum() == true_ligand_mask.sum()
        candidate_ligands[j] = lig_mask_j

    log_dict["num_sym_pocket"] = len(candidate_pockets)
    log_dict["num_sym_ligand"] = len(candidate_ligands)
    log_dict["has_sym_chain"] = len(candidate_ligands) + len(candidate_pockets) > 2

    # Enumerate over the batch dimension of pred_coord
    # to find the optimal chain assignment for each sample.

    def _find_protein_ligand_chains_for_one_sample(
        coord: torch.Tensor,
    ):
        best_results = {}
        unpermuted_results = {}
        for poc_asym_id, pocket_mask in candidate_pockets.items():
            # Align pocket_i to true pocket
            rot, trans = get_optimal_transform(
                src_atoms=coord[pocket_mask].clone(),
                tgt_atoms=true_coord[true_pocket_mask],
                mask=true_coord_mask[true_pocket_mask],
            )
            # Transform predicted coordinates according to the aligment results
            aligned_pred_coord = apply_transform(coord.clone(), rot=rot, trans=trans)

            # Find the best ligand
            ordered_lig_asym_ids = [i for i in candidate_ligands]
            orderd_lig_masks = [candidate_ligands[i] for i in ordered_lig_asym_ids]
            aligned_lig_coords = torch.stack(
                [aligned_pred_coord[m] for m in orderd_lig_masks], dim=0
            )  # [N_lig, N_lig_atom, 3]

            if use_center_rmsd:
                mask = true_coord_mask[true_ligand_mask].bool()  # [N_lig_atom]
                aligned_lig_center = aligned_lig_coords[:, mask, :].mean(
                    dim=-2, keepdim=True
                )  # [N_lig, 1, 3]
                true_coord_center = true_coord[true_ligand_mask][mask, :].mean(
                    dim=-2, keepdim=True
                )  # [1, 3]
                per_lig_rmsd = rmsd(
                    aligned_lig_center,  # [N_lig, 1, 3]
                    expand_at_dim(
                        true_coord_center,
                        dim=0,
                        n=aligned_lig_coords.size(0),
                    ),
                    reduce=False,
                )  # [N_lig]
            else:
                per_lig_rmsd = rmsd(
                    aligned_lig_coords,
                    expand_at_dim(
                        true_coord[true_ligand_mask],
                        dim=0,
                        n=aligned_lig_coords.size(0),
                    ),
                    mask=true_coord_mask[true_ligand_mask],
                    reduce=False,
                )  # [N_lig]
            lig_rmsd, idx = per_lig_rmsd.min(dim=0)
            lig_asym_id = ordered_lig_asym_ids[idx]

            if lig_rmsd < best_results.get("rmsd", torch.inf):
                best_results = {
                    "rmsd": lig_rmsd,
                    "pocket_asym_id": poc_asym_id,
                    "ligand_asym_id": lig_asym_id,
                    "aligned_pred_coord": aligned_pred_coord,
                }
            if poc_asym_id == pocket_asym_id:
                # record the unpermuted result
                i = ordered_lig_asym_ids.index(ligand_asym_id)
                unpermuted_lig_rmsd = per_lig_rmsd[i].item()
                unpermuted_results = {
                    "rmsd": unpermuted_lig_rmsd,
                    "aligned_pred_coord": aligned_pred_coord,
                }

        # record stats
        per_sample_log_dict = {
            "is_permuted": best_results["pocket_asym_id"] != pocket_asym_id
            or best_results["ligand_asym_id"] != ligand_asym_id,
            "is_permuted_pocket": best_results["pocket_asym_id"] != pocket_asym_id,
            "is_permuted_ligand": best_results["ligand_asym_id"] != ligand_asym_id,
            "algo:no_permute": best_results["pocket_asym_id"] == pocket_asym_id
            and best_results["ligand_asym_id"] == ligand_asym_id,
        }
        improved_rmsd = (unpermuted_results["rmsd"] - best_results["rmsd"]).item()
        if improved_rmsd >= 1e-12:
            # better
            per_sample_log_dict.update(
                {
                    "algo:equivalent_permute": False,
                    "algo:worse_permute": False,
                    "algo:better_permute": True,
                    "algo:better_rmsd": improved_rmsd,
                }
            )
        elif improved_rmsd < 0:
            # worse
            per_sample_log_dict.update(
                {
                    "algo:equivalent_permute": False,
                    "algo:worse_permute": True,
                    "algo:better_permute": False,
                    "algo:worse_rmsd": -improved_rmsd,
                }
            )
        elif per_sample_log_dict["is_permuted"]:
            # equivalent
            per_sample_log_dict.update(
                {
                    "algo:equivalent_permute": True,
                    "algo:worse_permute": False,
                    "algo:better_permute": False,
                }
            )
        # atom indices to permute coordinates
        N_atom = aligned_pred_coord.size(-2)
        device = aligned_pred_coord.device
        atom_indices = torch.arange(N_atom, device=device)

        permute_asym_pair = [
            (best_results["pocket_asym_id"], pocket_asym_id),
            (best_results["ligand_asym_id"], ligand_asym_id),
        ]
        for asym_new, asym_old in permute_asym_pair:
            if asym_new == asym_old:
                continue
            # switch two chains
            ori_indices = atom_indices[atom_asym_id == asym_old]
            new_indices = atom_indices[atom_asym_id == asym_new]
            atom_indices[ori_indices.tolist()] = new_indices.clone()
            atom_indices[new_indices.tolist()] = ori_indices.clone()

        aligned_pred_coord = best_results.pop("aligned_pred_coord")[atom_indices, :]
        per_sample_log_dict["rmsd"] = best_results["rmsd"].item()

        return atom_indices, aligned_pred_coord, per_sample_log_dict

    N_sample = pred_coord.size(0)
    permute_pred_indices = []
    permuted_aligned_pred_coord = []
    sample_log_dicts = []
    for i in range(N_sample):
        atom_indices, aligned_pred_coord, per_sample_log_dict = (
            _find_protein_ligand_chains_for_one_sample(pred_coord[i])
        )
        permute_pred_indices.append(atom_indices)
        permuted_aligned_pred_coord.append(aligned_pred_coord)
        sample_log_dicts.append(per_sample_log_dict)

    permuted_aligned_pred_coord = torch.stack(permuted_aligned_pred_coord, dim=0)

    log_dict.update(
        traverse_and_aggregate(
            sample_log_dicts, aggregation_func=lambda x_list: sum(x_list) / N_sample
        )
    )

    # rmsd variance
    all_sample_rmsd = torch.tensor([x["rmsd"] for x in sample_log_dicts]).float()
    log_dict.update(
        {
            "rmsd_sample_std": all_sample_rmsd.std().item(),
            "rmsd_sample_gap": (all_sample_rmsd.max() - all_sample_rmsd.min()).item(),
        }
    )

    return permute_pred_indices, permuted_aligned_pred_coord, log_dict
