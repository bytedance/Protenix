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

import os

import torch

from protenix.utils.permutation import atom_permutation, chain_permutation


class SymmetricPermutation(object):
    """
    A symmetric permutation class for chain and atom permutations.

    Attributes:
        configs: Configuration settings for the permutation process.
        error_dir (str, optional): Directory to save error data. Defaults to None.
    """

    def __init__(self, configs, error_dir: str = None):
        self.configs = configs
        if error_dir is not None:
            self.chain_error_dir = os.path.join(error_dir, "chain_permutation")
            self.atom_error_dir = os.path.join(error_dir, "atom_permutation")
        else:
            self.chain_error_dir = None
            self.atom_error_dir = None

    def permute_label_to_match_mini_rollout(
        self,
        mini_coord: torch.Tensor,
        input_feature_dict: dict,
        label_dict: dict,
        label_full_dict: dict,
    ):
        """
        Apply permutation to label structure to match the predicted structure.
        This is mainly used to align label structure to the mini-rollout structure during training.

        Args:
            mini_coord (torch.Tensor): Coordinates of the predicted mini-rollout structure.
            input_feature_dict (dict): Input feature dictionary.
            label_dict (dict): Label dictionary.
            label_full_dict (dict): Full label dictionary.
        """

        assert mini_coord.dim() == 3

        log_dict = {}
        # 1. ChainPermutation: permute ground-truth chains to match mini-rollout prediction
        permuted_label_dict, chain_perm_log_dict, _, _ = chain_permutation.run(
            mini_coord[0],  # Only accepts a single structure
            input_feature_dict,
            label_full_dict,
            permute_label=True,
            error_dir=self.chain_error_dir,
            **self.configs.chain_permutation.configs,
        )
        if self.configs.chain_permutation.train.mini_rollout:
            label_dict.update(permuted_label_dict)
            log_dict.update(
                {
                    f"minirollout_perm/Chain-{k}": v
                    for k, v in chain_perm_log_dict.items()
                }
            )
        else:
            # Log only, not update the label_dict
            log_dict.update(
                {
                    f"minirollout_perm/Chain.F-{k}": v
                    for k, v in chain_perm_log_dict.items()
                }
            )

        # 2. AtomPermutation: permute ground-truth atoms to match mini-rollout prediction
        permuted_label_dict, atom_perm_log_dict, _ = atom_permutation.run(
            pred_coord=mini_coord[0],
            true_coord=label_dict["coordinate"],
            true_coord_mask=label_dict["coordinate_mask"],
            ref_space_uid=input_feature_dict["ref_space_uid"],
            atom_perm_list=input_feature_dict["atom_perm_list"],
            permute_label=True,
            error_dir=self.atom_error_dir,
            global_align_wo_symmetric_atom=self.configs.atom_permutation.global_align_wo_symmetric_atom,
        )

        if self.configs.atom_permutation.train.mini_rollout:
            label_dict.update(permuted_label_dict)
            log_dict.update(
                {f"minirollout_perm/Atom-{k}": v for k, v in atom_perm_log_dict.items()}
            )
        else:
            # Log only, not update the label_dict
            log_dict.update(
                {
                    f"minirollout_perm/Atom.F-{k}": v
                    for k, v in atom_perm_log_dict.items()
                }
            )

        return label_dict, log_dict

    def permute_diffusion_sample_to_match_label(
        self,
        input_feature_dict: dict,
        pred_dict: dict,
        label_dict: dict,
        stage: str,
        permute_by_pocket: bool = False,
    ):
        """
        Apply per-sample permutation to predicted structures to correct symmetries.
        Permutations are performed independently for each diffusion sample.

        Args:
            input_feature_dict (dict): Input feature dictionary.
            pred_dict (dict): Prediction dictionary.
            label_dict (dict): Label dictionary.
            stage (str): Current stage of the diffusion process, in ['train', 'test'].
            permute_by_pocket (bool): Whether to permute by pocket (for PoseBusters dataset). Defaults to False.
        """

        assert pred_dict["coordinate"].size(-2) == label_dict["coordinate"].size(
            -2
        ), "Cannot perform per-sample permutation on predicted structures if the label structure has more atoms."

        log_dict = {}
        permute_pred_indices, permute_label_indices = [], []
        if (
            stage != "train"
        ):  # During training stage, the label_dict is cropped after mini-rollout permutation.
            # In this case, chain permutation is not handled.

            # ChainPermutation: permute predicted chains to match label structure.

            (
                permuted_pred_dict,
                chain_perm_log_dict,
                permute_pred_indices,
                _,
            ) = chain_permutation.run(
                pred_dict["coordinate"],
                input_feature_dict,
                label_full_dict=label_dict,
                max_num_chains=-1,
                permute_label=False,
                permute_by_pocket=permute_by_pocket
                and self.configs.chain_permutation.permute_by_pocket,
                error_dir=self.chain_error_dir,
                **self.configs.chain_permutation.configs,
            )
            if self.configs.chain_permutation.get(stage).diffusion_sample:
                pred_dict.update(permuted_pred_dict)
                log_dict.update(
                    {
                        f"sample_perm/Chain-{k}": v
                        for k, v in chain_perm_log_dict.items()
                    }
                )
            else:
                # Log only, not update the pred_dict.
                log_dict.update(
                    {
                        f"sample_perm/Chain.F-{k}": v
                        for k, v in chain_perm_log_dict.items()
                    }
                )

        # AtomPermutation: permute predicted atoms to match label structure.
        # Permutations are performed independently for each diffusion sample.
        if permute_by_pocket and self.configs.atom_permutation.permute_by_pocket:
            if label_dict["pocket_mask"].dim() == 2:
                # the 0-the pocket is assumed to be the `main` pocket
                pocket_mask = label_dict["pocket_mask"][0]
                ligand_mask = label_dict["interested_ligand_mask"][0]
            else:
                pocket_mask = label_dict["pocket_mask"]
                ligand_mask = label_dict["interested_ligand_mask"]
            chain_mask = self.get_chain_mask_from_atom_mask(
                pocket_mask + ligand_mask,
                atom_to_token_idx=input_feature_dict["atom_to_token_idx"],
                token_asym_id=input_feature_dict["asym_id"],
            )
            alignment_mask = pocket_mask
        else:
            chain_mask = 1
            alignment_mask = None

        permuted_pred_dict, atom_perm_log_dict, atom_perm_pred_indices = (
            atom_permutation.run(
                pred_coord=pred_dict["coordinate"],
                true_coord=label_dict["coordinate"],
                true_coord_mask=label_dict["coordinate_mask"] * chain_mask,
                ref_space_uid=input_feature_dict["ref_space_uid"],
                atom_perm_list=input_feature_dict["atom_perm_list"],
                permute_label=False,
                alignment_mask=alignment_mask,
                error_dir=self.atom_error_dir,
                global_align_wo_symmetric_atom=self.configs.atom_permutation.global_align_wo_symmetric_atom,
            )
        )
        if permute_pred_indices:
            # Update `permute_pred_indices' according to the results of atom permutation
            updated_permute_pred_indices = []
            assert len(permute_pred_indices) == len(atom_perm_pred_indices)
            for chain_perm_indices, atom_perm_indices in zip(
                permute_pred_indices, atom_perm_pred_indices
            ):
                updated_permute_pred_indices.append(
                    chain_perm_indices[atom_perm_indices]
                )
            permute_pred_indices = updated_permute_pred_indices
        elif atom_perm_pred_indices is not None:
            permute_pred_indices = [
                atom_perm_indices for atom_perm_indices in atom_perm_pred_indices
            ]

        if self.configs.atom_permutation.get(stage).diffusion_sample:
            pred_dict.update(permuted_pred_dict)
            log_dict.update(
                {f"sample_perm/Atom-{k}": v for k, v in atom_perm_log_dict.items()}
            )
        else:
            # Log only, not update the pred_dict.
            log_dict.update(
                {f"sample_perm/Atom.F-{k}": v for k, v in atom_perm_log_dict.items()}
            )

        return pred_dict, log_dict, permute_pred_indices, permute_label_indices

    @staticmethod
    def get_chain_mask_from_atom_mask(
        atom_mask: torch.Tensor,
        atom_to_token_idx: torch.Tensor,
        token_asym_id: torch.Tensor,
    ):
        """
        Generate a chain mask from an atom mask.

        This method maps atoms to their corresponding token indices and then to their asym IDs. It then filters these asym IDs based on the atom mask and returns a mask indicating which atoms belong to the filtered chains.

        Args:
            atom_mask (torch.Tensor): A boolean atom mask. Shape: [N_atom].
            atom_to_token_idx (torch.Tensor): A tensor mapping each atom to its corresponding token index. Shape: [N_atom].
            token_asym_id (torch.Tensor): A tensor containing the asym ID for each token. Shape: [N_token].

        Returns:
            torch.Tensor: Chain mask. Shape: [N_atom].

        """

        atom_asym_id = token_asym_id[atom_to_token_idx.long()].long()
        assert atom_asym_id.size(0) == atom_mask.size(0)
        masked_asym_id = torch.unique(atom_asym_id[atom_mask.bool()])
        return torch.isin(atom_asym_id, masked_asym_id)

    @staticmethod
    def get_asym_id_match(
        permute_indices: torch.Tensor,
        atom_to_token_idx: torch.Tensor,
        token_asym_id: torch.Tensor,
    ) -> dict[int, int]:
        """Function to match asym IDs between original and permuted structure.

        Args:
            permute_indices (torch.Tensor): indices that specify the permuted ordering of atoms.
                [N_atom]
            atom_to_token_idx (torch.Tensor):  each entry maps an atom to its corresponding token index.
                [N_atom]
            token_asym_id (torch.Tensor): contains the asym ID for each token.
                [N_token]
        Returns:
            asym_id_match (Dict[int])
                A dictionary where the key is the original asym ID and the value is the permuted asym ID.
        """
        token_asym_id = token_asym_id.long()
        atom_to_token_idx = atom_to_token_idx.long()

        # Get the asym IDs for the original atoms
        original_atom_asym_id = token_asym_id[atom_to_token_idx]

        # Permute these IDs using the provided indices
        permuted_atom_asym_id = original_atom_asym_id[permute_indices]
        unique_asym_ids = torch.unique(original_atom_asym_id)

        asym_id_match = {}
        for ori_aid in unique_asym_ids:
            ori_aid = ori_aid.item()
            asym_mask = original_atom_asym_id == ori_aid
            perm_aid = permuted_atom_asym_id[asym_mask]

            assert (
                len(torch.unique(perm_aid)) == 1
            ), "Permuted asym ID must be unique for each original ID."

            asym_id_match[ori_aid] = perm_aid[0].item()

        return asym_id_match

    @staticmethod
    def permute_summary_confidence(
        summary_confidence_list: list[dict],
        permute_pred_indices: list[torch.Tensor],  # [N_atom]
        atom_to_token_idx: torch.Tensor,  # [N_atom]
        token_asym_id: torch.Tensor,  # [N_token]
        chain_keys: list[str] = ["chain_ptm", "chain_iptm", "chain_plddt"],
        chain_pair_keys: list[str] = [
            "chain_pair_iptm",
            "chain_pair_iptm_global",
            "chain_pair_plddt",
        ],
    ):
        """
        Permute summary confidence based on predicted indices.

        Args:
            summary_confidence_list (list[dict]): List of summary confidence dictionaries.
            permute_pred_indices (list[torch.Tensor]): List of predicted indices for permutation.
            atom_to_token_idx (torch.Tensor): Mapping from atoms to token indices.
            token_asym_id (torch.Tensor): Asym ID for each token.
            chain_keys (list[str], optional): Keys for chain-level confidence metrics. Defaults to ["chain_ptm", "chain_iptm", "chain_plddt"].
            chain_pair_keys (list[str], optional): Keys for chain pair-level confidence metrics. Defaults to ["chain_pair_iptm", "chain_pair_iptm_global", "chain_pair_plddt"].
        """

        assert len(summary_confidence_list) == len(permute_pred_indices)

        def _permute_one_sample(summary_confidence, permute_indices):
            # asym_id_match : {ori_asym_id: permuted_asym_id}
            asym_id_match = SymmetricPermutation.get_asym_id_match(
                permute_indices=permute_indices,
                atom_to_token_idx=atom_to_token_idx,
                token_asym_id=token_asym_id,
            )
            id_indices = torch.arange(len(asym_id_match), device=permute_indices.device)
            for i, j in asym_id_match.items():
                id_indices[j] = i

            # fix chain_id (asym_id) in summary_confidence
            for key in chain_keys:
                assert summary_confidence[key].dim() == 1
                summary_confidence[key] = summary_confidence[key][id_indices]
            for key in chain_pair_keys:
                assert summary_confidence[key].dim() == 2
                summary_confidence[key] = summary_confidence[key][:, id_indices]
                summary_confidence[key] = summary_confidence[key][id_indices, :]
            return summary_confidence, asym_id_match

        asym_id_match_list = []
        permuted_summary_confidence_list = []
        for i, (summary_confidence, perm_indices) in enumerate(
            zip(summary_confidence_list, permute_pred_indices)
        ):
            summary_confidence, asym_id_match = _permute_one_sample(
                summary_confidence, perm_indices
            )
            permuted_summary_confidence_list.append(summary_confidence)
            asym_id_match_list.append(asym_id_match)

        return permuted_summary_confidence_list, asym_id_match_list

    def permute_heads(
        self,
        pred_dict: dict,
        permute_pred_indices: list,
        atom_to_token_idx: torch.Tensor,
        rep_atom_mask: torch.Tensor,
    ):
        """
        Permute heads based on predicted indices.


        Args:
            pred_dict (dict): A dictionary containing the predicted components.
            permute_pred_indices (list): A list of tensors, each containing the predicted indices for the permutation of a diffusion sample.
            atom_to_token_idx (torch.Tensor): A tensor mapping each atom to its corresponding token index. Shape: [N_atom].
            rep_atom_mask (torch.Tensor): A boolean mask indicating which atoms are representative. Shape: [N_atom].

        Returns:
            dict: The updated `pred_dict`
        """

        for i, perm_indices in enumerate(permute_pred_indices):
            # permute atoms at dim=-2
            for key in ["plddt", "resolved"]:
                if key in pred_dict:
                    assert pred_dict[key].size(-2) == len(perm_indices)
                    pred_dict[key][..., i, :, :] = pred_dict[key][
                        ..., i, perm_indices, :
                    ]

            # permute tokens at dim=-2 and -3
            perm_atom_to_token_idx = atom_to_token_idx[perm_indices]
            perm_rep_atom_mask = rep_atom_mask[perm_indices]
            perm_token_indices = perm_atom_to_token_idx[perm_rep_atom_mask]
            for key in ["pae", "pde"]:
                if key in pred_dict:
                    assert (
                        pred_dict[key].size(-2)
                        == pred_dict[key].size(-3)
                        == len(perm_token_indices)
                    )
                    pred_dict[key] = pred_dict[key].to(perm_token_indices.device)
                    assert pred_dict[key].device == perm_token_indices.device
                    pred_dict[key][..., i, :, :, :] = pred_dict[key][
                        ..., i, perm_token_indices, :, :
                    ]
                    pred_dict[key][..., i, :, :, :] = pred_dict[key][
                        ..., i, :, perm_token_indices, :
                    ]

            # contact_probs
            if "contact_probs" in pred_dict:
                contact_probs_i = pred_dict["contact_probs"].clone()
                assert (
                    contact_probs_i.size(-1)
                    == contact_probs_i.size(-2)
                    == len(perm_token_indices)
                )
                contact_probs_i = contact_probs_i[..., perm_token_indices, :][
                    ..., perm_token_indices
                ]  # [N_token, N_token]
                pred_dict.setdefault("per_sample_contact_probs", []).append(
                    contact_probs_i
                )

        if "per_sample_contact_probs" in pred_dict:
            pred_dict["per_sample_contact_probs"] = torch.stack(
                pred_dict["per_sample_contact_probs"], dim=0
            )  # [N_sample, N_token, N_token]

        return pred_dict

    def permute_inference_pred_dict(
        self,
        input_feature_dict: dict,
        pred_dict: dict,
        label_dict: dict,
        permute_by_pocket: bool = False,
    ):
        """
        Permute predicted coordinates during inference.

        Args:
            input_feature_dict (dict): Input features dictionary.
            pred_dict (dict): Predicted dictionary.
            label_dict (dict): Label dictionary.
            permute_by_pocket (bool, optional): Whether to permute by pocket. Defaults to False.
        """
        # 1. Permute predicted coordinates
        pred_dict, log_dict, permute_pred_indices, _ = (
            self.permute_diffusion_sample_to_match_label(
                input_feature_dict,
                pred_dict=pred_dict,
                label_dict=label_dict,
                stage="test",
                permute_by_pocket=permute_by_pocket,
            )
        )

        if permute_pred_indices:
            # 2. Permute confidence logits
            pred_dict = self.permute_heads(
                pred_dict,
                permute_pred_indices=permute_pred_indices,
                atom_to_token_idx=input_feature_dict["atom_to_token_idx"],
                rep_atom_mask=input_feature_dict["pae_rep_atom_mask"].bool(),
            )

        return pred_dict, log_dict
