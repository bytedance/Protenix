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

import random

import torch

from protenix.metrics.rmsd import rmsd, self_aligned_rmsd
from protenix.utils.logger import get_logger
from protenix.utils.permutation.chain_permutation.utils import (
    apply_transform,
    get_optimal_transform,
    num_unique_matches,
)
from protenix.utils.permutation.utils import Checker

logger = get_logger(__name__)

ExtraLabelKeys = [
    "pocket_mask",
    "interested_ligand_mask",
    "chain_1_mask",
    "chain_2_mask",
    "entity_mol_id",
    "mol_id",
    "mol_atom_index",
    "pae_rep_atom_mask",
]


def correct_symmetric_chains(
    pred_dict: dict,
    label_full_dict: dict,
    extra_label_keys: list[str] = ExtraLabelKeys,
    max_num_chains: int = 20,
    permute_label: bool = True,
    **kwargs,
):
    """Inputs

    Args:
        pred_dict (dict[str, torch.Tensor]): A dictionary containing:
            - coordinate: pred_dict["coordinate"]
                shape = [N_cropped_atom, 3] or [Batch, N_cropped_atom, 3].
            - other keys: entity_mol_id, mol_id, mol_atom_index, pae_rep_atom_mask, is_ligand.
                shape = [N_cropped_atom]
        label_full_dict (dict[str, torch.Tensor]): A dictionary containing
            - coordinate: label_full_dict["coordinate"] and label_full_dict["coordinate_mask"]
                shape = [N_atom, 3] and [N_atom] (for coordinate_mask)
            - other keys: entity_mol_id, mol_id, mol_atom_index, pae_rep_atom_mask.
                shape = [N_atom]
            - extra keys: keys specified by extra_feature_keys.
        extra_label_keys (list[str]):
            - Additional features in label_full_dict that should be returned along with the permuted coordinates.
        max_num_chains (int): if the number of chains is more than this number, than skip permutation to
            avoid expensive computations.
        permute_label (bool): if true, permute the groundtruth chains, otherwise premute the prediction chains

    Return:
        output_dict:
            If permute_label=True, this is a dictionary containing
            - coordinate
            - coordinate_mask
            - features specified by extra_label_keys.
            If permute_label=False, this is a dictionary containing
            - coordinate.

        log_dict: statistics.

        permute_pred_indices / permute_label_indices:
            If batch_mode, this is a list of LongTensor. Otherwise, this is a LongTensor.
            The LongTensor gives the indices to permute either prediction or label.
    """

    assert pred_dict["coordinate"].dim() in [2, 3]
    batch_mode = pred_dict["coordinate"].dim() > 2

    if not batch_mode:
        (
            best_match,
            permute_pred_indices,
            permute_label_indices,
            output_dict,
            log_dict,
        ) = _correct_symmetric_chains_for_one_sample(
            pred_dict,
            label_full_dict,
            max_num_chains,
            permute_label,
            extra_label_keys=extra_label_keys,
            **kwargs,
        )
        return output_dict, log_dict, permute_pred_indices, permute_label_indices
    else:
        assert not permute_label, "Only supports prediction permutations in batch mode."
        pred_coord = []
        log_dict = {}
        best_matches = []
        permute_pred_indices = []
        permute_label_indices = []
        # Loop over all samples to find best matches one by one
        for i, pred_coord_i in enumerate(pred_dict["coordinate"]):
            (
                best_match_i,
                permute_pred_indices_i,
                permute_label_indices_i,
                pred_dict_i,
                log_dict_i,
            ) = _correct_symmetric_chains_for_one_sample(
                {**pred_dict, "coordinate": pred_coord_i},
                label_full_dict,
                max_num_chains,
                permute_label=False,
                extra_label_keys=[],
                **kwargs,
            )

            best_matches.append(best_match_i)
            permute_pred_indices.append(permute_pred_indices_i)
            permute_label_indices.append(permute_label_indices_i)
            pred_coord.append(pred_dict_i["coordinate"])
            for key, value in log_dict_i.items():
                log_dict.setdefault(key, []).append(value)

        output_dict = {"coordinate": torch.stack(pred_coord, dim=0)}

        log_dict = {key: sum(value) / len(value) for key, value in log_dict.items()}
        log_dict["N_unique_perm"] = num_unique_matches(best_matches)

        return output_dict, log_dict, permute_pred_indices, permute_label_indices


def _correct_symmetric_chains_for_one_sample(
    pred_dict: dict,
    label_full_dict: dict,
    max_num_chains: int = 20,
    permute_label: bool = False,
    extra_label_keys: list[str] = [],
    **kwargs,
):
    """
    Correct symmetric chains for a single sample by permuting either the predicted or the ground truth coordinates.
    """

    if not permute_label:
        """
        Permutation will act on the predicted coordinate.
        In this case, predicted structures and true structure need to have
        the same number of atoms.
        """
        assert pred_dict["coordinate"].size(-2) == label_full_dict["coordinate"].size(
            -2
        )

    with torch.no_grad():
        # Do not compute gradient while optimizing the permutation
        (
            best_match,
            permute_pred_indices,
            permute_label_indices,
            log_dict,
        ) = MultiChainPermutation(**kwargs)(
            pred_dict=pred_dict,
            label_full_dict=label_full_dict,
            max_num_chains=max_num_chains,
        )

    if permute_label:
        # Permute groundtruth coord and coord mask with indices, along the first dimension.
        indices = permute_label_indices.tolist()
        output_dict = {
            "coordinate": label_full_dict["coordinate"][indices, :],
            "coordinate_mask": label_full_dict["coordinate_mask"][indices],
        }
        # Permute extra label features, along the last dimension.
        output_dict.update(
            {
                k: label_full_dict[k][..., indices]
                for k in extra_label_keys
                if k in label_full_dict
            }
        )

    else:
        # Permute the predicted coord with permuted_indices
        indices = permute_pred_indices.tolist()
        output_dict = {
            "coordinate": pred_dict["coordinate"][indices, :],
        }

    return (
        best_match,
        permute_pred_indices,
        permute_label_indices,
        output_dict,
        log_dict,
    )


class MultiChainPermutation(object):
    """Anchor-based heuristic method.
    Find the best match that maps predicted chains to chains in the true complex.
    Here the predicted chains could be cropped, which could be fewer and shorter than
    those in the true complex.
    """

    def __init__(
        self, use_center_rmsd, find_gt_anchor_first, accept_it_as_it_is, *args, **kwargs
    ):
        self.use_center_rmsd = use_center_rmsd
        self.find_gt_anchor_first = find_gt_anchor_first
        self.accept_it_as_it_is = accept_it_as_it_is

    @staticmethod
    def dict_of_interested_keys(
        input_dict: dict,
        keys: list = [
            "mol_id",
            "entity_mol_id",
            "mol_atom_index",
            "pae_rep_atom_mask",
            "coordinate",
            "coordinate_mask",
            "is_ligand",
        ],
    ):
        """
        Extract a subset of keys from the input dictionary from the list `keys`.
        """

        return {k: input_dict[k] for k in keys if k in input_dict}

    def process_input(
        self,
        pred_dict: dict[str, torch.Tensor],
        label_full_dict: dict[str, torch.Tensor],
        max_num_chains: int = 20,
    ):
        """Process the input dicts

        Args:
            pred_dict (dict[str, torch.Tensor]): A dictionary containing
                entity_mol_id, mol_id, mol_atom_index, pae_rep_atom_mask, coordinate, is_ligand.
                All Tensors have shape = [N_cropped_atom]
            label_full_dict (dict[str, torch.Tensor]): A dictionary containing
                entity_mol_id, mol_id, mol_atom_index, pae_rep_atom_mask, coordinate, coordinate_mask.
                All Tensors have shape = [N_atom]
            max_num_chains (int): if the number of chains is more than this number, than skip permutation to
                avoid expensive computations.
            permute_label (bool): if true, permute the groundtruth chains, otherwise premute the prediction chains
        """

        log_dict = {}

        for key in ["entity_mol_id", "mol_id", "mol_atom_index"]:
            pred_dict[key] = pred_dict[key].long()
            label_full_dict[key] = label_full_dict[key].long()

        # get original unpermuted match
        pred_mol_id = set(torch.unique(pred_dict["mol_id"]).tolist())
        label_mol_id = set(torch.unique(label_full_dict["mol_id"]).tolist())
        if pred_mol_id.intersection(label_mol_id) != pred_mol_id:
            # if the mol_id in predicted structure is not a subset of label structure,
            # assert they contain the same number of atoms.
            assert pred_dict["coordinate"].size(-2) == label_full_dict[
                "coordinate"
            ].size(-2)
            self.unpermuted_match = self.check_pattern_and_create_mapping(
                pred_dict["mol_id"], label_full_dict["mol_id"]
            )
        else:
            self.unpermuted_match = {
                i: i for i in torch.unique(pred_dict["mol_id"]).tolist()
            }

        if len(torch.unique(label_full_dict["entity_mol_id"])) == len(
            torch.unique(label_full_dict["mol_id"])
        ):
            # No permutation is needed
            has_sym_chain = False
            return self.unpermuted_match, has_sym_chain
        else:
            has_sym_chain = True

        n_label_chain = len(torch.unique(label_full_dict["mol_id"]))
        if n_label_chain > 20:
            logger.warning(f"The label_full_dict contains {n_label_chain} asym chains.")

        if max_num_chains > 0 and n_label_chain > max_num_chains:
            logger.warning(
                f"The label_full_dict contains {n_label_chain} asym chains (max_num_chains: {max_num_chains}). Will skip chain permutation and keep the original chain assignment."
            )
            return self.unpermuted_match, has_sym_chain

        # parse features to token-level
        self.label_token_dict, self.label_asym_dict = self._parse_atom_feature_dict(
            self.dict_of_interested_keys(label_full_dict),
            rep_atom_mask=label_full_dict["pae_rep_atom_mask"],
        )
        self.pred_token_dict, self.pred_asym_dict = self._parse_atom_feature_dict(
            self.dict_of_interested_keys(pred_dict),
            rep_atom_mask=pred_dict["pae_rep_atom_mask"],
        )

        # get mapping between entity_id and asym_id
        self.label_token_dict.update(
            self._get_entity_asym_mapping(
                self.label_token_dict["entity_mol_id"], self.label_token_dict["mol_id"]
            )
        )
        self.pred_token_dict.update(
            self._get_entity_asym_mapping(
                self.pred_token_dict["entity_mol_id"], self.pred_token_dict["mol_id"]
            )
        )
        return None, has_sym_chain

    @staticmethod
    def check_pattern_and_create_mapping(mol_id1: torch.Tensor, mol_id2: torch.Tensor):
        """
        Check if the patterns between two mol_id tensors match and create a mapping between them.

        Args:
            mol_id1 (torch.Tensor): A tensor of mol IDs from the first set.
            mol_id2 (torch.Tensor): A tensor of mol IDs from the second set.

        Returns:
            dict: A dictionary mapping mol IDs from mol_id1 to mol_id2.
        """
        if mol_id1.shape != mol_id2.shape:
            raise ValueError("mol_id1 and mol_id2 must have the same shape")

        pattern_mapping = {}
        for id1, id2 in zip(mol_id1.tolist(), mol_id2.tolist()):
            if id1 in pattern_mapping:
                if pattern_mapping[id1] != id2:
                    raise ValueError(
                        f"Inconsistent pattern: {id1} mapped to different values in mol_id2"
                    )
            else:
                if id2 in pattern_mapping.values():
                    raise ValueError(
                        f"Value {id2} in mol_id2 already mapped to another value"
                    )
                pattern_mapping[id1] = id2

        return pattern_mapping

    def _parse_atom_feature_dict(
        self, atom_features: dict, rep_atom_mask: torch.Tensor
    ):
        """
        Parse the atom feature dictionary and convert it to token features and per-asym token features.

        Args:
            atom_features (dict): A dictionary containing atom features.
            rep_atom_mask (torch.Tensor): The rep atom mask.

        Returns:
            tuple: A tuple containing:
                - token_dict (dict): A dictionary containing the token features corresponding to the rep atoms.
                - asym_token_dict (dict): A dictionary where keys are asym IDs and values are dictionaries of features corresponding to each asym ID.
        """

        # Atom features --> Token features
        token_dict = self._convert_to_token_dict(
            atom_dict=atom_features,
            rep_atom_mask=rep_atom_mask.bool(),
        )

        # Token features --> per asym token features
        asym_token_dict = self._convert_to_per_asym_feature_dict(
            asym_id=token_dict["mol_id"],
            feature_dict={
                "coordinate": token_dict.get("coordinate"),
                "coordinate_mask": token_dict.get("coordinate_mask", None),
                "mol_atom_index": token_dict.get("mol_atom_index"),
            },
        )

        return token_dict, asym_token_dict

    @staticmethod
    def _convert_to_token_dict(
        atom_dict: dict[str, torch.Tensor], rep_atom_mask: torch.Tensor
    ) -> dict[str, torch.Tensor]:
        """
        Convert the atom feature dictionary to a token feature dictionary based on the rep atom mask.

        Args:
            atom_dict (dict[str, torch.Tensor]): A dictionary containing atom features.
            rep_atom_mask (torch.Tensor): The rep atom mask.

        Returns:
            dict[str, torch.Tensor]: A dictionary containing the token features corresponding to the rep atoms.
        """

        rep_atom_mask = rep_atom_mask.bool()
        return {k: v[rep_atom_mask] for k, v in atom_dict.items() if v is not None}

    @staticmethod
    def _convert_to_per_asym_feature_dict(asym_id: torch.Tensor, feature_dict: dict):
        """
        Convert the feature dictionary to a dictionary where keys are asym IDs and values are dictionaries of features corresponding to each asym ID.

        Args:
            asym_id (torch.Tensor): A tensor of asym IDs.
            feature_dict (dict): A dictionary containing features for all atoms.

        Returns:
            dict: A dictionary where keys are asym IDs and values are dictionaries of features corresponding to each asym.
        """
        out = {}

        for aid in torch.unique(asym_id):
            mask = asym_id == aid
            out[aid.item()] = {
                k: v[mask] for k, v in feature_dict.items() if v is not None
            }
        return out

    @staticmethod
    def _get_entity_asym_mapping(
        entity_id: torch.Tensor, asym_id: torch.Tensor
    ) -> tuple[dict]:
        """
        Generate mappings between entity IDs and asym IDs.

        Args:
            entity_id (torch.Tensor): A tensor of entity IDs.
            asym_id (torch.Tensor): A tensor of asym IDs.

        Returns:
            tuple[dict]: A tuple containing two dictionaries:
                - entity_to_asym: A dictionary mapping entity IDs to their corresponding asym IDs.
                - asym_to_entity: A dictionary mapping asym IDs to their corresponding entity IDs.
        """

        entity_to_asym = {}
        asym_to_entity = {}
        for ein in torch.unique(entity_id):
            ein = ein.item()
            asyms = torch.unique(asym_id[entity_id == ein])
            entity_to_asym[ein] = asyms
            asym_to_entity.update({a.item(): ein for a in asyms})

        return {"entity_to_asym": entity_to_asym, "asym_to_entity": asym_to_entity}

    def find_anchor_asym_chain_in_predictions(self) -> tuple[int]:
        """
        Find anchor chains in the prediction.

        Ref: AlphaFold3 SI Chapter 4.2. -> AlphaFold Multimer Chapter 7.3.1
        In the alignment phase, we pick a pair of anchor asyms to align,
        one in the ground truth and one in the prediction.
        The ground truth anchor asym a_gt is chosen to be the least ambiguous possible,
        for example in an A3B2 complex an arbitrary B asym is chosen.
        In the event of a tie e.g. A2B2 stoichiometry, the longest asym is chosen,
        with the hope that in general the longer asyms are likely to have higher confident predictions.
        The prediction anchor asym is chosen from the set {a^pred_m} of all prediction asyms
        with the same sequence as the ground truth anchor asym.

        Return:
            anchor_pred_asym_id (int): selected asym chain.
        """

        # Do not consider asym with fewer than 4 tokens in Prediction
        asym_to_asym_length = {
            asym_id: len(asym_dict["coordinate"])
            for asym_id, asym_dict in self.pred_asym_dict.items()
        }
        valid_asyms = [asym_id for asym_id, l in asym_to_asym_length.items() if l >= 4]

        # Do not consider entities with fewer than 4 resolved tokens in GT
        valid_entities = []
        for ent, asyms in self.label_token_dict["entity_to_asym"].items():
            if any(
                self.label_asym_dict[asym.item()]["coordinate_mask"].sum().item() >= 4
                for asym in asyms
            ):
                valid_entities.append(ent)

        valid_entity_asym = [
            (ent, asym.item())
            for ent in valid_entities
            for asym in self.pred_token_dict["entity_to_asym"][ent]
            if asym.item() in valid_asyms
        ]

        candidate_entities = set(ent for ent, _ in valid_entity_asym)

        # Find polymer chains in the prediction
        pred_polymer_entity_id = []
        for ent_id in candidate_entities:
            mask = self.pred_token_dict["entity_mol_id"] == ent_id
            is_ligand = self.pred_token_dict["is_ligand"][mask]
            if (
                torch.sum(is_ligand) <= is_ligand.shape[0] / 2
                and is_ligand.shape[0]
                >= 12  # do not prioritize asym with too few tokens
            ):
                pred_polymer_entity_id.append(ent_id)

        # Prioritize polymer
        if len(pred_polymer_entity_id) > 0:
            candidate_entities = pred_polymer_entity_id

        # Choose entities with fewest asyms in GT
        entity_to_asym_count = {
            k: len(self.label_token_dict["entity_to_asym"][k])
            for k in candidate_entities
        }
        min_asym_count = min(list(entity_to_asym_count.values()))
        candidate_entities = [
            ent
            for ent, count in entity_to_asym_count.items()
            if count == min_asym_count
        ]

        # Choose longest asyms in Prediction
        candidate_asyms = [
            asym_id for ent, asym_id in valid_entity_asym if ent in candidate_entities
        ]
        max_asym_length = max(
            asym_to_asym_length[asym_id] for asym_id in candidate_asyms
        )
        candidate_asyms = [
            asym_id
            for asym_id in candidate_asyms
            if asym_to_asym_length[asym_id] == max_asym_length
        ]

        # If multiple asym chains remain, return a random one.
        anchor_pred_asym_id = random.choice(candidate_asyms)

        return anchor_pred_asym_id

    @staticmethod
    def _select_atoms_by_mol_atom_index(input_dict: dict, mol_atom_index: torch.Tensor):
        """
        Select atoms from the input dictionary based on the specified mol_atom_index.

        Args:
            input_dict (dict): Input dict.
            mol_atom_index (torch.Tensor): A tensor of atom indices.

        Returns:
            dict: A dictionary containing the selected atom features.
        """
        mask = torch.isin(input_dict["mol_atom_index"], mol_atom_index)
        out_dict = {k: v[mask] for k, v in input_dict.items()}
        assert (out_dict["mol_atom_index"] == mol_atom_index).all()
        return out_dict

    def compute_best_match_heuristic(self):
        """
        Compute the best chain permutation between prediction and groundtruth.


        Returns:
            dict[int, int]: A dictionary mapping pred chain IDs to those of the groundtruth.
        """

        # Find anchor asym chain in predictions
        anchor_pred_asym_id = self.find_anchor_asym_chain_in_predictions()
        anchor_entity_id = self.pred_token_dict["asym_to_entity"][anchor_pred_asym_id]

        if self.find_gt_anchor_first:
            # Randomly sample a groundtruth asym chain using this entity id
            anchor_gt_asym_id = self.label_token_dict["entity_to_asym"][
                anchor_entity_id
            ].tolist()
            anchor_gt_asym_id = random.choice(anchor_gt_asym_id)

            # The candidate anchors to be matched are from prediction
            candidate_anchors = self.pred_token_dict["entity_to_asym"][anchor_entity_id]
        else:

            # The candidate anchors to be matched are from groundtruth
            candidate_anchors = self.label_token_dict["entity_to_asym"][
                anchor_entity_id
            ]

        # Find best match
        best_rmsd = torch.inf
        best_match = None

        for anchor_k in candidate_anchors:
            anchor_k = anchor_k.item()

            if self.find_gt_anchor_first:
                gt_anchor, pred_anchor = anchor_gt_asym_id, anchor_k
            else:
                gt_anchor, pred_anchor = anchor_k, anchor_pred_asym_id

            # Find atoms in GT chain to match atoms in predicted chain (which could be cropped)
            gt_anchor_dict = MultiChainPermutation._select_atoms_by_mol_atom_index(
                self.label_asym_dict[gt_anchor],
                mol_atom_index=self.pred_asym_dict[pred_anchor]["mol_atom_index"],
            )

            # Align GT Anchor to Pred Anchor
            mask = gt_anchor_dict["coordinate_mask"].bool()  # use GT coordinate_mask
            if not mask.any():
                continue
            rot, trans = get_optimal_transform(
                gt_anchor_dict["coordinate"][mask],
                self.pred_asym_dict[pred_anchor]["coordinate"][mask],
            )

            # Transform all GT coordinates according to the aligment results
            aligned_coordinate = apply_transform(
                self.label_token_dict["coordinate"], rot, trans
            )
            for asym_id in self.label_asym_dict:
                self.label_asym_dict[asym_id]["aligned_coordinate"] = (
                    aligned_coordinate[self.label_token_dict["mol_id"] == asym_id]
                )

            # Greedily matches all remaining chains
            matched_asym = {pred_anchor: gt_anchor}
            to_be_matched = [k for k in self.pred_asym_dict if k != pred_anchor]
            candidate_gt_asym_id = [k for k in self.label_asym_dict if k != gt_anchor]

            # Sort the remaining chains by their length, so that longer chain chooses its match first.
            to_be_matched = sorted(
                to_be_matched,
                key=lambda k: -self.pred_asym_dict[k]["coordinate"].size(-2),
            )

            while len(to_be_matched) > 0:
                cur_pred_asym_id = to_be_matched.pop(0)
                cur_entity_id = self.pred_token_dict["asym_to_entity"][cur_pred_asym_id]
                cur_gt_asym_ids = self.label_token_dict["entity_to_asym"][
                    cur_entity_id
                ].tolist()
                matched_gt_asym_id, _ = self.match_pred_asym_to_gt_asym(
                    cur_pred_asym_id,
                    [asym for asym in cur_gt_asym_ids if asym in candidate_gt_asym_id],
                )
                matched_asym[cur_pred_asym_id] = matched_gt_asym_id
                candidate_gt_asym_id.remove(matched_gt_asym_id)

            assert len(matched_asym) == len(self.pred_asym_dict)

            # Calculate RMSD
            total_rmsd = self.calculate_rmsd(matched_asym)

            if total_rmsd < best_rmsd:
                best_rmsd = total_rmsd
                best_match = matched_asym

        assert best_match is not None

        return best_match

    def calculate_rmsd(self, asym_match: dict):
        """
        Calculate the RMSD given a match.
        """

        return sum(self._calculate_rmsd(a, b) for a, b in asym_match.items()) / len(
            asym_match
        )

    def _calculate_rmsd(self, pred_asym_id: int, gt_asym_id: int):
        """
        Calculate the RMSD between the predicted and ground truth chains, either using the average of the representative atoms or all of them.

        Args:
            pred_asym_id (int): The ID of the predicted asymmetric chain.
            gt_asym_id (int): The ID of the ground truth asymmetric chain.

        Returns:
            float: The calculated RMSD.
        """

        pred_asym_dict = self.pred_asym_dict[pred_asym_id]
        label_asym_dict = MultiChainPermutation._select_atoms_by_mol_atom_index(
            self.label_asym_dict[gt_asym_id], pred_asym_dict["mol_atom_index"]
        )
        mask = label_asym_dict["coordinate_mask"].bool()
        if not mask.any():
            return 0.0
        elif self.use_center_rmsd:
            return rmsd(
                pred_asym_dict["coordinate"][mask].mean(dim=-2, keepdim=True),
                label_asym_dict["aligned_coordinate"][mask].mean(dim=-2, keepdim=True),
            ).item()
        else:
            return rmsd(
                pred_asym_dict["coordinate"][mask],
                label_asym_dict["aligned_coordinate"][mask],
            ).item()

    def match_pred_asym_to_gt_asym(self, pred_asym_id: int, gt_asym_ids: list):
        """
        Match a predicted  chain to the groundtruth chain based on the average of the representative atoms.

        Args:
            pred_asym_id (int): The ID of the predicted asymmetric chain.
            gt_asym_ids (list[int]): A list or tensor of ground truth asymmetric chain IDs.

        Returns:
            tuple: A tuple containing:
                - best_gt_asym_id (int): The ID of the best matched ground truth asymmetric chain.
                - best_error (float): The distance error between the centers of mass of the best matched chains.
        """

        pred_asym_dict = self.pred_asym_dict[pred_asym_id]

        best_error = torch.inf
        best_gt_asym_id = None
        unresolved_gt_asym_id = []
        for gt_asym_id in gt_asym_ids:
            if isinstance(gt_asym_id, torch.Tensor):
                gt_asym_id = gt_asym_id.item()

            # Select cropped atoms by comparing to mol_atom_index in prediction
            label_asym_dict = MultiChainPermutation._select_atoms_by_mol_atom_index(
                self.label_asym_dict[gt_asym_id], pred_asym_dict["mol_atom_index"]
            )
            mask = label_asym_dict["coordinate_mask"].bool()

            if not mask.any():
                # Skip unresolved ones
                unresolved_gt_asym_id.append(gt_asym_id)
                continue

            gt_center = label_asym_dict["aligned_coordinate"][mask].mean(dim=0)
            pred_center = pred_asym_dict["coordinate"][mask].mean(dim=0)

            delta = torch.norm(gt_center - pred_center)

            if delta < best_error:
                best_error = delta
                best_gt_asym_id = gt_asym_id

        if best_gt_asym_id is None:
            # If only unresolved ones remains, return the first one
            assert len(unresolved_gt_asym_id) > 0
            best_gt_asym_id, best_error = gt_asym_ids[0], 0

        return best_gt_asym_id, best_error

    @staticmethod
    def build_permuted_indice(
        pred_dict: dict, label_full_dict: dict, best_match: dict[int, int]
    ):
        """
        Build permutation indices from the pred-gt chain mapping.
        Args:
            pred_dict (dict): A dictionary containing the predicted coordinates.
            label_full_dict (dict): A dictionary containing the true coordinates and their masks.
            best_match (dict[int, int]): {pred_mol_id: gt_mol_id} best match between pred asym chains and gt asym chains

        Returns:
            indices (torch.Tensor): Permutation indices.
        """

        # Get the number of predicted (cropped) atoms
        N_pred_atom = pred_dict["mol_id"].size(0)
        N_label_atom = label_full_dict["mol_id"].size(0)
        indices = pred_dict["mol_id"].new_zeros(size=(N_pred_atom,))
        full_indices = torch.arange(N_label_atom, device=indices.device)

        for pred_asym_id, gt_asym_id in best_match.items():
            # Create a mask for the predicted asym_id
            mask = pred_dict["mol_id"] == pred_asym_id
            mol_atom_index = pred_dict["mol_atom_index"][mask]

            # Creat a mask for the matched gt asym_id
            gt_mask = label_full_dict["mol_id"] == gt_asym_id
            # Extract indices according to 'mol_atom_index'
            gt_asym_dict = MultiChainPermutation._select_atoms_by_mol_atom_index(
                {
                    "mol_atom_index": label_full_dict["mol_atom_index"][gt_mask],
                    "indices": full_indices[gt_mask],
                },
                mol_atom_index,
            )
            indices[mask] = gt_asym_dict["indices"].clone()

        assert len(torch.unique(indices)) == len(indices)
        return indices

    @staticmethod
    def aligned_rmsd(
        pred_dict: dict,
        label_full_dict: dict,
        indices: torch.Tensor,
        reduce: bool = True,
        eps: float = 1e-8,
    ):
        """
        Calculate the global aligned RMSD between predicted and true coordinates.

        Args:
            pred_dict (dict): A dictionary containing the predicted coordinates.
            label_full_dict (dict): A dictionary containing the true coordinates and their masks.
            indices (torch.Tensor): Indices to select from the true coordinates.
            reduce (bool): If True, reduce the RMSD over the batch dimension.
            eps (float): A small value to avoid division by zero.

        Returns:
            float: The aligned RMSD value.
        """

        with torch.cuda.amp.autocast(enabled=False):
            aligned_rmsd, _, _, _ = self_aligned_rmsd(
                pred_pose=pred_dict["coordinate"].to(torch.float32),
                true_pose=label_full_dict["coordinate"][indices, :].to(torch.float32),
                atom_mask=label_full_dict["coordinate_mask"][indices],
                allowing_reflection=False,
                reduce=reduce,
                eps=eps,
            )
        return aligned_rmsd.item()

    def __call__(
        self,
        pred_dict: dict[str, torch.Tensor],
        label_full_dict: dict[str, torch.Tensor],
        max_num_chains: int = 20,
    ):
        """
        Call function for the class

        Args:
            pred_dict (dict): A dictionary containing the predicted coordinates.
            label_full_dict (dict): A dictionary containing the groundtruth and its attributes.
            max_num_chains (int): Maximum number of chains allowed.

        Returns:
            tuple: A tuple containing:
                - best_match (dict[int, int]): The best match between predicted and groundtruth chains.
                - permute_pred_indices (torch.Tensor or None): Indices to permute the predicted coordinates.
                - permuted_indices (torch.Tensor): Indices to permute the groundtruth coordinates.
                - log_dict (dict): A dictionary detailing the permutation information.
        """
        match, has_sym_chain = self.process_input(
            pred_dict, label_full_dict, max_num_chains
        )

        if match is not None:
            """
            Either the structure does not contain symmetric chains, or
            there are too many chains so that the algorithm gives up.
            """
            indices = self.build_permuted_indice(pred_dict, label_full_dict, match)
            pred_indices = torch.argsort(indices)
            return match, pred_indices, indices, {"has_sym_chain": False}

        # Core step: get best mol_id match

        best_match = self.compute_best_match_heuristic()

        permuted_indices = self.build_permuted_indice(
            pred_dict, label_full_dict, best_match
        )

        log_dict = {
            "has_sym_chain": True,
            "is_permuted": num_unique_matches([best_match, self.unpermuted_match]) > 1,
            "algo:no_permute": num_unique_matches([best_match, self.unpermuted_match])
            == 1,
        }

        if log_dict["algo:no_permute"]:
            # No permutation, return now
            pred_indices = torch.argsort(permuted_indices)
            return best_match, pred_indices, permuted_indices, log_dict

        # Compare rmsd before/after permutation
        unpermuted_indices = self.build_permuted_indice(
            pred_dict, label_full_dict, self.unpermuted_match
        )

        permuted_rmsd = self.aligned_rmsd(pred_dict, label_full_dict, permuted_indices)
        unpermuted_rmsd = self.aligned_rmsd(
            pred_dict, label_full_dict, unpermuted_indices
        )
        improved_rmsd = unpermuted_rmsd - permuted_rmsd
        if improved_rmsd >= 1e-12:
            # Case with better permutation
            log_dict.update(
                {
                    "algo:equivalent_permute": False,
                    "algo:worse_permute": False,
                    "algo:better_permute": True,
                    "algo:better_rmsd": improved_rmsd,
                }
            )
        elif improved_rmsd < 0:
            # Case with worse permutation
            log_dict.update(
                {
                    "algo:equivalent_permute": False,
                    "algo:worse_permute": True,
                    "algo:better_permute": False,
                    "algo:worse_rmsd": -improved_rmsd,
                }
            )
        elif not log_dict["algo:no_permute"]:
            # Case with equivalent permutation
            log_dict.update(
                {
                    "algo:equivalent_permute": True,
                    "algo:worse_permute": False,
                    "algo:better_permute": False,
                }
            )
        else:
            # No permutation
            log_dict["debug:zero_rmsd"] = improved_rmsd

        # Revert worse/equivalent permute to original chain assignment
        if (not self.accept_it_as_it_is) and (
            log_dict["algo:equivalent_permute"] or log_dict["algo:worse_permute"]
        ):
            # Revert to original chain assignment
            best_match = self.unpermuted_match
            permuted_indices = unpermuted_indices
            log_dict["is_permuted"] = False

        if pred_dict["coordinate"].size(-2) == label_full_dict["coordinate"].size(-2):
            Checker.is_permutation(permuted_indices)  # indices to permute/crop label
            permute_pred_indices = torch.argsort(
                permuted_indices
            )  # Indices to permute pred
        else:
            # Hard to `define` permute_pred_indices in this case
            permute_pred_indices = None

        return best_match, permute_pred_indices, permuted_indices, log_dict
