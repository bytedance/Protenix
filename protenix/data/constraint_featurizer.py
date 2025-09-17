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

import copy
import hashlib
from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from biotite.structure import AtomArray
from scipy.spatial.distance import cdist

from protenix.data.constants import ELEMS, STD_RESIDUES
from protenix.data.tokenizer import Token, TokenArray
from protenix.data.utils import get_atom_mask_by_name
from protenix.utils.logger import get_logger

logger = get_logger(__name__)


class ConstraintFeatureGenerator:
    """Generates constraint features by coordinating different featurizers"""

    def __init__(
        self, constraint_config: dict[str, Any], ab_top2_clusters: set[str] = None
    ):
        """
        Args:
            constraint_config: Configuration dict for constraints
            ab_top2_clusters: Optional antibody cluster info
        """
        self.constraint = constraint_config
        self.ab_top2_clusters = ab_top2_clusters

    # function set for inference
    @staticmethod
    def _canonicalize_contact_format(
        sequences: list[dict[str, Any]], pair: dict[str, Any]
    ) -> dict[str, Any]:
        pair = copy.deepcopy(pair)
        _pair = {}

        for id_num in ["1", "2"]:
            res_info = []
            for key in ["entity", "copy", "position", "atom"]:
                identifier_value = pair.get(f"{key}{id_num}", None)
                if identifier_value is None:
                    assert (
                        key == "atom"
                    ), "contact should have at least 3 identifiers('entity', 'copy', 'position')"
                if key == "atom" and (identifier_value is not None):
                    if isinstance(identifier_value, int):
                        entity_dict = list(
                            sequences[
                                int(pair.get(f"entity{id_num}", None) - 1)
                            ].values()
                        )[0]
                        assert "atom_map_to_atom_name" in entity_dict
                        identifier_value = entity_dict["atom_map_to_atom_name"][
                            identifier_value
                        ]
                res_info.append(identifier_value)
            _pair[f"id{id_num}"] = res_info

        if hash(tuple(_pair["id1"][:2])) == hash(tuple(_pair["id2"][:2])):
            raise ValueError("A contact pair can not be specified on the same chain")

        _pair["min_distance"] = float(pair.get("min_distance", 0))
        _pair["max_distance"] = float(pair["max_distance"])
        if _pair["max_distance"] < _pair["min_distance"]:
            raise ValueError("max_distance must be greater than min_distance")
        if "atom1" in pair and "atom2" in pair:
            _pair["contact_type"] = "atom_contact"
        else:
            _pair["contact_type"] = "token_contact"
        return _pair

    @staticmethod
    def _canonicalize_pocket_res_format(binder: dict, pocket_pos: dict) -> dict:
        assert len(pocket_pos) == 3
        if hash(tuple([binder["entity"], binder["copy"]])) == hash(
            tuple([pocket_pos["entity"], pocket_pos["copy"]])
        ):
            raise ValueError("Pockets can not be the same chain with the binder")
        return pocket_pos

    @staticmethod
    def _log_constraint_feature(
        atom_array: AtomArray, token_array: TokenArray, constraint_feature: dict
    ):

        atom_feature = constraint_feature["contact_atom"]
        if atom_feature.sum() > 0:
            # logging contact feature
            token_idx_1, token_idx_2 = torch.nonzero(
                torch.triu(atom_feature[..., 1]), as_tuple=True
            )

            atom1_index = token_array[token_idx_1].get_annotation("centre_atom_index")
            atom2_index = token_array[token_idx_2].get_annotation("centre_atom_index")

            res_name_1 = atom_array.res_name[atom1_index]
            res_name_2 = atom_array.res_name[atom2_index]

            atom_name_1 = atom_array.atom_name[atom1_index]
            atom_name_2 = atom_array.atom_name[atom2_index]

            chain_id_1 = atom_array.chain_id[atom1_index]
            chain_id_2 = atom_array.chain_id[atom2_index]

            max_distance = atom_feature[token_idx_1, token_idx_2, 1]
            min_distance = atom_feature[token_idx_1, token_idx_2, 0]

            contact_info = {
                "chain_id": np.stack([chain_id_1, chain_id_2]).T.tolist(),
                "res_name": np.stack([res_name_1, res_name_2]).T.tolist(),
                "atom_name": np.stack([atom_name_1, atom_name_2]).T.tolist(),
                "max_distance": max_distance.tolist(),
                "min_distance": min_distance.tolist(),
            }
            logger.info(f"loaded atom contact info:{contact_info}")

        token_feature = constraint_feature["contact"]
        if token_feature.sum() > 0:
            # logging contact feature
            token_idx_1, token_idx_2 = torch.nonzero(
                torch.triu(token_feature[..., 1]), as_tuple=True
            )

            atom1_index = token_array[token_idx_1].get_annotation("centre_atom_index")
            atom2_index = token_array[token_idx_2].get_annotation("centre_atom_index")

            res_name_1 = atom_array.res_name[atom1_index]
            res_name_2 = atom_array.res_name[atom2_index]

            atom_name_1 = atom_array.atom_name[atom1_index]
            atom_name_2 = atom_array.atom_name[atom2_index]

            chain_id_1 = atom_array.chain_id[atom1_index]
            chain_id_2 = atom_array.chain_id[atom2_index]

            max_distance = token_feature[token_idx_1, token_idx_2, 1]
            min_distance = token_feature[token_idx_1, token_idx_2, 0]

            contact_info = {
                "chain_id": np.stack([chain_id_1, chain_id_2]).T.tolist(),
                "res_name": np.stack([res_name_1, res_name_2]).T.tolist(),
                "atom_name": np.stack([atom_name_1, atom_name_2]).T.tolist(),
                "max_distance": max_distance.tolist(),
                "min_distance": min_distance.tolist(),
            }
            logger.info(f"loaded contact info:{contact_info}")

        pocket_feature = constraint_feature["pocket"]
        if pocket_feature.sum() > 0:
            # logging contact feature
            binder_idx, pocket_idx = torch.nonzero(
                pocket_feature[..., 0], as_tuple=True
            )
            atom1_index = token_array[binder_idx].get_annotation("centre_atom_index")
            atom2_index = token_array[pocket_idx].get_annotation("centre_atom_index")

            res_name_1 = atom_array.res_name[atom1_index]
            res_name_2 = atom_array.res_name[atom2_index]

            atom_name_1 = atom_array.atom_name[atom1_index]
            atom_name_2 = atom_array.atom_name[atom2_index]

            chain_id_1 = atom_array.chain_id[atom1_index]
            chain_id_2 = atom_array.chain_id[atom2_index]

            distance = pocket_feature[binder_idx, pocket_idx, 0]

            pocket_info = {
                "binder_chain_id": np.unique(chain_id_1).tolist(),
                "binder_res_name": np.unique(res_name_1).tolist(),
                "binder_atom_name": np.unique(atom_name_1).tolist(),
                "pocket_chain_id": np.unique(chain_id_2).tolist(),
                "pocket_res_name": np.unique(res_name_2).tolist(),
                "pocket_atom_name": np.unique(atom_name_2).tolist(),
                "distance": distance.unique().item(),
            }
            logger.info(f"loaded pocket info:{pocket_info}")

    @staticmethod
    def generate_from_json(
        token_array: TokenArray,
        atom_array: AtomArray,
        sequences: list[dict[str, Any]],
        constraint_param: dict,
    ) -> tuple[dict[str, Any], TokenArray, AtomArray]:
        feature_dict = {}

        # build atom-level contact features
        contact_inputs = [
            ConstraintFeatureGenerator._canonicalize_contact_format(sequences, pair)
            for pair in constraint_param.get("contact", [])
        ]
        atom_contact_inputs = list(
            filter(lambda d: d["contact_type"] == "atom_contact", contact_inputs)
        )
        token_contact_inputs = list(
            filter(lambda d: d["contact_type"] == "token_contact", contact_inputs)
        )

        atom_to_token_idx_dict = {}
        for idx, token in enumerate(token_array.tokens):
            for atom_idx in token.atom_indices:
                atom_to_token_idx_dict[atom_idx] = idx

        atom_contact_specifics = []
        for i, pair in enumerate(atom_contact_inputs):
            atom_mask1 = get_atom_mask_by_name(
                atom_array=atom_array,
                entity_id=pair["id1"][0],
                copy_id=pair["id1"][1],
                position=pair["id1"][2],
                atom_name=pair["id1"][3],
            )
            atom_mask2 = get_atom_mask_by_name(
                atom_array=atom_array,
                entity_id=pair["id2"][0],
                copy_id=pair["id2"][1],
                position=pair["id2"][2],
                atom_name=pair["id2"][3],
            )
            atom_list_1 = np.nonzero(atom_mask1)[0]
            atom_list_2 = np.nonzero(atom_mask2)[0]
            if np.size(atom_list_1) != 1 or np.size(atom_list_2) != 1:
                logger.info(f"Atom contact {i} not found for the input")
                continue
            atom_contact_specifics.append(
                (
                    atom_list_1.item(),
                    atom_list_2.item(),
                    pair["max_distance"],
                    pair["min_distance"],
                )
            )

        token_array, atom_array, atom_contact_specifics, _, _ = (
            ConstraintFeatureGenerator.expand_tokens(
                token_array,
                atom_array,
                atom_contact_specifics,
                atom_to_token_idx_dict,
                None,
            )
        )
        atom_contact_featurizer = ContactAtomFeaturizer(
            token_array=token_array, atom_array=atom_array
        )
        feature_dict["contact_atom"], _ = (
            atom_contact_featurizer.generate_spec_constraint(
                atom_contact_specifics,
                feature_type="continuous",
                shape=(len(token_array), len(token_array), 2),
            )
        )

        # build token-level contact
        atom_to_token_idx_dict = {}
        for idx, token in enumerate(token_array.tokens):
            for atom_idx in token.atom_indices:
                atom_to_token_idx_dict[atom_idx] = idx
        atom_to_token_idx = np.array(
            [atom_to_token_idx_dict[atom_idx] for atom_idx in range(len(atom_array))]
        )
        token_contact_specifics = []
        for i, pair in enumerate(token_contact_inputs):
            atom_mask1 = get_atom_mask_by_name(
                atom_array=atom_array,
                entity_id=pair["id1"][0],
                copy_id=pair["id1"][1],
                position=pair["id1"][2],
                atom_name=pair["id1"][3],
            )
            atom_mask2 = get_atom_mask_by_name(
                atom_array=atom_array,
                entity_id=pair["id2"][0],
                copy_id=pair["id2"][1],
                position=pair["id2"][2],
                atom_name=pair["id2"][3],
            )
            token_list_1 = atom_to_token_idx[atom_mask1]
            token_list_2 = atom_to_token_idx[atom_mask2]
            if np.size(token_list_1) == 0 or np.size(token_list_2) == 0:
                logger.info(f"Contact {i} not found for the input")
                continue
            token_contact_specifics.append(
                (token_list_1, token_list_2, pair["max_distance"], 0)
            )  # default min_distance=0

        contact_featurizer = ContactFeaturizer(
            token_array=token_array, atom_array=atom_array
        )
        feature_dict["contact"], _ = contact_featurizer.generate_spec_constraint(
            token_contact_specifics, feature_type="continuous"
        )

        # build pocket features
        pocket_specifics = []
        if pocket := constraint_param.get("pocket", {}):
            distance = pocket["max_distance"]
            binder = pocket["binder_chain"]

            assert len(binder) == 2

            atom_mask_binder = get_atom_mask_by_name(
                atom_array=atom_array,
                entity_id=binder["entity"],
                copy_id=binder["copy"],
            )

            binder_asym_id = torch.tensor(
                atom_array.asym_id_int[atom_mask_binder], dtype=torch.long
            )
            num_binder = binder_asym_id.unique().numel()
            if num_binder == 0:
                logger.info(f"Binder does not exist. {i},{num_binder}")
            elif num_binder > 1:
                logger.info(f"#Binders is more than 1. {i},{num_binder}")
            else:
                binder_token_list = atom_to_token_idx[atom_mask_binder]

                for j, pocket_res in enumerate(pocket["contact_residues"]):
                    pocket_res = (
                        ConstraintFeatureGenerator._canonicalize_pocket_res_format(
                            binder, pocket_res
                        )
                    )

                    atom_mask_pocket = get_atom_mask_by_name(
                        atom_array=atom_array,
                        entity_id=pocket_res["entity"],
                        copy_id=pocket_res["copy"],
                        position=pocket_res["position"],
                    )
                    pocket_token_list = atom_to_token_idx[atom_mask_pocket]

                    if np.size(pocket_token_list) == 0:
                        logger.info(f"Pocket not found: {i}:{j}")
                        continue

                    pocket_specifics.append(
                        (binder_token_list, pocket_token_list, distance)
                    )
                logger.info(f"#pocket:{len(pocket_specifics)}")

        pocket_featurizer = PocketFeaturizer(
            token_array=token_array, atom_array=atom_array
        )
        feature_dict["pocket"], _ = pocket_featurizer.generate_spec_constraint(
            pocket_specifics,
            feature_type="continuous",
        )

        # build substructure features
        substructure_specifics = {
            "token_indices": [],
            "token_coords": [],
        }
        if substructure := constraint_param.get("structure", {}):
            # TODO parse substructure specifics
            pass
        substructure_featurizer = SubStructureFeaturizer(
            token_array=token_array, atom_array=atom_array
        )
        feature_dict["substructure"] = substructure_featurizer.generate_spec_constraint(
            substructure_specifics, feature_type="one_hot"
        )

        logger.info(
            f"Loaded constraint feature: #atom contact:{len(atom_contact_specifics)} #contact:{len(token_contact_specifics)} #pocket:{len(pocket_specifics)}"
        )
        ConstraintFeatureGenerator._log_constraint_feature(
            atom_array, token_array, feature_dict
        )

        return feature_dict, token_array, atom_array

    # function set for training
    def generate(
        self,
        atom_array: AtomArray,
        token_array: TokenArray,
        sample_indice: pd.core.series.Series,
        pdb_indice: pd.core.series.Series,
        msa_features: dict[str, np.ndarray],
        max_entity_mol_id: int,
        full_atom_array: AtomArray,
    ) -> tuple[
        TokenArray,
        AtomArray,
        dict[str, np.ndarray],
        dict[str, torch.Tensor],
        dict[str, Any],
        torch.Tensor,
        AtomArray,
    ]:
        """Generate all constraint features

        Args:
            idx: Data index
            atom_array: Atom array data
            token_array: Token array data
            sample_indice: Sample index information
            pdb_indice: PDB index information
            msa_features: MSA features
            max_entity_mol_id: Maximum entity mol id
            full_atom_array: Full atom array data
        Returns:
            Dictionary of constraint features
        """
        # Setup constraint generator
        constraint_generator = self._setup_constraint_generator(sample_indice)

        # Get base features for constraint featurizer
        features_dict = self._get_base_features(token_array, atom_array)

        # Generate contact atom features
        (_, contact_atom_featurizer, contact_pairs, tokens_w_atom_contact) = (
            self._generate_contact_atom_features(
                atom_array,
                token_array,
                constraint_generator,
                features_dict,
                pdb_indice,
            )
        )

        # Expand token according to atom-contact pairs
        if len(contact_pairs) > 0:
            logger.info("Expanding tokens for contact atom constraint feature")
            token_array, atom_array, contact_pairs, token_map, full_atom_array = (
                ConstraintFeatureGenerator.expand_tokens(
                    token_array,
                    atom_array,
                    contact_pairs,
                    features_dict["atom_to_token_dict"],
                    full_atom_array,
                )
            )
            features_dict = self._get_base_features(token_array, atom_array)
        else:
            token_map = {}
        # make atom contact feature
        contact_atom_constraint_feature, tokens_w_atom_contact = (
            contact_atom_featurizer.generate_spec_constraint(
                contact_pairs,
                feature_type=self.constraint.get("contact_atom", {}).get(
                    "feature_type", "continuous"
                ),
                shape=(len(token_array), len(token_array), 2),
            )
        )

        # Expand MSA features
        if len(msa_features) > 0 and len(contact_pairs) > 0:
            msa_features = self.expand_msa_features(msa_features, token_map)

        # Generate pocket features
        pocket_constraint_feature, _, tokens_w_pocket = self._generate_pocket_features(
            atom_array,
            token_array,
            constraint_generator,
            sample_indice,
        )

        # Generate contact features
        contact_constraint_feature, _, _, tokens_w_contact = (
            self._generate_contact_features(
                atom_array,
                token_array,
                constraint_generator,
                features_dict,
                pdb_indice,
            )
        )

        # Generate substructure features
        (
            substructure_constraint_feature,
            substructure_featurizer,
            tokens_w_substructure,
        ) = self._generate_substructure_features(
            atom_array,
            token_array,
            constraint_generator,
            sample_indice,
        )

        # Combine features
        constraint_feature = {
            "contact": contact_constraint_feature,
            "pocket": pocket_constraint_feature,
            "contact_atom": contact_atom_constraint_feature,
            "substructure": substructure_constraint_feature,
        }

        # change entity_mol_id in case of permutation of constraint pairs
        featured_tokens = (
            tokens_w_contact
            | tokens_w_atom_contact
            | tokens_w_pocket
            | tokens_w_substructure
        )

        if max_entity_mol_id is not None:
            atom_array, full_atom_array = self.change_entity_mol_id(
                token_array,
                atom_array,
                max_entity_mol_id,
                full_atom_array,
                featured_tokens,
            )

        # Log feature statistics
        feature_info = self._get_feature_statistics(
            constraint_feature,
            atom_array,
            token_array,
            features_dict,
            contact_atom_featurizer,
            substructure_featurizer,
        )
        log_constraint = self._log_feature_statistics(feature_info)
        return (
            token_array,
            atom_array,
            msa_features,
            constraint_feature,
            feature_info,
            log_constraint,
            full_atom_array,
        )

    @staticmethod
    def expand_tokens(
        token_array: TokenArray,
        atom_array: AtomArray,
        contact_pairs: list[tuple[int, int, float, float]],
        atom_to_token: dict[int, int],
        full_atom_array: AtomArray,
    ) -> tuple[TokenArray, AtomArray, list[tuple[int, int, float, float]]]:
        """
        Expand selected tokens into atom-level tokens and update related arrays.

        Args:
            token_array: Original token array
            atom_array: Original atom array
            contact_pairs: Original contact pairs
            atom_to_token: Atom to token mapping
            full_atom_array: Full atom array
        Returns:
            Updated token array, atom array and transformed constraint pairs
        """
        # Update token array
        tokens_to_expand = set()
        for atom_i, atom_j, _, _ in contact_pairs:
            tokens_to_expand.add(atom_to_token[atom_i])
            tokens_to_expand.add(atom_to_token[atom_j])

        if len(tokens_to_expand) == 0:
            return token_array, atom_array, contact_pairs, {}, full_atom_array

        new_tokens = []
        # Maps old token idx to list of new token indices
        token_map = {}
        curr_token_idx = 0

        for old_token_idx, token in enumerate(token_array):
            if old_token_idx in tokens_to_expand:
                # Check if token represents standard residue
                centre_atom = atom_array[token.centre_atom_index]
                if (
                    centre_atom.res_name in STD_RESIDUES
                    and centre_atom.mol_type != "ligand"
                ):
                    # Expand token into atom-level tokens
                    atom_tokens = []
                    for atom_idx, atom_name in zip(
                        token.atom_indices, token.atom_names
                    ):
                        atom = atom_array[atom_idx]
                        atom_token = Token(ELEMS[atom.element])
                        atom_token.atom_indices = [atom_idx]
                        atom_token.atom_names = [atom_name]
                        atom_token.centre_atom_index = atom_idx
                        atom_tokens.append(atom_token)
                    new_tokens.extend(atom_tokens)
                    token_map[old_token_idx] = list(
                        range(curr_token_idx, curr_token_idx + len(atom_tokens))
                    )
                    curr_token_idx += len(atom_tokens)
                else:
                    new_tokens.append(token)
                    token_map[old_token_idx] = [curr_token_idx]
                    curr_token_idx += 1
            else:
                new_tokens.append(token)
                token_map[old_token_idx] = [curr_token_idx]
                curr_token_idx += 1

        updated_token_array = TokenArray(new_tokens)

        # Create atom_idx to new_token_idx mapping for expanded tokens
        atom_to_new_token = {}
        for new_token_idx, token in enumerate(updated_token_array):
            for atom_idx in token.atom_indices:
                atom_to_new_token[atom_idx] = new_token_idx

        # Update atom array annotations
        atom_array.centre_atom_mask = np.zeros(len(atom_array), dtype=bool)
        for token in updated_token_array:
            atom_array.centre_atom_mask[token.centre_atom_index] = True

        # Update tokatom_idx and distogram_rep_atom_mask
        expanded_atoms = set()

        for tokens in token_map.values():
            if len(tokens) > 1:  # Expanded tokens
                for token_idx in tokens:
                    token = updated_token_array[token_idx]
                    expanded_atoms.update(token.atom_indices)

        atom_array.tokatom_idx = np.array(
            [
                0 if i in expanded_atoms else idx
                for i, idx in enumerate(atom_array.tokatom_idx)
            ]
        )

        atom_array.distogram_rep_atom_mask = np.array(
            [
                1 if i in expanded_atoms else mask
                for i, mask in enumerate(atom_array.distogram_rep_atom_mask)
            ]
        )
        if len(expanded_atoms) > 0:
            logger.info(f"Expanded atoms: {expanded_atoms}")

        # Create mapping between atom_array and full_atom_array atoms
        if full_atom_array is not None:
            expanded_atom_keys = set()
            for atom_idx in expanded_atoms:
                atom = atom_array[atom_idx]
                # Create unique key using chain_id, res_id and atom_name
                atom_key = (atom.chain_id, atom.res_id, atom.atom_name)
                expanded_atom_keys.add(atom_key)

            # Update full_atom_array centre_atom_mask using the mapping
            full_atom_array.centre_atom_mask = np.array(
                [
                    (
                        1
                        if (atom.chain_id, atom.res_id, atom.atom_name)
                        in expanded_atom_keys
                        else mask
                    )
                    for atom, mask in zip(
                        full_atom_array, full_atom_array.centre_atom_mask
                    )
                ]
            )
        # Transform constraint pairs using atom_to_new_token mapping
        transformed_pairs = []
        for atom_i, atom_j, min_dist, max_dist in contact_pairs:
            new_token_i = atom_to_new_token[atom_i]
            new_token_j = atom_to_new_token[atom_j]
            transformed_pairs.append((new_token_i, new_token_j, min_dist, max_dist))

        return (
            updated_token_array,
            atom_array,
            transformed_pairs,
            token_map,
            full_atom_array,
        )

    def expand_msa_features(
        self,
        msa_features: dict[str, np.ndarray],
        token_map: dict[
            int, list[int]
        ],  # Maps old token idx to list of new token indices
    ) -> dict[str, np.ndarray]:
        """
        Expand MSA features for expanded tokens.

        Args:
            msa_features: Original MSA features
            token_map: Mapping from old token indices to new token indices

        Returns:
            Updated MSA features with expanded tokens
        """
        # Calculate new number of tokens
        num_new_tokens = max(max(new_idxs) for new_idxs in token_map.values()) + 1

        old_indices = []
        new_indices = []
        for old_idx, new_idxs in token_map.items():
            old_indices.extend([old_idx] * len(new_idxs))
            new_indices.extend(new_idxs)
        old_indices = np.array(old_indices, dtype=int)
        new_indices = np.array(new_indices, dtype=int)

        # For sequence-based features (msa, has_deletion, deletion_value)
        for feat_name in ["msa", "has_deletion", "deletion_value"]:
            if feat_name not in msa_features:
                continue

            feat = msa_features[feat_name]
            num_seqs = feat.shape[0]  # Number of sequences in MSA

            # Create new feature array
            new_feat = np.zeros((num_seqs, num_new_tokens), dtype=feat.dtype)

            # Copy features according to token mapping
            new_feat = np.zeros((num_seqs, num_new_tokens), dtype=feat.dtype)
            new_feat[:, new_indices] = feat[:, old_indices]
            msa_features[feat_name] = new_feat

        # Handle deletion_mean (1D array)
        if "deletion_mean" in msa_features:
            feat = msa_features["deletion_mean"]
            new_feat = np.zeros(num_new_tokens, dtype=feat.dtype)
            new_feat[new_indices] = feat[old_indices]
            msa_features["deletion_mean"] = new_feat

        # Handle profile (2D array: tokens x channels)
        if "profile" in msa_features:
            feat = msa_features["profile"]
            num_channels = feat.shape[1]
            new_feat = np.zeros((num_new_tokens, num_channels), dtype=feat.dtype)
            new_feat[new_indices, :] = feat[old_indices, :]
            msa_features["profile"] = new_feat

        return msa_features

    def change_entity_mol_id(
        self,
        token_array: TokenArray,
        atom_array: AtomArray,
        max_entity_mol_id: int,
        full_atom_array: AtomArray,
        featured_tokens: set[int],
    ) -> tuple[AtomArray, AtomArray]:
        """Update entity_mol_id for atoms involved in constraints"""
        if max_entity_mol_id is None or len(featured_tokens) == 0:
            return atom_array, full_atom_array

        # Get atom indices for all constrained tokens
        constrained_atom_indices = set()
        centre_atom_indices = token_array.get_annotation("centre_atom_index")

        for token_idx in featured_tokens:
            constrained_atom_indices.add(centre_atom_indices[token_idx])

        # Get mol_ids for constrained atoms
        constrained_mol_ids = set(atom_array.mol_id[list(constrained_atom_indices)])

        # Create mapping for new entity_mol_ids
        new_id = max_entity_mol_id + 1
        id_mapping = {}
        for old_id in constrained_mol_ids:
            id_mapping[old_id] = new_id
            new_id += 1

        # Update entity_mol_ids in atom_array
        for old_id, new_id in id_mapping.items():
            mask = atom_array.mol_id == old_id
            atom_array.entity_mol_id[mask] = new_id
            mask = full_atom_array.mol_id == old_id
            full_atom_array.entity_mol_id[mask] = new_id

        return atom_array, full_atom_array

    def _get_chain_interface_mask(
        self,
        pdb_indice: pd.core.series.Series,
        atom_array_chain_id: np.array,
    ) -> tuple[list[torch.Tensor], list[torch.Tensor]]:
        df = pdb_indice.copy()

        def get_atom_mask(row):
            chain_1_mask = atom_array_chain_id == row["chain_1_id"]
            if row["type"] == "chain":
                chain_2_mask = chain_1_mask
            else:
                chain_2_mask = atom_array_chain_id == row["chain_2_id"]
            chain_1_mask = torch.tensor(chain_1_mask).bool()
            chain_2_mask = torch.tensor(chain_2_mask).bool()
            if chain_1_mask.sum() == 0 or chain_2_mask.sum() == 0:
                return None, None
            return chain_1_mask, chain_2_mask

        df["chain_1_mask"], df["chain_2_mask"] = zip(*df.apply(get_atom_mask, axis=1))
        df = df[df["chain_1_mask"].notna()]  # drop NaN

        chain_1_mask = df["chain_1_mask"].tolist()  # [N_eval, N_atom]
        chain_2_mask = df["chain_2_mask"].tolist()  # [N_eval, N_atom]

        return chain_1_mask, chain_2_mask

    def _setup_constraint_generator(self, sample_indice: pd.core.series.Series) -> Any:
        """Setup random generator with optional fixed seed"""
        if not self.constraint.get("fix_seed", False):
            return None

        constraint_seed = int(
            hashlib.sha256(sample_indice.pdb_id.encode("utf-8")).hexdigest(), 16
        ) % (2**32)
        logger.info(
            f"[constraint_seed seed]: {constraint_seed} for pdb {sample_indice.pdb_id}"
        )
        return torch.Generator().manual_seed(constraint_seed)

    def _get_base_features(
        self, token_array: TokenArray, atom_array: AtomArray
    ) -> dict[str, Any]:
        """Get base features for constraint featurizer"""
        features_dict = {}
        centre_atoms_indices = token_array.get_annotation("centre_atom_index")
        centre_atoms = atom_array[centre_atoms_indices]
        features_dict["asym_id"] = torch.Tensor(centre_atoms.asym_id_int).long()

        atom_to_token_idx_dict = {}
        for idx, token in enumerate(token_array.tokens):
            for atom_idx in token.atom_indices:
                atom_to_token_idx_dict[atom_idx] = idx

        # ensure the order of the atom_to_token_idx is the same as the atom_array
        atom_to_token_idx = [
            atom_to_token_idx_dict[atom_idx] for atom_idx in range(len(atom_array))
        ]
        features_dict["atom_to_token_idx"] = torch.Tensor(atom_to_token_idx).long()

        features_dict["is_dna"] = torch.Tensor(atom_array.is_dna).bool()
        features_dict["is_rna"] = torch.Tensor(atom_array.is_rna).bool()
        features_dict["is_ligand"] = torch.Tensor(atom_array.is_ligand).bool()
        features_dict["atom_to_token_dict"] = atom_to_token_idx_dict
        return features_dict

    def _generate_pocket_features(
        self,
        atom_array: AtomArray,
        token_array: TokenArray,
        generator: torch.Generator,
        sample_indice: pd.core.series.Series,
    ) -> tuple[torch.Tensor, Any, set[int]]:
        """Generate pocket constraint features"""
        # pocket feature
        pocket_featurizer = PocketFeaturizer(
            token_array=token_array,
            atom_array=atom_array,
            generator=generator,
        )
        use_random_pocket = np.random.rand() < self.constraint.get("pocket", {}).get(
            "prob", 0
        )
        if use_random_pocket:
            binder_asym_id = None
            if self.constraint["pocket"].get("spec_binder_chain", False):
                # find ligand binder
                if sample_indice.mol_1_type == "ligand":
                    binder_chain_id = str(sample_indice.chain_1_id)
                elif sample_indice.mol_2_type == "ligand":
                    binder_chain_id = str(sample_indice.chain_2_id)
                # find antibody binder
                elif (
                    f"{sample_indice.pdb_id.lower()}_{sample_indice.entity_1_id}"
                    in self.ab_top2_clusters
                ):
                    binder_chain_id = str(sample_indice.chain_1_id)
                elif (
                    f"{sample_indice.pdb_id.lower()}_{sample_indice.entity_2_id}"
                    in self.ab_top2_clusters
                ):
                    binder_chain_id = str(sample_indice.chain_2_id)
                else:
                    binder_chain_id = -1
                    logger.info(f"No binder found!")
                binder_asym_id = atom_array.asym_id_int[
                    atom_array.chain_id == binder_chain_id
                ]
                num_unique = len(np.unique(binder_asym_id))
                assert num_unique <= 1
                logger.info(f"found binder chains: {num_unique}")
                binder_asym_id = binder_asym_id[0] if num_unique == 1 else -1
            pocket_constraint_feature, constrained_tokens = pocket_featurizer.generate(
                size=self.constraint["pocket"].get("size", 0.0),
                feature_type=self.constraint["pocket"].get(
                    "feature_type", "continuous"
                ),
                max_distance_range=self.constraint["pocket"].get(
                    "max_distance_range", {}
                ),
                sample_group=self.constraint["pocket"].get("group", "complex"),
                distance_type=self.constraint["pocket"].get(
                    "distance_type", "center_atom"
                ),
                spec_binder_asym_id=binder_asym_id,
            )
        else:
            pocket_constraint_feature, constrained_tokens = pocket_featurizer.generate(
                {}, "continuous", size=0, distance_type=None
            )
        return pocket_constraint_feature, pocket_featurizer, constrained_tokens

    def _generate_contact_features(
        self,
        atom_array: AtomArray,
        token_array: TokenArray,
        generator: torch.Generator,
        features_dict: dict[str, Any],
        pdb_indice: pd.core.series.Series,
    ) -> tuple[torch.Tensor, Any, list[tuple[int, int]], set[int]]:
        """Generate contact constraint features"""
        contact_featurizer = ContactFeaturizer(
            token_array=token_array,
            atom_array=atom_array,
            generator=generator,
        )
        # contact_feature_prob =
        use_random_contact = np.random.rand() < self.constraint.get("contact", {}).get(
            "prob", 0
        )
        selected_pairs = []
        if use_random_contact:
            # collecte all interfaces in cropped_tokens
            interface_asym_pairs = []
            if self.constraint["contact"].get("group", "complex") == "interface":

                def _get_asym_id(atom_mask):
                    token_idx = features_dict["atom_to_token_idx"][
                        atom_mask.bool()
                    ].long()
                    asym_ids = features_dict["asym_id"][token_idx]
                    assert len(torch.unique(asym_ids)) == 1
                    return asym_ids[0].item()

                chain_1_mask, chain_2_mask = self._get_chain_interface_mask(
                    pdb_indice=pdb_indice, atom_array_chain_id=atom_array.chain_id
                )

                interface_asym_pairs = [
                    (_get_asym_id(c1), _get_asym_id(c2))
                    for c1, c2 in zip(chain_1_mask, chain_2_mask)
                    if not torch.equal(c1, c2)
                ]

                if len(interface_asym_pairs) == 0:
                    logger.info("Interface for constraint feature is not found")

            contact_constraint_feature, selected_pairs, constrained_tokens = (
                contact_featurizer.generate(
                    size=self.constraint["contact"].get("size", 0.0),
                    chain_pairs=interface_asym_pairs,
                    feature_type=self.constraint["contact"].get(
                        "feature_type", "continuous"
                    ),
                    max_distance_range=self.constraint["contact"].get(
                        "max_distance_range", {}
                    ),
                    min_distance_threshold=self.constraint["contact"].get(
                        "min_distance", 0.0
                    ),
                    sample_group=self.constraint["contact"].get("group", "complex"),
                    distance_type=self.constraint["contact"].get(
                        "distance_type", "center_atom"
                    ),
                )
            )
        else:
            contact_constraint_feature, selected_pairs, constrained_tokens = (
                contact_featurizer.generate(
                    {}, None, "continuous", size=0, distance_type=None
                )
            )
        return (
            contact_constraint_feature,
            contact_featurizer,
            selected_pairs,
            constrained_tokens,
        )

    def _generate_contact_atom_features(
        self,
        atom_array: AtomArray,
        token_array: TokenArray,
        generator: torch.Generator,
        features_dict: dict[str, Any],
        pdb_indice: pd.core.series.Series,
    ) -> tuple[torch.Tensor, Any, list[tuple[int, int]], set[int]]:
        """Generate contact atom constraint features"""
        contact_atom_featurizer = ContactAtomFeaturizer(
            token_array=token_array, atom_array=atom_array, generator=generator
        )
        use_random_contact_atom = np.random.rand() < self.constraint.get(
            "contact_atom", {}
        ).get("prob", 0)
        selected_pairs = []
        if use_random_contact_atom:
            # collecte all interfaces in cropped_tokens
            interface_asym_pairs = []
            if self.constraint["contact_atom"].get("group", "complex") == "interface":

                def _get_asym_id(atom_mask):
                    token_idx = features_dict["atom_to_token_idx"][
                        atom_mask.bool()
                    ].long()
                    asym_ids = features_dict["asym_id"][token_idx]
                    assert len(torch.unique(asym_ids)) == 1
                    return asym_ids[0].item()

                chain_1_mask, chain_2_mask = self._get_chain_interface_mask(
                    pdb_indice=pdb_indice, atom_array_chain_id=atom_array.chain_id
                )

                interface_asym_pairs = [
                    (_get_asym_id(c1), _get_asym_id(c2))
                    for c1, c2 in zip(chain_1_mask, chain_2_mask)
                    if not torch.equal(c1, c2)
                ]

                if len(interface_asym_pairs) == 0:
                    logger.info("Interface for constraint feature is not found")

            contact_atom_constraint_feature, selected_pairs, constrained_tokens = (
                contact_atom_featurizer.generate(
                    size=self.constraint["contact_atom"].get("size", 0.0),
                    chain_pairs=interface_asym_pairs,
                    feature_type=self.constraint["contact_atom"].get(
                        "feature_type", "continuous"
                    ),
                    max_distance_range=self.constraint["contact_atom"].get(
                        "max_distance_range", {}
                    ),
                    min_distance_threshold=self.constraint["contact_atom"].get(
                        "min_distance", 0.0
                    ),
                    sample_group=self.constraint["contact_atom"].get(
                        "group", "complex"
                    ),
                    distance_type=self.constraint["contact_atom"].get(
                        "distance_type", "atom"
                    ),
                )
            )
        else:
            contact_atom_constraint_feature, selected_pairs, constrained_tokens = (
                contact_atom_featurizer.generate(
                    {}, None, "continuous", size=0, distance_type=None
                )
            )
        return (
            contact_atom_constraint_feature,
            contact_atom_featurizer,
            selected_pairs,
            constrained_tokens,
        )

    def _generate_substructure_features(
        self,
        atom_array: AtomArray,
        token_array: TokenArray,
        generator: torch.Generator,
        sample_indice: pd.core.series.Series,
    ) -> tuple[torch.Tensor, Any, set[int]]:
        """Generate substructure constraint features"""
        substructure_featurizer = SubStructureFeaturizer(
            token_array=token_array,
            atom_array=atom_array,
            generator=generator,
        )
        use_random_substructure = np.random.rand() < self.constraint.get(
            "substructure", {}
        ).get("prob", 0)
        if use_random_substructure:
            spec_asym_id = None
            if self.constraint["substructure"].get("spec_asym_id", False):
                # find ligand chain
                if sample_indice.mol_1_type == "ligand":
                    binder_chain_id = str(sample_indice.chain_1_id)
                    target_chain_id = str(sample_indice.chain_2_id)
                elif sample_indice.mol_2_type == "ligand":
                    binder_chain_id = str(sample_indice.chain_2_id)
                    target_chain_id = str(sample_indice.chain_1_id)
                # find antibody chain
                elif (
                    f"{sample_indice.pdb_id.lower()}_{sample_indice.entity_1_id}"
                    in self.ab_top2_clusters
                ):
                    binder_chain_id = str(sample_indice.chain_1_id)
                    target_chain_id = str(sample_indice.chain_2_id)
                elif (
                    f"{sample_indice.pdb_id.lower()}_{sample_indice.entity_2_id}"
                    in self.ab_top2_clusters
                ):
                    binder_chain_id = str(sample_indice.chain_2_id)
                    target_chain_id = str(sample_indice.chain_1_id)
                else:
                    target_chain_id = binder_chain_id = -1
                    logger.info(f"No specific chain found!")

                chain_choice = self.constraint["substructure"].get(
                    "spec_chain_type", "binder"
                )
                if chain_choice == "binder":
                    spec_chain_id = binder_chain_id
                elif chain_choice == "target":
                    spec_chain_id = target_chain_id
                else:
                    raise ValueError(
                        f"Invalid spec_chain_type: {self.constraint['substructure'].get('spec_chain_type', 'binder')}"
                    )

                if spec_chain_id != -1:
                    spec_asym_id = atom_array.asym_id_int[
                        atom_array.chain_id == spec_chain_id
                    ]
                    num_unique = len(np.unique(spec_asym_id))
                    assert num_unique <= 1
                    logger.info(f"found {chain_choice} chain: {num_unique}")
                    spec_asym_id = spec_asym_id[0] if num_unique == 1 else -1

            substructure_constraint_feature, constrained_tokens = (
                substructure_featurizer.generate(
                    mol_type_pairs=self.constraint["substructure"].get(
                        "mol_type_pairs", {}
                    ),
                    feature_type=self.constraint["substructure"].get(
                        "feature_type", "one_hot"
                    ),
                    size=self.constraint["substructure"].get("size", 0),
                    ratios=self.constraint["substructure"].get(
                        "ratios", {"full": [0.0, 0.5, 1.0], "partial": 0.3}
                    ),
                    coord_noise_scale=self.constraint["substructure"].get(
                        "coord_noise_scale", 0.05
                    ),
                    spec_asym_id=spec_asym_id,
                )
            )
        else:
            substructure_constraint_feature, constrained_tokens = (
                substructure_featurizer.generate(
                    mol_type_pairs={},
                    feature_type="one_hot",
                    size=0,
                    ratios={"full": [0.0, 0.5, 1.0], "partial": 0.3},
                    coord_noise_scale=0.05,
                    spec_asym_id=None,
                )
            )
        return (
            substructure_constraint_feature,
            substructure_featurizer,
            constrained_tokens,
        )

    def _get_feature_statistics(
        self,
        constraint_feature: dict[str, torch.Tensor],
        atom_array: AtomArray,
        token_array: TokenArray,
        features_dict: dict[str, Any],
        contact_atom_featurizer: Any,
        substructure_featurizer: Any,
    ) -> dict[str, Any]:
        """Log statistics about generated features"""
        token_idx_1, token_idx_2 = torch.nonzero(
            torch.triu(constraint_feature["contact"][..., 1]), as_tuple=True
        )
        asym_id_1 = features_dict["asym_id"][token_idx_1]
        asym_id_2 = features_dict["asym_id"][token_idx_2]

        res_id_1 = atom_array.res_id[
            token_array[token_idx_1].get_annotation("centre_atom_index")
        ]
        res_id_2 = atom_array.res_id[
            token_array[token_idx_2].get_annotation("centre_atom_index")
        ]
        contact_distance = constraint_feature["contact"][token_idx_1, token_idx_2, 1]

        # logging contact atom feature
        atom_idx_1, atom_idx_2 = torch.nonzero(
            torch.triu(constraint_feature["contact_atom"][..., 1]), as_tuple=True
        )
        contact_atom_real_distance = (
            contact_atom_featurizer.get_real_distance(atom_idx_1, atom_idx_2)
            if len(atom_idx_1) > 0
            else None
        )
        contact_atom_max_distance = constraint_feature["contact_atom"][
            atom_idx_1, atom_idx_2, 1
        ]
        contact_atom_min_distance = constraint_feature["contact_atom"][
            atom_idx_1, atom_idx_2, 0
        ]
        num_contact_atom = contact_atom_max_distance.shape[0]

        # logging pocket feature
        binder_idx, pocket_idx = torch.nonzero(
            constraint_feature["pocket"].squeeze(-1), as_tuple=True
        )
        binder_idx = binder_idx.unique()
        pocket_idx = pocket_idx.unique()
        if binder_idx.numel() > 1:
            pocket_distance = constraint_feature["pocket"][binder_idx[0], pocket_idx, 0]
        else:
            pocket_distance = -1

        # logging substructure feature
        sub_structure_info = substructure_featurizer.analyze_features(
            constraint_feature["substructure"]
        )

        feature_info = {
            "contact_info": {
                "asym_id": torch.stack([asym_id_1, asym_id_2]),
                "token_id": torch.stack([token_idx_1, token_idx_2]),
                "res_id": torch.tensor(np.stack([res_id_1, res_id_2])),
                "distance": contact_distance,
            },
            "contact_atom_info": {
                "distance": contact_atom_real_distance,
                "max_distance": contact_atom_max_distance,
                "min_distance": contact_atom_min_distance,
                "num_contact_atom": num_contact_atom,
            },
            "pocket_info": {
                "binder_tokenid": binder_idx,
                "pocket_tokenid": pocket_idx,
                "distance": pocket_distance,
            },
            "substructure_info": sub_structure_info,
        }

        return feature_info

    def _log_feature_statistics(self, feature_info: dict[str, Any]) -> dict[str, Any]:
        log_constraint = {}
        if feature_info.get("pocket_info", None) is not None:
            binder_tokens, pocket_tokens, distance = feature_info[
                "pocket_info"
            ].values()
            pocket_msg = ";".join(
                [
                    ",".join(map(str, binder_tokens.flatten().tolist())),
                    ",".join(map(str, pocket_tokens.flatten().tolist())),
                ]
            )
            log_constraint["pocket_msg"] = pocket_msg
            log_constraint["pocket_N_binder"] = torch.tensor(
                feature_info["pocket_info"]["binder_tokenid"].shape[0]
            )
            log_constraint["pocket_N_pocket"] = torch.tensor(
                feature_info["pocket_info"]["pocket_tokenid"].shape[0]
            )
        if feature_info.get("contact_info", None) is not None:
            asym_id, token_id, res_id, distance = feature_info["contact_info"].values()
            N_contact = asym_id.shape[-1]
            contact_msg = ";".join(
                [
                    ",".join(map(str, asym_id.flatten().tolist())),
                    ",".join(map(str, token_id.flatten().tolist())),
                    ",".join(map(str, distance.flatten().tolist())),
                ]
            )
            log_constraint["contact_N_pair"] = torch.tensor(N_contact)
            log_constraint["contact_msg"] = contact_msg
        if feature_info.get("contact_atom_info", None) is not None:
            distance, max_distance, min_distance, num_contact_atom = feature_info[
                "contact_atom_info"
            ].values()
            N_contact = num_contact_atom
            log_constraint["contact_atom_N_pair"] = torch.tensor(N_contact)
            log_constraint["contact_atom_distance"] = distance
            log_constraint["contact_atom_max_distance"] = max_distance
            log_constraint["contact_atom_min_distance"] = min_distance

        if feature_info.get("substructure_info", None) is not None:
            log_constraint["substructure_active_tokens"] = torch.tensor(
                feature_info["substructure_info"]["num_active_tokens"]
            )
            log_constraint["substructure_active_token_ratio"] = torch.tensor(
                feature_info["substructure_info"]["active_token_ratio"]
            )
            log_constraint["substructure_bin0_cnt"] = torch.tensor(
                feature_info["substructure_info"]["distance_distribution"][
                    "bin_0_count"
                ]
            )
            log_constraint["substructure_bin1_cnt"] = torch.tensor(
                feature_info["substructure_info"]["distance_distribution"][
                    "bin_1_count"
                ]
            )
            log_constraint["substructure_bin2_cnt"] = torch.tensor(
                feature_info["substructure_info"]["distance_distribution"][
                    "bin_2_count"
                ]
            )
            log_constraint["substructure_bin3_cnt"] = torch.tensor(
                feature_info["substructure_info"]["distance_distribution"][
                    "bin_3_count"
                ]
            )

        return log_constraint


class ConstraintFeaturizer(object):
    def __init__(
        self,
        token_array: TokenArray,
        atom_array: AtomArray,
        pad_value: float = 0,
        generator=None,
    ):
        self.token_array = token_array
        self.atom_array = atom_array
        self.pad_value = pad_value
        self.generator = generator
        self._get_base_info()

    @staticmethod
    def one_hot_encoder(feature: torch.Tensor, num_classes: int):
        # Create mask for padding values (-1)
        pad_mask = feature == -1

        # Replace -1 with 0 temporarily for F.one_hot
        feature = torch.where(pad_mask, torch.zeros_like(feature), feature)

        # Convert to one-hot
        one_hot = F.one_hot(feature, num_classes=num_classes).float()

        # Zero out the one-hot vectors for padding positions
        one_hot[pad_mask] = 0.0

        return one_hot

    def encode(self, feature: torch.Tensor, feature_type: str, **kwargs):
        if feature_type == "one_hot":
            return ConstraintFeaturizer.one_hot_encoder(
                feature, num_classes=kwargs.get("num_classes", -1)
            )
        elif feature_type == "continuous":
            return feature
        else:
            raise RuntimeError(f"Invalid feature_type: {feature_type}")

    def _get_base_info(self):
        token_centre_atom_indices = self.token_array.get_annotation("centre_atom_index")
        centre_atoms = self.atom_array[token_centre_atom_indices]
        self.asymid = torch.tensor(centre_atoms.asym_id_int, dtype=torch.long)
        self.is_ligand = torch.tensor(centre_atoms.is_ligand, dtype=torch.bool)
        self.is_protein = torch.tensor(centre_atoms.is_protein, dtype=torch.bool)
        self.entity_type_dict = {"P": self.is_protein, "L": self.is_ligand}

    def _get_generation_basics(self, distance_type: str = "center_atom"):
        token_centre_atom_indices = self.token_array.get_annotation("centre_atom_index")
        centre_atoms = self.atom_array[token_centre_atom_indices]

        # is_resolved mask
        self.token_resolved_mask = torch.tensor(
            centre_atoms.is_resolved, dtype=torch.bool
        )
        self.token_resolved_maskmat = (
            self.token_resolved_mask[:, None] * self.token_resolved_mask[None, :]
        )

        # distance matrix
        if distance_type == "center_atom":
            # center atom distance
            self.token_distance = torch.tensor(
                cdist(centre_atoms.coord, centre_atoms.coord), dtype=torch.float64
            )
        elif distance_type == "any_atom":
            # any atom distance
            all_atom_resolved_mask = (
                self.atom_array.is_resolved[:, None]
                * self.atom_array.is_resolved[None, :]
            )
            all_atom_distance = cdist(self.atom_array.coord, self.atom_array.coord)
            all_atom_distance[~all_atom_resolved_mask] = np.inf

            token_atoms_num = [
                len(_atoms)
                for _atoms in self.token_array.get_annotation("atom_indices")
            ]
            atom_token_num = np.repeat(
                np.arange(len(self.token_array)), token_atoms_num
            )

            self.token_distance = torch.zeros(
                (len(centre_atoms), len(centre_atoms)), dtype=torch.float64
            )
            for i, j in np.ndindex(self.token_distance.shape):
                atom_pairs_mask = np.ix_(atom_token_num == i, atom_token_num == j)
                self.token_distance[i, j] = np.min(all_atom_distance[atom_pairs_mask])
        elif distance_type == "atom":
            raise ValueError(
                "Not implement in this class, please use ContactAtomFeaturizer"
            )
        else:
            raise ValueError(f"Not recognized distance_type: {distance_type}")

    def generate(self):
        pass

    def generate_spec_constraint(self):
        pass


class ContactFeaturizer(ConstraintFeaturizer):
    def __init__(self, **kargs):
        super().__init__(**kargs)

    def get_valid_contact_feature(
        self,
        valid_contact_type: str,
        max_distance_threshold: float,
        min_distance_threshold: float,
    ) -> torch.Tensor:
        """
        Find all valid pairs of entities that satisfy the given contact distance requirements.

        Parameters:
         - valid_contact_type: A two-charactor type to represent entity pairs under consideration.
            e.g. PP, PL, P_
         - max_distance_threshold: The maximum allowable distance for a contact to be considered valid.
         - min_distance_threshold: The minimum allowable distance for a contact to be considered valid.
        Returns:
         - shape=(N_token, N_token), A matrix generate contact within the specified distance range.
        """
        # get valid contact type
        query_type, key_type = valid_contact_type

        if key_type == "_":
            # intra chain contact
            assert query_type == "P", "only support intra-protein contact for now"

            valid_chain_mask = self.asymid[:, None] == self.asymid[None, :]

            # skip closest squence pairs
            n_neighbor = 20  # TODO: tune this parameter
            valid_chain_mask_right = torch.triu(valid_chain_mask, diagonal=n_neighbor)
            valid_chain_mask_left = torch.tril(valid_chain_mask, diagonal=-n_neighbor)
            valid_chain_mask = valid_chain_mask_right | valid_chain_mask_left

        else:
            # inter chain contact
            assert (
                query_type in "PL" and key_type in "PL"
            ), f"[Error]contact type not support: {valid_contact_type}"

            valid_type_mask = (
                self.entity_type_dict[query_type][None, :]
                & self.entity_type_dict[key_type][:, None]
            )
            valid_type_mask |= (
                self.entity_type_dict[key_type][None, :]
                & self.entity_type_dict[query_type][:, None]
            )
            # get different chain mask
            valid_chain_mask = self.asymid[:, None] != self.asymid[None, :]
            valid_chain_mask &= valid_type_mask

        # get min & max distance threshold
        if min_distance_threshold == -1:  # default false, hard to determine
            # random select a min threshold in [0, max_distance_threshold), but may lead to zero contact pairs
            min_distance_threshold = (
                torch.zeros(1).uniform_(to=max_distance_threshold).item()
            )
        valid_dist_mask = (self.token_distance <= max_distance_threshold) & (
            self.token_distance >= min_distance_threshold
        )

        # make feature
        contact_valid_mask = (
            valid_chain_mask & self.token_resolved_maskmat & valid_dist_mask
        )
        return contact_valid_mask

    def _get_constraint_size(self, group: str, size: Any) -> list[int]:
        """
        If size is not fixed, then we generate it randomly for each group

        Args:
            - group: groups types
        """
        if group == "complex":
            k = 1
        elif group == "interface":
            N_asym = torch.unique(self.asymid).shape[0]
            k = (N_asym * N_asym - N_asym) // 2
        if size < 1 and size > 0:
            samples = torch.zeros(k).geometric_(size).int().tolist()
            return samples
        elif size >= 1 and isinstance(size, int):
            return [size] * k
        else:
            raise NotImplementedError

    def _sample_contacts(
        self,
        contact_valid_mask: torch.Tensor,
        contact_distance_mat: torch.Tensor,
        size: list[int],
        sample_group: str,
        chain_pairs: list[tuple[int, int]],
    ) -> tuple[torch.Tensor, list[tuple[int, int, float, float]]]:
        """
        Randomly select contact from all valid contact pairs

        Args:
            - contact_valid_mask: shape=(N_token, N_token), all valid contact pairs, bool
            - size: how many contacts to sample
            - sample_group: support sample group by 1. whole complex 2. interface
        """
        result_mat = torch.full(
            (contact_valid_mask.shape[0], contact_valid_mask.shape[1], 2),
            fill_value=self.pad_value,
            dtype=torch.float32,
        )
        selected_pairs = []

        if not contact_valid_mask.any().item():
            return result_mat, selected_pairs

        def _sample(valid_mask, cur_size):
            nonlocal selected_pairs

            valid_indices = torch.nonzero(torch.triu(valid_mask))
            selected_indices = valid_indices[
                torch.randperm(valid_indices.shape[0], generator=self.generator)[
                    :cur_size
                ]
            ]

            selected_indices = tuple(zip(*selected_indices))
            if len(selected_indices) == 0:
                return

            # Convert to pairs format
            for i, j in zip(*selected_indices):
                min_dist, max_dist = contact_distance_mat[i, j]
                selected_pairs.append(
                    (i.item(), j.item(), min_dist.item(), max_dist.item())
                )
                selected_pairs.append(
                    (j.item(), i.item(), min_dist.item(), max_dist.item())
                )  # Add symmetric pair

            # add symmetry indices
            selected_indices = (
                selected_indices[0] + selected_indices[1],
                selected_indices[1] + selected_indices[0],
            )
            result_mat[selected_indices] = contact_distance_mat[selected_indices]
            return

        # sample contacts by complex
        if sample_group == "complex":
            _sample(contact_valid_mask, size[0])
        # if group by interface, get all unique interfaces, and sample contact from each interface
        elif sample_group == "interface":
            # sample contacts from each interface iteratively

            idx = 0
            for asym1, asym2 in chain_pairs:
                asym1_mask, asym2_mask = self.asymid == asym1, self.asymid == asym2

                cur_interface_mask = asym1_mask[..., None, :] & asym2_mask[..., :, None]
                cur_interface_mask |= (
                    asym1_mask[..., :, None] & asym2_mask[..., None, :]
                )
                valid_contacts_interface = contact_valid_mask & cur_interface_mask
                _sample(valid_contacts_interface, size[idx])
                idx += 1
        else:
            raise NotImplementedError
        return result_mat, selected_pairs

    def generate(
        self,
        max_distance_range: dict[str, tuple[float, float]],
        sample_group: str,
        feature_type: str,
        size: Any,
        distance_type: str,
        min_distance_threshold: float = 0.0,
        chain_pairs: list[tuple[int, int]] = [],
        **kwargs,
    ) -> tuple[torch.Tensor, list[tuple[int, int, float, float]], set[int]]:
        """
        training & evaluation
        """
        constrained_tokens = set()
        if size == 0:
            return (
                self.encode(
                    torch.full(
                        (self.asymid.shape[0], self.asymid.shape[0], 2),
                        fill_value=self.pad_value,
                        dtype=torch.float32,
                    ),
                    feature_type=feature_type,
                ),
                [],
                constrained_tokens,
            )
        self._get_generation_basics(distance_type=distance_type)

        # contact mask
        n_token = len(self.asymid)
        contact_valid_mask = torch.zeros((n_token, n_token), dtype=torch.bool)
        contact_distance_mat = torch.zeros((n_token, n_token, 2), dtype=torch.float32)

        for contact_type, max_d_range in max_distance_range.items():
            # generate max_distance_mask for different contact type
            max_distance_threshold = torch.zeros(1).uniform_(*max_d_range).item()

            # get all valid contact pairs
            contact_mask = self.get_valid_contact_feature(
                contact_type,
                max_distance_threshold=max_distance_threshold,
                min_distance_threshold=min_distance_threshold,
            )
            contact_valid_mask |= contact_mask
            contact_distance_mat[contact_mask] = torch.tensor(
                [0, max_distance_threshold], dtype=torch.float32
            )

        # random select contact
        size = self._get_constraint_size(sample_group, size)
        sampled_contact_feature, selected_pairs = self._sample_contacts(
            contact_valid_mask,
            contact_distance_mat,
            size,
            sample_group,
            chain_pairs,
        )

        # encode the feature
        contact_feature = self.encode(
            feature=sampled_contact_feature, feature_type=feature_type
        )
        # Track constrained tokens
        constrained_tokens = set()
        for token_i, token_j, _, _ in selected_pairs:
            constrained_tokens.add(token_i)
            constrained_tokens.add(token_j)

        return contact_feature, selected_pairs, constrained_tokens

    def generate_spec_constraint(
        self,
        contact_specifics: list[tuple[int, int, float, float]],
        feature_type: str,
    ) -> tuple[torch.Tensor, set[int]]:
        """
        parse constraint from user specification
        """

        contact_feature = torch.full(
            (self.asymid.shape[0], self.asymid.shape[0], 2),
            fill_value=self.pad_value,
            dtype=torch.float32,
        )
        constrained_tokens = set()

        for token_list_1, token_list_2, max_distance, min_distance in contact_specifics:
            token_id_1 = token_list_1[
                torch.randint(
                    high=token_list_1.shape[0], size=(1,), generator=self.generator
                ).item()
            ]
            token_id_2 = token_list_2[
                torch.randint(
                    high=token_list_2.shape[0], size=(1,), generator=self.generator
                ).item()
            ]

            contact_feature[token_id_1, token_id_2, 1] = max_distance
            contact_feature[token_id_2, token_id_1, 1] = max_distance
            contact_feature[token_id_1, token_id_2, 0] = min_distance
            contact_feature[token_id_2, token_id_1, 0] = min_distance

            constrained_tokens.add(token_id_1)
            constrained_tokens.add(token_id_2)

        contact_feature = self.encode(
            feature=contact_feature, feature_type=feature_type
        )
        return contact_feature, constrained_tokens


class PocketFeaturizer(ConstraintFeaturizer):
    def __init__(self, **kargs):
        super().__init__(**kargs)

    def get_valid_pocket_feature(
        self,
        binder_pocket_type: str,
        max_distance_threshold: float,
        asym_list: list[int],
    ) -> torch.Tensor:
        """
        Parameters:
            - binder_pocket_type: PP, LP
            - max_distance_threshold
            - asym list
        Returns:
            - binder_pocket_valid_masks, shape=(N_asym, N_token), A matrix of all pocket tokens within the specified distance range.
            - max_distance_value, shape=(N_asym, N_token), A matrix of the max_distance_threshold to determine the pocket token for each binder chain
        """
        # get valid binder-pocket type, not symmetry
        query_type, key_type = binder_pocket_type
        valid_type_mask = (
            self.entity_type_dict[query_type][:, None]
            & self.entity_type_dict[key_type][None, :]
        )

        # get different chain mask
        diff_chain_mask = self.asymid[:, None] != self.asymid[None, :]

        # get distance mask
        # Note: To simplify the implementation, we only consider token-token distance, instead of any heavy atom distance for each residue for now.
        dist_mask = (
            (self.token_distance <= max_distance_threshold)
            & self.token_resolved_maskmat
            & diff_chain_mask
            & valid_type_mask
        )

        # pocket mask
        n_token = len(self.asymid)
        n_asym = len(asym_list)
        binder_pocket_valid_masks = torch.zeros((n_asym, n_token), dtype=torch.bool)

        for idx, asym_id in enumerate(asym_list):
            cur_chain_dist_mask = (
                self.asymid[:, None].expand(-1, self.asymid.shape[0]) == asym_id
            )
            if cur_chain_dist_mask[:, 0].sum() > 1:  # ensure num of binder tokens > 1
                cur_chain_dist_mask = cur_chain_dist_mask & dist_mask
                cur_valid_mask = cur_chain_dist_mask.any(dim=0)
                binder_pocket_valid_masks[idx] = cur_valid_mask
        return binder_pocket_valid_masks

    def _sample_pocket(
        self,
        binder_pocket_valid_masks: torch.Tensor,
        size: Any,
        asym_list: list[int],
        max_distance_value: torch.Tensor,
        spec_binder_asym_id: Any = None,
        **_,  # do not use
    ) -> torch.Tensor:
        pocket_dist_feature = torch.full(
            (self.asymid.shape[0], self.asymid.shape[0]),
            fill_value=self.pad_value,
            dtype=torch.float32,
        )
        if spec_binder_asym_id is None:
            # random select 1 binder with valid pocket
            binders_with_valid_pocket = torch.nonzero(
                binder_pocket_valid_masks.any(-1), as_tuple=True
            )[0]
            if len(binders_with_valid_pocket) == 0:
                return pocket_dist_feature
            selected_binder = binders_with_valid_pocket[
                torch.randperm(
                    binders_with_valid_pocket.shape[0], generator=self.generator
                )[0]
            ]
        else:
            selected_binder = torch.nonzero(
                asym_list == spec_binder_asym_id, as_tuple=True
            )[0].squeeze()
        selected_all_pocket_res = binder_pocket_valid_masks[selected_binder]
        binder_asym_id = (
            asym_list[selected_binder]
            if spec_binder_asym_id is None
            else spec_binder_asym_id
        )

        # random select k residues within the selected pocket
        selected_pocket_res_mask = torch.zeros(self.asymid.shape[0], dtype=torch.bool)
        valid_pockets = torch.nonzero(selected_all_pocket_res, as_tuple=True)[0]
        selected_pocket_res_indices = torch.randperm(
            valid_pockets.shape[0], generator=self.generator
        )[:size]
        selected_pocket_res_mask[valid_pockets[selected_pocket_res_indices]] = True

        binder_pocket_mask_mat = (self.asymid == binder_asym_id)[
            :, None
        ] * selected_pocket_res_mask[None, :]

        # set distance threshold to selected pocket
        # [binder_mask, pocket_mask]
        binder_dist_mat = max_distance_value[selected_binder][selected_pocket_res_mask][
            None, :
        ].expand((self.asymid == binder_asym_id).sum(), -1)

        pocket_dist_feature[binder_pocket_mask_mat] = binder_dist_mat.reshape(-1)

        return pocket_dist_feature

    def generate(
        self,
        max_distance_range: dict[str, tuple[float, float]],
        feature_type: str,
        size: Any,
        distance_type: Any,
        spec_binder_asym_id: int = None,
        **_,  # do not use
    ) -> tuple[torch.Tensor, set[int]]:
        """
        trainint & evaluation
        """
        constrained_tokens = set()
        if size == 0 or spec_binder_asym_id == -1:
            pocket_dist_feature = self.encode(
                torch.full(
                    (self.asymid.shape[0], self.asymid.shape[0], 1),
                    fill_value=self.pad_value,
                    dtype=torch.float32,
                ),
                feature_type=feature_type,
            )  # [..., N_token, 1]

            return pocket_dist_feature, constrained_tokens
        self._get_generation_basics(distance_type=distance_type)

        # get all binder-pocket pairs & masks
        n_token = len(self.asymid)
        asym_list = torch.unique(self.asymid, sorted=True)
        n_asym = len(asym_list)
        binder_pocket_valid_masks = torch.zeros((n_asym, n_token), dtype=torch.bool)
        max_distance_value = torch.zeros((n_asym, n_token), dtype=torch.float32)
        for binder_pocket_type, max_d_range in max_distance_range.items():
            # generate max_distance_mask for different binder_pocket type
            max_distance_threshold = torch.zeros(1).uniform_(*max_d_range).item()

            # get all valid binder-pocket pairs
            valid_pocket_token_mask = self.get_valid_pocket_feature(
                binder_pocket_type,
                max_distance_threshold=max_distance_threshold,
                asym_list=asym_list,
            )
            binder_pocket_valid_masks |= valid_pocket_token_mask
            max_distance_value[valid_pocket_token_mask] = max_distance_threshold

        # random select k residues from pocket, only consider one binder
        size = self._get_constraint_size(size)
        sampled_pocket_feature = self._sample_pocket(
            binder_pocket_valid_masks,
            size,
            asym_list,
            max_distance_value,
            spec_binder_asym_id,
        ).unsqueeze(-1)

        # encode the feature
        pocket_dist_feature = self.encode(
            feature=sampled_pocket_feature, feature_type=feature_type
        )
        # Track constrained tokens
        nonzero_indices = torch.nonzero(pocket_dist_feature)
        for i, j, _ in nonzero_indices:
            constrained_tokens.add(i.item())
            constrained_tokens.add(j.item())

        return pocket_dist_feature, constrained_tokens

    def generate_spec_constraint(
        self, pocket_specifics, feature_type: str
    ) -> tuple[torch.Tensor, set[int]]:
        """
        parse constraint from user specification
        """
        pocket_dist_mat = torch.full(
            (self.asymid.shape[0], self.asymid.shape[0], 1),
            fill_value=self.pad_value,
            dtype=torch.float32,
        )

        for binder_token_list, pocket_token_list, max_distance in pocket_specifics:
            pocket_token_id = pocket_token_list[
                torch.randint(
                    high=pocket_token_list.shape[0], size=(1,), generator=self.generator
                ).item()
            ]

            binder_token_idx = torch.tensor(binder_token_list)[:, None]
            pocket_dist_mat[binder_token_idx, pocket_token_id, 0] = max_distance

        pocket_dist_feature = self.encode(
            feature=pocket_dist_mat, feature_type=feature_type
        )

        constrained_tokens = set()
        nonzero_indices = torch.nonzero(pocket_dist_feature)
        for binder_token_id, pocket_token_id, _ in nonzero_indices:
            constrained_tokens.add(binder_token_id)
            constrained_tokens.add(pocket_token_id)
        return pocket_dist_feature, constrained_tokens

    def _get_constraint_size(self, size) -> list[int]:
        """
        If size is not fixed, then we generate it randomly for each group
        """
        if size < 1 and size > 0:
            # TODO: to be determined!
            samples = torch.zeros(1).geometric_(size).int().tolist()[0]
            return samples
        elif size >= 1 and isinstance(size, int):
            return size
        else:
            raise NotImplementedError


class ContactAtomFeaturizer(ContactFeaturizer):
    def __init__(
        self,
        token_array: TokenArray,
        atom_array: AtomArray,
        pad_value: float = 0,
        generator=None,
    ):
        self.token_array = token_array
        self.atom_array = atom_array
        self.pad_value = pad_value
        self.generator = generator
        self._get_base_info()

    def _get_base_info(self):
        self.asymid = torch.tensor(self.atom_array.asym_id_int, dtype=torch.long)
        self.is_ligand = torch.tensor(self.atom_array.is_ligand, dtype=torch.bool)
        self.is_protein = torch.tensor(self.atom_array.is_protein, dtype=torch.bool)
        self.entity_type_dict = {"P": self.is_protein, "L": self.is_ligand}

    def _get_generation_basics(self, distance_type="atom"):
        # is_resolved mask
        self.token_resolved_mask = torch.tensor(
            self.atom_array.is_resolved, dtype=torch.bool
        )
        self.token_resolved_maskmat = (
            self.token_resolved_mask[:, None] * self.token_resolved_mask[None, :]
        )

        # distance matrix
        self.token_distance = torch.tensor(
            cdist(self.atom_array.coord, self.atom_array.coord), dtype=torch.float64
        )

    def generate_spec_constraint(
        self,
        contact_specifics: list[tuple[int, int, float, float]],
        feature_type: str,
        shape: tuple[int, int, int],
    ) -> tuple[torch.Tensor, set[int]]:
        """
        parse constraint from user specification
        """

        contact_feature = torch.full(
            shape, fill_value=self.pad_value, dtype=torch.float32
        )
        for token_id_1, token_id_2, max_distance, min_distance in contact_specifics:
            contact_feature[token_id_1, token_id_2, 1] = max_distance
            contact_feature[token_id_2, token_id_1, 1] = max_distance
            contact_feature[token_id_1, token_id_2, 0] = min_distance
            contact_feature[token_id_2, token_id_1, 0] = min_distance

        contact_feature = self.encode(
            feature=contact_feature, feature_type=feature_type
        )
        constrained_tokens = set()
        for token_id_1, token_id_2, _, _ in contact_specifics:
            constrained_tokens.add(token_id_1)
            constrained_tokens.add(token_id_2)
        return contact_feature, constrained_tokens

    def get_real_distance(self, atom_idx_1: int, atom_idx_2: int) -> float:
        return self.token_distance[atom_idx_1, atom_idx_2]


class SubStructureFeaturizer(ConstraintFeaturizer):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        # Default distance bins
        self.distance_bins = torch.tensor([0, 4, 8, 16, torch.inf])

    def _add_coordinate_noise(self, coord_noise_scale: float = 0.05) -> torch.Tensor:
        """Add Gaussian noise to coordinates"""
        # Convert coordinates to tensor first
        coords = torch.tensor(self.atom_array.coord)
        if coord_noise_scale <= 0:
            return coords
        noisy_coords = (
            coords
            + torch.randn(
                coords.shape,
                generator=self.generator,
            )
            * coord_noise_scale
        )

        return noisy_coords

    def _get_distance_feature(
        self,
        selected_token_mask: torch.Tensor,
        coord_noise_scale: float = 0.05,
        feature_type: str = "one_hot",
    ) -> torch.Tensor:
        """
        Encode pairwise distances between selected tokens into binned one-hot features

        Parameters:
            - selected_token_mask: Boolean mask of selected tokens
            - coord_noise_scale: Scale of Gaussian noise to add to coordinates
            - feature_type: "one_hot" or "continuous"

        Returns:
            - distance_feature: For one_hot: [..., N_token, N_token] tensor with bin indices
                              For continuous: [..., N_token, N_token] tensor with distance values
        """
        n_tokens = len(self.asymid)

        # Initialize output tensor
        distance_feature = torch.full(
            (n_tokens, n_tokens),
            fill_value=-1 if feature_type == "one_hot" else self.pad_value,
            dtype=torch.long if feature_type == "one_hot" else torch.float32,
        )

        # Get selected token indices
        selected_tokens = torch.nonzero(selected_token_mask).squeeze(-1)
        if len(selected_tokens) <= 1:
            return distance_feature

        # Add noise to coordinates and calculate distances
        noisy_coords = self._add_coordinate_noise(coord_noise_scale)

        # Get token center atom indices
        token_centre_atom_indices = torch.tensor(
            self.token_array.get_annotation("centre_atom_index"), dtype=torch.long
        )
        selected_indices = token_centre_atom_indices[selected_tokens.long()]
        selected_coords = noisy_coords[selected_indices]

        # Calculate pairwise distances between selected tokens
        pairwise_distances = torch.cdist(selected_coords, selected_coords)

        if feature_type == "one_hot":
            # Digitize distances into bins
            binned_distances = (
                torch.bucketize(pairwise_distances, self.distance_bins, right=True) - 1
            )

            # Create mask for valid bins
            valid_bins = binned_distances > 0

            # Create indices for the full matrix
            rows = selected_tokens.repeat_interleave(len(selected_tokens))
            cols = selected_tokens.repeat(len(selected_tokens))

            # Get valid indices and their corresponding bin values
            valid_mask = valid_bins.flatten()
            valid_rows = rows[valid_mask]
            valid_cols = cols[valid_mask]
            valid_bins = binned_distances.flatten()[valid_mask]

            # Fill the distance feature matrix
            distance_feature[valid_rows, valid_cols] = valid_bins
            distance_feature[valid_cols, valid_rows] = valid_bins  # symmetric

        else:  # continuous
            # Create indices for the full matrix
            rows = selected_tokens.repeat_interleave(len(selected_tokens))
            cols = selected_tokens.repeat(len(selected_tokens))

            # Fill the distance feature matrix
            distance_feature[rows, cols] = pairwise_distances.flatten()
            distance_feature[cols, rows] = pairwise_distances.flatten()  # symmetric

        return distance_feature

    def get_valid_substructure_feature(
        self, mol_type_pairs: dict[str, float]
    ) -> torch.Tensor:
        """
        Find valid chains that form interfaces based on mol_type_pairs

        Parameters:
            - mol_type_pairs: Dict of type pairs (e.g. {'PP': threshold, 'LP': threshold})
                First type in pair is the one to be selected as substructure, threshold is the distance threshold to determine the interface
        Returns:
            - valid_asym_ids: List of valid asym_ids that can form interfaces
        """
        valid_asym_ids = []

        for type_pair, threshold in mol_type_pairs.items():
            # Parse type pair (e.g. 'PP' -> ('P','P'))
            query_type, key_type = type_pair

            # Get type masks
            valid_type_mask = (
                self.entity_type_dict[query_type][:, None]
                & self.entity_type_dict[key_type][None, :]
            )

            # Get different chain mask
            diff_chain_mask = self.asymid[:, None] != self.asymid[None, :]

            # Get distance mask
            dist_mask = (
                (self.token_distance <= threshold)
                & self.token_resolved_maskmat
                & diff_chain_mask
                & valid_type_mask
            )

            # Find chains that form interfaces
            for asym_id in torch.unique(self.asymid):
                # Only consider chains of query_type
                if not self.entity_type_dict[query_type][self.asymid == asym_id].any():
                    continue

                # Check if this chain forms interface with any chain of key_type
                cur_chain_mask = self.asymid == asym_id
                other_chains_mask = ~cur_chain_mask

                has_interface = (dist_mask[cur_chain_mask][:, other_chains_mask]).any()

                if has_interface:
                    valid_asym_ids.append(asym_id)

        return torch.tensor(valid_asym_ids)

    def _get_constraint_size(self, size: Any) -> int:
        """
        If size is not fixed, then we generate it randomly

        Args:
            - size: If >=1, randomly select chains; if <1, select size chains
        Returns:
            - size: Number of chains/proportion of tokens to select
        """
        if size < 1 and size > 0:
            # For size < 1, use as proportion directly
            return torch.rand(1).geometric_(size, generator=self.generator).int().item()
        elif size >= 1 and isinstance(size, (int, float)):
            # For size >= 1, return the size directly
            return int(size)
        else:
            raise NotImplementedError(f"Invalid size: {size}")

    def _sample_substructure(
        self,
        valid_asym_ids: torch.Tensor,
        size: Any,
        ratios: dict[str, list[float]],
        spec_asym_id: Any = None,
    ) -> torch.Tensor:
        """
        Sample substructure based on size and ratios

        Parameters:
            - valid_asym_ids: List of valid asym_ids
            - size: Total number of chains to select
            - ratios: Dict containing:
                - full: List of possible proportions for full chain selection
                - partial: Proportion of tokens to select for partial chains [0,1]
            - spec_asym_id: If provided, select from this specific chain
        """
        selected_token_mask = torch.zeros(len(self.asymid), dtype=torch.bool)

        if len(valid_asym_ids) == 0:
            return selected_token_mask

        # Handle spec_asym_id case
        if spec_asym_id is not None:
            if spec_asym_id not in valid_asym_ids:
                return selected_token_mask
            # Use partial_ratio for spec_asym_id
            chain_mask = (self.asymid == spec_asym_id) & self.token_resolved_mask
            chain_tokens = torch.nonzero(chain_mask).squeeze()
            if len(chain_tokens) == 0:  # Skip if no resolved tokens
                return selected_token_mask
            num_tokens = max(1, int(len(chain_tokens) * ratios["partial"]))
            selected_tokens = chain_tokens[
                torch.randperm(len(chain_tokens), generator=self.generator)[:num_tokens]
            ]
            selected_token_mask[selected_tokens] = True
            return selected_token_mask

        # Regular case: sample based on size and ratios
        if size == 0:
            return selected_token_mask

        # Randomly select full chain ratio from the list
        full_ratio_idx = torch.randint(
            len(ratios["full"]), (1,), generator=self.generator
        ).item()
        full_ratio = ratios["full"][full_ratio_idx]

        # Calculate number of chains for full and partial selection
        num_full_chains = min(int(size * full_ratio), len(valid_asym_ids))
        num_partial_chains = min(
            size - num_full_chains, len(valid_asym_ids) - num_full_chains
        )

        # Randomly shuffle and split valid_asym_ids
        shuffled_indices = torch.randperm(len(valid_asym_ids), generator=self.generator)
        full_chain_ids = valid_asym_ids[shuffled_indices[:num_full_chains]]
        partial_chain_ids = valid_asym_ids[
            shuffled_indices[num_full_chains : num_full_chains + num_partial_chains]
        ]

        # Select full chains
        for asym_id in full_chain_ids:
            chain_mask = (self.asymid == asym_id) & self.token_resolved_mask
            selected_token_mask |= chain_mask
        # Select partial chains
        for asym_id in partial_chain_ids:
            chain_mask = (self.asymid == asym_id) & self.token_resolved_mask
            chain_tokens = torch.nonzero(chain_mask).squeeze(-1)
            if len(chain_tokens) == 0:  # Skip if no resolved tokens
                continue
            num_tokens = max(
                1,
                int(len(chain_tokens) * torch.rand(1, generator=self.generator).item()),
            )
            selected_tokens = chain_tokens[
                torch.randperm(len(chain_tokens), generator=self.generator)[:num_tokens]
            ]
            selected_token_mask[selected_tokens] = True

        return selected_token_mask

    def generate(
        self,
        mol_type_pairs: dict[str, float],
        feature_type: str,
        size: Any,
        ratios: dict[str, list[float]],
        coord_noise_scale: float,
        spec_asym_id: int = None,
    ) -> tuple[torch.Tensor, set[int]]:
        """
        Generate substructure features

        Parameters:
            - mol_type_pairs: Dict of type pairs and their distance thresholds
            - feature_type: Type of feature encoding
            - size: Number of chains to select
            - ratios: Dict containing:
                - full: List of possible proportions for full chain selection
                - partial: Proportion of tokens to select for partial chains [0,1]
            - coord_noise_scale: Scale of Gaussian noise to add to coordinates
            - spec_asym_id: Specific chain to select from
        """
        constrained_tokens = set()
        if size == 0 or spec_asym_id == -1:
            distance_feature = torch.full(
                (self.asymid.shape[0], self.asymid.shape[0]),
                fill_value=-1 if feature_type == "one_hot" else self.pad_value,
                dtype=torch.long if feature_type == "one_hot" else torch.float32,
            )
            return (
                self.encode(
                    feature=distance_feature,
                    feature_type=feature_type,
                    num_classes=len(self.distance_bins) - 1,
                ),
                constrained_tokens,
            )

        self._get_generation_basics()

        # Get valid asym_ids that form interfaces
        valid_asym_ids = self.get_valid_substructure_feature(mol_type_pairs)

        if len(valid_asym_ids) == 0:
            distance_feature = torch.full(
                (self.asymid.shape[0], self.asymid.shape[0]),
                fill_value=-1 if feature_type == "one_hot" else self.pad_value,
                dtype=torch.long if feature_type == "one_hot" else torch.float32,
            )
            return (
                self.encode(
                    feature=distance_feature,
                    feature_type=feature_type,
                    num_classes=len(self.distance_bins) - 1,
                ),
                constrained_tokens,
            )

        size = self._get_constraint_size(size)

        # Sample tokens based on ratio and spec_asym_id
        selected_token_mask = self._sample_substructure(
            valid_asym_ids, size, ratios, spec_asym_id
        )

        # Get distance features (bin indices)
        distance_feature = self._get_distance_feature(
            selected_token_mask, coord_noise_scale, feature_type
        )

        # Track constrained tokens
        constrained_tokens = set(torch.nonzero(selected_token_mask).flatten().tolist())

        # Encode using base class method
        return (
            self.encode(
                feature=distance_feature,
                feature_type=feature_type,
                num_classes=len(self.distance_bins) - 1,
            ),
            constrained_tokens,
        )

    def analyze_features(self, feature_tensor: torch.Tensor) -> dict[str, Any]:
        """
        Analyze the features generated by the generate method
        """
        is_one_hot = len(feature_tensor.shape) == 3
        n_tokens = feature_tensor.shape[0]
        if is_one_hot:
            # For one-hot features
            # A token is active if it has any non-zero one-hot vector
            has_valid_distance = torch.any(feature_tensor, dim=-1)
            active_tokens = torch.any(has_valid_distance, dim=1)

            # Get distribution of distance bins (excluding zero vectors)
            valid_distances = feature_tensor[has_valid_distance]
            bin_counts = torch.sum(valid_distances, dim=0)  # Sum over all valid pairs
            distance_stats = {
                f"bin_{i}_count": count.item() for i, count in enumerate(bin_counts)
            }

        else:
            # For continuous features
            # A token is active if it has any non-zero distance to other tokens
            has_valid_distance = feature_tensor != 0
            active_tokens = torch.any(has_valid_distance, dim=1)

            # Get distribution of actual distances (excluding zeros)
            valid_distances = feature_tensor[has_valid_distance]
            if len(valid_distances) > 0:
                distance_stats = {
                    "mean": valid_distances.mean().item(),
                    "std": valid_distances.std().item(),
                    "min": valid_distances.min().item(),
                    "max": valid_distances.max().item(),
                }
            else:
                distance_stats = {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0}

        stats = {
            "num_active_tokens": active_tokens.sum().item(),
            "active_token_ratio": (active_tokens.sum().item() / n_tokens),
            "distance_distribution": distance_stats,
        }

        return stats

    def generate_spec_constraint(self, substructure_specifics, feature_type):
        """parse constraint from user specification

        Args:
            substructure_specifics (Dict): dictionary speicifing fixed tokens
                token_indices: List[int]
                token_coords: List[List[float]]
            feature_type: 'ont hot' by default
        """
        distance_feature_mat = torch.full(
            (self.asymid.shape[0], self.asymid.shape[0]),
            fill_value=-1 if feature_type == "one_hot" else self.pad_value,
            dtype=torch.long if feature_type == "one_hot" else torch.float32,
        )

        if len(substructure_specifics["token_indices"]) > 0:
            token_indices = torch.tensor(substructure_specifics["token_indices"])
            coords = torch.tensor(substructure_specifics["token_coords"])

            distance_mat = torch.cdist(coords, coords)
            distance_feature_mat[token_indices[:, None], token_indices[None, :]] = (
                distance_mat
            )

        distance_feature_mat = self.encode(
            distance_feature_mat, feature_type, num_classes=len(self.distance_bins) - 1
        )

        return distance_feature_mat
