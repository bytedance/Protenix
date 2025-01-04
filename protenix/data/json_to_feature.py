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
import logging

import numpy as np
import torch
from biotite.structure import AtomArray

from protenix.data.featurizer import Featurizer
from protenix.data.json_parser import add_entity_atom_array, remove_leaving_atoms
from protenix.data.parser import AddAtomArrayAnnot
from protenix.data.tokenizer import AtomArrayTokenizer, TokenArray
from protenix.data.utils import int_to_letters

logger = logging.getLogger(__name__)


class SampleDictToFeatures:
    def __init__(self, single_sample_dict):
        self.single_sample_dict = single_sample_dict
        self.input_dict = add_entity_atom_array(single_sample_dict)
        self.entity_poly_type = self.get_entity_poly_type()

    def get_entity_poly_type(self) -> dict[str, str]:
        """
        Get the entity type for each entity.

        Allowed Value for "_entity_poly.type":
        · cyclic-pseudo-peptide
        · other
        · peptide nucleic acid
        · polydeoxyribonucleotide
        · polydeoxyribonucleotide/polyribonucleotide hybrid
        · polypeptide(D)
        · polypeptide(L)
        · polyribonucleotide

        Returns:
            dict[str, str]: a dict of polymer entity id to entity type.
        """
        entity_type_mapping_dict = {
            "proteinChain": "polypeptide(L)",
            "dnaSequence": "polydeoxyribonucleotide",
            "rnaSequence": "polyribonucleotide",
        }
        entity_poly_type = {}
        for idx, type2entity_dict in enumerate(self.input_dict["sequences"]):
            assert len(type2entity_dict) == 1, "Only one entity type is allowed."
            for entity_type, entity in type2entity_dict.items():
                if "sequence" in entity:
                    assert entity_type in [
                        "proteinChain",
                        "dnaSequence",
                        "rnaSequence",
                    ], 'The "sequences" field accepts only these entity types: ["proteinChain", "dnaSequence", "rnaSequence"].'
                    entity_poly_type[str(idx + 1)] = entity_type_mapping_dict[
                        entity_type
                    ]
        return entity_poly_type

    def build_full_atom_array(self) -> AtomArray:
        """
        By assembling the AtomArray of each entity, a complete AtomArray is created.

        Returns:
            AtomArray: Biotite Atom array.
        """
        atom_array = None
        asym_chain_idx = 0
        for idx, type2entity_dict in enumerate(self.input_dict["sequences"]):
            for entity_type, entity in type2entity_dict.items():
                entity_id = str(idx + 1)

                entity_atom_array = None
                for asym_chain_count in range(1, entity["count"] + 1):
                    asym_id_str = int_to_letters(asym_chain_idx + 1)
                    asym_chain = copy.deepcopy(entity["atom_array"])
                    chain_id = [asym_id_str] * len(asym_chain)
                    copy_id = [asym_chain_count] * len(asym_chain)
                    asym_chain.set_annotation("label_asym_id", chain_id)
                    asym_chain.set_annotation("auth_asym_id", chain_id)
                    asym_chain.set_annotation("chain_id", chain_id)
                    asym_chain.set_annotation("label_seq_id", asym_chain.res_id)
                    asym_chain.set_annotation("copy_id", copy_id)
                    if entity_atom_array is None:
                        entity_atom_array = asym_chain
                    else:
                        entity_atom_array += asym_chain
                    asym_chain_idx += 1

                entity_atom_array.set_annotation(
                    "label_entity_id", [entity_id] * len(entity_atom_array)
                )

                if entity_type in ["proteinChain", "dnaSequence", "rnaSequence"]:
                    entity_atom_array.hetero[:] = False
                else:
                    entity_atom_array.hetero[:] = True

                if atom_array is None:
                    atom_array = entity_atom_array
                else:
                    atom_array += entity_atom_array
        return atom_array

    @staticmethod
    def get_a_bond_atom(
        atom_array: AtomArray,
        entity_id: int,
        position: int,
        atom_name: str,
        copy_id: int = None,
    ) -> np.ndarray:
        """
        Get the atom index of a bond atom.

        Args:
            atom_array (AtomArray): Biotite Atom array.
            entity_id (int): Entity id.
            position (int): Residue index of the atom.
            atom_name (str): Atom name.
            copy_id (copy_id): A asym chain id in N copies of an entity.

        Returns:
            np.ndarray: Array of indices for specified atoms on each asym chain.
        """
        entity_mask = atom_array.label_entity_id == str(entity_id)
        position_mask = atom_array.res_id == int(position)
        atom_name_mask = atom_array.atom_name == str(atom_name)

        if copy_id is not None:
            copy_mask = atom_array.copy_id == int(copy_id)
            mask = entity_mask & position_mask & atom_name_mask & copy_mask
        else:
            mask = entity_mask & position_mask & atom_name_mask
        atom_indices = np.where(mask)[0]
        return atom_indices

    def add_bonds_between_entities(self, atom_array: AtomArray) -> AtomArray:
        """
        Based on the information in the "covalent_bonds",
        add a bond between specified atoms on each pair of asymmetric chains of the two entities.
        Note that this requires the number of asymmetric chains in both entities to be equal.

        Args:
            atom_array (AtomArray): Biotite Atom array.

        Returns:
            AtomArray: Biotite Atom array with bonds added.
        """
        if "covalent_bonds" not in self.input_dict:
            return atom_array

        bond_count = {}
        for bond_info_dict in self.input_dict["covalent_bonds"]:
            bond_atoms = []
            for idx, i in enumerate(["left", "right"]):
                entity_id = int(
                    bond_info_dict.get(
                        f"{i}_entity", bond_info_dict.get(f"entity{idx+1}")
                    )
                )
                copy_id = bond_info_dict.get(
                    f"{i}_copy", bond_info_dict.get(f"copy{idx+1}")
                )
                position = int(
                    bond_info_dict.get(
                        f"{i}_position", bond_info_dict.get(f"position{idx+1}")
                    )
                )
                atom_name = bond_info_dict.get(
                    f"{i}_atom", bond_info_dict.get(f"atom{idx+1}")
                )

                if copy_id is not None:
                    copy_id = int(copy_id)

                if isinstance(atom_name, str):
                    if atom_name.isdigit():
                        # Convert SMILES atom index to int
                        atom_name = int(atom_name)

                if isinstance(atom_name, int):
                    # Convert AtomMap in SMILES to atom name in AtomArray
                    entity_dict = self.input_dict["sequences"][
                        int(entity_id - 1)
                    ].values()
                    assert "atom_map_to_atom_name" in entity_dict
                    atom_name = entity_dict["atom_map_to_atom_name"][atom_name]

                # Get bond atoms by entity_id, position, atom_name
                atom_indices = self.get_a_bond_atom(
                    atom_array, entity_id, position, atom_name, copy_id
                )
                assert (
                    atom_indices.size > 0
                ), f"No atom found for {atom_name} in entity {entity_id} at position {position}."
                bond_atoms.append(atom_indices)
            assert len(bond_atoms[0]) == len(bond_atoms[1]), (
                'Can not create bonds because the "count" of entity1 '
                f'({bond_info_dict.get("left_entity", bond_info_dict.get("entity1"))}) '
                f'and entity2 ({bond_info_dict.get("right_entity", bond_info_dict.get("entity2"))}) are not equal. '
            )

            # Create bond between each asym chain pair
            for atom_idx1, atom_idx2 in zip(bond_atoms[0], bond_atoms[1]):
                atom_array.bonds.add_bond(atom_idx1, atom_idx2, 1)
                bond_count[atom_idx1] = bond_count.get(atom_idx1, 0) + 1
                bond_count[atom_idx2] = bond_count.get(atom_idx2, 0) + 1

        atom_array = remove_leaving_atoms(atom_array, bond_count)

        return atom_array

    @staticmethod
    def add_atom_array_attributes(
        atom_array: AtomArray, entity_poly_type: dict[str, str]
    ) -> AtomArray:
        """
        Add attributes to the Biotite AtomArray.

        Args:
            atom_array (AtomArray): Biotite Atom array.
            entity_poly_type (dict[str, str]): a dict of polymer entity id to entity type.

        Returns:
            AtomArray: Biotite Atom array with attributes added.
        """
        atom_array = AddAtomArrayAnnot.add_token_mol_type(atom_array, entity_poly_type)
        atom_array = AddAtomArrayAnnot.add_centre_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_atom_mol_type_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_distogram_rep_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_plddt_m_rep_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_cano_seq_resname(atom_array)
        atom_array = AddAtomArrayAnnot.add_tokatom_idx(atom_array)
        atom_array = AddAtomArrayAnnot.add_modified_res_mask(atom_array)
        atom_array = AddAtomArrayAnnot.unique_chain_and_add_ids(atom_array)
        atom_array = AddAtomArrayAnnot.find_equiv_mol_and_assign_ids(
            atom_array, check_final_equiv=False
        )
        atom_array = AddAtomArrayAnnot.add_ref_space_uid(atom_array)
        return atom_array

    @staticmethod
    def mse_to_met(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI chapter 2.1
        MSE residues are converted to MET residues.

        Args:
            atom_array (AtomArray): Biotite AtomArray object.

        Returns:
            AtomArray: Biotite AtomArray object after converted MSE to MET.
        """
        mse = atom_array.res_name == "MSE"
        se = mse & (atom_array.atom_name == "SE")
        atom_array.atom_name[se] = "SD"
        atom_array.element[se] = "S"
        atom_array.res_name[mse] = "MET"
        atom_array.hetero[mse] = False
        return atom_array

    def get_atom_array(self) -> AtomArray:
        """
        Create a Biotite AtomArray and add attributes from the input dict.

        Returns:
            AtomArray: Biotite Atom array.
        """
        atom_array = self.build_full_atom_array()
        atom_array = self.add_bonds_between_entities(atom_array)
        atom_array = self.mse_to_met(atom_array)
        atom_array = self.add_atom_array_attributes(atom_array, self.entity_poly_type)
        return atom_array

    def get_feature_dict(self) -> tuple[dict[str, torch.Tensor], AtomArray, TokenArray]:
        """
        Generates a feature dictionary from the input sample dictionary.

        Returns:
            A tuple containing:
                - A dictionary of features.
                - An AtomArray object.
                - A TokenArray object.
        """
        atom_array = self.get_atom_array()

        aa_tokenizer = AtomArrayTokenizer(atom_array)
        token_array = aa_tokenizer.get_token_array()

        featurizer = Featurizer(token_array, atom_array)
        feature_dict = featurizer.get_all_input_features()

        token_array_with_frame = featurizer.get_token_frame(
            token_array=token_array,
            atom_array=atom_array,
            ref_pos=feature_dict["ref_pos"],
            ref_mask=feature_dict["ref_mask"],
        )

        # [N_token]
        feature_dict["has_frame"] = torch.Tensor(
            token_array_with_frame.get_annotation("has_frame")
        ).long()

        # [N_token, 3]
        feature_dict["frame_atom_index"] = torch.Tensor(
            token_array_with_frame.get_annotation("frame_atom_index")
        ).long()
        return feature_dict, atom_array, token_array
