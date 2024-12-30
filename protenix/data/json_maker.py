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


import argparse
import copy
import json
import os
from collections import defaultdict

import numpy as np
from biotite.structure import AtomArray, get_chain_starts, get_residue_starts

from protenix.data.constants import STD_RESIDUES
from protenix.data.filter import Filter
from protenix.data.parser import AddAtomArrayAnnot, MMCIFParser
from protenix.data.utils import get_lig_lig_bonds, get_ligand_polymer_bond_mask


def merge_covalent_bonds(
    covalent_bonds: list[dict], all_entity_counts: dict[str, int]
) -> list[dict]:
    """
    Merge covalent bonds with same entity and position.

    Args:
        covalent_bonds (list[dict]): A list of covalent bond dicts.
        all_entity_counts (dict[str, int]): A dict of entity id to chain count.

    Returns:
        list[dict]: A list of merged covalent bond dicts.
    """
    bonds_recorder = defaultdict(list)
    bonds_entity_counts = {}
    for bond_dict in covalent_bonds:
        bond_unique_string = []
        entity_counts = (
            all_entity_counts[str(bond_dict["entity1"])],
            all_entity_counts[str(bond_dict["entity2"])],
        )
        for i in range(2):
            for j in ["entity", "position", "atom"]:
                k = f"{j}{i+1}"
                bond_unique_string.append(str(bond_dict[k]))
        bond_unique_string = "_".join(bond_unique_string)
        bonds_recorder[bond_unique_string].append(bond_dict)
        bonds_entity_counts[bond_unique_string] = entity_counts

    merged_covalent_bonds = []
    for k, v in bonds_recorder.items():
        counts1 = bonds_entity_counts[k][0]
        counts2 = bonds_entity_counts[k][1]
        if counts1 == counts2 == len(v):
            bond_dict_copy = copy.deepcopy(v[0])
            del bond_dict_copy["copy1"]
            del bond_dict_copy["copy2"]
            merged_covalent_bonds.append(bond_dict_copy)
        else:
            merged_covalent_bonds.extend(v)
    return merged_covalent_bonds


def atom_array_to_input_json(
    atom_array: AtomArray,
    parser: MMCIFParser,
    assembly_id: str = None,
    output_json: str = None,
    sample_name=None,
    save_entity_and_asym_id=False,
) -> dict:
    """
    Convert a Biotite AtomArray to a dict that can be used as input to the model.

    Args:
        atom_array (AtomArray): Biotite Atom array.
        parser (MMCIFParser): Instantiated Protenix MMCIFParer.
        assembly_id (str, optional): Assembly ID. Defaults to None.
        output_json (str, optional): Output json file path. Defaults to None.
        sample_name (_type_, optional): The "name" filed in json file. Defaults to None.
        save_entity_and_asym_id (bool, optional): Whether to save entity and asym ids to json.
                                                  Defaults to False.

    Returns:
        dict: Protenix input json dict.
    """
    # get sequences after modified AtomArray
    entity_seq = parser.get_sequences(atom_array)

    # add unique chain id
    atom_array = AddAtomArrayAnnot.unique_chain_and_add_ids(atom_array)

    # get lig entity sequences and position
    label_entity_id_to_sequences = {}
    lig_chain_ids = []  # record chain_id of the first asym chain
    for label_entity_id in np.unique(atom_array.label_entity_id):
        if label_entity_id not in parser.entity_poly_type:
            current_lig_chain_ids = np.unique(
                atom_array.chain_id[atom_array.label_entity_id == label_entity_id]
            ).tolist()
            lig_chain_ids += current_lig_chain_ids
            for chain_id in current_lig_chain_ids:
                lig_atom_array = atom_array[atom_array.chain_id == chain_id]
                starts = get_residue_starts(lig_atom_array, add_exclusive_stop=True)
                seq = lig_atom_array.res_name[starts[:-1]].tolist()
                label_entity_id_to_sequences[label_entity_id] = seq

    # find polymer modifications
    entity_id_to_mod_list = {}
    for entity_id, res_names in parser.get_poly_res_names(atom_array).items():
        modifications_list = []
        for idx, res_name in enumerate(res_names):
            if res_name not in STD_RESIDUES:
                position = idx + 1
                modifications_list.append([position, f"CCD_{res_name}"])
        if modifications_list:
            entity_id_to_mod_list[entity_id] = modifications_list

    chain_starts = get_chain_starts(atom_array, add_exclusive_stop=False)
    chain_starts_atom_array = atom_array[chain_starts]

    json_dict = {
        "sequences": [],
    }
    if assembly_id is not None:
        json_dict["assembly_id"] = assembly_id

    unique_label_entity_id = np.unique(atom_array.label_entity_id)
    label_entity_id_to_entity_id_in_json = {}
    chain_id_to_copy_id_dict = {}
    for idx, label_entity_id in enumerate(unique_label_entity_id):
        entity_id_in_json = str(idx + 1)
        label_entity_id_to_entity_id_in_json[label_entity_id] = entity_id_in_json
        chain_ids_in_entity = chain_starts_atom_array.chain_id[
            chain_starts_atom_array.label_entity_id == label_entity_id
        ]
        for chain_count, chain_id in enumerate(chain_ids_in_entity):
            chain_id_to_copy_id_dict[chain_id] = chain_count + 1
    copy_id = np.vectorize(chain_id_to_copy_id_dict.get)(atom_array.chain_id)
    atom_array.set_annotation("copy_id", copy_id)

    all_entity_counts = {}
    skipped_entity_id = []
    for label_entity_id in unique_label_entity_id:
        entity_dict = {}
        asym_chains = chain_starts_atom_array[
            chain_starts_atom_array.label_entity_id == label_entity_id
        ]
        entity_type = parser.entity_poly_type.get(label_entity_id, "ligand")
        if entity_type != "ligand":
            if entity_type == "polypeptide(L)":
                entity_type = "proteinChain"
            elif entity_type == "polydeoxyribonucleotide":
                entity_type = "dnaSequence"
            elif entity_type == "polyribonucleotide":
                entity_type = "rnaSequence"
            else:
                # DNA/RNA hybrid, polypeptide(D), etc.
                skipped_entity_id.append(label_entity_id)
                continue

            sequence = entity_seq.get(label_entity_id)
            entity_dict["sequence"] = sequence
        else:
            # ligand
            lig_ccd = "_".join(label_entity_id_to_sequences[label_entity_id])
            entity_dict["ligand"] = f"CCD_{lig_ccd}"
        entity_dict["count"] = len(asym_chains)
        all_entity_counts[label_entity_id_to_entity_id_in_json[label_entity_id]] = len(
            asym_chains
        )
        if save_entity_and_asym_id:
            entity_dict["label_entity_id"] = str(label_entity_id)
            entity_dict["label_asym_id"] = asym_chains.label_asym_id.tolist()

        # add PTM info
        if label_entity_id in entity_id_to_mod_list:
            modifications = entity_id_to_mod_list[label_entity_id]
            if entity_type == "proteinChain":
                entity_dict["modifications"] = [
                    {"ptmPosition": position, "ptmType": mod_ccd_code}
                    for position, mod_ccd_code in modifications
                ]
            else:
                entity_dict["modifications"] = [
                    {"basePosition": position, "modificationType": mod_ccd_code}
                    for position, mod_ccd_code in modifications
                ]

        json_dict["sequences"].append({entity_type: entity_dict})

    # skip some uncommon entities
    atom_array = atom_array[~np.isin(atom_array.label_entity_id, skipped_entity_id)]

    # add covalent bonds
    atom_array = AddAtomArrayAnnot.add_token_mol_type(
        atom_array, parser.entity_poly_type
    )
    lig_polymer_bonds = get_ligand_polymer_bond_mask(atom_array, lig_include_ions=False)
    lig_lig_bonds = get_lig_lig_bonds(atom_array, lig_include_ions=False)
    inter_entity_bonds = np.vstack((lig_polymer_bonds, lig_lig_bonds))

    lig_indices = np.where(np.isin(atom_array.chain_id, lig_chain_ids))[0]
    lig_bond_mask = np.any(np.isin(inter_entity_bonds[:, :2], lig_indices), axis=1)
    inter_entity_bonds = inter_entity_bonds[lig_bond_mask]  # select bonds of ligands
    if inter_entity_bonds.size != 0:
        covalent_bonds = []
        for atoms in inter_entity_bonds[:, :2]:
            bond_dict = {}
            for i in range(2):
                atom = atom_array[atoms[i]]
                positon = atom.res_id
                bond_dict[f"entity{i+1}"] = int(
                    label_entity_id_to_entity_id_in_json[atom.label_entity_id]
                )
                bond_dict[f"position{i+1}"] = int(positon)
                bond_dict[f"atom{i+1}"] = atom.atom_name
                bond_dict[f"copy{i+1}"] = int(atom.copy_id)

            covalent_bonds.append(bond_dict)

        # merge covalent_bonds for same entity
        merged_covalent_bonds = merge_covalent_bonds(covalent_bonds, all_entity_counts)
        json_dict["covalent_bonds"] = merged_covalent_bonds

    json_dict["name"] = sample_name

    if output_json is not None:
        with open(output_json, "w") as f:
            json.dump([json_dict], f, indent=4)
    return json_dict


def cif_to_input_json(
    mmcif_file: str,
    assembly_id: str = None,
    altloc="first",
    output_json: str = None,
    sample_name=None,
    save_entity_and_asym_id=False,
) -> dict:
    """
    Convert mmcif file to Protenix input json file.

    Args:
        mmcif_file (str): mmCIF file path.
        assembly_id (str, optional): Assembly ID. Defaults to None.
        altloc (str, optional): Altloc selection. Defaults to "first".
        output_json (str, optional): Output json file path. Defaults to None.
        sample_name (_type_, optional): The "name" filed in json file. Defaults to None.
        save_entity_and_asym_id (bool, optional): Whether to save entity and asym ids to json.
                                                  Defaults to False.

    Returns:
        dict: Protenix input json dict.
    """
    parser = MMCIFParser(mmcif_file)
    atom_array = parser.get_structure(altloc, model=1, bond_lenth_threshold=None)

    # remove HOH from entities
    atom_array = Filter.remove_water(atom_array)
    atom_array = Filter.remove_hydrogens(atom_array)
    atom_array = parser.mse_to_met(atom_array)
    atom_array = Filter.remove_element_X(atom_array)

    # remove crystallization_aids
    if any(["DIFFRACTION" in m for m in parser.methods]):
        atom_array = Filter.remove_crystallization_aids(
            atom_array, parser.entity_poly_type
        )

    if assembly_id is not None:
        # expand created AtomArray by expand bioassembly
        atom_array = parser.expand_assembly(atom_array, assembly_id)

    if sample_name is None:
        sample_name = os.path.basename(mmcif_file).split(".")[0]

    json_dict = atom_array_to_input_json(
        atom_array,
        parser,
        assembly_id,
        output_json,
        sample_name,
        save_entity_and_asym_id=save_entity_and_asym_id,
    )
    return json_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cif_file", type=str, required=True, help="The cif file to parse"
    )
    parser.add_argument(
        "--json_file",
        type=str,
        required=False,
        default=None,
        help="The json file path to generate",
    )
    args = parser.parse_args()
    print(cif_to_input_json(args.cif_file, output_json=args.json_file))
