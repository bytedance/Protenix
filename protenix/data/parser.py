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
import functools
import gzip
import logging
import random
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Optional, Union

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
from biotite.structure import AtomArray, get_chain_starts, get_residue_starts
from biotite.structure.io.pdbx import convert as pdbx_convert
from biotite.structure.molecules import get_molecule_indices

from protenix.data import ccd
from protenix.data.ccd import get_ccd_ref_info
from protenix.data.constants import (
    CRYSTALLIZATION_METHODS,
    DNA_STD_RESIDUES,
    GLYCANS,
    LIGAND_EXCLUSION,
    PRO_STD_RESIDUES,
    PROT_STD_RESIDUES_ONE_TO_THREE,
    RES_ATOMS_DICT,
    RNA_STD_RESIDUES,
    STD_RESIDUES,
)
from protenix.data.filter import Filter
from protenix.data.utils import (
    atom_select,
    get_inter_residue_bonds,
    get_ligand_polymer_bond_mask,
    get_starts_by,
    parse_pdb_cluster_file_to_dict,
)

logger = logging.getLogger(__name__)

# Ignore inter residue metal coordinate bonds in mmcif _struct_conn
if "metalc" in pdbx_convert.PDBX_COVALENT_TYPES:  # for reload
    pdbx_convert.PDBX_COVALENT_TYPES.remove("metalc")


class MMCIFParser:
    """
    Parsing and extracting information from mmCIF files.
    """

    def __init__(self, mmcif_file: Union[str, Path]):
        self.cif = self._parse(mmcif_file=mmcif_file)

    def _parse(self, mmcif_file: Union[str, Path]) -> pdbx.CIFFile:
        mmcif_file = Path(mmcif_file)
        if mmcif_file.suffix == ".gz":
            with gzip.open(mmcif_file, "rt") as f:
                cif_file = pdbx.CIFFile.read(f)
        else:
            with open(mmcif_file, "rt") as f:
                cif_file = pdbx.CIFFile.read(f)
        return cif_file

    def get_category_table(self, name: str) -> Union[pd.DataFrame, None]:
        """
        Retrieve a category table from the CIF block and return it as a pandas DataFrame.

        Args:
            name (str): The name of the category to retrieve from the CIF block.

        Returns:
            Union[pd.DataFrame, None]: A pandas DataFrame containing the category data if the category exists,
                                       otherwise None.
        """
        if name not in self.cif.block:
            return None
        category = self.cif.block[name]
        category_dict = {k: column.as_array() for k, column in category.items()}
        return pd.DataFrame(category_dict, dtype=str)

    @functools.cached_property
    def pdb_id(self) -> str:
        """
        Extracts and returns the PDB ID from the CIF block.

        Returns:
            str: The PDB ID in lowercase if present, otherwise an empty string.
        """

        if "entry" not in self.cif.block:
            return ""
        else:
            return self.cif.block["entry"]["id"].as_item().lower()

    def num_assembly_polymer_chains(self, assembly_id: str = "1") -> int:
        """
        Calculate the number of polymer chains in a specified assembly.

        Args:
            assembly_id (str): The ID of the assembly to count polymer chains for.
                               Defaults to "1". If "all", counts chains for all assemblies.

        Returns:
            int: The total number of polymer chains in the specified assembly.
                 If the oligomeric count is invalid (e.g., '?'), the function returns None.
        """
        chain_count = 0
        for _assembly_id, _chain_count in zip(
            self.cif.block["pdbx_struct_assembly"]["id"].as_array(),
            self.cif.block["pdbx_struct_assembly"]["oligomeric_count"].as_array(),
        ):
            if assembly_id == "all" or _assembly_id == assembly_id:
                try:
                    chain_count += int(_chain_count)
                except ValueError:
                    # oligomeric_count == '?'.  e.g. 1hya.cif
                    return
        return chain_count

    @functools.cached_property
    def resolution(self) -> float:
        """
        Get resolution for X-ray and cryoEM.
        Some methods don't have resolution, set as -1.0

        Returns:
            float: resolution (set to -1.0 if not found)
        """
        block = self.cif.block
        resolution_names = [
            "refine.ls_d_res_high",
            "em_3d_reconstruction.resolution",
            "reflns.d_resolution_high",
        ]
        for category_item in resolution_names:
            category, item = category_item.split(".")
            if category in block and item in block[category]:
                try:
                    resolution = block[category][item].as_array(float)[0]
                    # "." will be converted to 0.0, but it is not a valid resolution.
                    if resolution == 0.0:
                        continue
                    return resolution
                except ValueError:
                    # in some cases, resolution_str is "?"
                    continue
        return -1.0

    @functools.cached_property
    def release_date(self) -> str:
        """
        Get first release date.

        Returns:
            str: yyyy-mm-dd
        """

        def _is_valid_date_format(date_string):
            try:
                datetime.strptime(date_string, "%Y-%m-%d")
                return True
            except ValueError:
                return False

        if "pdbx_audit_revision_history" in self.cif.block:
            history = self.cif.block["pdbx_audit_revision_history"]
            # np.str_ is inherit from str, so return is str
            date = history["revision_date"].as_array()[0]
        else:
            # no release date
            date = "9999-12-31"

        valid_date = _is_valid_date_format(date)
        assert (
            valid_date
        ), f"Invalid date format: {date}, it should be yyyy-mm-dd format"
        return date

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

    @staticmethod
    def fix_arginine(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI chapter 2.1
        Arginine naming ambiguities are fixed (ensuring NH1 is always closer to CD than NH2).

        Args:
            atom_array (AtomArray): Biotite AtomArray object.

        Returns:
            AtomArray: Biotite AtomArray object after fix arginine .
        """

        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start_i, stop_i in zip(starts[:-1], starts[1:]):
            if atom_array[start_i].res_name != "ARG":
                continue
            cd_idx, nh1_idx, nh2_idx = None, None, None
            for idx in range(start_i, stop_i):
                if atom_array.atom_name[idx] == "CD":
                    cd_idx = idx
                if atom_array.atom_name[idx] == "NH1":
                    nh1_idx = idx
                if atom_array.atom_name[idx] == "NH2":
                    nh2_idx = idx
            if cd_idx and nh1_idx and nh2_idx:  # all not None
                cd_nh1 = atom_array.coord[nh1_idx] - atom_array.coord[cd_idx]
                d2_cd_nh1 = np.sum(cd_nh1**2)
                cd_nh2 = atom_array.coord[nh2_idx] - atom_array.coord[cd_idx]
                d2_cd_nh2 = np.sum(cd_nh2**2)
                if d2_cd_nh2 < d2_cd_nh1:
                    atom_array.coord[[nh1_idx, nh2_idx]] = atom_array.coord[
                        [nh2_idx, nh1_idx]
                    ]
        return atom_array

    @functools.cached_property
    def methods(self) -> list[str]:
        """the methods to get the structure

        most of the time, methods only has one method, such as 'X-RAY DIFFRACTION',
        but about 233 entries have multi methods, such as ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'].

        Allowed Values:
        https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl.method.html

        Returns:
            list[str]: such as ['X-RAY DIFFRACTION'], ['ELECTRON MICROSCOPY'], ['SOLUTION NMR', 'THEORETICAL MODEL'],
                ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'], ['ELECTRON MICROSCOPY', 'SOLUTION NMR'], etc.
        """
        if "exptl" not in self.cif.block:
            return []
        else:
            methods = self.cif.block["exptl"]["method"]
            return methods.as_array()

    def get_poly_res_names(
        self, atom_array: Optional[AtomArray] = None
    ) -> dict[str, list[str]]:
        """get 3-letter residue names by combining mmcif._entity_poly_seq and atom_array

        if ref_atom_array is None: keep first altloc residue of the same res_id based in mmcif._entity_poly_seq
        if ref_atom_array is provided: keep same residue of ref_atom_array.

        Returns
            dict[str, list[str]]: label_entity_id --> [res_ids, res_names]
        """
        entity_res_names = {}
        if atom_array is not None:
            # build entity_id -> res_id -> res_name for input atom array
            res_starts = struc.get_residue_starts(atom_array, add_exclusive_stop=False)
            for start in res_starts:
                entity_id = atom_array.label_entity_id[start]
                res_id = atom_array.res_id[start]
                res_name = atom_array.res_name[start]
                if entity_id in entity_res_names:
                    entity_res_names[entity_id][res_id] = res_name
                else:
                    entity_res_names[entity_id] = {res_id: res_name}

        # build reference entity atom array, including missing residues
        entity_poly_seq = self.get_category_table("entity_poly_seq")
        if entity_poly_seq is None:
            return {}

        poly_res_names = {}
        for entity_id, poly_type in self.entity_poly_type.items():
            chain_mask = entity_poly_seq.entity_id == entity_id
            seq_mon_ids = entity_poly_seq.mon_id[chain_mask].to_numpy(dtype=str)

            # replace all MSE to MET in _entity_poly_seq.mon_id
            seq_mon_ids[seq_mon_ids == "MSE"] = "MET"

            seq_nums = entity_poly_seq.num[chain_mask].to_numpy(dtype=int)

            if np.unique(seq_nums).size == seq_nums.size:
                # no altloc residues
                poly_res_names[entity_id] = seq_mon_ids
                continue

            # filter altloc residues, eg: 181 ALA (altloc A); 181 GLY (altloc B)
            select_mask = np.zeros(len(seq_nums), dtype=bool)
            matching_res_id = seq_nums[0]
            for i, res_id in enumerate(seq_nums):
                if res_id != matching_res_id:
                    continue

                res_name_in_atom_array = entity_res_names.get(entity_id, {}).get(res_id)
                if res_name_in_atom_array is None:
                    # res_name is mssing in atom_array,
                    # keep first altloc residue of the same res_id
                    select_mask[i] = True
                else:
                    # keep match residue to atom_array
                    if res_name_in_atom_array == seq_mon_ids[i]:
                        select_mask[i] = True

                if select_mask[i]:
                    matching_res_id += 1

            seq_mon_ids = seq_mon_ids[select_mask]
            seq_nums = seq_nums[select_mask]
            assert len(seq_nums) == max(seq_nums)
            poly_res_names[entity_id] = seq_mon_ids
        return poly_res_names

    def get_sequences(self, atom_array=None) -> dict:
        """get sequence by combining mmcif._entity_poly_seq and atom_array

        if ref_atom_array is None: keep first altloc residue of the same res_id based in mmcif._entity_poly_seq
        if ref_atom_array is provided: keep same residue of atom_array.

        Return
            Dict{str:str}: label_entity_id --> canonical_sequence
        """
        sequences = {}
        for entity_id, res_names in self.get_poly_res_names(atom_array).items():
            seq = ccd.res_names_to_sequence(res_names)
            sequences[entity_id] = seq
        return sequences

    @functools.cached_property
    def entity_poly_type(self) -> dict[str, str]:
        """
        Ref: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity_poly.type.html
        Map entity_id to entity_poly_type.

        Allowed Value:
        · cyclic-pseudo-peptide
        · other
        · peptide nucleic acid
        · polydeoxyribonucleotide
        · polydeoxyribonucleotide/polyribonucleotide hybrid
        · polypeptide(D)
        · polypeptide(L)
        · polyribonucleotide

        Returns:
            Dict: a dict of label_entity_id --> entity_poly_type.
        """
        entity_poly = self.get_category_table("entity_poly")
        if entity_poly is None:
            return {}

        return {i: t for i, t in zip(entity_poly.entity_id, entity_poly.type)}

    def filter_altloc(self, atom_array: AtomArray, altloc: str = "first") -> AtomArray:
        """
        Filter alternate conformations (altloc) of a given AtomArray based on the specified criteria.
        For example, in 2PXS, there are two res_name (XYG|DYG) at res_id 63.

        Args:
            atom_array : AtomArray
                The array of atoms to filter.
            altloc : str, optional
                The criteria for filtering alternate conformations. Possible values are:
                - "first": Keep the first alternate conformation.
                - "all": Keep all alternate conformations.
                - "A", "B", etc.: Keep the specified alternate conformation.
                - "global_largest": Keep the alternate conformation with the largest average occupancy.

        Returns:
            AtomArray
                The filtered AtomArray based on the specified altloc criteria.
        """
        if altloc == "all":
            return atom_array

        elif altloc == "first":
            letter_altloc_ids = np.unique(atom_array.label_alt_id)
            if len(letter_altloc_ids) == 1 and letter_altloc_ids[0] == ".":
                return atom_array
            letter_altloc_ids = letter_altloc_ids[letter_altloc_ids != "."]
            altloc_id = np.sort(letter_altloc_ids)[0]
            return atom_array[np.isin(atom_array.label_alt_id, [altloc_id, "."])]

        elif altloc == "global_largest":
            occ_dict = defaultdict(list)
            res_altloc = defaultdict(list)

            res_starts = get_residue_starts(atom_array, add_exclusive_stop=True)
            for res_start, _res_end in zip(res_starts[:-1], res_starts[1:]):
                altloc_char = atom_array.label_alt_id[res_start]
                if altloc_char == ".":
                    continue

                occupency = atom_array.occupancy[res_start]
                occ_dict[altloc_char].append(occupency)

                chain_id = atom_array.chain_id[res_start]
                res_id = atom_array.res_id[res_start]
                res_altloc[(chain_id, res_id)].append(altloc_char)

            alt_and_avg_occ = [
                (altloc_char, np.mean(occ_list))
                for altloc_char, occ_list in occ_dict.items()
            ]
            sorted_altloc_chars = [
                i[0] for i in sorted(alt_and_avg_occ, key=lambda x: x[1], reverse=True)
            ]

            selected_mask = np.zeros(len(atom_array), dtype=bool)
            for res_start, res_end in zip(res_starts[:-1], res_starts[1:]):
                chain_id = atom_array.chain_id[res_start]
                res_id = atom_array.res_id[res_start]
                altloc_char = atom_array.label_alt_id[res_start]

                if altloc_char == ".":
                    selected_mask[res_start:res_end] = True
                else:
                    res_sorted_altloc = [
                        i
                        for i in sorted_altloc_chars
                        if i in res_altloc[(chain_id, res_id)]
                    ]
                    selected_altloc = res_sorted_altloc[0]
                    if altloc_char == selected_altloc:
                        selected_mask[res_start:res_end] = True
            return atom_array[selected_mask]

        else:
            return atom_array[np.isin(atom_array.label_alt_id, [altloc, "."])]

    @staticmethod
    def replace_auth_with_label(atom_array: AtomArray) -> AtomArray:
        """
        Replace the author-provided chain ID with the label asym ID in the given AtomArray.

        This function addresses the issue described in https://github.com/biotite-dev/biotite/issues/553.
        It updates the `chain_id` of the `atom_array` to match the `label_asym_id` and resets the ligand
        residue IDs (`res_id`) for chains where the `label_seq_id` is ".". The residue IDs are reset
        sequentially starting from 1 within each chain.

        Args:
            atom_array (AtomArray): The input AtomArray object to be modified.

        Returns:
            AtomArray: The modified AtomArray with updated chain IDs and residue IDs.
        """
        atom_array.chain_id = atom_array.label_asym_id

        # reset ligand res_id
        res_id = copy.deepcopy(atom_array.label_seq_id)
        chain_starts = get_chain_starts(atom_array, add_exclusive_stop=True)
        for chain_start, chain_stop in zip(chain_starts[:-1], chain_starts[1:]):
            if atom_array.label_seq_id[chain_start] != ".":
                continue
            else:
                res_starts = get_residue_starts(
                    atom_array[chain_start:chain_stop], add_exclusive_stop=True
                )
                num = 1
                for res_start, res_stop in zip(res_starts[:-1], res_starts[1:]):
                    res_id[chain_start:chain_stop][res_start:res_stop] = num
                    num += 1

        atom_array.res_id = res_id.astype(int)
        return atom_array

    def get_structure(
        self,
        altloc: str = "first",
        model: int = 1,
        bond_lenth_threshold: Union[float, None] = 2.4,
    ) -> AtomArray:
        """
        Get an AtomArray created by bioassembly of MMCIF.

        altloc: "first", "all", "A", "B", etc
        model: the model number of the structure.
        bond_lenth_threshold: the threshold of bond length. If None, no filter will be applied.
                              Default is 2.4 Angstroms.

        Returns:
            AtomArray: Biotite AtomArray object created by bioassembly of MMCIF.
        """
        use_author_fields = True
        extra_fields = ["label_asym_id", "label_entity_id", "auth_asym_id"]  # chain
        extra_fields += ["label_seq_id", "auth_seq_id"]  # residue
        atom_site_fields = {
            "occupancy": "occupancy",
            "pdbx_formal_charge": "charge",
            "B_iso_or_equiv": "b_factor",
            "label_alt_id": "label_alt_id",
        }  # atom
        for atom_site_name, alt_name in atom_site_fields.items():
            if atom_site_name in self.cif.block["atom_site"]:
                extra_fields.append(alt_name)

        block = self.cif.block

        extra_fields = set(extra_fields)

        atom_site = block.get("atom_site")

        models = atom_site["pdbx_PDB_model_num"].as_array(np.int32)
        model_starts = pdbx_convert._get_model_starts(models)
        model_count = len(model_starts)

        if model == 0:
            raise ValueError("The model index must not be 0")
        # Negative models mean model indexing starting from last model

        model = model_count + model + 1 if model < 0 else model
        if model > model_count:
            raise ValueError(
                f"The file has {model_count} models, "
                f"the given model {model} does not exist"
            )

        model_atom_site = pdbx_convert._filter_model(atom_site, model_starts, model)
        # Any field of the category would work here to get the length
        model_length = model_atom_site.row_count
        atoms = AtomArray(model_length)

        atoms.coord[:, 0] = model_atom_site["Cartn_x"].as_array(np.float32)
        atoms.coord[:, 1] = model_atom_site["Cartn_y"].as_array(np.float32)
        atoms.coord[:, 2] = model_atom_site["Cartn_z"].as_array(np.float32)

        atoms.box = pdbx_convert._get_box(block)

        # The below part is the same for both, AtomArray and AtomArrayStack
        pdbx_convert._fill_annotations(
            atoms, model_atom_site, extra_fields, use_author_fields
        )

        bonds = struc.connect_via_residue_names(atoms, inter_residue=False)
        if "struct_conn" in block:
            conn_bonds = pdbx_convert._parse_inter_residue_bonds(
                model_atom_site, block["struct_conn"]
            )
            coord1 = atoms.coord[conn_bonds._bonds[:, 0]]
            coord2 = atoms.coord[conn_bonds._bonds[:, 1]]
            dist = np.linalg.norm(coord1 - coord2, axis=1)
            if bond_lenth_threshold is not None:
                conn_bonds._bonds = conn_bonds._bonds[dist < bond_lenth_threshold]
            bonds = bonds.merge(conn_bonds)
        atoms.bonds = bonds

        atom_array = self.filter_altloc(atoms, altloc=altloc)

        # inference inter residue bonds based on res_id (auth_seq_id) and label_asym_id.
        atom_array = ccd.add_inter_residue_bonds(
            atom_array,
            exclude_struct_conn_pairs=True,
            remove_far_inter_chain_pairs=True,
        )

        # use label_seq_id to match seq and structure
        atom_array = self.replace_auth_with_label(atom_array)

        # inference inter residue bonds based on new res_id (label_seq_id).
        # the auth_seq_id is not reliable, some are discontinuous (8bvh), some with insertion codes (6ydy).
        atom_array = ccd.add_inter_residue_bonds(
            atom_array, exclude_struct_conn_pairs=True
        )
        return atom_array

    def expand_assembly(
        self, structure: AtomArray, assembly_id: str = "1"
    ) -> AtomArray:
        """
        Expand the given assembly to all chains
        copy from biotite.structure.io.pdbx.get_assembly

        Args:
            structure (AtomArray): The AtomArray of the structure to expand.
            assembly_id (str, optional): The assembly ID in mmCIF file. Defaults to "1".
                                         If assembly_id is "all", all assemblies will be returned.

        Returns:
            AtomArray: The assembly AtomArray.
        """
        block = self.cif.block

        try:
            assembly_gen_category = block["pdbx_struct_assembly_gen"]
        except KeyError:
            logging.info(
                "File has no 'pdbx_struct_assembly_gen' category, return original structure."
            )
            return structure

        try:
            struct_oper_category = block["pdbx_struct_oper_list"]
        except KeyError:
            logging.info(
                "File has no 'pdbx_struct_oper_list' category, return original structure."
            )
            return structure

        assembly_ids = assembly_gen_category["assembly_id"].as_array(str)

        if assembly_id != "all":
            if assembly_id is None:
                assembly_id = assembly_ids[0]
            elif assembly_id not in assembly_ids:
                raise KeyError(f"File has no Assembly ID '{assembly_id}'")

        ### Calculate all possible transformations
        transformations = pdbx_convert._get_transformations(struct_oper_category)

        ### Get transformations and apply them to the affected asym IDs
        assembly = None
        assembly_1_mask = []
        for id, op_expr, asym_id_expr in zip(
            assembly_gen_category["assembly_id"].as_array(str),
            assembly_gen_category["oper_expression"].as_array(str),
            assembly_gen_category["asym_id_list"].as_array(str),
        ):
            # Find the operation expressions for given assembly ID
            # We already asserted that the ID is actually present
            if assembly_id == "all" or id == assembly_id:
                operations = pdbx_convert._parse_operation_expression(op_expr)
                asym_ids = asym_id_expr.split(",")
                # Filter affected asym IDs
                sub_structure = copy.deepcopy(
                    structure[..., np.isin(structure.label_asym_id, asym_ids)]
                )
                sub_assembly = pdbx_convert._apply_transformations(
                    sub_structure, transformations, operations
                )
                # Merge the chains with asym IDs for this operation
                # with chains from other operations
                if assembly is None:
                    assembly = sub_assembly
                else:
                    assembly += sub_assembly

                if id == "1":
                    assembly_1_mask.extend([True] * len(sub_assembly))
                else:
                    assembly_1_mask.extend([False] * len(sub_assembly))

        if assembly_id == "1" or assembly_id == "all":
            assembly.set_annotation("assembly_1", np.array(assembly_1_mask))
        return assembly

    def _get_core_indices(self, atom_array):
        if "assembly_1" in atom_array._annot:
            core_indices = np.where(atom_array.assembly_1)[0]
        else:
            core_indices = None
        return core_indices

    def get_bioassembly(
        self,
        assembly_id: str = "1",
        max_assembly_chains: int = 1000,
    ) -> dict[str, Any]:
        """
        Build the given biological assembly.

        Args:
            assembly_id (str, optional): Assembly ID. Defaults to "1".
            max_assembly_chains (int, optional): Max allowed chains in the assembly. Defaults to 1000.

        Returns:
            dict[str, Any]: A dictionary containing basic Bioassembly information, including:
                - "pdb_id": The PDB ID.
                - "sequences": The sequences associated with the assembly.
                - "release_date": The release date of the structure.
                - "assembly_id": The assembly ID.
                - "num_assembly_polymer_chains": The number of polymer chains in the assembly.
                - "num_prot_chains": The number of protein chains in the assembly.
                - "entity_poly_type": The type of polymer entities.
                - "resolution": The resolution of the structure. Set to -1.0 if resolution not found.
                - "atom_array": The AtomArray object representing the structure.
                - "num_tokens": The number of tokens in the AtomArray.
        """
        num_assembly_polymer_chains = self.num_assembly_polymer_chains(assembly_id)
        bioassembly_dict = {
            "pdb_id": self.pdb_id,
            "sequences": self.get_sequences(),  # label_entity_id --> canonical_sequence
            "release_date": self.release_date,
            "assembly_id": assembly_id,
            "num_assembly_polymer_chains": num_assembly_polymer_chains,
            "num_prot_chains": -1,
            "entity_poly_type": self.entity_poly_type,
            "resolution": self.resolution,
            "atom_array": None,
        }
        if (not num_assembly_polymer_chains) or (
            num_assembly_polymer_chains > max_assembly_chains
        ):
            return bioassembly_dict

        # created AtomArray of first model from mmcif atom_site (Asymmetric Unit)
        atom_array = self.get_structure()

        # update sequences: keep same altloc residue with atom_array
        bioassembly_dict["sequences"] = self.get_sequences(atom_array)

        pipeline_functions = [
            Filter.remove_water,
            Filter.remove_hydrogens,
            lambda aa: Filter.remove_polymer_chains_all_residues_unknown(
                aa, self.entity_poly_type
            ),
            # Note: Filter.remove_polymer_chains_too_short not being used
            lambda aa: Filter.remove_polymer_chains_with_consecutive_c_alpha_too_far_away(
                aa, self.entity_poly_type
            ),
            self.fix_arginine,
            self.add_missing_atoms_and_residues,  # and add annotation is_resolved (False for missing atoms)
            self.mse_to_met,  # do mse_to_met() after add_missing_atoms_and_residues()
            Filter.remove_element_X,  # remove X element (including ASX->ASP, GLX->GLU) after add_missing_atoms_and_residues()
        ]

        if set(self.methods) & CRYSTALLIZATION_METHODS:
            # AF3 SI 2.5.4 Crystallization aids are removed if the mmCIF method information indicates that crystallography was used.
            pipeline_functions.append(
                lambda aa: Filter.remove_crystallization_aids(aa, self.entity_poly_type)
            )

        for func in pipeline_functions:
            atom_array = func(atom_array)
            if len(atom_array) == 0:
                # no atoms left
                return bioassembly_dict

        atom_array = AddAtomArrayAnnot.add_token_mol_type(
            atom_array, self.entity_poly_type
        )
        atom_array = AddAtomArrayAnnot.add_centre_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_atom_mol_type_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_distogram_rep_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_plddt_m_rep_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_cano_seq_resname(atom_array)
        atom_array = AddAtomArrayAnnot.add_tokatom_idx(atom_array)
        atom_array = AddAtomArrayAnnot.add_modified_res_mask(atom_array)
        assert (
            atom_array.centre_atom_mask.sum()
            == atom_array.distogram_rep_atom_mask.sum()
        )

        # expand created AtomArray by expand bioassembly
        atom_array = self.expand_assembly(atom_array, assembly_id)

        if len(atom_array) == 0:
            # If no chains corresponding to the assembly_id remain in the AtomArray
            # expand_assembly will return an empty AtomArray.
            return bioassembly_dict

        # reset the coords after expand assembly
        atom_array.coord[~atom_array.is_resolved, :] = 0.0

        # rename chain_ids from A A B to A0 A1 B0 and add asym_id_int, entity_id_int, sym_id_int
        atom_array = AddAtomArrayAnnot.unique_chain_and_add_ids(atom_array)

        # get chain id before remove chains
        core_indices = self._get_core_indices(atom_array)
        if core_indices is not None:
            ori_chain_ids = np.unique(atom_array.chain_id[core_indices])
        else:
            ori_chain_ids = np.unique(atom_array.chain_id)

        atom_array = AddAtomArrayAnnot.add_mol_id(atom_array)
        atom_array = Filter.remove_unresolved_mols(atom_array)

        # update core indices after remove unresolved mols
        core_indices = np.where(np.isin(atom_array.chain_id, ori_chain_ids))[0]

        # If the number of chains has already reached `max_chains_num`, but the token count hasn't reached `max_tokens_num`,
        # chains will continue to be added until `max_tokens_num` is exceeded.
        atom_array, _input_chains_num = Filter.too_many_chains_filter(
            atom_array,
            core_indices=core_indices,
            max_chains_num=20,
            max_tokens_num=5120,
        )

        if atom_array is None:
            # The distance between the central atoms in any two chains is greater than 15 angstroms.
            return bioassembly_dict

        # update core indices after too_many_chains_filter
        core_indices = np.where(np.isin(atom_array.chain_id, ori_chain_ids))[0]

        atom_array, _removed_chain_ids = Filter.remove_clashing_chains(
            atom_array, core_indices=core_indices
        )

        # remove asymmetric polymer ligand bonds (including protein-protein bond, like disulfide bond)
        # apply to assembly atom array
        atom_array = Filter.remove_asymmetric_polymer_ligand_bonds(
            atom_array, self.entity_poly_type
        )

        # add_mol_id before applying the two filters below to ensure that covalent components are not removed as individual chains.
        atom_array = AddAtomArrayAnnot.find_equiv_mol_and_assign_ids(
            atom_array, self.entity_poly_type
        )

        # numerical encoding of (chain id, residue index)
        atom_array = AddAtomArrayAnnot.add_ref_space_uid(atom_array)
        atom_array = AddAtomArrayAnnot.add_ref_info_and_res_perm(atom_array)

        # the number of protein chains in the assembly
        prot_label_entity_ids = [
            k for k, v in self.entity_poly_type.items() if "polypeptide" in v
        ]
        num_prot_chains = len(
            np.unique(
                atom_array.chain_id[
                    np.isin(atom_array.label_entity_id, prot_label_entity_ids)
                ]
            )
        )
        bioassembly_dict["num_prot_chains"] = num_prot_chains

        bioassembly_dict["atom_array"] = atom_array
        bioassembly_dict["num_tokens"] = atom_array.centre_atom_mask.sum()
        return bioassembly_dict

    @staticmethod
    def create_empty_annotation_like(
        source_array: AtomArray, target_array: AtomArray
    ) -> AtomArray:
        """create empty annotation like source_array"""
        # create empty annotation, atom array addition only keep common annotation
        for k, v in source_array._annot.items():
            if k not in target_array._annot:
                target_array._annot[k] = np.zeros(len(target_array), dtype=v.dtype)
        return target_array

    @staticmethod
    def find_non_ccd_leaving_atoms(
        atom_array: AtomArray,
        select_dict: dict[str, Any],
        component: AtomArray,
    ) -> list[str]:
        """ "
        handle mismatch bettween CCD and mmcif
        some residue has bond in non-central atom (without leaving atoms in CCD)
        and its neighbors should be removed like atom_array from mmcif.

        Args:
            atom_array (AtomArray): Biotite AtomArray object from mmcif.
            select_dict dict[str, Any]: entity_id, res_id, atom_name,... of central atom in atom_array.
            component (AtomArray): CCD component AtomArray object.

        Returns:
            list[str]: list of atom_name to be removed.
        """
        # find non-CCD central atoms in atom_array
        indices_in_atom_array = atom_select(atom_array, select_dict)

        if len(indices_in_atom_array) == 0:
            return []

        if component.bonds is None:
            return []

        # atom_name not in CCD component, return []
        atom_name = select_dict["atom_name"]
        idx_in_comp = np.where(component.atom_name == atom_name)[0]
        if len(idx_in_comp) == 0:
            return []
        idx_in_comp = idx_in_comp[0]

        # find non-CCD leaving atoms in atom_array
        remove_atom_names = []
        for idx in indices_in_atom_array:
            neighbor_idx, types = atom_array.bonds.get_bonds(idx)
            ref_neighbor_idx, types = component.bonds.get_bonds(idx_in_comp)
            # neighbor_atom only bond to central atom in CCD component
            ref_neighbor_idx = [
                i for i in ref_neighbor_idx if len(component.bonds.get_bonds(i)[0]) == 1
            ]
            removed_mask = ~np.isin(
                component.atom_name[ref_neighbor_idx],
                atom_array.atom_name[neighbor_idx],
            )
            remove_atom_names.append(
                component.atom_name[ref_neighbor_idx][removed_mask].tolist()
            )
        max_id = np.argmax(map(len, remove_atom_names))
        return remove_atom_names[max_id]

    def build_ref_chain_with_atom_array(self, atom_array: AtomArray) -> AtomArray:
        """
        build ref chain with atom_array and poly_res_names
        """
        # count inter residue bonds of each potential central atom for removing leaving atoms later
        central_bond_count = Counter()  # (entity_id,res_id,atom_name) -> bond_count

        # build reference entity atom array, including missing residues
        poly_res_names = self.get_poly_res_names(atom_array)
        entity_atom_array = {}
        for entity_id, poly_type in self.entity_poly_type.items():
            chain = struc.AtomArray(0)
            for res_id, res_name in enumerate(poly_res_names[entity_id]):
                # keep all leaving atoms, will remove leaving atoms later in this function
                residue = ccd.get_component_atom_array(
                    res_name, keep_leaving_atoms=True, keep_hydrogens=False
                )
                residue.res_id[:] = res_id + 1
                chain += residue
            res_starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
            inter_bonds = ccd._connect_inter_residue(chain, res_starts)

            # filter out non-std polymer bonds
            bond_mask = np.ones(len(inter_bonds._bonds), dtype=bool)
            for b_idx, (atom_i, atom_j, b_type) in enumerate(inter_bonds._bonds):
                idx_i = atom_select(
                    atom_array,
                    {
                        "label_entity_id": entity_id,
                        "res_id": chain.res_id[atom_i],
                        "atom_name": chain.atom_name[atom_i],
                    },
                )
                idx_j = atom_select(
                    atom_array,
                    {
                        "label_entity_id": entity_id,
                        "res_id": chain.res_id[atom_j],
                        "atom_name": chain.atom_name[atom_j],
                    },
                )
                for i in idx_i:
                    for j in idx_j:
                        # both i, j exist in same chain but not bond in atom_array, non-std polymer bonds, remove from chain
                        if atom_array.chain_id[i] == atom_array.chain_id[j]:
                            bonds, types = atom_array.bonds.get_bonds(i)
                            if j not in bonds:
                                bond_mask[b_idx] = False
                                break

                if bond_mask[b_idx]:
                    # keep this bond, add to central_bond_count
                    central_atom_idx = (
                        atom_i if chain.atom_name[atom_i] in ("C", "P") else atom_j
                    )
                    atom_key = (
                        entity_id,
                        chain.res_id[central_atom_idx],
                        chain.atom_name[central_atom_idx],
                    )
                    # use ref chain bond count if no inter bond in atom_array.
                    central_bond_count[atom_key] = 1

            inter_bonds._bonds = inter_bonds._bonds[bond_mask]
            chain.bonds = chain.bonds.merge(inter_bonds)

            chain.hetero[:] = False
            entity_atom_array[entity_id] = chain

        # remove leaving atoms of residues based on atom_array

        # count inter residue bonds from atom_array for removing leaving atoms later
        inter_residue_bonds = get_inter_residue_bonds(atom_array)
        for i in inter_residue_bonds.flat:
            bonds, types = atom_array.bonds.get_bonds(i)
            bond_count = (
                (atom_array.res_id[bonds] != atom_array.res_id[i])
                | (atom_array.chain_id[bonds] != atom_array.chain_id[i])
            ).sum()
            atom_key = (
                atom_array.label_entity_id[i],
                atom_array.res_id[i],
                atom_array.atom_name[i],
            )
            # remove leaving atoms if central atom has inter residue bond in any copy of a entity
            central_bond_count[atom_key] = max(central_bond_count[atom_key], bond_count)

        # remove leaving atoms for each central atom based in atom_array info
        # so the residue in reference chain can be used directly.
        for entity_id, chain in entity_atom_array.items():
            keep_atom_mask = np.ones(len(chain), dtype=bool)
            starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
            for start, stop in zip(starts[:-1], starts[1:]):
                res_name = chain.res_name[start]
                remove_atom_names = []
                for i in range(start, stop):
                    central_atom_name = chain.atom_name[i]
                    atom_key = (entity_id, chain.res_id[i], central_atom_name)
                    inter_bond_count = central_bond_count[atom_key]

                    if inter_bond_count == 0:
                        continue

                    # num of remove leaving groups equals to num of inter residue bonds (inter_bond_count)
                    component = ccd.get_component_atom_array(
                        res_name, keep_leaving_atoms=True
                    )

                    if component.central_to_leaving_groups is None:
                        # The leaving atoms might be labeled wrongly. The residue remains as it is.
                        break

                    # central_to_leaving_groups:dict[str, list[list[str]]], central atom name to leaving atom groups (atom names).
                    if central_atom_name in component.central_to_leaving_groups:
                        leaving_groups = component.central_to_leaving_groups[
                            central_atom_name
                        ]
                        # removed only when there are leaving atoms.
                        if inter_bond_count >= len(leaving_groups):
                            remove_groups = leaving_groups
                        else:
                            # subsample leaving atoms, keep resolved leaving atoms first
                            exist_group = []
                            not_exist_group = []
                            for group in leaving_groups:
                                for leaving_atom_name in group:
                                    atom_idx = atom_select(
                                        atom_array,
                                        select_dict={
                                            "label_entity_id": entity_id,
                                            "res_id": chain.res_id[i],
                                            "atom_name": leaving_atom_name,
                                        },
                                    )
                                    if len(atom_idx) > 0:  # resolved
                                        exist_group.append(group)
                                        break
                                else:
                                    not_exist_group.append(group)
                            if inter_bond_count <= len(not_exist_group):
                                remove_groups = random.sample(
                                    not_exist_group, inter_bond_count
                                )
                            else:
                                remove_groups = not_exist_group + random.sample(
                                    exist_group, inter_bond_count - len(not_exist_group)
                                )
                        names = [name for group in remove_groups for name in group]
                        remove_atom_names.extend(names)

                    else:
                        # may has non-std leaving atom
                        non_std_leaving_atoms = self.find_non_ccd_leaving_atoms(
                            atom_array=atom_array,
                            select_dict={
                                "label_entity_id": entity_id,
                                "res_id": chain.res_id[i],
                                "atom_name": chain.atom_name[i],
                            },
                            component=component,
                        )
                        if len(non_std_leaving_atoms) > 0:
                            remove_atom_names.extend(non_std_leaving_atoms)

                # remove leaving atoms of this residue
                remove_mask = np.isin(chain.atom_name[start:stop], remove_atom_names)
                keep_atom_mask[np.arange(start, stop)[remove_mask]] = False

            entity_atom_array[entity_id] = chain[keep_atom_mask]
        return entity_atom_array

    @staticmethod
    def make_new_residue(
        atom_array, res_start, res_stop, ref_chain=None
    ) -> tuple[AtomArray, dict[int, int]]:
        """
        make new residue from atom_array[res_start:res_stop], ref_chain is the reference chain.
        1. only remove leavning atom when central atom covalent to other residue.
        2. if ref_chain is provided, remove all atoms not match the residue in ref_chain.
        """
        res_id = atom_array.res_id[res_start]
        res_name = atom_array.res_name[res_start]
        ref_residue = ccd.get_component_atom_array(
            res_name,
            keep_leaving_atoms=True,
            keep_hydrogens=False,
        )
        if ref_residue is None:  # only https://www.rcsb.org/ligand/UNL
            return atom_array[res_start:res_stop]

        if ref_residue.central_to_leaving_groups is None:
            # ambiguous: one leaving group bond to more than one central atom, keep same atoms with PDB entry.
            return atom_array[res_start:res_stop]

        if ref_chain is not None:
            return ref_chain[ref_chain.res_id == res_id]

        keep_atom_mask = np.ones(len(ref_residue), dtype=bool)

        # remove leavning atoms when covalent to other residue
        for i in range(res_start, res_stop):
            central_name = atom_array.atom_name[i]
            old_atom_names = atom_array.atom_name[res_start:res_stop]
            idx = np.where(old_atom_names == central_name)[0]
            if len(idx) == 0:
                # central atom is not resolved in atom_array, not remove leaving atoms
                continue
            idx = idx[0] + res_start
            bonds, types = atom_array.bonds.get_bonds(idx)
            bond_count = (res_id != atom_array.res_id[bonds]).sum()
            if bond_count == 0:
                # central atom is not covalent to other residue, not remove leaving atoms
                continue

            if central_name in ref_residue.central_to_leaving_groups:
                leaving_groups = ref_residue.central_to_leaving_groups[central_name]
                # removed only when there are leaving atoms.
                if bond_count >= len(leaving_groups):
                    remove_groups = leaving_groups
                else:
                    # subsample leaving atoms, remove unresolved leaving atoms first
                    exist_group = []
                    not_exist_group = []
                    for group in leaving_groups:
                        for leaving_atom_name in group:
                            atom_idx = atom_select(
                                atom_array,
                                select_dict={
                                    "chain_id": atom_array.chain_id[i],
                                    "res_id": atom_array.res_id[i],
                                    "atom_name": leaving_atom_name,
                                },
                            )
                            if len(atom_idx) > 0:  # resolved
                                exist_group.append(group)
                                break
                        else:
                            not_exist_group.append(group)

                    # not remove leaving atoms of B and BE, if all leaving atoms is exist in atom_array
                    if central_name in ["B", "BE"]:
                        if not not_exist_group:
                            continue

                    if bond_count <= len(not_exist_group):
                        remove_groups = random.sample(not_exist_group, bond_count)
                    else:
                        remove_groups = not_exist_group + random.sample(
                            exist_group, bond_count - len(not_exist_group)
                        )
            else:
                leaving_atoms = MMCIFParser.find_non_ccd_leaving_atoms(
                    atom_array=atom_array,
                    select_dict={
                        "chain_id": atom_array.chain_id[i],
                        "res_id": atom_array.res_id[i],
                        "atom_name": atom_array.atom_name[i],
                    },
                    component=ref_residue,
                )
                remove_groups = [leaving_atoms]

            names = [name for group in remove_groups for name in group]
            remove_mask = np.isin(ref_residue.atom_name, names)
            keep_atom_mask &= ~remove_mask

        return ref_residue[keep_atom_mask]

    def add_missing_atoms_and_residues(self, atom_array: AtomArray) -> AtomArray:
        """add missing atoms and residues based on CCD and mmcif info.

        Args:
            atom_array (AtomArray): structure with missing residues and atoms, from PDB entry.

        Returns:
            AtomArray: structure added missing residues and atoms (label atom_array.is_resolved as False).
        """
        # build reference entity atom array, including missing residues
        entity_atom_array = self.build_ref_chain_with_atom_array(atom_array)

        # build new atom array and copy info from input atom array to it (new_array).
        new_array = None
        new_global_start = 0
        o2n_amap = {}  # old to new atom map
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        res_starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for c_start, c_stop in zip(chain_starts[:-1], chain_starts[1:]):
            # get reference chain atom array
            entity_id = atom_array.label_entity_id[c_start]
            has_ref_chain = False
            if entity_id in entity_atom_array:
                has_ref_chain = True
                ref_chain_array = entity_atom_array[entity_id].copy()
                ref_chain_array = self.create_empty_annotation_like(
                    atom_array, ref_chain_array
                )

            chain_array = None
            c_res_starts = res_starts[(c_start <= res_starts) & (res_starts <= c_stop)]

            # add missing residues
            prev_res_id = 0
            for r_start, r_stop in zip(c_res_starts[:-1], c_res_starts[1:]):
                curr_res_id = atom_array.res_id[r_start]
                if has_ref_chain and curr_res_id - prev_res_id > 1:
                    # missing residue in head or middle, res_id is 1-based int.
                    segment = ref_chain_array[
                        (prev_res_id < ref_chain_array.res_id)
                        & (ref_chain_array.res_id < curr_res_id)
                    ]
                    if chain_array is None:
                        chain_array = segment
                    else:
                        chain_array += segment

                new_global_start = 0 if new_array is None else len(new_array)
                new_global_start += 0 if chain_array is None else len(chain_array)

                # add missing atoms of existing residue
                ref_chain = ref_chain_array if has_ref_chain else None
                new_residue = self.make_new_residue(
                    atom_array, r_start, r_stop, ref_chain
                )

                new_residue = self.create_empty_annotation_like(atom_array, new_residue)

                # copy residue level info
                residue_fields = ["res_id", "hetero", "label_seq_id", "auth_seq_id"]
                for k in residue_fields:
                    v = atom_array._annot[k][r_start]
                    new_residue._annot[k][:] = v

                # make o2n_amap: old to new atom map
                name_to_index_new = {
                    name: idx for idx, name in enumerate(new_residue.atom_name)
                }
                res_o2n_amap = {}
                res_mismatch_idx = []
                for old_idx in range(r_start, r_stop):
                    old_name = atom_array.atom_name[old_idx]
                    if old_name not in name_to_index_new:
                        # AF3 SI 2.5.4 Filtering
                        # For residues or small molecules with CCD codes, atoms outside of the CCD code’s defined set of atom names are removed.
                        res_mismatch_idx.append(old_idx)
                    else:
                        new_idx = name_to_index_new[old_name]
                        res_o2n_amap[old_idx] = new_global_start + new_idx
                if len(res_o2n_amap) > len(res_mismatch_idx):
                    # Match residues only if more than half of their resolved atoms are matched.
                    # e.g. 1gbt GBS shows 2/12 match, not add to o2n_amap, all atoms are marked as is_resolved=False.
                    o2n_amap.update(res_o2n_amap)

                if chain_array is None:
                    chain_array = new_residue
                else:
                    chain_array += new_residue

                prev_res_id = curr_res_id

            # missing residue in tail
            if has_ref_chain:
                last_res_id = ref_chain_array.res_id[-1]
                if last_res_id > curr_res_id:
                    chain_array += ref_chain_array[ref_chain_array.res_id > curr_res_id]

            # copy chain level info
            chain_fields = [
                "chain_id",
                "label_asym_id",
                "label_entity_id",
                "auth_asym_id",
                # "asym_id_int",
                # "entity_id_int",
                # "sym_id_int",
            ]
            for k in chain_fields:
                chain_array._annot[k][:] = atom_array._annot[k][c_start]

            if new_array is None:
                new_array = chain_array
            else:
                new_array += chain_array

        # copy atom level info
        old_idx = list(o2n_amap.keys())
        new_idx = list(o2n_amap.values())
        atom_fields = ["b_factor", "occupancy", "charge"]
        for k in atom_fields:
            if k not in atom_array._annot:
                continue
            new_array._annot[k][new_idx] = atom_array._annot[k][old_idx]

        # add is_resolved annotation
        is_resolved = np.zeros(len(new_array), dtype=bool)
        is_resolved[new_idx] = True
        new_array.set_annotation("is_resolved", is_resolved)

        # copy coord
        new_array.coord[:] = 0.0
        new_array.coord[new_idx] = atom_array.coord[old_idx]
        # copy bonds
        old_bonds = atom_array.bonds.as_array()  # *n x 3* np.ndarray (i,j,bond_type)

        # some non-leaving atoms are not in the new_array for atom name mismatch, e.g. 4msw TYF
        # only keep bonds of matching atoms
        old_bonds = old_bonds[
            np.isin(old_bonds[:, 0], old_idx) & np.isin(old_bonds[:, 1], old_idx)
        ]

        old_bonds[:, 0] = [o2n_amap[i] for i in old_bonds[:, 0]]
        old_bonds[:, 1] = [o2n_amap[i] for i in old_bonds[:, 1]]
        new_bonds = struc.BondList(len(new_array), old_bonds)
        if new_array.bonds is None:
            new_array.bonds = new_bonds
        else:
            new_array.bonds = new_array.bonds.merge(new_bonds)

        # add peptide bonds and nucleic acid bonds based on CCD type
        new_array = ccd.add_inter_residue_bonds(
            new_array, exclude_struct_conn_pairs=True, remove_far_inter_chain_pairs=True
        )
        return new_array

    def make_chain_indices(
        self, atom_array: AtomArray, pdb_cluster_file: Union[str, Path] = None
    ) -> list:
        """
        Make chain indices.

        Args:
            atom_array (AtomArray): Biotite AtomArray object.
            pdb_cluster_file (Union[str, Path]): cluster info txt file.
        """
        if pdb_cluster_file is None:
            pdb_cluster_dict = {}
        else:
            pdb_cluster_dict = parse_pdb_cluster_file_to_dict(pdb_cluster_file)
        poly_res_names = self.get_poly_res_names(atom_array)
        starts = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        chain_indices_list = []

        is_centre_atom_and_is_resolved = (
            atom_array.is_resolved & atom_array.centre_atom_mask.astype(bool)
        )
        for start, stop in zip(starts[:-1], starts[1:]):
            chain_id = atom_array.chain_id[start]
            entity_id = atom_array.label_entity_id[start]

            # skip if centre atoms within a chain are all unresolved, e.g. 1zc8
            if ~np.any(is_centre_atom_and_is_resolved[start:stop]):
                continue

            # AF3 SI 2.5.1 Weighted PDB dataset
            entity_type = self.entity_poly_type.get(entity_id, "non-poly")

            res_names = poly_res_names.get(entity_id, None)
            if res_names is None:
                chain_atoms = atom_array[start:stop]
                res_ids, res_names = struc.get_residues(chain_atoms)

            if "polypeptide" in entity_type:
                mol_type = "prot"
                sequence = ccd.res_names_to_sequence(res_names)
                if len(sequence) < 10:
                    cluster_id = sequence
                else:
                    pdb_entity = f"{self.pdb_id}_{entity_id}"
                    if pdb_entity in pdb_cluster_dict:
                        cluster_id, _ = pdb_cluster_dict[pdb_entity]
                    elif entity_type == "polypeptide(D)":
                        cluster_id = sequence
                    elif sequence == "X" * len(sequence):
                        chain_atoms = atom_array[start:stop]
                        res_ids, res_names = struc.get_residues(chain_atoms)
                        if np.all(res_names == "UNK"):
                            cluster_id = "poly_UNK"
                        else:
                            cluster_id = "_".join(res_names)
                    else:
                        cluster_id = "NotInClusterTxt"

            elif "ribonucleotide" in entity_type:
                mol_type = "nuc"
                cluster_id = ccd.res_names_to_sequence(res_names)
            else:
                mol_type = "ligand"
                cluster_id = "_".join(res_names)

            chain_dict = {
                "entity_id": entity_id,  # str
                "chain_id": chain_id,
                "mol_type": mol_type,
                "cluster_id": cluster_id,
            }
            chain_indices_list.append(chain_dict)
        return chain_indices_list

    def make_interface_indices(
        self, atom_array: AtomArray, chain_indices_list: list
    ) -> list:
        """make interface indices
        As described in SI 2.5.1, interfaces defined as pairs of chains with minimum heavy atom
        (i.e. non-hydrogen) separation less than 5 Å
        Args:
            atom_array (AtomArray): _description_
            chain_indices_list (List): _description_
        """

        chain_indices_dict = {i["chain_id"]: i for i in chain_indices_list}
        interface_indices_dict = {}

        cell_list = struc.CellList(
            atom_array, cell_size=5, selection=atom_array.is_resolved
        )
        for chain_i, chain_i_dict in chain_indices_dict.items():
            chain_mask = atom_array.chain_id == chain_i
            coord = atom_array.coord[chain_mask & atom_array.is_resolved]
            neighbors_indices_2d = cell_list.get_atoms(
                coord, radius=5
            )  # shape:(n_coord, max_n_neighbors), padding with -1
            neighbors_indices = np.unique(neighbors_indices_2d)
            neighbors_indices = neighbors_indices[neighbors_indices != -1]

            chain_j_list = np.unique(atom_array.chain_id[neighbors_indices])
            for chain_j in chain_j_list:
                if chain_i == chain_j:
                    continue

                # skip if centre atoms within a chain are all unresolved, e.g. 1zc8
                if chain_j not in chain_indices_dict:
                    continue

                interface_id = "_".join(sorted([chain_i, chain_j]))
                if interface_id in interface_indices_dict:
                    continue
                chain_j_dict = chain_indices_dict[chain_j]
                interface_dict = {}
                # chain_id --> chain_1_id
                # mol_type --> mol_1_type
                # entity_id --> entity_1_id
                # cluster_id --> cluster_1_id
                interface_dict.update(
                    {k.replace("_", "_1_"): v for k, v in chain_i_dict.items()}
                )
                interface_dict.update(
                    {k.replace("_", "_2_"): v for k, v in chain_j_dict.items()}
                )
                interface_indices_dict[interface_id] = interface_dict
        return list(interface_indices_dict.values())

    @staticmethod
    def add_sub_mol_type(
        atom_array: AtomArray,
        indices_dict: dict[str, Any],
    ) -> dict[str, Any]:
        """
        Add a "sub_mol_[i]_type" field to indices_dict.
        It includes the following mol_types and sub_mol_types:

        prot
            - prot
            - glycosylation_prot
            - modified_prot

        nuc
            - dna
            - rna
            - modified_dna
            - modified_rna
            - dna_rna_hybrid

        ligand
            - bonded_ligand
            - non_bonded_ligand

        excluded_ligand
            - excluded_ligand

        glycans
            - glycans

        ions
            - ions

        Args:
            atom_array (AtomArray): Biotite AtomArray object of bioassembly.
            indices_dict (dict[str, Any]): A dict of chain or interface indices info.

        Returns:
            dict[str, Any]: A dict of chain or interface indices info with "sub_mol_[i]_type" field.
        """
        polymer_lig_bonds = get_ligand_polymer_bond_mask(atom_array)
        if len(polymer_lig_bonds) == 0:
            lig_polymer_bond_chain_id = []
        else:
            lig_polymer_bond_chain_id = atom_array.chain_id[
                np.unique(polymer_lig_bonds[:, :2])
            ]

        for i in ["1", "2"]:
            if indices_dict[f"entity_{i}_id"] == "":
                indices_dict[f"sub_mol_{i}_type"] = ""
                continue
            entity_type = indices_dict[f"mol_{i}_type"]
            mol_id = atom_array.mol_id[
                atom_array.label_entity_id == indices_dict[f"entity_{i}_id"]
            ][0]
            mol_all_res_name = atom_array.res_name[atom_array.mol_id == mol_id]
            chain_all_mol_type = atom_array.mol_type[
                atom_array.chain_id == indices_dict[f"chain_{i}_id"]
            ]
            chain_all_res_name = atom_array.res_name[
                atom_array.chain_id == indices_dict[f"chain_{i}_id"]
            ]

            if entity_type == "ligand":
                ccd_code = indices_dict[f"cluster_{i}_id"]
                if ccd_code in GLYCANS:
                    indices_dict[f"sub_mol_{i}_type"] = "glycans"

                elif ccd_code in LIGAND_EXCLUSION:
                    indices_dict[f"sub_mol_{i}_type"] = "excluded_ligand"

                elif indices_dict[f"chain_{i}_id"] in lig_polymer_bond_chain_id:
                    indices_dict[f"sub_mol_{i}_type"] = "bonded_ligand"
                else:
                    indices_dict[f"sub_mol_{i}_type"] = "non_bonded_ligand"

            elif entity_type == "prot":
                # glycosylation
                if np.any(np.isin([mol_all_res_name], list(GLYCANS))):
                    indices_dict[f"sub_mol_{i}_type"] = "glycosylation_prot"

                if ~np.all(np.isin(chain_all_res_name, list(PRO_STD_RESIDUES.keys()))):
                    indices_dict[f"sub_mol_{i}_type"] = "modified_prot"

            elif entity_type == "nuc":
                if np.all(chain_all_mol_type == "dna"):
                    if np.any(
                        np.isin(chain_all_res_name, list(DNA_STD_RESIDUES.keys()))
                    ):
                        indices_dict[f"sub_mol_{i}_type"] = "dna"
                    else:
                        indices_dict[f"sub_mol_{i}_type"] = "modified_dna"

                elif np.all(chain_all_mol_type == "rna"):
                    if np.any(
                        np.isin(chain_all_res_name, list(RNA_STD_RESIDUES.keys()))
                    ):
                        indices_dict[f"sub_mol_{i}_type"] = "rna"
                    else:
                        indices_dict[f"sub_mol_{i}_type"] = "modified_rna"
                else:
                    indices_dict[f"sub_mol_{i}_type"] = "dna_rna_hybrid"

            else:
                indices_dict[f"sub_mol_{i}_type"] = [f"mol_{i}_type"]

            if indices_dict.get(f"sub_mol_{i}_type") is None:
                indices_dict[f"sub_mol_{i}_type"] = indices_dict[f"mol_{i}_type"]
        return indices_dict

    @staticmethod
    def add_eval_type(indices_dict: dict[str, Any]) -> dict[str, Any]:
        """
        Differentiate DNA and RNA from the nucleus.

        Args:
            indices_dict (dict[str, Any]): A dict of chain or interface indices info.

        Returns:
            dict[str, Any]: A dict of chain or interface indices info with "eval_type" field.
        """
        if indices_dict["mol_type_group"] not in ["intra_nuc", "nuc_prot"]:
            eval_type = indices_dict["mol_type_group"]
        elif "dna_rna_hybrid" in [
            indices_dict["sub_mol_1_type"],
            indices_dict["sub_mol_2_type"],
        ]:
            eval_type = indices_dict["mol_type_group"]
        else:
            if indices_dict["mol_type_group"] == "intra_nuc":
                nuc_type = str(indices_dict["sub_mol_1_type"]).split("_")[-1]
                eval_type = f"intra_{nuc_type}"
            else:
                nuc_type1 = str(indices_dict["sub_mol_1_type"]).split("_")[-1]
                nuc_type2 = str(indices_dict["sub_mol_2_type"]).split("_")[-1]
                if "dna" in [nuc_type1, nuc_type2]:
                    eval_type = "dna_prot"
                else:
                    eval_type = "rna_prot"
        indices_dict["eval_type"] = eval_type
        return indices_dict

    def make_indices(
        self,
        bioassembly_dict: dict[str, Any],
        pdb_cluster_file: Union[str, Path] = None,
    ) -> list:
        """generate indices of chains and interfaces for sampling data

        Args:
            bioassembly_dict (dict): dict from MMCIFParser.get_bioassembly().
            cluster_file (str): PDB cluster file. Defaults to None.
        Return:
            List(Dict(str, str)): sample_indices_list
        """
        atom_array = bioassembly_dict["atom_array"]
        if atom_array is None:
            print(
                f"Warning: make_indices() input atom_array is None, return empty list (PDB Code:{bioassembly_dict['pdb_id']})"
            )
            return []
        chain_indices_list = self.make_chain_indices(atom_array, pdb_cluster_file)
        interface_indices_list = self.make_interface_indices(
            atom_array, chain_indices_list
        )
        meta_dict = {
            "pdb_id": bioassembly_dict["pdb_id"],
            "assembly_id": bioassembly_dict["assembly_id"],
            "release_date": self.release_date,
            "num_tokens": bioassembly_dict["num_tokens"],
            "num_prot_chains": bioassembly_dict["num_prot_chains"],
            "resolution": self.resolution,
        }
        sample_indices_list = []
        for chain_dict in chain_indices_list:
            chain_dict_out = {k.replace("_", "_1_"): v for k, v in chain_dict.items()}
            chain_dict_out.update(
                {k.replace("_", "_2_"): "" for k, v in chain_dict.items()}
            )
            chain_dict_out["cluster_id"] = chain_dict["cluster_id"]
            chain_dict_out.update(meta_dict)
            chain_dict_out["type"] = "chain"
            sample_indices_list.append(chain_dict_out)

        for interface_dict in interface_indices_list:
            cluster_ids = [
                interface_dict["cluster_1_id"],
                interface_dict["cluster_2_id"],
            ]
            interface_dict["cluster_id"] = ":".join(sorted(cluster_ids))
            interface_dict.update(meta_dict)
            interface_dict["type"] = "interface"
            sample_indices_list.append(interface_dict)

        for indices in sample_indices_list:
            for i in ["1", "2"]:
                chain_id = indices[f"chain_{i}_id"]
                if chain_id == "":
                    continue
                chain_atom_num = np.sum([atom_array.chain_id == chain_id])
                if chain_atom_num == 1:
                    indices[f"mol_{i}_type"] = "ions"

            if indices["type"] == "chain":
                indices["mol_type_group"] = f'intra_{indices["mol_1_type"]}'
            else:
                indices["mol_type_group"] = "_".join(
                    sorted([indices["mol_1_type"], indices["mol_2_type"]])
                )
            indices = self.add_sub_mol_type(atom_array, indices)
            indices = self.add_eval_type(indices)
        return sample_indices_list


class DistillationMMCIFParser(MMCIFParser):

    def get_structure_dict(self) -> dict[str, Any]:
        """
        Get an AtomArray from a CIF file of distillation data.

        Returns:
            Dict[str, Any]: a dict of asymmetric unit structure info.
        """
        # created AtomArray of first model from mmcif atom_site (Asymmetric Unit)
        atom_array = self.get_structure()

        structure_dict = {
            "pdb_id": self.pdb_id,
            "atom_array": None,
            "assembly_id": None,
            "sequences": self.get_sequences(atom_array),
            "entity_poly_type": self.entity_poly_type,
            "num_tokens": -1,
            "num_prot_chains": -1,
        }

        pipeline_functions = [
            self.fix_arginine,
            self.add_missing_atoms_and_residues,  # add UNK
            self.mse_to_met,  # do mse_to_met() after add_missing_atoms_and_residues()
        ]

        for func in pipeline_functions:
            atom_array = func(atom_array)
            if len(atom_array) == 0:
                # no atoms left
                return structure_dict

        atom_array = AddAtomArrayAnnot.add_token_mol_type(
            atom_array, self.entity_poly_type
        )
        atom_array = AddAtomArrayAnnot.add_centre_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_atom_mol_type_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_distogram_rep_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_plddt_m_rep_atom_mask(atom_array)
        atom_array = AddAtomArrayAnnot.add_cano_seq_resname(atom_array)
        atom_array = AddAtomArrayAnnot.add_tokatom_idx(atom_array)
        atom_array = AddAtomArrayAnnot.add_modified_res_mask(atom_array)
        assert (
            atom_array.centre_atom_mask.sum()
            == atom_array.distogram_rep_atom_mask.sum()
        )

        # rename chain_ids from A A B to A0 A1 B0 and add asym_id_int, entity_id_int, sym_id_int
        atom_array = AddAtomArrayAnnot.unique_chain_and_add_ids(atom_array)
        atom_array = AddAtomArrayAnnot.find_equiv_mol_and_assign_ids(
            atom_array, self.entity_poly_type
        )

        # numerical encoding of (chain id, residue index)
        atom_array = AddAtomArrayAnnot.add_ref_space_uid(atom_array)
        atom_array = AddAtomArrayAnnot.add_ref_info_and_res_perm(atom_array)

        # the number of protein chains in the structure
        prot_label_entity_ids = [
            k for k, v in self.entity_poly_type.items() if "polypeptide" in v
        ]
        num_prot_chains = len(
            np.unique(
                atom_array.chain_id[
                    np.isin(atom_array.label_entity_id, prot_label_entity_ids)
                ]
            )
        )
        structure_dict["num_prot_chains"] = num_prot_chains
        structure_dict["atom_array"] = atom_array
        structure_dict["num_tokens"] = atom_array.centre_atom_mask.sum()
        return structure_dict


class AddAtomArrayAnnot(object):
    """
    The methods in this class are all designed to add annotations to an AtomArray
    without altering the information in the original AtomArray.
    """

    @staticmethod
    def add_token_mol_type(
        atom_array: AtomArray, sequences: dict[str, str]
    ) -> AtomArray:
        """
        Add molecule types in atom_arry.mol_type based on ccd pdbx_type.

        Args:
            atom_array (AtomArray): Biotite AtomArray object.
            sequences (dict[str, str]): A dict of label_entity_id --> canonical_sequence

        Return
            AtomArray: add atom_arry.mol_type = "protein" | "rna" | "dna" | "ligand"
        """
        mol_types = np.zeros(len(atom_array), dtype="U7")
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            entity_id = atom_array.label_entity_id[start]
            if entity_id not in sequences:
                # non-poly is ligand
                mol_types[start:stop] = "ligand"
                continue
            res_name = atom_array.res_name[start]

            mol_types[start:stop] = ccd.get_mol_type(res_name)

        atom_array.set_annotation("mol_type", mol_types)
        return atom_array

    @staticmethod
    def add_atom_mol_type_mask(atom_array: AtomArray) -> AtomArray:
        """
        Mask indicates is_protein / rna / dna / ligand.
        It is atom-level which is different with paper (token-level).
        The type of each atom is determined based on the most frequently
        occurring type in the chain to which it belongs.

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with
                       "is_ligand", "is_dna", "is_rna", "is_protein" annotation added.
        """
        # it should be called after mmcif_parser.add_token_mol_type
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        chain_mol_type = []
        for start, end in zip(chain_starts[:-1], chain_starts[1:]):
            mol_types = atom_array.mol_type[start:end]
            mol_type_count = Counter(mol_types)
            most_freq_mol_type = max(mol_type_count, key=mol_type_count.get)
            chain_mol_type.extend([most_freq_mol_type] * (end - start))
        atom_array.set_annotation("chain_mol_type", chain_mol_type)

        for type_str in ["ligand", "dna", "rna", "protein"]:
            mask = (atom_array.chain_mol_type == type_str).astype(int)
            atom_array.set_annotation(f"is_{type_str}", mask)
        return atom_array

    @staticmethod
    def add_modified_res_mask(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI Chapter 5.9.3

        Determine if an atom belongs to a modified residue,
        which is used to calculate the Modified Residue Scores in sample ranking:
        Modified residue scores are ranked according to the average pLDDT of the modified residue.

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with
                       "modified_res_mask" annotation added.
        """
        modified_res_mask = []
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            res_name = atom_array.res_name[start]
            mol_type = atom_array.mol_type[start]
            res_atom_nums = stop - start
            if res_name not in STD_RESIDUES and mol_type != "ligand":
                modified_res_mask.extend([1] * res_atom_nums)
            else:
                modified_res_mask.extend([0] * res_atom_nums)
        atom_array.set_annotation("modified_res_mask", modified_res_mask)
        return atom_array

    @staticmethod
    def add_centre_atom_mask(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI Chapter 2.6
            • A standard amino acid residue (Table 13) is represented as a single token.
            • A standard nucleotide residue (Table 13) is represented as a single token.
            • A modified amino acid or nucleotide residue is tokenized per-atom (i.e. N tokens for an N-atom residue)
            • All ligands are tokenized per-atom
        For each token we also designate a token centre atom, used in various places below:
            • Cα for standard amino acids
            • C1′ for standard nucleotides
            • For other cases take the first and only atom as they are tokenized per-atom.

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with "centre_atom_mask" annotation added.
        """
        res_name = list(STD_RESIDUES.keys())
        std_res = np.isin(atom_array.res_name, res_name) & (
            atom_array.mol_type != "ligand"
        )
        prot_res = np.char.str_len(atom_array.res_name) == 3
        prot_centre_atom = prot_res & (atom_array.atom_name == "CA")
        nuc_centre_atom = (~prot_res) & (atom_array.atom_name == r"C1'")
        not_std_res = ~std_res
        centre_atom_mask = (
            std_res & (prot_centre_atom | nuc_centre_atom)
        ) | not_std_res
        centre_atom_mask = centre_atom_mask.astype(int)
        atom_array.set_annotation("centre_atom_mask", centre_atom_mask)
        return atom_array

    @staticmethod
    def add_distogram_rep_atom_mask(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI Chapter 4.4
        the representative atom mask for each token for distogram head
        • Cβ for protein residues (Cα for glycine),
        • C4 for purines and C2 for pyrimidines.
        • All ligands already have a single atom per token.

        Due to the lack of explanation regarding the handling of "N" and "DN" in the article,
        it is impossible to determine the representative atom based on whether it is a purine or pyrimidine.
        Therefore, C1' is chosen as the representative atom for both "N" and "DN".

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with "distogram_rep_atom_mask" annotation added.
        """
        std_res = np.isin(atom_array.res_name, list(STD_RESIDUES.keys())) & (
            atom_array.mol_type != "ligand"
        )

        # for protein std res
        std_prot_res = std_res & (np.char.str_len(atom_array.res_name) == 3)
        gly = atom_array.res_name == "GLY"
        prot_cb = std_prot_res & (~gly) & (atom_array.atom_name == "CB")
        prot_gly_ca = gly & (atom_array.atom_name == "CA")

        # for nucleotide std res
        purines_c4 = np.isin(atom_array.res_name, ["DA", "DG", "A", "G"]) & (
            atom_array.atom_name == "C4"
        )
        pyrimidines_c2 = np.isin(atom_array.res_name, ["DC", "DT", "C", "U"]) & (
            atom_array.atom_name == "C2"
        )

        # for nucleotide unk res
        unk_nuc = np.isin(atom_array.res_name, ["DN", "N"]) & (
            atom_array.atom_name == r"C1'"
        )

        distogram_rep_atom_mask = (
            prot_cb | prot_gly_ca | purines_c4 | pyrimidines_c2 | unk_nuc
        ) | (~std_res)
        distogram_rep_atom_mask = distogram_rep_atom_mask.astype(int)

        atom_array.set_annotation("distogram_rep_atom_mask", distogram_rep_atom_mask)

        assert np.sum(atom_array.distogram_rep_atom_mask) == np.sum(
            atom_array.centre_atom_mask
        )

        return atom_array

    @staticmethod
    def add_plddt_m_rep_atom_mask(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI Chapter 4.3.1
        the representative atom for plddt loss
        • Atoms such that the distance in the ground truth between atom l and atom m is less than 15 Å
            if m is a protein atom or less than 30 Å if m is a nucleic acid atom.
        • Only atoms in polymer chains.
        • One atom per token - Cα for standard protein residues
            and C1′ for standard nucleic acid residues.

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with "plddt_m_rep_atom_mask" annotation added.
        """
        std_res = np.isin(atom_array.res_name, list(STD_RESIDUES.keys())) & (
            atom_array.mol_type != "ligand"
        )
        ca_or_c1 = (atom_array.atom_name == "CA") | (atom_array.atom_name == r"C1'")
        plddt_m_rep_atom_mask = (std_res & ca_or_c1).astype(int)
        atom_array.set_annotation("plddt_m_rep_atom_mask", plddt_m_rep_atom_mask)
        return atom_array

    @staticmethod
    def add_ref_space_uid(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI Chapter 2.8 Table 5
        Numerical encoding of the chain id and residue index associated with this reference conformer.
        Each (chain id, residue index) tuple is assigned an integer on first appearance.

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with "ref_space_uid" annotation added.
        """
        # [N_atom, 2]
        chain_res_id = np.vstack((atom_array.asym_id_int, atom_array.res_id)).T
        unique_id = np.unique(chain_res_id, axis=0)

        mapping_dict = {}
        for idx, chain_res_id_pair in enumerate(unique_id):
            asym_id_int, res_id = chain_res_id_pair
            mapping_dict[(asym_id_int, res_id)] = idx

        ref_space_uid = [
            mapping_dict[(asym_id_int, res_id)] for asym_id_int, res_id in chain_res_id
        ]
        atom_array.set_annotation("ref_space_uid", ref_space_uid)
        return atom_array

    @staticmethod
    def add_cano_seq_resname(atom_array: AtomArray) -> AtomArray:
        """
        Assign to each atom the three-letter residue name (resname)
        corresponding to its place in the canonical sequences.
        Non-standard residues are mapped to standard ones.
        Residues that cannot be mapped to standard residues and ligands are all labeled as "UNK".

        Note: Some CCD Codes in the canonical sequence are mapped to three letters. It is labeled as one "UNK".

        Args:
            atom_array (AtomArray): Biotite AtomArray object

        Returns:
            AtomArray: Biotite AtomArray object with "cano_seq_resname" annotation added.
        """
        cano_seq_resname = []
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            res_atom_nums = stop - start
            mol_type = atom_array.mol_type[start]
            resname = atom_array.res_name[start]

            one_letter_code = ccd.get_one_letter_code(resname)
            if one_letter_code is None or len(one_letter_code) != 1:
                # Some non-standard residues cannot be mapped back to one standard residue.
                one_letter_code = "X" if mol_type == "protein" else "N"

            if mol_type == "protein":
                res_name_in_cano_seq = PROT_STD_RESIDUES_ONE_TO_THREE.get(
                    one_letter_code, "UNK"
                )
            elif mol_type == "dna":
                res_name_in_cano_seq = "D" + one_letter_code
                if res_name_in_cano_seq not in DNA_STD_RESIDUES:
                    res_name_in_cano_seq = "DN"
            elif mol_type == "rna":
                res_name_in_cano_seq = one_letter_code
                if res_name_in_cano_seq not in RNA_STD_RESIDUES:
                    res_name_in_cano_seq = "N"
            else:
                # some molecules attached to a polymer like ATP-RNA. e.g.
                res_name_in_cano_seq = "UNK"

            cano_seq_resname.extend([res_name_in_cano_seq] * res_atom_nums)

        atom_array.set_annotation("cano_seq_resname", cano_seq_resname)
        return atom_array

    @staticmethod
    def remove_bonds_between_polymer_chains(
        atom_array: AtomArray, entity_poly_type: dict[str, str]
    ) -> struc.BondList:
        """
        Remove bonds between polymer chains based on entity_poly_type

        Args:
            atom_array (AtomArray): Biotite AtomArray object
            entity_poly_type (dict[str, str]): entity_id to poly_type

        Returns:
            BondList: Biotite BondList object (copy) with bonds between polymer chains removed
        """
        copy = atom_array.bonds.copy()
        polymer_mask = np.isin(
            atom_array.label_entity_id, list(entity_poly_type.keys())
        )
        i = copy._bonds[:, 0]
        j = copy._bonds[:, 1]
        pp_bond_mask = polymer_mask[i] & polymer_mask[j]
        diff_chain_mask = atom_array.chain_id[i] != atom_array.chain_id[j]
        pp_bond_mask = pp_bond_mask & diff_chain_mask
        copy._bonds = copy._bonds[~pp_bond_mask]

        # post-process after modified bonds manually
        # due to the extraction of bonds using a mask, the lower one of the two atom indices is still in the first
        copy._remove_redundant_bonds()
        copy._max_bonds_per_atom = copy._get_max_bonds_per_atom()
        return copy

    @staticmethod
    def find_equiv_mol_and_assign_ids(
        atom_array: AtomArray,
        entity_poly_type: Optional[dict[str, str]] = None,
        check_final_equiv: bool = True,
    ) -> AtomArray:
        """
        Assign a unique integer to each molecule in the structure.
        All atoms connected by covalent bonds are considered as a molecule, with unique mol_id (int).
        different copies of same molecule will assign same entity_mol_id (int).
        for each mol, assign mol_atom_index starting from 0.

        Args:
            atom_array (AtomArray): Biotite AtomArray object
            entity_poly_type (Optional[dict[str, str]]): label_entity_id to entity.poly_type.
                              Defaults to None.
            check_final_equiv (bool, optional): check if the final mol_ids of same entity_mol_id are all equivalent.

        Returns:
            AtomArray: Biotite AtomArray object with new annotations
            - mol_id: atoms with covalent bonds connected, 0-based int
            - entity_mol_id: equivalent molecules will assign same entity_mol_id, 0-based int
            - mol_residue_index: mol_atom_index for each mol, 0-based int
        """
        # Re-assign mol_id to AtomArray after break asym bonds
        if entity_poly_type is None:
            mol_indices: list[np.ndarray] = get_molecule_indices(atom_array)
        else:
            bonds_filtered = AddAtomArrayAnnot.remove_bonds_between_polymer_chains(
                atom_array, entity_poly_type
            )
            mol_indices: list[np.ndarray] = get_molecule_indices(bonds_filtered)

        # assign mol_id
        mol_ids = np.array([-1] * len(atom_array), dtype=np.int32)
        for mol_id, atom_indices in enumerate(mol_indices):
            mol_ids[atom_indices] = mol_id
        atom_array.set_annotation("mol_id", mol_ids)

        assert ~np.isin(-1, atom_array.mol_id), "Some mol_id is not assigned."
        assert len(np.unique(atom_array.mol_id)) == len(
            mol_indices
        ), "Some mol_id is duplicated."

        # assign entity_mol_id
        # --------------------
        # first atom of mol with infos in attrubites, eg: info.num_atoms, info.bonds, ...
        ref_mol_infos = []
        # perm for keep multiple chains in one mol are together and in same chain order
        new_atom_perm = []
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=False)
        entity_mol_ids = np.zeros_like(mol_ids)
        for mol_id, atom_indices in enumerate(mol_indices):
            atom_indices = np.sort(atom_indices)
            # keep multiple chains-mol has same chain order in different copies
            chain_perm = np.argsort(
                atom_array.label_entity_id[atom_indices], kind="stable"
            )
            atom_indices = atom_indices[chain_perm]
            # save indices for finally re-ordering atom_array
            new_atom_perm.extend(atom_indices)

            # check mol equal, keep chain order consistent with atom_indices
            mol_chain_mask = np.isin(atom_indices, chain_starts)
            entity_ids = atom_array.label_entity_id[atom_indices][
                mol_chain_mask
            ].tolist()

            match_entity_mol_id = None
            for entity_mol_id, mol_info in enumerate(ref_mol_infos):
                # check mol equal
                # same entity_ids and same atom name will assign same entity_mol_id
                if entity_ids != mol_info.entity_ids:
                    continue

                if len(atom_indices) != len(mol_info.atom_name):
                    continue

                atom_name_not_equal = (
                    atom_array.atom_name[atom_indices] != mol_info.atom_name
                )
                if np.any(atom_name_not_equal):
                    diff_indices = np.where(atom_name_not_equal)[0]
                    query_atom = atom_array[atom_indices[diff_indices[0]]]
                    ref_atom = atom_array[mol_info.atom_indices[diff_indices[0]]]
                    logger.warning(
                        f"Two mols have entity_ids and same number of atoms, but diff atom name:\n{query_atom=}\n{  ref_atom=}"
                    )
                    continue

                # pass all checks, it is a match
                match_entity_mol_id = entity_mol_id
                break

            if match_entity_mol_id is None:  # not found match mol
                # use first atom as a placeholder for mol info.
                mol_info = atom_array[atom_indices[0]]
                mol_info.atom_indices = atom_indices
                mol_info.entity_ids = entity_ids
                mol_info.atom_name = atom_array.atom_name[atom_indices]
                mol_info.entity_mol_id = len(ref_mol_infos)
                ref_mol_infos.append(mol_info)
                match_entity_mol_id = mol_info.entity_mol_id

            entity_mol_ids[atom_indices] = match_entity_mol_id

        atom_array.set_annotation("entity_mol_id", entity_mol_ids)

        # re-order atom_array to make atoms with same mol_id together.
        atom_array = atom_array[new_atom_perm]

        # assign mol_atom_index
        mol_starts = get_starts_by(
            atom_array, by_annot="mol_id", add_exclusive_stop=True
        )
        mol_atom_index = np.zeros_like(atom_array.mol_id, dtype=np.int32)
        for start, stop in zip(mol_starts[:-1], mol_starts[1:]):
            mol_atom_index[start:stop] = np.arange(stop - start)
        atom_array.set_annotation("mol_atom_index", mol_atom_index)

        # check mol equivalence again
        if check_final_equiv:
            num_mols = len(mol_starts) - 1
            for i in range(num_mols):
                for j in range(i + 1, num_mols):
                    start_i, stop_i = mol_starts[i], mol_starts[i + 1]
                    start_j, stop_j = mol_starts[j], mol_starts[j + 1]
                    if (
                        atom_array.entity_mol_id[start_i]
                        != atom_array.entity_mol_id[start_j]
                    ):
                        continue
                    for key in ["res_name", "atom_name", "mol_atom_index"]:
                        # not check res_id for ligand may have different res_id
                        annot = getattr(atom_array, key)
                        assert np.all(
                            annot[start_i:stop_i] == annot[start_j:stop_j]
                        ), f"not equal {key} when find_equiv_mol_and_assign_ids()"

        return atom_array

    @staticmethod
    def add_tokatom_idx(atom_array: AtomArray) -> AtomArray:
        """
        Add a tokatom_idx corresponding to the residue and atom name for each atom.
        For non-standard residues or ligands, the tokatom_idx should be set to 0.

        Parameters:
        atom_array (AtomArray): The AtomArray object to which the annotation will be added.

        Returns:
        AtomArray: The AtomArray object with the 'tokatom_idx' annotation added.
        """
        # pre-defined atom name order for tokatom_idx
        tokatom_idx_list = []
        for atom in atom_array:
            atom_name_position = RES_ATOMS_DICT.get(atom.res_name, None)
            if atom.mol_type == "ligand" or atom_name_position is None:
                tokatom_idx = 0
            else:
                tokatom_idx = atom_name_position[atom.atom_name]
            tokatom_idx_list.append(tokatom_idx)
        atom_array.set_annotation("tokatom_idx", tokatom_idx_list)
        return atom_array

    @staticmethod
    def add_mol_id(atom_array: AtomArray) -> AtomArray:
        """
        Assign a unique integer to each molecule in the structure.

        Args:
            atom_array (AtomArray): Biotite AtomArray object
        Returns:
            AtomArray: Biotite AtomArray object with new annotations
            - mol_id: atoms with covalent bonds connected, 0-based int
        """
        mol_indices = get_molecule_indices(atom_array)

        # assign mol_id
        mol_ids = np.array([-1] * len(atom_array), dtype=np.int32)
        for mol_id, atom_indices in enumerate(mol_indices):
            mol_ids[atom_indices] = mol_id
        atom_array.set_annotation("mol_id", mol_ids)
        return atom_array

    @staticmethod
    def unique_chain_and_add_ids(atom_array: AtomArray) -> AtomArray:
        """
        Unique chain ID and add asym_id, entity_id, sym_id.
        Adds a number to the chain ID to make chain IDs in the assembly unique.
        Example: [A, B, A, B, C] ==> [A0, B0, A1, B1, C0]

        Args:
            atom_array (AtomArray): Biotite AtomArray object.

        Returns:
            AtomArray: Biotite AtomArray object with new annotations:
                - asym_id_int: np.array(int)
                - entity_id_int: np.array(int)
                - sym_id_int: np.array(int)
        """
        entity_id_uniq = np.sort(np.unique(atom_array.label_entity_id))
        entity_id_dict = {e: i for i, e in enumerate(entity_id_uniq)}
        asym_ids = np.zeros(len(atom_array), dtype=int)
        entity_ids = np.zeros(len(atom_array), dtype=int)
        sym_ids = np.zeros(len(atom_array), dtype=int)
        chain_ids = np.zeros(len(atom_array), dtype="U4")
        counter = Counter()
        start_indices = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        for i in range(len(start_indices) - 1):
            start_i = start_indices[i]
            stop_i = start_indices[i + 1]
            asym_ids[start_i:stop_i] = i

            entity_id = atom_array.label_entity_id[start_i]
            entity_ids[start_i:stop_i] = entity_id_dict[entity_id]

            sym_ids[start_i:stop_i] = counter[entity_id]
            counter[entity_id] += 1
            new_chain_id = f"{atom_array.chain_id[start_i]}{sym_ids[start_i]}"
            chain_ids[start_i:stop_i] = new_chain_id

        atom_array.set_annotation("asym_id_int", asym_ids)
        atom_array.set_annotation("entity_id_int", entity_ids)
        atom_array.set_annotation("sym_id_int", sym_ids)
        atom_array.chain_id = chain_ids
        return atom_array

    @staticmethod
    def add_int_id(atom_array):
        """
        Unique chain ID and add asym_id, entity_id, sym_id.
        Adds a number to the chain ID to make chain IDs in the assembly unique.
        Example: [A, B, A, B, C] ==> [A0, B0, A1, B1, C0]

        Args:
            atom_array (AtomArray): Biotite AtomArray object.

        Returns:
            AtomArray: Biotite AtomArray object with new annotations:
                - asym_id_int: np.array(int)
                - entity_id_int: np.array(int)
                - sym_id_int: np.array(int)
        """
        entity_id_uniq = np.sort(np.unique(atom_array.label_entity_id))
        entity_id_dict = {e: i for i, e in enumerate(entity_id_uniq)}
        asym_ids = np.zeros(len(atom_array), dtype=int)
        entity_ids = np.zeros(len(atom_array), dtype=int)
        sym_ids = np.zeros(len(atom_array), dtype=int)
        counter = Counter()
        start_indices = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        for i in range(len(start_indices) - 1):
            start_i = start_indices[i]
            stop_i = start_indices[i + 1]
            asym_ids[start_i:stop_i] = i

            entity_id = atom_array.label_entity_id[start_i]
            entity_ids[start_i:stop_i] = entity_id_dict[entity_id]

            sym_ids[start_i:stop_i] = counter[entity_id]
            counter[entity_id] += 1

        atom_array.set_annotation("asym_id_int", asym_ids)
        atom_array.set_annotation("entity_id_int", entity_ids)
        atom_array.set_annotation("sym_id_int", sym_ids)
        return atom_array

    @staticmethod
    def add_ref_feat_info(
        atom_array: AtomArray,
    ) -> tuple[np.ndarray, np.ndarray, list[int]]:
        """
        Get info of reference structure of atoms based on the atom array.

        Args:
            atom_array (AtomArray): The atom array.

        Returns:
            tuple:
                ref_pos (numpy.ndarray): Atom positions in the reference conformer,
                                         with a random rotation and translation applied.
                                         Atom positions are given in Å. Shape=(num_atom, 3).
                ref_charge (numpy.ndarray): Charge for each atom in the reference conformer. Shape=(num_atom）
                ref_mask ((numpy.ndarray): Mask indicating which atom slots are used in the reference conformer. Shape=(num_atom）
        """
        info_dict = {}
        for ccd_id in np.unique(atom_array.res_name):
            # create ref conformer for each CCD ID
            ref_result = get_ccd_ref_info(ccd_id)
            if ref_result:
                for space_uid in np.unique(
                    atom_array[atom_array.res_name == ccd_id].ref_space_uid
                ):
                    if ref_result:
                        info_dict[space_uid] = [
                            ref_result["atom_map"],
                            ref_result["coord"],
                            ref_result["charge"],
                            ref_result["mask"],
                        ]
            else:
                # get conformer failed will result in an empty dictionary
                continue

        ref_mask = []  # [N_atom]
        ref_pos = []  # [N_atom, 3]
        ref_charge = []  # [N_atom]
        for atom in atom_array:
            ref_result = info_dict.get(atom.ref_space_uid)
            if ref_result is None:
                # get conformer failed
                ref_mask.append(0)
                ref_pos.append([0.0, 0.0, 0.0])
                ref_charge.append(0)

            else:
                atom_map, coord, charge, mask = ref_result
                atom_sub_idx = atom_map[atom.atom_name]
                ref_mask.append(mask[atom_sub_idx])
                ref_pos.append(coord[atom_sub_idx])
                ref_charge.append(charge[atom_sub_idx])

        ref_pos = np.array(ref_pos)
        ref_charge = np.array(ref_charge).astype(int)
        ref_mask = np.array(ref_mask).astype(int)
        return ref_pos, ref_charge, ref_mask

    @staticmethod
    def add_res_perm(
        atom_array: AtomArray,
    ) -> tuple[np.ndarray, np.ndarray, list[int]]:
        """
        Get permutations of each atom within the residue.

        Args:
            atom_array (AtomArray): biotite AtomArray object.

        Returns:
            list[list[int]]: 2D list of (N_atom, N_perm)
        """
        starts = get_residue_starts(atom_array, add_exclusive_stop=True)
        res_perm = []
        for start, stop in zip(starts[:-1], starts[1:]):
            res_atom = atom_array[start:stop]
            curr_res_atom_idx = list(range(len(res_atom)))

            res_dict = get_ccd_ref_info(ccd_code=res_atom.res_name[0])
            if not res_dict:
                res_perm.extend([[i] for i in curr_res_atom_idx])
                continue

            perm_array = res_dict["perm"]  # [N_atoms, N_perm]
            perm_atom_idx_in_res_order = [
                res_dict["atom_map"][i] for i in res_atom.atom_name
            ]
            perm_idx_to_present_atom_idx = dict(
                zip(perm_atom_idx_in_res_order, curr_res_atom_idx)
            )

            precent_row_mask = np.isin(perm_array[:, 0], perm_atom_idx_in_res_order)
            perm_array_row_filtered = perm_array[precent_row_mask]

            precent_col_mask = np.isin(
                perm_array_row_filtered, perm_atom_idx_in_res_order
            ).all(axis=0)
            perm_array_filtered = perm_array_row_filtered[:, precent_col_mask]

            # replace the elem in new_perm_array according to the perm_idx_to_present_atom_idx dict
            new_perm_array = np.vectorize(perm_idx_to_present_atom_idx.get)(
                perm_array_filtered
            )

            assert (
                new_perm_array.shape[1] <= 1000
                and new_perm_array.shape[1] <= perm_array.shape[1]
            )
            res_perm.extend(new_perm_array.tolist())
        return res_perm

    @staticmethod
    def add_ref_info_and_res_perm(atom_array: AtomArray) -> AtomArray:
        """
        Add info of reference structure of atoms to the atom array.

        Args:
            atom_array (AtomArray): The atom array.

        Returns:
            AtomArray: The atom array with the 'ref_pos', 'ref_charge', 'ref_mask', 'res_perm' annotations added.
        """
        ref_pos, ref_charge, ref_mask = AddAtomArrayAnnot.add_ref_feat_info(atom_array)
        res_perm = AddAtomArrayAnnot.add_res_perm(atom_array)

        str_res_perm = []  # encode [N_atom, N_perm] -> list[str]
        for i in res_perm:
            str_res_perm.append("_".join([str(j) for j in i]))

        assert (
            len(atom_array)
            == len(ref_pos)
            == len(ref_charge)
            == len(ref_mask)
            == len(res_perm)
        ), f"{len(atom_array)=}, {len(ref_pos)=}, {len(ref_charge)=}, {len(ref_mask)=}, {len(str_res_perm)=}"

        atom_array.set_annotation("ref_pos", ref_pos)
        atom_array.set_annotation("ref_charge", ref_charge)
        atom_array.set_annotation("ref_mask", ref_mask)
        atom_array.set_annotation("res_perm", str_res_perm)
        return atom_array
