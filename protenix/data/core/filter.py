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

from collections import defaultdict

import biotite.structure as struc
import networkx as nx
import numpy as np
from biotite.structure import AtomArray, get_residue_starts
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist

from protenix.data.constants import CRYSTALLIZATION_AIDS
from protenix.data.utils import get_inter_residue_bonds


class Filter(object):
    """
    Ref: AlphaFold3 SI Chapter 2.5.4
    """

    @staticmethod
    def remove_hydrogens(atom_array: AtomArray) -> AtomArray:
        """remove hydrogens and deuteriums"""
        return atom_array[~np.isin(atom_array.element, ["H", "D"])]

    @staticmethod
    def remove_water(atom_array: AtomArray) -> AtomArray:
        """remove water (HOH) and deuterated water (DOD)"""
        return atom_array[~np.isin(atom_array.res_name, ["HOH", "DOD"])]

    @staticmethod
    def remove_element_X(atom_array: AtomArray) -> AtomArray:
        """
        remove element X
        following residues have element X:
        - UNX: unknown one atom or ion
        - UNL: unknown ligand, some atoms are marked as X
        - ASX: ASP/ASN ambiguous, two ambiguous atoms are marked as X, 6 entries in the PDB
        - GLX: GLU/GLN ambiguous, two ambiguous atoms are marked as X, 5 entries in the PDB
        """
        X_mask = np.zeros(len(atom_array), dtype=bool)
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            res_name = atom_array.res_name[start]
            if res_name in ["UNX", "UNL"]:
                X_mask[start:stop] = True
        atom_array = atom_array[~X_mask]

        # map ASX to ASP, as ASP is more symmetric than ASN
        mask = atom_array.res_name == "ASX"
        atom_array.res_name[mask] = "ASP"
        atom_array.atom_name[mask & (atom_array.atom_name == "XD1")] = "OD1"
        atom_array.atom_name[mask & (atom_array.atom_name == "XD2")] = "OD2"
        atom_array.element[mask & (atom_array.element == "X")] = "O"

        # map GLX to GLU, as GLU is more symmetric than GLN
        mask = atom_array.res_name == "GLX"
        atom_array.res_name[mask] = "GLU"
        atom_array.atom_name[mask & (atom_array.atom_name == "XE1")] = "OE1"
        atom_array.atom_name[mask & (atom_array.atom_name == "XE2")] = "OE2"
        atom_array.element[mask & (atom_array.element == "X")] = "O"
        return atom_array

    @staticmethod
    def remove_crystallization_aids(
        atom_array: AtomArray, entity_poly_type: dict
    ) -> AtomArray:
        """remove crystallization aids, eg: SO4, GOL, etc.

        Only remove crystallization aids if the chain is not polymer.

        Ref: AlphaFold3 SI Chapter 2.5.4
        """
        non_aids_mask = ~np.isin(atom_array.res_name, CRYSTALLIZATION_AIDS)
        poly_mask = np.isin(atom_array.label_entity_id, list(entity_poly_type.keys()))
        return atom_array[poly_mask | non_aids_mask]

    @staticmethod
    def remove_polymer_chains_all_residues_unknown(
        atom_array: AtomArray,
        entity_poly_type: dict,
    ) -> AtomArray:
        """remove chains with all residues unknown"""
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        invalid_chains = []  # list of [start, end)
        for index in range(len(chain_starts) - 1):
            start, end = chain_starts[index], chain_starts[index + 1]
            entity_id = atom_array.label_entity_id[start]
            if (
                entity_poly_type.get(entity_id, "non-poly") == "polypeptide(L)"
                and np.all(atom_array.res_name[start:end] == "UNK")
            ) or (
                entity_poly_type.get(entity_id, "non-poly")
                in (
                    "polyribonucleotide",
                    "polydeoxyribonucleotide",
                )
                and np.all(atom_array.res_name[start:end] == "N")
            ):
                invalid_chains.append((start, end))
        mask = np.ones(len(atom_array), dtype=bool)
        for start, end in invalid_chains:
            mask[start:end] = False
        atom_array = atom_array[mask]
        return atom_array

    @staticmethod
    def remove_polymer_chains_too_short(
        atom_array: AtomArray, entity_poly_type: dict
    ) -> AtomArray:
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        invalid_chains = []  # list of [start, end)
        for index in range(len(chain_starts) - 1):
            start, end = chain_starts[index], chain_starts[index + 1]
            entity_id = atom_array.label_entity_id[start]
            num_residue_ids = len(set(atom_array.label_seq_id[start:end]))
            if (
                entity_poly_type.get(entity_id, "non-poly")
                in (
                    "polypeptide(L)",  # TODO: how to handle polypeptide(D)?
                    "polyribonucleotide",
                    "polydeoxyribonucleotide",
                )
                and num_residue_ids < 4
            ):
                invalid_chains.append((start, end))
        mask = np.ones(len(atom_array), dtype=bool)
        for start, end in invalid_chains:
            mask[start:end] = False
        atom_array = atom_array[mask]
        return atom_array

    @staticmethod
    def remove_polymer_chains_with_consecutive_c_alpha_too_far_away(
        atom_array: AtomArray, entity_poly_type: dict, max_distance: float = 10.0
    ) -> AtomArray:
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        invalid_chains = []  # list of [start, end)
        for index in range(len(chain_starts) - 1):
            start, end = chain_starts[index], chain_starts[index + 1]
            entity_id = atom_array.label_entity_id[start]
            if entity_poly_type.get(entity_id, "non-poly") == "polypeptide(L)":
                peptide_atoms = atom_array[start:end]
                ca_atoms = peptide_atoms[peptide_atoms.atom_name == "CA"]
                seq_ids = ca_atoms.label_seq_id
                seq_ids[seq_ids == "."] = "-100"
                seq_ids = seq_ids.astype(np.int64)
                dist_square = np.sum(
                    (ca_atoms[:-1].coord - ca_atoms[1:].coord) ** 2, axis=-1
                )
                invalid_neighbor_mask = (dist_square > max_distance**2) & (
                    seq_ids[:-1] + 1 == seq_ids[1:]
                )
                if np.any(invalid_neighbor_mask):
                    invalid_chains.append((start, end))
        mask = np.ones(len(atom_array), dtype=bool)
        for start, end in invalid_chains:
            mask[start:end] = False
        atom_array = atom_array[mask]
        return atom_array

    @staticmethod
    def too_many_chains_filter(
        atom_array: AtomArray,
        interface_radius: int = 15,
        max_chains_num: int = 20,
        core_indices: list[int] = None,
        max_tokens_num: int = None,
    ) -> tuple[AtomArray, int]:
        """
        Ref: AlphaFold3 SI Chapter 2.5.4

        For bioassemblies with greater than 20 chains, we select a random interface token
        (with a centre atom <15 Å to the centre atom of a token in another chain)
        and select the closest 20 chains to this token based on
        minimum distance between any tokens centre atom.

        Note: due to the presence of covalent small molecules,
        treat the covalent small molecule and the polymer it is attached to
        as a single chain to avoid inadvertently removing the covalent small molecules.
        Use the mol_id added to the AtomArray to differentiate between the various
        parts of the structure composed of covalent bonds.

        Args:
            atom_array (AtomArray): Biotite AtomArray Object of a Bioassembly.
            interface_radius (int, optional): Atoms within this distance of the central atom are considered interface atoms.
                                            Defaults to 15.
            max_chains_num (int, optional): The maximum number of chains permitted in a bioassembly.
                                            Filtration will be applied if exceeds this value. Defaults to 20.
            core_indices (list[int], optional): A list of indices to be used as chose the central atom.
                                                     And corresponding chains in the list will be selected proriority.
                                                     If None, a random index from whole AtomArray will be selected. Defaults to None.
            max_tokens_num (int, optional): The maximum number of tokens permitted in a bioassembly.
                                            If not None,  after more than max_chains_num, if the max_tokens_num is not reached,
                                            it will continue to append the chains.

        Returns:
            tuple:
                - atom_array (AtomArray): An AtomArray that has been processed through this filter.
                - input_chains_num (int): The number of chain in the input AtomArray.
                                          This is to log whether the filter has been utilized.
        """
        # each mol is a so called "chain" in the context of this filter.
        input_chains_num = len(np.unique(atom_array.mol_id))
        if input_chains_num <= max_chains_num:
            # no change
            return atom_array, input_chains_num

        is_resolved_centre_atom = (
            atom_array.centre_atom_mask == 1
        ) & atom_array.is_resolved

        cell_list = struc.CellList(
            atom_array, cell_size=interface_radius, selection=is_resolved_centre_atom
        )
        resolved_centre_atom = atom_array[is_resolved_centre_atom]

        assert resolved_centre_atom, "There is no resolved central atom."

        # random pick centre atom
        if core_indices is None:
            index_shuf = np.random.default_rng(seed=42).permutation(
                len(resolved_centre_atom)
            )
        else:
            index_shuf = np.array(core_indices)
            resolved_centre_atom_indices = np.nonzero(is_resolved_centre_atom)[0]

            # get indices of resolved_centre_atom
            index_shuf = np.array(
                [
                    np.where(resolved_centre_atom_indices == idx)[0][0]
                    for idx in index_shuf
                    if idx in resolved_centre_atom_indices
                ]
            )
            np.random.default_rng(seed=42).shuffle(index_shuf)

        chosen_centre_atom = None
        for idx in index_shuf:
            centre_atom = resolved_centre_atom[idx]
            neighbors_indices = cell_list.get_atoms(
                centre_atom.coord, radius=interface_radius
            )
            neighbors_indices = neighbors_indices[neighbors_indices != -1]

            neighbors_chain_ids = np.unique(atom_array.mol_id[neighbors_indices])
            # neighbors include centre atom itself
            if len(neighbors_chain_ids) > 1:
                chosen_centre_atom = centre_atom
                break

        # The distance between the central atoms in any two chains is greater than 15 angstroms.
        if chosen_centre_atom is None:
            return None, input_chains_num

        dist_mat = cdist(centre_atom.coord.reshape((1, -1)), resolved_centre_atom.coord)
        sorted_chain_id = np.array(
            [
                chain_id
                for chain_id, _dist in sorted(
                    zip(resolved_centre_atom.mol_id, dist_mat[0]),
                    key=lambda pair: pair[1],
                )
            ]
        )

        if core_indices is not None:
            # select core proriority
            core_mol_id = np.unique(atom_array.mol_id[core_indices])
            in_core_mask = np.isin(sorted_chain_id, core_mol_id)
            sorted_chain_id = np.concatenate(
                (sorted_chain_id[in_core_mask], sorted_chain_id[~in_core_mask])
            )

        closest_chain_id = set()
        chain_ids_to_token_num = {}
        if max_tokens_num is None:
            max_tokens_num = 0

        tokens = 0
        for chain_id in sorted_chain_id:
            # get token num
            if chain_id not in chain_ids_to_token_num:
                chain_ids_to_token_num[chain_id] = atom_array.centre_atom_mask[
                    atom_array.mol_id == chain_id
                ].sum()
            chain_token_num = chain_ids_to_token_num[chain_id]

            if len(closest_chain_id) >= max_chains_num:
                if tokens + chain_token_num > max_tokens_num:
                    break

            closest_chain_id.add(chain_id)
            tokens += chain_token_num

        atom_array = atom_array[np.isin(atom_array.mol_id, list(closest_chain_id))]
        output_chains_num = len(np.unique(atom_array.mol_id))
        assert (
            output_chains_num == max_chains_num
            or atom_array.centre_atom_mask.sum() <= max_tokens_num
        )
        return atom_array, input_chains_num

    @staticmethod
    def _get_clash_mask(
        atom_array: AtomArray,
        kdtree_query_result: np.ndarray,
        removed_chain_ids: set[str],
    ) -> np.ndarray:
        """
        Identify atomic clashes between chains using KDTree query results.

        Implements clash detection criteria from AlphaFold3 (SI Chapter 2.5.4):
        - Atoms within 1.7 Å considered clashing (but we use 1.6 Å here)
        - Excludes intra-chain clashes and covalent bonds
        - Respects chain removal priority from previous iterations

        Args:
            atom_array (AtomArray): Full atomic structure data
            kdtree_query_result (np.ndarray): Precomputed neighbor indices from KDTree.query_ball_point
            removed_chain_ids (set[str]): Chains already scheduled for removal

        Returns:
            np.ndarray: Boolean array where True indicates atoms involved in inter-chain clashes
            meeting removal criteria
        """
        clash_mask = np.zeros(len(atom_array), dtype=bool)
        for atom_idx, clashes in enumerate(kdtree_query_result):
            chain_id_i = atom_array.chain_id[atom_idx]
            if chain_id_i in removed_chain_ids:
                continue

            if len(clashes) == 1:
                # only clash with itself
                continue

            bonded_atoms, _bond_types = atom_array.bonds.get_bonds(atom_idx)
            for clashed_atom_idx in clashes:
                chain_id_j = atom_array.chain_id[clashed_atom_idx]
                if chain_id_j in removed_chain_ids:
                    continue

                # clash with itself
                if chain_id_i == chain_id_j:
                    continue

                # if two atoms are covalent, they are not considered as clashing.
                if clashed_atom_idx in bonded_atoms:
                    continue
                clash_mask[atom_idx] = True
                break
        return clash_mask

    @staticmethod
    def remove_clashing_chains(
        atom_array: AtomArray,
        clash_radius=1.6,
        clash_ratio_threshold=0.3,
    ) -> list[str]:
        """
        Implements AlphaFold3's chain clash removal protocol (SI Chapter 2.5.4).

        Algorithm Steps:
        1. Identify resolved center atoms for clash detection
        2. Build spatial index (KDTree) for efficient neighbor lookup
        3. Iteratively detect and remove worst-offending chains until:
           - No chains exceed clash ratio threshold
           - All remaining chains meet quality criteria

        Removal Priority Hierarchy:
        1. Core chains (from core_indices) get retention priority
        2. Chains with higher clash ratios removed first
        3. For equal clash ratios: smaller chains removed first
        4. For equal size: lexicographically larger chain IDs removed first

        Args:
            atom_array: Input bioassembly structure
            clash_radius: Distance threshold (Å) for atomic clashes (1.6Å ~ 0.4Å less than vdW radii)
            clash_ratio_threshold: Minimum fraction of clashing atoms (0.3 = 30%) to trigger removal

        Returns:
            List of removed chain IDs for auditing purposes
        """
        is_resolved_centre_atom = (
            atom_array.centre_atom_mask == 1
        ) & atom_array.is_resolved

        # eg: 1qzb
        assert is_resolved_centre_atom.sum() > 0, "No resolved center atoms found"

        resolved_centre_atom_array = atom_array[is_resolved_centre_atom]

        kdtree = KDTree(resolved_centre_atom_array.coord)
        query_result = kdtree.query_ball_point(
            resolved_centre_atom_array.coord,
            r=clash_radius,
        )

        # record basic clash info
        chain_id_to_clashed_atom_index = defaultdict(set)
        atom_index_to_num_clashed_chains = {}
        clash_mask = np.zeros(len(resolved_centre_atom_array), dtype=bool)
        for atom_idx, clashes in enumerate(query_result):
            chain_id = resolved_centre_atom_array.chain_id[atom_idx]
            bonded_atoms, _bond_types = resolved_centre_atom_array.bonds.get_bonds(
                atom_idx
            )
            filtered_clashes = set(clashes) - set(bonded_atoms) - set([atom_idx])
            filtered_clashes = set(
                [
                    clash_atom_idx
                    for clash_atom_idx in filtered_clashes
                    if resolved_centre_atom_array.chain_id[clash_atom_idx] != chain_id
                ]
            )
            chain_id_to_clashed_atom_index[chain_id] = (
                chain_id_to_clashed_atom_index[chain_id] | filtered_clashes
            )
            atom_index_to_num_clashed_chains[atom_idx] = len(
                np.unique(resolved_centre_atom_array.chain_id[list(filtered_clashes)])
            )
            clash_mask[atom_idx] = atom_index_to_num_clashed_chains[atom_idx] > 0

        chain_ids = np.unique(resolved_centre_atom_array.chain_id)
        sorted_chain_indices = np.argsort(chain_ids)
        avg_occ_per_chain = []
        clash_ratio_per_chain = []
        atom_num_per_chain = []
        clash_num_per_chain = []

        # Some structure has multiple altloc chains in same space, e.g. 3ok4
        first_altloc_per_chain = [
            sorted(
                np.unique(
                    resolved_centre_atom_array.label_alt_id[
                        resolved_centre_atom_array.chain_id == chain_id
                    ]
                )
            )[0]
            for chain_id in chain_ids
        ]

        # ['A', '.', '.', 'C', 'B'] -> [1 0 0 3 2]
        ranked_altloc_per_chain = np.searchsorted(
            np.unique(first_altloc_per_chain), first_altloc_per_chain
        )

        chain_id_to_index = {}
        for i, chain_id in enumerate(chain_ids):
            chain_mask = resolved_centre_atom_array.chain_id == chain_id
            atom_num = np.sum(chain_mask)
            clash_num = np.sum(clash_mask[chain_mask])
            clash_ratio = clash_num / atom_num
            avg_occ = np.mean(resolved_centre_atom_array.occupancy[chain_mask])

            avg_occ_per_chain.append(avg_occ)
            atom_num_per_chain.append(atom_num)
            clash_num_per_chain.append(clash_num)
            clash_ratio_per_chain.append(clash_ratio)
            chain_id_to_index[chain_id] = i

        in_core_clashing_chain_ids = [
            1 if "." not in chain_id else 0 for chain_id in chain_ids
        ]

        removed_chain_ids = set()
        last_removed_chain_id = None

        for _cycle in range(len(chain_ids)):
            if last_removed_chain_id is not None:
                chain_i = chain_id_to_index[last_removed_chain_id]
                clashing_atoms = chain_id_to_clashed_atom_index[last_removed_chain_id]

                for atom_idx in clashing_atoms:
                    chain_id_j = resolved_centre_atom_array.chain_id[atom_idx]

                    if chain_id_j in removed_chain_ids:
                        continue

                    chain_j = chain_id_to_index[chain_id_j]
                    assert clash_ratio_per_chain[chain_j] > 0

                    # update atom_index_to_num_clashed_chains after removing a chain
                    atom_index_to_num_clashed_chains[atom_idx] -= 1
                    assert atom_index_to_num_clashed_chains[atom_idx] >= 0

                    if atom_index_to_num_clashed_chains[atom_idx] == 0:
                        # when num_clashed_chains for some atom reach zero, substract 1 from clash_num_per_chain
                        clash_num_per_chain[chain_j] -= 1
                        assert clash_num_per_chain[chain_j] >= 0

                        # update clash_ratio_per_chain[chain_j]
                        clash_ratio_per_chain[chain_j] = (
                            clash_num_per_chain[chain_j] / atom_num_per_chain[chain_j]
                        )

                # set clash_ratio_per_chain for last_removed_chain to 0
                clash_ratio_per_chain[chain_i] = 0

            has_clash_chain = np.any(
                [
                    clash_ratio > clash_ratio_threshold
                    for clash_ratio in clash_ratio_per_chain
                ]
            )

            if not has_clash_chain:
                # no clash over threshold in all chains
                break

            # filter out non clashing chains and removed chains
            filtered_chains = [
                item
                for item in zip(
                    in_core_clashing_chain_ids,
                    avg_occ_per_chain,
                    ranked_altloc_per_chain,
                    clash_ratio_per_chain,
                    atom_num_per_chain,
                    sorted_chain_indices,
                    chain_ids,
                )
                if item[3] > clash_ratio_threshold
            ]

            # Sort removal candidates by priority rules
            remove_priority = sorted(
                filtered_chains,
                key=lambda x: (
                    x[0],
                    x[1],
                    -x[2],
                    -x[3],
                    x[4],
                    -x[5],
                ),
            )

            # Remove highest priority candidate
            last_removed_chain_id = remove_priority[0][6]
            removed_chain_ids.add(last_removed_chain_id)
        else:
            raise ValueError(
                "Find the number of clashes in chain incidents matching the total number of chains, "
                "indicating a likely record issue with the removed_chain_ids."
            )

        return list(removed_chain_ids)

    @staticmethod
    def remove_unresolved_chains(atom_array: AtomArray) -> AtomArray:
        """
        Remove chains from a bioassembly object which all atoms are not resolved.

        Args:
            atom_array (AtomArray): Biotite AtomArray Object of a bioassembly.

        Returns:
            AtomArray: An AtomArray object with unresolved chains removed.
        """
        valid_chain_id = []
        for chain_id in np.unique(atom_array.chain_id):
            resolved = atom_array.is_resolved[atom_array.chain_id == chain_id]
            if np.any(resolved):
                valid_chain_id.append(chain_id)

        atom_array = atom_array[np.isin(atom_array.chain_id, valid_chain_id)]
        return atom_array

    @staticmethod
    def remove_asymmetric_polymer_ligand_bonds(
        atom_array: AtomArray, entity_poly_type: dict[str, str]
    ) -> AtomArray:
        """
        Remove asymmetric polymer ligand bonds (including protein-protein bond, like disulfide bond).

        AF3 SI 5.1 Structure filters
        Bonds for structures with homomeric subcomplexes lacking the corresponding homomeric symmetry are also removed
        e.g. 4MZ in 1keq

        - if a certain bonded ligand only exists for some of the symmetric copies, but not for all,
        we remove the corresponding bond information from the input.
        In consequence the model has to learn to infer these bonds by itself.

        Args:
            atom_array (AtomArray): input atom array
            entity_poly_type (dict): entity_poly_type dict from MMCIFParser

        Returns:
            AtomArray: output atom array with asymmetric polymer ligand bonds removed.
        """
        # get inter chain bonds
        bonds = atom_array.bonds.as_array()[:, :2]
        inter_chain_bonds_all = bonds[
            atom_array.chain_id[bonds[:, 0]] != atom_array.chain_id[bonds[:, 1]]
        ]

        i_is_polymer_mask = np.isin(
            atom_array.label_entity_id[inter_chain_bonds_all[:, 0]],
            list(entity_poly_type.keys()),
        )
        j_is_polymer_mask = np.isin(
            atom_array.label_entity_id[inter_chain_bonds_all[:, 1]],
            list(entity_poly_type.keys()),
        )

        # filter to at least one polymer
        inter_chain_bonds = inter_chain_bonds_all[i_is_polymer_mask | j_is_polymer_mask]

        if len(inter_chain_bonds) == 0:
            # no inter chain bonds found
            return atom_array

        # (sorted_entity_key, entity_id): chain_ids
        bonded_entity_to_chains = defaultdict(set)

        # (bond_site, bonded_entity_id): chain_ids
        bond_site_and_paired_entity_to_chains = defaultdict(set)
        for i, j in inter_chain_bonds:
            chain_i = atom_array.chain_id[i]
            chain_j = atom_array.chain_id[j]
            entity_i = atom_array.label_entity_id[i]
            entity_j = atom_array.label_entity_id[j]
            sorted_entity_key = tuple(sorted([entity_i, entity_j]))
            bonded_entity_to_chains[(sorted_entity_key, entity_i)].add(chain_i)
            bonded_entity_to_chains[(sorted_entity_key, entity_j)].add(chain_j)

            bond_site_i = (
                atom_array.label_entity_id[i],
                atom_array.res_id[i],
                atom_array.res_name[i],
                atom_array.atom_name[i],
            )
            bond_site_j = (
                atom_array.label_entity_id[j],
                atom_array.res_id[j],
                atom_array.res_name[j],
                atom_array.atom_name[j],
            )
            bond_site_and_paired_entity_to_chains[(bond_site_i, entity_j)].add(chain_j)
            bond_site_and_paired_entity_to_chains[(bond_site_j, entity_i)].add(chain_i)

        # (sorted_entity_key, entity_id): num_copies
        bonded_entity_to_num_copies = {
            entity_key: len(chain_ids)
            for entity_key, chain_ids in bonded_entity_to_chains.items()
        }

        # (bond_site, bonded_entity_id):: num_copies of bonded
        bond_site_and_paired_entity_to_num_copies = {
            bond_site_and_entity_key: len(chain_ids)
            for bond_site_and_entity_key, chain_ids in bond_site_and_paired_entity_to_chains.items()
        }

        # find asymmetric bonds
        asymmetric_bonds = set()
        for i, j in inter_chain_bonds:
            entity_i = atom_array.label_entity_id[i]
            entity_j = atom_array.label_entity_id[j]
            i_copies = bonded_entity_to_num_copies[
                (tuple(sorted([entity_i, entity_j])), entity_i)
            ]
            j_copies = bonded_entity_to_num_copies[
                (tuple(sorted([entity_i, entity_j])), entity_j)
            ]
            num_copies = min(i_copies, j_copies)

            if entity_i not in entity_poly_type:
                # entity_i is not polymer
                polymer_atom = j
                ligand_atom = i
            else:
                polymer_atom = i
                ligand_atom = j

            polymer_bond_site = (
                atom_array.label_entity_id[polymer_atom],
                atom_array.res_id[polymer_atom],
                atom_array.res_name[polymer_atom],
                atom_array.atom_name[polymer_atom],
            )
            bonded_ligand_num = bond_site_and_paired_entity_to_num_copies[
                (polymer_bond_site, atom_array.label_entity_id[ligand_atom])
            ]

            if bonded_ligand_num != num_copies:
                # asymmetric bond
                asymmetric_bonds.add((i, j))

        for i, j in asymmetric_bonds:
            atom_array.bonds.remove_bond(i, j)
        return atom_array

    @staticmethod
    def remove_ligand_absent_atoms(atom_array: AtomArray) -> AtomArray:
        """
        Remove ligand atoms absent from the input atom array.
        For each ligand type (by res_name), retain only atoms present in at least one ligand copy in the input.
        Args:
            atom_array (AtomArray): Input atom array.
        Returns:
            AtomArray: Output atom array with absent ligand atoms removed.
        """
        res_to_atom_name = defaultdict(set)
        lig_indices = np.flatnonzero(atom_array.label_seq_id == ".")
        lig_res_name = atom_array.res_name[lig_indices]
        lig_atom_name = atom_array.atom_name[lig_indices]
        lig_is_resolved = atom_array.is_resolved[lig_indices]

        for res_name, atom_name in zip(
            lig_res_name[lig_is_resolved], lig_atom_name[lig_is_resolved]
        ):
            res_to_atom_name[res_name].add(atom_name)

        remove_indices = []
        for idx, res_name, atom_name in zip(lig_indices, lig_res_name, lig_atom_name):
            if atom_name not in res_to_atom_name[res_name]:
                remove_indices.append(idx)

        keep_mask = np.ones(len(atom_array), dtype=bool)
        keep_mask[remove_indices] = False

        return atom_array[keep_mask]

    @staticmethod
    def remove_ligand_unresolved_leaving_atoms(atom_array: AtomArray) -> AtomArray:
        """
        For a ligand involved in covalent bonding, remove the unresolved leaving atoms from the covalently-bonded central atom,
        regardless of whether they are marked as leaving atoms in the CCD (Chemical Component Dictionary).

        Args:
            atom_array (AtomArray): Input atom array.
        Returns:
            AtomArray: Output atom array with unresolved leaving atoms removed.
        """
        inter_residue_bonds = get_inter_residue_bonds(atom_array)
        bonds = atom_array.bonds.copy()

        removed_indices = []
        for centre_idx in np.unique(inter_residue_bonds):
            if atom_array.label_seq_id[centre_idx] != ".":  # not ligand
                continue
            neighbors, _ = bonds.get_bonds(centre_idx)
            bonds.remove_bonds_to(centre_idx)
            if atom_array.is_resolved[neighbors].sum() == 1:
                # only one neighbor is resolved, means one atom is resolved, eg: 4MZ in 1keq
                # do not remove leaving atoms
                continue

            for n in neighbors:
                if not atom_array.is_resolved[n]:
                    group_idx = struc.find_connected(bonds, n)
                    if atom_array.is_resolved[group_idx].any():
                        continue
                    if (
                        atom_array.chain_id[group_idx]
                        != atom_array.chain_id[centre_idx]
                    ).any():
                        continue
                    if (
                        atom_array.res_id[group_idx] != atom_array.res_id[centre_idx]
                    ).any():
                        continue
                    # only remove unresolved atoms in same residue with centre atom
                    removed_indices.extend(group_idx)

        keep_mask = np.ones(len(atom_array), dtype=bool)
        keep_mask[removed_indices] = False

        return atom_array[keep_mask]

    @staticmethod
    def get_rep_chains_for_too_many_chains(
        atom_array: AtomArray, radius: float = 5.0
    ) -> list[str]:
        """
        Get the representative chains from the given atom array.

        Args:
            atom_array (AtomArray): Biotite AtomArray Object of a Bioassembly.
            radius (float): radius to find interface chains.

        Returns:
            list: A list of chain IDs that are considered representative chains.
        """
        core_chains = set()
        sele_entities = set()

        all_entities = np.unique(atom_array.label_entity_id)

        chains, chain_atom_counts = np.unique(atom_array.chain_id, return_counts=True)
        sort_idx_by_atom_cnt = np.argsort(chain_atom_counts)[::-1]
        sorted_chains = chains[sort_idx_by_atom_cnt]

        cell_list = struc.CellList(
            atom_array, cell_size=radius, selection=atom_array.is_resolved
        )

        for chain in sorted_chains:
            chain_mask = (atom_array.chain_id == chain) & atom_array.is_resolved
            if not np.any(chain_mask):
                # none of atom are resolved
                continue
            entity_id = atom_array.label_entity_id[chain_mask][0]
            if entity_id in sele_entities:
                continue

            coords = atom_array.coord[chain_mask]

            # Get atom indices from the current cell and the eight surrounding cells.
            neighbors_ids_2d = cell_list.get_atoms_in_cells(
                coords,
                cell_radius=1,
            )
            neighbors_ids = np.unique(neighbors_ids_2d)

            neighbors_entity_ids = np.unique(atom_array.label_entity_id[neighbors_ids])
            sele_entities |= set(neighbors_entity_ids)

            neighbors_chain_ids = set(atom_array.chain_id[neighbors_ids])
            core_chains |= neighbors_chain_ids

            if sele_entities == set(all_entities):
                # Stop when already sele all entities
                break
        return list(core_chains)

    @staticmethod
    def _filter_altloc_by_local_largest(atom_array: AtomArray) -> np.ndarray:
        """
        Select alternate conformations with highest average occupancy in contiguous
        residue groups. Handles multi-residue conformation blocks that must share
        consistent altloc selection.

        Algorithm:
        1. Group consecutive residues with altlocs into contiguous blocks
        2. For each block, calculate average occupancy per altloc character across all
           residues in the group
        3. Select altloc with highest average occupancy that exists in all group residues
        4. Apply selection consistently across the entire residue group

        Args:
            atom_array: Input structure with potential alternate conformations

        Returns:
            Boolean mask selecting atoms with either:
            - The chosen altloc character from contiguous groups
            - No altloc specification ('.')
        """
        res_starts = get_residue_starts(atom_array, add_exclusive_stop=True)
        res_groups_list = []
        last_res_has_alt = False
        last_chain_id = None
        last_res_id = None

        for start, end in zip(res_starts[:-1], res_starts[1:]):
            chain_id = atom_array.chain_id[start]
            res_id = atom_array.res_id[start]
            altloc = atom_array.label_alt_id[start:end]
            if np.all(altloc == "."):
                last_res_has_alt = False
            elif (chain_id == last_chain_id) and (res_id == last_res_id):
                # They are same res, but get_residue_starts return 2 starts
                # e.g. res_id 168 in 3i0v
                assert last_res_has_alt
                res_groups_list[-1][-1] = (res_groups_list[-1][-1][0], end)
            else:
                if (chain_id == last_chain_id) and last_res_has_alt:
                    res_groups_list[-1].append((start, end))
                else:
                    res_groups_list.append([(start, end)])
                last_res_has_alt = True
            last_chain_id = chain_id
            last_res_id = res_id

        selected_mask = np.ones(len(atom_array), dtype=bool)
        for group in res_groups_list:
            occ_dict = defaultdict(list)
            chars_set = None
            for start, end in group:
                label_alt_id = atom_array.label_alt_id[start:end]

                altloc_chars, char_indices = np.unique(label_alt_id, return_index=True)
                if chars_set is None:
                    chars_set = set(altloc_chars)
                else:
                    chars_set &= set(altloc_chars)

                for altloc_char, idx in zip(altloc_chars, char_indices):
                    # count occ once for each char of a res
                    occupancy = atom_array.occupancy[start:end][idx]
                    occ_dict[altloc_char].append(occupancy)

            alt_and_avg_occ = [
                (altloc_char, np.mean(occ_list))
                for altloc_char, occ_list in occ_dict.items()
            ]
            sorted_altloc_chars = [
                i[0] for i in sorted(alt_and_avg_occ, key=lambda x: x[1], reverse=True)
            ]
            for start, end in group:
                label_alt_id = atom_array.label_alt_id[start:end]
                chosen_char = None
                for i in sorted_altloc_chars:
                    if i == ".":
                        continue
                    elif i not in label_alt_id:
                        continue
                    elif chars_set and i not in chars_set:
                        continue
                    else:
                        chosen_char = i
                        break

                selected_mask[start:end] = (label_alt_id == chosen_char) + (
                    label_alt_id == "."
                )
        return selected_mask

    @staticmethod
    def filter_altloc(
        atom_array: AtomArray, altloc: str = "local_largest"
    ) -> AtomArray:
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
                - "local_largest": Keep the alternate conformation with the largest average occupancy of
                                   a contiguous set of residues with alternate locations.

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

        elif altloc == "local_largest":
            selected_mask = Filter._filter_altloc_by_local_largest(atom_array)
            return atom_array[selected_mask]

        else:
            return atom_array[np.isin(atom_array.label_alt_id, [altloc, "."])]

    @staticmethod
    def remove_dissociation(
        atom_array: AtomArray, all_chain_pairs: list[tuple[str, str]]
    ) -> AtomArray:
        """
        Remove dissociation chains from the given atom array.

        Args:
            atom_array (AtomArray): Biotite AtomArray Object of a Bioassembly.
            all_chain_pairs (list[tuple[str, str]]): A list of chain pairs that represent interfaces.

        Returns:
            AtomArray: An AtomArray object with dissociation chains removed.
        """
        interface_graph = nx.Graph()
        interface_graph.add_edges_from(all_chain_pairs)

        if len(all_chain_pairs) == 0:
            # no interface
            unique_chain_ids = np.unique(atom_array.chain_id)
            if len(unique_chain_ids) == 1:
                return atom_array
            else:
                chain_info = []
                for chain_id in unique_chain_ids:
                    num_resolved_atoms = np.sum(
                        (atom_array.chain_id == chain_id) & atom_array.is_resolved
                    )
                    chain_info.append((chain_id, num_resolved_atoms))

                sorted_chain_info = sorted(chain_info, key=lambda x: x[1], reverse=True)
                largest_chain = sorted_chain_info[0][0]
                return atom_array[atom_array.chain_id == largest_chain]

        else:
            largest_cc = max(nx.connected_components(interface_graph), key=len)
            atom_array = atom_array[np.isin(atom_array.chain_id, list(largest_cc))]
            return atom_array

    @staticmethod
    def remove_too_far_away_chains(
        atom_array: AtomArray, dist_thres: float = 30.0
    ) -> AtomArray:
        """
        Filter out symmetry copies that are too far from the asymmetric unit.

        Args:
            atom_array: Input structure containing both main and symmetry copies
            dist_thres: Maximum allowed distance (Å) from main unit's bounding box

        Returns:
            AtomArray: Filtered structure with distant symmetry copies removed
        """
        resolved_asym_unit_mask = (
            np.char.find(atom_array.chain_id, ".") == -1
        ) & atom_array.is_resolved

        other_resolved_mask = (
            np.char.find(atom_array.chain_id, ".") != -1
        ) & atom_array.is_resolved

        asym_unit_coord = atom_array.coord[resolved_asym_unit_mask]
        xyz_min = np.min(asym_unit_coord, axis=0)
        xyz_max = np.max(asym_unit_coord, axis=0)

        inner_chain_id = list(np.unique(atom_array.chain_id[resolved_asym_unit_mask]))
        for chain_id in np.unique(atom_array.chain_id[other_resolved_mask]):
            chain_mask = (atom_array.chain_id == chain_id) & atom_array.is_resolved
            coords = atom_array.coord[chain_mask]
            if np.any(
                (coords - xyz_min > -dist_thres) | (coords - xyz_max < dist_thres)
            ):
                inner_chain_id.append(chain_id)

        atom_array = atom_array[np.isin(atom_array.chain_id, inner_chain_id)]

        resolved_asym_unit_mask = (
            np.char.find(atom_array.chain_id, ".") == -1
        ) & atom_array.is_resolved

        other_resolved_mask = (
            np.char.find(atom_array.chain_id, ".") != -1
        ) & atom_array.is_resolved

        kdtree = KDTree(atom_array.coord[resolved_asym_unit_mask])
        query_result = kdtree.query_ball_point(
            atom_array.coord[other_resolved_mask],
            r=dist_thres,
        )
        in_dist_mask = [True if q else False for q in query_result]
        selected_chain_ids = np.unique(
            atom_array.chain_id[other_resolved_mask][in_dist_mask]
        )

        atom_array = atom_array[
            (np.isin(atom_array.chain_id, selected_chain_ids))
            | (np.char.find(atom_array.chain_id, ".") == -1)
        ]
        return atom_array
