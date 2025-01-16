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

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray, get_molecule_indices
from scipy.spatial.distance import cdist

from protenix.data.constants import CRYSTALLIZATION_AIDS


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
    def _get_clashing_chains(
        atom_array: AtomArray, chain_ids: list[str]
    ) -> tuple[np.ndarray, list[int]]:
        """
        Calculate the number of atoms clashing with other chains for each chain
        and return a matrix that records the count of clashing atoms.

        Note: if two chains are covalent, they are not considered as clashing.

        Args:
            atom_array (AtomArray): All atoms, including those not resolved.
            chain_ids (list[str]): Unique chain indices of resolved atoms.

        Returns:
            tuple:
                clash_records (numpy.ndarray): Matrix of clashing atom num.
                                               (i, j) means the ratio of i's atom clashed with j's atoms.
                                               Note: (i, j) != (j, i).
                chain_resolved_atom_nums (list[int]): The number of resolved atoms corresponding to each chain ID.
        """
        is_resolved_centre_atom = (
            atom_array.centre_atom_mask == 1
        ) & atom_array.is_resolved
        cell_list = struc.CellList(
            atom_array, cell_size=1.7, selection=is_resolved_centre_atom
        )

        # (i, j) means the ratio of i's atom clashed with j's atoms
        clash_records = np.zeros((len(chain_ids), len(chain_ids)))

        # record the number of resolved atoms for each chain
        chain_resolved_atom_nums = []

        # record covalent relationship between chains
        chains_covalent_dict = {}
        for idx, chain_id_i in enumerate(chain_ids):
            for chain_id_j in chain_ids[idx + 1 :]:
                mol_indices = get_molecule_indices(
                    atom_array[np.isin(atom_array.chain_id, [chain_id_i, chain_id_j])]
                )
                if len(mol_indices) == 1:
                    covalent = 1
                else:
                    covalent = 0
                chains_covalent_dict[(chain_id_i, chain_id_j)] = covalent
                chains_covalent_dict[(chain_id_j, chain_id_i)] = covalent

        for i, chain_id in enumerate(chain_ids):
            coords = atom_array.coord[
                (atom_array.chain_id == chain_id) & is_resolved_centre_atom
            ]
            chain_resolved_atom_nums.append(len(coords))
            chain_atom_ids = np.where(atom_array.chain_id == chain_id)[0]
            chain_atom_ids_set = set(chain_atom_ids) | {-1}

            # Get atom indices from the current cell and the eight surrounding cells.
            neighbors_ids_2d = cell_list.get_atoms_in_cells(coords, cell_radius=1)
            neighbors_ids = np.unique(neighbors_ids_2d)

            # Remove the atom indices of the current chain.
            other_chain_atom_ids = list(set(neighbors_ids) - chain_atom_ids_set)

            if not other_chain_atom_ids:
                continue
            else:
                # Calculate the distance matrix with neighboring atoms.
                other_chain_atom_coords = atom_array.coord[other_chain_atom_ids]
                dist_mat = cdist(coords, other_chain_atom_coords, metric="euclidean")
                clash_mat = dist_mat < 1.6  # change 1.7 to 1.6 for more compatibility
                if np.any(clash_mat):
                    clashed_other_chain_ids = atom_array.chain_id[other_chain_atom_ids]

                    for other_chain_id in set(clashed_other_chain_ids):

                        # two chains covalent with each other
                        if chains_covalent_dict[(chain_id, other_chain_id)]:
                            continue

                        cols = np.where(clashed_other_chain_ids == other_chain_id)[0]

                        # how many i's atoms clashed with j
                        any_atom_clashed = np.any(
                            clash_mat[:, cols].astype(int), axis=1
                        )
                        clashed_atom_num = np.sum(any_atom_clashed.astype(int))

                        if clashed_atom_num > 0:
                            j = chain_ids.index(other_chain_id)
                            clash_records[i][j] += clashed_atom_num
        return clash_records, chain_resolved_atom_nums

    @staticmethod
    def _get_removed_clash_chain_ids(
        clash_records: np.ndarray,
        chain_ids: list[str],
        chain_resolved_atom_nums: list[int],
        core_chain_id: np.ndarray = [],
    ) -> list[str]:
        """
        Perform pairwise comparisons on the chains, and select the chain IDs
        to be deleted according to the clahsing chain rules.

        Args:
            clash_records (numpy.ndarray): Matrix of clashing atom num.
                                           (i, j) means the ratio of i's atom clashed with j's atoms.
                                           Note: (i, j) != (j, i).
            chain_ids (list[str]): Unique chain indices of resolved atoms.
            chain_resolved_atom_nums (list[int]): The number of resolved atoms corresponding to each chain ID.
            core_chain_id (np.ndarray): The chain ID of the core chain.

        Returns:
            list[str]: A list of chain IDs that have been determined for deletion.
        """
        removed_chain_ids = []
        for i in range(len(chain_ids)):
            atom_num_i = chain_resolved_atom_nums[i]
            chain_idx_i = chain_ids[i]

            if chain_idx_i in removed_chain_ids:
                continue

            for j in range(i + 1, len(chain_ids)):
                atom_num_j = chain_resolved_atom_nums[j]
                chain_idx_j = chain_ids[j]

                if chain_idx_j in removed_chain_ids:
                    continue

                clash_num_ij, clash_num_ji = (
                    clash_records[i][j],
                    clash_records[j][i],
                )

                clash_ratio_ij = clash_num_ij / atom_num_i
                clash_ratio_ji = clash_num_ji / atom_num_j

                if clash_ratio_ij <= 0.3 and clash_ratio_ji <= 0.3:
                    # not reaches the threshold
                    continue
                else:
                    # clashing chains
                    if (
                        chain_idx_i in core_chain_id
                        and chain_idx_j not in core_chain_id
                    ):
                        removed_chain_idx = chain_idx_j
                    elif (
                        chain_idx_i not in core_chain_id
                        and chain_idx_j in core_chain_id
                    ):
                        removed_chain_idx = chain_idx_i

                    elif clash_ratio_ij > clash_ratio_ji:
                        removed_chain_idx = chain_idx_i
                    elif clash_ratio_ij < clash_ratio_ji:
                        removed_chain_idx = chain_idx_j
                    else:
                        if atom_num_i < atom_num_j:
                            removed_chain_idx = chain_idx_i
                        elif atom_num_i > atom_num_j:
                            removed_chain_idx = chain_idx_j
                        else:
                            removed_chain_idx = sorted([chain_idx_i, chain_idx_j])[1]

                    removed_chain_ids.append(removed_chain_idx)

                    if removed_chain_idx == chain_idx_i:
                        # chain i already removed
                        break
        return removed_chain_ids

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
            entity_id = atom_array[start].label_entity_id
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
            entity_id = atom_array[start].label_entity_id
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
    def remove_clashing_chains(
        atom_array: AtomArray,
        core_indices: list[int] = None,
    ) -> AtomArray:
        """
        Ref: AlphaFold3 SI Chapter 2.5.4

        Clashing chains are removed.
        Clashing chains are defined as those with >30% of atoms within 1.7 Å of an atom in another chain.
        If two chains are clashing with each other, the chain with the greater percentage of clashing atoms will be removed.
        If the same fraction of atoms are clashing, the chain with fewer total atoms is removed.
        If the chains have the same number of atoms, then the chain with the larger chain id is removed.

        Note: if two chains are covalent, they are not considered as clashing.

        Args:
            atom_array (AtomArray): Biotite AtomArray Object of a Bioassembly.
            core_indices (list[int]): A list of indices for core structures,
                                      where these indices correspond to structures that will be preferentially
                                      retained when pairwise clash chain assessments are performed.

        Returns:
            atom_array (AtomArray): An AtomArray that has been processed through this filter.
            removed_chain_ids (list[str]): A list of chain IDs that have been determined for deletion.
                                           This is to log whether the filter has been utilized.
        """
        chain_ids = np.unique(atom_array.chain_id[atom_array.is_resolved]).tolist()

        if core_indices is not None:
            core_chain_id = np.unique(atom_array.chain_id[core_indices])
        else:
            core_chain_id = np.array([])

        clash_records, chain_resolved_atom_nums = Filter._get_clashing_chains(
            atom_array, chain_ids
        )
        removed_chain_ids = Filter._get_removed_clash_chain_ids(
            clash_records,
            chain_ids,
            chain_resolved_atom_nums,
            core_chain_id=core_chain_id,
        )

        atom_array = atom_array[~np.isin(atom_array.chain_id, removed_chain_ids)]
        return atom_array, removed_chain_ids

    @staticmethod
    def remove_unresolved_mols(atom_array: AtomArray) -> AtomArray:
        """
        Remove molecules from a bioassembly object which all atoms are not resolved.

        Args:
            atom_array (AtomArray): Biotite AtomArray Object of a bioassembly.

        Returns:
            AtomArray: An AtomArray object with unresolved molecules removed.
        """
        valid_mol_id = []
        for mol_id in np.unique(atom_array.mol_id):
            resolved = atom_array.is_resolved[atom_array.mol_id == mol_id]
            if np.any(resolved):
                valid_mol_id.append(mol_id)

        atom_array = atom_array[np.isin(atom_array.mol_id, valid_mol_id)]
        return atom_array

    @staticmethod
    def remove_asymmetric_polymer_ligand_bonds(
        atom_array: AtomArray, entity_poly_type: dict
    ) -> AtomArray:
        """remove asymmetric polymer ligand bonds (including protein-protein bond, like disulfide bond).

        AF3 SI 5.1 Structure filters
        Bonds for structures with homomeric subcomplexes lacking the corresponding homomeric symmetry are also removed
        - e.g. if a certain bonded ligand only exists for some of the symmetric copies, but not for all,
        we remove the corresponding bond information from the input.
        In consequence the model has to learn to infer these bonds by itself.

        Args:
            atom_array (AtomArray): input atom array

        Returns:
            AtomArray: output atom array with asymmetric polymer ligand bonds removed.
        """
        # get inter chain bonds
        inter_chain_bonds = set()
        for i, j, b in atom_array.bonds.as_array():
            if atom_array.chain_id[i] != atom_array.chain_id[j]:
                inter_chain_bonds.add((i, j))

        # get asymmetric polymer ligand bonds
        asymmetric_bonds = set()
        chain_starts = struc.get_chain_starts(atom_array, add_exclusive_stop=False)
        for bond in inter_chain_bonds:

            if bond in asymmetric_bonds:
                continue

            i, j = bond
            atom_i = atom_array[i]
            atom_j = atom_array[j]
            i_is_polymer = atom_i.label_entity_id in entity_poly_type
            j_is_polymer = atom_j.label_entity_id in entity_poly_type
            if i_is_polymer:
                pass
            elif j_is_polymer:
                i, j = j, i
                atom_i, atom_j = atom_j, atom_i
                i_is_polymer, j_is_polymer = j_is_polymer, i_is_polymer
            else:
                # both entity is not polymer
                continue

            # get atom i mask from all entity i copies
            entity_mask_i = atom_array.label_entity_id == atom_i.label_entity_id
            num_copies = np.isin(chain_starts, np.flatnonzero(entity_mask_i)).sum()
            mask_i = (
                entity_mask_i
                & (atom_array.res_id == atom_i.res_id)
                & (atom_array.atom_name == atom_i.atom_name)
            )
            indices_i = np.flatnonzero(mask_i)

            if len(indices_i) != num_copies:
                # not every copy of entity i has atom i.
                asymmetric_bonds.add(bond)
                continue

            # check all atom i in entity i bond to an atom j in entity j.
            target_bonds = []
            for ii in indices_i:
                ii_bonds = [b for b in inter_chain_bonds if ii in b]
                for bond in ii_bonds:
                    jj = bond[1] if ii == bond[0] else bond[0]
                    atom_jj = atom_array[jj]
                    if atom_jj.label_entity_id != atom_j.label_entity_id:
                        continue
                    if atom_jj.res_name != atom_j.res_name:
                        continue
                    if atom_jj.atom_name != atom_j.atom_name:
                        continue
                    if j_is_polymer and atom_jj.res_id != atom_j.res_id:
                        # only for polymer, check res_id
                        continue
                    # found bond (ii, jj) with same enity_id, res_name, atom_name to bond (i,j)
                    target_bonds.append((min(ii, jj), max(ii, jj)))
                    break
            if len(target_bonds) != num_copies:
                asymmetric_bonds |= set(target_bonds)

        for bond in asymmetric_bonds:
            atom_array.bonds.remove_bond(bond[0], bond[1])
        return atom_array
