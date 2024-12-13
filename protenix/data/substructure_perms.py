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

import itertools
from collections import defaultdict

import numpy as np
from rdkit import Chem


def neutralize_atoms(mol: Chem.Mol):
    pattern = Chem.MolFromSmarts(
        "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!#4!#5!$([*]~[+1,+2,+3,+4])]"
    )
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def recursive_permutation(atom_inds, permutation_list, res):
    def _permute_atom_ind(atom_inds, permutation):
        # atom_inds: list of atom (positional) indices
        # permutation: values to be permutated in the given order
        permute_inds = [i for i, a in enumerate(atom_inds) if a in permutation]
        for i, perm_ind in enumerate(permute_inds):
            atom_inds[perm_ind] = permutation[i]
        return atom_inds

    if len(permutation_list) == 0:
        res.append(atom_inds)
    else:
        current_permutation_list = permutation_list.copy()
        for permutation in current_permutation_list.pop(0):
            atom_inds_permed = _permute_atom_ind(atom_inds.copy(), permutation)
            recursive_permutation(atom_inds_permed, current_permutation_list, res)


def augment_atom_maps_with_conjugate_terminal_groups(
    original_maps, atomic_number_mapping, terminal_group_tuples, MaxMatches=1e6
):
    """
    Augment atom maps from GetSubstructMatches with extra symmetry from confjugated terminal groups.
    Parameters
    --------------
    original_maps:  Tuple(Tuples), all possible atom index mappings, note we require that the mappings should range from 0 to n_heavy_atom-1 (a.k.a. no gap in indexing)
    atomic_number_mapping: dict, mapping from atom (positional) indices to its atomic numbers, for splitting/removing different types of atoms in each terminal group
    terminal_group_tuples: Tuple(Tuples), a group of pair of atoms whose bonds match the SMARTS string. Ex: ((0, 1), (2, 1), (10, 9), (11, 9), (12, 9), (14, 13), (15, 13))
    MaxMatches: int, cutoff for total number of matches (n_original_perm * n_conjugate perm)

    Returns
    --------------
    augmented_maps: Tuple(Tuples) , original_maps augmented by muliplying the permutations induced by terminal_group_tuples.
    """

    def _terminal_atom_cluster_from_pairs(edges):
        graph = defaultdict(set)
        for u, v in edges:
            graph[u].add(v)
            graph[v].add(u)
        return graph

    def _split_sets_by_mapped_values(list_of_sets, mapping):
        result = []
        for s in list_of_sets:
            mapped_sets = {}
            for elem in s:
                mapped_value = mapping.get(elem)
                if mapped_value not in mapped_sets:
                    mapped_sets[mapped_value] = set()
                mapped_sets[mapped_value].add(elem)
            result.extend(mapped_sets.values())
        return result

    # group terminal group tuples with common atom_indices: [{0, 2}, {10, 11, 12}, {14, 15}]
    terminal_atom_clusters = _terminal_atom_cluster_from_pairs(terminal_group_tuples)
    MaxTerminalGroups = max(
        1, int(np.ceil(np.emath.logn(3, MaxMatches / len(original_maps))))
    )
    # if MaxTerminalGroups is less than the total number terminal groups, sample the first {MaxTerminalGroups} groups (to remove randomness)

    perm_groups = sorted(
        [
            atom_inds
            for common_id, atom_inds in terminal_atom_clusters.items()
            if len(atom_inds) > 1
        ]
    )[: min(MaxTerminalGroups, len(terminal_atom_clusters))]

    # within each terminal group, if there are different atom types, split by atom type (if only one left, discard)
    perm_groups = _split_sets_by_mapped_values(perm_groups, atomic_number_mapping)
    perm_groups = [p for p in perm_groups if len(p) > 1]

    # all permutations according to symmetric conjugate terminal atoms: [[(0, 2), (2, 0)], [(10, 11, 12), (10, 12, 11), (11, 10, 12), (11, 12, 10), (12, 10, 11), (12, 11, 10)], [(14, 15), (15, 14)]]
    perm_groups = [sorted(list(itertools.permutations(g))) for g in perm_groups]

    # recursively permute the original mappings
    augmented_maps = []
    for initial_mapping in original_maps:
        recursive_permutation(list(initial_mapping), perm_groups, augmented_maps)

    # Convert to the same data type as in original_maps
    augmented_maps = tuple(tuple(a) for a in augmented_maps)
    # Remove duplicates: original_maps might have already permutated some of the conjugate_terminal group indices
    return tuple(set(augmented_maps))


def _get_substructure_perms(
    mol: Chem.Mol,
    Neutralize: bool = False,
    CheckStereochem: bool = True,
    SymmetrizeConjugatedTerminal: bool = True,
    MaxMatches: int = 512,
) -> np.ndarray:
    """
    Args:
        CheckStereochem: whether to assure stereochem does not change after permutation
        Neutralize: if true, neutralize the mol before computing the permutations
        SymmetrizeConjugatedTerminal: if true, consider symmetrization of conjugated terminal groups
        MaxMatches: int, cutoff for total number of matches

    return shape=[num_perms, num_atoms]
    """
    ori_idx_w_h = []
    for atom in mol.GetAtoms():
        atom.SetProp("ori_idx_w_h", str(atom.GetIdx()))
        ori_idx_w_h.append(atom.GetIdx())

    # Attention !!!
    # Remove Hs; Otherwise, there will be too many matches.
    mol = Chem.RemoveHs(mol)
    if Neutralize:
        mol = neutralize_atoms(mol)

    # Get substructure matches
    base_perms = np.array(
        mol.GetSubstructMatches(mol, uniquify=False, maxMatches=MaxMatches)
    )
    assert len(base_perms) > 0, "no matches found, error"
    # Check stereochem
    if CheckStereochem:
        chem_order = np.array(
            list(Chem.rdmolfiles.CanonicalRankAtoms(mol, breakTies=False))
        )
        perms_mask = (chem_order[base_perms] == chem_order[None]).sum(
            -1
        ) == mol.GetNumAtoms()
        base_perms = base_perms[perms_mask]

    # Add terminal conjugate groups
    sma = "[O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*]"
    patt = Chem.MolFromSmarts(sma)
    terminal_group_tuples = mol.GetSubstructMatches(patt)
    if (
        len(terminal_group_tuples) > 0 and SymmetrizeConjugatedTerminal
    ):  # Only augment if there exist conjugate pairs or if user sets to
        atomic_number_mapping = {
            i: atom.GetAtomicNum() for i, atom in enumerate(mol.GetAtoms())
        }
        base_perms = augment_atom_maps_with_conjugate_terminal_groups(
            tuple(tuple(a) for a in base_perms),
            atomic_number_mapping,
            terminal_group_tuples,
            MaxMatches,
        )
        base_perms = np.array(base_perms)

    if len(base_perms) > MaxMatches:
        base_perms = base_perms[:MaxMatches]

    new_to_ori_idx_map = {}
    ori_to_new_idx_map = {}
    for atom in mol.GetAtoms():
        ori_idx = int(atom.GetProp("ori_idx_w_h"))
        new_idx = atom.GetIdx()
        new_to_ori_idx_map[new_idx] = ori_idx
        ori_to_new_idx_map[ori_idx] = new_idx

    base_perms = np.vectorize(new_to_ori_idx_map.get)(base_perms)
    perms = np.zeros(shape=(base_perms.shape[0], len(ori_idx_w_h)))
    for i in range(len(ori_idx_w_h)):
        if i in ori_to_new_idx_map:
            perms[:, i] = base_perms[:, ori_to_new_idx_map[i]]
        else:
            # The position of the H atom will not be exchanged.
            perms[:, i] = i
    return perms


def get_substructure_perms(
    mol: Chem.Mol,
    CheckStereochem: bool = True,
    SymmetrizeConjugatedTerminal: bool = True,
    MaxMatches: int = 512,
    KeepProtonation: bool = False,
) -> np.ndarray:
    kwargs = {
        "CheckStereochem": CheckStereochem,
        "SymmetrizeConjugatedTerminal": SymmetrizeConjugatedTerminal,
        "MaxMatches": MaxMatches,
    }

    if KeepProtonation:
        perms = _get_substructure_perms(mol, Neutralize=False, **kwargs)
    else:
        # Have to deuplicate permutations across the two protonation states
        perms = np.unique(
            np.row_stack(
                (
                    _get_substructure_perms(mol, Neutralize=False, **kwargs),
                    _get_substructure_perms(mol, Neutralize=True, **kwargs),
                )
            ),
            axis=0,
        )

    nperm = len(perms)
    if nperm > MaxMatches:
        perms = perms[np.random.choice(range(nperm), MaxMatches, replace=False)]
    return perms


def test():
    testcases = [
        "C1=CC=CC=C1",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    ]
    for smiles in testcases:
        print(smiles)
        molecule = Chem.MolFromSmiles(smiles)
        perms = get_substructure_perms(molecule)
        print(perms.shape)
        print(perms.T)


if __name__ == "__main__":
    test()
