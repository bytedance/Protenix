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

"""Geometry feature extraction utilities.

This module extracts geometry and stereochemistry constraints used in training
free guidance from reference chemical descriptions (CCD + RDKit) and aligns
them with the atom indices used by the BioTite `AtomArray`.

"""

import copy
import functools
from collections import defaultdict
from itertools import chain, combinations
from typing import List, Optional, Tuple

import numpy as np
import torch
from biotite.structure import AtomArray, get_residue_starts
from rdkit import Chem
from rdkit.Chem.rdchem import BondStereo
from rdkit.Chem.rdDistGeom import GetExperimentalTorsions, GetMoleculeBoundsMatrix
from rdkit.Chem.rdMolTransforms import GetDihedralRad

from protenix.data.constants import STD_RESIDUES
from protenix.data.core.ccd import get_ccd_ref_info, get_component_rdkit_mol
from protenix.utils.logger import get_logger

logger = get_logger(__name__)

RDKIT_GEOMETRY_FEATURES = [
    "pairwise_distance_index",
    "pairwise_distance_upper_bound",
    "pairwise_distance_lower_bound",
    "pairwise_distance_is_bond",
    "pairwise_distance_is_angle",
    "experimental_torsion_index",
    "experimental_torsion_force_constant",
    "experimental_torsion_sign",
    "linear_triple_bond_index",
    "chiral_index",
    "chiral_orientation",
    "stereo_bond_index",
    "stereo_bond_orientation",
    "planar_improper_index",
    "planar_improper_is_carbonyl",
]

GEOMETRY_FEATURES = [
    "interchain_bond_index",
    "symmetric_chain_index",
] + RDKIT_GEOMETRY_FEATURES
METAL_ATOMIC_NUMBERS = frozenset(
    chain(
        (3, 4),  # Li, Be
        range(11, 14),  # Na..Al
        range(19, 32),  # K..Ga
        range(37, 52),  # Rb..Sb
        range(55, 85),  # Cs..Po
        range(87, 119),  # Fr..Og
    )
)

INT_DTYPES_TORCH = [
    torch.uint8,
    torch.int8,
    torch.uint16,
    torch.int16,
    torch.uint32,
    torch.int32,
    torch.uint64,
    torch.int64,
]

# Expected number of atoms per "*_index" feature. Used only to normalize empty
# index tensors to a stable 2D shape (K, 0).
_INDEX_WIDTHS = {
    "interchain_bond_index": 2,
    "symmetric_chain_index": 2,
    "pairwise_distance_index": 2,
    "linear_triple_bond_index": 3,
    "experimental_torsion_index": 4,
    "chiral_index": 4,
    "stereo_bond_index": 4,
    "planar_improper_index": 4,
}


class GeometryFeaturizer:
    """Extract geometry-related features for an AtomArray.

    This class combines two sources of information:
    1) Chain-level topology from the AtomArray (inter-chain bonds and
       symmetry between chains), and
    2) Per-residue geometry constraints computed from CCD reference structures
       via RDKit.

    The per-residue features are first computed in *local RDKit atom indices*
    and then mapped onto the *global atom indices* of atom_array.

    Args:
        atom_array: AtomArray that holds atom records and bond graph.
        exclude_std_residue: If True, skip standard residues (AA/NA) and only
            featurize non-standard components.
        ccd_mols: Optional pre-built RDKit molecules keyed by CCD code.
            When provided, they will be used instead of loading from CCD cache.
    """

    def __init__(
        self,
        atom_array: AtomArray,
        exclude_std_residue: bool = False,
        ccd_mols: Optional[dict[str, Chem.Mol]] = None,
    ):
        self.atom_array = atom_array
        self.exclude_std_residue = exclude_std_residue
        # `lru_cache` requires hashable arguments; convert dict to tuple for caching.
        self.ccd_mols = tuple(ccd_mols.items()) if ccd_mols else None

    def get_features(self):
        """Compute and return all geometry features as torch tensors."""
        features = self.get_residue_geometry_features()
        features.update(self.get_chain_topology_features())
        return self.dict_to_tensor(features)

    def get_chain_topology_features(self):
        """
        Extract chain-level topology features from AtomArray.

        Returned features:
            - `interchain_bond_index`: list of bonded atom index pairs where the
              two atoms belong to different `asym_id_int` (different chains).
            - `symmetric_chain_index`: list of chain-id pairs that share the same
              `entity_id_int`.
        """
        bonds = self.atom_array.bonds
        asym_id_int = self.atom_array.asym_id_int
        entity_id_int = self.atom_array.entity_id_int
        if bonds is None or bonds.get_bond_count() == 0:
            interchain_bond_index = []
        else:
            bond_index = bonds.as_array()[:, :2]
            interchain_mask = (
                asym_id_int[bond_index[:, 0]] != asym_id_int[bond_index[:, 1]]
            )
            interchain_bond_index = bond_index[interchain_mask].tolist()

        symmetric_chain_index = []
        unique_asym_ids, index = np.unique(asym_id_int, return_index=True)
        for asym_id_1, idx_1 in zip(unique_asym_ids, index):
            for asym_id_2, idx_2 in zip(unique_asym_ids, index):
                if asym_id_2 <= asym_id_1:
                    continue
                if entity_id_int[idx_1] == entity_id_int[idx_2]:
                    symmetric_chain_index.append([asym_id_1, asym_id_2])

        return {
            "interchain_bond_index": interchain_bond_index,
            "symmetric_chain_index": symmetric_chain_index,
        }

    @staticmethod
    def dict_to_tensor(feature_dict):
        """Convert a python feature dict into torch tensors.

        Conventions:
        - Keys ending with `_index` are treated as index arrays.
          If the tensor is 2D, it is transposed to shape `(2|3|4, N)`.
        - Floating tensors are cast to `float32`; integer tensors are cast to `int64`.
        """
        for k, v in feature_dict.items():
            if k.endswith("_index"):
                # Index tensors are always int64 and conventionally shaped as (K, N).
                v = torch.as_tensor(v, dtype=torch.int64)
                if v.ndim == 2:
                    v = v.T
                # Normalize empty lists to a stable 2D shape to avoid downstream
                # shape assumptions like idx.shape[1] failing on 1D empty tensors.
                if v.numel() == 0:
                    width = _INDEX_WIDTHS.get(k)
                    if width is not None:
                        v = v.new_empty((width, 0))
            else:
                v = torch.as_tensor(v)
                if v.is_floating_point():
                    v = v.to(torch.float32)
                elif v.dtype in INT_DTYPES_TORCH:
                    v = v.to(torch.int64)
            feature_dict[k] = v
        return feature_dict

    def get_residue_geometry_features(self):
        """Compute residue-level geometry features and map to global atom indices.

        Returns:
            A dict whose keys are `RDKIT_GEOMETRY_FEATURES`.
            All values are Python lists.
        """
        all_features = {feature_name: [] for feature_name in RDKIT_GEOMETRY_FEATURES}
        starts = get_residue_starts(self.atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            res_name = self.atom_array.res_name[start]
            is_hetero = self.atom_array.hetero[start]
            if self.exclude_std_residue and res_name in STD_RESIDUES and not is_hetero:
                continue
            ref_info = get_ccd_ref_info(
                res_name, ccd_mols=self.ccd_mols, return_atomic_number=True
            )
            if not ref_info:
                continue
            # Skip residues containing metal atoms
            if any(
                at_num in METAL_ATOMIC_NUMBERS for at_num in ref_info["atomic_number"]
            ):
                continue
            features = copy.deepcopy(
                get_ccd_geometry_features(res_name, ccd_mols=self.ccd_mols)
            )
            # Map from CCD atom name -> RDKit atom index -> AtomArray global index.
            idx_mol_to_atom_array = {
                ref_info["atom_map"][atom_name]: idx.item()
                for atom_name, idx in zip(
                    self.atom_array.atom_name[start:stop], np.arange(start, stop)
                )
            }
            for feature_name in [
                "pairwise_distance",
                "experimental_torsion",
                "linear_triple_bond",
                "chiral",
                "stereo_bond",
                "planar_improper",
            ]:
                keep_feature_index = []
                atom_index_key = f"{feature_name}_index"
                for feature_index, atom_indices in enumerate(features[atom_index_key]):
                    # Keep only constraints whose atoms exist in this residue.
                    if all([i in idx_mol_to_atom_array for i in atom_indices]):
                        atom_indices = [idx_mol_to_atom_array[i] for i in atom_indices]
                        keep_feature_index.append(feature_index)
                        all_features[atom_index_key].append(atom_indices)
                value_keys = [
                    f
                    for f in RDKIT_GEOMETRY_FEATURES
                    if feature_name in f and f != f"{feature_name}_index"
                ]
                for value_key in value_keys:
                    all_features[value_key].extend(
                        [features[value_key][i] for i in keep_feature_index]
                    )
        return all_features


def build_angle_triples_from_bonds(bond_index: np.ndarray) -> np.ndarray:
    """
    Build all angle triples (i, j, k) from a bond list.

    An angle is defined by two bonds sharing the center atom `j`: (i-j) and (j-k).
    Returns an int array with shape `(n_angle, 3)`.
    """
    if bond_index.size == 0:
        return np.empty((0, 3), dtype=int)

    num_atoms = bond_index.max() + 1
    neighbors = [[] for _ in range(num_atoms)]
    for u, v in bond_index:
        if u == v:
            continue
        neighbors[u].append(v)
        neighbors[v].append(u)

    angles = []
    for j, neighbor in enumerate(neighbors):
        if len(neighbor) < 2:
            continue
        neighbor = sorted(neighbor)
        for i, k in combinations(neighbor, 2):
            angles.append([i, j, k])

    if not angles:
        return np.empty((0, 3), dtype=int)

    return np.asarray(angles, dtype=int)


def extract_pairwise_distance_bounds_from_mol(mol: Chem.Mol):
    """
    Extract RDKit bounds-matrix distance constraints.

    This function flattens constraints into a list of atom pairs.
    """
    # RDKit returns an NxN bounds matrix where:
    # - upper bounds are stored in the upper triangle `bm[i, j]` (i < j)
    # - lower bounds are stored in the lower triangle `bm[j, i]` (i < j)
    bm = GetMoleculeBoundsMatrix(mol)
    n = mol.GetNumAtoms()
    pairwise_distance_index = []
    pairwise_distance_upper_bound = []
    pairwise_distance_lower_bound = []
    pairwise_distance_is_bond = []
    pairwise_distance_is_angle = []

    # Mark whether a pair corresponds to a chemical bond or a bond angle.
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    if len(bonds) == 0:
        bond_index = np.empty((0, 2), dtype=int)
    else:
        bond_index = np.asarray(bonds, dtype=int)

    bond_pairs = {tuple(sorted((int(i), int(j)))) for i, j in bond_index.tolist()}

    angle_triples = build_angle_triples_from_bonds(bond_index)
    angle_pairs = (
        {tuple(sorted((int(i), int(k)))) for i, _, k in angle_triples.tolist()}
        if angle_triples.size > 0
        else set()
    )
    for i in range(n - 1):
        for j in range(i + 1, n):
            upper_bound = float(bm[i, j])
            lower_bound = float(bm[j, i])
            pairwise_distance_index.append([i, j])
            pairwise_distance_upper_bound.append(upper_bound)
            pairwise_distance_lower_bound.append(lower_bound)
            pairwise_distance_is_bond.append(int((i, j) in bond_pairs))
            pairwise_distance_is_angle.append(int((i, j) in angle_pairs))
    return {
        "pairwise_distance_index": pairwise_distance_index,
        "pairwise_distance_upper_bound": pairwise_distance_upper_bound,
        "pairwise_distance_lower_bound": pairwise_distance_lower_bound,
        "pairwise_distance_is_bond": pairwise_distance_is_bond,
        "pairwise_distance_is_angle": pairwise_distance_is_angle,
    }


def extract_experimental_torsion_from_mol(mol: Chem.Mol):
    """Extract RDKit experimental torsions.

    RDKit provides a set of torsions with (atomIndices, V, signs). We flatten
    these into per-torsion lists so downstream code can treat them as batched
    constraints.
    """
    et = GetExperimentalTorsions(mol, useSmallRingTorsions=True)
    index = []
    force_constant = []
    sign = []
    marked_bonds = set()
    for torsion in et:
        atom_indices = list(torsion["atomIndices"])
        index.append(atom_indices)
        force_constant.append(list(torsion["V"]))
        sign.append(list(torsion["signs"]))
        marked_bonds.add(tuple(sorted((atom_indices[1], atom_indices[2]))))
    # Add extra defined torsions for sp2 atoms in rings of size 4, 5, 6
    ring_info = mol.GetRingInfo()
    if ring_info is None:
        return {
            "experimental_torsion_index": index,
            "experimental_torsion_force_constant": force_constant,
            "experimental_torsion_sign": sign,
        }
    for ring in ring_info.AtomRings():
        n = len(ring)
        if 3 < n < 7:
            for a in range(n):
                idx4 = [
                    ring[a],
                    ring[(a + 1) % n],
                    ring[(a + 2) % n],
                    ring[(a + 3) % n],
                ]
                bond_jk = tuple(sorted((idx4[1], idx4[2])))
                if bond_jk in marked_bonds:
                    continue
                atoms = [mol.GetAtomWithIdx(i) for i in idx4]
                if all(
                    atom.GetHybridization() == Chem.HybridizationType.SP2
                    for atom in atoms
                ):
                    index.append(idx4)
                    force_constant.append([0.0, 100.0, 0.0, 0.0, 0.0, 0.0])
                    sign.append([1, -1, 1, 1, 1, 1])
                    marked_bonds.add(bond_jk)

    return {
        "experimental_torsion_index": index,
        "experimental_torsion_force_constant": force_constant,
        "experimental_torsion_sign": sign,
    }


def extract_chiral_dihedral_from_mol(mol: Chem.Mol):
    """Extract chiral-center dihedral signs from the reference conformer.

    For each tetrahedral chiral center, we take combinations of 3 neighbors and
    compute the dihedral angle sign w.r.t. the center atom. This yields a simple
    +/- orientation signal.
    """
    chiral_index = []
    chiral_orientation = []
    assert mol.GetNumConformers() > 0, "mol does not have ref pos"
    conf = mol.GetConformer(0)
    for atom in mol.GetAtoms():
        chiral_tag = atom.GetChiralTag()
        if (chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW) or (
            chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
        ):
            center_idx = atom.GetIdx()
            neighbors = sorted([neighbor.GetIdx() for neighbor in atom.GetNeighbors()])
            if len(neighbors) > 4 or len(neighbors) < 3:
                continue

            for idx in combinations(neighbors, 3):
                dihedral_idx = [*idx, center_idx]
                chiral_index.append(dihedral_idx)
                phi = GetDihedralRad(conf, *dihedral_idx)
                chiral_orientation.append(1.0 if phi >= 0.0 else -1.0)
    return {
        "chiral_index": chiral_index,
        "chiral_orientation": chiral_orientation,
    }


def extract_linear_triple_bond_from_mol(mol: Chem.Mol):
    """Extract triple-bond angle triplets.

    For non-aromatic SP-SP triple bonds, we create angle triplets to encourage
    linear geometry at both ends.
    """
    triple_bond_index = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.TRIPLE:
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        begin_idx = begin.GetIdx()
        end_idx = end.GetIdx()
        if bond.GetIsAromatic() or begin.GetIsAromatic() or end.GetIsAromatic():
            continue
        if not (
            begin.GetHybridization() == Chem.HybridizationType.SP
            and end.GetHybridization() == Chem.HybridizationType.SP
        ):
            continue
        begin_neighbors = sorted(
            [
                neighbor.GetIdx()
                for neighbor in begin.GetNeighbors()
                if neighbor.GetIdx() != end_idx
            ]
        )
        end_neighbors = sorted(
            [
                neighbor.GetIdx()
                for neighbor in end.GetNeighbors()
                if neighbor.GetIdx() != begin_idx
            ]
        )

        for neighbor_idx in begin_neighbors:
            triple_bond_index.append([neighbor_idx, begin_idx, end_idx])
        for neighbor_idx in end_neighbors:
            triple_bond_index.append([begin_idx, end_idx, neighbor_idx])
    return {
        "linear_triple_bond_index": triple_bond_index,
    }


def compute_planar_dihedral_orientation(conf, dihedral_idx):
    """Compute a coarse planarity indicator for a dihedral.

    Returns:
        0.0 if the absolute dihedral angle is close to 0 (planar),
        1.0 if it is closer to pi (still planar but flipped).

    Notes:
        This function intentionally uses a simple threshold at pi/2.
    """
    phi = np.abs(GetDihedralRad(conf, *dihedral_idx))
    return 0.0 if phi < np.pi / 2 else 1.0


def extract_planar_improper_from_mol(mol: Chem.Mol):
    """Extract improper dihedrals to encourage planarity around sp2 centers (C, N, O)."""
    planar_index = []
    is_carbonyl = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ["C", "N", "O"]:
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP2:
            continue

        neighbors = atom.GetNeighbors()
        if len(neighbors) != 3:
            continue

        center_idx = atom.GetIdx()
        nb_indices = sorted([nb.GetIdx() for nb in neighbors])
        nb1, nb2, nb3 = nb_indices

        # Record three quadruples: [nb1, nb2, center, nb3], [nb3, nb1, center, nb2], [nb2, nb3, center, nb1]
        planar_index.append([nb1, nb2, center_idx, nb3])
        planar_index.append([nb3, nb1, center_idx, nb2])
        planar_index.append([nb2, nb3, center_idx, nb1])

        # Check for C=O where both are sp2
        has_carbonyl = False
        for nb in neighbors:
            if (atom.GetSymbol() == "C" and nb.GetSymbol() == "O") and (
                nb.GetHybridization() == Chem.HybridizationType.SP2
            ):
                has_carbonyl = True
                break

        is_carbonyl.extend([float(has_carbonyl)] * 3)

    return {
        "planar_improper_index": planar_index,
        "planar_improper_is_carbonyl": is_carbonyl,
    }


def extract_stereo_bond_from_mol(mol):
    """Extract E/Z stereo bond dihedrals from the reference conformer.

    For bonds with RDKit stereo tags (E or Z), we build one or two dihedral
    definitions using substituents on each side and compute a coarse orientation.
    """
    stereo_bond_index = []
    stereo_bond_orientation = []
    assert mol.GetNumConformers() > 0, "mol does not have ref pos"
    conf = mol.GetConformer(0)

    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if (stereo != BondStereo.STEREOE) and (stereo != BondStereo.STEREOZ):
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        begin_idx = begin.GetIdx()
        end_idx = end.GetIdx()
        begin_neighbors = sorted(
            [
                neighbor.GetIdx()
                for neighbor in begin.GetNeighbors()
                if neighbor.GetIdx() != end_idx
            ]
        )
        end_neighbors = sorted(
            [
                neighbor.GetIdx()
                for neighbor in end.GetNeighbors()
                if neighbor.GetIdx() != begin_idx
            ]
        )

        if len(begin_neighbors) == 0 or len(end_neighbors) == 0:
            continue
        dihedral_idx = [
            begin_neighbors[0],
            begin.GetIdx(),
            end.GetIdx(),
            end_neighbors[0],
        ]
        stereo_bond_index.append(dihedral_idx)
        stereo_bond_orientation.append(
            compute_planar_dihedral_orientation(conf, dihedral_idx)
        )
        if len(begin_neighbors) == 2 and len(end_neighbors) == 2:
            dihedral_idx = [
                begin_neighbors[1],
                begin_idx,
                end_idx,
                end_neighbors[1],
            ]
            stereo_bond_index.append(dihedral_idx)
            stereo_bond_orientation.append(
                compute_planar_dihedral_orientation(conf, dihedral_idx)
            )
    return {
        "stereo_bond_index": stereo_bond_index,
        "stereo_bond_orientation": stereo_bond_orientation,
    }


@functools.lru_cache
def get_ccd_geometry_features(
    ccd_code: str, ccd_mols: Optional[tuple[tuple[str, Chem.Mol]]] = None
):
    """
    Get geometry features for a CCD component.

    The returned dict uses local RDKit atom indices.
    """
    if ccd_mols is None:
        ccd_mols_dict = {}
    else:
        ccd_mols_dict = {k: v for k, v in ccd_mols}

    # `ccd_mols` is passed as a tuple so this function can be cached.
    if ccd_code in ccd_mols_dict:
        mol = copy.deepcopy(ccd_mols_dict[ccd_code])
    else:
        mol = copy.deepcopy(get_component_rdkit_mol(ccd_code))

    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        Chem.GetSymmSSSR(mol)
        features = {
            **extract_pairwise_distance_bounds_from_mol(mol),
            **extract_experimental_torsion_from_mol(mol),
            **extract_chiral_dihedral_from_mol(mol),
            **extract_linear_triple_bond_from_mol(mol),
            **extract_stereo_bond_from_mol(mol),
            **extract_planar_improper_from_mol(mol),
        }
    except Exception:
        logger.warning(
            f"Warning: mol {ccd_code} compute geometry feature from mol failed, return empty features"
        )
        features = {feature_name: [] for feature_name in RDKIT_GEOMETRY_FEATURES}

    return features


def update_features_by_perm(features: dict, perm: np.ndarray) -> dict:
    """
    Resolve feature conflicts caused by symmetric-atom permutations.

    This function loosens or removes constraints that would otherwise over-constrain equivalent atoms.
    """
    # `perm` encodes atom permutations for symmetry-equivalent atoms.
    # For simplicity, we divided symmetry-equivalent atoms into subgroups
    # and then defensively relax constraints across each class.
    groups = compute_equivalence_groups_from_permutations(perm.T.astype(int))
    perm_groups = [g for g in groups if len(g) > 1]
    eq_atoms = set().union(*perm_groups) if perm_groups else set()
    # remove chiral index related to equivalent atoms
    chiral_index = features["chiral_index"]
    if len(chiral_index) > 0:
        remove_centers = set(
            [
                c[3]
                for c in chiral_index
                if any([sum(i in g for i in c[:3]) > 1 for g in perm_groups])
            ]
        )
        keep = [j for j, c in enumerate(chiral_index) if c[3] not in remove_centers]
        features["chiral_index"] = [features["chiral_index"][j] for j in keep]
        features["chiral_orientation"] = [
            features["chiral_orientation"][j] for j in keep
        ]

    # update pairwise distance bound related to equivalent atoms
    pairwise_distance_index = features["pairwise_distance_index"]
    if len(pairwise_distance_index) > 0:
        atom_i, atom_j = (list(col) for col in zip(*pairwise_distance_index))
        (
            features["pairwise_distance_upper_bound"],
            features["pairwise_distance_lower_bound"],
        ) = update_bounds_by_equivalence_groups(
            atom_i,
            atom_j,
            features["pairwise_distance_upper_bound"],
            features["pairwise_distance_lower_bound"],
            groups,
        )

    # remove experimental torsions related to equivalent atoms
    experimental_torsion_index = features["experimental_torsion_index"]
    if len(experimental_torsion_index) > 0:
        keep = [
            j
            for j, c in enumerate(experimental_torsion_index)
            if sum(i in eq_atoms for i in c) == 0
        ]
        features["experimental_torsion_index"] = [
            features["experimental_torsion_index"][j] for j in keep
        ]
        features["experimental_torsion_force_constant"] = [
            features["experimental_torsion_force_constant"][j] for j in keep
        ]
        features["experimental_torsion_sign"] = [
            features["experimental_torsion_sign"][j] for j in keep
        ]

    return features


def update_bounds_by_equivalence_groups(
    atom_idx_i: List[int],
    atom_idx_j: List[int],
    upper_A: List[float],
    lower_B: List[float],
    groups: List[List[int]],
    undirected: bool = True,
) -> Tuple[List[float], List[float]]:
    """Relax bounds for constraints involving symmetry-equivalent atoms.

    For each constraint (i, j), we map atoms to their equivalence-group IDs. All
    constraints sharing the same (group(i), group(j)) are aggregated by taking:
    - upper bound = max over the group-pair
    - lower bound = min over the group-pair
    The aggregated bounds are then assigned back to each original constraint.
    """
    n = len(atom_idx_i)
    if not (len(atom_idx_j) == len(upper_A) == len(lower_B) == n):
        raise ValueError(
            "`atom_idx_i`, `atom_idx_j`, `upper_A`, and `lower_B` must have the same length."
        )

    reps = {atom: group_idx for group_idx, group in enumerate(groups) for atom in group}

    max_upper, min_lower = {}, {}

    for i, j, ua, lb in zip(atom_idx_i, atom_idx_j, upper_A, lower_B):
        ri, rj = reps[int(i)], reps[int(j)]
        if undirected and ri > rj:
            ri, rj = rj, ri
        key = (ri, rj)
        max_upper[key] = ua if key not in max_upper else max(max_upper[key], ua)
        min_lower[key] = lb if key not in min_lower else min(min_lower[key], lb)

    new_upper, new_lower = [], []
    for i, j in zip(atom_idx_i, atom_idx_j):
        ri, rj = reps[int(i)], reps[int(j)]
        if undirected and ri > rj:
            ri, rj = rj, ri
        key = (ri, rj)
        new_upper.append(max_upper[key])
        new_lower.append(min_lower[key])

    return new_upper, new_lower


class DSU:
    """Disjoint Set Union (Union-Find) helper.

    Used by `compute_equivalence_groups_from_permutations()` to compute
    connected components induced by permutation generators.
    """

    def __init__(self, n):
        self.p = list(range(n))
        self.r = [0] * n

    def find(self, x):
        """Find the representative of x with path compression."""
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        """Union the sets containing a and b (union by rank)."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1


def compute_equivalence_groups_from_permutations(perms, n=None, one_based=False):
    """
    Compute equivalence classes induced by a set of permutations.

    Args:
        perms: A list/array of permutations. Each permutation is a length-n sequence
            representing index mapping. Example (0-based): `[2, 0, 1, 3]` means
            0→2, 1→0, 2→1, 3→3.
        n: Number of indices. Defaults to `len(perms[0])`.
        one_based: If True, interpret permutations as 1..n and return 1-based groups.

    Returns:
        A list of sorted equivalence groups. Two indices are equivalent if one can be
        mapped to the other by composing the provided permutations.
    """
    if len(perms) == 0:
        if n is None:
            return []
        return [[i + 1] if one_based else [i] for i in range(n)]

    if n is None:
        n = len(perms[0])

    # Normalize permutations to 0-based and validate bijectivity.
    norm_perms = []
    for pi in perms:
        if len(pi) != n:
            raise ValueError("All permutations must have the same length n.")
        if one_based:
            pi0 = [x - 1 for x in pi]
            if set(pi0) != set(range(n)):
                raise ValueError("Each permutation must be a bijection of 1..n.")
            norm_perms.append(pi0)
        else:
            if set(pi) != set(range(n)):
                raise ValueError("Each permutation must be a bijection of 0..n-1.")
            norm_perms.append(list(pi))

    dsu = DSU(n)
    # Union all edges (i, σ(i)) for each generator permutation σ.
    for pi in norm_perms:
        for i, j in enumerate(pi):
            dsu.union(i, j)

    buckets = defaultdict(list)
    for i in range(n):
        buckets[dsu.find(i)].append(i)

    groups = []
    for comp in buckets.values():
        comp.sort()
        if one_based:
            comp = [x + 1 for x in comp]
        groups.append(comp)
    groups.sort(key=lambda g: (len(g), g))
    return groups
