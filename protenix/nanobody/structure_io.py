"""Small structure IO helpers for nanobody evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from biotite.structure import AtomArray, AtomArrayStack
import biotite.structure.io as strucio


@dataclass(frozen=True)
class CAResidue:
    chain_id: str
    res_id: int
    res_name: str
    coord: np.ndarray


def load_atom_array(path: str | Path) -> AtomArray:
    """Load a PDB/mmCIF structure as a single Biotite AtomArray."""

    atom_array = strucio.load_structure(str(path), model=1)
    if isinstance(atom_array, AtomArrayStack):
        return atom_array[0]
    return atom_array


def chain_ids_with_ca(atom_array: AtomArray) -> list[str]:
    mask = atom_array.atom_name == "CA"
    ids = [str(chain) for chain in np.unique(atom_array.chain_id[mask])]
    return ids


def extract_ca_residues(
    path: str | Path,
    chain_id: str | None = None,
    fallback_first_chain: bool = True,
) -> list[CAResidue]:
    atom_array = load_atom_array(path)
    ca = atom_array[atom_array.atom_name == "CA"]
    if chain_id and np.any(ca.chain_id == chain_id):
        ca = ca[ca.chain_id == chain_id]
    elif chain_id and not fallback_first_chain:
        return []
    elif chain_id and fallback_first_chain:
        chains = chain_ids_with_ca(atom_array)
        if not chains:
            return []
        ca = ca[ca.chain_id == chains[0]]

    residues: list[CAResidue] = []
    for atom in ca:
        residues.append(
            CAResidue(
                chain_id=str(atom.chain_id),
                res_id=int(atom.res_id),
                res_name=str(atom.res_name),
                coord=np.asarray(atom.coord, dtype=np.float64),
            )
        )
    return residues


def match_ca_by_residue_or_order(
    pred: list[CAResidue],
    ref: list[CAResidue],
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Return matched coordinates, preferring shared residue IDs."""

    pred_by_res = {res.res_id: res for res in pred}
    ref_by_res = {res.res_id: res for res in ref}
    common = sorted(set(pred_by_res) & set(ref_by_res))
    if common:
        pred_coords = np.stack([pred_by_res[idx].coord for idx in common])
        ref_coords = np.stack([ref_by_res[idx].coord for idx in common])
        return pred_coords, ref_coords, common

    n = min(len(pred), len(ref))
    if n == 0:
        return np.empty((0, 3)), np.empty((0, 3)), []
    pred_coords = np.stack([res.coord for res in pred[:n]])
    ref_coords = np.stack([res.coord for res in ref[:n]])
    return pred_coords, ref_coords, [res.res_id for res in ref[:n]]
