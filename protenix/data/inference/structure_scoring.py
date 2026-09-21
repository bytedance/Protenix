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

import fnmatch
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

import numpy as np
from biotite.structure import AtomArray, get_chain_starts

from protenix.data.core.filter import Filter
from protenix.data.core.parser import MMCIFParser
from protenix.data.inference.json_maker import atom_array_to_input_json
from protenix.data.utils import pdb_to_cif

SUPPORTED_STRUCTURE_SUFFIXES = {".cif", ".pdb"}
MISSING_ATOM_POLICIES = {"error", "reference", "zero"}


def sanitize_sample_name(name: str, max_length: int = 60) -> str:
    """
    Convert a file stem into a stable sample/output name.
    """
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return (safe or "sample")[:max_length]


def collect_structure_files(
    input_path: str | os.PathLike[str],
    recursive: bool = False,
    globs: Optional[Sequence[str]] = None,
) -> list[Path]:
    """
    Collect PDB/mmCIF structures from a file or directory.
    """
    root = Path(input_path)
    if not root.exists():
        raise RuntimeError(f"input path {input_path} does not exist")
    if root.is_file():
        paths = [root]
    else:
        paths = [path for path in (root.rglob("*") if recursive else root.glob("*"))]

    patterns = [pattern.strip() for pattern in globs or [] if pattern.strip()]
    files = []
    seen = set()
    for path in paths:
        if not path.is_file():
            continue
        if path.suffix.lower() not in SUPPORTED_STRUCTURE_SUFFIXES:
            continue
        if patterns and not any(
            fnmatch.fnmatch(path.name, pattern)
            or fnmatch.fnmatch(path.name.lower(), pattern.lower())
            for pattern in patterns
        ):
            continue
        resolved = path.resolve()
        if resolved not in seen:
            files.append(resolved)
            seen.add(resolved)
    return sorted(files)


def copy_label_asym_ids_to_entity_ids(input_json: dict[str, Any]) -> bool:
    """
    Preserve source chain IDs when every chain has a unique label_asym_id.

    The inference JSON reader already accepts an optional per-entity ``id`` list.
    ``json_maker`` records source ``label_asym_id`` values but does not copy them
    to ``id``; scoring needs that preservation so supplied coordinates can be
    mapped back to the featurized AtomArray.
    """
    pending_ids = []
    all_ids = []
    for _, entity in iter_entity_dicts(input_json):
        label_ids = entity.get("label_asym_id")
        count = int(entity.get("count", 1))
        if label_ids is None or len(label_ids) != count:
            return False
        label_ids = [str(label_id) for label_id in label_ids]
        if len(set(label_ids)) != len(label_ids):
            return False
        pending_ids.append((entity, label_ids))
        all_ids.extend(label_ids)

    if not pending_ids or len(set(all_ids)) != len(all_ids):
        return False

    for entity, label_ids in pending_ids:
        entity["id"] = label_ids
    return True


def iter_entity_dicts(
    input_json: dict[str, Any],
) -> Iterable[tuple[str, dict[str, Any]]]:
    """
    Yield ``(entity_type, entity_dict)`` pairs from a Protenix input JSON sample.
    """
    for sequence in input_json.get("sequences", []):
        if len(sequence) != 1:
            raise RuntimeError(f"invalid sequence entity: {sequence}")
        entity_type, entity = next(iter(sequence.items()))
        yield entity_type, entity


def structure_to_score_input_json(
    input_file: str | os.PathLike[str],
    sample_name: Optional[str] = None,
    assembly_id: Optional[str] = None,
    altloc: str = "first",
    include_discont_poly_poly_bonds: bool = True,
) -> tuple[dict[str, Any], AtomArray, bool]:
    """
    Parse a PDB/mmCIF structure into Protenix JSON plus the normalized AtomArray.
    """
    input_file = Path(input_file)
    if sample_name is None:
        sample_name = sanitize_sample_name(input_file.stem)

    if input_file.suffix.lower() == ".pdb":
        with tempfile.NamedTemporaryFile(suffix=".cif") as tmp:
            pdb_to_cif(str(input_file), tmp.name, entry_id=sample_name)
            return cif_to_score_input_json(
                tmp.name,
                sample_name=sample_name,
                assembly_id=assembly_id,
                altloc=altloc,
                include_discont_poly_poly_bonds=include_discont_poly_poly_bonds,
            )
    if input_file.suffix.lower() == ".cif":
        return cif_to_score_input_json(
            input_file,
            sample_name=sample_name,
            assembly_id=assembly_id,
            altloc=altloc,
            include_discont_poly_poly_bonds=include_discont_poly_poly_bonds,
        )
    raise RuntimeError(f"unsupported structure format: {input_file}")


def cif_to_score_input_json(
    mmcif_file: str | os.PathLike[str],
    sample_name: str,
    assembly_id: Optional[str] = None,
    altloc: str = "first",
    include_discont_poly_poly_bonds: bool = True,
) -> tuple[dict[str, Any], AtomArray, bool]:
    """
    Convert mmCIF to scoring JSON while keeping the source AtomArray in sync.
    """
    parser = MMCIFParser(str(mmcif_file))
    atom_array = parser.get_structure(altloc, model=1, bond_lenth_threshold=None)
    if atom_array is None:
        raise RuntimeError(f"failed to parse atom_site records from {mmcif_file}")

    atom_array = Filter.remove_water(atom_array)
    atom_array = Filter.remove_hydrogens(atom_array)
    atom_array = parser.mse_to_met(atom_array)
    atom_array = Filter.remove_element_X(atom_array)

    if any("DIFFRACTION" in method for method in parser.methods):
        atom_array = Filter.remove_crystallization_aids(
            atom_array, parser.entity_poly_type
        )

    if assembly_id is not None:
        atom_array = parser.expand_assembly(atom_array, assembly_id)

    input_json = atom_array_to_input_json(
        atom_array=atom_array,
        parser=parser,
        assembly_id=assembly_id,
        output_json=None,
        sample_name=sample_name,
        save_entity_and_asym_id=True,
        include_discont_poly_poly_bonds=include_discont_poly_poly_bonds,
    )
    chain_ids_preserved = copy_label_asym_ids_to_entity_ids(input_json)
    atom_array = _filter_atom_array_to_json_entities(atom_array, input_json)
    return input_json, atom_array, chain_ids_preserved


def build_chain_id_map(
    source_atom_array: AtomArray,
    internal_atom_array: AtomArray,
    prefer_identity: bool = True,
) -> dict[str, Any]:
    """
    Map featurized chain IDs back to chain IDs in the parsed source AtomArray.
    """
    source_order = _chain_order(source_atom_array)
    internal_order = _chain_order(internal_atom_array)
    if len(source_order) != len(internal_order):
        raise RuntimeError(
            "source and featurized structures have different chain counts: "
            f"{len(source_order)} vs {len(internal_order)}"
        )
    if len(set(source_order)) != len(source_order):
        raise RuntimeError(f"duplicate source chain IDs: {source_order}")
    if len(set(internal_order)) != len(internal_order):
        raise RuntimeError(f"duplicate internal chain IDs: {internal_order}")

    if prefer_identity and set(source_order) == set(internal_order):
        internal_to_source = {chain_id: chain_id for chain_id in internal_order}
        source_to_internal = {chain_id: chain_id for chain_id in source_order}
    else:
        internal_to_source = dict(zip(internal_order, source_order))
        source_to_internal = dict(zip(source_order, internal_order))

    return {
        "source_chain_order": source_order,
        "internal_chain_order": internal_order,
        "internal_to_source": internal_to_source,
        "source_to_internal": source_to_internal,
    }


def build_source_coord_map(
    source_atom_array: AtomArray,
) -> dict[tuple[str, int, str], np.ndarray]:
    """
    Build a coordinate lookup keyed by normalized chain, residue ID, and atom name.
    """
    source_coord_map: dict[tuple[str, int, str], np.ndarray] = {}
    is_resolved = (
        np.asarray(source_atom_array.is_resolved, dtype=bool)
        if hasattr(source_atom_array, "is_resolved")
        else np.ones(len(source_atom_array), dtype=bool)
    )
    for idx in range(len(source_atom_array)):
        if not is_resolved[idx]:
            continue
        key = atom_lookup_key(source_atom_array, idx)
        if key in source_coord_map:
            raise RuntimeError(
                "duplicate atom key while building score coordinate map: "
                f"{key}. Select a single model/altloc before scoring."
            )
        source_coord_map[key] = np.asarray(
            source_atom_array.coord[idx], dtype=np.float32
        )
    return source_coord_map


def atom_lookup_key(atom_array: AtomArray, idx: int) -> tuple[str, int, str]:
    """
    Return the normalized coordinate lookup key for one atom.
    """
    return (
        str(atom_array.chain_id[idx]),
        int(atom_array.res_id[idx]),
        str(atom_array.atom_name[idx]).strip(),
    )


def map_source_coords_to_internal(
    source_coord_map: dict[tuple[str, int, str], np.ndarray],
    internal_atom_array: AtomArray,
    chain_id_map: dict[str, Any],
    missing_atom_policy: str = "error",
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """
    Reorder source coordinates into the featurized AtomArray atom order.
    """
    if missing_atom_policy not in MISSING_ATOM_POLICIES:
        raise ValueError(
            f"missing_atom_policy must be one of {sorted(MISSING_ATOM_POLICIES)}, "
            f"got {missing_atom_policy}"
        )

    internal_to_source = chain_id_map["internal_to_source"]
    coords = np.zeros((len(internal_atom_array), 3), dtype=np.float32)
    missing_atoms = []
    for idx in range(len(internal_atom_array)):
        internal_chain_id = str(internal_atom_array.chain_id[idx])
        source_chain_id = internal_to_source.get(internal_chain_id, internal_chain_id)
        source_key = (
            source_chain_id,
            int(internal_atom_array.res_id[idx]),
            str(internal_atom_array.atom_name[idx]).strip(),
        )
        if source_key in source_coord_map:
            coords[idx] = source_coord_map[source_key]
            continue

        missing_atoms.append(
            {
                "internal_chain_id": internal_chain_id,
                "source_chain_id": source_chain_id,
                "res_id": int(internal_atom_array.res_id[idx]),
                "res_name": str(internal_atom_array.res_name[idx]).strip(),
                "atom_name": str(internal_atom_array.atom_name[idx]).strip(),
            }
        )
        if missing_atom_policy == "reference":
            coords[idx] = np.asarray(internal_atom_array.coord[idx], dtype=np.float32)
        elif missing_atom_policy == "zero":
            coords[idx] = 0.0

    if missing_atoms and missing_atom_policy == "error":
        preview = ", ".join(
            f"{atom['source_chain_id']}:{atom['res_id']}:{atom['atom_name']}"
            for atom in missing_atoms[:5]
        )
        raise RuntimeError(
            f"{len(missing_atoms)} featurized atoms are missing source coordinates. "
            f"First missing atoms: {preview}. Use missing_atom_policy='reference' "
            "or 'zero' only when fallback coordinates are acceptable."
        )

    return coords, missing_atoms


def _filter_atom_array_to_json_entities(
    atom_array: AtomArray,
    input_json: dict[str, Any],
) -> AtomArray:
    label_entity_ids = {
        str(entity["label_entity_id"])
        for _, entity in iter_entity_dicts(input_json)
        if "label_entity_id" in entity
    }
    if not label_entity_ids:
        return atom_array
    return atom_array[np.isin(atom_array.label_entity_id, list(label_entity_ids))]


def _chain_order(atom_array: AtomArray) -> list[str]:
    chain_starts = get_chain_starts(atom_array, add_exclusive_stop=False)
    return [str(atom_array.chain_id[idx]) for idx in chain_starts]
