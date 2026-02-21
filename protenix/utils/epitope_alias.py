from __future__ import annotations

import copy
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Tuple, Union


def _require_int(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"`{field}` must be an int, got {type(value).__name__}")
    return value


def _require_int_list(values: Any, field: str) -> List[int]:
    if not isinstance(values, list) or len(values) == 0:
        raise ValueError(f"`{field}` must be a non-empty list of ints")
    out: List[int] = []
    for i, v in enumerate(values):
        out.append(_require_int(v, f"{field}[{i}]"))
    return out


def _require_mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"`{field}` must be an object/dict")
    return value


def _parse_chain_ref(obj: Any, field: str) -> Tuple[int, int]:
    m = _require_mapping(obj, field)
    entity = _require_int(m.get("entity"), f"{field}.entity")
    copy_id = _require_int(m.get("copy"), f"{field}.copy")
    return entity, copy_id


def _build_pocket_from_epitopes(epitopes: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Convert a user-friendly `epitopes` payload into Protenix's `constraint.pocket` format.

    Supported shapes:
    1) Direct pocket-like:
       {
         "binder_chain": {"entity": 2, "copy": 1},
         "contact_residues": [{"entity": 1, "copy": 1, "position": 45}, ...],
         "max_distance": 8
       }
    2) Antigen+positions shorthand:
       {
         "binder_chain": {"entity": 2, "copy": 1},
         "antigen": {"entity": 1, "copy": 1},
         "positions": [45, 46, 47],
         "max_distance": 8
       }
    """
    binder_entity, binder_copy = _parse_chain_ref(
        epitopes.get("binder_chain"), "epitopes.binder_chain"
    )

    max_distance = epitopes.get("max_distance", 8)
    if isinstance(max_distance, bool) or not isinstance(max_distance, (int, float)):
        raise ValueError("`epitopes.max_distance` must be a number")
    max_distance = float(max_distance)

    if "contact_residues" in epitopes:
        contact_residues_raw = epitopes.get("contact_residues")
        if not isinstance(contact_residues_raw, list) or len(contact_residues_raw) == 0:
            raise ValueError("`epitopes.contact_residues` must be a non-empty list")
        contact_residues: List[Dict[str, Any]] = []
        for i, r in enumerate(contact_residues_raw):
            r_m = _require_mapping(r, f"epitopes.contact_residues[{i}]")
            entity = _require_int(r_m.get("entity"), f"epitopes.contact_residues[{i}].entity")
            copy_id = _require_int(r_m.get("copy"), f"epitopes.contact_residues[{i}].copy")
            position = _require_int(r_m.get("position"), f"epitopes.contact_residues[{i}].position")
            contact_residues.append({"entity": entity, "copy": copy_id, "position": position})
    else:
        antigen_entity, antigen_copy = _parse_chain_ref(
            epitopes.get("antigen"), "epitopes.antigen"
        )
        positions = _require_int_list(epitopes.get("positions"), "epitopes.positions")
        contact_residues = [
            {"entity": antigen_entity, "copy": antigen_copy, "position": p} for p in positions
        ]

    if (binder_entity, binder_copy) in {
        (r["entity"], r["copy"]) for r in contact_residues
    }:
        raise ValueError("epitopes: binder_chain and contact_residues must be on different chains")

    return {
        "binder_chain": {"entity": binder_entity, "copy": binder_copy},
        "contact_residues": contact_residues,
        "max_distance": max_distance,
    }


def apply_epitopes_alias(sample: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Accept an `epitopes` field as an alias and translate it into `constraint.pocket`.

    - If `constraint.pocket` already exists, it is left unchanged.
    - If `epitopes` exists and `constraint.pocket` does not, it is converted to `constraint.pocket`.
    - The `epitopes` field is removed from the returned dict after conversion to avoid ambiguity.
    """
    out: Dict[str, Any] = copy.deepcopy(dict(sample))

    epitopes = out.get("epitopes")
    if epitopes is None:
        return out
    if not isinstance(epitopes, Mapping):
        raise ValueError("`epitopes` must be an object/dict")

    constraint = out.get("constraint")
    if constraint is None:
        constraint = {}
        out["constraint"] = constraint
    if not isinstance(constraint, MutableMapping):
        raise ValueError("`constraint` must be an object/dict when provided")

    # Respect explicit pocket constraints over the alias.
    if isinstance(constraint.get("pocket"), Mapping):
        out.pop("epitopes", None)
        return out

    constraint["pocket"] = _build_pocket_from_epitopes(epitopes)
    out.pop("epitopes", None)
    return out

