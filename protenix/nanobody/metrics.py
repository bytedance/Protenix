"""Nanobody evaluation metrics."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np


def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Align mobile coordinates onto target coordinates with the Kabsch algorithm."""

    mobile = np.asarray(mobile, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if mobile.shape != target.shape:
        raise ValueError(f"Coordinate shape mismatch: {mobile.shape} vs {target.shape}")
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(f"Expected coordinates with shape [N, 3], got {mobile.shape}")
    if len(mobile) == 0:
        return mobile.copy()

    mobile_center = mobile.mean(axis=0, keepdims=True)
    target_center = target.mean(axis=0, keepdims=True)
    mobile0 = mobile - mobile_center
    target0 = target - target_center
    cov = mobile0.T @ target0
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    return mobile0 @ rot + target_center


def rmsd(mobile: np.ndarray, target: np.ndarray) -> float:
    mobile = np.asarray(mobile, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if mobile.shape != target.shape:
        raise ValueError(f"Coordinate shape mismatch: {mobile.shape} vs {target.shape}")
    if len(mobile) == 0:
        return math.nan
    return float(np.sqrt(np.mean(np.sum((mobile - target) ** 2, axis=-1))))


def aligned_rmsd(mobile: np.ndarray, target: np.ndarray) -> float:
    if len(mobile) == 0:
        return math.nan
    return rmsd(kabsch_align(mobile, target), target)


def region_rmsd_after_alignment(
    mobile: np.ndarray,
    target: np.ndarray,
    region_mask: np.ndarray,
    align_mask: np.ndarray | None = None,
) -> float:
    if align_mask is None:
        align_mask = np.ones(len(mobile), dtype=bool)
    if np.sum(region_mask) == 0 or np.sum(align_mask) == 0:
        return math.nan
    aligned = kabsch_align(mobile[align_mask], target[align_mask])

    # Recompute the same transform for all atoms.
    mobile_center = mobile[align_mask].mean(axis=0, keepdims=True)
    target_center = target[align_mask].mean(axis=0, keepdims=True)
    cov = (mobile[align_mask] - mobile_center).T @ (target[align_mask] - target_center)
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    transformed = (mobile - mobile_center) @ rot + target_center
    _ = aligned  # keep the direct alignment path exercised for shape validation
    return rmsd(transformed[region_mask], target[region_mask])


def evaluate_single_chain_prediction(
    pred_cif: str | Path,
    reference_cif: str | Path,
    nanobody_chain: str | None,
    cdrh3_start: Any = None,
    cdrh3_end: Any = None,
) -> dict[str, Any]:
    """Compute whole-chain and CDR-H3 CA RMSD for one nanobody target."""

    from protenix.nanobody.structure_io import (
        extract_ca_residues,
        match_ca_by_residue_or_order,
    )

    pred_residues = extract_ca_residues(pred_cif, nanobody_chain)
    ref_residues = extract_ca_residues(reference_cif, nanobody_chain)
    pred_coords, ref_coords, residue_ids = match_ca_by_residue_or_order(
        pred_residues, ref_residues
    )

    result: dict[str, Any] = {
        "nanobody_rmsd_ca": math.nan,
        "cdrh3_rmsd_ca": math.nan,
        "cdrh3_status": "missing",
        "metric_status": "ok",
        "notes": "",
    }
    if len(pred_coords) < 3:
        result["metric_status"] = "insufficient_ca"
        result["notes"] = "fewer than three matched CA atoms"
        return result

    result["nanobody_rmsd_ca"] = aligned_rmsd(pred_coords, ref_coords)

    try:
        start = int(float(str(cdrh3_start)))
        end = int(float(str(cdrh3_end)))
    except (TypeError, ValueError):
        return result

    region_mask = np.asarray([start <= int(res_id) <= end for res_id in residue_ids])
    if np.sum(region_mask) == 0:
        result["cdrh3_status"] = "no_matched_cdrh3_ca"
        return result
    result["cdrh3_rmsd_ca"] = region_rmsd_after_alignment(
        pred_coords,
        ref_coords,
        region_mask=region_mask,
        align_mask=np.ones(len(pred_coords), dtype=bool),
    )
    result["cdrh3_status"] = "ok"
    return result


def _finite_values(rows: Iterable[Mapping[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = row.get(key)
        if value is None or value == "":
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(numeric):
            values.append(numeric)
    return values


def mean_or_none(values: list[float]) -> float | None:
    return float(np.mean(values)) if values else None


def median_or_none(values: list[float]) -> float | None:
    return float(np.median(values)) if values else None


def summarize_rows(rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    dockq = _finite_values(rows, "dockq")
    cdrh3 = _finite_values(rows, "cdrh3_rmsd_ca")
    nb = _finite_values(rows, "nanobody_rmsd_ca")
    successful = [row for row in rows if row.get("metric_status") in {"ok", "dockq_unavailable"}]
    dockq_success = [
        row
        for row in rows
        if str(row.get("dockq_success_0_23", "")).lower() in {"true", "1"}
    ]
    return {
        "n_total": len(rows),
        "n_successful": len(successful),
        "n_missing_reference": sum(row.get("metric_status") == "missing_reference" for row in rows),
        "n_missing_prediction": sum(row.get("metric_status") == "missing_prediction" for row in rows),
        "n_dockq_available": len(dockq),
        "mean_dockq": mean_or_none(dockq),
        "median_dockq": median_or_none(dockq),
        "dockq_success_rate_0_23": (
            len(dockq_success) / len(dockq) if dockq else None
        ),
        "mean_cdrh3_rmsd_ca": mean_or_none(cdrh3),
        "median_cdrh3_rmsd_ca": median_or_none(cdrh3),
        "mean_nanobody_rmsd_ca": mean_or_none(nb),
        "median_nanobody_rmsd_ca": median_or_none(nb),
    }
