#!/usr/bin/env python3
"""Evaluate Protenix nanobody predictions."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


CSV_COLUMNS = [
    "pdb_id",
    "split",
    "task_type",
    "nanobody_chain",
    "antigen_chains",
    "pred_cif",
    "reference_cif",
    "sample_id",
    "sample_mode",
    "ranking_score",
    "plddt",
    "gpde",
    "ptm",
    "dockq",
    "dockq_success_0_23",
    "irmsd",
    "lrmsd",
    "fnat",
    "cdrh3_rmsd_ca",
    "nanobody_rmsd_ca",
    "lddt_ca",
    "metric_status",
    "notes",
]

pd = None


def _load_pandas() -> None:
    global pd
    if pd is None:
        import pandas as pd_module

        pd = pd_module


def _safe_get(row: Any, key: str) -> str:
    value = row.get(key, "")
    if pd.isna(value):
        return ""
    return str(value)


def _sample_id(path: Path) -> str:
    match = re.search(r"sample_([0-9]+)", path.name)
    return match.group(1) if match else ""


def find_prediction_files(pred_dir: Path, pdb_id: str, task_type: str, sample_mode: str) -> list[Path]:
    names = [f"{pdb_id}_{task_type}", pdb_id, pdb_id.lower(), pdb_id.upper()]
    seen: set[Path] = set()
    files: list[Path] = []
    for name in names:
        for path in pred_dir.glob(f"**/{name}_sample_*.cif"):
            if path not in seen:
                seen.add(path)
                files.append(path)
    if not files:
        for path in pred_dir.glob(f"**/*{pdb_id}*sample_*.cif"):
            if path not in seen:
                files.append(path)
    files = sorted(files)
    if sample_mode == "rank0":
        rank0 = [path for path in files if "_sample_0.cif" in path.name]
        return rank0[:1] if rank0 else files[:1]
    return files


def read_confidence(pred_cif: Path) -> dict[str, Any]:
    conf = pred_cif.with_name(pred_cif.name.replace(".cif", ".json").replace("_sample_", "_summary_confidence_sample_"))
    if not conf.exists():
        return {}
    with open(conf, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _base_row(row: pd.Series, pred_cif: Path | None, sample_mode: str) -> dict[str, Any]:
    reference = _safe_get(row, "prepared_mmcif_path") or _safe_get(row, "structure_path")
    out = {column: "" for column in CSV_COLUMNS}
    out.update(
        {
            "pdb_id": _safe_get(row, "pdb_id"),
            "split": _safe_get(row, "split"),
            "task_type": _safe_get(row, "task_type"),
            "nanobody_chain": _safe_get(row, "nanobody_chain"),
            "antigen_chains": _safe_get(row, "antigen_chains"),
            "pred_cif": str(pred_cif) if pred_cif else "",
            "reference_cif": reference,
            "sample_id": _sample_id(pred_cif) if pred_cif else "",
            "sample_mode": sample_mode,
        }
    )
    return out


def _merge_confidence(out: dict[str, Any], confidence: dict[str, Any]) -> None:
    for key in ("ranking_score", "plddt", "gpde", "ptm"):
        value = confidence.get(key, "")
        if isinstance(value, list):
            value = value[0] if value else ""
        out[key] = value


def evaluate_one(row: Any, pred_cif: Path, sample_mode: str, dockq_timeout: int) -> dict[str, Any]:
    out = _base_row(row, pred_cif, sample_mode)
    reference = Path(out["reference_cif"]).expanduser()
    if not reference.exists():
        out["metric_status"] = "missing_reference"
        out["notes"] = f"reference not found: {reference}"
        return out
    if not pred_cif.exists():
        out["metric_status"] = "missing_prediction"
        out["notes"] = f"prediction not found: {pred_cif}"
        return out

    _merge_confidence(out, read_confidence(pred_cif))
    task_type = out["task_type"]
    if task_type == "complex":
        from protenix.nanobody.dockq import compute_dockq

        chains = ";".join(part for part in [out["nanobody_chain"], out["antigen_chains"]] if part)
        result = compute_dockq(
            pred_cif,
            reference,
            model_chains=chains or None,
            native_chains=chains or None,
            timeout=dockq_timeout,
        )
        for key in ("dockq", "dockq_success_0_23", "irmsd", "lrmsd", "fnat", "metric_status", "notes"):
            out[key] = result.get(key, "")
    else:
        from protenix.nanobody.metrics import evaluate_single_chain_prediction

        result = evaluate_single_chain_prediction(
            pred_cif,
            reference,
            out["nanobody_chain"] or None,
            row.get("cdrh3_start"),
            row.get("cdrh3_end"),
        )
        out["cdrh3_rmsd_ca"] = result.get("cdrh3_rmsd_ca", math.nan)
        out["nanobody_rmsd_ca"] = result.get("nanobody_rmsd_ca", math.nan)
        out["metric_status"] = result.get("metric_status", "ok")
        out["notes"] = "; ".join(
            part
            for part in [result.get("notes", ""), f"cdrh3_status={result.get('cdrh3_status', '')}"]
            if part
        )
    return out


def _oracle_select(rows: list[dict[str, Any]], task_type: str) -> list[dict[str, Any]]:
    if not rows:
        return rows
    if task_type == "complex":
        scored = [row for row in rows if row.get("dockq") not in ("", None)]
        if scored:
            return [max(scored, key=lambda row: float(row["dockq"]))]
    scored = []
    for row in rows:
        try:
            value = float(row.get("cdrh3_rmsd_ca", "nan"))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            scored.append((value, row))
    if scored:
        return [min(scored, key=lambda item: item[0])[1]]
    return rows[:1]


def evaluate_manifest(args: argparse.Namespace) -> list[dict[str, Any]]:
    _load_pandas()
    manifest = pd.read_csv(args.manifest)
    rows: list[dict[str, Any]] = []
    for _, manifest_row in manifest.iterrows():
        pdb_id = _safe_get(manifest_row, "pdb_id")
        task_type = _safe_get(manifest_row, "task_type")
        pred_files = find_prediction_files(args.pred_dir, pdb_id, task_type, args.sample_mode)
        if not pred_files:
            missing = _base_row(manifest_row, None, args.sample_mode)
            missing["metric_status"] = "missing_prediction"
            missing["notes"] = f"no prediction CIF found below {args.pred_dir}"
            rows.append(missing)
            continue
        per_target = [
            evaluate_one(manifest_row, pred_cif, args.sample_mode, args.dockq_timeout)
            for pred_cif in pred_files
        ]
        if args.sample_mode == "oracle":
            per_target = _oracle_select(per_target, task_type)
            for row in per_target:
                row["sample_mode"] = "oracle"
        rows.extend(per_target)
    return rows


def write_outputs(rows: list[dict[str, Any]], out_csv: Path, out_json: Path) -> None:
    from protenix.nanobody.metrics import summarize_rows

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in CSV_COLUMNS})
    summary = summarize_rows(rows)
    with open(out_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate nanobody predictions from Protenix inference output.")
    parser.add_argument("--manifest", required=True, type=Path, help="Canonical nanobody manifest CSV.")
    parser.add_argument("--pred-dir", required=True, type=Path, help="Prediction root produced by runner/inference.py.")
    parser.add_argument("--out-csv", required=True, type=Path, help="Per-target metrics CSV.")
    parser.add_argument("--out-json", required=True, type=Path, help="Aggregate summary JSON.")
    parser.add_argument("--sample-mode", choices=["rank0", "all", "oracle"], default="rank0", help="Prediction sample selection mode.")
    parser.add_argument("--dockq-timeout", type=int, default=300, help="DockQ CLI timeout in seconds.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    rows = evaluate_manifest(args)
    write_outputs(rows, args.out_csv, args.out_json)


if __name__ == "__main__":
    main()
