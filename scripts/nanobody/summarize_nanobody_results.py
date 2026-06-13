#!/usr/bin/env python3
"""Summarize nanobody experiment summary JSON files."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


SUMMARY_FIELDS = [
    "experiment_id",
    "team",
    "model_name",
    "base_checkpoint",
    "output_checkpoint",
    "train_sets",
    "validation_sets",
    "template_mode",
    "mean_dockq",
    "dockq_success_rate_0_23",
    "median_cdrh3_rmsd_ca",
    "mean_cdrh3_rmsd_ca",
    "notes",
]


def _read(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect nanobody experiment summaries into CSV/JSON.")
    parser.add_argument("--input-root", required=True, type=Path, help="Root containing experiment_summary.json files.")
    parser.add_argument("--out-csv", required=True, type=Path, help="Summary CSV path.")
    parser.add_argument("--out-json", required=True, type=Path, help="Summary JSON path.")
    args = parser.parse_args()

    rows = []
    for path in sorted(args.input_root.glob("**/experiment_summary.json")):
        summary = _read(path)
        rows.append({field: summary.get(field, "") for field in SUMMARY_FIELDS})

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)


if __name__ == "__main__":
    main()
