#!/usr/bin/env python3
"""Prepare a local OpenFold-style structural subset for Protenix."""

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


def _structure_files(input_dir: Path) -> list[Path]:
    patterns = ["*.cif", "*.mmcif", "*.cif.gz", "*.mmcif.gz"]
    files: list[Path] = []
    for pattern in patterns:
        files.extend(input_dir.glob(pattern))
    return sorted(files)


def _manifest_files(args: argparse.Namespace) -> list[Path]:
    if args.manifest is None:
        return _structure_files(args.input_dir)
    import pandas as pd

    df = pd.read_csv(args.manifest)
    if args.structure_column not in df.columns:
        raise ValueError(f"Manifest is missing --structure-column {args.structure_column!r}")
    if args.max_resolution is not None and "resolution" in df.columns:
        df = df[pd.to_numeric(df["resolution"], errors="coerce") <= args.max_resolution]
    if args.max_release_date and "release_date" in df.columns:
        dates = pd.to_datetime(df["release_date"], errors="coerce")
        df = df[dates < pd.to_datetime(args.max_release_date)]
    if args.min_length is not None and "sequence_length" in df.columns:
        df = df[pd.to_numeric(df["sequence_length"], errors="coerce") >= args.min_length]
    if args.max_length is not None and "sequence_length" in df.columns:
        df = df[pd.to_numeric(df["sequence_length"], errors="coerce") <= args.max_length]
    if args.cluster_ids and "cluster_id" in df.columns:
        allowed = set(args.cluster_ids.split(","))
        df = df[df["cluster_id"].astype(str).isin(allowed)]
    paths = []
    for value in df[args.structure_column].dropna().tolist():
        path = Path(str(value)).expanduser()
        if not path.is_absolute() and args.input_dir is not None:
            path = args.input_dir / path
        paths.append(path)
    return paths


def prepare_subset(args: argparse.Namespace) -> dict[str, Any]:
    import pandas as pd
    from tqdm import tqdm

    from protenix.data.pipeline.data_pipeline import DataPipeline
    from protenix.nanobody.manifest import ensure_index_columns
    from protenix.utils.file_io import dump_gzip_pickle

    out_root = args.out_root.expanduser()
    index_path = out_root / "indices" / "openfold_subset_train.csv"
    bio_dir = out_root / "bioassembly" / "openfold_subset"
    index_path.parent.mkdir(parents=True, exist_ok=True)
    bio_dir.mkdir(parents=True, exist_ok=True)

    dataset = "Distillation" if args.distillation else "WeightedPDB"
    all_rows: list[pd.DataFrame] = []
    failures: list[str] = []
    for mmcif in tqdm(_manifest_files(args), desc="Preparing OpenFold subset"):
        if not mmcif.exists():
            failures.append(f"missing:{mmcif}")
            continue
        sample_indices, bioassembly = DataPipeline.get_data_from_mmcif(
            mmcif,
            args.cluster_file,
            dataset,
        )
        if not sample_indices or not bioassembly:
            failures.append(f"parse_failed:{mmcif}")
            continue
        pdb_id = str(bioassembly["pdb_id"])
        dump_gzip_pickle(bioassembly, bio_dir / f"{pdb_id}.pkl.gz")
        rows = ensure_index_columns(pd.DataFrame(sample_indices))
        if args.max_tokens is not None and "num_tokens" in rows.columns:
            rows = rows[pd.to_numeric(rows["num_tokens"], errors="coerce") <= args.max_tokens]
        all_rows.append(rows)

    if all_rows:
        index_df = ensure_index_columns(pd.concat(all_rows, ignore_index=True))
    else:
        index_df = pd.DataFrame(columns=ensure_index_columns(pd.DataFrame()).columns)
    index_df.to_csv(index_path, index=False, quoting=csv.QUOTE_NONNUMERIC)
    summary = {
        "indices": str(index_path),
        "bioassembly_dir": str(bio_dir),
        "n_rows": int(len(index_df)),
        "distillation": bool(args.distillation),
        "failures": failures,
    }
    with open(out_root / "indices" / "openfold_subset_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare OpenFold-style mmCIF structures for Protenix continued pre-training.")
    parser.add_argument("--input-dir", type=Path, default=Path("."), help="Directory containing local mmCIF files.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional CSV manifest listing structures.")
    parser.add_argument("--structure-column", default="structure_path", help="Manifest column containing structure paths.")
    parser.add_argument("--out-root", type=Path, required=True, help="Output root, usually ${PROTENIX_ROOT_DIR}/nanobody.")
    parser.add_argument("--cluster-file", type=Path, default=None, help="Optional cluster file passed to DataPipeline.")
    parser.add_argument("--distillation", action="store_true", help="Use Protenix DistillationMMCIFParser mode.")
    parser.add_argument("--max-tokens", type=int, default=None, help="Filter generated index rows by num_tokens.")
    parser.add_argument("--max-resolution", type=float, default=None, help="Optional manifest resolution filter.")
    parser.add_argument("--max-release-date", default=None, help="Optional manifest release-date cutoff.")
    parser.add_argument("--min-length", type=int, default=None, help="Optional manifest sequence length lower bound.")
    parser.add_argument("--max-length", type=int, default=None, help="Optional manifest sequence length upper bound.")
    parser.add_argument("--cluster-ids", default=None, help="Comma-separated cluster IDs to keep when manifest has cluster_id.")
    parser.add_argument("--n-cpu", type=int, default=1, help="Reserved for CLI symmetry; conversion is deterministic.")
    return parser


def main() -> None:
    summary = prepare_subset(build_arg_parser().parse_args())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
