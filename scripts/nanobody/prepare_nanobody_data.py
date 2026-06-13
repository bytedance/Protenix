#!/usr/bin/env python3
"""Prepare SAbDab-style nanobody data for Protenix training/evaluation."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import warnings
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


TASK_SPLITS = [
    ("single", "train_pre2022"),
    ("single", "val_2022h1"),
    ("single", "test_2022h2"),
    ("complex", "train_pre2022"),
    ("complex", "val_2022h1"),
    ("complex", "test_2022h2"),
]


def _load_runtime_deps() -> None:
    global pd
    global tqdm
    global annotate_cdrh3
    global CANONICAL_MANIFEST_COLUMNS, PROTENIX_INDEX_COLUMNS
    global ensure_index_columns, filter_indices_for_manifest_row
    global normalize_metadata_csv, split_manifest_filename, write_split_manifests

    import pandas as pd_module
    from tqdm import tqdm as tqdm_func
    from protenix.nanobody.cdr import annotate_cdrh3 as annotate_cdrh3_func
    from protenix.nanobody.manifest import (
        CANONICAL_MANIFEST_COLUMNS as CANONICAL_MANIFEST_COLUMNS_VALUE,
        PROTENIX_INDEX_COLUMNS as PROTENIX_INDEX_COLUMNS_VALUE,
        ensure_index_columns as ensure_index_columns_func,
        filter_indices_for_manifest_row as filter_indices_for_manifest_row_func,
        normalize_metadata_csv as normalize_metadata_csv_func,
        split_manifest_filename as split_manifest_filename_func,
        write_split_manifests as write_split_manifests_func,
    )

    pd = pd_module
    tqdm = tqdm_func
    annotate_cdrh3 = annotate_cdrh3_func
    CANONICAL_MANIFEST_COLUMNS = CANONICAL_MANIFEST_COLUMNS_VALUE
    PROTENIX_INDEX_COLUMNS = PROTENIX_INDEX_COLUMNS_VALUE
    ensure_index_columns = ensure_index_columns_func
    filter_indices_for_manifest_row = filter_indices_for_manifest_row_func
    normalize_metadata_csv = normalize_metadata_csv_func
    split_manifest_filename = split_manifest_filename_func
    write_split_manifests = write_split_manifests_func


def _materialize_mmcif(source: str, target_dir: Path, pdb_id: str) -> tuple[str, str]:
    """Copy mmCIF input or best-effort convert PDB to mmCIF."""

    target_dir.mkdir(parents=True, exist_ok=True)
    if not source:
        return "", "missing structure_path"
    src = Path(source).expanduser()
    if not src.exists():
        return str(src), "structure file not found"

    suffixes = "".join(src.suffixes).lower()
    if suffixes.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        target = target_dir / f"{pdb_id.lower()}{''.join(src.suffixes[-2:]) if suffixes.endswith('.gz') else '.cif'}"
        if src.resolve() != target.resolve():
            shutil.copy2(src, target)
        return str(target), ""

    if suffixes.endswith((".pdb", ".ent")):
        from Bio.PDB import MMCIFIO, PDBParser

        target = target_dir / f"{pdb_id.lower()}.cif"
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(src))
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(target))
        return str(target), "converted from PDB with Bio.PDB; rich mmCIF metadata may be limited"

    return str(src), f"unsupported structure suffix: {suffixes}"


def materialize_structures(manifest: pd.DataFrame, out_root: Path) -> pd.DataFrame:
    manifest = manifest.copy()
    for idx, row in manifest.iterrows():
        task_type = row["task_type"]
        path, note = _materialize_mmcif(
            row["structure_path"],
            out_root / "mmcif" / task_type,
            row["pdb_id"],
        )
        manifest.at[idx, "prepared_mmcif_path"] = path
        if note:
            manifest.at[idx, "notes"] = "; ".join(
                part for part in [str(manifest.at[idx, "notes"]), note] if part
            )
    return manifest


def fill_cdr_annotations(manifest: pd.DataFrame) -> pd.DataFrame:
    manifest = manifest.copy()
    for idx, row in manifest.iterrows():
        annotation = annotate_cdrh3(
            metadata_start=row.get("cdrh3_start"),
            metadata_end=row.get("cdrh3_end"),
            use_optional_numbering=False,
        )
        manifest.at[idx, "cdrh3_start"] = annotation.start or ""
        manifest.at[idx, "cdrh3_end"] = annotation.end or ""
        manifest.at[idx, "cdrh3_source"] = annotation.source
        if annotation.notes and annotation.source != "metadata":
            manifest.at[idx, "notes"] = "; ".join(
                part for part in [str(manifest.at[idx, "notes"]), annotation.notes] if part
            )
    return manifest


def _write_empty_index(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=PROTENIX_INDEX_COLUMNS).to_csv(
        path, index=False, quoting=csv.QUOTE_NONNUMERIC
    )


def _convert_task_split(
    subset: pd.DataFrame,
    task_type: str,
    split_name: str,
    out_root: Path,
    cluster_file: Path | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    from protenix.data.pipeline.data_pipeline import DataPipeline
    from protenix.utils.file_io import dump_gzip_pickle

    output_index = out_root / "indices" / split_manifest_filename(task_type, split_name)
    bio_dir = out_root / "bioassembly" / task_type
    bio_dir.mkdir(parents=True, exist_ok=True)
    output_index.parent.mkdir(parents=True, exist_ok=True)

    if subset.empty:
        _write_empty_index(output_index)
        return pd.DataFrame(columns=PROTENIX_INDEX_COLUMNS), {"path": str(output_index), "n_rows": 0}

    generated_rows: list[pd.DataFrame] = []
    token_counts: dict[str, int] = {}
    for mmcif_path, rows_for_file in subset.groupby("prepared_mmcif_path"):
        if not mmcif_path or not Path(str(mmcif_path)).exists():
            warnings.warn(f"Skipping missing mmCIF for {task_type}/{split_name}: {mmcif_path}")
            continue
        sample_indices, bioassembly = DataPipeline.get_data_from_mmcif(
            mmcif_path,
            cluster_file,
            "WeightedPDB",
        )
        if not sample_indices or not bioassembly:
            warnings.warn(f"Protenix parser produced no rows for {mmcif_path}")
            continue

        pdb_id = str(bioassembly["pdb_id"])
        dump_gzip_pickle(bioassembly, bio_dir / f"{pdb_id}.pkl.gz")
        token_array = bioassembly.get("token_array")
        if token_array is not None:
            try:
                token_counts[pdb_id.upper()] = len(token_array)
            except TypeError:
                pass

        indices_df = pd.DataFrame(sample_indices)
        for _, manifest_row in rows_for_file.iterrows():
            filtered = filter_indices_for_manifest_row(indices_df, manifest_row)
            if filtered.empty:
                warnings.warn(
                    f"No task-relevant Protenix index rows for {manifest_row['pdb_id']} "
                    f"{task_type} {manifest_row['nanobody_chain']} {manifest_row['antigen_chains']}"
                )
                continue
            filtered["release_date"] = manifest_row["release_date"]
            generated_rows.append(filtered)

    if generated_rows:
        index_df = ensure_index_columns(pd.concat(generated_rows, ignore_index=True))
    else:
        index_df = pd.DataFrame(columns=PROTENIX_INDEX_COLUMNS)
    index_df.to_csv(output_index, index=False, quoting=csv.QUOTE_NONNUMERIC)
    return index_df, {
        "path": str(output_index),
        "n_rows": int(len(index_df)),
        "token_counts": token_counts,
    }


def convert_to_protenix_format(
    manifest: pd.DataFrame,
    out_root: Path,
    cluster_file: Path | None,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    manifest = manifest.copy()
    summaries: list[dict[str, Any]] = []
    for task_type, split_name in tqdm(TASK_SPLITS, desc="Converting task splits"):
        subset = manifest[
            (manifest["task_type"] == task_type) & (manifest["split_name"] == split_name)
        ]
        _, summary = _convert_task_split(
            subset=subset,
            task_type=task_type,
            split_name=split_name,
            out_root=out_root,
            cluster_file=cluster_file,
        )
        summary.update({"task_type": task_type, "split_name": split_name})
        summaries.append(summary)
        for pdb_id, n_token in summary.get("token_counts", {}).items():
            mask = manifest["pdb_id"].astype(str).str.upper() == pdb_id
            manifest.loc[mask, "num_tokens"] = n_token
    return manifest, summaries


def make_inference_jsons(manifest: pd.DataFrame, out_root: Path) -> list[Path]:
    from protenix.data.inference.json_maker import cif_to_input_json

    json_dir = out_root / "inference_jsons"
    json_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for task_type, split_name in [
        ("single", "val_2022h1"),
        ("complex", "val_2022h1"),
        ("single", "test_2022h2"),
        ("complex", "test_2022h2"),
    ]:
        subset = manifest[
            (manifest["task_type"] == task_type) & (manifest["split_name"] == split_name)
        ]
        output_name = f"{task_type}_{split_name.replace('val_', 'val_').replace('test_', 'test_')}.json"
        output_name = output_name.replace("val_2022h1", "val_2022h1").replace("test_2022h2", "test_2022h2")
        output_path = json_dir / output_name
        entries: list[dict[str, Any]] = []
        for _, row in subset.iterrows():
            mmcif_path = row["prepared_mmcif_path"]
            if not mmcif_path or not Path(str(mmcif_path)).exists():
                continue
            try:
                entry = cif_to_input_json(
                    str(mmcif_path),
                    sample_name=f"{row['pdb_id']}_{task_type}",
                    save_entity_and_asym_id=True,
                )
                entries.append(entry)
            except Exception as exc:
                warnings.warn(f"Could not create inference JSON for {mmcif_path}: {exc}")
        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(entries, handle, indent=2)
        written.append(output_path)
    return written


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize SAbDab-like nanobody metadata, split by release date, "
            "and convert local structures to Protenix bioassembly/index files."
        )
    )
    parser.add_argument("--metadata", required=True, type=Path, help="SAbDab-style metadata CSV.")
    parser.add_argument("--structure-dir", required=True, type=Path, help="Directory with local PDB/mmCIF files.")
    parser.add_argument("--out-root", required=True, type=Path, help="Output root, usually ${PROTENIX_ROOT_DIR}/nanobody.")
    parser.add_argument("--cluster-file", type=Path, default=None, help="Optional Protenix/PDB cluster file.")
    parser.add_argument("--n-cpu", type=int, default=1, help="Reserved for compatibility; conversion is parser-bound and deterministic.")
    parser.add_argument("--train-cutoff", default="2022-01-01", help="Train split uses release_date before this date.")
    parser.add_argument("--val-cutoff", default="2022-07-01", help="Validation split upper date.")
    parser.add_argument("--test-cutoff", default="2023-01-01", help="Test split upper date.")
    parser.add_argument(
        "--strict-protenix-cutoff",
        action="store_true",
        help="Use 2021-09-30 as the train cutoff and exclude 2021-Q4 records from train.",
    )
    parser.add_argument(
        "--fetch-missing",
        action="store_true",
        help="Accepted for workflow compatibility; fetching is not required and is not performed by default.",
    )
    parser.add_argument(
        "--skip-protenix-conversion",
        action="store_true",
        help="Only write normalized manifests and inference JSONs; skip DataPipeline conversion.",
    )
    parser.add_argument(
        "--skip-inference-jsons",
        action="store_true",
        help="Skip inference JSON generation; useful for manifest-only smoke checks.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    _load_runtime_deps()
    if args.fetch_missing:
        warnings.warn("--fetch-missing is currently a no-op; provide local structures under --structure-dir.")

    out_root = args.out_root.expanduser()
    out_root.mkdir(parents=True, exist_ok=True)

    manifest = normalize_metadata_csv(
        args.metadata,
        structure_dir=args.structure_dir,
        train_cutoff=args.train_cutoff,
        val_cutoff=args.val_cutoff,
        test_cutoff=args.test_cutoff,
        strict_protenix_cutoff=args.strict_protenix_cutoff,
    )
    manifest = materialize_structures(manifest, out_root)
    manifest = fill_cdr_annotations(manifest)

    conversion_summaries: list[dict[str, Any]] = []
    if not args.skip_protenix_conversion:
        manifest, conversion_summaries = convert_to_protenix_format(
            manifest,
            out_root=out_root,
            cluster_file=args.cluster_file,
        )
    else:
        for task_type, split_name in TASK_SPLITS:
            _write_empty_index(out_root / "indices" / split_manifest_filename(task_type, split_name))

    manifest = manifest[CANONICAL_MANIFEST_COLUMNS]
    write_split_manifests(manifest, out_root)
    manifest.to_csv(out_root / "manifests" / "nanobody_all.csv", index=False)
    inference_jsons = [] if args.skip_inference_jsons else make_inference_jsons(manifest, out_root)

    summary = {
        "out_root": str(out_root),
        "n_manifest_rows": int(len(manifest)),
        "conversion": conversion_summaries,
        "inference_jsons": [str(path) for path in inference_jsons],
    }
    with open(out_root / "manifests" / "prepare_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)


if __name__ == "__main__":
    main()
