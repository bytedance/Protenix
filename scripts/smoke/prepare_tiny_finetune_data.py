#!/usr/bin/env python
# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0

"""Prepare a tiny local Protenix fine-tuning dataset for smoke tests."""

import argparse
import csv
import gzip
import json
import os
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from protenix.data.pipeline.data_pipeline import DataPipeline  # noqa: E402
from protenix.utils.file_io import dump_gzip_pickle  # noqa: E402
from protenix.web_service.dependency_url import URL  # noqa: E402


DEFAULT_CIF_CANDIDATES = [
    "examples/2lwu.cif",
    "scripts/msa/data/mmcif/102m.cif",
    "scripts/msa/data/mmcif/1k1a.cif",
    "scripts/msa/data/mmcif/5zyh.cif",
    "scripts/msa/data/mmcif/7zzx.cif",
    "scripts/msa/data/mmcif/8f9h.cif",
]

COMMON_FILES = {
    "ccd_components_file": "components.cif",
    "ccd_components_rdkit_mol_file": "components.cif.rdkit_mol.pkl",
    "pdb_cluster_file": "clusters-by-entity-40.txt",
    "obsolete_release_data_csv": "obsolete_release_date.csv",
}


def default_out_root() -> Path:
    protenix_root = Path(os.environ.get("PROTENIX_ROOT_DIR", str(Path.home())))
    return protenix_root / "tiny_finetune"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--out-root", type=Path, default=default_out_root())
    parser.add_argument("--max-structures", type=int, default=6)
    parser.add_argument("--n-train", type=int, default=4)
    parser.add_argument("--n-val", type=int, default=2)
    parser.add_argument("--max-n-token", type=int, default=256)
    parser.add_argument("--n-cpu", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--extra-cif", action="append", type=Path, default=[])
    return parser.parse_args()


def as_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def local_common_candidates(filename: str) -> List[Path]:
    env_source = os.environ.get("PROTENIX_COMMON_SOURCE_DIR")
    candidates = []
    if env_source:
        candidates.append(Path(env_source) / filename)
    candidates.append(REPO_ROOT / "data" / "protenix_root" / "common" / filename)
    return candidates


def download_if_missing(url_key: str, output_path: Path) -> None:
    if output_path.exists() and output_path.stat().st_size > 0:
        return
    if output_path.exists():
        output_path.unlink()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    for local_path in local_common_candidates(output_path.name):
        if local_path.exists() and local_path.stat().st_size > 0:
            print(f"[tiny-data] copying {local_path} -> {output_path}")
            shutil.copy2(local_path, output_path)
            return
    url = URL[url_key]
    print(f"[tiny-data] downloading {url} -> {output_path}")
    try:
        urllib.request.urlretrieve(url, output_path)
    except Exception as exc:
        curl_cmd = f"curl -L -o '{output_path}' '{url}'"
        raise RuntimeError(
            f"Could not download {output_path} from {url}: {exc}\n"
            f"Download it manually with:\n  {curl_cmd}"
        ) from exc


def ensure_common_files(protenix_root: Path) -> Dict[str, str]:
    common_dir = protenix_root / "common"
    resolved = {}
    for url_key, filename in COMMON_FILES.items():
        path = common_dir / filename
        download_if_missing(url_key, path)
        resolved[url_key] = str(path)
    return resolved


def unique_existing_paths(paths: Iterable[Path]) -> List[Path]:
    out = []
    seen = set()
    for path in paths:
        abs_path = path.resolve()
        if abs_path in seen:
            continue
        seen.add(abs_path)
        out.append(path)
    return out


def candidate_cifs(extra_cifs: Sequence[Path]) -> Tuple[List[Path], List[Dict[str, str]]]:
    candidates = [REPO_ROOT / rel for rel in DEFAULT_CIF_CANDIDATES]
    candidates.extend(extra_cifs)
    existing = []
    skipped = []
    for path in unique_existing_paths(candidates):
        if path.exists():
            existing.append(path)
        else:
            skipped.append({"path": str(path), "error": "missing"})
    return existing, skipped


def copy_cif_to_dataset(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as in_f, open(dst, "wb") as out_f:
            shutil.copyfileobj(in_f, out_f)
    else:
        shutil.copy2(src, dst)


def row_token_count(row: Dict[str, Any]) -> Optional[int]:
    if "num_tokens" not in row or row["num_tokens"] in ("", None):
        return None
    try:
        return int(float(row["num_tokens"]))
    except (TypeError, ValueError):
        return None


def normalize_rows(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, str]]:
    normalized = []
    for row in rows:
        normalized.append({key: as_text(value) for key, value in row.items()})
    return normalized


def write_csv(path: Path, rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    header = list(rows[0].keys())
    for row in rows:
        for key in row.keys():
            if key not in header:
                header.append(key)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        writer.writerows(rows)


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def validate_dataset(out_root: Path, require_summary: bool = True) -> Tuple[bool, str]:
    train_csv = out_root / "indices" / "train.csv"
    val_csv = out_root / "indices" / "val.csv"
    summary = out_root / "prepare_summary.json"
    required_paths = [train_csv, val_csv]
    if require_summary:
        required_paths.append(summary)
    for path in required_paths:
        if not path.exists():
            return False, f"missing {path}"
    try:
        rows = read_csv_rows(train_csv) + read_csv_rows(val_csv)
    except Exception as exc:
        return False, f"could not read CSVs: {exc}"
    if not rows:
        return False, "train/val CSVs are empty"
    for row in rows:
        pdb_id = row.get("pdb_id", "")
        if not pdb_id:
            return False, "row without pdb_id"
        if not (out_root / "bioassembly" / f"{pdb_id}.pkl.gz").exists():
            return False, f"missing bioassembly for {pdb_id}"
    return True, "ok"


def split_rows(
    rows: Sequence[Dict[str, Any]],
    n_train: int,
    n_val: int,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[str], List[str], bool]:
    pdb_ids = []
    for row in rows:
        pdb_id = as_text(row["pdb_id"])
        if pdb_id not in pdb_ids:
            pdb_ids.append(pdb_id)
    if not pdb_ids:
        raise ValueError("No PDB IDs available for splitting.")

    if len(pdb_ids) >= n_train + n_val:
        train_pdbs = pdb_ids[:n_train]
        val_pdbs = pdb_ids[n_train : n_train + n_val]
    elif len(pdb_ids) >= 2:
        train_count = min(max(1, n_train), len(pdb_ids) - 1)
        train_pdbs = pdb_ids[:train_count]
        val_pdbs = pdb_ids[train_count : train_count + n_val]
        if not val_pdbs:
            val_pdbs = [pdb_ids[-1]]
        if len(val_pdbs) < n_val:
            for pdb_id in reversed(train_pdbs):
                if len(val_pdbs) >= n_val:
                    break
                val_pdbs.append(pdb_id)
    else:
        train_pdbs = [pdb_ids[0]]
        val_pdbs = [pdb_ids[0]]

    train_set = set(train_pdbs)
    val_set = set(val_pdbs)
    train_rows = [row for row in rows if as_text(row["pdb_id"]) in train_set]
    val_rows = [row for row in rows if as_text(row["pdb_id"]) in val_set]
    overlap = bool(train_set.intersection(val_set))
    return train_rows, val_rows, train_pdbs, val_pdbs, overlap


def build_dataset(args: argparse.Namespace) -> Dict[str, Any]:
    out_root = args.out_root.resolve()
    protenix_root = Path(
        os.environ.get("PROTENIX_ROOT_DIR", str(out_root.parent))
    ).resolve()

    if args.force and out_root.exists():
        shutil.rmtree(out_root)
    if out_root.exists() and not args.force:
        is_valid, reason = validate_dataset(out_root)
        if is_valid:
            print(f"[tiny-data] existing dataset is valid: {out_root}")
            with open(out_root / "prepare_summary.json") as f:
                return json.load(f)
        print(f"[tiny-data] rebuilding invalid existing dataset: {reason}")
        shutil.rmtree(out_root)

    mmcif_dir = out_root / "mmcif"
    bioassembly_dir = out_root / "bioassembly"
    indices_dir = out_root / "indices"
    for path in [mmcif_dir, bioassembly_dir, indices_dir]:
        path.mkdir(parents=True, exist_ok=True)

    common_files = ensure_common_files(protenix_root)
    cluster_file = Path(common_files["pdb_cluster_file"])

    candidates, skipped_cifs = candidate_cifs(args.extra_cif)
    if not candidates:
        raise RuntimeError("No local mmCIF candidates exist.")

    all_rows = []
    used_cifs = []
    num_bioassemblies = 0
    for cif_path in candidates:
        if num_bioassemblies >= args.max_structures:
            break
        print(f"[tiny-data] parsing {cif_path}")
        try:
            sample_indices_list, bioassembly_dict = DataPipeline.get_data_from_mmcif(
                cif_path,
                cluster_file,
                "WeightedPDB",
            )
        except Exception as exc:
            skipped_cifs.append({"path": str(cif_path), "error": repr(exc)})
            continue
        if not sample_indices_list or not bioassembly_dict:
            skipped_cifs.append(
                {"path": str(cif_path), "error": "no usable rows returned"}
            )
            continue

        pdb_id = as_text(bioassembly_dict["pdb_id"]).lower()
        bioassembly_dict["pdb_id"] = pdb_id
        dump_gzip_pickle(bioassembly_dict, bioassembly_dir / f"{pdb_id}.pkl.gz")
        copy_cif_to_dataset(cif_path, mmcif_dir / f"{pdb_id}.cif")

        for row in sample_indices_list:
            row["pdb_id"] = pdb_id
        all_rows.extend(sample_indices_list)
        used_cifs.append(str(cif_path))
        num_bioassemblies += 1

    if not all_rows:
        raise RuntimeError(f"No usable rows generated. Skipped CIFs: {skipped_cifs}")

    filtered_rows = [
        row
        for row in all_rows
        if row_token_count(row) is None or row_token_count(row) <= args.max_n_token
    ]
    if not filtered_rows:
        raise RuntimeError(
            "All generated rows were filtered by "
            f"--max-n-token={args.max_n_token}. Increase --max-n-token."
        )

    train_rows, val_rows, train_pdbs, val_pdbs, overlap = split_rows(
        filtered_rows,
        args.n_train,
        args.n_val,
    )
    if not train_rows or not val_rows:
        raise RuntimeError("Tiny dataset split produced an empty train or val split.")

    train_csv = indices_dir / "train.csv"
    val_csv = indices_dir / "val.csv"
    write_csv(train_csv, normalize_rows(train_rows))
    write_csv(val_csv, normalize_rows(val_rows))

    is_valid, reason = validate_dataset(out_root, require_summary=False)
    if not is_valid:
        raise RuntimeError(f"Generated dataset failed validation: {reason}")

    summary = {
        "out_root": str(out_root),
        "train_csv": str(train_csv),
        "val_csv": str(val_csv),
        "num_train_rows": len(train_rows),
        "num_val_rows": len(val_rows),
        "num_bioassemblies": num_bioassemblies,
        "used_cifs": used_cifs,
        "skipped_cifs": skipped_cifs,
        "max_n_token": args.max_n_token,
        "train_pdb_ids": train_pdbs,
        "val_pdb_ids": val_pdbs,
        "train_val_overlap": overlap,
        "common_files": common_files,
    }
    summary_path = out_root / "prepare_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(f"[tiny-data] wrote {summary_path}")
    return summary


def main() -> None:
    args = parse_args()
    if args.max_structures <= 0:
        raise ValueError("--max-structures must be positive")
    if args.n_train <= 0 or args.n_val <= 0:
        raise ValueError("--n-train and --n-val must be positive")
    if args.n_cpu != 1:
        print("[tiny-data] --n-cpu is accepted for compatibility; processing is serial.")
    summary = build_dataset(args)
    print(
        "[tiny-data] tiny data OK: "
        f"{summary['num_train_rows']} train rows, "
        f"{summary['num_val_rows']} val rows"
    )


if __name__ == "__main__":
    main()
