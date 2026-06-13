#!/usr/bin/env python3
"""Verify the SAbDab-nano temporal dataset bundle and write the dataset card."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set

from sabdab_nano_temporal_common import (
    COMPLEX_TEST_EXPECTED,
    COMPLEX_TRAIN_PAPER_REFERENCE,
    COMPLEX_VAL_PAPER_REFERENCE,
    MANIFEST_FIELDS,
    RCSB_ENTRY_API_TEMPLATE,
    SABDAB_ARCHIVE_URL,
    SABDAB_DOWNLOADER_URL,
    SABDAB_MAIN_URL,
    SABDAB_NANO_URL,
    SABDAB_SUMMARY_URL,
    SINGLE_TEST_EXPECTED,
    TFOLD_GDRIVE_DIRECT_URL,
    TFOLD_GDRIVE_VIEW_URL,
    TFOLD_PAPER_URL,
    TFOLD_REPO_URL,
    TFOLD_WEIYUN_URL,
    clean_value,
    dataset_paths,
    manifest_path,
    read_csv_rows,
    relpath,
    sha256_file,
    validate_summary_tsv,
    validate_tar_gz,
    validate_zip,
    write_json,
)


def _print_step(message: str) -> None:
    print(f"[sabdab-nano-temporal] {message}", flush=True)


def _checksum(path: Path) -> str:
    return sha256_file(path) if path.exists() and path.stat().st_size > 0 else ""


def _load_json(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _manifest_rows(paths, split_name: str, warnings: List[str], errors: List[str]) -> List[Dict[str, str]]:
    path = manifest_path(paths, split_name)
    if not path.exists():
        errors.append(f"Missing split manifest: {relpath(path, paths.output_root)}")
        return []
    rows = read_csv_rows(path)
    fieldnames = set(rows[0].keys()) if rows else set(MANIFEST_FIELDS)
    missing_fields = [field for field in MANIFEST_FIELDS if field not in fieldnames]
    if missing_fields:
        errors.append(f"{path.name} missing required columns: {', '.join(missing_fields)}")
    return rows


def _path_exists(output_root: Path, value: str) -> bool:
    if not value:
        return True
    path = Path(value)
    if not path.is_absolute():
        path = output_root / path
    return path.exists()


def _validate_manifest_rows(paths, split_name: str, rows: Sequence[Mapping[str, str]], warnings: List[str], errors: List[str]) -> None:
    seen: Set[str] = set()
    for idx, row in enumerate(rows, start=2):
        target_id = clean_value(row.get("target_id", ""))
        if not target_id:
            errors.append(f"{split_name}: row {idx} missing target_id")
            continue
        if target_id in seen:
            errors.append(f"{split_name}: duplicate target_id {target_id}")
        seen.add(target_id)

        if clean_value(row.get("light_chain", "")) not in {"", "NA"}:
            errors.append(f"{split_name}: {target_id} has non-empty light_chain")
        if row.get("task") == "complex_nanoag" and not clean_value(row.get("antigen_chains", "")):
            errors.append(f"{split_name}: {target_id} complex row has empty antigen_chains")

        for column in ["raw_pdb_path", "chothia_pdb_path", "imgt_pdb_path", "native_cif_path", "fasta_path", "json_path", "source_summary_path", "source_archive"]:
            value = clean_value(row.get(column, ""))
            if value and not _path_exists(paths.output_root, value):
                errors.append(f"{split_name}: {target_id} points to missing {column}: {value}")

        structure_values = [
            clean_value(row.get("raw_pdb_path", "")),
            clean_value(row.get("chothia_pdb_path", "")),
            clean_value(row.get("imgt_pdb_path", "")),
            clean_value(row.get("native_cif_path", "")),
        ]
        if not any(value and _path_exists(paths.output_root, value) for value in structure_values):
            errors.append(f"{split_name}: {target_id} has no existing structure file path")


def _target_ids(rows: Sequence[Mapping[str, str]]) -> Set[str]:
    return {clean_value(row.get("target_id", "")) for row in rows if clean_value(row.get("target_id", ""))}


def _date_in_range(value: str, start: Optional[dt.date], end: dt.date) -> bool:
    if not value:
        return False
    date_value = dt.date.fromisoformat(value)
    if start is not None and date_value < start:
        return False
    return date_value < end


def _validate_temporal_bounds(split_name: str, rows: Sequence[Mapping[str, str]], warnings: List[str], errors: List[str]) -> None:
    for row in rows:
        target_id = row.get("target_id", "")
        release_date = clean_value(row.get("release_date", ""))
        if split_name.endswith("train_pre2022"):
            ok = _date_in_range(release_date, None, dt.date(2022, 1, 1))
        elif split_name.endswith("val_2022h1"):
            ok = _date_in_range(release_date, dt.date(2022, 1, 1), dt.date(2022, 7, 1))
        elif split_name.endswith("test_2022h2"):
            ok = _date_in_range(release_date, dt.date(2022, 7, 1), dt.date(2023, 1, 1))
        else:
            ok = True
        if not ok:
            warnings.append(f"{split_name}: {target_id} release_date {release_date or 'missing'} outside expected window")


def _write_dataset_card(paths, report: Mapping[str, object]) -> None:
    splits = report.get("splits", {})
    checksums = report.get("checksums", {})
    warnings = report.get("warnings", [])
    counts_lines = []
    if isinstance(splits, dict):
        for split_name in [
            "single_nano_train_pre2022",
            "single_nano_val_2022h1",
            "single_nano_test_2022h2",
            "complex_nanoag_train_pre2022",
            "complex_nanoag_val_2022h1",
            "complex_nanoag_test_2022h2",
        ]:
            item = splits.get(split_name, {})
            if isinstance(item, dict):
                counts_lines.append(f"- `{split_name}`: {item.get('actual_targets', 0)} targets ({item.get('status', 'unknown')})")

    lines = [
        "# SAbDab-nano temporal structure prediction dataset",
        "",
        "## Purpose",
        "This local bundle supports fine-tuning and evaluating nanobody structure prediction workflows without training models or downloading model weights.",
        "",
        "## Tasks",
        "- `single_nano`: one single-domain antibody/VHH chain per target. Bound structures are retained as single-chain nanobody targets and annotated as bound.",
        "- `complex_nanoag`: one nanobody chain plus one or more bound antigen chains per target.",
        "",
        "## Splits",
        "- `train_pre2022`: RCSB initial release date before 2022-01-01.",
        "- `val_2022h1`: RCSB initial release date from 2022-01-01 through 2022-06-30.",
        "- `test_2022h2`: official tFold H2 archive targets from 2022-07-01 through 2022-12-31.",
        "",
        "## Sources",
        f"- SAbDab: {SABDAB_MAIN_URL}",
        f"- SAbDab-nano: {SABDAB_NANO_URL}",
        f"- SAbDab all-structures archive: {SABDAB_ARCHIVE_URL}",
        f"- SAbDab all-summary TSV: {SABDAB_SUMMARY_URL}",
        f"- tFold repository: {TFOLD_REPO_URL}",
        f"- tFold official test archive: {TFOLD_GDRIVE_VIEW_URL}",
        f"- tFold paper: {TFOLD_PAPER_URL}",
        f"- RCSB entry API: {RCSB_ENTRY_API_TEMPLATE}",
        "",
        f"Download/verification date UTC: {report.get('download_date_utc', '')}",
        "",
        "## SHA256",
        f"- `all_structures.zip`: `{checksums.get('all_structures_zip_sha256', '') if isinstance(checksums, dict) else ''}`",
        f"- `sabdab_summary_all.tsv`: `{checksums.get('sabdab_summary_all_sha256', '') if isinstance(checksums, dict) else ''}`",
        f"- `tFold_test_set.tar.gz`: `{checksums.get('tfold_test_set_tar_gz_sha256', '') if isinstance(checksums, dict) else ''}`",
        "",
        "## Counts",
        *counts_lines,
        "",
        "## Caveats",
        "- SAbDab is live and updated regularly, so train/validation counts can vary from historical snapshots.",
        "- H2 test data is taken from the official tFold archive rather than reconstructed from current SAbDab.",
        "- SAbDab summary `date` is deposition metadata; RCSB `initial_release_date` is preferred for temporal splitting.",
        "- Complex train/validation counts may differ from the tFold paper unless the same historical SAbDab snapshot and every author-side filter are reproduced.",
        "",
        "## Suggested Citations",
        f"- SAbDab/SAbDab-nano: {SABDAB_MAIN_URL}",
        f"- tFold: {TFOLD_PAPER_URL}",
    ]
    if warnings:
        lines.extend(["", "## Verification Warnings"])
        lines.extend(f"- {warning}" for warning in warnings)
    (paths.manifests / "DATASET_CARD.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def verify_dataset(output_root: Path) -> Dict[str, object]:
    paths = dataset_paths(output_root)
    paths.manifests.mkdir(parents=True, exist_ok=True)
    warnings: List[str] = []
    errors: List[str] = []

    if not validate_zip(paths.sabdab_archive):
        errors.append(f"SAbDab all-structures archive is missing or invalid: {relpath(paths.sabdab_archive, paths.output_root)}")
    if not validate_summary_tsv(paths.sabdab_summary):
        errors.append(f"SAbDab summary is missing, invalid, or lacks required columns: {relpath(paths.sabdab_summary, paths.output_root)}")
    if not validate_tar_gz(paths.tfold_archive):
        errors.append(f"Official tFold archive is missing or invalid: {relpath(paths.tfold_archive, paths.output_root)}")

    build_summary = _load_json(paths.manifests / "build_summary.json")
    leakage_report = _load_json(paths.manifests / "leakage_report.json")
    split_rows: Dict[str, List[Dict[str, str]]] = {}
    for split_name in [
        "single_nano_train_pre2022",
        "single_nano_val_2022h1",
        "single_nano_test_2022h2",
        "complex_nanoag_train_pre2022",
        "complex_nanoag_val_2022h1",
        "complex_nanoag_test_2022h2",
    ]:
        rows = _manifest_rows(paths, split_name, warnings, errors)
        split_rows[split_name] = rows
        _validate_manifest_rows(paths, split_name, rows, warnings, errors)
        _validate_temporal_bounds(split_name, rows, warnings, errors)

    single_train_ids = _target_ids(split_rows["single_nano_train_pre2022"])
    single_val_ids = _target_ids(split_rows["single_nano_val_2022h1"])
    single_test_ids = _target_ids(split_rows["single_nano_test_2022h2"])
    complex_train_ids = _target_ids(split_rows["complex_nanoag_train_pre2022"])
    complex_val_ids = _target_ids(split_rows["complex_nanoag_val_2022h1"])
    complex_test_ids = _target_ids(split_rows["complex_nanoag_test_2022h2"])
    train_ids = single_train_ids | complex_train_ids
    val_ids = single_val_ids | complex_val_ids
    test_ids = single_test_ids | complex_test_ids

    train_test_overlap = len(train_ids & test_ids)
    val_test_overlap = len(val_ids & test_ids)
    train_val_overlap = len(train_ids & val_ids)
    if train_test_overlap:
        errors.append(f"Train/test target overlap detected: {train_test_overlap}")
    if val_test_overlap:
        errors.append(f"Validation/test target overlap detected: {val_test_overlap}")
    if train_val_overlap:
        errors.append(f"Train/validation target overlap detected: {train_val_overlap}")

    if len(single_test_ids) != SINGLE_TEST_EXPECTED:
        errors.append(f"single_nano_test_2022h2 has {len(single_test_ids)} targets, expected {SINGLE_TEST_EXPECTED}")
    if len(complex_test_ids) != COMPLEX_TEST_EXPECTED:
        errors.append(f"complex_nanoag_test_2022h2 has {len(complex_test_ids)} targets, expected {COMPLEX_TEST_EXPECTED}")

    complex_train_status = "ok"
    complex_val_status = "ok"
    if len(complex_train_ids) != COMPLEX_TRAIN_PAPER_REFERENCE:
        complex_train_status = "warn"
        warnings.append(
            f"complex_nanoag_train_pre2022 has {len(complex_train_ids)} targets; "
            f"paper reference is {COMPLEX_TRAIN_PAPER_REFERENCE}. Current live SAbDab or local filters may differ."
        )
    if len(complex_val_ids) != COMPLEX_VAL_PAPER_REFERENCE:
        complex_val_status = "warn"
        warnings.append(
            f"complex_nanoag_val_2022h1 has {len(complex_val_ids)} targets; "
            f"paper reference is {COMPLEX_VAL_PAPER_REFERENCE}. Current live SAbDab or local filters may differ."
        )

    if (paths.manifests / "complex_nanoag_test_archive_filtered_out.csv").exists():
        filtered = read_csv_rows(paths.manifests / "complex_nanoag_test_archive_filtered_out.csv")
        if filtered:
            warnings.append(
                f"Official tFold NanoAg archive contained {len(complex_test_ids) + len(filtered)} candidates; "
                f"{len(filtered)} pre-filter candidates were excluded and documented in complex_nanoag_test_archive_filtered_out.csv."
            )

    fallback_dates = []
    if paths.rcsb_release_dates_csv.exists():
        for row in read_csv_rows(paths.rcsb_release_dates_csv):
            if row.get("source") == "sabdab_date_fallback":
                fallback_dates.append(row.get("pdb_id", ""))
    if fallback_dates:
        warnings.append(f"Used SAbDab date fallback for {len(fallback_dates)} PDB IDs after RCSB release-date lookup failed.")

    if not shutil_which("cd-hit"):
        warnings.append("CD-HIT unavailable; optional approximate 95% identity clustering was not run.")

    splits = {
        "single_nano_train_pre2022": {
            "actual_targets": len(single_train_ids),
            "status": "ok",
            "manifest_path": relpath(manifest_path(paths, "single_nano_train_pre2022"), paths.output_root),
        },
        "single_nano_val_2022h1": {
            "actual_targets": len(single_val_ids),
            "status": "ok",
            "manifest_path": relpath(manifest_path(paths, "single_nano_val_2022h1"), paths.output_root),
        },
        "single_nano_test_2022h2": {
            "expected_targets": SINGLE_TEST_EXPECTED,
            "actual_targets": len(single_test_ids),
            "status": "ok" if len(single_test_ids) == SINGLE_TEST_EXPECTED else "failed",
            "manifest_path": relpath(manifest_path(paths, "single_nano_test_2022h2"), paths.output_root),
        },
        "complex_nanoag_train_pre2022": {
            "paper_reference_targets": COMPLEX_TRAIN_PAPER_REFERENCE,
            "actual_targets": len(complex_train_ids),
            "status": complex_train_status,
            "manifest_path": relpath(manifest_path(paths, "complex_nanoag_train_pre2022"), paths.output_root),
        },
        "complex_nanoag_val_2022h1": {
            "paper_reference_targets": COMPLEX_VAL_PAPER_REFERENCE,
            "actual_targets": len(complex_val_ids),
            "status": complex_val_status,
            "manifest_path": relpath(manifest_path(paths, "complex_nanoag_val_2022h1"), paths.output_root),
        },
        "complex_nanoag_test_2022h2": {
            "expected_targets": COMPLEX_TEST_EXPECTED,
            "actual_targets": len(complex_test_ids),
            "status": "ok" if len(complex_test_ids) == COMPLEX_TEST_EXPECTED else "failed",
            "manifest_path": relpath(manifest_path(paths, "complex_nanoag_test_2022h2"), paths.output_root),
        },
    }

    leakage_checks = {
        "train_test_target_overlap": train_test_overlap,
        "val_test_target_overlap": val_test_overlap,
        "train_val_target_overlap": train_val_overlap,
        "exact_nanobody_sequence_leaks_removed": int(leakage_report.get("exact_nanobody_sequence_leaks_removed", 0) or 0),
        "exact_antigen_sequence_leaks_removed": int(leakage_report.get("exact_antigen_sequence_leaks_removed", 0) or 0),
    }

    report: Dict[str, object] = {
        "dataset_root": relpath(paths.dataset_root, paths.output_root),
        "download_date_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "source_urls": {
            "sabdab": SABDAB_MAIN_URL,
            "sabdab_nano": SABDAB_NANO_URL,
            "sabdab_archive_all": SABDAB_ARCHIVE_URL,
            "sabdab_summary_all": SABDAB_SUMMARY_URL,
            "sabdab_downloader": SABDAB_DOWNLOADER_URL,
            "tfold_repo": TFOLD_REPO_URL,
            "tfold_google_drive": TFOLD_GDRIVE_VIEW_URL,
            "tfold_google_drive_direct": TFOLD_GDRIVE_DIRECT_URL,
            "tfold_weiyun": TFOLD_WEIYUN_URL,
            "tfold_paper": TFOLD_PAPER_URL,
            "rcsb_entry_api": RCSB_ENTRY_API_TEMPLATE,
        },
        "checksums": {
            "all_structures_zip_sha256": _checksum(paths.sabdab_archive),
            "sabdab_summary_all_sha256": _checksum(paths.sabdab_summary),
            "tfold_test_set_tar_gz_sha256": _checksum(paths.tfold_archive),
        },
        "splits": splits,
        "leakage_checks": leakage_checks,
        "native_structures_present_in_tfold_archive": bool(
            leakage_report.get(
                "native_structures_present_in_tfold_archive",
                build_summary.get("native_structures_present_in_tfold_archive", False),
            )
        ),
        "warnings": warnings,
        "errors": errors,
    }
    report_path = paths.manifests / "verify_report.json"
    write_json(report_path, report)
    _write_dataset_card(paths, report)

    if errors:
        raise RuntimeError(f"Verification failed with {len(errors)} hard error(s); see {report_path}")
    return report


def shutil_which(name: str) -> Optional[str]:
    import shutil

    return shutil.which(name)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify SAbDab-nano temporal dataset bundle.")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("."),
        help="Directory under which data/sabdab_nano_temporal will be read/written.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        report = verify_dataset(args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    for split_name, payload in report["splits"].items():  # type: ignore[index,union-attr]
        print(f"{split_name}: {payload['actual_targets']} targets ({payload['status']})")
    print(f"Verification report: {dataset_paths(args.output_root).manifests / 'verify_report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
