#!/usr/bin/env python3
"""Verify and manifest the tFold SAbDab-22H2 nanobody benchmarks."""

from __future__ import annotations

import argparse
import csv
import datetime as _datetime
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from download_tfold_sabdab22h2 import (
    MANIFEST_FIELDS,
    SOURCE_URLS,
    SUBSET_EXPECTED_COUNTS,
    DatasetPaths,
    FilterResult,
    dataset_paths,
    final_target_filter,
    locate_subset_dir,
    parse_target_id,
    read_prot_ids,
    relpath,
    sha256_file,
    validate_archive,
)


def _path_or_blank(path: Path, output_root: Path) -> str:
    return relpath(path, output_root) if path.exists() else ""


def build_manifest_rows(paths: DatasetPaths, subset_name: str) -> List[Dict[str, str]]:
    subset_dir = paths.subsets / subset_name
    target_ids = read_prot_ids(subset_dir)
    rows: List[Dict[str, str]] = []
    for target_id in target_ids:
        info = parse_target_id(target_id)
        fasta = subset_dir / "fasta.files" / f"{target_id}.fasta"
        native_pdb = subset_dir / "pdb.files.native" / f"{target_id}.pdb"
        native_cif = subset_dir / "native_cif" / f"{target_id}.cif"
        json_path = subset_dir / "json.files" / f"{target_id}.json"

        msa_paths = []
        for chain in info.antigen_chains:
            msa = subset_dir / "msa.files" / f"{info.pdb_id}_{chain}.a3m"
            if msa.exists():
                msa_paths.append(relpath(msa, paths.output_root))

        native_note = "native PDB included in official tFold archive" if native_pdb.exists() else ""
        rows.append(
            {
                "subset": subset_name,
                "target_id": target_id,
                "pdb_id": info.pdb_id,
                "nanobody_chain": ",".join(info.nanobody_chains),
                "heavy_chain": ",".join(info.nanobody_chains),
                "light_chain": "",
                "antigen_chains": ",".join(info.antigen_chains),
                "fasta_path": _path_or_blank(fasta, paths.output_root),
                "json_path": _path_or_blank(json_path, paths.output_root),
                "native_pdb_path": _path_or_blank(native_pdb, paths.output_root),
                "native_cif_path": _path_or_blank(native_cif, paths.output_root),
                "msa_path": ";".join(msa_paths),
                "notes": native_note,
            }
        )
    return rows


def write_manifest(paths: DatasetPaths, subset_name: str, rows: List[Dict[str, str]]) -> Path:
    path = paths.manifests / f"{subset_name}_manifest.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def summarize_filter_results(results: Iterable[FilterResult]) -> List[Dict[str, object]]:
    output: List[Dict[str, object]] = []
    for result in results:
        if result.include:
            continue
        output.append(
            {
                "target_id": result.target_id,
                "reasons": list(result.reasons),
                "nanobody_sequence_residues": result.nanobody_sequence_residues,
                "nanobody_structure_ca_residues": result.nanobody_structure_ca_residues,
                "antigen_sequence_residues": result.antigen_sequence_residues,
                "interfacial_ca_contacts_lt_10a": result.interfacial_ca_contacts_lt_10a,
            }
        )
    return output


def write_excluded_targets(paths: DatasetPaths, subset_name: str, results: List[FilterResult]) -> Path:
    path = paths.manifests / f"{subset_name}_excluded_targets.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "subset",
            "target_id",
            "reasons",
            "nanobody_sequence_residues",
            "nanobody_structure_ca_residues",
            "antigen_sequence_residues",
            "interfacial_ca_contacts_lt_10a",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            if result.include:
                continue
            writer.writerow(
                {
                    "subset": subset_name,
                    "target_id": result.target_id,
                    "reasons": ";".join(result.reasons),
                    "nanobody_sequence_residues": result.nanobody_sequence_residues,
                    "nanobody_structure_ca_residues": result.nanobody_structure_ca_residues,
                    "antigen_sequence_residues": result.antigen_sequence_residues,
                    "interfacial_ca_contacts_lt_10a": result.interfacial_ca_contacts_lt_10a,
                }
            )
    return path


def write_dataset_card(
    paths: DatasetPaths,
    source_sha256: str,
    report: Dict[str, object],
    manifest_paths: Dict[str, Path],
) -> None:
    today = _datetime.date.today().isoformat()
    lines = [
        "# tFold SAbDab-22H2 nanobody benchmark subsets",
        "",
        f"Download date: {today}",
        f"Source archive SHA256: `{source_sha256}`",
        "",
        "## Datasets",
    ]
    for subset_name, expected in SUBSET_EXPECTED_COUNTS.items():
        subset_report = report["subsets"][subset_name]  # type: ignore[index]
        lines.extend(
            [
                f"- `{subset_name}`: expected {expected}, actual {subset_report['actual_targets']} target rows.",
                f"  Local subset path: `{relpath(paths.subsets / subset_name, paths.output_root)}`",
                f"  Manifest: `{relpath(manifest_paths[subset_name], paths.output_root)}`",
            ]
        )

    lines.extend(
        [
            "",
            "## Source URLs",
            *[f"- {url}" for url in SOURCE_URLS],
            "",
            "## Local paths",
            f"- Archive: `{relpath(paths.archive, paths.output_root)}`",
            f"- Raw provenance: `{relpath(paths.raw, paths.output_root)}`",
            f"- Extracted official archive: `{relpath(paths.extracted, paths.output_root)}`",
            f"- Canonical subsets: `{relpath(paths.subsets, paths.output_root)}`",
            f"- Verification report: `{relpath(paths.manifests / 'verify_report.json', paths.output_root)}`",
            "",
            "## Caveats and discrepancies",
        ]
    )

    notes = report.get("notes", [])
    if notes:
        lines.extend([f"- {note}" for note in notes])  # type: ignore[union-attr]
    else:
        lines.append("- No discrepancies were detected.")

    lines.extend(
        [
            "",
            "Native structures are included as PDB files in the official tFold archive for all retained targets. No RCSB native files were fetched by this workflow.",
            "",
            "## Provenance and citation",
            "These are tFold benchmark subsets constructed from SAbDab/PDB structures. Cite the tFold paper and SAbDab when using them.",
            "",
            "## How to verify",
            "",
            "```bash",
            "python scripts/nanobody/datasets/verify_tfold_sabdab22h2.py",
            "```",
            "",
            "## One-command download",
            "",
            "```bash",
            "bash scripts/nanobody/datasets/download_tfold_sabdab22h2.sh",
            "```",
            "",
        ]
    )
    (paths.manifests / "DATASET_CARD.md").write_text("\n".join(lines), encoding="utf-8")


def verify_dataset(output_root: Path) -> Dict[str, object]:
    paths = dataset_paths(output_root)
    paths.manifests.mkdir(parents=True, exist_ok=True)
    notes: List[str] = []
    status_failed = False

    if not validate_archive(paths.archive):
        raise RuntimeError(f"Missing or invalid source archive: {paths.archive}")
    source_sha256 = sha256_file(paths.archive)

    report: Dict[str, object] = {
        "source_archive": relpath(paths.archive, paths.output_root),
        "source_archive_sha256": source_sha256,
        "source_urls": SOURCE_URLS,
        "subsets": {},
        "notes": notes,
    }
    manifest_paths: Dict[str, Path] = {}

    for subset_name, expected in SUBSET_EXPECTED_COUNTS.items():
        source_subset = locate_subset_dir(paths.extracted, subset_name)
        archive_ids = read_prot_ids(source_subset)
        expected_ids, filter_results = final_target_filter(subset_name, source_subset)
        excluded = summarize_filter_results(filter_results)
        if excluded:
            write_excluded_targets(paths, subset_name, filter_results)
            notes.append(
                f"{subset_name}: official archive contains {len(archive_ids)} targets; "
                f"published final-paper filters exclude {len(excluded)} targets "
                f"({', '.join(item['target_id'] for item in excluded)}), yielding {len(expected_ids)}."
            )

        rows = build_manifest_rows(paths, subset_name)
        manifest_path = write_manifest(paths, subset_name, rows)
        manifest_paths[subset_name] = manifest_path

        actual_ids = [row["target_id"] for row in rows]
        unique_actual_ids = sorted(set(actual_ids))
        missing = sorted(set(expected_ids) - set(unique_actual_ids))
        extra = sorted(set(unique_actual_ids) - set(expected_ids))
        duplicate_count = len(actual_ids) - len(unique_actual_ids)
        subset_status = "ok"
        subset_notes: List[str] = []
        if len(unique_actual_ids) != expected:
            subset_status = "failed"
            subset_notes.append(
                f"Expected {expected} unique targets, found {len(unique_actual_ids)}."
            )
        if missing:
            subset_status = "failed"
            subset_notes.append(f"Missing expected targets: {', '.join(missing)}")
        if extra:
            subset_status = "failed"
            subset_notes.append(f"Unexpected targets: {', '.join(extra)}")
        if duplicate_count:
            subset_status = "failed"
            subset_notes.append(f"Duplicate target_id rows: {duplicate_count}")
        if subset_status != "ok":
            status_failed = True
            notes.extend(f"{subset_name}: {note}" for note in subset_notes)

        report["subsets"][subset_name] = {  # type: ignore[index]
            "expected_targets": expected,
            "actual_targets": len(unique_actual_ids),
            "archive_targets": len(archive_ids),
            "status": subset_status,
            "subset_path": relpath(paths.subsets / subset_name, paths.output_root),
            "manifest_path": relpath(manifest_path, paths.output_root),
            "excluded_targets": excluded,
            "notes": subset_notes,
        }

    report_path = paths.manifests / "verify_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_dataset_card(paths, source_sha256, report, manifest_paths)

    if status_failed:
        raise RuntimeError(f"Verification failed; see {report_path}")
    return report


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify tFold SAbDab-22H2 nanobody benchmark subsets."
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("."),
        help="Directory under which data/tfold_sabdab22h2 will be read. Defaults to '.'.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        report = verify_dataset(args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    for subset_name, subset_report in report["subsets"].items():  # type: ignore[union-attr]
        print(
            f"{subset_name}: {subset_report['actual_targets']} targets "
            f"({subset_report['status']})"
        )
    print(f"Verification report: {dataset_paths(args.output_root).manifests / 'verify_report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
