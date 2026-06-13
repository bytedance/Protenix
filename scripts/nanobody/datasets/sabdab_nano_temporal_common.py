#!/usr/bin/env python3
"""Shared helpers for the SAbDab-nano temporal dataset pipeline."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence


DATASET_REL = Path("data") / "sabdab_nano_temporal"

SABDAB_ARCHIVE_URL = "https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/archive/all/"
SABDAB_SUMMARY_URL = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/"
)
SABDAB_DOWNLOADER_URL = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/"
    "sabdab_downloader.py"
)
SABDAB_MAIN_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab"
SABDAB_NANO_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/nano/"

TFOLD_REPO_URL = "https://github.com/TencentAI4S/tfold"
TFOLD_GIT_URL = "https://github.com/TencentAI4S/tfold.git"
TFOLD_GDRIVE_ID = "1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ"
TFOLD_GDRIVE_VIEW_URL = (
    "https://drive.google.com/file/d/1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ/view"
    "?usp=drive_link"
)
TFOLD_GDRIVE_DIRECT_URL = (
    "https://drive.google.com/uc?export=download&id="
    "1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ"
)
TFOLD_WEIYUN_URL = "https://share.weiyun.com/zycZDrfA"
TFOLD_PAPER_URL = "https://www.nature.com/articles/s41467-025-67361-9"
RCSB_ENTRY_API_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_CIF_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"

SABDAB_ARCHIVE_NAME = "all_structures.zip"
SABDAB_SUMMARY_NAME = "sabdab_summary_all.tsv"
TFOLD_ARCHIVE_NAME = "tFold_test_set.tar.gz"

SINGLE_TEST_EXPECTED = 73
COMPLEX_TEST_EXPECTED = 41
COMPLEX_TRAIN_PAPER_REFERENCE = 1319
COMPLEX_VAL_PAPER_REFERENCE = 40

SPLIT_NAMES = [
    "single_nano_train_pre2022",
    "single_nano_val_2022h1",
    "single_nano_test_2022h2",
    "complex_nanoag_train_pre2022",
    "complex_nanoag_val_2022h1",
    "complex_nanoag_test_2022h2",
]

MANIFEST_FIELDS = [
    "task",
    "split",
    "subset",
    "target_id",
    "pdb_id",
    "release_date",
    "release_date_source",
    "nanobody_chain",
    "heavy_chain",
    "light_chain",
    "antigen_chains",
    "antigen_type",
    "antigen_name",
    "method",
    "resolution",
    "r_factor",
    "r_free",
    "raw_pdb_path",
    "chothia_pdb_path",
    "imgt_pdb_path",
    "native_cif_path",
    "fasta_path",
    "json_path",
    "source_summary_path",
    "source_archive",
    "source_url",
    "was_bound_for_single_task",
    "passed_quality_filters",
    "filter_notes",
    "notes",
]


@dataclass(frozen=True)
class DatasetPaths:
    output_root: Path
    dataset_root: Path
    raw: Path
    raw_sabdab: Path
    sabdab_archive: Path
    sabdab_structures: Path
    sabdab_summary: Path
    sabdab_downloader: Path
    sabdab_sha256sums: Path
    sabdab_source_urls: Path
    raw_tfold: Path
    tfold_archive: Path
    tfold_extracted: Path
    tfold_head: Path
    tfold_sha256sums: Path
    tfold_source_urls: Path
    raw_rcsb: Path
    rcsb_release_dates_csv: Path
    rcsb_release_dates_cache: Path
    rcsb_downloaded_cif: Path
    rcsb_downloaded_pdb: Path
    splits: Path
    manifests: Path
    archive_file_trees: Path


def dataset_paths(output_root: Path) -> DatasetPaths:
    root = output_root.resolve()
    dataset = root / DATASET_REL
    raw = dataset / "raw"
    raw_sabdab = raw / "sabdab"
    raw_tfold = raw / "tfold"
    raw_rcsb = raw / "rcsb"
    return DatasetPaths(
        output_root=root,
        dataset_root=dataset,
        raw=raw,
        raw_sabdab=raw_sabdab,
        sabdab_archive=raw_sabdab / SABDAB_ARCHIVE_NAME,
        sabdab_structures=raw_sabdab / "all_structures",
        sabdab_summary=raw_sabdab / SABDAB_SUMMARY_NAME,
        sabdab_downloader=raw_sabdab / "sabdab_downloader.py",
        sabdab_sha256sums=raw_sabdab / "SHA256SUMS.txt",
        sabdab_source_urls=raw_sabdab / "source_urls.txt",
        raw_tfold=raw_tfold,
        tfold_archive=raw_tfold / TFOLD_ARCHIVE_NAME,
        tfold_extracted=raw_tfold / "extracted",
        tfold_head=raw_tfold / "tfold_repo_HEAD.txt",
        tfold_sha256sums=raw_tfold / "SHA256SUMS.txt",
        tfold_source_urls=raw_tfold / "source_urls.txt",
        raw_rcsb=raw_rcsb,
        rcsb_release_dates_csv=raw_rcsb / "pdb_release_dates.csv",
        rcsb_release_dates_cache=raw_rcsb / "pdb_release_dates_cache.json",
        rcsb_downloaded_cif=raw_rcsb / "downloaded_cif",
        rcsb_downloaded_pdb=raw_rcsb / "downloaded_pdb",
        splits=dataset / "splits",
        manifests=dataset / "manifests",
        archive_file_trees=dataset / "manifests" / "archive_file_trees.txt",
    )


def ensure_base_layout(paths: DatasetPaths) -> None:
    for path in [
        paths.raw_sabdab,
        paths.sabdab_structures,
        paths.raw_tfold,
        paths.tfold_extracted,
        paths.raw_rcsb,
        paths.rcsb_downloaded_cif,
        paths.rcsb_downloaded_pdb,
        paths.splits / "single_nano" / "train_pre2022" / "structures",
        paths.splits / "single_nano" / "train_pre2022" / "metadata",
        paths.splits / "single_nano" / "val_2022h1" / "structures",
        paths.splits / "single_nano" / "val_2022h1" / "metadata",
        paths.splits / "single_nano" / "test_2022h2" / "structures",
        paths.splits / "single_nano" / "test_2022h2" / "metadata",
        paths.splits / "complex_nanoag" / "train_pre2022" / "structures",
        paths.splits / "complex_nanoag" / "train_pre2022" / "metadata",
        paths.splits / "complex_nanoag" / "val_2022h1" / "structures",
        paths.splits / "complex_nanoag" / "val_2022h1" / "metadata",
        paths.splits / "complex_nanoag" / "test_2022h2" / "structures",
        paths.splits / "complex_nanoag" / "test_2022h2" / "metadata",
        paths.manifests,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def relpath(path: Path, output_root: Path) -> str:
    if not path:
        return ""
    absolute_path = Path(os.path.abspath(path))
    absolute_root = output_root.resolve()
    try:
        return absolute_path.relative_to(absolute_root).as_posix()
    except ValueError:
        return absolute_path.as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_json(path: Path, default: object) -> object:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def is_missing(value: object) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"na", "n/a", "none", "null", "-", "nan"}


def clean_value(value: object) -> str:
    return "" if is_missing(value) else str(value).strip()


def normalize_pdb_id(value: object) -> str:
    return clean_value(value).lower()


def split_chains(value: object) -> List[str]:
    text = clean_value(value)
    if not text:
        return []
    chains = []
    for item in re.split(r"[,;|/\s]+", text):
        item = item.strip()
        if item and item.upper() not in {"NA", "NONE", "NULL"}:
            chains.append(item)
    return sorted(dict.fromkeys(chains))


def joined_chains(chains: Sequence[str]) -> str:
    return ",".join(chains)


def target_chains_for_id(chains: Sequence[str]) -> str:
    return "_".join(chains)


def write_source_urls(paths: DatasetPaths) -> None:
    paths.sabdab_source_urls.write_text(
        "\n".join(
            [
                SABDAB_MAIN_URL,
                SABDAB_NANO_URL,
                SABDAB_ARCHIVE_URL,
                SABDAB_SUMMARY_URL,
                SABDAB_DOWNLOADER_URL,
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    paths.tfold_source_urls.write_text(
        "\n".join(
            [
                TFOLD_REPO_URL,
                TFOLD_GDRIVE_VIEW_URL,
                TFOLD_GDRIVE_DIRECT_URL,
                TFOLD_WEIYUN_URL,
                TFOLD_PAPER_URL,
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def command_exists(name: str) -> bool:
    return shutil.which(name) is not None


def run_command(cmd: Sequence[str], timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        list(cmd),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout,
    )


def download_with_wget_or_curl(url: str, output: Path, timeout: Optional[int] = None) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + ".partial")
    if tmp.exists():
        tmp.unlink()
    errors: List[str] = []
    if command_exists("wget"):
        cmd = ["wget", "-q", "-O", str(tmp), url]
        try:
            subprocess.run(cmd, check=True, timeout=timeout)
            tmp.replace(output)
            return
        except Exception as exc:
            errors.append(f"{' '.join(cmd)}: {exc}")
            if tmp.exists():
                tmp.unlink()
    if command_exists("curl"):
        cmd = ["curl", "-L", "--fail", "--silent", "--show-error", url, "-o", str(tmp)]
        try:
            subprocess.run(cmd, check=True, timeout=timeout)
            tmp.replace(output)
            return
        except Exception as exc:
            errors.append(f"{' '.join(cmd)}: {exc}")
            if tmp.exists():
                tmp.unlink()
    raise RuntimeError(f"Could not download {url}: {'; '.join(errors)}")


def looks_like_html(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return True
    with path.open("rb") as handle:
        sample = handle.read(4096).lstrip().lower()
    return sample.startswith(b"<!doctype html") or sample.startswith(b"<html")


def validate_zip(path: Path) -> bool:
    return path.exists() and not looks_like_html(path) and zipfile.is_zipfile(path)


def validate_tar_gz(path: Path) -> bool:
    if not path.exists() or looks_like_html(path):
        return False
    with path.open("rb") as handle:
        if handle.read(2) != b"\x1f\x8b":
            return False
    try:
        with tarfile.open(path, "r:gz") as tar:
            tar.next()
        return True
    except tarfile.TarError:
        return False


def validate_summary_tsv(path: Path) -> bool:
    if not path.exists() or looks_like_html(path):
        return False
    required = {"pdb", "hchain", "lchain", "antigen_chain", "antigen_type", "date"}
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
    except Exception:
        return False
    present = {name.strip().lower() for name in header}
    return required.issubset(present)


def _safe_member_path(name: str) -> Path:
    member = Path(name)
    if member.is_absolute() or ".." in member.parts:
        raise RuntimeError(f"Unsafe archive member path: {name}")
    return member


def safe_extract_zip(zip_path: Path, dest: Path) -> None:
    if not validate_zip(zip_path):
        raise RuntimeError(f"Not a valid zip archive: {zip_path}")
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        for info in archive.infolist():
            _safe_member_path(info.filename)
        archive.extractall(dest)


def safe_extract_tar_gz(tar_path: Path, dest: Path) -> None:
    if not validate_tar_gz(tar_path):
        raise RuntimeError(f"Not a valid gzip tar archive: {tar_path}")
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as archive:
        members = archive.getmembers()
        for member in members:
            _safe_member_path(member.name)
            if member.issym() or member.islnk():
                _safe_member_path(member.linkname)
        archive.extractall(dest, members=members)


def normalize_sabdab_structure_root(structures_root: Path) -> None:
    """Move a nested all_structures folder up so raw/chothia/imgt are direct children."""
    nested = structures_root / "all_structures"
    if nested.is_dir() and any((nested / name).exists() for name in ["raw", "chothia", "imgt"]):
        for child in list(nested.iterdir()):
            target = structures_root / child.name
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink()
            shutil.move(str(child), str(target))
        nested.rmdir()


def write_archive_file_trees(paths: DatasetPaths) -> None:
    entries: List[str] = []
    for label, root in [
        ("sabdab", paths.sabdab_structures),
        ("tfold", paths.tfold_extracted),
    ]:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.is_file() or path.is_symlink():
                entries.append(f"{label}\t{relpath(path, paths.output_root)}")
    paths.archive_file_trees.parent.mkdir(parents=True, exist_ok=True)
    paths.archive_file_trees.write_text("\n".join(sorted(entries)) + "\n", encoding="utf-8")


def link_or_copy(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        relative = os.path.relpath(src.resolve(), start=dst.parent.resolve())
        dst.symlink_to(relative)
    except OSError:
        shutil.copy2(src, dst)
    return True


def maybe_rel(path: Optional[Path], output_root: Path) -> str:
    if path is None or not path.exists():
        return ""
    return relpath(path, output_root)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def manifest_path(paths: DatasetPaths, split_name: str) -> Path:
    return paths.manifests / f"{split_name}_manifest.csv"
