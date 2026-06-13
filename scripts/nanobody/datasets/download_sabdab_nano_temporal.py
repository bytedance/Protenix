#!/usr/bin/env python3
"""Download and safely extract raw inputs for the SAbDab-nano temporal bundle."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from sabdab_nano_temporal_common import (
    SABDAB_ARCHIVE_NAME,
    SABDAB_ARCHIVE_URL,
    SABDAB_DOWNLOADER_URL,
    SABDAB_SUMMARY_NAME,
    SABDAB_SUMMARY_URL,
    TFOLD_ARCHIVE_NAME,
    TFOLD_GDRIVE_DIRECT_URL,
    TFOLD_GDRIVE_ID,
    TFOLD_GDRIVE_VIEW_URL,
    TFOLD_GIT_URL,
    TFOLD_WEIYUN_URL,
    command_exists,
    dataset_paths,
    download_with_wget_or_curl,
    ensure_base_layout,
    normalize_sabdab_structure_root,
    safe_extract_tar_gz,
    safe_extract_zip,
    sha256_file,
    validate_summary_tsv,
    validate_tar_gz,
    validate_zip,
    write_archive_file_trees,
    write_source_urls,
)


def _print_step(message: str) -> None:
    print(f"[sabdab-nano-temporal] {message}", flush=True)


def _has_extracted_sabdab_structures(root: Path) -> bool:
    return any((root / name).is_dir() and any((root / name).glob("*.pdb")) for name in ["raw", "chothia", "imgt"])


def _has_extracted_tfold(paths_root: Path) -> bool:
    return any(path.is_dir() and path.name in {"SAbDab-22H2-Nano", "SAbDab-22H2-NanoAg"} for path in paths_root.rglob("*"))


def _write_sabdab_sha256(paths) -> None:
    rows: List[str] = []
    for path in [paths.sabdab_archive, paths.sabdab_summary, paths.sabdab_downloader]:
        if path.exists() and path.stat().st_size > 0:
            rows.append(f"{sha256_file(path)}  {path.name}")
    paths.sabdab_sha256sums.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def _write_tfold_sha256(paths) -> None:
    rows: List[str] = []
    if paths.tfold_archive.exists() and paths.tfold_archive.stat().st_size > 0:
        rows.append(f"{sha256_file(paths.tfold_archive)}  {TFOLD_ARCHIVE_NAME}")
    paths.tfold_sha256sums.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def _download_required(url: str, output: Path, validator, force: bool, label: str) -> None:
    if not force and validator(output):
        _print_step(f"Reusing valid {label}: {output}")
        return
    if output.exists():
        output.unlink()
    _print_step(f"Downloading {label}: {url}")
    download_with_wget_or_curl(url, output)
    if not validator(output):
        raise RuntimeError(f"Downloaded {label} is invalid or appears to be HTML: {output}")


def _download_optional(url: str, output: Path, force: bool, label: str) -> None:
    if not force and output.exists() and output.stat().st_size > 0:
        _print_step(f"Reusing {label}: {output}")
        return
    try:
        _print_step(f"Downloading optional {label}: {url}")
        download_with_wget_or_curl(url, output, timeout=120)
    except Exception as exc:
        output.write_text(
            f"Could not download optional provenance helper from {url}\nReason: {exc}\n",
            encoding="utf-8",
        )
        _print_step(f"Optional {label} unavailable; wrote note to {output}")


def _record_tfold_head(paths) -> None:
    cmd = ["git", "ls-remote", TFOLD_GIT_URL, "HEAD"]
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=60,
        )
        paths.tfold_head.write_text(result.stdout, encoding="utf-8")
    except Exception as exc:
        paths.tfold_head.write_text(
            "Could not record tFold repository HEAD with:\n"
            f"{' '.join(cmd)}\n\nReason: {exc}\n",
            encoding="utf-8",
        )


def _run_download_command(cmd: Sequence[str], tmp_path: Path) -> bool:
    if tmp_path.exists():
        tmp_path.unlink()
    _print_step(f"Trying tFold download command: {' '.join(cmd)}")
    try:
        subprocess.run(list(cmd), check=True, timeout=None)
    except Exception as exc:
        _print_step(f"Command failed: {exc}")
        return False
    if validate_tar_gz(tmp_path):
        return True
    _print_step("Downloaded tFold file was not a valid gzip tar archive.")
    return False


def _copy_legacy_tfold_archive(paths, force: bool) -> bool:
    if force or validate_tar_gz(paths.tfold_archive):
        return False
    legacy = paths.output_root / "data" / "tfold_sabdab22h2" / "raw" / TFOLD_ARCHIVE_NAME
    if legacy.exists() and validate_tar_gz(legacy):
        _print_step(f"Reusing existing legacy tFold archive: {legacy}")
        shutil.copy2(legacy, paths.tfold_archive)
        return True
    return False


def _download_tfold_archive(paths, force: bool) -> None:
    if not force and validate_tar_gz(paths.tfold_archive):
        _print_step(f"Reusing valid tFold archive: {paths.tfold_archive}")
        return

    if _copy_legacy_tfold_archive(paths, force) and validate_tar_gz(paths.tfold_archive):
        return

    if paths.tfold_archive.exists():
        paths.tfold_archive.unlink()
    tmp = paths.tfold_archive.with_name(paths.tfold_archive.name + ".partial")
    commands: List[List[str]] = []

    gdown = shutil.which("gdown")
    if gdown:
        commands.extend(
            [
                [gdown, "--id", TFOLD_GDRIVE_ID, "-O", str(tmp)],
                [gdown, TFOLD_GDRIVE_ID, "-O", str(tmp)],
            ]
        )

    commands.extend(
        [
            [sys.executable, "-m", "gdown", "--id", TFOLD_GDRIVE_ID, "-O", str(tmp)],
            [sys.executable, "-m", "gdown", TFOLD_GDRIVE_ID, "-O", str(tmp)],
        ]
    )

    if command_exists("pip") or command_exists("pip3"):
        try:
            _print_step("Installing gdown with pip if it is not already available.")
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "-q", "--user", "gdown"],
                check=False,
                timeout=300,
            )
            commands.extend(
                [
                    [sys.executable, "-m", "gdown", "--id", TFOLD_GDRIVE_ID, "-O", str(tmp)],
                    [sys.executable, "-m", "gdown", TFOLD_GDRIVE_ID, "-O", str(tmp)],
                ]
            )
        except Exception:
            pass

    if shutil.which("curl"):
        commands.append(["curl", "-L", TFOLD_GDRIVE_DIRECT_URL, "-o", str(tmp)])

    for cmd in commands:
        if _run_download_command(cmd, tmp):
            tmp.replace(paths.tfold_archive)
            return

    instructions = (
        "Failed to download a valid tFold test-set archive from Google Drive.\n"
        "Manual fallback:\n"
        f"1. Open the official Google Drive URL: {TFOLD_GDRIVE_VIEW_URL}\n"
        f"2. If Google Drive is unavailable, open the Tencent Weiyun mirror: {TFOLD_WEIYUN_URL}\n"
        f"3. Download {TFOLD_ARCHIVE_NAME}\n"
        f"4. Place it at: {paths.tfold_archive}\n"
        "5. Re-run bash scripts/nanobody/datasets/download_sabdab_nano_temporal.sh\n"
    )
    raise RuntimeError(instructions)


def _extract_archives(paths, force: bool) -> None:
    if force or not _has_extracted_sabdab_structures(paths.sabdab_structures):
        _print_step(f"Extracting SAbDab archive to {paths.sabdab_structures}")
        safe_extract_zip(paths.sabdab_archive, paths.sabdab_structures)
        normalize_sabdab_structure_root(paths.sabdab_structures)
    else:
        _print_step(f"Reusing extracted SAbDab structures: {paths.sabdab_structures}")

    if force or not _has_extracted_tfold(paths.tfold_extracted):
        _print_step(f"Extracting tFold archive to {paths.tfold_extracted}")
        safe_extract_tar_gz(paths.tfold_archive, paths.tfold_extracted)
    else:
        _print_step(f"Reusing extracted tFold archive: {paths.tfold_extracted}")


def run_download(output_root: Path, force: bool = False) -> None:
    paths = dataset_paths(output_root)
    ensure_base_layout(paths)
    write_source_urls(paths)

    _download_required(
        SABDAB_ARCHIVE_URL,
        paths.sabdab_archive,
        validate_zip,
        force,
        "SAbDab all-structures archive",
    )
    _download_required(
        SABDAB_SUMMARY_URL,
        paths.sabdab_summary,
        validate_summary_tsv,
        force,
        "SAbDab all-summary TSV",
    )
    _download_optional(
        SABDAB_DOWNLOADER_URL,
        paths.sabdab_downloader,
        force,
        "SAbDab downloader provenance script",
    )
    _download_tfold_archive(paths, force)
    _record_tfold_head(paths)

    if not validate_zip(paths.sabdab_archive):
        raise RuntimeError(f"SAbDab archive is invalid: {paths.sabdab_archive}")
    if not validate_summary_tsv(paths.sabdab_summary):
        raise RuntimeError(f"SAbDab summary is invalid or missing required columns: {paths.sabdab_summary}")
    if not validate_tar_gz(paths.tfold_archive):
        raise RuntimeError(f"tFold archive is invalid: {paths.tfold_archive}")

    _write_sabdab_sha256(paths)
    _write_tfold_sha256(paths)
    _extract_archives(paths, force)
    write_archive_file_trees(paths)

    _print_step(f"Raw SAbDab data: {paths.raw_sabdab}")
    _print_step(f"Raw tFold data: {paths.raw_tfold}")
    _print_step(f"Archive file tree: {paths.archive_file_trees}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download raw inputs for the SAbDab-nano temporal dataset bundle."
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("."),
        help="Directory under which data/sabdab_nano_temporal will be created.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download and re-extract raw archives even if valid files exist.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        run_download(args.output_root, args.force)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
