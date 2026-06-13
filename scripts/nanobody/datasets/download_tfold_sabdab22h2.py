#!/usr/bin/env python3
"""Download and canonicalize the tFold SAbDab-22H2 nanobody benchmarks."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tarfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DATASET_REL = Path("data") / "tfold_sabdab22h2"
ARCHIVE_NAME = "tFold_test_set.tar.gz"
GOOGLE_DRIVE_ID = "1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ"
GOOGLE_DRIVE_DIRECT_URL = (
    "https://drive.google.com/uc?export=download&id="
    "1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ"
)
TFOLD_GIT_URL = "https://github.com/TencentAI4S/tfold.git"

SUBSET_EXPECTED_COUNTS = {
    "SAbDab-22H2-Nano": 73,
    "SAbDab-22H2-NanoAg": 41,
}

SOURCE_URLS = [
    "https://github.com/TencentAI4S/tfold",
    "https://drive.google.com/file/d/1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ/view?usp=drive_link",
    GOOGLE_DRIVE_DIRECT_URL,
    "https://share.weiyun.com/zycZDrfA",
    "https://www.nature.com/articles/s41467-025-67361-9",
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab",
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/nano/",
]

MANIFEST_FIELDS = [
    "subset",
    "target_id",
    "pdb_id",
    "nanobody_chain",
    "heavy_chain",
    "light_chain",
    "antigen_chains",
    "fasta_path",
    "json_path",
    "native_pdb_path",
    "native_cif_path",
    "msa_path",
    "notes",
]


@dataclass(frozen=True)
class DatasetPaths:
    output_root: Path
    dataset: Path
    raw: Path
    archive: Path
    sha256sums: Path
    source_urls: Path
    tfold_head: Path
    extracted: Path
    subsets: Path
    manifests: Path
    archive_file_tree: Path


@dataclass(frozen=True)
class TargetInfo:
    target_id: str
    pdb_id: str
    nanobody_chains: Tuple[str, ...]
    antigen_chains: Tuple[str, ...]


@dataclass(frozen=True)
class FilterResult:
    target_id: str
    include: bool
    reasons: Tuple[str, ...]
    nanobody_sequence_residues: Optional[int] = None
    nanobody_structure_ca_residues: Optional[int] = None
    antigen_sequence_residues: Optional[int] = None
    interfacial_ca_contacts_lt_10a: Optional[int] = None


def dataset_paths(output_root: Path) -> DatasetPaths:
    root = output_root.resolve()
    dataset = root / DATASET_REL
    raw = dataset / "raw"
    manifests = dataset / "manifests"
    return DatasetPaths(
        output_root=root,
        dataset=dataset,
        raw=raw,
        archive=raw / ARCHIVE_NAME,
        sha256sums=raw / "SHA256SUMS.txt",
        source_urls=raw / "source_urls.txt",
        tfold_head=raw / "tfold_repo_HEAD.txt",
        extracted=dataset / "extracted",
        subsets=dataset / "subsets",
        manifests=manifests,
        archive_file_tree=manifests / "archive_file_tree.txt",
    )


def relpath(path: Path, output_root: Path) -> str:
    try:
        return path.resolve().relative_to(output_root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def ensure_layout(paths: DatasetPaths) -> None:
    for path in [paths.raw, paths.extracted, paths.subsets, paths.manifests]:
        path.mkdir(parents=True, exist_ok=True)


def write_source_urls(paths: DatasetPaths) -> None:
    paths.source_urls.write_text("\n".join(SOURCE_URLS) + "\n", encoding="utf-8")


def record_tfold_head(paths: DatasetPaths) -> None:
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
            f"{' '.join(cmd)}\n\n"
            f"Reason: {exc}\n",
            encoding="utf-8",
        )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_archive(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
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


def _run_download_command(cmd: Sequence[str], tmp_path: Path) -> bool:
    if tmp_path.exists():
        tmp_path.unlink()
    print(f"Trying: {' '.join(cmd)}", flush=True)
    try:
        subprocess.run(cmd, check=True, timeout=None)
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"Command failed: {exc}", flush=True)
        return False
    if validate_archive(tmp_path):
        return True
    print("Downloaded file was not a valid gzip/tar archive.", flush=True)
    return False


def download_archive(paths: DatasetPaths) -> None:
    if validate_archive(paths.archive):
        print(f"Reusing valid archive: {paths.archive}", flush=True)
        write_sha256(paths)
        return

    paths.raw.mkdir(parents=True, exist_ok=True)
    tmp_path = paths.archive.with_suffix(paths.archive.suffix + ".partial")
    commands: List[List[str]] = []

    gdown = shutil.which("gdown")
    if gdown:
        commands.extend(
            [
                [gdown, "--id", GOOGLE_DRIVE_ID, "-O", str(tmp_path)],
                [gdown, GOOGLE_DRIVE_ID, "-O", str(tmp_path)],
            ]
        )

    commands.extend(
        [
            [sys.executable, "-m", "gdown", "--id", GOOGLE_DRIVE_ID, "-O", str(tmp_path)],
            [sys.executable, "-m", "gdown", GOOGLE_DRIVE_ID, "-O", str(tmp_path)],
        ]
    )

    pip = [sys.executable, "-m", "pip"]
    try:
        subprocess.run(
            pip + ["--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
        )
        print("Installing gdown with pip if needed.", flush=True)
        subprocess.run(pip + ["install", "-q", "--user", "gdown"], check=False)
        commands.extend(
            [
                [sys.executable, "-m", "gdown", "--id", GOOGLE_DRIVE_ID, "-O", str(tmp_path)],
                [sys.executable, "-m", "gdown", GOOGLE_DRIVE_ID, "-O", str(tmp_path)],
            ]
        )
    except Exception:
        pass

    uvx = shutil.which("uvx")
    if uvx:
        commands.append([uvx, "gdown", GOOGLE_DRIVE_ID, "-O", str(tmp_path)])

    curl = shutil.which("curl")
    if curl:
        commands.append([curl, "-L", GOOGLE_DRIVE_DIRECT_URL, "-o", str(tmp_path)])

    for cmd in commands:
        if _run_download_command(cmd, tmp_path):
            tmp_path.replace(paths.archive)
            write_sha256(paths)
            return

    instructions = (
        "Failed to download a valid tFold test-set archive from Google Drive.\n"
        "Manual fallback:\n"
        "1. Open the official Tencent Weiyun mirror: https://share.weiyun.com/zycZDrfA\n"
        f"2. Download {ARCHIVE_NAME}\n"
        f"3. Place it at: {paths.archive}\n"
        "4. Re-run this script.\n"
    )
    raise RuntimeError(instructions)


def write_sha256(paths: DatasetPaths) -> str:
    digest = sha256_file(paths.archive)
    paths.sha256sums.write_text(f"{digest}  {ARCHIVE_NAME}\n", encoding="utf-8")
    return digest


def safe_extract_archive(paths: DatasetPaths) -> None:
    if not validate_archive(paths.archive):
        raise RuntimeError(f"Not a valid gzip/tar archive: {paths.archive}")

    if paths.extracted.exists():
        shutil.rmtree(paths.extracted)
    paths.extracted.mkdir(parents=True, exist_ok=True)

    with tarfile.open(paths.archive, "r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            name = Path(member.name)
            if name.is_absolute() or ".." in name.parts:
                raise RuntimeError(f"Unsafe archive member path: {member.name}")
            if member.issym() or member.islnk():
                link = Path(member.linkname)
                if link.is_absolute() or ".." in link.parts:
                    raise RuntimeError(
                        f"Unsafe archive link target: {member.name} -> {member.linkname}"
                    )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The default behavior of tarfile extraction has been changed.*",
                category=RuntimeWarning,
            )
            tar.extractall(paths.extracted, members=members)


def write_archive_file_tree(paths: DatasetPaths) -> None:
    entries: List[str] = []
    for path in paths.extracted.rglob("*"):
        if not path.is_file():
            continue
        depth = len(path.relative_to(paths.extracted).parts)
        if depth <= 8:
            entries.append(relpath(path, paths.output_root))
    paths.archive_file_tree.write_text("\n".join(sorted(entries)) + "\n", encoding="utf-8")


def read_prot_ids(subset_dir: Path) -> List[str]:
    prot_ids = subset_dir / "prot_ids.txt"
    if prot_ids.exists():
        return [line.strip() for line in prot_ids.read_text(encoding="utf-8").splitlines() if line.strip()]

    ids = set()
    for suffix in [".fasta", ".pdb", ".cif", ".json"]:
        for path in subset_dir.rglob(f"*{suffix}"):
            ids.add(path.stem)
    return sorted(ids)


def parse_target_id(target_id: str) -> TargetInfo:
    parts = target_id.split("_")
    pdb_id = parts[0].lower() if parts else ""
    if "NA" in parts:
        idx = parts.index("NA")
        nanobody_chains = tuple(parts[1:idx])
        antigen_chains = tuple(parts[idx + 1 :])
    else:
        nanobody_chains = tuple(parts[1:2])
        antigen_chains = ()
    return TargetInfo(target_id, pdb_id, nanobody_chains, antigen_chains)


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def parse_pdb_ca(path: Path) -> Dict[str, Dict[Tuple[str, str], Tuple[float, float, float]]]:
    residues: Dict[str, Dict[Tuple[str, str], Tuple[float, float, float]]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain = line[21].strip()
        residue_id = (line[22:26].strip(), line[26].strip())
        try:
            xyz = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue
        residues.setdefault(chain, {})[residue_id] = xyz
    return residues


def _distance_lt_10a(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> bool:
    return (
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    ) < 100.0


def evaluate_nanoag_target(subset_dir: Path, target_id: str) -> FilterResult:
    fasta = subset_dir / "fasta.files" / f"{target_id}.fasta"
    native_pdb = subset_dir / "pdb.files.native" / f"{target_id}.pdb"
    reasons: List[str] = []
    if not fasta.exists():
        return FilterResult(target_id, False, ("missing_fasta",))
    if not native_pdb.exists():
        return FilterResult(target_id, False, ("missing_native_pdb",))

    records = read_fasta(fasta)
    if len(records) < 2:
        return FilterResult(target_id, False, ("fasta_missing_antigen",))

    nanobody_sequence_residues = len(records[0][1])
    antigen_sequence_residues = sum(len(seq) for _, seq in records[1:])
    ca = parse_pdb_ca(native_pdb)
    nanobody_structure_ca_residues = len(ca.get("H", {}))
    antigen_ca = [xyz for chain, residues in ca.items() if chain != "H" for xyz in residues.values()]
    nanobody_ca = list(ca.get("H", {}).values())

    contacts = 0
    for nb_xyz in nanobody_ca:
        for ag_xyz in antigen_ca:
            if _distance_lt_10a(nb_xyz, ag_xyz):
                contacts += 1

    if nanobody_sequence_residues == 0:
        reasons.append("empty_nanobody_sequence")
    else:
        missing = 1.0 - (nanobody_structure_ca_residues / nanobody_sequence_residues)
        if missing > 0.5:
            reasons.append("nanobody_structure_missing_residues_gt_50pct")
    if antigen_sequence_residues > 600:
        reasons.append("antigen_sequence_gt_600_residues")
    if contacts == 0:
        reasons.append("no_interfacial_ca_contacts_lt_10a")

    return FilterResult(
        target_id=target_id,
        include=not reasons,
        reasons=tuple(reasons),
        nanobody_sequence_residues=nanobody_sequence_residues,
        nanobody_structure_ca_residues=nanobody_structure_ca_residues,
        antigen_sequence_residues=antigen_sequence_residues,
        interfacial_ca_contacts_lt_10a=contacts,
    )


def final_target_filter(subset_name: str, subset_dir: Path) -> Tuple[List[str], List[FilterResult]]:
    target_ids = read_prot_ids(subset_dir)
    if subset_name != "SAbDab-22H2-NanoAg":
        return target_ids, [FilterResult(target_id, True, ()) for target_id in target_ids]

    results = [evaluate_nanoag_target(subset_dir, target_id) for target_id in target_ids]
    included = [result.target_id for result in results if result.include]
    return included, results


def locate_subset_dir(extracted: Path, subset_name: str) -> Path:
    exact = [path for path in extracted.rglob(subset_name) if path.is_dir()]
    if exact:
        return sorted(exact, key=lambda p: (len(p.parts), p.as_posix()))[0]

    tokens = {
        "SAbDab-22H2-Nano": ["22H2", "Nano"],
        "SAbDab-22H2-NanoAg": ["22H2", "NanoAg"],
    }[subset_name]
    candidates = []
    for path in extracted.rglob("*"):
        if not path.is_dir():
            continue
        lower = path.name.lower()
        if all(token.lower() in lower for token in tokens):
            candidates.append(path)
    if candidates:
        return sorted(candidates, key=lambda p: (len(p.parts), p.as_posix()))[0]

    raise RuntimeError(f"Could not locate official archive subset directory: {subset_name}")


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def canonicalize_subset(source_dir: Path, dest_dir: Path, subset_name: str) -> List[FilterResult]:
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    included, filter_results = final_target_filter(subset_name, source_dir)
    expected = SUBSET_EXPECTED_COUNTS[subset_name]
    if len(included) != expected:
        excluded = [result for result in filter_results if not result.include]
        reasons = "; ".join(
            f"{result.target_id}: {','.join(result.reasons)}" for result in excluded
        )
        raise RuntimeError(
            f"{subset_name} resolved to {len(included)} targets, expected {expected}. "
            f"Excluded targets: {reasons or 'none'}"
        )

    (dest_dir / "prot_ids.txt").write_text("\n".join(included) + "\n", encoding="utf-8")

    for target_id in included:
        info = parse_target_id(target_id)
        for folder, suffix in [
            ("fasta.files", ".fasta"),
            ("pdb.files.native", ".pdb"),
            ("json.files", ".json"),
            ("native_cif", ".cif"),
        ]:
            copy_if_exists(
                source_dir / folder / f"{target_id}{suffix}",
                dest_dir / folder / f"{target_id}{suffix}",
            )

        for chain in info.antigen_chains:
            copy_if_exists(
                source_dir / "msa.files" / f"{info.pdb_id}_{chain}.a3m",
                dest_dir / "msa.files" / f"{info.pdb_id}_{chain}.a3m",
            )

        for path in source_dir.rglob(f"{target_id}.*"):
            if path.is_file() and not (dest_dir / path.relative_to(source_dir)).exists():
                copy_if_exists(path, dest_dir / path.relative_to(source_dir))

    return filter_results


def canonicalize_subsets(paths: DatasetPaths) -> Dict[str, List[FilterResult]]:
    results: Dict[str, List[FilterResult]] = {}
    for subset_name in SUBSET_EXPECTED_COUNTS:
        source_dir = locate_subset_dir(paths.extracted, subset_name)
        dest_dir = paths.subsets / subset_name
        results[subset_name] = canonicalize_subset(source_dir, dest_dir, subset_name)
        print(f"Canonicalized {subset_name}: {dest_dir}", flush=True)
    return results


def run_download(output_root: Path) -> DatasetPaths:
    paths = dataset_paths(output_root)
    ensure_layout(paths)
    write_source_urls(paths)
    record_tfold_head(paths)
    download_archive(paths)
    safe_extract_archive(paths)
    write_archive_file_tree(paths)
    canonicalize_subsets(paths)
    print(f"Archive: {paths.archive}", flush=True)
    print(f"Subsets: {paths.subsets}", flush=True)
    return paths


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download and canonicalize tFold SAbDab-22H2 nanobody benchmarks."
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("."),
        help="Directory under which data/tfold_sabdab22h2 will be created. Defaults to '.'.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        run_download(args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
