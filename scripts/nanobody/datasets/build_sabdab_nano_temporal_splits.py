#!/usr/bin/env python3
"""Build temporal SAbDab-nano train/val/test splits from raw downloads."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from sabdab_nano_temporal_common import (
    COMPLEX_TEST_EXPECTED,
    MANIFEST_FIELDS,
    RCSB_CIF_URL_TEMPLATE,
    RCSB_ENTRY_API_TEMPLATE,
    SABDAB_ARCHIVE_URL,
    SABDAB_SUMMARY_URL,
    SINGLE_TEST_EXPECTED,
    TFOLD_GDRIVE_VIEW_URL,
    clean_value,
    dataset_paths,
    ensure_base_layout,
    is_missing,
    joined_chains,
    link_or_copy,
    manifest_path,
    maybe_rel,
    normalize_pdb_id,
    read_csv_rows,
    relpath,
    split_chains,
    target_chains_for_id,
    validate_summary_tsv,
    validate_tar_gz,
    validate_zip,
    write_archive_file_trees,
    write_csv,
    write_json,
)


TRAIN_END = dt.date(2022, 1, 1)
VAL_END = dt.date(2022, 7, 1)
TEST_END = dt.date(2023, 1, 1)
CONTACT_CUTOFF_A = 5.0
RCSB_WORKERS = 12

ACCEPTED_COMPLEX_ANTIGEN_TYPES = {"protein", "peptide"}

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
}
STANDARD_AA1 = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class PdbParsed:
    seqres: Dict[str, str] = field(default_factory=dict)
    atom_seq: Dict[str, str] = field(default_factory=dict)
    atom_residue_names: Dict[str, Set[str]] = field(default_factory=dict)
    heavy_atoms: Dict[str, List[Tuple[float, float, float]]] = field(default_factory=dict)


@dataclass
class TargetRecord:
    task: str
    split: str
    subset: str
    target_id: str
    pdb_id: str
    release_date: str
    release_date_source: str
    nanobody_chain: str
    heavy_chain: str
    light_chain: str = ""
    antigen_chains: str = ""
    antigen_type: str = ""
    antigen_name: str = ""
    method: str = ""
    resolution: str = ""
    r_factor: str = ""
    r_free: str = ""
    source_summary_path: str = ""
    source_archive: str = ""
    source_url: str = ""
    was_bound_for_single_task: str = ""
    passed_quality_filters: str = "true"
    filter_notes: str = ""
    notes: str = ""
    raw_src: Optional[Path] = None
    chothia_src: Optional[Path] = None
    imgt_src: Optional[Path] = None
    native_cif_src: Optional[Path] = None
    fasta_src: Optional[Path] = None
    json_src: Optional[Path] = None
    sequence: str = ""
    antigen_sequence: str = ""
    metadata: Dict[str, object] = field(default_factory=dict)


def _print_step(message: str) -> None:
    print(f"[sabdab-nano-temporal] {message}", flush=True)


def read_summary_rows(summary_path: Path) -> Tuple[List[Dict[str, str]], Dict[str, str]]:
    if not validate_summary_tsv(summary_path):
        raise RuntimeError(f"Invalid SAbDab summary TSV: {summary_path}")
    with summary_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"SAbDab summary has no header: {summary_path}")
        lower_to_original = {name.strip().lower(): name for name in reader.fieldnames}
        required = {"pdb", "hchain", "lchain", "model", "antigen_chain", "antigen_type", "date"}
        missing = sorted(required - set(lower_to_original))
        if missing:
            raise RuntimeError(f"SAbDab summary missing required columns: {', '.join(missing)}")
        rows: List[Dict[str, str]] = []
        for row in reader:
            normalized = {key: clean_value(row.get(original, "")) for key, original in lower_to_original.items()}
            rows.append(normalized)
    return rows, lower_to_original


def write_sabdab_raw_manifest(paths, rows: Sequence[Mapping[str, str]], release_info: Mapping[str, Mapping[str, str]]) -> None:
    fields = [
        "pdb",
        "Hchain",
        "Lchain",
        "model",
        "antigen_chain",
        "antigen_type",
        "antigen_name",
        "date",
        "rcsb_initial_release_date",
        "release_date_source",
        "is_sabdab_nano_candidate",
        "is_nanoag_candidate",
    ]
    manifest_rows: List[Dict[str, object]] = []
    for row in rows:
        pdb_id = normalize_pdb_id(row.get("pdb", ""))
        hchain = clean_value(row.get("hchain", ""))
        lchain = clean_value(row.get("lchain", ""))
        antigen_chains = split_chains(row.get("antigen_chain", ""))
        info = release_info.get(pdb_id, {})
        manifest_rows.append(
            {
                "pdb": pdb_id,
                "Hchain": hchain,
                "Lchain": lchain,
                "model": row.get("model", ""),
                "antigen_chain": joined_chains(antigen_chains),
                "antigen_type": row.get("antigen_type", ""),
                "antigen_name": row.get("antigen_name", ""),
                "date": row.get("date", ""),
                "rcsb_initial_release_date": info.get("initial_release_date", ""),
                "release_date_source": info.get("source", ""),
                "is_sabdab_nano_candidate": str(bool(hchain and is_missing(lchain))).lower(),
                "is_nanoag_candidate": str(bool(hchain and is_missing(lchain) and antigen_chains)).lower(),
            }
        )
    write_csv(paths.manifests / "sabdab_all_raw_manifest.csv", manifest_rows, fields)


def build_structure_index(structures_root: Path) -> Dict[str, Dict[str, Path]]:
    index: Dict[str, Dict[str, Path]] = {"raw": {}, "chothia": {}, "imgt": {}}
    for variant in index:
        root = structures_root / variant
        if not root.exists():
            continue
        for path in root.rglob("*.pdb"):
            index[variant].setdefault(path.stem.lower(), path)
    return index


def parse_sabdab_date(value: object) -> str:
    text = clean_value(value)
    if not text:
        return ""
    for fmt in ("%m/%d/%y", "%Y-%m-%d", "%m/%d/%Y"):
        try:
            return dt.datetime.strptime(text, fmt).date().isoformat()
        except ValueError:
            pass
    return ""


def _fetch_rcsb_release_date(pdb_id: str) -> Dict[str, str]:
    url = RCSB_ENTRY_API_TEMPLATE.format(pdb_id=pdb_id.upper())
    last_error = ""
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sabdab-nano-temporal/1.0"})
            with urllib.request.urlopen(req, timeout=30) as response:
                payload = json.loads(response.read().decode("utf-8"))
            date_value = (
                payload.get("rcsb_accession_info", {})
                .get("initial_release_date", "")
            )
            if not date_value:
                raise RuntimeError("rcsb_accession_info.initial_release_date missing")
            return {
                "pdb_id": pdb_id,
                "initial_release_date": str(date_value).split("T", 1)[0],
                "source": "rcsb_accession_info.initial_release_date",
                "status": "ok",
                "error": "",
            }
        except Exception as exc:
            last_error = str(exc)
            time.sleep(0.5 * (attempt + 1))
    return {
        "pdb_id": pdb_id,
        "initial_release_date": "",
        "source": "",
        "status": "error",
        "error": last_error,
    }


def resolve_release_dates(paths, pdb_ids: Iterable[str], sabdab_dates: Mapping[str, str]) -> Dict[str, Dict[str, str]]:
    unique_ids = sorted({pdb_id.lower() for pdb_id in pdb_ids if pdb_id})
    cache: Dict[str, Dict[str, str]] = {}
    if paths.rcsb_release_dates_cache.exists():
        cache = json.loads(paths.rcsb_release_dates_cache.read_text(encoding="utf-8"))

    missing = [pdb_id for pdb_id in unique_ids if cache.get(pdb_id, {}).get("status") != "ok"]
    if missing:
        _print_step(f"Resolving RCSB release dates for {len(missing)} PDB IDs.")
        with ThreadPoolExecutor(max_workers=RCSB_WORKERS) as executor:
            futures = {executor.submit(_fetch_rcsb_release_date, pdb_id): pdb_id for pdb_id in missing}
            done = 0
            for future in as_completed(futures):
                pdb_id = futures[future]
                cache[pdb_id] = future.result()
                done += 1
                if done % 250 == 0:
                    _print_step(f"Resolved {done}/{len(missing)} RCSB release-date requests.")

    resolved: Dict[str, Dict[str, str]] = {}
    for pdb_id in unique_ids:
        info = dict(cache.get(pdb_id, {}))
        if info.get("status") != "ok" or not info.get("initial_release_date"):
            fallback = parse_sabdab_date(sabdab_dates.get(pdb_id, ""))
            if fallback:
                info = {
                    "pdb_id": pdb_id,
                    "initial_release_date": fallback,
                    "source": "sabdab_date_fallback",
                    "status": "fallback",
                    "error": info.get("error", "rcsb release date unavailable"),
                }
            else:
                info = {
                    "pdb_id": pdb_id,
                    "initial_release_date": "",
                    "source": "",
                    "status": "error",
                    "error": info.get("error", "no release date available"),
                }
        info["sabdab_date"] = sabdab_dates.get(pdb_id, "")
        resolved[pdb_id] = info
        cache[pdb_id] = {key: value for key, value in info.items() if key != "sabdab_date"}

    paths.rcsb_release_dates_cache.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(
        paths.rcsb_release_dates_csv,
        [resolved[pdb_id] for pdb_id in sorted(resolved)],
        ["pdb_id", "initial_release_date", "source", "status", "error", "sabdab_date"],
    )
    return resolved


def release_split(date_text: str) -> str:
    if not date_text:
        return ""
    date_value = dt.date.fromisoformat(date_text)
    if date_value < TRAIN_END:
        return "train_pre2022"
    if TRAIN_END <= date_value < VAL_END:
        return "val_2022h1"
    return ""


def parse_pdb(path: Optional[Path]) -> PdbParsed:
    parsed = PdbParsed()
    if path is None or not path.exists():
        return parsed

    seqres_chunks: Dict[str, List[str]] = {}
    atom_order: Dict[str, List[Tuple[str, str, str]]] = {}
    seen_residues: Dict[str, Set[Tuple[str, str, str]]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("SEQRES"):
                chain = line[11].strip() or "_"
                residues = [AA3_TO_AA1.get(res.upper(), "X") for res in line[19:70].split()]
                seqres_chunks.setdefault(chain, []).extend(residues)
                continue
            is_atom = line.startswith("ATOM")
            is_hetatm = line.startswith("HETATM")
            if not (is_atom or is_hetatm):
                continue
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip().upper()
            if is_hetatm and residue_name not in AA3_TO_AA1:
                continue
            chain = line[21].strip() or "_"
            residue_key = (line[22:26].strip(), line[26].strip(), residue_name)
            aa = AA3_TO_AA1.get(residue_name, "X")
            parsed.atom_residue_names.setdefault(chain, set()).add(residue_name)
            if residue_key not in seen_residues.setdefault(chain, set()):
                seen_residues[chain].add(residue_key)
                atom_order.setdefault(chain, []).append(residue_key)
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = atom_name[:1].upper()
            if element == "H" or atom_name.upper().startswith("H"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            parsed.heavy_atoms.setdefault(chain, []).append(xyz)

    parsed.seqres = {chain: "".join(seq) for chain, seq in seqres_chunks.items()}
    parsed.atom_seq = {
        chain: "".join(AA3_TO_AA1.get(residue_name.upper(), "X") for _, _, residue_name in residues)
        for chain, residues in atom_order.items()
    }
    return parsed


def sequence_for_chain(parsed: PdbParsed, chain: str) -> str:
    return parsed.seqres.get(chain) or parsed.atom_seq.get(chain, "")


def sequence_for_chains(parsed: PdbParsed, chains: Sequence[str]) -> str:
    return "|".join(sequence_for_chain(parsed, chain) for chain in chains if sequence_for_chain(parsed, chain))


def count_interfacial_contacts(parsed: PdbParsed, nanobody_chain: str, antigen_chains: Sequence[str]) -> int:
    nb_atoms = parsed.heavy_atoms.get(nanobody_chain, [])
    antigen_atoms: List[Tuple[float, float, float]] = []
    for chain in antigen_chains:
        antigen_atoms.extend(parsed.heavy_atoms.get(chain, []))
    if not nb_atoms or not antigen_atoms:
        return 0

    cutoff2 = CONTACT_CUTOFF_A * CONTACT_CUTOFF_A
    cell_size = CONTACT_CUTOFF_A
    grid: Dict[Tuple[int, int, int], List[Tuple[float, float, float]]] = {}
    for xyz in antigen_atoms:
        cell = tuple(int(math.floor(coord / cell_size)) for coord in xyz)
        grid.setdefault(cell, []).append(xyz)

    contacts = 0
    for x, y, z in nb_atoms:
        base = (int(math.floor(x / cell_size)), int(math.floor(y / cell_size)), int(math.floor(z / cell_size)))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for ax, ay, az in grid.get((base[0] + dx, base[1] + dy, base[2] + dz), []):
                        if (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2 <= cutoff2:
                            contacts += 1
                            break
                    if contacts:
                        break
                if contacts:
                    break
            if contacts:
                break
    return contacts


def tsv_join(values: Iterable[str]) -> str:
    return ";".join(sorted({clean_value(value) for value in values if clean_value(value)}))


def row_antigen_types(row: Mapping[str, str]) -> Set[str]:
    text = clean_value(row.get("antigen_type", ""))
    if not text:
        return set()
    return {item.strip().lower() for item in text.replace("|", ",").replace("/", ",").split(",") if item.strip()}


def locate_tfold_subset(extracted: Path, subset_name: str) -> Path:
    exact = [path for path in extracted.rglob(subset_name) if path.is_dir()]
    if exact:
        return sorted(exact, key=lambda path: (len(path.parts), path.as_posix()))[0]
    tokens = ["22H2", "NanoAg"] if subset_name.endswith("NanoAg") else ["22H2", "Nano"]
    candidates = []
    for path in extracted.rglob("*"):
        if path.is_dir() and all(token.lower() in path.name.lower() for token in tokens):
            candidates.append(path)
    if candidates:
        return sorted(candidates, key=lambda path: (len(path.parts), path.as_posix()))[0]
    raise RuntimeError(f"Could not locate official tFold subset directory: {subset_name}")


def read_prot_ids(subset_dir: Path) -> List[str]:
    prot_ids = subset_dir / "prot_ids.txt"
    if prot_ids.exists():
        return [line.strip() for line in prot_ids.read_text(encoding="utf-8").splitlines() if line.strip()]
    ids: Set[str] = set()
    for suffix in [".fasta", ".pdb", ".cif", ".json"]:
        for path in subset_dir.rglob(f"*{suffix}"):
            ids.add(path.stem)
    return sorted(ids)


def parse_target_id(target_id: str) -> Tuple[str, List[str], List[str]]:
    parts = target_id.split("_")
    pdb_id = parts[0].lower() if parts else ""
    if "NA" in parts:
        idx = parts.index("NA")
        return pdb_id, parts[1:idx], parts[idx + 1 :]
    return pdb_id, parts[1:2], []


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    if not path.exists():
        return []
    records: List[Tuple[str, str]] = []
    name = ""
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if name:
        records.append((name, "".join(chunks)))
    return records


def write_fasta(path: Path, records: Sequence[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    for name, seq in records:
        if not seq:
            continue
        lines.append(f">{name}")
        for idx in range(0, len(seq), 80):
            lines.append(seq[idx : idx + 80])
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def best_structure_sources(index: Mapping[str, Mapping[str, Path]], pdb_id: str) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    key = pdb_id.lower()
    return (
        index.get("raw", {}).get(key),
        index.get("chothia", {}).get(key),
        index.get("imgt", {}).get(key),
    )


def build_single_records(
    paths,
    rows: Sequence[Mapping[str, str]],
    release_info: Mapping[str, Mapping[str, str]],
    structure_index: Mapping[str, Mapping[str, Path]],
) -> List[TargetRecord]:
    by_target: Dict[str, TargetRecord] = {}
    parsed_cache: Dict[Path, PdbParsed] = {}
    for row in rows:
        pdb_id = normalize_pdb_id(row.get("pdb", ""))
        hchain = clean_value(row.get("hchain", ""))
        lchain = clean_value(row.get("lchain", ""))
        if not pdb_id or not hchain or not is_missing(lchain):
            continue
        info = release_info.get(pdb_id, {})
        release_date = info.get("initial_release_date", "")
        split = release_split(release_date)
        if split not in {"train_pre2022", "val_2022h1"}:
            continue
        target_id = f"{pdb_id}_{hchain}"
        raw_src, chothia_src, imgt_src = best_structure_sources(structure_index, pdb_id)
        parse_src = raw_src or chothia_src or imgt_src
        if parse_src and parse_src not in parsed_cache:
            parsed_cache[parse_src] = parse_pdb(parse_src)
        sequence = sequence_for_chain(parsed_cache.get(parse_src, PdbParsed()), hchain) if parse_src else ""
        existing = by_target.get(target_id)
        antigen_chains = split_chains(row.get("antigen_chain", ""))
        if existing:
            merged_chains = split_chains(existing.antigen_chains) + antigen_chains
            existing.antigen_chains = joined_chains(merged_chains)
            existing.antigen_type = tsv_join([existing.antigen_type, row.get("antigen_type", "")])
            existing.antigen_name = tsv_join([existing.antigen_name, row.get("antigen_name", "")])
            existing.was_bound_for_single_task = str(bool(split_chains(existing.antigen_chains))).lower()
            if sequence and not existing.sequence:
                existing.sequence = sequence
            continue
        by_target[target_id] = TargetRecord(
            task="single_nano",
            split=split,
            subset=f"single_nano/{split}",
            target_id=target_id,
            pdb_id=pdb_id,
            release_date=release_date,
            release_date_source=info.get("source", ""),
            nanobody_chain=hchain,
            heavy_chain=hchain,
            light_chain="",
            antigen_chains=joined_chains(antigen_chains),
            antigen_type=clean_value(row.get("antigen_type", "")),
            antigen_name=clean_value(row.get("antigen_name", "")),
            method=clean_value(row.get("method", "")),
            resolution=clean_value(row.get("resolution", "")),
            r_factor=clean_value(row.get("r_factor", "")),
            r_free=clean_value(row.get("r_free", "")),
            source_summary_path=relpath(paths.sabdab_summary, paths.output_root),
            source_archive=relpath(paths.sabdab_archive, paths.output_root),
            source_url=SABDAB_SUMMARY_URL,
            was_bound_for_single_task=str(bool(antigen_chains)).lower(),
            passed_quality_filters="true",
            raw_src=raw_src,
            chothia_src=chothia_src,
            imgt_src=imgt_src,
            sequence=sequence,
            metadata={"sabdab_row": dict(row)},
        )
    return sorted(by_target.values(), key=lambda item: item.target_id)


def evaluate_complex_quality(record: TargetRecord) -> Tuple[bool, List[str], Dict[str, object]]:
    notes: List[str] = []
    metadata: Dict[str, object] = {}
    parsed = parse_pdb(record.raw_src or record.chothia_src or record.imgt_src)
    if not (record.raw_src or record.chothia_src or record.imgt_src):
        return False, ["missing_sabdab_structure_file"], metadata

    nb_chain = record.nanobody_chain
    antigen_chains = split_chains(record.antigen_chains)
    nb_sequence = sequence_for_chain(parsed, nb_chain)
    antigen_sequence = sequence_for_chains(parsed, antigen_chains)
    record.sequence = nb_sequence
    record.antigen_sequence = antigen_sequence
    metadata["nanobody_sequence_length"] = len(nb_sequence)
    metadata["antigen_sequence_length"] = sum(len(sequence_for_chain(parsed, chain)) for chain in antigen_chains)

    nonstandard = sorted(
        residue
        for residue in parsed.atom_residue_names.get(nb_chain, set())
        if AA3_TO_AA1.get(residue, "X") not in STANDARD_AA1
    )
    if nonstandard:
        notes.append(f"nonstandard_nanobody_residues:{','.join(nonstandard)}")

    contacts = count_interfacial_contacts(parsed, nb_chain, antigen_chains)
    metadata["heavy_atom_contacts_le_5a"] = contacts
    if contacts == 0:
        notes.append("no_interfacial_heavy_atom_contacts_le_5a")

    antigen_length = metadata["antigen_sequence_length"]
    if isinstance(antigen_length, int) and antigen_length > 600:
        notes.append("antigen_sequence_gt_600_residues")

    seqres_nb = parsed.seqres.get(nb_chain, "")
    atom_nb = parsed.atom_seq.get(nb_chain, "")
    if seqres_nb:
        observed = len(atom_nb)
        missing_fraction = 1.0 - (observed / len(seqres_nb)) if seqres_nb else 0.0
        metadata["missing_nanobody_residue_fraction"] = round(missing_fraction, 4)
        metadata["missing_fraction_status"] = "computed"
        if missing_fraction > 0.5:
            notes.append("nanobody_missing_observed_residues_gt_50pct")
    else:
        metadata["missing_fraction_status"] = "unknown_no_seqres"

    return not notes, notes, metadata


def build_complex_records(
    paths,
    rows: Sequence[Mapping[str, str]],
    release_info: Mapping[str, Mapping[str, str]],
    structure_index: Mapping[str, Mapping[str, Path]],
) -> Tuple[List[TargetRecord], List[Dict[str, object]]]:
    records: List[TargetRecord] = []
    filtered: List[Dict[str, object]] = []
    seen_targets: Set[str] = set()

    for row in rows:
        pdb_id = normalize_pdb_id(row.get("pdb", ""))
        hchain = clean_value(row.get("hchain", ""))
        lchain = clean_value(row.get("lchain", ""))
        antigen_chains = split_chains(row.get("antigen_chain", ""))
        antigen_types = row_antigen_types(row)
        if not pdb_id or not hchain or not is_missing(lchain) or not antigen_chains:
            continue
        target_id = f"{pdb_id}_{hchain}_NA_{target_chains_for_id(antigen_chains)}"
        base_filtered = {
            "target_id": target_id,
            "pdb_id": pdb_id,
            "nanobody_chain": hchain,
            "antigen_chains": joined_chains(antigen_chains),
            "antigen_type": clean_value(row.get("antigen_type", "")),
        }
        if not antigen_types or not antigen_types.issubset(ACCEPTED_COMPLEX_ANTIGEN_TYPES):
            filtered.append({**base_filtered, "reason": "unsupported_or_missing_antigen_type"})
            continue
        if target_id in seen_targets:
            filtered.append({**base_filtered, "reason": "duplicate_target_id"})
            continue
        seen_targets.add(target_id)

        info = release_info.get(pdb_id, {})
        release_date = info.get("initial_release_date", "")
        split = release_split(release_date)
        if split not in {"train_pre2022", "val_2022h1"}:
            if release_date:
                filtered.append({**base_filtered, "reason": f"outside_train_val_temporal_window:{release_date}"})
            else:
                filtered.append({**base_filtered, "reason": "missing_release_date"})
            continue

        raw_src, chothia_src, imgt_src = best_structure_sources(structure_index, pdb_id)
        record = TargetRecord(
            task="complex_nanoag",
            split=split,
            subset=f"complex_nanoag/{split}",
            target_id=target_id,
            pdb_id=pdb_id,
            release_date=release_date,
            release_date_source=info.get("source", ""),
            nanobody_chain=hchain,
            heavy_chain=hchain,
            light_chain="",
            antigen_chains=joined_chains(antigen_chains),
            antigen_type=clean_value(row.get("antigen_type", "")),
            antigen_name=clean_value(row.get("antigen_name", "")),
            method=clean_value(row.get("method", "")),
            resolution=clean_value(row.get("resolution", "")),
            r_factor=clean_value(row.get("r_factor", "")),
            r_free=clean_value(row.get("r_free", "")),
            source_summary_path=relpath(paths.sabdab_summary, paths.output_root),
            source_archive=relpath(paths.sabdab_archive, paths.output_root),
            source_url=SABDAB_SUMMARY_URL,
            was_bound_for_single_task="",
            raw_src=raw_src,
            chothia_src=chothia_src,
            imgt_src=imgt_src,
            metadata={"sabdab_row": dict(row)},
        )
        passed, notes, metadata = evaluate_complex_quality(record)
        record.metadata.update(metadata)
        record.passed_quality_filters = str(passed).lower()
        record.filter_notes = ";".join(notes)
        if passed:
            records.append(record)
        else:
            filtered.append({**base_filtered, **metadata, "reason": record.filter_notes})
    return sorted(records, key=lambda item: item.target_id), filtered


def evaluate_tfold_nanoag_filter(subset_dir: Path, target_id: str) -> Tuple[bool, List[str], Dict[str, object]]:
    fasta = subset_dir / "fasta.files" / f"{target_id}.fasta"
    native_pdb = subset_dir / "pdb.files.native" / f"{target_id}.pdb"
    records = read_fasta(fasta)
    parsed = parse_pdb(native_pdb)
    metadata: Dict[str, object] = {}
    reasons: List[str] = []
    if len(records) < 2:
        reasons.append("fasta_missing_antigen")
    nb_len = len(records[0][1]) if records else 0
    ag_len = sum(len(seq) for _, seq in records[1:])
    metadata["nanobody_sequence_length"] = nb_len
    metadata["antigen_sequence_length"] = ag_len
    nb_atom_len = len(parsed.atom_seq.get("H", ""))
    metadata["nanobody_observed_residue_count_native_pdb_chain_H"] = nb_atom_len
    if nb_len:
        missing_fraction = 1.0 - (nb_atom_len / nb_len)
        metadata["missing_nanobody_residue_fraction"] = round(missing_fraction, 4)
        if missing_fraction > 0.5:
            reasons.append("nanobody_structure_missing_residues_gt_50pct")
    if ag_len > 600:
        reasons.append("antigen_sequence_gt_600_residues")
    return not reasons, reasons, metadata


def build_tfold_test_records(paths, release_info: Mapping[str, Mapping[str, str]]) -> Tuple[List[TargetRecord], List[TargetRecord], List[Dict[str, object]], bool]:
    single_dir = locate_tfold_subset(paths.tfold_extracted, "SAbDab-22H2-Nano")
    complex_dir = locate_tfold_subset(paths.tfold_extracted, "SAbDab-22H2-NanoAg")
    single_ids = read_prot_ids(single_dir)
    complex_ids_archive = read_prot_ids(complex_dir)
    filtered: List[Dict[str, object]] = []

    if len(single_ids) != SINGLE_TEST_EXPECTED:
        raise RuntimeError(f"Official SAbDab-22H2-Nano target count is {len(single_ids)}, expected {SINGLE_TEST_EXPECTED}")

    complex_ids: List[str] = []
    for target_id in complex_ids_archive:
        passed, reasons, metadata = evaluate_tfold_nanoag_filter(complex_dir, target_id)
        if passed:
            complex_ids.append(target_id)
        else:
            filtered.append({"target_id": target_id, "reason": ";".join(reasons), **metadata})
    if len(complex_ids) != COMPLEX_TEST_EXPECTED:
        raise RuntimeError(
            f"Official SAbDab-22H2-NanoAg resolved to {len(complex_ids)} targets, expected {COMPLEX_TEST_EXPECTED}. "
            f"Archive target count: {len(complex_ids_archive)}"
        )

    native_structures_present = True
    single_records: List[TargetRecord] = []
    complex_records: List[TargetRecord] = []
    for task, split, subset, subset_dir, ids, output in [
        ("single_nano", "test_2022h2", "SAbDab-22H2-Nano", single_dir, single_ids, single_records),
        ("complex_nanoag", "test_2022h2", "SAbDab-22H2-NanoAg", complex_dir, complex_ids, complex_records),
    ]:
        for target_id in sorted(ids):
            pdb_id, nb_chains, antigen_chains = parse_target_id(target_id)
            info = release_info.get(pdb_id, {})
            fasta_src = subset_dir / "fasta.files" / f"{target_id}.fasta"
            native_pdb = subset_dir / "pdb.files.native" / f"{target_id}.pdb"
            native_cif = subset_dir / "native_cif" / f"{target_id}.cif"
            json_src = subset_dir / "json.files" / f"{target_id}.json"
            if not native_pdb.exists() and not native_cif.exists():
                native_structures_present = False
                fetched = fetch_rcsb_cif(paths, pdb_id)
                native_cif = fetched
            records = read_fasta(fasta_src)
            nb_seq = records[0][1] if records else ""
            ag_seq = "|".join(seq for _, seq in records[1:])
            output.append(
                TargetRecord(
                    task=task,
                    split=split,
                    subset=subset,
                    target_id=target_id,
                    pdb_id=pdb_id,
                    release_date=info.get("initial_release_date", ""),
                    release_date_source=info.get("source", ""),
                    nanobody_chain=joined_chains(nb_chains),
                    heavy_chain=joined_chains(nb_chains),
                    light_chain="",
                    antigen_chains=joined_chains(antigen_chains),
                    antigen_type="",
                    antigen_name="",
                    raw_src=native_pdb if native_pdb.exists() else None,
                    native_cif_src=native_cif if native_cif.exists() else None,
                    fasta_src=fasta_src if fasta_src.exists() else None,
                    json_src=json_src if json_src.exists() else None,
                    source_archive=relpath(paths.tfold_archive, paths.output_root),
                    source_url=TFOLD_GDRIVE_VIEW_URL,
                    passed_quality_filters="true",
                    notes="official tFold H2 test target",
                    sequence=nb_seq,
                    antigen_sequence=ag_seq,
                    metadata={"source_tfold_subset": str(subset_dir), "fasta_records": [name for name, _ in records]},
                )
            )
    return single_records, complex_records, filtered, native_structures_present


def fetch_rcsb_cif(paths, pdb_id: str) -> Path:
    output = paths.rcsb_downloaded_cif / f"{pdb_id.lower()}.cif"
    if output.exists() and output.stat().st_size > 0:
        return output
    url = RCSB_CIF_URL_TEMPLATE.format(pdb_id=pdb_id.upper())
    req = urllib.request.Request(url, headers={"User-Agent": "sabdab-nano-temporal/1.0"})
    with urllib.request.urlopen(req, timeout=120) as response:
        payload = response.read()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    return output


def remove_leakage(
    single_train: List[TargetRecord],
    single_val: List[TargetRecord],
    single_test: List[TargetRecord],
    complex_train: List[TargetRecord],
    complex_val: List[TargetRecord],
    complex_test: List[TargetRecord],
) -> Tuple[List[TargetRecord], List[TargetRecord], List[TargetRecord], List[TargetRecord], Dict[str, object]]:
    report: Dict[str, object] = {
        "train_test_target_overlap": 0,
        "val_test_target_overlap": 0,
        "train_val_target_overlap": 0,
        "exact_nanobody_sequence_leaks_removed": 0,
        "exact_antigen_sequence_leaks_removed": 0,
        "removed_targets": [],
    }

    test_ids = {record.target_id for record in single_test + complex_test}
    val_ids = {record.target_id for record in single_val + complex_val}
    train_ids = {record.target_id for record in single_train + complex_train}
    report["train_test_target_overlap"] = len(train_ids & test_ids)
    report["val_test_target_overlap"] = len(val_ids & test_ids)
    report["train_val_target_overlap"] = len(train_ids & val_ids)

    def remove_target_id_overlaps(records: List[TargetRecord], forbidden: Set[str], reason: str) -> List[TargetRecord]:
        kept = []
        for record in records:
            if record.target_id in forbidden:
                report["removed_targets"].append({"target_id": record.target_id, "reason": reason})
            else:
                kept.append(record)
        return kept

    single_train = remove_target_id_overlaps(single_train, test_ids, "target_id_overlap_with_test")
    single_val = remove_target_id_overlaps(single_val, test_ids, "target_id_overlap_with_test")
    complex_train = remove_target_id_overlaps(complex_train, test_ids, "target_id_overlap_with_test")
    complex_val = remove_target_id_overlaps(complex_val, test_ids, "target_id_overlap_with_test")

    def seq_set(records: Sequence[TargetRecord]) -> Set[str]:
        return {record.sequence for record in records if record.sequence}

    single_test_seq = seq_set(single_test)
    single_val = [
        record for record in single_val
        if not _remove_if_sequence(record, single_test_seq, report, "nanobody_sequence_overlap_with_single_test")
    ]
    single_val_seq = seq_set(single_val)
    single_train = [
        record for record in single_train
        if not _remove_if_sequence(record, single_test_seq | single_val_seq, report, "nanobody_sequence_overlap_with_single_val_or_test")
    ]

    complex_test_seq = seq_set(complex_test)
    complex_val = [
        record for record in complex_val
        if not _remove_if_sequence(record, complex_test_seq, report, "nanobody_sequence_overlap_with_complex_test")
    ]
    complex_val_seq = seq_set(complex_val)
    complex_train = [
        record for record in complex_train
        if not _remove_if_sequence(record, complex_test_seq | complex_val_seq, report, "nanobody_sequence_overlap_with_complex_val_or_test")
    ]

    train_antigen_sequences = {record.antigen_sequence for record in complex_train if record.antigen_sequence}
    filtered_complex_val = []
    for record in complex_val:
        if record.antigen_sequence and record.antigen_sequence in train_antigen_sequences:
            report["exact_antigen_sequence_leaks_removed"] += 1
            report["removed_targets"].append({"target_id": record.target_id, "reason": "complex_validation_antigen_sequence_overlap_with_train"})
        else:
            filtered_complex_val.append(record)
    complex_val = filtered_complex_val

    return single_train, single_val, complex_train, complex_val, report


def _remove_if_sequence(record: TargetRecord, forbidden_sequences: Set[str], report: Dict[str, object], reason: str) -> bool:
    if record.sequence and record.sequence in forbidden_sequences:
        report["exact_nanobody_sequence_leaks_removed"] += 1
        report["removed_targets"].append({"target_id": record.target_id, "reason": reason})
        return True
    return False


def reset_split_dirs(paths) -> None:
    for task, split in [
        ("single_nano", "train_pre2022"),
        ("single_nano", "val_2022h1"),
        ("single_nano", "test_2022h2"),
        ("complex_nanoag", "train_pre2022"),
        ("complex_nanoag", "val_2022h1"),
        ("complex_nanoag", "test_2022h2"),
    ]:
        root = paths.splits / task / split
        if root.exists():
            shutil.rmtree(root)
        (root / "structures").mkdir(parents=True, exist_ok=True)
        (root / "metadata").mkdir(parents=True, exist_ok=True)


def materialize_records(paths, records: Sequence[TargetRecord]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for record in records:
        split_root = paths.splits / record.task / record.split
        structure_dir = split_root / "structures" / record.target_id
        metadata_dir = split_root / "metadata"
        raw_dst = structure_dir / "raw.pdb"
        chothia_dst = structure_dir / "chothia.pdb"
        imgt_dst = structure_dir / "imgt.pdb"
        native_cif_dst = structure_dir / "native.cif"
        fasta_dst = metadata_dir / f"{record.target_id}.fasta"
        json_dst = metadata_dir / f"{record.target_id}.json"

        linked_raw = link_or_copy(record.raw_src, raw_dst) if record.raw_src else False
        linked_chothia = link_or_copy(record.chothia_src, chothia_dst) if record.chothia_src else False
        linked_imgt = link_or_copy(record.imgt_src, imgt_dst) if record.imgt_src else False
        linked_cif = link_or_copy(record.native_cif_src, native_cif_dst) if record.native_cif_src else False

        if record.fasta_src and record.fasta_src.exists():
            link_or_copy(record.fasta_src, fasta_dst)
        elif record.sequence or record.antigen_sequence:
            fasta_records = [(record.nanobody_chain or record.heavy_chain or "nanobody", record.sequence)]
            antigen_chains = split_chains(record.antigen_chains)
            antigen_sequences = record.antigen_sequence.split("|") if record.antigen_sequence else []
            for chain, seq in zip(antigen_chains, antigen_sequences):
                fasta_records.append((chain, seq))
            write_fasta(fasta_dst, fasta_records)

        if record.json_src and record.json_src.exists():
            link_or_copy(record.json_src, json_dst)
        else:
            metadata = {
                "target_id": record.target_id,
                "task": record.task,
                "split": record.split,
                "pdb_id": record.pdb_id,
                "nanobody_chain": record.nanobody_chain,
                "antigen_chains": record.antigen_chains,
                "release_date": record.release_date,
                "release_date_source": record.release_date_source,
                **record.metadata,
            }
            write_json(json_dst, metadata)

        row = {
            "task": record.task,
            "split": record.split,
            "subset": record.subset,
            "target_id": record.target_id,
            "pdb_id": record.pdb_id,
            "release_date": record.release_date,
            "release_date_source": record.release_date_source,
            "nanobody_chain": record.nanobody_chain,
            "heavy_chain": record.heavy_chain,
            "light_chain": record.light_chain,
            "antigen_chains": record.antigen_chains,
            "antigen_type": record.antigen_type,
            "antigen_name": record.antigen_name,
            "method": record.method,
            "resolution": record.resolution,
            "r_factor": record.r_factor,
            "r_free": record.r_free,
            "raw_pdb_path": relpath(raw_dst, paths.output_root) if linked_raw else "",
            "chothia_pdb_path": relpath(chothia_dst, paths.output_root) if linked_chothia else "",
            "imgt_pdb_path": relpath(imgt_dst, paths.output_root) if linked_imgt else "",
            "native_cif_path": relpath(native_cif_dst, paths.output_root) if linked_cif else "",
            "fasta_path": relpath(fasta_dst, paths.output_root) if fasta_dst.exists() else "",
            "json_path": relpath(json_dst, paths.output_root) if json_dst.exists() else "",
            "source_summary_path": record.source_summary_path,
            "source_archive": record.source_archive,
            "source_url": record.source_url,
            "was_bound_for_single_task": record.was_bound_for_single_task,
            "passed_quality_filters": record.passed_quality_filters,
            "filter_notes": record.filter_notes,
            "notes": record.notes,
        }
        rows.append(row)
    return rows


def write_split_manifests(paths, split_rows: Mapping[str, Sequence[Mapping[str, object]]]) -> None:
    all_rows: List[Dict[str, object]] = []
    for split_name, rows in split_rows.items():
        write_csv(manifest_path(paths, split_name), rows, MANIFEST_FIELDS)
        for row in rows:
            all_rows.append(dict(row))
    all_fields = ["dataset_row_id"] + MANIFEST_FIELDS
    all_with_ids = []
    for idx, row in enumerate(all_rows, start=1):
        all_with_ids.append({"dataset_row_id": f"sabdab_nano_temporal_{idx:07d}", **row})
    write_csv(paths.manifests / "all_splits_manifest.csv", all_with_ids, all_fields)


def build_splits(output_root: Path) -> Dict[str, object]:
    paths = dataset_paths(output_root)
    ensure_base_layout(paths)
    if not validate_zip(paths.sabdab_archive):
        raise RuntimeError(f"Missing or invalid SAbDab archive: {paths.sabdab_archive}")
    if not validate_summary_tsv(paths.sabdab_summary):
        raise RuntimeError(f"Missing or invalid SAbDab summary: {paths.sabdab_summary}")
    if not validate_tar_gz(paths.tfold_archive):
        raise RuntimeError(f"Missing or invalid tFold archive: {paths.tfold_archive}")

    rows, _ = read_summary_rows(paths.sabdab_summary)
    structure_index = build_structure_index(paths.sabdab_structures)
    if not any(structure_index.values()):
        raise RuntimeError(f"No extracted SAbDab PDB files found under {paths.sabdab_structures}")

    single_dir = locate_tfold_subset(paths.tfold_extracted, "SAbDab-22H2-Nano")
    complex_dir = locate_tfold_subset(paths.tfold_extracted, "SAbDab-22H2-NanoAg")
    test_pdb_ids = [parse_target_id(target_id)[0] for target_id in read_prot_ids(single_dir) + read_prot_ids(complex_dir)]

    candidate_rows = [
        row for row in rows
        if normalize_pdb_id(row.get("pdb", "")) and clean_value(row.get("hchain", "")) and is_missing(row.get("lchain", ""))
    ]
    sabdab_dates = {
        normalize_pdb_id(row.get("pdb", "")): clean_value(row.get("date", ""))
        for row in rows
        if normalize_pdb_id(row.get("pdb", ""))
    }
    release_info = resolve_release_dates(
        paths,
        [normalize_pdb_id(row.get("pdb", "")) for row in candidate_rows] + test_pdb_ids,
        sabdab_dates,
    )
    write_sabdab_raw_manifest(paths, rows, release_info)

    _print_step("Building SAbDab single-nanobody train/val targets.")
    single_records = build_single_records(paths, rows, release_info, structure_index)
    single_train = [record for record in single_records if record.split == "train_pre2022"]
    single_val = [record for record in single_records if record.split == "val_2022h1"]

    _print_step("Building SAbDab nanobody-antigen train/val targets with quality filters.")
    complex_records, complex_filtered = build_complex_records(paths, rows, release_info, structure_index)
    complex_train = [record for record in complex_records if record.split == "train_pre2022"]
    complex_val = [record for record in complex_records if record.split == "val_2022h1"]

    _print_step("Extracting official tFold H2 test subsets.")
    single_test, complex_test, tfold_filtered, native_structures_present = build_tfold_test_records(paths, release_info)

    single_train, single_val, complex_train, complex_val, leakage_report = remove_leakage(
        single_train,
        single_val,
        single_test,
        complex_train,
        complex_val,
        complex_test,
    )

    reset_split_dirs(paths)
    split_rows = {
        "single_nano_train_pre2022": materialize_records(paths, single_train),
        "single_nano_val_2022h1": materialize_records(paths, single_val),
        "single_nano_test_2022h2": materialize_records(paths, single_test),
        "complex_nanoag_train_pre2022": materialize_records(paths, complex_train),
        "complex_nanoag_val_2022h1": materialize_records(paths, complex_val),
        "complex_nanoag_test_2022h2": materialize_records(paths, complex_test),
    }
    write_split_manifests(paths, split_rows)

    write_csv(
        paths.manifests / "complex_nanoag_filtered_out.csv",
        complex_filtered,
        [
            "target_id",
            "pdb_id",
            "nanobody_chain",
            "antigen_chains",
            "antigen_type",
            "reason",
            "nanobody_sequence_length",
            "antigen_sequence_length",
            "heavy_atom_contacts_le_5a",
            "missing_nanobody_residue_fraction",
            "missing_fraction_status",
        ],
    )
    write_csv(
        paths.manifests / "complex_nanoag_test_archive_filtered_out.csv",
        tfold_filtered,
        [
            "target_id",
            "reason",
            "nanobody_sequence_length",
            "antigen_sequence_length",
            "nanobody_observed_residue_count_native_pdb_chain_H",
            "missing_nanobody_residue_fraction",
        ],
    )
    write_json(
        paths.manifests / "leakage_report.json",
        {
            **leakage_report,
            "native_structures_present_in_tfold_archive": native_structures_present,
        },
    )
    write_archive_file_trees(paths)

    summary = {
        "counts": {split: len(rows) for split, rows in split_rows.items()},
        "native_structures_present_in_tfold_archive": native_structures_present,
        "leakage": leakage_report,
        "complex_filtered_out": len(complex_filtered),
        "tfold_test_archive_filtered_out": len(tfold_filtered),
    }
    write_json(paths.manifests / "build_summary.json", summary)
    for split, count in summary["counts"].items():
        _print_step(f"{split}: {count} targets")
    return summary


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build SAbDab-nano temporal train/val/test splits.")
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
        build_splits(args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
