"""Manifest helpers for nanobody data preparation."""

from __future__ import annotations

import math
import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any


CANONICAL_MANIFEST_COLUMNS = [
    "pdb_id",
    "release_date",
    "split",
    "split_name",
    "task_type",
    "nanobody_chain",
    "antigen_chains",
    "structure_path",
    "prepared_mmcif_path",
    "cdrh3_start",
    "cdrh3_end",
    "cdrh3_source",
    "num_tokens",
    "notes",
]

PROTENIX_INDEX_COLUMNS = [
    "type",
    "pdb_id",
    "cluster_id",
    "assembly_id",
    "release_date",
    "resolution",
    "num_tokens",
    "num_prot_chains",
    "eval_type",
    "entity_1_id",
    "chain_1_id",
    "mol_1_type",
    "sub_mol_1_type",
    "cluster_1_id",
    "entity_2_id",
    "chain_2_id",
    "mol_2_type",
    "sub_mol_2_type",
    "cluster_2_id",
]

COLUMN_ALIASES = {
    "pdb_id": ["pdb", "pdb_id", "PDB"],
    "release_date": ["release_date", "date", "deposition_date"],
    "nanobody_chain": [
        "nanobody_chain",
        "heavy_chain",
        "Hchain",
        "H_chain",
        "vh_chain",
    ],
    "antigen_chains": [
        "antigen_chain",
        "antigen_chains",
        "antigen_chain_ids",
        "antigen",
    ],
    "structure_path": ["structure_path", "pdb_path", "mmcif_path", "cif_path"],
    "cdrh3_start": ["cdrh3_start", "cdr3_start"],
    "cdrh3_end": ["cdrh3_end", "cdr3_end"],
}


@dataclass(frozen=True)
class SplitPolicy:
    """Date cutoffs used to assign nanobody records to train/val/test splits."""

    train_cutoff: date = date(2022, 1, 1)
    val_cutoff: date = date(2022, 7, 1)
    test_cutoff: date = date(2023, 1, 1)
    strict_protenix_cutoff: bool = False

    @property
    def effective_train_cutoff(self) -> date:
        if self.strict_protenix_cutoff:
            return date(2021, 9, 30)
        return self.train_cutoff


def _clean_key(key: str) -> str:
    return re.sub(r"[^a-z0-9]", "", key.lower())


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "null", "-"}


def _find_column(columns: Iterable[str], canonical: str) -> str | None:
    columns = list(columns)
    by_exact = {column.lower(): column for column in columns}
    by_clean = {_clean_key(column): column for column in columns}
    for alias in COLUMN_ALIASES[canonical]:
        if alias.lower() in by_exact:
            return by_exact[alias.lower()]
        cleaned = _clean_key(alias)
        if cleaned in by_clean:
            return by_clean[cleaned]
    return None


def _string(value: Any) -> str:
    if _is_missing(value):
        return ""
    return str(value).strip()


def normalize_pdb_id(value: Any) -> str:
    """Return an uppercase four-character-style PDB identifier when possible."""

    return _string(value).upper()


def normalize_chain_id(value: Any) -> str:
    return _string(value)


def normalize_antigen_chains(value: Any) -> str:
    """Normalize antigen chain separators to semicolons and remove duplicates."""

    if _is_missing(value):
        return ""
    parts = re.split(r"[;,\|/\s]+", str(value).strip())
    seen: set[str] = set()
    normalized: list[str] = []
    for part in parts:
        part = part.strip()
        if not part or part in seen:
            continue
        seen.add(part)
        normalized.append(part)
    return ";".join(normalized)


def parse_date(value: Any) -> date | None:
    if _is_missing(value):
        return None
    import pandas as pd

    parsed = pd.to_datetime(value, errors="coerce")
    if pd.isna(parsed):
        return None
    return parsed.date()


def assign_split(
    release_date: Any,
    train_cutoff: date | str = date(2022, 1, 1),
    val_cutoff: date | str = date(2022, 7, 1),
    test_cutoff: date | str = date(2023, 1, 1),
    strict_protenix_cutoff: bool = False,
) -> tuple[str, str]:
    """Assign a release date to the canonical nanobody split names."""

    policy = SplitPolicy(
        train_cutoff=parse_date(train_cutoff) or train_cutoff,
        val_cutoff=parse_date(val_cutoff) or val_cutoff,
        test_cutoff=parse_date(test_cutoff) or test_cutoff,
        strict_protenix_cutoff=strict_protenix_cutoff,
    )
    parsed = parse_date(release_date)
    if parsed is None:
        return "unknown", "unknown"
    if parsed < policy.effective_train_cutoff:
        return "train", "train_pre2022"
    if strict_protenix_cutoff and parsed < policy.train_cutoff:
        return "gap", "excluded_2021q4"
    if parsed < policy.val_cutoff:
        return "val", "val_2022h1"
    if parsed < policy.test_cutoff:
        return "test", "test_2022h2"
    return "holdout", "holdout_post2022"


def resolve_structure_path(
    pdb_id: str,
    explicit_path: Any,
    structure_dir: str | os.PathLike[str] | None,
    metadata_dir: str | os.PathLike[str] | None = None,
) -> str:
    """Resolve metadata-provided or conventional local structure paths."""

    if not _is_missing(explicit_path):
        candidate = Path(str(explicit_path)).expanduser()
        if not candidate.is_absolute() and metadata_dir is not None:
            candidate = Path(metadata_dir) / candidate
        return str(candidate)
    if structure_dir is None:
        return ""

    root = Path(structure_dir).expanduser()
    names = [
        pdb_id,
        pdb_id.lower(),
        pdb_id.upper(),
    ]
    suffixes = [".cif", ".cif.gz", ".mmcif", ".mmcif.gz", ".pdb", ".pdb.gz"]
    for name in names:
        for suffix in suffixes:
            candidate = root / f"{name}{suffix}"
            if candidate.exists():
                return str(candidate)
    return str(root / f"{pdb_id.lower()}.cif")


def _get(row: Mapping[str, Any], column: str | None) -> Any:
    if column is None:
        return ""
    return row.get(column, "")


def _coerce_optional_int(value: Any) -> str:
    if _is_missing(value):
        return ""
    try:
        return str(int(float(str(value).strip())))
    except ValueError:
        return str(value).strip()


def normalize_metadata_table(
    metadata: Any,
    structure_dir: str | os.PathLike[str] | None = None,
    metadata_dir: str | os.PathLike[str] | None = None,
    train_cutoff: date | str = date(2022, 1, 1),
    val_cutoff: date | str = date(2022, 7, 1),
    test_cutoff: date | str = date(2023, 1, 1),
    strict_protenix_cutoff: bool = False,
) -> Any:
    """Normalize SAbDab-like metadata into task-specific manifest rows."""

    columns = metadata.columns.tolist()
    mapped = {key: _find_column(columns, key) for key in COLUMN_ALIASES}
    required = ["pdb_id", "release_date", "nanobody_chain"]
    missing = [key for key in required if mapped[key] is None]
    if missing:
        raise ValueError(
            "Metadata is missing required columns: "
            + ", ".join(missing)
            + ". Accepted aliases are defined in protenix.nanobody.manifest.COLUMN_ALIASES."
        )

    records: list[dict[str, Any]] = []
    for raw in metadata.to_dict(orient="records"):
        pdb_id = normalize_pdb_id(_get(raw, mapped["pdb_id"]))
        release = parse_date(_get(raw, mapped["release_date"]))
        release_text = release.isoformat() if release is not None else ""
        split, split_name = assign_split(
            release,
            train_cutoff=train_cutoff,
            val_cutoff=val_cutoff,
            test_cutoff=test_cutoff,
            strict_protenix_cutoff=strict_protenix_cutoff,
        )
        nanobody_chain = normalize_chain_id(_get(raw, mapped["nanobody_chain"]))
        antigen_chains = normalize_antigen_chains(_get(raw, mapped["antigen_chains"]))
        structure_path = resolve_structure_path(
            pdb_id=pdb_id,
            explicit_path=_get(raw, mapped["structure_path"]),
            structure_dir=structure_dir,
            metadata_dir=metadata_dir,
        )
        cdrh3_start = _coerce_optional_int(_get(raw, mapped["cdrh3_start"]))
        cdrh3_end = _coerce_optional_int(_get(raw, mapped["cdrh3_end"]))
        cdrh3_source = "metadata" if cdrh3_start and cdrh3_end else "missing"
        base = {
            "pdb_id": pdb_id,
            "release_date": release_text,
            "split": split,
            "split_name": split_name,
            "nanobody_chain": nanobody_chain,
            "structure_path": structure_path,
            "prepared_mmcif_path": "",
            "cdrh3_start": cdrh3_start,
            "cdrh3_end": cdrh3_end,
            "cdrh3_source": cdrh3_source,
            "num_tokens": "",
            "notes": "",
        }
        single = dict(base)
        single["task_type"] = "single"
        single["antigen_chains"] = ""
        records.append(single)
        if antigen_chains:
            complex_row = dict(base)
            complex_row["task_type"] = "complex"
            complex_row["antigen_chains"] = antigen_chains
            records.append(complex_row)

    import pandas as pd

    return pd.DataFrame(records, columns=CANONICAL_MANIFEST_COLUMNS)


def read_metadata_csv(path: str | os.PathLike[str]) -> Any:
    import pandas as pd

    return pd.read_csv(path)


def normalize_metadata_csv(
    metadata_path: str | os.PathLike[str],
    structure_dir: str | os.PathLike[str] | None = None,
    **kwargs: Any,
) -> Any:
    path = Path(metadata_path)
    return normalize_metadata_table(
        read_metadata_csv(path),
        structure_dir=structure_dir,
        metadata_dir=path.parent,
        **kwargs,
    )


def split_manifest_filename(task_type: str, split_name: str) -> str:
    return f"nanobody_{task_type}_{split_name}.csv"


def write_split_manifests(manifest: Any, out_root: str | os.PathLike[str]) -> list[Path]:
    """Write canonical task/split manifests and return the paths written."""

    manifest_dir = Path(out_root) / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for (task_type, split_name), subset in manifest.groupby(["task_type", "split_name"]):
        if subset.empty or split_name in {"unknown", "excluded_2021q4", "holdout_post2022"}:
            continue
        path = manifest_dir / split_manifest_filename(task_type, split_name)
        subset.to_csv(path, index=False)
        written.append(path)

    cdr_cols = ["pdb_id", "nanobody_chain", "cdrh3_start", "cdrh3_end", "cdrh3_source"]
    cdr_path = manifest_dir / "cdrh3_annotations.csv"
    manifest[cdr_cols].drop_duplicates().to_csv(cdr_path, index=False)
    written.append(cdr_path)
    return written


def ensure_index_columns(df: Any) -> Any:
    """Make generated Protenix index rows explicit and stable across parser variants."""

    out = df.copy()
    for column in PROTENIX_INDEX_COLUMNS:
        if column not in out.columns:
            out[column] = ""
    return out[PROTENIX_INDEX_COLUMNS + [c for c in out.columns if c not in PROTENIX_INDEX_COLUMNS]]


def filter_indices_for_manifest_row(indices: Any, row: Mapping[str, Any]) -> Any:
    """Keep Protenix parser rows relevant to a single nanobody task row."""

    if indices.empty:
        return ensure_index_columns(indices)
    df = ensure_index_columns(indices)
    task_type = str(row["task_type"])
    nanobody_chain = str(row["nanobody_chain"])
    if task_type == "single":
        mask = (df["type"].astype(str) == "chain") & (
            df["chain_1_id"].astype(str) == nanobody_chain
        )
        return df[mask].copy()

    antigens = set(str(row.get("antigen_chains", "")).split(";")) - {""}
    chain_1 = df["chain_1_id"].astype(str)
    chain_2 = df["chain_2_id"].astype(str)
    mask = (df["type"].astype(str) == "interface") & (
        ((chain_1 == nanobody_chain) & chain_2.isin(antigens))
        | ((chain_2 == nanobody_chain) & chain_1.isin(antigens))
    )
    return df[mask].copy()
