"""CDR-H3 annotation helpers with optional antibody-numbering support."""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class CDRH3Annotation:
    start: int | None
    end: int | None
    source: str
    notes: str = ""

    @property
    def is_available(self) -> bool:
        return self.start is not None and self.end is not None


def _to_int(value: Any) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null", "-"}:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _annotate_with_optional_numbering(sequence: str) -> CDRH3Annotation | None:
    """Try optional ANARCI/anarcii-style packages without hard dependency."""

    try:
        from anarci import anarci  # type: ignore
    except Exception:
        return None

    try:
        numbering, _, _ = anarci([("nanobody", sequence)], scheme="imgt", output=False)
        domains = numbering[0] if numbering else []
        for domain in domains or []:
            numbered_residues = domain[0]
            seq_pos = 0
            cdr_positions: list[int] = []
            for (imgt_pos, _insertion), aa in numbered_residues:
                if aa == "-":
                    continue
                seq_pos += 1
                if 105 <= int(imgt_pos) <= 117:
                    cdr_positions.append(seq_pos)
            if cdr_positions:
                return CDRH3Annotation(
                    min(cdr_positions),
                    max(cdr_positions),
                    "anarci",
                    "IMGT positions 105-117",
                )
    except Exception as exc:
        warnings.warn(
            f"ANARCI CDR-H3 numbering failed ({exc}); falling back to heuristic.",
            RuntimeWarning,
            stacklevel=2,
        )
    return None


def _heuristic_cdrh3(sequence: str) -> CDRH3Annotation:
    """Conservatively find the loop between the conserved Cys and FG/WG motif."""

    seq = re.sub(r"[^A-Za-z]", "", sequence).upper()
    if len(seq) < 80:
        return CDRH3Annotation(None, None, "missing", "sequence too short")

    motif_match = None
    for match in re.finditer(r"[WF]G.G", seq):
        if match.start() >= max(70, len(seq) // 2):
            motif_match = match
            break
    if motif_match is None:
        return CDRH3Annotation(None, None, "missing", "no downstream FG/WG motif")

    cysteine_pos = seq.rfind("C", 0, motif_match.start())
    if cysteine_pos < 0:
        return CDRH3Annotation(None, None, "missing", "no upstream cysteine")

    start = cysteine_pos + 2
    end = motif_match.start()
    if end < start:
        return CDRH3Annotation(None, None, "missing", "empty heuristic loop")
    return CDRH3Annotation(start, end, "heuristic", "1-based sequence positions")


def annotate_cdrh3(
    sequence: str | None = None,
    metadata_start: Any = None,
    metadata_end: Any = None,
    use_optional_numbering: bool = True,
) -> CDRH3Annotation:
    """Annotate CDR-H3 bounds, preferring explicit metadata."""

    start = _to_int(metadata_start)
    end = _to_int(metadata_end)
    if start is not None and end is not None and start <= end:
        return CDRH3Annotation(start, end, "metadata", "metadata-provided bounds")

    if not sequence:
        return CDRH3Annotation(None, None, "missing", "no metadata bounds or sequence")

    if use_optional_numbering:
        numbered = _annotate_with_optional_numbering(sequence)
        if numbered is not None and numbered.is_available:
            return numbered

    return _heuristic_cdrh3(sequence)
