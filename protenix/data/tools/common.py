# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import contextlib
import shutil
import tempfile
from typing import Any, Iterable, List, Optional, Sequence, Tuple


@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None):
    """Context manager that deletes a temporary directory on exit."""
    tmpdir = tempfile.mkdtemp(dir=base_dir)
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def to_a3m(sequences: Sequence[str]) -> str:
    """Converts sequences to an a3m file."""
    names = [f"sequence {i}" for i in range(1, len(sequences) + 1)]
    a3m = []
    for sequence, name in zip(sequences, names):
        a3m.append(">" + name + "\n")
        a3m.append(sequence + "\n")
    return "".join(a3m)


def parse_fasta(data: str) -> Tuple[List[str], List[str]]:
    """
    Parse FASTA/A3M string into sequences and descriptions.

    Args:
        data: The input FASTA or A3M formatted string.

    Returns:
        A tuple of (sequences, descriptions).
    """
    sequences = []
    descriptions = []
    index = -1
    for line in data.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])
            sequences.append("")
        else:
            if index >= 0:
                sequences[index] += line
    return sequences, descriptions


def parse_kalign_a3m(a3m_string: str) -> List[str]:
    """Parses sequences from an A3M format alignment."""
    sequences, _ = parse_fasta(a3m_string)
    return sequences


def lazy_fasta_parse(data: str) -> Iterable[Tuple[str, str]]:
    """
    Memory-friendly FASTA parsing yielding (sequence, description) pairs.

    Args:
        data: The input FASTA or A3M formatted string.

    Yields:
        Tuples of (sequence, description).
    """
    curr_desc, curr_seq = None, []
    for line in data.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if curr_desc:
                yield "".join(curr_seq), curr_desc
            curr_desc, curr_seq = line[1:].strip(), []
        else:
            curr_seq.append(line)
    if curr_desc:
        yield "".join(curr_seq), curr_desc


def a3m_to_sto_list(a3m_seqs: List[str]) -> List[str]:
    """
    Convert a list of A3M sequences to Stockholm format aligned sequences.

    Args:
        a3m_seqs: List of sequences in A3M format.

    Returns:
        List of aligned sequences in Stockholm format.

    Raises:
        ValueError: If A3M rows have inconsistent lengths.
    """
    if not a3m_seqs:
        return []
    cursors = [list(s) for s in a3m_seqs]
    out = [""] * len(a3m_seqs)
    while any(cursors):
        if any(c and c[0].islower() for c in cursors):
            for i, c in enumerate(cursors):
                if c and c[0].islower():
                    out[i] += c.pop(0).upper()
                else:
                    out[i] += "-"
        else:
            for i, c in enumerate(cursors):
                if not c:
                    raise ValueError("Inconsistent row lengths in A3M")
                out[i] += c.pop(0)
    return out


def align_to_query(seq: str, query: str) -> str:
    """
    Align sequence to a gapless query sequence.

    Args:
        seq: Sequence to align.
        query: Reference query sequence.

    Returns:
        Aligned sequence string.

    Raises:
        ValueError: If lengths do not match.
    """
    if len(seq) != len(query):
        raise ValueError("Length mismatch")
    return "".join(
        s if q != "-" else (s.lower() if s != "-" else "") for s, q in zip(seq, query)
    )


def convert_a3m_to_stockholm(a3m: str, max_seqs: Optional[int] = None) -> str:
    """
    Convert A3M format string to Stockholm format string.

    Args:
        a3m: Input A3M format string.
        max_seqs: Maximum number of sequences to include.

    Returns:
        Stockholm format string.
    """
    seqs, descs = parse_fasta(a3m)
    if max_seqs:
        seqs, descs = seqs[:max_seqs], descs[:max_seqs]

    names = [f"{d.split()[0]}_{i}" for i, d in enumerate(descs)]
    gs_lines = [
        f"#=GS {n} DE {d.partition(' ')[2].strip() or '<EMPTY>'}"
        for n, d in zip(names, descs)
    ]
    aligned = a3m_to_sto_list(seqs)

    w = max(len(n) for n in names + ["#=GC RF"])
    msa_lines = [f"{n:<{w}} {s}" for n, s in zip(names, aligned)]
    rf = "".join("." if c == "-" else "x" for c in aligned[0])

    return "\n".join(
        ["# STOCKHOLM 1.0", ""]
        + gs_lines
        + [""]
        + msa_lines
        + [f'{"#=GC RF":<{w}} {rf}', "//"]
    )


def convert_stockholm_to_a3m(
    sto_io: Any,
    max_seqs: Optional[int] = None,
    remove_gaps: bool = True,
    linewidth: Optional[int] = None,
) -> str:
    """
    Convert Stockholm format from an IO object to A3M format string.

    Args:
        sto_io: Iterable of Stockholm lines.
        max_seqs: Maximum sequences to parse.
        remove_gaps: Whether to remove gaps in the first row.
        linewidth: Optional line width for the output sequence.

    Returns:
        A3M format string.
    """
    seqs, descs = {}, {}
    for line in sto_io:
        line = line.strip()
        if not line or line.startswith(("//")):
            continue
        if line.startswith("#=GS"):
            parts = line.split(maxsplit=3)
            if len(parts) >= 3 and parts[2] == "DE":
                descs[parts[1]] = parts[3] if len(parts) == 4 else ""
        elif not line.startswith("#"):
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                name, seq = parts
                if name not in seqs and max_seqs and len(seqs) >= max_seqs:
                    continue
                seqs[name] = seqs.get(name, "") + seq

    if not seqs:
        return ""
    query_seq = next(iter(seqs.values()))
    chunks = []
    for name, s in seqs.items():
        chunks.append(f">{name} {descs.get(name, '')}")
        final_s = (
            align_to_query(s, query_seq).replace(".", "")
            if remove_gaps
            else s.replace(".", "")
        )
        if linewidth:
            chunks.extend(
                final_s[i : i + linewidth] for i in range(0, len(final_s), linewidth)
            )
        else:
            chunks.append(final_s)
    return "\n".join(chunks) + "\n"
