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

from collections.abc import Sequence
from dataclasses import dataclass

from protenix.data.constants import MSA_PROTEIN_SEQ_TO_ID
from protenix.data.msa.msa_utils import extract_species_id
from protenix.data.tools.common import parse_fasta


@dataclass(frozen=True)
class MonomerMsaInput:
    """Logical split of one monomer A3M into paired and unpaired MSA streams."""

    paired_a3m: str
    unpaired_a3m: str
    total_rows: int
    pairable_rows: int
    unpairable_rows: int
    invalid_rows: int


def _aligned_len(seq: str) -> int:
    """Count A3M aligned columns using the same protein alphabet as MSA encoding."""
    return sum(1 for c in seq if c in MSA_PROTEIN_SEQ_TO_ID)


def _to_a3m(rows: Sequence[tuple[str, str]]) -> str:
    return "".join(f">{description}\n{sequence}\n" for description, sequence in rows)


def split_monomer_msa(
    query_sequence: str,
    a3m: str,
    source_name: str = "<monomerMsaPath>",
) -> MonomerMsaInput:
    """
    Split one monomer A3M into Protenix paired and unpaired MSA streams.

    All valid non-query rows are kept as unpaired MSA signal. Rows with an
    extractable species/taxonomy identifier are also exposed to the existing
    taxonomy-based pairing path.
    """
    sequences, descriptions = parse_fasta(a3m)
    if not sequences:
        raise ValueError(f"{source_name} does not contain any A3M/FASTA rows.")

    query_aligned_len = _aligned_len(query_sequence)
    first_aligned_len = _aligned_len(sequences[0])
    if first_aligned_len != query_aligned_len:
        raise ValueError(
            f"{source_name} query row aligned length ({first_aligned_len}) does not "
            f"match protein sequence length ({query_aligned_len})."
        )

    query_row = ("query", query_sequence)
    paired_rows = [query_row]
    unpaired_rows = [query_row]
    pairable_rows = 0
    unpairable_rows = 0
    invalid_rows = 0

    for seq, desc in zip(sequences[1:], descriptions[1:]):
        if _aligned_len(seq) != query_aligned_len:
            invalid_rows += 1
            continue

        unpaired_rows.append((desc, seq))
        species_id = extract_species_id(desc)
        if species_id:
            paired_rows.append((desc, seq))
            pairable_rows += 1
        else:
            unpairable_rows += 1

    return MonomerMsaInput(
        paired_a3m=_to_a3m(paired_rows),
        unpaired_a3m=_to_a3m(unpaired_rows),
        total_rows=len(sequences),
        pairable_rows=pairable_rows,
        unpairable_rows=unpairable_rows,
        invalid_rows=invalid_rows,
    )
