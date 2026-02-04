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
# Copyright 2021 AlQuraishi Laboratory

import collections
import dataclasses
import functools
import io
import re
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from Bio import PDB
from Bio.Data import PDBData

from protenix.data.constants import (
    DNA_CHAIN,
    PROTEIN_CHAIN,
    RNA_CHAIN,
    TEMPLATE_DNA_SEQ_TO_ID,
    TEMPLATE_PROTEIN_SEQ_TO_ID,
    TEMPLATE_RNA_SEQ_TO_ID,
)
from protenix.data.tools.common import parse_fasta
from protenix.utils.logger import get_logger

logger = get_logger(__name__)

# Type aliases for clarity.
MmCIFDict = Mapping[str, Sequence[str]]
PdbHeader = Mapping[str, Any]
ChainId = str
SeqRes = str


# --- Exceptions ---
"""Adapted from openfold.data"""


class TemplateError(Exception):
    """Base class for exceptions in this module."""


class ParsingError(TemplateError):
    """Raised when parsing mmCIF file fails or contains invalid data."""


class NoChainsError(ParsingError):
    """Raised when template mmCIF doesn't have any chains."""


class NoAtomDataInTemplateError(ParsingError):
    """Raised when template mmCIF doesn't contain atom positions."""


class TemplateAtomMaskAllZerosError(ParsingError):
    """Raised when template mmCIF had all atom positions masked."""


class MultipleChainsError(ParsingError):
    """Raised when multiple chains are found for a given ID where only one was expected."""


class AlignmentError(TemplateError):
    """Raised when alignment between query and template fails."""


class SequenceNotInTemplateError(AlignmentError):
    """Raised when template mmCIF doesn't contain the expected sequence."""


class QueryToTemplateAlignError(AlignmentError):
    """Raised when the query cannot be aligned to the template."""


class CaDistanceError(AlignmentError):
    """Raised when a CA atom distance exceeds a threshold."""


class PrefilterError(TemplateError):
    """Base class for template prefilter exceptions."""


class DateError(PrefilterError):
    """Raised when the hit date is after the max allowed date."""


class AlignRatioError(PrefilterError):
    """Raised when the hit align ratio to the query is too small."""


class DuplicateError(PrefilterError):
    """Raised when the hit is an exact subsequence of the query."""


class LengthError(PrefilterError):
    """Raised when the hit is too short."""


# --- Data Classes ---


@dataclasses.dataclass(frozen=True)
class ResiduePosition:
    """Represents a residue position in a structure."""

    chain_id: str
    residue_number: int
    insertion_code: str


@dataclasses.dataclass(frozen=True)
class ResidueAtPosition:
    """Represents a residue at a specific position, possibly missing."""

    position: Optional[ResiduePosition]
    name: str
    is_missing: bool
    hetflag: str


@dataclasses.dataclass(frozen=True)
class MmcifObject:
    """Representation of a parsed mmCIF file."""

    file_id: str
    header: PdbHeader
    structure: PDB.Structure.Structure
    chain_to_seqres: Mapping[ChainId, SeqRes]
    seqres_to_structure: Mapping[ChainId, Mapping[int, ResidueAtPosition]]
    raw_string: Any


@dataclasses.dataclass(frozen=True)
class ParsingResult:
    """Result of parsing an mmCIF file."""

    mmcif_object: Optional[MmcifObject]
    errors: Mapping[Tuple[str, str], Any]


@dataclasses.dataclass(frozen=True)
class Monomer:
    """Represents a monomer in a polymer chain."""

    id: str
    num: int


@dataclasses.dataclass(frozen=True)
class AtomSite:
    """Represents an atom site from mmCIF _atom_site loop."""

    residue_name: str
    author_chain_id: str
    mmcif_chain_id: str
    author_seq_num: str
    mmcif_seq_num: int
    insertion_code: str
    hetatm_atom: str
    model_num: str


@dataclasses.dataclass(frozen=True)
class PrefilterResult:
    """Result of pre-filtering a template hit."""

    valid: bool
    error: Optional[str] = None
    warning: Optional[str] = None


@dataclasses.dataclass(frozen=True)
class SingleHitResult:
    """Result of processing a single template hit."""

    hit: Optional["TemplateHit"]
    features: Optional[Mapping[str, Any]]
    error: Optional[str]
    warning: Optional[str]


@dataclasses.dataclass(frozen=True)
class TemplateSearchResult:
    """Result of a template search across multiple hits."""

    features: Sequence[Mapping[str, Any]]
    hits: Sequence["TemplateHit"]
    errors: Sequence[str]
    warnings: Sequence[str]


# --- mmCIF Parsing Helpers ---


def encode_template_restype(chain_type: str, sequence: str) -> List[int]:
    """Encodes a sequence of residues into integer IDs based on chain type."""
    if chain_type == PROTEIN_CHAIN:
        return [
            TEMPLATE_PROTEIN_SEQ_TO_ID.get(r, TEMPLATE_PROTEIN_SEQ_TO_ID["X"])
            for r in sequence
        ]
    if chain_type == RNA_CHAIN:
        return [
            TEMPLATE_RNA_SEQ_TO_ID.get(r, TEMPLATE_RNA_SEQ_TO_ID["N"]) for r in sequence
        ]
    if chain_type == DNA_CHAIN:
        return [
            TEMPLATE_DNA_SEQ_TO_ID.get(r, TEMPLATE_DNA_SEQ_TO_ID["N"]) for r in sequence
        ]
    raise NotImplementedError(f"Unsupported chain type: {chain_type}")


def get_pdb_id_and_chain(hit: "TemplateHit") -> Tuple[str, str]:
    """Extracts PDB ID and chain from hit name."""
    match = re.match(r"([a-zA-Z\d]{4})_([a-zA-Z0-9.]+)", hit.name)
    if not match:
        raise ValueError(f"Invalid hit name format: {hit.name}")
    return match.group(1).lower(), match.group(2)


class TemplateParser:
    """Class to parse mmCIF files into MmcifObject."""

    @staticmethod
    def _mmcif_loop_to_list(
        prefix: str, parsed_info: MmCIFDict
    ) -> List[Dict[str, str]]:
        """Extracts an mmCIF loop into a list of dictionaries."""
        cols = []
        data = []
        for key, value in parsed_info.items():
            if key.startswith(prefix):
                cols.append(key)
                data.append(value)

        if not data:
            return []

        assert all(
            len(xs) == len(data[0]) for xs in data
        ), f"mmCIF error: Not all loops are the same length for prefix {prefix}"

        return [dict(zip(cols, xs)) for xs in zip(*data)]

    @staticmethod
    def _mmcif_loop_to_dict(
        prefix: str, index_key: str, parsed_info: MmCIFDict
    ) -> Dict[str, Dict[str, str]]:
        """Extracts an mmCIF loop into a dictionary keyed by a specific column."""
        entries = TemplateParser._mmcif_loop_to_list(prefix, parsed_info)
        return {entry[index_key]: entry for entry in entries}

    @staticmethod
    def _get_release_date(parsed_info: MmCIFDict) -> str:
        """Returns the oldest revision date from mmCIF data."""
        revision_dates = parsed_info.get("_pdbx_audit_revision_history.revision_date")
        if revision_dates:
            return min(revision_dates)
        return ""

    @staticmethod
    def _get_header(parsed_info: MmCIFDict) -> PdbHeader:
        """Extracts basic header information (method, release date, resolution)."""
        header = {}
        experiments = TemplateParser._mmcif_loop_to_list("_exptl.", parsed_info)
        header["structure_method"] = ",".join(
            [exp.get("_exptl.method", "").lower() for exp in experiments]
        )

        release_date = TemplateParser._get_release_date(parsed_info)
        if release_date:
            header["release_date"] = release_date
        else:
            logger.warning(
                "Could not determine release_date for %s", parsed_info.get("_entry.id")
            )

        header["resolution"] = 0.00
        for res_key in (
            "_refine.ls_d_res_high",
            "_em_3d_reconstruction.resolution",
            "_reflns.d_resolution_high",
        ):
            if res_key in parsed_info:
                try:
                    header["resolution"] = float(parsed_info[res_key][0])
                    break
                except (ValueError, IndexError):
                    continue
        return header

    @staticmethod
    def _get_first_model(structure: PDB.Structure.Structure) -> PDB.Model.Model:
        """Returns the first model in a Biopython structure."""
        return next(structure.get_models())

    @staticmethod
    def _get_atom_site_list(parsed_info: MmCIFDict) -> List[AtomSite]:
        """Returns a list of AtomSite objects from mmCIF data."""
        return [
            AtomSite(*site)
            for site in zip(
                parsed_info["_atom_site.label_comp_id"],
                parsed_info["_atom_site.auth_asym_id"],
                parsed_info["_atom_site.label_asym_id"],
                parsed_info["_atom_site.auth_seq_id"],
                parsed_info["_atom_site.label_seq_id"],
                parsed_info["_atom_site.pdbx_PDB_ins_code"],
                parsed_info["_atom_site.group_PDB"],
                parsed_info["_atom_site.pdbx_PDB_model_num"],
            )
        ]

    @staticmethod
    def _get_protein_chains(parsed_info: MmCIFDict) -> Dict[ChainId, Sequence[Monomer]]:
        """Extracts valid protein chains and their sequences."""
        entity_poly_seqs = TemplateParser._mmcif_loop_to_list(
            "_entity_poly_seq.", parsed_info
        )
        polymers = collections.defaultdict(list)
        for entry in entity_poly_seqs:
            polymers[entry["_entity_poly_seq.entity_id"]].append(
                Monomer(
                    id=entry["_entity_poly_seq.mon_id"],
                    num=int(entry["_entity_poly_seq.num"]),
                )
            )

        chem_comps = TemplateParser._mmcif_loop_to_dict(
            "_chem_comp.", "_chem_comp.id", parsed_info
        )
        struct_asyms = TemplateParser._mmcif_loop_to_list("_struct_asym.", parsed_info)

        entity_to_mmcif_chains = collections.defaultdict(list)
        for asym in struct_asyms:
            entity_to_mmcif_chains[asym["_struct_asym.entity_id"]].append(
                asym["_struct_asym.id"]
            )

        valid_chains = {}
        for entity_id, seq_info in polymers.items():
            # Check if any component is peptide-like.
            is_protein = any(
                "peptide" in chem_comps.get(m.id, {}).get("_chem_comp.type", "")
                for m in seq_info
            )
            if is_protein:
                for chain_id in entity_to_mmcif_chains[entity_id]:
                    valid_chains[chain_id] = seq_info
        return valid_chains

    @staticmethod
    def _is_set(data: str) -> bool:
        """Checks if mmCIF data is set (not '.' or '?')."""
        return data not in (".", "?")

    @staticmethod
    @functools.lru_cache(maxsize=16)
    def parse(
        *,
        file_id: str,
        mmcif_string: str,
        auth_chain_id: Optional[str] = None,
        catch_all_errors: bool = True,
    ) -> ParsingResult:
        """Parses an mmCIF string into an MmcifObject."""
        errors = {}
        try:
            parser = PDB.MMCIFParser(QUIET=True)
            structure = parser.get_structure("", io.StringIO(mmcif_string))
            first_model = TemplateParser._get_first_model(structure)
            parsed_info = parser._mmcif_dict  # pylint: disable=protected-access

            # Ensure all values are lists.
            for k, v in parsed_info.items():
                if not isinstance(v, list):
                    parsed_info[k] = [v]

            header = TemplateParser._get_header(parsed_info)
            valid_chains = TemplateParser._get_protein_chains(parsed_info)
            if not valid_chains:
                return ParsingResult(None, {(file_id, ""): "No protein chains found."})

            seq_start_num = {
                cid: min(m.num for m in seq) for cid, seq in valid_chains.items()
            }
            mmcif_to_author_chain_id = {}
            seq_to_structure_mappings = collections.defaultdict(dict)

            for atom in TemplateParser._get_atom_site_list(parsed_info):
                if atom.model_num != "1":
                    continue
                if auth_chain_id is not None and atom.author_chain_id != auth_chain_id:
                    continue

                mmcif_to_author_chain_id[atom.mmcif_chain_id] = atom.author_chain_id

                if atom.mmcif_chain_id in valid_chains:
                    hetflag = " "
                    if atom.hetatm_atom == "HETATM":
                        hetflag = (
                            "W"
                            if atom.residue_name in ("HOH", "WAT")
                            else f"H_{atom.residue_name}"
                        )

                    ins_code = (
                        atom.insertion_code
                        if TemplateParser._is_set(atom.insertion_code)
                        else " "
                    )
                    pos = ResiduePosition(
                        chain_id=atom.author_chain_id,
                        residue_number=int(atom.author_seq_num),
                        insertion_code=ins_code,
                    )
                    seq_idx = (
                        int(atom.mmcif_seq_num) - seq_start_num[atom.mmcif_chain_id]
                    )
                    seq_to_structure_mappings[atom.author_chain_id][
                        seq_idx
                    ] = ResidueAtPosition(
                        position=pos,
                        name=atom.residue_name,
                        is_missing=False,
                        hetflag=hetflag,
                    )

            # Fill in missing residues.
            for cid, seq_info in valid_chains.items():
                if cid not in mmcif_to_author_chain_id:
                    continue
                auth_cid = mmcif_to_author_chain_id[cid]
                mapping = seq_to_structure_mappings[auth_cid]
                for idx, monomer in enumerate(seq_info):
                    if idx not in mapping:
                        mapping[idx] = ResidueAtPosition(
                            position=None, name=monomer.id, is_missing=True, hetflag=" "
                        )

            # Convert 3-letter to 1-letter sequences.
            auth_chain_to_seq = {}
            for cid, seq_info in valid_chains.items():
                if cid not in mmcif_to_author_chain_id:
                    continue
                auth_cid = mmcif_to_author_chain_id[cid]
                seq = "".join(
                    [
                        PDBData.protein_letters_3to1.get(monomer.id, "X")
                        for monomer in seq_info
                    ]
                )
                # Ensure it's 1-letter.
                seq = "".join([c if len(c) == 1 else "X" for c in seq])
                auth_chain_to_seq[auth_cid] = seq

            mmcif_obj = MmcifObject(
                file_id=file_id,
                header=header,
                structure=first_model,
                chain_to_seqres=auth_chain_to_seq,
                seqres_to_structure=seq_to_structure_mappings,
                raw_string=parsed_info,
            )
            return ParsingResult(mmcif_object=mmcif_obj, errors=errors)
        except Exception as e:
            errors[(file_id, "")] = e
            if not catch_all_errors:
                raise
            return ParsingResult(None, errors=errors)


@dataclasses.dataclass(frozen=True)
class TemplateHit:
    """Represents a template hit from a search tool (e.g., HHSearch, hmmsearch)."""

    index: int
    name: str
    aligned_cols: int
    sum_probs: Optional[float]
    query: str
    hit_sequence: str
    indices_query: List[int]
    indices_hit: List[int]

    @functools.cached_property
    def query_to_hit_mapping(self) -> Mapping[int, int]:
        """Maps 0-based query indices to 0-based hit indices."""
        mapping = {}
        for q_idx, h_idx in zip(self.indices_query, self.indices_hit):
            if (q_idx != -1) and (h_idx != -1):
                mapping[q_idx] = h_idx
        return mapping


@dataclasses.dataclass(frozen=True)
class HitMetadata:
    """Metadata parsed from an hmmsearch A3M description line."""

    pdb_id: str
    chain: str
    start: int
    end: int
    length: int
    text: str


class HHRParser:
    """Class to parse HHR files from HHSearch."""

    @staticmethod
    def parse(hhr_string: str) -> List[TemplateHit]:
        """
        Parses an entire HHR file content.

        Args:
            hhr_string: The content of the HHR file.

        Returns:
            A list of TemplateHit objects.
        """
        lines = hhr_string.splitlines()
        block_starts = [i for i, line in enumerate(lines) if line.startswith("No ")]
        hits = []
        if block_starts:
            block_starts.append(len(lines))
            for i in range(len(block_starts) - 1):
                hits.append(
                    HHRParser._parse_hit(lines[block_starts[i] : block_starts[i + 1]])
                )
        return hits

    @staticmethod
    def _parse_hit(lines: Sequence[str]) -> TemplateHit:
        """Parses a single hit block from an HHR file."""
        hit_num = int(lines[0].split()[-1])
        hit_name = lines[1][1:].strip()
        summary = lines[2]
        match = re.search(r"Aligned_cols=(\d+).*Sum_probs=([0-9.]+)", summary)
        if not match:
            raise RuntimeError(f"Could not parse HHR summary: {summary}")
        cols, probs = int(match.group(1)), float(match.group(2))

        query, hit_seq = "", ""
        idx_q, idx_h = [], []
        for line in lines[3:]:
            if line.startswith("Q ") and not any(
                line.startswith(x) for x in ("Q ss_", "Q Consensus")
            ):
                match = re.search(r"\s+(\d+)\s+([A-Z-]+)\s+(\d+)", line[17:])
                if match:
                    start = int(match.group(1)) - 1
                    seq = match.group(2)
                    query += seq
                    HHRParser._update_residue_indices(seq, start, idx_q)
            elif line.startswith("T ") and not any(
                line.startswith(x) for x in ("T ss_", "T Consensus")
            ):
                match = re.search(r"\s+(\d+)\s+([A-Z-]+)", line[17:])
                if match:
                    start = int(match.group(1)) - 1
                    seq = match.group(2)
                    hit_seq += seq
                    HHRParser._update_residue_indices(seq, start, idx_h)

        return TemplateHit(hit_num, hit_name, cols, probs, query, hit_seq, idx_q, idx_h)

    @staticmethod
    def _update_residue_indices(seq: str, start: int, indices: List[int]):
        """Updates the list of residue indices for a sequence segment."""
        curr = start
        for char in seq:
            if char == "-":
                indices.append(-1)
            else:
                indices.append(curr)
                curr += 1


class HmmsearchA3MParser:
    """Class to parse A3M files from hmmsearch."""

    @staticmethod
    def parse(
        query_seq: str, a3m_str: str, skip_first: bool = True
    ) -> List[TemplateHit]:
        """
        Parses an A3M string from hmmsearch.

        Args:
            query_seq: The query sequence.
            a3m_str: The content of the A3M file.
            skip_first: Whether to skip the first sequence (usually the query).

        Returns:
            A list of TemplateHit objects.
        """
        seqs, descs = parse_fasta(a3m_str)
        parsed = list(zip(seqs, descs))
        if skip_first:
            parsed = parsed[1:]

        idx_q = HmmsearchA3MParser._get_indices(query_seq, 0)
        hits = []
        for i, (h_seq, h_desc) in enumerate(parsed, start=1):
            if "mol:protein" not in h_desc:
                continue
            meta = HmmsearchA3MParser._parse_description(h_desc)
            cols = sum(1 for r in h_seq if r.isupper() and r != "-")
            idx_h = HmmsearchA3MParser._get_indices(h_seq, meta.start - 1)
            hits.append(
                TemplateHit(
                    i,
                    f"{meta.pdb_id}_{meta.chain}",
                    cols,
                    None,
                    query_seq,
                    h_seq.upper(),
                    idx_q,
                    idx_h,
                )
            )
        return hits

    @staticmethod
    def _get_indices(seq: str, start: int) -> List[int]:
        """Calculates residue indices for a sequence with gaps or insertions."""
        indices = []
        curr = start
        for char in seq:
            if char == "-":
                indices.append(-1)
            elif char.islower():
                curr += 1
            else:
                indices.append(curr)
                curr += 1
        return indices

    @staticmethod
    def _parse_description(desc: str) -> HitMetadata:
        """Parses the description line from an hmmsearch hit."""
        pattern = (
            r"^>?([a-z0-9]+)_(\w+)/([0-9]+)-([0-9]+).*protein length:([0-9]+) *(.*)$"
        )
        match = re.match(pattern, desc.strip())
        if not match:
            raise ValueError(f"Could not parse hmmsearch description: {desc}")
        return HitMetadata(
            match[1], match[2], int(match[3]), int(match[4]), int(match[5]), match[6]
        )
