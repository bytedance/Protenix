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
import re
import string
from collections.abc import Mapping, MutableMapping, Sequence
from typing import Any, List, Optional, Set, Tuple

import numpy as np

from protenix.data.constants import (
    DNA_CHAIN,
    MSA_DNA_SEQ_TO_ID,
    MSA_PROTEIN_SEQ_TO_ID,
    MSA_RNA_SEQ_TO_ID,
    PROTEIN_CHAIN,
    RNA_CHAIN,
    STANDARD_POLYMER_CHAIN_TYPES,
    STD_RESIDUES_WITH_GAP,
)
from protenix.data.tools.common import parse_fasta
from protenix.utils.logger import get_logger

logger = get_logger(__name__)

# Type Aliases
FeatureDict = Mapping[str, np.ndarray]
MutableFeatureDict = MutableMapping[str, np.ndarray]

# Regex for species identification
_UNIPROT_REGEX = re.compile(
    r"(?:tr|sp)\|[A-Z0-9]{6,10}(?:_\d+)?\|(?:[A-Z0-9]{1,10}_)(?P<SpeciesId>[A-Z0-9]{1,5})"
)
_UNIREF_REGEX = re.compile(r"^UniRef100_[^_]+_([^_/]+)")

MSA_GAP_IDX = STD_RESIDUES_WITH_GAP.get("-")
NUM_SEQ_NUM_RES_MSA_FEATURES = ("msa", "msa_mask", "deletion_matrix")
NUM_SEQ_MSA_FEATURES = ("msa_species_identifiers",)
MSA_PAD_VALUES = {"msa": MSA_GAP_IDX, "msa_mask": 1, "deletion_matrix": 0}


class MSACore:
    """Basic MSA parsing and numerical conversion operations."""

    @staticmethod
    def parse_fasta(fasta_str: str) -> Tuple[List[str], List[str]]:
        """
        Parses a FASTA/A3M format string into sequences and descriptions.

        Args:
            fasta_str: The input string in FASTA or A3M format.

        Returns:
            A tuple of (list of sequences, list of descriptions).
        """
        return parse_fasta(fasta_str)

    @staticmethod
    def sequences_to_array(
        sequences: Sequence[str],
        chain_type: str,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Converts sequence strings into numerical arrays and deletion matrices.

        Args:
            sequences: A list of aligned sequence strings.
            chain_type: Type of the chain (PROTEIN_CHAIN, DNA_CHAIN, RNA_CHAIN).

        Returns:
            A tuple of (msa_array, deletion_array).

        Raises:
            ValueError: If the chain type is unsupported.
        """
        char_maps = {
            PROTEIN_CHAIN: MSA_PROTEIN_SEQ_TO_ID,
            DNA_CHAIN: MSA_DNA_SEQ_TO_ID,
            RNA_CHAIN: MSA_RNA_SEQ_TO_ID,
        }
        if chain_type not in char_maps:
            raise ValueError(f"Unsupported chain type: {chain_type}")

        cmap = char_maps[chain_type]
        if not sequences:
            return np.zeros((0, 0), dtype=np.int32), np.zeros((0, 0), dtype=np.int32)

        rows, cols = len(sequences), sum(1 for c in sequences[0] if c in cmap)
        msa_arr = np.zeros((rows, cols), dtype=np.int32)
        del_arr = np.zeros((rows, cols), dtype=np.int32)

        for i, seq in enumerate(sequences):
            d_count, c_idx = 0, 0
            for char in seq:
                val = cmap.get(char, -1)
                if val == -1:
                    d_count += 1
                else:
                    if c_idx < cols:
                        msa_arr[i, c_idx] = val
                        del_arr[i, c_idx] = d_count
                    d_count, c_idx = 0, c_idx + 1
        return msa_arr, del_arr


class RawMsa:
    """
    Represents a single MSA and handles deduplication and featurization logic.

    Args:
        query: The query sequence string.
        chain_type: Type of the chain.
        sequences: List of sequence strings in the alignment.
        descriptions: List of description strings for the sequences.
        deduplicate: Whether to deduplicate sequences.
    """

    def __init__(
        self,
        query: str,
        chain_type: str,
        sequences: Sequence[str],
        descriptions: Sequence[str],
        deduplicate: bool = True,
    ) -> None:
        self.query = query
        self.chain_type = chain_type
        if deduplicate:
            self.seqs, self.descs = self._deduplicate_sequences(sequences, descriptions)
        else:
            self.seqs, self.descs = list(sequences), list(descriptions)

        # Make sure the MSA always has at least the query.
        self.seqs = self.seqs or [query]
        self.descs = self.descs or ["Original query"]

        if not self._verify_query():
            raise ValueError(
                f"MSA query/size mismatch for {chain_type}, falling back to query only."
            )

    def _verify_query(self) -> bool:
        """Ensures the first sequence in the MSA is equivalent to the query sequence."""
        if self.chain_type not in STANDARD_POLYMER_CHAIN_TYPES:
            return self.seqs[0] == self.query
        q_arr, _ = MSACore.sequences_to_array([self.query], self.chain_type)
        m_arr, _ = MSACore.sequences_to_array([self.seqs[0]], self.chain_type)
        # Only check the shape in case the query (the first sequence in the MSA)
        # is not explicitly included in the A3M file provided by the user.
        return np.array_equal(q_arr.shape, m_arr.shape)

    @staticmethod
    def _deduplicate_sequences(
        seqs: Sequence[str], descs: Sequence[str]
    ) -> Tuple[List[str], List[str]]:
        """Removes duplicate sequences, ignoring insertions (lowercase characters)."""
        u_seqs, u_descs, seen = [], [], set()
        table = str.maketrans("", "", string.ascii_lowercase)
        for s, d in zip(seqs, descs):
            stripped = s.translate(table)
            if stripped not in seen:
                seen.add(stripped)
                u_seqs.append(s)
                u_descs.append(d)
        return u_seqs, u_descs

    @classmethod
    def from_a3m(
        cls,
        query: str,
        ctype: str,
        a3m: str,
        depth_limit: Optional[int] = None,
        dedup: bool = True,
    ) -> "RawMsa":
        """
        Constructs a RawMsa instance from A3M content.

        Args:
            query: The query sequence.
            ctype: The chain type.
            a3m: The A3M format string.
            depth_limit: Optional maximum number of sequences to keep.
            dedup: Whether to deduplicate sequences.

        Returns:
            An instance of RawMsa.
        """
        s, d = MSACore.parse_fasta(a3m)
        if depth_limit and 0 < depth_limit < len(s):
            s, d = s[:depth_limit], d[:depth_limit]
        return cls(query, ctype, s, d, dedup)

    @classmethod
    def merge(cls, msas: Sequence["RawMsa"], deduplicate: bool = True) -> "RawMsa":
        """
        Merges multiple RawMsa objects into one.

        Args:
            msas: A sequence of RawMsa objects.
            deduplicate: Whether to deduplicate the resulting merged MSA.

        Returns:
            A new RawMsa object representing the merged alignment.

        Raises:
            ValueError: If no MSAs are provided.
        """
        if not msas:
            raise ValueError("No MSAs provided for merging.")
        return cls(
            msas[0].query,
            msas[0].chain_type,
            [s for m in msas for s in m.seqs],
            [d for m in msas for d in m.descs],
            deduplicate,
        )

    @property
    def depth(self) -> int:
        """Returns the number of sequences in the MSA."""
        return len(self.seqs)

    def featurize(self) -> MutableFeatureDict:
        """
        Extracts numerical features from the MSA.

        Returns:
            A mutable dictionary of MSA features.
        """
        msa_arr, del_arr = MSACore.sequences_to_array(self.seqs, self.chain_type)
        return {
            "msa": msa_arr,
            "deletion_matrix": del_arr,
            "msa_species_identifiers": np.array(
                MSAPairingEngine.get_species_ids(self.descs), dtype=object
            ),
            "num_alignments": np.array(len(self.seqs), dtype=np.int32),
        }

    def to_a3m(self) -> str:
        """Serializes the MSA back into A3M string format."""
        return "\n".join([f">{d}\n{s}" for d, s in zip(self.descs, self.seqs)]) + "\n"


class MSAPairingEngine:
    """Manages MSA pairing by species and post-pairing processing."""

    @staticmethod
    def get_species_ids(descriptions: Sequence[str]) -> List[str]:
        """Extracts species identifiers from UniProt or UniRef description lines."""
        ids = []
        for d in descriptions:
            d = d.strip()
            m = _UNIPROT_REGEX.match(d) or _UNIREF_REGEX.match(d)
            if m:
                ids.append(
                    m.group("SpeciesId") if "SpeciesId" in m.groupdict() else m.group(1)
                )
            else:
                ids.append("")
        return ids

    @staticmethod
    def _align_species(
        all_species: Sequence[str],
        chain_species_map: Sequence[Mapping[str, np.ndarray]],
        species_min_hits: Mapping[str, int],
    ) -> np.ndarray:
        """Aligns MSA row indices based on species."""
        species_blocks = []
        for species in all_species:
            chain_row_indices = []
            for species_to_rows in chain_species_map:
                min_msa_size = species_min_hits[species]
                if species not in species_to_rows:
                    row_indices = np.full(min_msa_size, fill_value=-1, dtype=np.int32)
                else:
                    row_indices = species_to_rows[species][:min_msa_size]
                chain_row_indices.append(row_indices)
            species_block = np.stack(chain_row_indices, axis=1)
            species_blocks.append(species_block)
        return np.concatenate(species_blocks, axis=0)

    @staticmethod
    def pair_chains_by_species(
        chains: Sequence[MutableFeatureDict],
        max_paired: int,
        active_chains: Set[str],
        max_per_species: int,
    ) -> List[MutableFeatureDict]:
        """
        Aligns MSA rows across different chains based on species identifiers.

        Args:
            chains: List of chain feature dictionaries.
            max_paired: Maximum number of paired rows.
            active_chains: Set of chain IDs that should be paired.
            max_per_species: Maximum number of rows allowed per species.

        Returns:
            A list of feature dictionaries with paired MSA information.
        """
        chain_species_map = []
        all_species_counts = {}
        species_min_hits = {}

        for c in chains:
            ids = c.get("msa_species_identifiers_all_seq", np.array([]))
            if (
                ids.size == 0
                or (ids.size == 1 and not ids[0])
                or c["chain_id"] not in active_chains
            ):
                chain_species_map.append({})
                continue

            row_indices = np.arange(len(ids))
            sort_idxs = ids.argsort()
            ids = ids[sort_idxs]
            row_indices = row_indices[sort_idxs]

            species, unique_row_indices = np.unique(ids, return_index=True)
            grouped_row_indices = np.split(row_indices, unique_row_indices[1:])
            mapping = dict(zip(species, grouped_row_indices))
            chain_species_map.append(mapping)

            for s in species:
                all_species_counts[s] = all_species_counts.get(s, 0) + 1

            for s, idxs in mapping.items():
                species_min_hits[s] = min(
                    species_min_hits.get(s, max_per_species), len(idxs)
                )

        ranked_species = {}
        for s, count in all_species_counts.items():
            if not s or count <= 1:
                continue
            ranked_species.setdefault(count, []).append(s)

        pair_idxs = [np.zeros((1, len(chains)), dtype=np.int32)]  # Always keep Query
        current_total = 0
        for count in sorted(ranked_species.keys(), reverse=True):
            all_species = ranked_species[count]
            rows = MSAPairingEngine._align_species(
                all_species, chain_species_map, species_min_hits
            )
            rank = np.sum(np.log(np.abs(rows.astype(np.float32)) + 1e-10), axis=1)
            pair_idxs.append(rows[np.argsort(rank)])
            current_total += rows.shape[0]
            if current_total >= max_paired:
                break

        final_idxs = np.concatenate(pair_idxs, axis=0)[:max_paired]

        new_chains = []
        for i, c in enumerate(chains):
            nc = {k: v for k, v in c.items() if "all_seq" not in k}
            sel = final_idxs[:, i]
            for f in ["msa", "deletion_matrix"]:
                src = c[f"{f}_all_seq"]
                padded = np.concatenate(
                    [src, np.full((1, src.shape[1]), MSA_PAD_VALUES[f], src.dtype)],
                    axis=0,
                )
                nc[f"{f}_all_seq"] = padded[sel]
            nc["num_alignments_all_seq"] = np.array(nc["msa_all_seq"].shape[0])
            new_chains.append(nc)
        return new_chains

    @staticmethod
    def cleanup_unpaired_features(
        chains: List[MutableFeatureDict],
    ) -> List[MutableFeatureDict]:
        """Removes sequences from unpaired MSA that are already present in paired MSA."""
        for c in chains:
            paired_bytes = {row.tobytes() for row in c["msa_all_seq"].astype(np.int8)}
            keep = [
                i
                for i, row in enumerate(c["msa"].astype(np.int8))
                if row.tobytes() not in paired_bytes
            ]
            for k in ["msa", "deletion_matrix", "msa_species_identifiers"]:
                c[k] = c[k][keep]
            c["num_alignments"] = np.array(c["msa"].shape[0], dtype=np.int32)
        return chains

    @staticmethod
    def filter_all_gapped_rows(
        chains: List[MutableFeatureDict], active_asyms: Sequence[int]
    ) -> List[MutableFeatureDict]:
        """Removes rows from the paired MSA that consist entirely of gaps."""
        subset = [c["msa_all_seq"] for c in chains if c["asym_id"][0] in active_asyms]
        if not subset:
            return chains
        non_gap_mask = np.any(np.concatenate(subset, axis=1) != MSA_GAP_IDX, axis=1)
        for c in chains:
            for k in ["msa_all_seq", "deletion_matrix_all_seq"]:
                if k in c:
                    c[k] = c[k][non_gap_mask]
            c["num_alignments_all_seq"] = np.array(np.sum(non_gap_mask))
        return chains

    @staticmethod
    def merge_chain_features(chains: List[FeatureDict], key: str) -> np.ndarray:
        """Merges MSA features across different chains, padding unpaired sequences where necessary."""
        if "_all_seq" in key:
            return np.concatenate([c[key] for c in chains], axis=1)
        max_d = max(c[key].shape[0] for c in chains)
        pads = [
            np.pad(
                c[key],
                ((0, max_d - c[key].shape[0]), (0, 0)),
                constant_values=MSA_PAD_VALUES.get(key, 0),
            )
            for c in chains
        ]
        return np.concatenate(pads, axis=1)


def map_to_standard(
    asym_ids: np.ndarray, res_ids: np.ndarray, meta: Mapping[int, Mapping[str, Any]]
) -> np.ndarray:
    """
    Maps residue indices to a standardized MSA coordinate system.

    Args:
        asym_ids: Array of asymmetric IDs.
        res_ids: Array of residue IDs.
        meta: Metadata dictionary containing chain info and sequences.

    Returns:
        Array of standardized indices.
    """
    uids = [f"{a}-{b}" for a, b in zip(asym_ids, res_ids)]
    std_uids = []
    for aid, info in meta.items():
        std_uids.extend([f"{aid}-{x}" for x in range(1, len(info["sequence"]) + 1)])

    lookup = {uid: i for i, uid in enumerate(std_uids)}
    return np.array(
        [lookup.get(u, lookup.get(f"{u.split('-')[0]}-1", -1)) for u in uids],
        dtype=np.int32,
    )
