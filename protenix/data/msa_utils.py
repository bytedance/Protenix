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

import logging
import os
import shutil
import string
import subprocess
import time
import uuid
from collections import OrderedDict, defaultdict
from os.path import exists as opexists
from typing import (
    Any,
    Callable,
    Iterable,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Union,
)

import numpy as np

from protenix.data.constants import (
    PRO_STD_RESIDUES,
    PROT_STD_RESIDUES_ONE_TO_THREE,
    RNA_ID_TO_NT,
    RNA_NT_TO_ID,
    RNA_STD_RESIDUES,
)
from protenix.openfold_local.data import parsers
from protenix.openfold_local.data.msa_identifiers import (
    Identifiers,
    _extract_sequence_identifier,
    _parse_sequence_identifier,
)
from protenix.openfold_local.data.msa_pairing import (
    CHAIN_FEATURES,
    MSA_FEATURES,
    MSA_GAP_IDX,
    MSA_PAD_VALUES,
    SEQ_FEATURES,
    block_diag,
    create_paired_features,
    deduplicate_unpaired_sequences,
)
from protenix.openfold_local.np import residue_constants
from protenix.utils.distributed import DIST_WRAPPER

# FeatureDict, make_dummy_msa_obj, convert_monomer_features
# These are modified from openfold: data/data_pipeline
try:
    from protenix.openfold_local.data.tools import jackhmmer
except ImportError:
    print(
        "Failed to import packages for searching MSA; can only run with precomputed MSA"
    )

logger = logging.getLogger(__name__)
FeatureDict = MutableMapping[str, np.ndarray]


SEQ_FEATURES = list(SEQ_FEATURES) + ["profile"]

HHBLITS_INDEX_TO_OUR_INDEX = {
    hhblits_index: (
        PRO_STD_RESIDUES[PROT_STD_RESIDUES_ONE_TO_THREE[hhblits_letter]]
        if hhblits_letter != "-"
        else 31
    )
    for hhblits_index, hhblits_letter in residue_constants.ID_TO_HHBLITS_AA.items()
}
NEW_ORDER_LIST = tuple(
    HHBLITS_INDEX_TO_OUR_INDEX[idx] for idx in range(len(HHBLITS_INDEX_TO_OUR_INDEX))
)

RNA_ID_TO_OUR_INDEX = {
    key: RNA_STD_RESIDUES[value] if value != "-" else 31
    for key, value in RNA_ID_TO_NT.items()
}
RNA_NEW_ORDER_LIST = tuple(
    RNA_ID_TO_OUR_INDEX[idx] for idx in range(len(RNA_ID_TO_OUR_INDEX))
)

RNA_MSA_GAP_IDX = RNA_NT_TO_ID["-"]

RNA_MSA_PAD_VALUES = {
    "msa_all_seq": RNA_MSA_GAP_IDX,
    "msa_mask_all_seq": 1,
    "deletion_matrix_all_seq": 0,
    "deletion_matrix_int_all_seq": 0,
    "msa": RNA_MSA_GAP_IDX,
    "msa_mask": 1,
    "deletion_matrix": 0,
    "deletion_matrix_int": 0,
}

REQUIRED_FEATURES = frozenset(
    {
        "asym_id",
        "entity_id",
        "sym_id",
        "has_deletion",
        "deletion_mean",
        "deletion_value",
        "msa",
        "profile",
        "num_alignments",
        "residue_index",
        "prot_pair_num_alignments",
        "prot_unpair_num_alignments",
        "rna_pair_num_alignments",
        "rna_unpair_num_alignments",
    }
)

PROT_TYPE_NAME = "proteinChain"  # inference protein name in json


def make_dummy_msa_obj(input_sequence) -> parsers.Msa:
    deletion_matrix = [[0 for _ in input_sequence]]
    return parsers.Msa(
        sequences=[input_sequence],
        deletion_matrix=deletion_matrix,
        descriptions=["dummy"],
    )


def convert_monomer_features(monomer_features: FeatureDict) -> FeatureDict:
    """Reshapes and modifies monomer features for multimer models."""
    converted = {}
    unnecessary_leading_dim_feats = {
        "sequence",
        "domain_name",
        "num_alignments",
        "seq_length",
    }
    for feature_name, feature in monomer_features.items():
        if feature_name in unnecessary_leading_dim_feats:
            # asarray ensures it's a np.ndarray.
            feature = np.asarray(feature[0], dtype=feature.dtype)
        elif feature_name == "aatype":
            # The multimer model performs the one-hot operation itself.
            feature = np.argmax(feature, axis=-1).astype(np.int32)
        converted[feature_name] = feature
    return converted


def make_sequence_features(
    sequence: str,
    num_res: int,
    mapping: dict = residue_constants.restype_order_with_x,
    x_token: str = "X",
) -> FeatureDict:
    """
    Construct a feature dict of sequence features

    Args:
        sequence (str): input sequence
        num_res (int): number of residues in the input sequence

    Returns:
        FeatureDict: basic features of the input sequence
    """
    features = {}
    features["aatype"] = residue_constants.sequence_to_onehot(
        sequence=sequence,
        mapping=mapping,
        map_unknown_to_x=True,
        x_token=x_token,
    )
    features["between_segment_residues"] = np.zeros((num_res,), dtype=np.int32)
    features["residue_index"] = np.array(range(num_res), dtype=np.int32) + 1
    features["seq_length"] = np.array([num_res] * num_res, dtype=np.int32)
    features["sequence"] = np.array([sequence.encode("utf-8")], dtype=object)
    return features


def make_msa_features(
    msas: Sequence[parsers.Msa],
    identifier_func: Callable,
    mapping: tuple[dict] = (
        residue_constants.HHBLITS_AA_TO_ID,
        residue_constants.ID_TO_HHBLITS_AA,
    ),
) -> FeatureDict:
    """
        Constructs a feature dict of MSA features

    Args:
        msas (Sequence[parsers.Msa]): input MSA arrays
        identifier_func (Callable): the function extracting species identifier from MSA

    Returns:
        FeatureDict: raw MSA features
    """
    if not msas:
        raise ValueError("At least one MSA must be provided.")

    int_msa = []
    deletion_matrix = []
    species_ids = []
    seen_sequences = set()
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f"MSA {msa_index} must contain at least one sequence.")
        for sequence_index, sequence in enumerate(msa.sequences):
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            int_msa.append([mapping[0][res] for res in sequence])
            deletion_matrix.append(msa.deletion_matrix[sequence_index])

            identifiers = identifier_func(
                msa.descriptions[sequence_index],
            )

            species_ids.append(identifiers.species_id.encode("utf-8"))

    # residue type from HHBLITS_AA_TO_ID
    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {}
    features["deletion_matrix_int"] = np.array(deletion_matrix, dtype=np.int32)
    features["msa"] = np.array(int_msa, dtype=np.int32)
    features["num_alignments"] = np.array([num_alignments] * num_res, dtype=np.int32)
    features["msa_species_identifiers"] = np.array(species_ids, dtype=np.object_)
    features["profile"] = _make_msa_profile(
        msa=features["msa"], dict_size=len(mapping[1])
    )  # [num_res, 27]
    return features


def _make_msa_profile(msa: np.ndarray, dict_size: int) -> np.ndarray:
    """
    Make MSA profile (distribution over residues)

    Args:
        msas (Sequence[parsers.Msa]): input MSA arrays
        dict_size (int): number of residue types

    Returns:
        np.array: MSA profile
    """
    num_seqs = msa.shape[0]
    all_res_types = np.arange(dict_size)
    res_type_hits = msa[..., None] == all_res_types[None, ...]
    res_type_counts = res_type_hits.sum(axis=0)
    profile = res_type_counts / num_seqs
    return profile


def parse_a3m(path: str, seq_limit: int) -> tuple[list[str], list[str]]:
    """
    Parse a .a3m file

    Args:
        path (str): file path
        seq_limit (int): the max number of MSA sequences read from the file
            seq_limit > 0: real limit
            seq_limit = 0: return empty results
            seq_limit < 0: no limit, return all results

    Returns:
        tuple[list[str], list[str]]: parsed MSA sequences and their corresponding descriptions
    """
    sequences, descriptions = [], []
    if seq_limit == 0:
        return sequences, descriptions

    index = -1
    with open(path, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith(">"):
                if seq_limit > 0 and len(sequences) > seq_limit:
                    break
                index += 1
                descriptions.append(line[1:])  # Remove the '>' at the beginning.
                sequences.append("")
                continue
            elif line.startswith("#"):
                continue
            elif not line:
                continue  # Skip blank lines.
            sequences[index] += line

    return sequences, descriptions


def calc_stockholm_RNA_msa(
    name_to_sequence: OrderedDict, query: Optional[str]
) -> dict[str, parsers.Msa]:
    """
    Parses sequences and deletion matrix from stockholm format alignment.

    Args:
      stockholm_string: The string contents of a stockholm file. The first
        sequence in the file should be the query sequence.

    Returns:
      A tuple of:
        * A list of sequences that have been aligned to the query. These
          might contain duplicates.
        * The deletion matrix for the alignment as a list of lists. The element
          at `deletion_matrix[i][j]` is the number of residues deleted from
          the aligned sequence i at residue position j.
        * The names of the targets matched, including the jackhmmer subsequence
          suffix.
    """
    msa = []
    deletion_matrix = []

    if query is not None:
        # Add query string to the alignment
        if len(name_to_sequence.keys()) == 0:
            # name_to_sequence = OrderedDict({"query": query})
            return OrderedDict()
        else:
            query = align_query_to_sto(query, list(name_to_sequence.values())[0])
            new_name_to_sequence = OrderedDict({"query": query})
            new_name_to_sequence.update(name_to_sequence)
            name_to_sequence = new_name_to_sequence

    keep_columns = []
    for seq_index, sequence in enumerate(name_to_sequence.values()):
        if seq_index == 0:
            # Gather the columns with gaps from the query
            query = sequence
            keep_columns = [i for i, res in enumerate(query) if res != "-"]

        if len(sequence) < len(query):
            sequence = sequence + "-" * (len(query) - len(sequence))

        # Remove the columns with gaps in the query from all sequences.
        aligned_sequence = "".join([sequence[c] for c in keep_columns])
        # Convert lower case letter to upper case
        aligned_sequence = aligned_sequence.upper()
        msa.append(aligned_sequence)

        # Count the number of deletions w.r.t. query.
        deletion_vec = []
        deletion_count = 0
        for seq_res, query_res in zip(sequence, query):
            if seq_res not in ["-", "."] or query_res != "-":
                if query_res == "-":
                    deletion_count += 1
                else:
                    deletion_vec.append(deletion_count)
                    deletion_count = 0
        deletion_matrix.append(deletion_vec)

    return parsers.Msa(
        sequences=msa,
        deletion_matrix=deletion_matrix,
        descriptions=list(name_to_sequence.keys()),
    )


def align_query_to_sto(query: str, sto_sequence: str):
    """
    Aligns the query sequence to a Stockholm sequence by inserting gaps where necessary.

    Args:
        query (str): The query sequence to be aligned.
        sto_sequence (str): The Stockholm sequence to which the query is aligned.

    Returns:
        str: The aligned query sequence.
    """
    query = query.strip()
    sto_sequence = sto_sequence.strip()

    query_chars = []
    j = 0

    for i in range(len(sto_sequence)):
        if sto_sequence[i].islower() or sto_sequence[i] == ".":
            query_chars.append("-")
        else:
            query_chars.append(query[j])
            j += 1

    aligned_query = "".join(query_chars)

    if j < len(query):
        aligned_query += query[-(len(query) - j) :]

    return aligned_query


def parse_sto(path: str) -> OrderedDict:
    """
    Parses a Stockholm file and returns an ordered dictionary mapping sequence names to their sequences.

    Args:
        path (str): The path to the Stockholm file.

    Returns:
        OrderedDict: An ordered dictionary where keys are sequence names and values are the corresponding sequences.
    """
    name_to_sequence = OrderedDict()
    with open(path, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith(("#", "//")):
                continue
            name, sequence = line.split()
            if name not in name_to_sequence:
                name_to_sequence[name] = ""
            name_to_sequence[name] += sequence
    return name_to_sequence


def parse_msa_data(
    raw_msa_paths: Sequence[str],
    seq_limits: Sequence[int],
    msa_entity_type: str,
    query: Optional[str] = None,
) -> dict[str, parsers.Msa]:
    """
    Parses MSA data based on the entity type (protein or RNA).

    Args:
        raw_msa_paths (Sequence[str]): Paths to the MSA files.
        seq_limits (Sequence[int]): Limits on the number of sequences to read from each file.
        msa_entity_type (str): Type of MSA entity, either "prot" for protein or "rna" for RNA.
        query (Optional[str]): The query sequence for RNA MSA parsing. Defaults to None.

    Returns:
        dict[str, parsers.Msa]: A dictionary containing the parsed MSA data.

    Raises:
        ValueError: If `msa_entity_type` is not "prot" or "rna".
    """
    if msa_entity_type == "prot":
        return parse_prot_msa_data(raw_msa_paths, seq_limits)

    if msa_entity_type == "rna":
        return parse_rna_msa_data(raw_msa_paths, seq_limits, query=query)
    return []


def parse_rna_msa_data(
    raw_msa_paths: Sequence[str],
    seq_limits: Sequence[int],
    query: Optional[str] = None,
) -> dict[str, parsers.Msa]:
    """
    Parse MSAs for a sequence

    Args:
        raw_msa_paths (Sequence[str]): Paths of MSA files
        seq_limits (Sequence[int]): The max number of MSA sequences read from each file

    Returns:
        Dict[str, parsers.Msa]: MSAs parsed from each file
    """
    msa_data = {}
    for path, seq_limit in zip(raw_msa_paths, seq_limits):
        name_to_sequence = parse_sto(path)
        # The sto file has been truncated to a maximum length of seq_limit
        msa = calc_stockholm_RNA_msa(
            name_to_sequence=name_to_sequence,
            query=query,
        )
        if len(msa) > 0:
            msa_data[path] = msa
    return msa_data


def parse_prot_msa_data(
    raw_msa_paths: Sequence[str],
    seq_limits: Sequence[int],
) -> dict[str, parsers.Msa]:
    """
    Parse MSAs for a sequence

    Args:
        raw_msa_paths (Sequence[str]): Paths of MSA files
        seq_limits (Sequence[int]): The max number of MSA sequences read from each file

    Returns:
        Dict[str, parsers.Msa]: MSAs parsed from each file
    """
    msa_data = {}
    for path, seq_limit in zip(raw_msa_paths, seq_limits):
        sequences, descriptions = parse_a3m(path, seq_limit)

        deletion_matrix = []
        for msa_sequence in sequences:
            deletion_vec = []
            deletion_count = 0
            for j in msa_sequence:
                if j.islower():
                    deletion_count += 1
                else:
                    deletion_vec.append(deletion_count)
                    deletion_count = 0
            deletion_matrix.append(deletion_vec)

        # Make the MSA matrix out of aligned (deletion-free) sequences.
        deletion_table = str.maketrans("", "", string.ascii_lowercase)
        aligned_sequences = [s.translate(deletion_table) for s in sequences]
        assert all([len(seq) == len(aligned_sequences[0]) for seq in aligned_sequences])
        assert all([len(vec) == len(deletion_matrix[0]) for vec in deletion_matrix])

        if len(aligned_sequences) > 0:
            # skip empty file
            msa = parsers.Msa(
                sequences=aligned_sequences,
                deletion_matrix=deletion_matrix,
                descriptions=descriptions,
            )
            msa_data[path] = msa

    return msa_data


def load_and_process_msa(
    pdb_name: str,
    msa_type: str,
    raw_msa_paths: Sequence[str],
    seq_limits: Sequence[int],
    identifier_func: Optional[Callable] = lambda x: Identifiers(),
    input_sequence: Optional[str] = None,
    handle_empty: str = "return_self",
    msa_entity_type: str = "prot",
) -> dict[str, Any]:
    """
    Load and process MSA features of a single sequence

    Args:
        pdb_name (str): f"{pdb_id}_{entity_id}" of the input entity
        msa_type (str): Type of MSA ("pairing" or "non_pairing")
        raw_msa_paths (Sequence[str]): Paths of MSA files
        identifier_func (Optional[Callable]): The function extracting species identifier from MSA
        input_sequence (str): The input sequence
        handle_empty (str): How to handle empty MSA ("return_self" or "raise_error")
        entity_type (str): rna or prot

    Returns:
        Dict[str, Any]: processed MSA features
    """
    msa_data = parse_msa_data(
        raw_msa_paths, seq_limits, msa_entity_type=msa_entity_type, query=input_sequence
    )
    if len(msa_data) == 0:
        if handle_empty == "return_self":
            msa_data["dummy"] = make_dummy_msa_obj(input_sequence)
        elif handle_empty == "raise_error":
            ValueError(f"No valid {msa_type} MSA for {pdb_name}")
        else:
            raise NotImplementedError(
                f"Unimplemented empty-handling method: {handle_empty}"
            )
    msas = list(msa_data.values())

    if msa_type == "non_pairing":
        return make_msa_features(
            msas=msas,
            identifier_func=identifier_func,
            mapping=(
                (residue_constants.HHBLITS_AA_TO_ID, residue_constants.ID_TO_HHBLITS_AA)
                if msa_entity_type == "prot"
                else (RNA_NT_TO_ID, RNA_ID_TO_NT)
            ),
        )
    elif msa_type == "pairing":
        all_seq_features = make_msa_features(
            msas=msas,
            identifier_func=identifier_func,
            mapping=(
                (residue_constants.HHBLITS_AA_TO_ID, residue_constants.ID_TO_HHBLITS_AA)
                if msa_entity_type == "prot"
                else (RNA_NT_TO_ID, RNA_ID_TO_NT)
            ),
        )
        valid_feats = MSA_FEATURES + ("msa_species_identifiers",)
        return {
            f"{k}_all_seq": v for k, v in all_seq_features.items() if k in valid_feats
        }


def add_assembly_features(
    pdb_id: str,
    all_chain_features: MutableMapping[str, FeatureDict],
    asym_to_entity_id: Mapping[int, str],
) -> dict[str, FeatureDict]:
    """
    Add features to distinguish between chains.

    Args:
        all_chain_features (MutableMapping[str, FeatureDict]): A dictionary which maps chain_id to a dictionary of features for each chain.
        asym_to_entity_id (Mapping[int, str]): A mapping from asym_id_int to entity_id

    Returns:
        all_chain_features (MutableMapping[str, FeatureDict]): all_chain_features with assembly features added
    """
    # Group the chains by entity
    grouped_chains = defaultdict(list)
    for asym_id_int, chain_features in all_chain_features.items():
        entity_id = asym_to_entity_id[asym_id_int]
        chain_features["asym_id"] = asym_id_int
        grouped_chains[entity_id].append(chain_features)

    new_all_chain_features = {}
    for entity_id, group_chain_features in grouped_chains.items():
        assert int(entity_id) >= 0
        for sym_id, chain_features in enumerate(group_chain_features, start=1):
            new_all_chain_features[f"{entity_id}_{sym_id}"] = chain_features
            seq_length = chain_features["seq_length"]
            chain_features["asym_id"] = (
                chain_features["asym_id"] * np.ones(seq_length)
            ).astype(np.int64)
            chain_features["sym_id"] = (sym_id * np.ones(seq_length)).astype(np.int64)
            chain_features["entity_id"] = (int(entity_id) * np.ones(seq_length)).astype(
                np.int64
            )
            chain_features["pdb_id"] = pdb_id
    return new_all_chain_features


def process_unmerged_features(
    all_chain_features: MutableMapping[str, Mapping[str, np.ndarray]]
):
    """
    Postprocessing stage for per-chain features before merging

    Args:
        all_chain_features (MutableMapping[str, Mapping[str, np.ndarray]]): MSA features of all chains

    Returns:
        post-processed per-chain features
    """
    for chain_features in all_chain_features.values():
        # Convert deletion matrices to float.
        chain_features["deletion_matrix"] = np.asarray(
            chain_features.pop("deletion_matrix_int"), dtype=np.float32
        )
        if "deletion_matrix_int_all_seq" in chain_features:
            chain_features["deletion_matrix_all_seq"] = np.asarray(
                chain_features.pop("deletion_matrix_int_all_seq"), dtype=np.float32
            )

        chain_features["deletion_mean"] = np.mean(
            chain_features["deletion_matrix"], axis=0
        )


def pair_and_merge(
    is_homomer_or_monomer: bool,
    all_chain_features: MutableMapping[str, Mapping[str, np.ndarray]],
    merge_method: str,
    msa_crop_size: int,
) -> dict[str, np.ndarray]:
    """
    Runs processing on features to augment, pair and merge

    Args:
        is_homomer_or_monomer (bool): True if the bioassembly is a homomer or a monomer
        all_chain_features (MutableMapping[str, Mapping[str, np.ndarray]]):
            A MutableMap of dictionaries of features for each chain.
        merge_method (str): How to merge unpaired MSA features

    Returns:
        Dict[str, np.ndarray]: A dictionary of features
    """

    process_unmerged_features(all_chain_features)

    np_chains_list = list(all_chain_features.values())

    pair_msa_sequences = not is_homomer_or_monomer

    if pair_msa_sequences:
        np_chains_list = create_paired_features(chains=np_chains_list)
        np_chains_list = deduplicate_unpaired_sequences(np_chains_list)

    np_chains_list = crop_chains(
        np_chains_list,
        msa_crop_size=msa_crop_size,
        pair_msa_sequences=pair_msa_sequences,
    )

    np_example = merge_chain_features(
        np_chains_list=np_chains_list,
        pair_msa_sequences=pair_msa_sequences,
        merge_method=merge_method,
        msa_entity_type="prot",
    )

    np_example = process_prot_final(np_example)

    return np_example


def rna_merge(
    all_chain_features: MutableMapping[str, Mapping[str, np.ndarray]],
    merge_method: str,
    msa_crop_size: int,
) -> dict[str, np.ndarray]:
    """
    Runs processing on features to augment and merge

    Args:
        all_chain_features (MutableMapping[str, Mapping[str, np.ndarray]]):
            A MutableMap of dictionaries of features for each chain.
        merge_method (str): how to merge unpaired MSA features

    Returns:
        Dict[str, np.ndarray]: A dictionary of features
    """
    process_unmerged_features(all_chain_features)
    np_chains_list = list(all_chain_features.values())
    np_chains_list = crop_chains(
        np_chains_list,
        msa_crop_size=msa_crop_size,
        pair_msa_sequences=False,
    )
    np_example = merge_chain_features(
        np_chains_list=np_chains_list,
        pair_msa_sequences=False,
        merge_method=merge_method,
        msa_entity_type="rna",
    )
    np_example = process_rna_final(np_example)

    return np_example


def merge_chain_features(
    np_chains_list: list[Mapping[str, np.ndarray]],
    pair_msa_sequences: bool,
    merge_method: str,
    msa_entity_type: str = "prot",
) -> Mapping[str, np.ndarray]:
    """
    Merges features for multiple chains to single FeatureDict

    Args:
        np_chains_list (List[Mapping[str, np.ndarray]]): List of FeatureDicts for each chain
        pair_msa_sequences (bool): Whether to merge paired MSAs
        merge_method (str): how to merge unpaired MSA features
        msa_entity_type (str): protein or rna

    Returns:
        Single FeatureDict for entire bioassembly
    """
    np_chains_list = _merge_homomers_dense_features(
        np_chains_list, merge_method, msa_entity_type=msa_entity_type
    )

    # Unpaired MSA features will be always block-diagonalised; paired MSA
    # features will be concatenated.
    np_example = _merge_features_from_multiple_chains(
        np_chains_list,
        pair_msa_sequences=False,
        merge_method=merge_method,
        msa_entity_type=msa_entity_type,
    )

    if pair_msa_sequences:
        np_example = _concatenate_paired_and_unpaired_features(np_example)

    np_example = _correct_post_merged_feats(
        np_example=np_example,
    )

    return np_example


def _merge_homomers_dense_features(
    chains: Iterable[Mapping[str, np.ndarray]],
    merge_method: str,
    msa_entity_type: str = "prot",
) -> list[dict[str, np.ndarray]]:
    """
    Merge all identical chains, making the resulting MSA dense

    Args:
        chains (Iterable[Mapping[str, np.ndarray]]): An iterable of features for each chain
        merge_method (str): how to merge unpaired MSA features
        msa_entity_type (str): protein or rna
    Returns:
        List[Dict[str, np.ndarray]]: A list of feature dictionaries. All features with the same entity_id
        will be merged - MSA features will be concatenated along the num_res dimension - making them dense.
    """
    entity_chains = defaultdict(list)
    for chain in chains:
        entity_id = chain["entity_id"][0]
        entity_chains[entity_id].append(chain)

    grouped_chains = []
    for entity_id in sorted(entity_chains):
        chains = entity_chains[entity_id]
        grouped_chains.append(chains)
    chains = [
        _merge_features_from_multiple_chains(
            chains,
            pair_msa_sequences=True,
            merge_method=merge_method,
            msa_entity_type=msa_entity_type,
        )
        for chains in grouped_chains
    ]
    return chains


def _merge_msa_features(
    *feats: np.ndarray,
    feature_name: str,
    merge_method: str,
    msa_entity_type: str = "prot",
) -> np.ndarray:
    """
    Merge unpaired MSA features

    Args:
        feats (np.ndarray): input features
        feature_name (str): feature name
        merge_method (str): how to merge unpaired MSA features
        msa_entity_type (str): protein or rna
    Returns:
        np.ndarray: merged feature
    """
    assert msa_entity_type in ["prot", "rna"]
    if msa_entity_type == "prot":
        mapping = MSA_PAD_VALUES
    elif msa_entity_type == "rna":
        mapping = RNA_MSA_PAD_VALUES
    if merge_method == "sparse":
        merged_feature = block_diag(*feats, pad_value=mapping[feature_name])
    elif merge_method in ["dense_min"]:
        merged_feature = truncate_at_min(*feats)
    elif merge_method in ["dense_max"]:
        merged_feature = pad_to_max(*feats, pad_value=mapping[feature_name])
    else:
        raise NotImplementedError(
            f"Unknown merge method {merge_method}! Allowed merged methods are: "
        )
    return merged_feature


def _merge_features_from_multiple_chains(
    chains: Sequence[Mapping[str, np.ndarray]],
    pair_msa_sequences: bool,
    merge_method: str,
    msa_entity_type: str = "prot",
) -> dict[str, np.ndarray]:
    """
    Merge features from multiple chains.

    Args:
        chains (Sequence[Mapping[str, np.ndarray]]):
            A list of feature dictionaries that we want to merge
        pair_msa_sequences (bool): Whether to concatenate MSA features along the
            num_res dimension (if True), or to block diagonalize them (if False)
        merge_method (str): how to merge unpaired MSA features
        msa_entity_type (str): protein or rna
    Returns:
        Dict[str, np.ndarray]: A feature dictionary for the merged example
    """
    merged_example = {}
    for feature_name in chains[0]:
        feats = [x[feature_name] for x in chains]
        feature_name_split = feature_name.split("_all_seq")[0]
        if feature_name_split in MSA_FEATURES:
            if pair_msa_sequences or "_all_seq" in feature_name:
                merged_example[feature_name] = np.concatenate(feats, axis=1)
            else:
                merged_example[feature_name] = _merge_msa_features(
                    *feats,
                    merge_method=merge_method,
                    feature_name=feature_name,
                    msa_entity_type=msa_entity_type,
                )
        elif feature_name_split in SEQ_FEATURES:
            merged_example[feature_name] = np.concatenate(feats, axis=0)
        elif feature_name_split in CHAIN_FEATURES:
            merged_example[feature_name] = np.sum([x for x in feats]).astype(np.int32)
        else:
            merged_example[feature_name] = feats[0]
    return merged_example


def merge_features_from_prot_rna(
    chains: Sequence[Mapping[str, np.ndarray]],
) -> dict[str, np.ndarray]:
    """
    Merge features from prot and rna chains.

    Args:
        chains (Sequence[Mapping[str, np.ndarray]]):
            A list of feature dictionaries that we want to merge

    Returns:
        Dict[str, np.ndarray]: A feature dictionary for the merged example
    """
    merged_example = {}
    if len(chains) == 1:  # only prot or rna msa exists
        return chains[0]
    final_msa_pad_values = {
        "msa": 31,
        "has_deletion": False,
        "deletion_value": 0,
    }
    for feature_name in set(chains[0].keys()).union(chains[1].keys()):
        feats = [x[feature_name] for x in chains if feature_name in x]
        if (
            feature_name in SEQ_FEATURES
        ):  # ["residue_index", "profile", "asym_id", "sym_id", "entity_id", "deletion_mean"]
            merged_example[feature_name] = np.concatenate(feats, axis=0)
        elif feature_name in ["msa", "has_deletion", "deletion_value"]:
            merged_example[feature_name] = pad_to_max(
                *feats, pad_value=final_msa_pad_values[feature_name]
            )
        elif feature_name in [
            "prot_pair_num_alignments",
            "prot_unpair_num_alignments",
            "rna_pair_num_alignments",
            "rna_unpair_num_alignments",
        ]:  # unmerged keys keep for tracking
            merged_example[feature_name] = feats[0]
        else:
            continue
    merged_example["num_alignments"] = np.asarray(
        merged_example["msa"].shape[0], dtype=np.int32
    )
    return merged_example


def _concatenate_paired_and_unpaired_features(
    np_example: Mapping[str, np.ndarray]
) -> dict[str, np.ndarray]:
    """
    Concatenate paired and unpaired features

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: features with paired and unpaired features concatenated
    """
    features = MSA_FEATURES
    for feature_name in features:
        if feature_name in np_example:
            feat = np_example[feature_name]
            feat_all_seq = np_example[feature_name + "_all_seq"]
            merged_feat = np.concatenate([feat_all_seq, feat], axis=0)
            np_example[feature_name] = merged_feat
    np_example["num_alignments"] = np.array(np_example["msa"].shape[0], dtype=np.int32)
    return np_example


def _correct_post_merged_feats(
    np_example: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """
    Adds features that need to be computed/recomputed post merging

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """

    np_example["seq_length"] = np.asarray(np_example["aatype"].shape[0], dtype=np.int32)
    np_example["num_alignments"] = np.asarray(
        np_example["msa"].shape[0], dtype=np.int32
    )
    return np_example


def _add_msa_num_alignment(
    np_example: Mapping[str, np.ndarray], msa_entity_type: str = "prot"
) -> dict[str, np.ndarray]:
    """
    Adds pair and unpair msa alignments num

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    assert msa_entity_type in ["prot", "rna"]
    if "msa_all_seq" in np_example:
        pair_num_alignments = np.asarray(
            np_example["msa_all_seq"].shape[0], dtype=np.int32
        )
    else:
        pair_num_alignments = np.asarray(0, dtype=np.int32)
    np_example[f"{msa_entity_type}_pair_num_alignments"] = pair_num_alignments
    np_example[f"{msa_entity_type}_unpair_num_alignments"] = (
        np_example["num_alignments"] - pair_num_alignments
    )
    return np_example


def process_prot_final(np_example: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Final processing steps in data pipeline, after merging and pairing

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    np_example = correct_msa_restypes(np_example)
    np_example = final_transform(np_example)
    np_example = _add_msa_num_alignment(np_example, msa_entity_type="prot")
    np_example = filter_features(np_example)

    return np_example


def correct_msa_restypes(np_example: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Correct MSA restype to have the same order as residue_constants

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    # remap msa
    np_example["msa"] = np.take(NEW_ORDER_LIST, np_example["msa"], axis=0)
    np_example["msa"] = np_example["msa"].astype(np.int32)

    seq_len, profile_dim = np_example["profile"].shape
    assert profile_dim == len(NEW_ORDER_LIST)
    profile = np.zeros((seq_len, 32))
    profile[:, np.array(NEW_ORDER_LIST)] = np_example["profile"]

    np_example["profile"] = profile
    return np_example


def process_rna_final(np_example: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Final processing steps in data pipeline, after merging and pairing

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    np_example = correct_rna_msa_restypes(np_example)
    np_example = final_transform(np_example)
    np_example = _add_msa_num_alignment(np_example, msa_entity_type="rna")
    np_example = filter_features(np_example)

    return np_example


def correct_rna_msa_restypes(
    np_example: Mapping[str, np.ndarray]
) -> dict[str, np.ndarray]:
    """
    Correct MSA restype to have the same order as residue_constants

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    # remap msa
    np_example["msa"] = np.take(RNA_NEW_ORDER_LIST, np_example["msa"], axis=0)
    np_example["msa"] = np_example["msa"].astype(np.int32)

    seq_len, profile_dim = np_example["profile"].shape
    assert profile_dim == len(RNA_NEW_ORDER_LIST)
    profile = np.zeros((seq_len, 32))
    profile[:, np.array(RNA_NEW_ORDER_LIST)] = np_example["profile"]

    np_example["profile"] = profile
    return np_example


def final_transform(np_example: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Performing some transformations related to deletion_matrix

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    deletion_mat = np_example.pop("deletion_matrix")
    np_example["has_deletion"] = np.clip(deletion_mat, a_min=0, a_max=1).astype(
        np.bool_
    )

    np_example["deletion_value"] = (2 / np.pi) * np.arctan(deletion_mat / 3)
    assert np.all(-1e-5 < np_example["deletion_value"]) and np.all(
        np_example["deletion_value"] < (1 + 1e-5)
    )
    return np_example


def filter_features(np_example: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Filters features of example to only those requested

    Args:
        np_example (Mapping[str, np.ndarray]): input features

    Returns:
        Dict[str, np.ndarray]: processed features
    """
    return {k: v for (k, v) in np_example.items() if k in REQUIRED_FEATURES}


def crop_chains(
    chains_list: Sequence[Mapping[str, np.ndarray]],
    msa_crop_size: int,
    pair_msa_sequences: bool,
) -> list[Mapping[str, np.ndarray]]:
    """
    Crops the MSAs for a set of chains.

    Args:
        chains_list (Sequence[Mapping[str, np.ndarray]]): A list of chains to be cropped.
        msa_crop_size (int): The total number of sequences to crop from the MSA.
        pair_msa_sequences (bool): Whether we are operating in sequence-pairing mode.

    Returns:
        List[Mapping[str, np.ndarray]]: The chains cropped
    """

    # Apply the cropping.
    cropped_chains = []
    for chain in chains_list:
        cropped_chain = _crop_single_chain(
            chain,
            msa_crop_size=msa_crop_size,
            pair_msa_sequences=pair_msa_sequences,
        )
        cropped_chains.append(cropped_chain)

    return cropped_chains


def _crop_single_chain(
    chain: Mapping[str, np.ndarray],
    msa_crop_size: int,  # 2048
    pair_msa_sequences: bool,
) -> dict[str, np.ndarray]:
    """
    Crops msa sequences to msa_crop_size

    Args:
        chain (Mapping[str, np.ndarray]): The chain to be cropped
        msa_crop_size (int): The total number of sequences to crop from the MSA
        pair_msa_sequences (bool): Whether we are operating in sequence-pairing mode

    Returns:
        Dict[str, np.ndarray]: The chains cropped
    """
    msa_size = chain["num_alignments"]

    if pair_msa_sequences:
        msa_size_all_seq = chain["num_alignments_all_seq"]
        msa_crop_size_all_seq = np.minimum(msa_size_all_seq, msa_crop_size // 2)

        # We reduce the number of un-paired sequences, by the number of times a
        # sequence from this chain's MSA is included in the paired MSA.  This keeps
        # the MSA size for each chain roughly constant.
        msa_all_seq = chain["msa_all_seq"][:msa_crop_size_all_seq, :]
        num_non_gapped_pairs = np.sum(np.any(msa_all_seq != MSA_GAP_IDX, axis=1))
        num_non_gapped_pairs = np.minimum(num_non_gapped_pairs, msa_crop_size_all_seq)

        # Restrict the unpaired crop size so that paired+unpaired sequences do not
        # exceed msa_seqs_per_chain for each chain.
        max_msa_crop_size = np.maximum(msa_crop_size - num_non_gapped_pairs, 0)
        msa_crop_size = np.minimum(msa_size, max_msa_crop_size)
    else:
        msa_crop_size = np.minimum(msa_size, msa_crop_size)

    for k in chain:
        k_split = k.split("_all_seq")[0]
        if k_split in MSA_FEATURES:
            if "_all_seq" in k and pair_msa_sequences:
                chain[k] = chain[k][:msa_crop_size_all_seq, :]
            else:
                chain[k] = chain[k][:msa_crop_size, :]

    chain["num_alignments"] = np.asarray(msa_crop_size, dtype=np.int32)
    if pair_msa_sequences:
        chain["num_alignments_all_seq"] = np.asarray(
            msa_crop_size_all_seq, dtype=np.int32
        )
    return chain


def truncate_at_min(*arrs: np.ndarray) -> np.ndarray:
    """
    Processing unpaired features by truncating at the min length

    Args:
        arrs (np.ndarray): input features

    Returns:
        np.ndarray: truncated features
    """
    min_num_msa = min([x.shape[0] for x in arrs])
    truncated_arrs = [x[:min_num_msa, :] for x in arrs]
    new_arrs = np.concatenate(truncated_arrs, axis=1)
    return new_arrs


def pad_to_max(*arrs: np.ndarray, pad_value: float = 0.0) -> np.ndarray:
    """
    Processing unpaired features by padding to the max length

    Args:
        arrs (np.ndarray): input features

    Returns:
        np.ndarray: padded features
    """
    max_num_msa = max([x.shape[0] for x in arrs])
    padded_arrs = [
        np.pad(
            x,
            ((0, (max_num_msa - x.shape[0])), (0, 0)),
            mode="constant",
            constant_values=pad_value,
        )
        for x in arrs
    ]
    new_arrs = np.concatenate(padded_arrs, axis=1)
    return new_arrs


def clip_msa(
    np_example: Mapping[str, np.ndarray], max_num_msa: int
) -> dict[str, np.ndarray]:
    """
    Clip MSA features to a maximum length

    Args:
        np_example (Mapping[str, np.ndarray]): input MSA features
        pad_value (float): pad value

    Returns:
        Dict[str, np.ndarray]: clipped MSA features
    """
    if np_example["msa"].shape[0] > max_num_msa:
        for k in ["msa", "has_deletion", "deletion_value"]:
            np_example[k] = np_example[k][:max_num_msa, :]
        np_example["num_alignments"] = max_num_msa
        assert np_example["num_alignments"] == np_example["msa"].shape[0]
    return np_example


def get_identifier_func(pairing_db: str) -> Callable:
    """
    Get the function the extracts species identifier from sequence descriptions

    Args:
        pairing_db (str): the database from which MSAs for pairing are searched

    Returns:
        Callable: the function the extracts species identifier from sequence descriptions
    """
    if pairing_db.startswith("uniprot"):

        def func(description: str) -> Identifiers:
            sequence_identifier = _extract_sequence_identifier(description)
            if sequence_identifier is None:
                return Identifiers()
            else:
                return _parse_sequence_identifier(sequence_identifier)

        return func

    elif pairing_db.startswith("uniref100"):

        def func(description: str) -> Identifiers:
            if (
                description.startswith("UniRef100")
                and "/" in description
                and (first_comp := description.split("/")[0]).count("_") == 2
            ):
                identifier = Identifiers(species_id=first_comp.split("_")[-1])
            else:
                identifier = Identifiers()
            return identifier

        return func
    else:
        raise NotImplementedError(
            f"Identifier func for {pairing_db} is not implemented"
        )


def run_msa_tool(
    msa_runner,
    fasta_path: str,
    msa_out_path: str,
    msa_format: str,
    max_sto_sequences: Optional[int] = None,
) -> Mapping[str, Any]:
    """Runs an MSA tool, checking if output already exists first."""
    if msa_format == "sto" and max_sto_sequences is not None:
        result = msa_runner.query(fasta_path, max_sto_sequences)[0]
    else:
        result = msa_runner.query(fasta_path)[0]

    assert msa_out_path.split(".")[-1] == msa_format
    with open(msa_out_path, "w") as f:
        f.write(result[msa_format])

    return result


def search_msa(sequence: str, db_fpath: str, res_fpath: str = ""):
    assert opexists(
        db_fpath
    ), f"Database path for MSA searching does not exists:\n{db_fpath}"
    seq_name = uuid.uuid4().hex
    db_name = os.path.basename(db_fpath)
    jackhmmer_binary_path = shutil.which("jackhmmer")
    msa_runner = jackhmmer.Jackhmmer(
        binary_path=jackhmmer_binary_path,
        database_path=db_fpath,
        n_cpu=2,
    )
    if res_fpath == "":
        tmp_dir = f"/tmp/{uuid.uuid4().hex}"
        res_fpath = os.path.join(tmp_dir, f"{seq_name}.a3m")
    else:
        tmp_dir = os.path.dirname(res_fpath)
        os.makedirs(tmp_dir, exist_ok=True)
    output_sto_path = os.path.join(tmp_dir, f"{seq_name}.sto")
    with open((tmp_fasta_path := f"{tmp_dir}/{seq_name}_{db_name}.fasta"), "w") as f:
        f.write(f">query\n")
        f.write(sequence)

    logger.info(f"Searching MSA for {seq_name}\n. Will be saved to {output_sto_path}")
    _ = run_msa_tool(msa_runner, tmp_fasta_path, output_sto_path, "sto")
    if not opexists(output_sto_path):
        logger.info(f"Failed to search MSA for {sequence} from the database {db_fpath}")
        return

    logger.info(f"Reformatting the MSA file. Will be saved to {res_fpath}")

    cmd = f"/opt/hhsuite/scripts/reformat.pl {output_sto_path} {res_fpath}"
    try:
        subprocess.check_call(cmd, shell=True, executable="/bin/bash")
    except Exception as e:
        logger.info(f"Reformatting failed:\n {e}\nRetry {cmd}...")
        time.sleep(1)
        subprocess.check_call(cmd, shell=True, executable="/bin/bash")
    if not os.path.exists(res_fpath):
        logger.info(
            f"Failed to reformat the MSA file. Please check the validity of the .sto file{output_sto_path}"
        )
        return


def search_msa_paired(
    sequence: str, pairing_db_fpath: str, non_pairing_db_fpath: str, idx: int = -1
) -> tuple[Union[str, None], int]:
    tmp_dir = f"/tmp/{uuid.uuid4().hex}_{str(time.time()).replace('.', '_')}_{DIST_WRAPPER.rank}_{idx}"
    os.makedirs(tmp_dir, exist_ok=True)
    pairing_file = os.path.join(tmp_dir, "pairing.a3m")
    search_msa(sequence, pairing_db_fpath, pairing_file)
    if not os.path.exists(pairing_file):
        return None, idx
    non_pairing_file = os.path.join(tmp_dir, "non_pairing.a3m")
    search_msa(sequence, non_pairing_db_fpath, non_pairing_file)
    if not os.path.exists(non_pairing_file):
        return None, idx
    else:
        return tmp_dir, idx


def msa_parallel(sequences: dict[int, tuple[str, str, str]]) -> dict[int, str]:
    from concurrent.futures import ThreadPoolExecutor

    num_threads = 4
    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(search_msa_paired, seq[0], seq[1], seq[2], idx)
            for idx, seq in sequences.items()
        ]
        # Wait for all threads to complete
        for future in futures:
            results.append(future.result())
    msa_res = {}
    for x in results:
        msa_res[x[1]] = x[0]
    return msa_res
