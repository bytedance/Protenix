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

import json
import os
import shutil
from abc import ABC, abstractmethod
from collections import defaultdict
from copy import deepcopy
from os.path import exists as opexists
from os.path import join as opjoin
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np
import torch
from biotite.structure import AtomArray

from protenix.data.constants import STD_RESIDUES, rna_order_with_x
from protenix.data.msa_utils import (
    PROT_TYPE_NAME,
    FeatureDict,
    add_assembly_features,
    clip_msa,
    convert_monomer_features,
    get_identifier_func,
    load_and_process_msa,
    make_sequence_features,
    merge_features_from_prot_rna,
    msa_parallel,
    pair_and_merge,
    rna_merge,
)
from protenix.data.tokenizer import TokenArray
from protenix.utils.logger import get_logger

logger = get_logger(__name__)

SEQ_LIMITS = {
    "uniref100": -1,
    "mmseqs_other": -1,
    "uniclust30": -1,
    "rfam": 10000,
    "rnacentral": 10000,
    "nucleotide": 10000,
}
MSA_MAX_SIZE = 16384


class BaseMSAFeaturizer(ABC):
    def __init__(
        self,
        indexing_method: str = "sequence",
        merge_method: str = "dense_max",
        seq_limits: Optional[dict[str, int]] = {},
        max_size: int = 16384,
        **kwargs,
    ):
        """
        Initializes the BaseMSAFeaturizer with the specified parameters.

        Args:
            indexing_method (str): The method used for indexing the MSA. Defaults to "sequence".
            merge_method (str): The method used for merging MSA features. Defaults to "dense_max".
            seq_limits (Optional[dict[str, int]]): Dictionary specifying sequence limits for different databases. Defaults to an empty dictionary.
            max_size (int): The maximum size of the MSA. Defaults to 16384.
            **kwargs: Additional keyword arguments.

        Raises:
            AssertionError: If the provided `merge_method` or `indexing_method` is not valid.
        """
        assert merge_method in ["dense_max", "dense_min", "sparse"]
        assert indexing_method in [
            "sequence",
            "pdb_id",
            "pdb_id_entity_id",
        ], f"Unknown indexing method: {indexing_method}"
        self.indexing_method = indexing_method
        self.merge_method = merge_method
        self.seq_limits = seq_limits
        self.max_size = max_size

    @abstractmethod
    def get_msa_path(self):
        pass

    @abstractmethod
    def process_single_sequence(self):
        pass

    def get_entity_ids(
        self, bioassembly_dict: Mapping[str, Any], msa_entity_type: str = "prot"
    ) -> set[str]:
        """
        Extracts the entity IDs that match the specified MSA entity type from the bioassembly dictionary.

        Args:
            bioassembly_dict (Mapping[str, Any]): The bioassembly dictionary containing entity information.
            msa_entity_type (str): The type of MSA entity to filter by. Defaults to "prot".

        Returns:
            set[str]: A set of entity IDs that match the specified MSA entity type.

        Raises:
            AssertionError: If the provided `msa_entity_type` is not "prot" or "rna".
        """
        assert msa_entity_type in ["prot", "rna"], "only protein and rna might have msa"
        poly_type_mapping = {
            "prot": "polypeptide",
            "rna": "polyribonucleotide",
            "dna": "polydeoxyribonucleotide",
        }
        entity_poly_type = bioassembly_dict["entity_poly_type"]

        entity_ids: set[str] = {
            entity_id
            for entity_id, poly_type in entity_poly_type.items()
            if poly_type_mapping[msa_entity_type] in poly_type
        }
        return entity_ids

    def get_selected_asym_ids(
        self,
        bioassembly_dict: Mapping[str, Any],
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
        selected_token_indices: Optional[torch.Tensor],
        entity_ids: set[str],
    ) -> tuple[set[int], set[int], dict[int, str], dict[int, str], dict[str, str]]:
        """
        Extracts the selected asym IDs based on the provided bioassembly dictionary and entity IDs.

        Args:
            bioassembly_dict (Mapping[str, Any]): The bioassembly dictionary containing entity information.
            entity_to_asym_id_int (Mapping[str, Sequence[int]]): Mapping from entity ID to asym ID integers.
            selected_token_indices (Optional[torch.Tensor]): Indices of selected tokens.
            entity_ids (set[str]): Set of entity IDs to consider.

        Returns:
            tuple: A tuple containing:
                - selected_asym_ids (set[int]): Set of selected asym IDs.
                - asym_id_ints (set[int]): Set of asym ID integers.
                - asym_to_entity_id (dict[int, str]): Mapping from asym ID integers to entity IDs.
                - asym_id_int_to_sequence (dict[int, str]): Mapping from asym ID integers to sequences.
                - entity_id_to_sequence (dict[str, str]): Mapping from entity IDs to sequences.
        """
        asym_to_entity_id: dict[int, str] = {}
        # Only count the selected Prot/RNA entities, many-to-one mapping
        for entity_id, asym_id_int_list in entity_to_asym_id_int.items():
            if entity_id in entity_ids:
                for asym_id_int in asym_id_int_list:
                    asym_to_entity_id[asym_id_int] = entity_id
        entity_id_to_sequence = {
            k: v
            for (k, v) in bioassembly_dict["sequences"].items()
            if k in entity_ids and k in entity_to_asym_id_int
        }
        asym_id_ints = set(
            [
                asym_id_int
                for (asym_id_int, entity_id) in asym_to_entity_id.items()
                if entity_id in entity_ids
            ]
        )
        # Only count Prot/RNA chains, many-to-one mapping
        asym_id_int_to_sequence = {
            asym_id_int: entity_id_to_sequence[entity_id]
            for (asym_id_int, entity_id) in asym_to_entity_id.items()
        }
        atom_array = bioassembly_dict["atom_array"]
        token_array = bioassembly_dict["token_array"]

        if selected_token_indices is None:
            selected_asym_ids = set(
                [
                    atom_array[idx].asym_id_int
                    for idx in token_array.get_annotation("centre_atom_index")
                ]
            )
        else:
            selected_asym_ids = set(
                [
                    atom_array[idx].asym_id_int
                    for idx in token_array[selected_token_indices].get_annotation(
                        "centre_atom_index"
                    )
                ]
            )
        return (
            selected_asym_ids,
            asym_id_ints,
            asym_to_entity_id,
            asym_id_int_to_sequence,
            entity_id_to_sequence,
        )

    def get_msa_pipeline(
        self,
        is_homomer_or_monomer: bool,
        selected_asym_ids: set[int],
        asym_to_entity_id: dict[int, str],
        asym_id_int_to_sequence: dict[int, str],
        entity_id_to_sequence: dict[str, str],
        bioassembly_dict: Mapping[str, Any],
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
        msa_entity_type="prot",
    ) -> Optional[dict[str, np.ndarray]]:
        """
        Processes the MSA pipeline for the given bioassembly dictionary and selected asym IDs.

        Args:
            is_homomer_or_monomer (bool): Indicates if the sequence is a homomer or monomer.
            selected_asym_ids (set[int]): Set of selected asym IDs.
            asym_to_entity_id (dict[int, str]): Mapping from asym ID integers to entity IDs.
            asym_id_int_to_sequence (dict[int, str]): Mapping from asym ID integers to sequences.
            entity_id_to_sequence (dict[str, str]): Mapping from entity IDs to sequences.
            bioassembly_dict (Mapping[str, Any]): The bioassembly dictionary containing entity information.
            entity_to_asym_id_int (Mapping[str, Sequence[int]]): Mapping from entity ID to asym ID integers.
            msa_entity_type (str): The type of MSA entity to process. Defaults to "prot".

        Returns:
            Optional[dict[str, np.ndarray]]: A dictionary containing the processed MSA features, or None if no features are processed.

        Raises:
            AssertionError: If `msa_entity_type` is "rna" and `is_homomer_or_monomer` is False.
        """
        if msa_entity_type == "rna":
            assert is_homomer_or_monomer, "RNA MSAs do not pairing"
        pdb_id = bioassembly_dict["pdb_id"]
        sequence_to_features: dict[str, dict[str, Any]] = {}

        for entity_id, sequence in entity_id_to_sequence.items():
            if sequence in sequence_to_features:
                # It is possible that different entity ids correspond to the same sequence
                continue

            if all(
                [
                    asym_id_int not in selected_asym_ids
                    for asym_id_int in entity_to_asym_id_int[entity_id]
                ]
            ):
                # All chains corresponding to this entity are not selected
                continue

            sequence_feat = self.process_single_sequence(
                pdb_name=f"{pdb_id}_{entity_id}",
                sequence=sequence,
                pdb_id=pdb_id,
                is_homomer_or_monomer=is_homomer_or_monomer,
            )
            sequence_feat = convert_monomer_features(sequence_feat)
            sequence_to_features[sequence] = sequence_feat

        all_chain_features = {
            asym_id_int: deepcopy(sequence_to_features[seq])
            for asym_id_int, seq in asym_id_int_to_sequence.items()
            if seq in sequence_to_features
        }
        del sequence_to_features

        if len(all_chain_features) == 0:
            return None
        np_example = merge_all_chain_features(
            pdb_id=pdb_id,
            all_chain_features=all_chain_features,
            asym_to_entity_id=asym_to_entity_id,
            is_homomer_or_monomer=is_homomer_or_monomer,
            merge_method=self.merge_method,
            max_size=self.max_size,
            msa_entity_type=msa_entity_type,
        )
        return np_example


class PROTMSAFeaturizer(BaseMSAFeaturizer):
    def __init__(
        self,
        dataset_name: str = "",
        seq_to_pdb_idx_path: str = "",
        distillation_index_file: str = None,
        indexing_method: str = "sequence",
        pairing_db: Optional[str] = "",
        non_pairing_db: str = "mmseqs_all",
        merge_method: str = "dense_max",
        seq_limits: Optional[dict[str, int]] = {},
        max_size: int = 16384,
        pdb_jackhmmer_dir: str = None,
        pdb_mmseqs_dir: str = None,
        distillation_mmseqs_dir: str = None,
        distillation_uniclust_dir: str = None,
        **kwargs,
    ):
        super().__init__(
            indexing_method=indexing_method,
            merge_method=merge_method,
            seq_limits=seq_limits,
            max_size=max_size,
            **kwargs,
        )
        self.dataset_name = dataset_name
        self.pdb_jackhmmer_dir = pdb_jackhmmer_dir
        self.pdb_mmseqs_dir = pdb_mmseqs_dir
        self.distillation_mmseqs_dir = distillation_mmseqs_dir
        self.distillation_uniclust_dir = distillation_uniclust_dir
        self.pairing_db = pairing_db if len(pairing_db) > 0 else None

        if non_pairing_db == "mmseqs_all":
            self.non_pairing_db = ["uniref100", "mmseqs_other"]
        else:
            self.non_pairing_db = [db_name for db_name in non_pairing_db.split(",")]

        with open(seq_to_pdb_idx_path, "r") as f:
            self.seq_to_pdb_idx = json.load(f)
        # If distillation data is avaiable
        if distillation_index_file is not None:
            with open(distillation_index_file, "r") as f:
                self.distillation_pdb_id_to_msa_dir = json.load(f)
        else:
            self.distillation_pdb_id_to_msa_dir = None

    def get_msa_path(self, db_name: str, sequence: str, pdb_id: str) -> str:
        """
        Get the path of an MSA file

        Args:
            db_name (str): name of genomics database
            sequence (str): input sequence
            pdb_id (str): pdb_id of input sequence

        Returns:
            str: file path
        """

        if self.indexing_method == "pdb_id" and self.distillation_pdb_id_to_msa_dir:
            rel_path = self.distillation_pdb_id_to_msa_dir[pdb_id]

            if db_name == "uniclust30":
                msa_dir_path = opjoin(self.distillation_uniclust_dir, rel_path)
            elif db_name in ["uniref100", "mmseqs_other"]:
                msa_dir_path = opjoin(self.distillation_mmseqs_dir, rel_path)
            else:
                raise ValueError(
                    f"Indexing with {self.indexing_method} is not supported for {db_name}"
                )

            if opexists(msa_path := opjoin(msa_dir_path, f"{db_name}_hits.a3m")):
                return msa_path
            else:
                return opjoin(msa_dir_path, f"{db_name}.a3m")
        else:
            # indexing_method == "sequence"
            pdb_index = self.seq_to_pdb_idx[sequence]
            if db_name in ["uniref100", "mmseqs_other"]:
                return opjoin(
                    self.pdb_mmseqs_dir, str(pdb_index), f"{db_name}_hits.a3m"
                )
            else:
                return opjoin(
                    self.pdb_jackhmmer_dir,
                    f"pdb_on_{db_name}",
                    "results",
                    f"{pdb_index}.a3m",
                )

    def process_single_sequence(
        self,
        pdb_name: str,
        sequence: str,
        pdb_id: str,
        is_homomer_or_monomer: bool,
    ) -> dict[str, np.ndarray]:
        """
        Get basic MSA features for a single sequence.

        Args:
            pdb_name (str): f"{pdb_id}_{entity_id}" of the input entity
            sequence (str): input sequnce
            pdb_id (str): pdb_id of input sequence
            is_homomer_or_monomer (bool): True if the input sequence is a homomer or a monomer

        Returns:
            Dict[str, np.ndarray]: the basic MSA features of the input sequence
        """

        raw_msa_paths, seq_limits = [], []
        for db_name in self.non_pairing_db:
            if opexists(
                path := self.get_msa_path(db_name, sequence, pdb_id)
            ) and path.endswith(".a3m"):
                raw_msa_paths.append(path)
                seq_limits.append(self.seq_limits.get(db_name, SEQ_LIMITS[db_name]))

        # Get sequence and non-pairing msa features
        sequence_features = process_single_sequence(
            pdb_name=pdb_name,
            sequence=sequence,
            raw_msa_paths=raw_msa_paths,
            seq_limits=seq_limits,
            msa_entity_type="prot",
            msa_type="non_pairing",
        )

        # Get pairing msa features
        if not is_homomer_or_monomer:
            # Separately process the MSA needed for pairing
            raw_msa_paths, seq_limits = [], []
            if opexists(
                path := self.get_msa_path(self.pairing_db, sequence, pdb_id)
            ) and path.endswith(".a3m"):
                raw_msa_paths = [
                    path,
                ]
                seq_limits.append(
                    self.seq_limits.get(self.pairing_db, SEQ_LIMITS[self.pairing_db])
                )

            if len(raw_msa_paths) == 0:
                raise ValueError(f"{pdb_name} does not have MSA for pairing")

            all_seq_msa_features = load_and_process_msa(
                pdb_name=pdb_name,
                msa_type="pairing",
                raw_msa_paths=raw_msa_paths,
                seq_limits=seq_limits,
                identifier_func=get_identifier_func(pairing_db=self.pairing_db),
                handle_empty="raise_error",
            )
            sequence_features.update(all_seq_msa_features)

        return sequence_features

    def get_msa_features_for_assembly(
        self,
        bioassembly_dict: Mapping[str, Any],
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
        selected_token_indices: Optional[torch.Tensor],
    ) -> dict[str, np.ndarray]:
        """
        Get MSA features for the bioassembly.

        Args:
            bioassembly_dict (Mapping[str, Any]): the bioassembly dict with sequence, atom_array and token_array.
            entity_to_asym_id_int (Mapping[str, Sequence[int]]): mapping from entity_id to asym_id_int.
            selected_token_indices (torch.Tensor): Cropped token indices.

        Returns:
            Dict[str, np.ndarray]: the basic MSA features of the bioassembly.
        """
        protein_entity_ids = self.get_entity_ids(
            bioassembly_dict, msa_entity_type="prot"
        )
        if len(protein_entity_ids) == 0:
            return None
        (
            selected_asym_ids,
            asym_id_ints,
            asym_to_entity_id,
            asym_id_int_to_sequence,
            entity_id_to_sequence,
        ) = self.get_selected_asym_ids(
            bioassembly_dict=bioassembly_dict,
            entity_to_asym_id_int=entity_to_asym_id_int,
            selected_token_indices=selected_token_indices,
            entity_ids=protein_entity_ids,
        )
        # No pairing_db specified (all proteins are treated as monomers) or only one sequence
        is_homomer_or_monomer = (self.pairing_db is None) or (
            len(
                set(
                    [
                        asym_id_int_to_sequence[asym_id_int]
                        for asym_id_int in selected_asym_ids
                        if asym_id_int in asym_id_ints
                    ]
                )
            )
            == 1
        )
        np_example = self.get_msa_pipeline(
            is_homomer_or_monomer=is_homomer_or_monomer,
            selected_asym_ids=selected_asym_ids,
            asym_to_entity_id=asym_to_entity_id,
            asym_id_int_to_sequence=asym_id_int_to_sequence,
            entity_id_to_sequence=entity_id_to_sequence,
            bioassembly_dict=bioassembly_dict,
            entity_to_asym_id_int=entity_to_asym_id_int,
            msa_entity_type="prot",
        )
        return np_example


class RNAMSAFeaturizer(BaseMSAFeaturizer):
    def __init__(
        self,
        seq_to_pdb_idx_path: str = "",
        indexing_method: str = "sequence",
        merge_method: str = "dense_max",
        seq_limits: Optional[dict[str, int]] = {},
        max_size: int = 16384,
        rna_msa_dir: str = None,
        **kwargs,
    ) -> None:
        super().__init__(
            indexing_method=indexing_method,
            merge_method=merge_method,
            seq_limits=seq_limits,
            max_size=max_size,
            **kwargs,
        )
        # By default, use all the database in paper
        self.rna_msa_dir = rna_msa_dir
        self.non_pairing_db = ["rfam", "rnacentral", "nucleotide"]
        with open(seq_to_pdb_idx_path, "r") as f:
            self.seq_to_pdb_idx = json.load(f)  # it's rna sequence to pdb list

    def get_msa_path(
        self, db_name: str, sequence: str, pdb_id_entity_id: str, reduced: bool = True
    ) -> str:
        """
        Get the path of an RNA MSA file

        Args:
            db_name (str): genetics databases for RNA chains
            sequence (str): input sequence
            pdb_id_entity_id (str): pdb_id_entity_id of input sequence
            reduced (bool): whether reduce the sto files to max 1w

        Returns:
            str: file path
        """
        assert self.indexing_method in [
            "pdb_id_entity_id",
            "sequence",
        ], "use the pdb_id_entity_id or sequence to search msa dir"
        if reduced:
            suffix = "_max_1w"
        else:
            suffix = ""
        if self.indexing_method == "sequence":
            # only the first pdb save the rna msa
            if sequence in self.seq_to_pdb_idx:
                pdb_id_entity_id = self.seq_to_pdb_idx[sequence][0]
            else:
                logger.info(f"{pdb_id_entity_id} not in seq_to_pdb_idx")
                pdb_id_entity_id = "not_exist"

        rel_path = f"{pdb_id_entity_id}/{db_name}.sto"
        msa_dir_path = opjoin(f"{self.rna_msa_dir}{suffix}", rel_path)
        return msa_dir_path

    def process_single_sequence(
        self,
        pdb_name: str,
        sequence: str,
        pdb_id: str,
        is_homomer_or_monomer: bool,
    ) -> dict[str, np.ndarray]:
        """
        Get basic MSA features for a single sequence.

        Args:
            pdb_name (str): f"{pdb_id}_{entity_id}" of the input entity
            sequence (str): input sequnce
            pdb_id (str): pdb_id of input sequence
            is_homomer_or_monomer (bool): True if the input sequence is a homomer or a monomer

        Returns:
            Dict[str, np.ndarray]: the basic MSA features of the input sequence
        """
        raw_msa_paths, seq_limits = [], []
        for db_name in self.non_pairing_db:
            if opexists(
                path := self.get_msa_path(db_name, sequence, pdb_name)
            ) and path.endswith(".sto"):
                raw_msa_paths.append(path)
                seq_limits.append(self.seq_limits.get(db_name, SEQ_LIMITS[db_name]))

        sequence_features = process_single_sequence(
            pdb_name=pdb_name,
            sequence=sequence,
            raw_msa_paths=raw_msa_paths,
            seq_limits=seq_limits,
            msa_entity_type="rna",
            msa_type="non_pairing",
        )

        return sequence_features

    def get_msa_features_for_assembly(
        self,
        bioassembly_dict: Mapping[str, Any],
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
        selected_token_indices: Optional[torch.Tensor],
    ) -> dict[str, np.ndarray]:
        """
        Get MSA features for the bioassembly.

        Args:
            bioassembly_dict (Mapping[str, Any]): the bioassembly dict with sequence, atom_array and token_array.
            entity_to_asym_id_int (Mapping[str, Sequence[int]]): mapping from entity_id to asym_id_int.
            selected_token_indices (torch.Tensor): Cropped token indices.

        Returns:
            Dict[str, np.ndarray]: the basic MSA features of the bioassembly.
        """
        rna_entity_ids = self.get_entity_ids(bioassembly_dict, msa_entity_type="rna")
        if len(rna_entity_ids) == 0:
            return None
        (
            selected_asym_ids,
            asym_id_ints,
            asym_to_entity_id,
            asym_id_int_to_sequence,
            entity_id_to_sequence,
        ) = self.get_selected_asym_ids(
            bioassembly_dict=bioassembly_dict,
            entity_to_asym_id_int=entity_to_asym_id_int,
            selected_token_indices=selected_token_indices,
            entity_ids=rna_entity_ids,
        )
        is_homomer_or_monomer = True
        np_example = self.get_msa_pipeline(
            is_homomer_or_monomer=is_homomer_or_monomer,
            selected_asym_ids=selected_asym_ids,
            asym_to_entity_id=asym_to_entity_id,
            asym_id_int_to_sequence=asym_id_int_to_sequence,
            entity_id_to_sequence=entity_id_to_sequence,
            bioassembly_dict=bioassembly_dict,
            entity_to_asym_id_int=entity_to_asym_id_int,
            msa_entity_type="rna",
        )
        return np_example


class MSAFeaturizer:
    def __init__(
        self,
        prot_msa_args: dict = {},
        rna_msa_args: dict = {},
        enable_rna_msa: bool = False,
    ):
        self.prot_msa_featurizer = PROTMSAFeaturizer(**prot_msa_args)
        self.enable_rna_msa = enable_rna_msa
        if self.enable_rna_msa:
            self.rna_msa_featurizer = RNAMSAFeaturizer(**rna_msa_args)

    def __call__(
        self,
        bioassembly_dict: dict[str, Any],
        selected_indices: np.ndarray,
        entity_to_asym_id_int: Mapping[str, int],
    ) -> Optional[dict[str, np.ndarray]]:
        """
        Processes the bioassembly dictionary to generate MSA features for both protein and RNA entities, if enabled.

        Args:
            bioassembly_dict (dict[str, Any]): The bioassembly dictionary containing entity information.
            selected_indices (np.ndarray): Indices of selected tokens.
            entity_to_asym_id_int (Mapping[str, int]): Mapping from entity ID to asym ID integers.

        Returns:
            Optional[dict[str, np.ndarray]]: A dictionary containing the merged MSA features for the bioassembly, or None if no features are generated.
        """
        prot_msa_feats = self.prot_msa_featurizer.get_msa_features_for_assembly(
            bioassembly_dict=bioassembly_dict,
            entity_to_asym_id_int=entity_to_asym_id_int,
            selected_token_indices=selected_indices,
        )
        if self.enable_rna_msa:
            rna_msa_feats = self.rna_msa_featurizer.get_msa_features_for_assembly(
                bioassembly_dict=bioassembly_dict,
                entity_to_asym_id_int=entity_to_asym_id_int,
                selected_token_indices=selected_indices,
            )
        else:
            rna_msa_feats = None
        np_chains_list = []
        if prot_msa_feats is not None:
            np_chains_list.append(prot_msa_feats)
        if rna_msa_feats is not None:
            np_chains_list.append(rna_msa_feats)
        if len(np_chains_list) == 0:
            return None

        msa_feats = merge_features_from_prot_rna(np_chains_list)
        msa_feats = self.tokenize(
            msa_feats=msa_feats,
            token_array=bioassembly_dict["token_array"],
            atom_array=bioassembly_dict["atom_array"],
        )

        return msa_feats

    def tokenize(
        self,
        msa_feats: Mapping[str, np.ndarray],
        token_array: TokenArray,
        atom_array: AtomArray,
    ) -> dict[str, np.ndarray]:
        """
        Tokenize raw MSA features.

        Args:
            msa_feats (Dict[str, np.ndarray]): raw MSA features.
            token_array (TokenArray): token array of this bioassembly
            atom_array (AtomArray): atom array of this bioassembly

        Returns:
            Dict[str, np.ndarray]: the tokenized MSA features of the bioassembly.
        """
        msa_feats = tokenize_msa(
            msa_feats=msa_feats, token_array=token_array, atom_array=atom_array
        )
        # Add to tracking for msa analysis
        msa_feats.update(
            {
                "prot_pair_num_alignments": msa_feats.get(
                    "prot_pair_num_alignments", np.asarray(0, dtype=np.int32)
                ),
                "prot_unpair_num_alignments": msa_feats.get(
                    "prot_unpair_num_alignments", np.asarray(0, dtype=np.int32)
                ),
                "rna_pair_num_alignments": msa_feats.get(
                    "rna_pair_num_alignments", np.asarray(0, dtype=np.int32)
                ),
                "rna_unpair_num_alignments": msa_feats.get(
                    "rna_unpair_num_alignments", np.asarray(0, dtype=np.int32)
                ),
            }
        )
        return {
            k: v
            for (k, v) in msa_feats.items()
            if k
            in ["msa", "has_deletion", "deletion_value", "deletion_mean", "profile"]
            + [
                "prot_pair_num_alignments",
                "prot_unpair_num_alignments",
                "rna_pair_num_alignments",
                "rna_unpair_num_alignments",
            ]
        }


# Common function for train and inference
def process_single_sequence(
    pdb_name: str,
    sequence: str,
    raw_msa_paths: Optional[list[str]],
    seq_limits: Optional[list[str]],
    msa_entity_type: str = "prot",
    msa_type: str = "non_pairing",
) -> FeatureDict:
    """
    Processes a single sequence to generate sequence and MSA features.

    Args:
        pdb_name (str): The name of the PDB entry.
        sequence (str): The input sequence.
        raw_msa_paths (Optional[list[str]]): List of paths to raw MSA files.
        seq_limits (Optional[list[str]]): List of sequence limits for different databases.
        msa_entity_type (str): The type of MSA entity, either "prot" or "rna". Defaults to "prot".
        msa_type (str): The type of MSA, either "non_pairing" or "pairing". Defaults to "non_pairing".

    Returns:
        FeatureDict: A dictionary containing the sequence and MSA features.

    Raises:
        AssertionError: If `msa_entity_type` is not "prot" or "rna".
    """
    assert msa_entity_type in ["prot", "rna"]
    num_res = len(sequence)

    if msa_entity_type == "prot":
        sequence_features = make_sequence_features(
            sequence=sequence,
            num_res=num_res,
        )
    elif msa_entity_type == "rna":
        sequence_features = make_sequence_features(
            sequence=sequence,
            num_res=num_res,
            mapping=rna_order_with_x,
            x_token="N",
        )

    msa_features = load_and_process_msa(
        pdb_name=pdb_name,
        msa_type=msa_type,
        raw_msa_paths=raw_msa_paths,
        seq_limits=seq_limits,
        input_sequence=sequence,
        msa_entity_type=msa_entity_type,
    )
    sequence_features.update(msa_features)
    return sequence_features


# Common function for train and inference
def tokenize_msa(
    msa_feats: Mapping[str, np.ndarray],
    token_array: TokenArray,
    atom_array: AtomArray,
) -> dict[str, np.ndarray]:
    """
    Tokenize raw MSA features.

    Args:
        msa_feats (Dict[str, np.ndarray]): raw MSA features.
        token_array (TokenArray): token array of this bioassembly
        atom_array (AtomArray): atom array of this bioassembly

    Returns:
        Dict[str, np.ndarray]: the tokenized MSA features of the bioassembly.
    """
    token_center_atom_idxs = token_array.get_annotation("centre_atom_index")
    # res_id: (asym_id, residue_index)
    # msa_idx refers to the column number of a residue in the msa array
    res_id_2_msa_idx = {
        (msa_feats["asym_id"][idx], msa_feats["residue_index"][idx]): idx
        for idx in range(msa_feats["msa"].shape[1])
    }

    restypes = []
    col_idxs_in_msa = []
    col_idxs_in_new_msa = []
    for token_idx, center_atom_idx in enumerate(token_center_atom_idxs):
        restypes.append(STD_RESIDUES[atom_array.cano_seq_resname[center_atom_idx]])
        if (
            res_id := (
                atom_array[center_atom_idx].asym_id_int,
                atom_array[center_atom_idx].res_id,
            )
        ) in res_id_2_msa_idx:
            col_idxs_in_msa.append(res_id_2_msa_idx[res_id])
            col_idxs_in_new_msa.append(token_idx)

    num_msa_seq, _ = msa_feats["msa"].shape
    num_tokens = len(token_center_atom_idxs)

    restypes = np.array(restypes)
    col_idxs_in_new_msa = np.array(col_idxs_in_new_msa)
    col_idxs_in_msa = np.array(col_idxs_in_msa)

    # msa
    # For non-amino acid tokens, copy the token itself
    feat_name = "msa"
    new_feat = np.repeat(restypes[None, ...], num_msa_seq, axis=0)
    new_feat[:, col_idxs_in_new_msa] = msa_feats[feat_name][:, col_idxs_in_msa]
    msa_feats[feat_name] = new_feat

    # has_deletion, deletion_value
    # Assign 0 to non-amino acid tokens
    for feat_name in ["has_deletion", "deletion_value"]:
        new_feat = np.zeros((num_msa_seq, num_tokens), dtype=msa_feats[feat_name].dtype)
        new_feat[:, col_idxs_in_new_msa] = msa_feats[feat_name][:, col_idxs_in_msa]
        msa_feats[feat_name] = new_feat

    # deletion_mean
    # Assign 0 to non-amino acid tokens
    feat_name = "deletion_mean"
    new_feat = np.zeros((num_tokens,))
    new_feat[col_idxs_in_new_msa] = msa_feats[feat_name][col_idxs_in_msa]
    msa_feats[feat_name] = new_feat

    # profile
    # Assign one-hot enbedding (one-hot distribution) to non-amino acid tokens corresponding to restype
    feat_name = "profile"
    new_feat = np.zeros((num_tokens, 32))
    new_feat[np.arange(num_tokens), restypes] = 1
    new_feat[col_idxs_in_new_msa, :] = msa_feats[feat_name][col_idxs_in_msa, :]
    msa_feats[feat_name] = new_feat
    return msa_feats


# Common function for train and inference
def merge_all_chain_features(
    pdb_id: str,
    all_chain_features: dict[str, FeatureDict],
    asym_to_entity_id: dict,
    is_homomer_or_monomer: bool = False,
    merge_method: str = "dense_max",
    max_size: int = 16384,
    msa_entity_type: str = "prot",
) -> dict[str, np.ndarray]:
    """
    Merges features from all chains in the bioassembly.

    Args:
        pdb_id (str): The PDB ID of the bioassembly.
        all_chain_features (dict[str, FeatureDict]): Features for each chain in the bioassembly.
        asym_to_entity_id (dict): Mapping from asym ID to entity ID.
        is_homomer_or_monomer (bool): Indicates if the bioassembly is a homomer or monomer. Defaults to False.
        merge_method (str): Method used for merging features. Defaults to "dense_max".
        max_size (int): Maximum size of the MSA. Defaults to 16384.
        msa_entity_type (str): Type of MSA entity, either "prot" or "rna". Defaults to "prot".

    Returns:
        dict[str, np.ndarray]: Merged features for the bioassembly.
    """
    all_chain_features = add_assembly_features(
        pdb_id,
        all_chain_features,
        asym_to_entity_id=asym_to_entity_id,
    )
    if msa_entity_type == "rna":
        np_example = rna_merge(
            is_homomer_or_monomer=is_homomer_or_monomer,
            all_chain_features=all_chain_features,
            merge_method=merge_method,
            msa_crop_size=max_size,
        )
    elif msa_entity_type == "prot":
        np_example = pair_and_merge(
            is_homomer_or_monomer=is_homomer_or_monomer,
            all_chain_features=all_chain_features,
            merge_method=merge_method,
            msa_crop_size=max_size,
        )
    np_example = clip_msa(np_example, max_num_msa=max_size)
    return np_example


class InferenceMSAFeaturizer(object):
    # Now we only support protein msa in inference

    @staticmethod
    def process_prot_single_sequence(
        sequence: str,
        description: str,
        is_homomer_or_monomer: bool,
        msa_dir: Union[str, None],
        pairing_db: str,
    ) -> FeatureDict:
        """
        Processes a single protein sequence to generate sequence and MSA features.

        Args:
            sequence (str): The input protein sequence.
            description (str): Description of the sequence, typically the PDB name.
            is_homomer_or_monomer (bool): Indicates if the sequence is a homomer or monomer.
            msa_dir (Union[str, None]): Directory containing the MSA files, or None if no pre-computed MSA is provided.
            pairing_db (str): Database used for pairing.

        Returns:
            FeatureDict: A dictionary containing the sequence and MSA features.

        Raises:
            AssertionError: If the pairing MSA file does not exist when `is_homomer_or_monomer` is False.
        """
        # For non-pairing MSA
        if msa_dir is None:
            # No pre-computed MSA was provided, and the MSA search failed
            raw_msa_paths = []
        else:
            raw_msa_paths = [opjoin(msa_dir, "non_pairing.a3m")]
        pdb_name = description

        sequence_features = process_single_sequence(
            pdb_name=pdb_name,
            sequence=sequence,
            raw_msa_paths=raw_msa_paths,
            seq_limits=[-1],
            msa_entity_type="prot",
            msa_type="non_pairing",
        )
        if not is_homomer_or_monomer:
            # Separately process the pairing MSA
            assert opexists(
                raw_msa_path := opjoin(msa_dir, "pairing.a3m")
            ), f"No pairing-MSA of {pdb_name} (please check {raw_msa_path})"

            all_seq_msa_features = load_and_process_msa(
                pdb_name=pdb_name,
                msa_type="pairing",
                raw_msa_paths=[raw_msa_path],
                seq_limits=[-1],
                identifier_func=get_identifier_func(
                    pairing_db=pairing_db,
                ),
                handle_empty="raise_error",
            )
            sequence_features.update(all_seq_msa_features)

        return sequence_features

    @staticmethod
    def get_inference_prot_msa_features_for_assembly(
        bioassembly: Sequence[Mapping[str, Mapping[str, Any]]],
        entity_to_asym_id: Mapping[str, set[int]],
    ) -> FeatureDict:
        """
        Processes the bioassembly to generate MSA features for protein entities in inference mode.

        Args:
            bioassembly (Sequence[Mapping[str, Mapping[str, Any]]]): The bioassembly containing entity information.
            entity_to_asym_id (Mapping[str, set[int]]): Mapping from entity ID to asym ID integers.

        Returns:
            FeatureDict: A dictionary containing the MSA features for the protein entities.

        Raises:
            AssertionError: If the provided precomputed MSA path does not exist.
        """
        entity_to_asym_id_int = dict(entity_to_asym_id)
        asym_to_entity_id = {}
        entity_id_to_sequence = {}
        # In inference mode, the keys in bioassembly is different from training
        # Only contains protein entity, many-to-one mapping
        entity_id_to_sequence = {}
        for i, entity_info_wrapper in enumerate(bioassembly):
            entity_id = str(i + 1)
            entity_type = list(entity_info_wrapper.keys())[0]
            entity_info = entity_info_wrapper[entity_type]

            if entity_type == PROT_TYPE_NAME:
                # Update entity_id_to_sequence
                entity_id_to_sequence[entity_id] = entity_info["sequence"]

                # Update asym_to_entity_id
                for asym_id_int in entity_to_asym_id_int[entity_id]:
                    asym_to_entity_id[asym_id_int] = entity_id
        if len(entity_id_to_sequence) == 0:
            # No protein entity
            return None
        is_homomer_or_monomer = (
            len(set(entity_id_to_sequence.values())) == 1
        )  # Only one protein sequence
        sequence_to_entity = defaultdict(list)
        for entity_id, seq in entity_id_to_sequence.items():
            sequence_to_entity[seq].append(entity_id)

        sequence_to_features: dict[str, dict[str, Any]] = {}
        msa_sequences = {}
        msa_dirs = {}
        for idx, (sequence, entity_id_list) in enumerate(sequence_to_entity.items()):
            msa_info = bioassembly[int(entity_id_list[0]) - 1][PROT_TYPE_NAME]["msa"]
            msa_dir = msa_info.get("precomputed_msa_dir", None)
            if msa_dir is not None:
                assert opexists(
                    msa_dir
                ), f"The provided precomputed MSA path of entities {entity_id_list} does not exists: \n{msa_dir}"
                msa_dirs[idx] = msa_dir
            else:
                pairing_db_fpath = msa_info.get("pairing_db_fpath", None)
                non_pairing_db_fpath = msa_info.get("non_pairing_db_fpath", None)
                assert (
                    pairing_db_fpath is not None
                ), "Path of pairing MSA database is not given."
                assert (
                    non_pairing_db_fpath is not None
                ), "Path of non-pairing MSA database is not given."
                assert msa_info["pairing_db"] in ["uniprot", "", None], (
                    f"Using {msa_info['pairing_db']} as the source for MSA pairing "
                    f"is not supported in online MSA searching."
                )

                msa_info["pairing_db"] = "uniprot"
                msa_sequences[idx] = (sequence, pairing_db_fpath, non_pairing_db_fpath)
        if len(msa_sequences) > 0:
            msa_dirs.update(msa_parallel(msa_sequences))

        for idx, (sequence, entity_id_list) in enumerate(sequence_to_entity.items()):

            if len(entity_id_list) > 1:
                logger.info(
                    f"Entities {entity_id_list} correspond to the same sequence."
                )
            msa_info = bioassembly[int(entity_id_list[0]) - 1][PROT_TYPE_NAME]["msa"]
            msa_dir = msa_dirs[idx]

            description = f"entity_{'_'.join(map(str, entity_id_list))}"
            sequence_feat = InferenceMSAFeaturizer.process_prot_single_sequence(
                sequence=sequence,
                description=description,
                is_homomer_or_monomer=is_homomer_or_monomer,
                msa_dir=msa_dir,
                pairing_db=msa_info["pairing_db"],
            )
            sequence_feat = convert_monomer_features(sequence_feat)
            sequence_to_features[sequence] = sequence_feat
            if msa_dir and opexists(msa_dir) and idx in msa_sequences.keys():
                if (msa_save_dir := msa_info.get("msa_save_dir", None)) is not None:
                    if opexists(dst_dir := opjoin(msa_save_dir, str(idx + 1))):
                        shutil.rmtree(dst_dir)
                    shutil.copytree(msa_dir, dst_dir)
                    for fname in os.listdir(dst_dir):
                        if not fname.endswith(".a3m"):
                            os.remove(opjoin(dst_dir, fname))
                else:
                    shutil.rmtree(msa_dir)

        all_chain_features = {
            asym_id_int: deepcopy(
                sequence_to_features[entity_id_to_sequence[entity_id]]
            )
            for asym_id_int, entity_id in asym_to_entity_id.items()
            if seq in sequence_to_features
        }
        if len(all_chain_features) == 0:
            return None

        np_example = merge_all_chain_features(
            pdb_id="test_assembly",
            all_chain_features=all_chain_features,
            asym_to_entity_id=asym_to_entity_id,
            is_homomer_or_monomer=is_homomer_or_monomer,
            merge_method="dense_max",
            max_size=MSA_MAX_SIZE,
            msa_entity_type="prot",
        )

        return np_example

    def make_msa_feature(
        bioassembly: Sequence[Mapping[str, Mapping[str, Any]]],
        entity_to_asym_id: Mapping[str, Sequence[str]],
        token_array: TokenArray,
        atom_array: AtomArray,
    ) -> Optional[dict[str, np.ndarray]]:
        """
        Processes the bioassembly to generate MSA features for protein entities in inference mode and tokenizes the features.

        Args:
            bioassembly (Sequence[Mapping[str, Mapping[str, Any]]]): The bioassembly containing entity information.
            entity_to_asym_id (Mapping[str, Sequence[str]]): Mapping from entity ID to asym ID strings.
            token_array (TokenArray): Token array of the bioassembly.
            atom_array (AtomArray): Atom array of the bioassembly.

        Returns:
            Optional[dict[str, np.ndarray]]: A dictionary containing the tokenized MSA features for the protein entities,
                or an empty dictionary if no features are generated.
        """
        msa_feats = InferenceMSAFeaturizer.get_inference_prot_msa_features_for_assembly(
            bioassembly=bioassembly,
            entity_to_asym_id=entity_to_asym_id,
        )

        if msa_feats is None:
            return {}

        msa_feats = tokenize_msa(
            msa_feats=msa_feats,
            token_array=token_array,
            atom_array=atom_array,
        )
        return {
            k: v
            for (k, v) in msa_feats.items()
            if k
            in ["msa", "has_deletion", "deletion_value", "deletion_mean", "profile"]
        }
