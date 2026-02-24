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
import logging
import os
import time
import traceback
import warnings
from typing import Any, Mapping

import torch
from biotite.structure import AtomArray
from torch.utils.data import DataLoader, Dataset, DistributedSampler

from protenix.data.esm.esm_featurizer import ESMFeaturizer
from protenix.data.inference.json_to_feature import SampleDictToFeatures
from protenix.data.msa.msa_featurizer import InferenceMSAFeaturizer
from protenix.data.template.template_featurizer import InferenceTemplateFeaturizer
from protenix.data.template.template_utils import TemplateHitFeaturizer
from protenix.data.utils import data_type_transform, make_dummy_feature
from protenix.utils.distributed import DIST_WRAPPER
from protenix.utils.torch_utils import collate_fn_identity, dict_to_tensor

logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore", module="biotite")


def get_inference_dataloader(configs: Any) -> DataLoader:
    """
    Creates and returns a DataLoader for inference using the InferenceDataset.

    Args:
        configs: A configuration object containing the necessary parameters for the DataLoader.

    Returns:
        A DataLoader object configured for inference.
    """
    inference_dataset = InferenceDataset(
        configs=configs,
    )
    sampler = DistributedSampler(
        dataset=inference_dataset,
        num_replicas=DIST_WRAPPER.world_size,
        rank=DIST_WRAPPER.rank,
        shuffle=False,
    )
    dataloader = DataLoader(
        dataset=inference_dataset,
        batch_size=1,
        sampler=sampler,
        collate_fn=collate_fn_identity,
        num_workers=configs.num_workers,
    )
    return dataloader


class InferenceDataset(Dataset):
    def __init__(
        self,
        configs,
    ) -> None:
        self.configs = configs

        self.input_json_path = configs.input_json_path
        self.dump_dir = configs.dump_dir
        self.use_msa = configs.use_msa
        self.msa_pair_as_unpair = configs.get("msa_pair_as_unpair", True)
        self.use_rna_msa = configs.get("use_rna_msa", True)
        self.use_template = configs.get("use_template", True)
        with open(self.input_json_path, "r") as f:
            self.inputs = json.load(f)
        json_task_name = os.path.basename(self.input_json_path).split(".")[0]
        if self.use_template:
            template_mmcif_dir = configs.data.template.prot_template_mmcif_dir
            fetch_remote = configs.data.template.get("fetch_remote", False)
            if not fetch_remote:
                assert template_mmcif_dir is not None and os.path.exists(
                    template_mmcif_dir
                ), (
                    "Inference with template depends on the mmcif directory.\n"
                    "The mmcif directory containing cif files should be placed under $PROTENIX_ROOT_DIR/mmcif.\n"
                    "You can download it from PDB https://www.wwpdb.org/ftp/pdb-ftp-sites or\n"
                    "refer to scripts/database/download_protenix_data.sh to download inference dependency files, "
                    "set use_template=false for inference, or set data.template.fetch_remote=true "
                    "to download mmCIF files on demand from PDBe."
                )
            else:
                if template_mmcif_dir:
                    os.makedirs(template_mmcif_dir, exist_ok=True)
            self.online_template_featurizer = TemplateHitFeaturizer(
                mmcif_dir=configs.data.template.prot_template_mmcif_dir,
                template_cache_dir=configs.data.template.prot_template_cache_dir,
                max_hits=4,
                kalign_binary_path=configs.data.template.kalign_binary_path,
                max_template_date="2021-09-30",
                release_dates_path=configs.data.template.release_dates_path,
                obsolete_pdbs_path=configs.data.template.obsolete_pdbs_path,
                _shuffle_top_k_prefiltered=None,
                _max_template_candidates_num=20,
                fetch_remote=fetch_remote,
            )
        else:
            self.online_template_featurizer = None
        esm_info = configs.get("esm", {})
        configs.esm.embedding_dir = f"./esm_embeddings/{configs.esm.model_name}"
        configs.esm.sequence_fpath = (
            f"./esm_embeddings/{json_task_name}_prot_sequences.csv"
        )
        self.esm_enable = esm_info.get("enable", False)
        if self.esm_enable:
            os.makedirs(configs.esm.embedding_dir, exist_ok=True)
            os.makedirs(os.path.dirname(configs.esm.sequence_fpath), exist_ok=True)
            ESMFeaturizer.precompute_esm_embedding(
                self.inputs,
                configs.esm.model_name,
                configs.esm.embedding_dir,
                configs.esm.sequence_fpath,
                configs.load_checkpoint_dir,
            )
            self.esm_featurizer = ESMFeaturizer(
                embedding_dir=esm_info.embedding_dir,
                sequence_fpath=esm_info.sequence_fpath,
                embedding_dim=esm_info.embedding_dim,
                error_dir="./esm_embeddings/",
            )

    def process_one(
        self,
        single_sample_dict: Mapping[str, Any],
    ) -> tuple[dict[str, torch.Tensor], AtomArray, dict[str, float]]:
        """
        Processes a single sample from the input JSON to generate features and statistics.

        Args:
            single_sample_dict: A dictionary containing the sample data.

        Returns:
            A tuple containing:
                - A dictionary of features.
                - An AtomArray object.
                - A dictionary of time tracking statistics.
        """
        # general features
        t0 = time.time()
        sample2feat = SampleDictToFeatures(
            single_sample_dict,
        )
        features_dict, atom_array, token_array = sample2feat.get_feature_dict()
        features_dict["distogram_rep_atom_mask"] = torch.Tensor(
            atom_array.distogram_rep_atom_mask
        ).long()
        entity_poly_type_and_seqs = (
            sample2feat.entity_poly_type_and_seqs
        )  # we include ligand as well
        t1 = time.time()
        msa_features = (
            InferenceMSAFeaturizer.make_msa_feature(
                bioassembly=single_sample_dict["sequences"],
                atom_array=atom_array,
                msa_pair_as_unpair=self.msa_pair_as_unpair,
                use_rna_msa=self.use_rna_msa,
            )
            if self.use_msa
            else {}
        )
        template_features = InferenceTemplateFeaturizer.make_template_feature(
            bioassembly=single_sample_dict["sequences"],
            atom_array=atom_array,
            use_template=self.use_template,
            online_template_featurizer=self.online_template_featurizer,
        )
        # Esm features
        if self.esm_enable:
            x_esm = self.esm_featurizer(
                token_array=token_array,
                atom_array=atom_array,
                bioassembly_dict=single_sample_dict,
                inference_mode=True,
            )
            features_dict["esm_token_embedding"] = x_esm

        # Make dummy features for not implemented features
        dummy_feats = []
        if len(template_features) == 0:
            dummy_feats.append("template")
        else:
            template_features = dict_to_tensor(template_features)
            features_dict.update(template_features)
        if len(msa_features) == 0:
            dummy_feats.append("msa")
        else:
            msa_features = dict_to_tensor(msa_features)
            features_dict.update(msa_features)
        features_dict = make_dummy_feature(
            features_dict=features_dict,
            dummy_feats=dummy_feats,
        )

        # Transform to right data type
        feat = data_type_transform(feat_or_label_dict=features_dict)

        t2 = time.time()

        data = {}
        data["input_feature_dict"] = feat
        # Add dimension related items
        N_token = feat["token_index"].shape[0]
        N_atom = feat["atom_to_token_idx"].shape[0]
        N_msa = feat["msa"].shape[0]
        stats = {}
        for mol_type in ["ligand", "protein", "dna", "rna"]:
            mol_type_mask = feat[f"is_{mol_type}"].bool()
            stats[f"{mol_type}/atom"] = int(mol_type_mask.sum(dim=-1).item())
            stats[f"{mol_type}/token"] = len(
                torch.unique(feat["atom_to_token_idx"][mol_type_mask])
            )
        N_asym = len(torch.unique(data["input_feature_dict"]["asym_id"]))
        data.update(
            {
                "N_asym": torch.tensor([N_asym]),
                "N_token": torch.tensor([N_token]),
                "N_atom": torch.tensor([N_atom]),
                "N_msa": torch.tensor([N_msa]),
            }
        )

        def formatted_key(key):
            type_, unit = key.split("/")
            if type_ == "protein":
                type_ = "prot"
            elif type_ == "ligand":
                type_ = "lig"
            else:
                pass
            return f"N_{type_}_{unit}"

        data.update(
            {
                formatted_key(k): torch.tensor([stats[k]])
                for k in [
                    "protein/atom",
                    "ligand/atom",
                    "dna/atom",
                    "rna/atom",
                    "protein/token",
                    "ligand/token",
                    "dna/token",
                    "rna/token",
                ]
            }
        )
        data.update({"entity_poly_type": entity_poly_type_and_seqs["entity_poly_type"]})
        t3 = time.time()
        time_tracker = {
            "crop": t1 - t0,
            "featurizer": t2 - t1,
            "added_feature": t3 - t2,
        }

        return data, atom_array, time_tracker

    def __len__(self) -> int:
        return len(self.inputs)

    def __getitem__(self, index: int) -> tuple[dict[str, torch.Tensor], AtomArray, str]:
        try:
            single_sample_dict = self.inputs[index]
            sample_name = single_sample_dict["name"]
            logger.info(f"Featurizing {sample_name}...")

            data, atom_array, _ = self.process_one(
                single_sample_dict=single_sample_dict
            )
            error_message = ""
        except Exception as e:
            data, atom_array = {}, None
            error_message = f"{e}:\n{traceback.format_exc()}"
        data["sample_name"] = single_sample_dict["name"]
        data["sample_index"] = index
        return data, atom_array, error_message
