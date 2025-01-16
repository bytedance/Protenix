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

# pylint: disable=C0114,C0301
import os
from copy import deepcopy

from protenix.config.extend_types import GlobalConfigValue, ListValue

default_test_configs = {
    "sampler_configs": {
        "sampler_type": "uniform",
    },
    "cropping_configs": {
        "method_weights": [
            0.0,  # ContiguousCropping
            0.0,  # SpatialCropping
            1.0,  # SpatialInterfaceCropping
        ],
        "crop_size": -1,
    },
    "lig_atom_rename": GlobalConfigValue("test_lig_atom_rename"),
    "shuffle_mols": GlobalConfigValue("test_shuffle_mols"),
    "shuffle_sym_ids": GlobalConfigValue("test_shuffle_sym_ids"),
}

default_weighted_pdb_configs = {
    "sampler_configs": {
        "sampler_type": "weighted",
        "beta_dict": {
            "chain": 0.5,
            "interface": 1,
        },
        "alpha_dict": {
            "prot": 3,
            "nuc": 3,
            "ligand": 1,
        },
        "force_recompute_weight": True,
    },
    "cropping_configs": {
        "method_weights": ListValue([0.2, 0.4, 0.4]),
        "crop_size": GlobalConfigValue("train_crop_size"),
    },
    "sample_weight": 0.5,
    "limits": -1,
    "lig_atom_rename": GlobalConfigValue("train_lig_atom_rename"),
    "shuffle_mols": GlobalConfigValue("train_shuffle_mols"),
    "shuffle_sym_ids": GlobalConfigValue("train_shuffle_sym_ids"),
}

DATA_ROOT_DIR = "/af3-dev/release_data/"

# Use CCD cache created by scripts/gen_ccd_cache.py priority. (without date in filename)
# See: docs/prepare_data.md
CCD_COMPONENTS_FILE_PATH = os.path.join(DATA_ROOT_DIR, "components.cif")
CCD_COMPONENTS_RDKIT_MOL_FILE_PATH = os.path.join(
    DATA_ROOT_DIR, "components.cif.rdkit_mol.pkl"
)

if (not os.path.exists(CCD_COMPONENTS_FILE_PATH)) or (
    not os.path.exists(CCD_COMPONENTS_RDKIT_MOL_FILE_PATH)
):
    CCD_COMPONENTS_FILE_PATH = os.path.join(DATA_ROOT_DIR, "components.v20240608.cif")
    CCD_COMPONENTS_RDKIT_MOL_FILE_PATH = os.path.join(
        DATA_ROOT_DIR, "components.v20240608.cif.rdkit_mol.pkl"
    )


# This is a patch in inference stage for users that do not have root permission.
# If you run
# ```
# bash inference_demo.sh
# ```
# or
# ```
# protenix predict --input examples/example.json --out_dir  ./output
# ````
# The checkpoint and the data cache will be downloaded to the current code directory.
if (not os.path.exists(CCD_COMPONENTS_FILE_PATH)) or (
    not os.path.exists(CCD_COMPONENTS_RDKIT_MOL_FILE_PATH)
):
    print("Try to find the ccd cache data in the code directory for inference.")
    current_file_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_file_path)
    code_directory = os.path.dirname(current_directory)

    data_cache_dir = os.path.join(code_directory, "release_data/ccd_cache")
    CCD_COMPONENTS_FILE_PATH = os.path.join(data_cache_dir, "components.cif")
    CCD_COMPONENTS_RDKIT_MOL_FILE_PATH = os.path.join(
        data_cache_dir, "components.cif.rdkit_mol.pkl"
    )
    if (not os.path.exists(CCD_COMPONENTS_FILE_PATH)) or (
        not os.path.exists(CCD_COMPONENTS_RDKIT_MOL_FILE_PATH)
    ):

        CCD_COMPONENTS_FILE_PATH = os.path.join(
            data_cache_dir, "components.v20240608.cif"
        )
        CCD_COMPONENTS_RDKIT_MOL_FILE_PATH = os.path.join(
            data_cache_dir, "components.v20240608.cif.rdkit_mol.pkl"
        )

data_configs = {
    "num_dl_workers": 16,
    "epoch_size": 10000,
    "train_ref_pos_augment": True,
    "test_ref_pos_augment": True,
    "train_sets": ListValue(["weightedPDB_before2109_wopb_nometalc_0925"]),
    "train_sampler": {
        "train_sample_weights": ListValue([1.0]),
        "sampler_type": "weighted",
    },
    "test_sets": ListValue(["recentPDB_1536_sample384_0925"]),
    "weightedPDB_before2109_wopb_nometalc_0925": {
        "base_info": {
            "mmcif_dir": os.path.join(DATA_ROOT_DIR, "mmcif"),
            "bioassembly_dict_dir": os.path.join(DATA_ROOT_DIR, "mmcif_bioassembly"),
            "indices_fpath": os.path.join(
                DATA_ROOT_DIR,
                "indices/weightedPDB_indices_before_2021-09-30_wo_posebusters_resolution_below_9.csv.gz",
            ),
            "pdb_list": "",
            "random_sample_if_failed": True,
            "max_n_token": -1,  # can be used for removing data with too many tokens.
            "use_reference_chains_only": False,
            "exclusion": {  # do not sample the data based on ions.
                "mol_1_type": ListValue(["ions"]),
                "mol_2_type": ListValue(["ions"]),
            },
        },
        **deepcopy(default_weighted_pdb_configs),
    },
    "recentPDB_1536_sample384_0925": {
        "base_info": {
            "mmcif_dir": os.path.join(DATA_ROOT_DIR, "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                DATA_ROOT_DIR, "recentPDB_bioassembly"
            ),
            "indices_fpath": os.path.join(
                DATA_ROOT_DIR, "indices/recentPDB_low_homology_maxtoken1536.csv"
            ),
            "pdb_list": os.path.join(
                DATA_ROOT_DIR,
                "indices/recentPDB_low_homology_maxtoken1024_sample384_pdb_id.txt",
            ),
            "max_n_token": GlobalConfigValue("test_max_n_token"),  # filter data
            "sort_by_n_token": False,
            "group_by_pdb_id": True,
            "find_eval_chain_interface": True,
        },
        **deepcopy(default_test_configs),
    },
    "posebusters_0925": {
        "base_info": {
            "mmcif_dir": os.path.join(DATA_ROOT_DIR, "posebusters_mmcif"),
            "bioassembly_dict_dir": os.path.join(
                DATA_ROOT_DIR, "posebusters_bioassembly"
            ),
            "indices_fpath": os.path.join(
                DATA_ROOT_DIR, "indices/posebusters_indices_mainchain_interface.csv"
            ),
            "pdb_list": "",
            "find_pocket": True,
            "find_all_pockets": False,
            "max_n_token": GlobalConfigValue("test_max_n_token"),  # filter data
        },
        **deepcopy(default_test_configs),
    },
    "msa": {
        "enable": True,
        "enable_rna_msa": False,
        "prot": {
            "pairing_db": "uniref100",
            "non_pairing_db": "mmseqs_other",
            "pdb_mmseqs_dir": os.path.join(DATA_ROOT_DIR, "mmcif_msa"),
            "seq_to_pdb_idx_path": os.path.join(DATA_ROOT_DIR, "seq_to_pdb_index.json"),
            "indexing_method": "sequence",
        },
        "rna": {
            "seq_to_pdb_idx_path": "",
            "rna_msa_dir": "",
            "indexing_method": "sequence",
        },
        "strategy": "random",
        "merge_method": "dense_max",
        "min_size": {
            "train": 1,
            "test": 2048,
        },
        "max_size": {
            "train": 16384,
            "test": 16384,
        },
        "sample_cutoff": {
            "train": 2048,
            "test": 2048,
        },
    },
    "template": {
        "enable": False,
    },
    "ccd_components_file": CCD_COMPONENTS_FILE_PATH,
    "ccd_components_rdkit_mol_file": CCD_COMPONENTS_RDKIT_MOL_FILE_PATH,
}
