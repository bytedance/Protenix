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
from pathlib import Path

from protenix.config.extend_types import GlobalConfigValue, ListValue

PROTENIX_ROOT_DIR = os.environ.get("PROTENIX_ROOT_DIR", str(Path.home()))

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
    "constraint": {
        "enable": False,
        "fix_seed": False,  # True means use use the same contact in each evaluation.
    },
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
    # If enabled, the training settings for different constraint types,
    # providing the model a certain proportion of constraints
    # that meet specific conditions.
    "constraint": {
        "enable": False,
        "fix_seed": False,
        "pocket": {
            "prob": 0.0,
            "size": 1 / 3,
            "spec_binder_chain": False,
            "max_distance_range": {"PP": ListValue([6, 20]), "LP": ListValue([6, 20])},
            "group": "complex",
            "distance_type": "center_atom",
        },
        "contact": {
            "prob": 0.0,
            "size": 1 / 3,
            "max_distance_range": {
                "PP": ListValue([6, 30]),
                "PL": ListValue([4, 10]),
            },
            "group": "complex",
            "distance_type": "center_atom",
        },
        "substructure": {
            "prob": 0.0,
            "size": 0.8,
            "mol_type_pairs": {
                "PP": 15,
                "PL": 10,
                "LP": 10,
            },
            "feature_type": "one_hot",
            "ratios": {
                "full": [
                    0.0,
                    0.5,
                    1.0,
                ],  # ratio options of full chain substructure constraint
                "partial": 0.3,  # ratio of partial chain substructure constraint
            },
            "coord_noise_scale": 0.05,
            "spec_asym_id": False,
        },
        "contact_atom": {
            "prob": 0.0,
            "size": 1 / 3,
            "max_distance_range": {
                "PP": ListValue([2, 12]),
                "PL": ListValue([2, 8]),
            },
            "min_distance": -1,
            "group": "complex",
            "distance_type": "atom",
            "feature_type": "continuous",
        },
    },
}


def _nanobody_base_info(task_type: str, split_name: str, is_train: bool) -> dict:
    base_info = {
        "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "nanobody", "mmcif", task_type),
        "bioassembly_dict_dir": os.path.join(
            PROTENIX_ROOT_DIR, "nanobody", "bioassembly", task_type
        ),
        "indices_fpath": os.path.join(
            PROTENIX_ROOT_DIR,
            "nanobody",
            "indices",
            f"nanobody_{task_type}_{split_name}.csv",
        ),
        "pdb_list": "",
        "max_n_token": -1 if is_train else GlobalConfigValue("test_max_n_token"),
        "use_reference_chains_only": True,
    }
    if is_train:
        base_info["random_sample_if_failed"] = True
    else:
        base_info["sort_by_n_token"] = True
        base_info["group_by_pdb_id"] = True
        base_info["find_eval_chain_interface"] = task_type == "complex"
    return base_info


def _nanobody_train_config(task_type: str, split_name: str) -> dict:
    # Generated by scripts/nanobody/prepare_nanobody_data.py.
    return {
        "base_info": _nanobody_base_info(task_type, split_name, is_train=True),
        **deepcopy(default_weighted_pdb_configs),
    }


def _nanobody_eval_config(task_type: str, split_name: str) -> dict:
    # Generated by scripts/nanobody/prepare_nanobody_data.py.
    return {
        "base_info": _nanobody_base_info(task_type, split_name, is_train=False),
        **deepcopy(default_test_configs),
    }


def _openfold_subset_train_config() -> dict:
    # Generated by scripts/nanobody/prepare_openfold_subset.py.
    return {
        "base_info": {
            "mmcif_dir": os.path.join(
                PROTENIX_ROOT_DIR, "nanobody", "mmcif", "openfold_subset"
            ),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "nanobody", "bioassembly", "openfold_subset"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR,
                "nanobody",
                "indices",
                "openfold_subset_train.csv",
            ),
            "pdb_list": "",
            "random_sample_if_failed": True,
            "max_n_token": -1,
            "use_reference_chains_only": False,
        },
        **deepcopy(default_weighted_pdb_configs),
    }


def _tiny_finetune_msa_template_overrides() -> dict:
    return {
        "msa": {
            "enable_prot_msa": False,
            "enable_rna_msa": False,
        },
        "template": {
            "enable_prot_template": False,
        },
    }


def _tiny_finetune_train_config() -> dict:
    return {
        "base_info": {
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "tiny_finetune", "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "tiny_finetune", "bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR, "tiny_finetune", "indices", "train.csv"
            ),
            "pdb_list": "",
            "random_sample_if_failed": True,
            "max_n_token": -1,
            "use_reference_chains_only": False,
        },
        **deepcopy(default_weighted_pdb_configs),
        **_tiny_finetune_msa_template_overrides(),
    }


def _tiny_finetune_val_config() -> dict:
    tiny_eval_configs = deepcopy(default_test_configs)
    tiny_eval_configs["cropping_configs"]["method_weights"] = ListValue(
        [1.0, 0.0, 0.0]
    )
    return {
        "base_info": {
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "tiny_finetune", "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "tiny_finetune", "bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR, "tiny_finetune", "indices", "val.csv"
            ),
            "pdb_list": "",
            "max_n_token": GlobalConfigValue("test_max_n_token"),
            "sort_by_n_token": True,
            "group_by_pdb_id": True,
            "find_eval_chain_interface": False,
            "use_reference_chains_only": False,
        },
        **tiny_eval_configs,
        **_tiny_finetune_msa_template_overrides(),
    }


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
    # NOTE:
    # `weightedPDB_before2109_wopb_nometalc_0925` is compatible with the `2024.05.22` data version
    # downloaded via `scripts/database/download_protenix_data.sh`.
    # `weightedPDB_before250701_v20260101` and `weightedPDB_before210930_v20260101` are compatible
    # with the `2026.01.01` data version.
    "weightedPDB_before2109_wopb_nometalc_0925": {
        "base_info": {
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "mmcif_bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR,
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
    "weightedPDB_before250701_v20260101": {
        "base_info": {
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "mmcif_bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR,
                "indices/indices_20260107-20chains_before_2025-07-01_res4.5.csv.gz",
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
    "weightedPDB_before210930_v20260101": {
        "base_info": {
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "mmcif_bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR,
                "indices/indices_20260107-20chains_before_2021-09-30_res4.5.csv.gz",
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
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "recentPDB_bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR, "indices/recentPDB_low_homology_maxtoken1536.csv"
            ),
            "pdb_list": os.path.join(
                PROTENIX_ROOT_DIR,
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
            "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "posebusters_mmcif"),
            "bioassembly_dict_dir": os.path.join(
                PROTENIX_ROOT_DIR, "posebusters_bioassembly"
            ),
            "indices_fpath": os.path.join(
                PROTENIX_ROOT_DIR, "indices/posebusters_indices_mainchain_interface.csv"
            ),
            "pdb_list": "",
            "find_pocket": True,
            "find_all_pockets": False,
            "max_n_token": GlobalConfigValue("test_max_n_token"),  # filter data
        },
        **deepcopy(default_test_configs),
    },
    "nanobody_single_train_pre2022": _nanobody_train_config(
        "single", "train_pre2022"
    ),
    "nanobody_complex_train_pre2022": _nanobody_train_config(
        "complex", "train_pre2022"
    ),
    "nanobody_single_val_2022h1": _nanobody_eval_config("single", "val_2022h1"),
    "nanobody_complex_val_2022h1": _nanobody_eval_config("complex", "val_2022h1"),
    "nanobody_single_test_2022h2": _nanobody_eval_config("single", "test_2022h2"),
    "nanobody_complex_test_2022h2": _nanobody_eval_config("complex", "test_2022h2"),
    "openfold_subset_train": _openfold_subset_train_config(),
    "tiny_protenix_finetune_train": _tiny_finetune_train_config(),
    "tiny_protenix_finetune_val": _tiny_finetune_val_config(),
    "msa": {
        "enable_prot_msa": True,
        "prot_seq_or_filename_to_msadir_jsons": ListValue(
            [os.path.join(PROTENIX_ROOT_DIR, "common/seq_to_pdb_index.json")]
        ),
        "prot_msadir_raw_paths": ListValue(
            [os.path.join(PROTENIX_ROOT_DIR, "mmcif_msa_template")]
        ),
        "prot_pairing_dbs": ListValue(["pairing"]),
        "prot_non_pairing_dbs": ListValue(
            ["pairing-non_pairing"]
        ),  # Separated by "-", "pairing-non_pairing" means both pairing and non_pairing are used as non_pairing,
        # with pairing used first.
        "prot_indexing_methods": ListValue(["sequence"]),
        "enable_rna_msa": True,  # enable rna msa
        "rna_seq_or_filename_to_msadir_jsons": ListValue(
            [os.path.join(PROTENIX_ROOT_DIR, "rna_msa/rna_sequence_to_pdb_chains.json")]
        ),
        "rna_msadir_raw_paths": ListValue(
            [os.path.join(PROTENIX_ROOT_DIR, "rna_msa/msas")]
        ),
        "rna_indexing_methods": ListValue(["sequence"]),
        "min_size": {
            "train": 1,
            "test": 1,
        },
        "max_size": {
            "train": 16384,
            "test": 16384,
        },
        "sample_cutoff": {
            "train": 16384,
            "test": 16384,
        },
    },
    "template": {
        "enable_prot_template": True,
        "template_dropout_rate": 0.0,
        "prot_template_mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
        "prot_template_cache_dir": "",
        "prot_template_raw_paths": ListValue(
            [os.path.join(PROTENIX_ROOT_DIR, "mmcif_msa_template")]
        ),
        "prot_seq_or_filename_to_templatedir_jsons": ListValue(
            [os.path.join(PROTENIX_ROOT_DIR, "common/seq_to_pdb_index.json")]
        ),
        "prot_indexing_methods": ListValue(["sequence"]),
        "release_dates_path": os.path.join(
            PROTENIX_ROOT_DIR, "common/release_date_cache.json"
        ),
        "obsolete_pdbs_path": os.path.join(
            PROTENIX_ROOT_DIR, "common/obsolete_to_successor.json"
        ),
        "kalign_binary_path": "/usr/bin/kalign",  # apt-get install kalign
    },
    "ccd_components_file": os.path.join(PROTENIX_ROOT_DIR, "common/components.cif"),
    "ccd_components_rdkit_mol_file": os.path.join(
        PROTENIX_ROOT_DIR, "common/components.cif.rdkit_mol.pkl"
    ),
    "obsolete_release_data_csv": os.path.join(
        PROTENIX_ROOT_DIR, "common/obsolete_release_date.csv"
    ),
    "pdb_cluster_file": os.path.join(
        PROTENIX_ROOT_DIR, "common/clusters-by-entity-40.txt"
    ),
}
