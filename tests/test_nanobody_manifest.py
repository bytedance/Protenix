import importlib

import pandas as pd

from protenix.nanobody.manifest import (
    assign_split,
    normalize_antigen_chains,
    normalize_metadata_table,
)


def test_metadata_column_normalization_and_antigen_chains(tmp_path):
    metadata = pd.DataFrame(
        [
            {
                "PDB": "1abc",
                "date": "2022-02-03",
                "Hchain": "H",
                "antigen_chain_ids": "A, B;A",
                "cdr3_start": 101,
                "cdr3_end": 112,
            }
        ]
    )
    manifest = normalize_metadata_table(metadata, structure_dir=tmp_path)

    assert set(manifest["task_type"]) == {"single", "complex"}
    complex_row = manifest[manifest["task_type"] == "complex"].iloc[0]
    assert complex_row["pdb_id"] == "1ABC"
    assert complex_row["antigen_chains"] == "A;B"
    assert complex_row["cdrh3_source"] == "metadata"
    assert complex_row["split"] == "val"
    assert complex_row["split_name"] == "val_2022h1"


def test_date_split_assignment_boundaries():
    assert assign_split("2021-12-31") == ("train", "train_pre2022")
    assert assign_split("2022-01-01") == ("val", "val_2022h1")
    assert assign_split("2022-06-30") == ("val", "val_2022h1")
    assert assign_split("2022-07-01") == ("test", "test_2022h2")
    assert assign_split("2023-01-01") == ("holdout", "holdout_post2022")
    assert assign_split("2021-10-01", strict_protenix_cutoff=True) == (
        "gap",
        "excluded_2021q4",
    )


def test_semicolon_normalization_for_antigen_chains():
    assert normalize_antigen_chains("A,B; C|A /D") == "A;B;C;D"
    assert normalize_antigen_chains("") == ""


def test_nanobody_dataset_config_paths_under_root(monkeypatch, tmp_path):
    monkeypatch.setenv("PROTENIX_ROOT_DIR", str(tmp_path))
    import configs.configs_data as configs_data

    configs_data = importlib.reload(configs_data)
    cfg = configs_data.data_configs["nanobody_complex_val_2022h1"]["base_info"]
    assert str(tmp_path) in cfg["mmcif_dir"]
    assert str(tmp_path) in cfg["bioassembly_dict_dir"]
    assert str(tmp_path) in cfg["indices_fpath"]
    assert cfg["find_eval_chain_interface"] is True
