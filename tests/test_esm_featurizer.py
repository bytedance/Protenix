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

import numpy as np
import pandas as pd
import pytest
import torch
from biotite.structure import AtomArray

from protenix.data.esm.esm_featurizer import ESMFeaturizer
from protenix.data.tokenizer import Token, TokenArray


def _inputs(inference_mode):
    atoms = AtomArray(4)
    atoms.set_annotation("label_entity_id", np.array(["1", "1", "2", "2"]))
    atoms.set_annotation("chain_mol_type", np.array(["protein"] * 4))
    atoms.res_id = np.array([1, 2, 1, 2])
    tokens = TokenArray([Token(0, centre_atom_index=i) for i in range(4)])
    if inference_mode:
        sample = {
            "name": "inference_case",
            "sequences": [
                {"proteinChain": {"sequence": "AC", "count": 1}},
                {"proteinChain": {"sequence": "DE", "count": 1}},
            ],
        }
    else:
        sample = {"pdb_id": "training_case", "sequences": {"1": "AC", "2": "DE"}}
    return tokens, atoms, sample


def _featurizer(tmp_path, failure=None, save_errors=True):
    cache = tmp_path / "cache"
    cache.mkdir()
    first = torch.arange(1, 7, dtype=torch.float32).reshape(2, 3)
    second = first + 10
    torch.save(first, cache / "first.pt")
    if failure != "missing_file":
        torch.save(
            second[:1] if failure == "wrong_length" else second, cache / "second.pt"
        )
    rows = [{"seq": "AC", "part_id": "cache", "seq_label": "first"}]
    if failure != "missing_sequence":
        rows.append({"seq": "DE", "part_id": "cache", "seq_label": "second"})
    sequences = tmp_path / "sequences.csv"
    pd.DataFrame(rows).to_csv(sequences, index=False)
    featurizer = ESMFeaturizer(
        embedding_dir=str(tmp_path),
        sequence_fpath=str(sequences),
        embedding_dim=3,
        error_dir=str(tmp_path) if save_errors else None,
    )
    return featurizer, first, second


@pytest.mark.parametrize("inference_mode", [True, False])
@pytest.mark.parametrize("save_errors", [True, False])
@pytest.mark.parametrize(
    "failure", ["missing_file", "wrong_length", "missing_sequence"]
)
def test_cached_embedding_failure_keeps_sample_diagnostics(
    tmp_path, caplog, inference_mode, save_errors, failure
):
    featurizer, first, _ = _featurizer(tmp_path, failure, save_errors)
    tokens, atoms, sample = _inputs(inference_mode)
    sample_id = sample["name" if inference_mode else "pdb_id"]
    with caplog.at_level(logging.WARNING):
        result = featurizer(tokens, atoms, sample, inference_mode=inference_mode)

    torch.testing.assert_close(result[:2], first)
    torch.testing.assert_close(result[2:], torch.zeros(2, 3))
    expected_error = {
        "missing_file": "FileNotFoundError",
        "wrong_length": "AssertionError",
        "missing_sequence": "KeyError: 'DE'",
    }[failure]
    assert f"[{sample_id}] ESM error:" in caplog.text
    assert expected_error in caplog.text
    assert "KeyError: 'pdb_id'" not in caplog.text
    if save_errors:
        errors = list((tmp_path / "esm_error").iterdir())
        assert [p.name for p in errors] == [f"{sample_id}_2.txt"]
        assert expected_error in errors[0].read_text()
    else:
        assert not (tmp_path / "esm_error").exists()


@pytest.mark.parametrize("inference_mode", [True, False])
def test_cached_embeddings_preserve_all_entity_values(tmp_path, caplog, inference_mode):
    featurizer, first, second = _featurizer(tmp_path)
    tokens, atoms, sample = _inputs(inference_mode)
    result = featurizer(tokens, atoms, sample, inference_mode=inference_mode)

    torch.testing.assert_close(result, torch.cat([first, second]))
    assert "ESM error:" not in caplog.text
    assert not list((tmp_path / "esm_error").iterdir())
