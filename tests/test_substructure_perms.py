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

import itertools
import logging

import numpy as np
import pytest

from protenix.data.core import substructure_perms
from rdkit import Chem


def test_invalid_optional_neutralization_preserves_alf_maps(caplog):
    mol = Chem.MolFromSmiles("[Al-](F)(F)(F)F")
    original = Chem.MolToMolBlock(mol)
    with pytest.raises(Chem.rdchem.AtomValenceException):
        substructure_perms._get_substructure_perms(mol, Neutralize=True)
    with caplog.at_level(logging.WARNING):
        actual = substructure_perms.get_substructure_perms(mol)
    expected = {(0, *p) for p in itertools.permutations(range(1, 5))}
    assert set(map(tuple, actual)) == expected
    assert Chem.MolToMolBlock(mol) == original
    assert "Skipping invalid neutralized symmetry state" in caplog.text


@pytest.mark.parametrize("smiles", ["c1ccccc1", "CC(=O)[O-]", "C[NH2+]C"])
def test_valid_states_retain_union_and_keep_protonation(smiles):
    mol = Chem.MolFromSmiles(smiles)
    original = substructure_perms._get_substructure_perms(mol, Neutralize=False)
    neutral = substructure_perms._get_substructure_perms(mol, Neutralize=True)
    expected = np.unique(np.vstack((original, neutral)), axis=0)
    np.testing.assert_array_equal(
        substructure_perms.get_substructure_perms(mol), expected
    )
    np.testing.assert_array_equal(
        substructure_perms.get_substructure_perms(mol, KeepProtonation=True), original
    )


def test_original_state_sanitization_failure_propagates():
    mol = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
    with pytest.raises(Chem.rdchem.AtomValenceException):
        substructure_perms.get_substructure_perms(mol)


def test_unexpected_optional_failure_propagates(monkeypatch):
    original = substructure_perms._get_substructure_perms

    def unexpected_failure(mol, Neutralize=False, **kwargs):
        if Neutralize:
            raise RuntimeError("unexpected implementation failure")
        return original(mol, Neutralize=False, **kwargs)

    monkeypatch.setattr(
        substructure_perms, "_get_substructure_perms", unexpected_failure
    )
    with pytest.raises(RuntimeError, match="unexpected implementation failure"):
        substructure_perms.get_substructure_perms(Chem.MolFromSmiles("CC"))


def test_failed_optional_state_still_honors_match_limit():
    actual = substructure_perms.get_substructure_perms(
        Chem.MolFromSmiles("[Al-](F)(F)(F)F"), MaxMatches=6
    )
    assert len(actual) <= 6
    assert all(sorted(row) == list(range(5)) for row in actual.tolist())
