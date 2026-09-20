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

"""Failed RDKit embedding must preserve the selected CCD fallback and its mask."""

from pathlib import Path
from types import SimpleNamespace

import gemmi
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from scripts import gen_ccd_cache


def _component(monkeypatch, valid_ideal):
    mol = Chem.MolFromSmiles("CCO")
    names = ["C1", "C2", "O1"]
    for atom, name in zip(mol.GetAtoms(), names):
        atom.SetProp("name", name)
    for offset in (0.0, 10.0):
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for index in range(mol.GetNumAtoms()):
            conformer.SetAtomPosition(index, (index + offset, offset, 0.0))
        mol.AddConformer(conformer, assignId=True)
    rows = "\n".join(
        f"{name} 0.0 {index if valid else '?'}"
        for index, (name, valid) in enumerate(zip(names, valid_ideal))
    )
    document = gemmi.cif.read_string(
        "data_TST\nloop_\n_chem_comp_atom.atom_id\n"
        "_chem_comp_atom.model_Cartn_x\n"
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal\n" + rows + "\n"
    )
    monkeypatch.setattr(gen_ccd_cache, "gemmi_load_ccd_cif", lambda _: document)
    monkeypatch.setattr(
        gen_ccd_cache.ccd_reader,
        "_parse_pdb_mmcif",
        lambda block, sanitize: SimpleNamespace(
            component=SimpleNamespace(mol=mol), sanitized=True
        ),
    )
    return mol


@pytest.mark.parametrize("valid_ideal", [[True, False, True], [False, False, False]])
def test_real_failed_embedding_preserves_ideal_metadata(monkeypatch, valid_ideal):
    mol = _component(monkeypatch, valid_ideal)
    original = [conf.GetPositions().copy() for conf in mol.GetConformers()]
    embed = AllChem.EmbedMolecule
    returned_ids = []

    def impossible_embedding(molecule, options):
        # Inconsistent triangle bounds force the actual RDKit failure sentinel.
        bounds = AllChem.GetMoleculeBoundsMatrix(molecule)
        bounds[0, 2] = 100.0
        bounds[2, 0] = 100.0
        options.SetBoundsMat(bounds)
        result = embed(molecule, options)
        returned_ids.append(result)
        return result

    monkeypatch.setattr(AllChem, "EmbedMolecule", impossible_embedding)
    result = gen_ccd_cache._get_component_rdkit_mol_processing(("TST", Path("unused")))
    assert returned_ids == [-1]
    assert result is mol
    assert result.ref_conf_id == 0
    assert result.ref_conf_type == "idea"
    np.testing.assert_array_equal(result.ref_mask, valid_ideal)
    assert result.GetNumConformers() == 2
    for conformer, expected in zip(result.GetConformers(), original):
        np.testing.assert_array_equal(conformer.GetPositions(), expected)
    assert result.atom_map == {"C1": 0, "C2": 1, "O1": 2}


def test_real_successful_embedding_promotes_generated_conformer(monkeypatch):
    mol = _component(monkeypatch, [False, False, False])
    original = mol.GetConformer(0).GetPositions().copy()
    embed = AllChem.EmbedMolecule

    def deterministic_embedding(molecule, options):
        options.randomSeed = 42
        return embed(molecule, options)

    monkeypatch.setattr(AllChem, "EmbedMolecule", deterministic_embedding)
    result = gen_ccd_cache._get_component_rdkit_mol_processing(("TST", Path("unused")))
    assert result.ref_conf_id == 2
    assert result.ref_conf_type == "rdkit"
    assert result.ref_mask.all()
    assert np.isfinite(result.GetConformer(result.ref_conf_id).GetPositions()).all()
    np.testing.assert_array_equal(result.GetConformer(0).GetPositions(), original)


def test_raised_embedding_failure_preserves_ideal_metadata(monkeypatch):
    _component(monkeypatch, [True, False, True])

    def failing_embedding(molecule, options):
        raise ValueError("embedding failed")

    monkeypatch.setattr(AllChem, "EmbedMolecule", failing_embedding)
    result = gen_ccd_cache._get_component_rdkit_mol_processing(("TST", Path("unused")))
    assert result.ref_conf_id == 0
    assert result.ref_conf_type == "idea"
    np.testing.assert_array_equal(result.ref_mask, [True, False, True])
