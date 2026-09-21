# ruff: noqa: E501
import textwrap

import numpy as np
import pytest

from protenix.data.inference.json_to_feature import SampleDictToFeatures


def _component_block(code: str, one_letter_code: str) -> str:
    return textwrap.dedent(
        f"""
        data_{code}
        #
        _chem_comp.id                                    {code}
        _chem_comp.name                                  '{code} TEST COMPONENT'
        _chem_comp.type                                  'L-PEPTIDE LINKING'
        _chem_comp.pdbx_type                             ATOMP
        _chem_comp.formula                               'C3 H7 N O2'
        _chem_comp.mon_nstd_parent_comp_id               ?
        _chem_comp.pdbx_synonyms                         ?
        _chem_comp.pdbx_formal_charge                    0
        _chem_comp.pdbx_initial_date                     2026-05-10
        _chem_comp.pdbx_modified_date                    2026-05-10
        _chem_comp.pdbx_ambiguous_flag                   N
        _chem_comp.pdbx_release_status                   REL
        _chem_comp.pdbx_replaced_by                      ?
        _chem_comp.pdbx_replaces                         ?
        _chem_comp.formula_weight                        89.094
        _chem_comp.one_letter_code                       {one_letter_code}
        _chem_comp.three_letter_code                     {code}
        _chem_comp.pdbx_model_coordinates_db_code        ?
        _chem_comp.pdbx_model_coordinates_details        ?
        _chem_comp.pdbx_ideal_coordinates_details        ?
        _chem_comp.pdbx_ideal_coordinates_missing_flag   N
        _chem_comp.pdbx_model_coordinates_missing_flag   N
        _chem_comp.pdbx_processing_site                  ?
        #
        loop_
        _chem_comp_atom.comp_id
        _chem_comp_atom.atom_id
        _chem_comp_atom.alt_atom_id
        _chem_comp_atom.type_symbol
        _chem_comp_atom.charge
        _chem_comp_atom.pdbx_align
        _chem_comp_atom.pdbx_aromatic_flag
        _chem_comp_atom.pdbx_leaving_atom_flag
        _chem_comp_atom.pdbx_stereo_config
        _chem_comp_atom.model_Cartn_x
        _chem_comp_atom.model_Cartn_y
        _chem_comp_atom.model_Cartn_z
        _chem_comp_atom.pdbx_model_Cartn_x_ideal
        _chem_comp_atom.pdbx_model_Cartn_y_ideal
        _chem_comp_atom.pdbx_model_Cartn_z_ideal
        _chem_comp_atom.pdbx_component_atom_id
        _chem_comp_atom.pdbx_component_comp_id
        _chem_comp_atom.pdbx_ordinal
        {code} N   N   N 0 1 N N N  0.000  0.000  0.000  0.000  0.000  0.000 N   {code} 1
        {code} CA  CA  C 0 1 N N N  1.450  0.000  0.000  1.450  0.000  0.000 CA  {code} 2
        {code} C   C   C 0 1 N N N  2.020  1.410  0.000  2.020  1.410  0.000 C   {code} 3
        {code} O   O   O 0 1 N N N  1.320  2.400  0.000  1.320  2.400  0.000 O   {code} 4
        {code} CB  CB  C 0 1 N N N  1.980 -0.780 -1.200  1.980 -0.780 -1.200 CB  {code} 5
        {code} OXT OXT O 0 1 N Y N  3.250  1.560  0.000  3.250  1.560  0.000 OXT {code} 6
        #
        loop_
        _chem_comp_bond.comp_id
        _chem_comp_bond.atom_id_1
        _chem_comp_bond.atom_id_2
        _chem_comp_bond.value_order
        _chem_comp_bond.pdbx_aromatic_flag
        _chem_comp_bond.pdbx_stereo_config
        _chem_comp_bond.pdbx_ordinal
        {code} N  CA  SING N N 1
        {code} CA C   SING N N 2
        {code} C  O   DOUB N N 3
        {code} CA CB  SING N N 4
        {code} C  OXT SING N N 5
        #
        """
    )


def _user_ccd_text() -> str:
    return "\n".join(
        [
            _component_block("ALA", "A"),
            _component_block("UAA", "K"),
            _component_block("THR", "T"),
        ]
    )


def _job_with_user_ccd(**extra):
    job = {
        "name": "custom_ptm",
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "AKT",
                    "count": 1,
                    "modifications": [{"ptmType": "CCD_UAA", "ptmPosition": 2}],
                }
            }
        ],
    }
    job.update(extra)
    return job


def _atom_index(atom_array, res_id: int, atom_name: str) -> int:
    indices = np.where(
        (atom_array.res_id == res_id) & (atom_array.atom_name == atom_name)
    )[0]
    assert len(indices) == 1
    return int(indices[0])


def _has_bond(atom_array, atom_i: int, atom_j: int) -> bool:
    bonds = atom_array.bonds.as_array()[:, :2]
    return any(set(pair) == {atom_i, atom_j} for pair in bonds)


def test_user_ccd_path_supports_custom_internal_protein_ptm(tmp_path):
    ccd_path = tmp_path / "components.cif"
    ccd_path.write_text(_user_ccd_text())
    sample = SampleDictToFeatures(
        _job_with_user_ccd(userCCDPath="components.cif"), input_json_dir=tmp_path
    )

    atom_array = sample.get_atom_array()

    custom_mask = atom_array.res_name == "UAA"
    assert custom_mask.any()
    assert set(atom_array.mol_type[custom_mask]) == {"protein"}
    assert not np.any(custom_mask & (atom_array.atom_name == "OXT"))
    assert np.all(atom_array.modified_res_mask[custom_mask] == 1)
    assert set(atom_array.cano_seq_resname[custom_mask]) == {"LYS"}

    assert _has_bond(
        atom_array,
        _atom_index(atom_array, 1, "C"),
        _atom_index(atom_array, 2, "N"),
    )
    assert _has_bond(
        atom_array,
        _atom_index(atom_array, 2, "C"),
        _atom_index(atom_array, 3, "N"),
    )
    assert np.all(atom_array.ref_mask[custom_mask] == 1)


def test_inline_user_ccd_matches_path_user_ccd(tmp_path):
    ccd_path = tmp_path / "components.cif"
    ccd_text = _user_ccd_text()
    ccd_path.write_text(ccd_text)

    from_path = SampleDictToFeatures(
        _job_with_user_ccd(userCCDPath="components.cif"), input_json_dir=tmp_path
    ).get_atom_array()
    from_inline = SampleDictToFeatures(
        _job_with_user_ccd(userCCD=ccd_text), input_json_dir=tmp_path
    ).get_atom_array()

    np.testing.assert_array_equal(from_path.res_name, from_inline.res_name)
    np.testing.assert_array_equal(from_path.atom_name, from_inline.atom_name)
    np.testing.assert_array_equal(from_path.mol_type, from_inline.mol_type)


def test_user_ccd_custom_ptm_builds_feature_dict_with_geometry(tmp_path):
    ccd_path = tmp_path / "components.cif"
    ccd_path.write_text(_user_ccd_text())
    sample = SampleDictToFeatures(
        _job_with_user_ccd(userCCDPath="components.cif"),
        extract_features_for_tfg=True,
        input_json_dir=tmp_path,
    )

    feature_dict, atom_array, token_array = sample.get_feature_dict()

    assert "UAA" in sample.input_dict["ccd_mols"]
    assert len(token_array) == int(np.sum(atom_array.centre_atom_mask))
    assert feature_dict["ref_pos"].shape[0] == len(atom_array)
    assert "pairwise_distance_index" in feature_dict


def test_user_ccd_rejects_mutually_exclusive_sources(tmp_path):
    ccd_path = tmp_path / "components.cif"
    ccd_path.write_text(_user_ccd_text())

    with pytest.raises(ValueError, match="Only one of"):
        SampleDictToFeatures(
            _job_with_user_ccd(userCCD=_user_ccd_text(), userCCDPath="components.cif"),
            input_json_dir=tmp_path,
        )


def test_user_ccd_rejects_missing_relative_path(tmp_path):
    with pytest.raises(FileNotFoundError, match="userCCDPath"):
        SampleDictToFeatures(
            _job_with_user_ccd(userCCDPath="missing.cif"), input_json_dir=tmp_path
        )


def test_user_ccd_rejects_empty_inline_ccd(tmp_path):
    with pytest.raises(ValueError, match="userCCD"):
        SampleDictToFeatures(_job_with_user_ccd(userCCD=""), input_json_dir=tmp_path)


def test_custom_ptm_position_is_validated(tmp_path):
    ccd_path = tmp_path / "components.cif"
    ccd_path.write_text(_user_ccd_text())
    job = _job_with_user_ccd(userCCDPath="components.cif")
    job["sequences"][0]["proteinChain"]["modifications"][0]["ptmPosition"] = 99

    with pytest.raises(ValueError, match="ptmPosition"):
        SampleDictToFeatures(job, input_json_dir=tmp_path)
