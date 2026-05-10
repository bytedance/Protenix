from pathlib import Path

import numpy as np
import pytest

from protenix.data.constants import ATOM37_NUM
from protenix.data.template import template_utils
from protenix.data.template.template_path import template_paths_need_mmcif_dir
from protenix.data.template.template_parser import TemplateParser
from protenix.data.template.template_utils import TemplateHitFeaturizer
from runner import template_search


AA3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
}

BACKBONE_ATOMS = (
    ("N", "N", 0.0, 0.0, 0.0),
    ("CA", "C", 1.3, 0.0, 0.0),
    ("C", "C", 2.6, 0.0, 0.0),
    ("O", "O", 3.2, 0.0, 0.0),
)


def _simple_cif(chains: dict[str, str]) -> str:
    lines = [
        "data_custom",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.id",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.occupancy",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.pdbx_PDB_model_num",
    ]
    atom_id = 1
    for entity_id, (chain_id, sequence) in enumerate(chains.items(), start=1):
        for res_idx, aa in enumerate(sequence, start=1):
            resname = AA3[aa]
            for atom_name, element, dx, dy, dz in BACKBONE_ATOMS:
                x = (res_idx - 1) * 3.8 + dx
                lines.append(
                    f"ATOM {element} {atom_name} . {resname} {chain_id} "
                    f"{entity_id} {res_idx} . {res_idx} {resname} "
                    f"{chain_id} {atom_name} {atom_id} 20.00 1.00 "
                    f"{x:.3f} {dy:.3f} {dz:.3f} 1"
                )
                atom_id += 1
    lines.append("#")
    return "\n".join(lines) + "\n"


def _write_cif(tmp_path: Path, chains: dict[str, str]) -> Path:
    path = tmp_path / "custom-template.cif"
    path.write_text(_simple_cif(chains))
    return path


def _featurizer(tmp_path: Path) -> TemplateHitFeaturizer:
    return TemplateHitFeaturizer(
        mmcif_dir=str(tmp_path / "mmcif"),
        template_cache_dir=None,
        kalign_binary_path=None,
        _zero_center_positions=True,
    )


def test_simple_cif_parser_keeps_all_protein_chains():
    result = TemplateParser.parse_simple_cif(
        file_id="custom",
        mmcif_string=_simple_cif({"A": "ACDEFG", "B": "GGGGG"}),
    )

    assert result.mmcif_object is not None
    assert result.mmcif_object.chain_to_seqres == {
        "A": "ACDEFG",
        "B": "GGGGG",
    }


def test_direct_cif_template_requires_chain_for_multichain(tmp_path):
    template_path = _write_cif(tmp_path, {"A": "ACDEFG", "B": "GGGGG"})

    result = _featurizer(tmp_path).parse_cif_template(
        str(template_path),
        query_sequence="ACDEFG",
    )

    assert result.features == []
    assert "multiple protein chains" in result.errors[0]
    assert "templateChainId" in result.errors[0]


def test_direct_cif_template_uses_selected_chain(tmp_path):
    template_path = _write_cif(tmp_path, {"A": "ACDEFG", "B": "GGGGG"})

    result = _featurizer(tmp_path).parse_cif_template(
        str(template_path),
        query_sequence="GGGGG",
        chain_id="B",
    )

    assert result.errors == []
    assert len(result.features) == 1
    feature = result.features[0]
    assert feature["template_all_atom_positions"].shape == (5, ATOM37_NUM, 3)
    assert feature["template_all_atom_masks"].shape == (5, ATOM37_NUM)
    assert feature["template_sequence"] == b"GGGGG"
    assert feature["template_release_date"].item() == b"1970-01-01"
    assert feature["template_domain_names"].item() == b"custom-template_B"
    assert np.sum(feature["template_all_atom_masks"]) > 0
    assert result.hits[0].name == "custom-template_B"


def test_direct_cif_template_reports_missing_chain(tmp_path):
    template_path = _write_cif(tmp_path, {"A": "ACDEFG"})

    result = _featurizer(tmp_path).parse_cif_template(
        str(template_path),
        query_sequence="ACDEFG",
        chain_id="B",
    )

    assert result.features == []
    assert "was not found" in result.errors[0]
    assert "Available chains" in result.errors[0]


def test_direct_cif_template_aligns_non_identical_query(tmp_path, monkeypatch):
    template_path = _write_cif(tmp_path, {"A": "ACDEFG"})

    class FakeKalign:
        def __init__(self, binary_path):
            self.binary_path = binary_path

        def align(self, sequences):
            assert sequences == ["AACDEFG", "ACDEFG"]
            return ["AACDEFG", "-ACDEFG"]

    monkeypatch.setattr(template_utils, "Kalign", FakeKalign)

    result = _featurizer(tmp_path).parse_cif_template(
        str(template_path),
        query_sequence="AACDEFG",
    )

    assert result.errors == []
    assert result.features[0]["template_sequence"] == b"-ACDEFG"
    assert result.hits[0].indices_hit[0] == -1
    assert result.hits[0].indices_hit[1] == 0
    assert np.sum(result.features[0]["template_all_atom_masks"][0]) == 0
    assert np.sum(result.features[0]["template_all_atom_masks"][1]) > 0


def test_template_info_respects_existing_explicit_template_path(tmp_path):
    template_path = _write_cif(tmp_path, {"A": "ACDEFG"})
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "pairing.a3m").write_text(">query\nACDEFG\n")
    data = [
        {
            "name": "custom",
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "ACDEFG",
                        "count": 1,
                        "pairedMsaPath": str(msa_dir / "pairing.a3m"),
                        "templatesPath": str(template_path),
                    }
                }
            ],
        }
    ]

    assert template_search.update_template_info(data) is False
    assert data[0]["sequences"][0]["proteinChain"]["templatesPath"] == str(
        template_path
    )


def test_template_info_raises_for_missing_explicit_template_path(tmp_path):
    data = [
        {
            "name": "custom",
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "ACDEFG",
                        "count": 1,
                        "templatesPath": str(tmp_path / "missing.cif"),
                    }
                }
            ],
        }
    ]

    with pytest.raises(FileNotFoundError, match="templatesPath"):
        template_search.update_template_info(data)


def test_template_info_autosearches_when_template_path_absent(tmp_path, monkeypatch):
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "pairing.a3m").write_text(">query\nACDEFG\n")
    calls = []

    def fake_run_template_search(
        msa_for_template_search_dir,
        msa_for_template_search_name,
        hmmsearch_binary_path=None,
        hmmbuild_binary_path=None,
        seqres_database_path=None,
    ):
        calls.append((msa_for_template_search_dir, msa_for_template_search_name))
        Path(msa_for_template_search_dir, "hmmsearch.a3m").write_text("")

    monkeypatch.setattr(template_search, "run_template_search", fake_run_template_search)
    data = [
        {
            "name": "custom",
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "ACDEFG",
                        "count": 1,
                        "pairedMsaPath": str(msa_dir / "pairing.a3m"),
                    }
                }
            ],
        }
    ]

    assert template_search.update_template_info(data) is True
    assert calls == [(str(msa_dir), "pairing")]
    assert data[0]["sequences"][0]["proteinChain"]["templatesPath"] == str(
        msa_dir / "hmmsearch.a3m"
    )


def test_template_mmcif_dir_gate_only_requires_search_hit_formats():
    direct_inputs = [
        {
            "sequences": [
                {"proteinChain": {"templatesPath": "/tmp/template.cif"}},
                {"proteinChain": {"templatesPath": "/tmp/template.mmcif"}},
                {"proteinChain": {"templatesPath": "/tmp/template.json"}},
            ]
        }
    ]
    search_inputs = [
        {
            "sequences": [
                {"proteinChain": {"templatesPath": "/tmp/hmmsearch.a3m"}},
                {"proteinChain": {"templatesPath": "/tmp/concat.hhr"}},
            ]
        }
    ]

    assert template_paths_need_mmcif_dir(direct_inputs) is False
    assert template_paths_need_mmcif_dir(search_inputs) is True
