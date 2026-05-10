import pytest

from protenix.data.msa.msa_input import split_monomer_msa
from protenix.data.msa.msa_utils import MSAPairingEngine
from protenix.data.tools.common import parse_fasta
from runner.msa_search import need_msa_search


def test_species_ids_support_explicit_taxonomy_tags():
    descriptions = [
        "UniRef100_HITA_9606/",
        "sp|P12345|ABC_HUMAN Some protein OX=9606",
        "hit TaxID=10090 description",
        "hit taxid=7227 description",
        "hit Tax=7955 description",
        "speciesless_hit_42",
    ]

    assert MSAPairingEngine.get_species_ids(descriptions) == [
        "9606",
        "HUMAN",
        "10090",
        "7227",
        "7955",
        "",
    ]


def test_split_monomer_msa_separates_pairable_and_unpairable_rows():
    split = split_monomer_msa(
        query_sequence="ACDE",
        a3m=(
            ">query\n"
            "ACDE\n"
            ">UniRef100_HITA_9606/\n"
            "ACDE\n"
            ">environmental_hit\n"
            "ACdDE\n"
            ">bad_length\n"
            "ACD\n"
        ),
        source_name="test.a3m",
    )

    paired_seqs, paired_descs = parse_fasta(split.paired_a3m)
    unpaired_seqs, unpaired_descs = parse_fasta(split.unpaired_a3m)

    assert split.total_rows == 4
    assert split.pairable_rows == 1
    assert split.unpairable_rows == 1
    assert split.invalid_rows == 1
    assert paired_descs == ["query", "UniRef100_HITA_9606/"]
    assert paired_seqs == ["ACDE", "ACDE"]
    assert unpaired_descs == ["query", "UniRef100_HITA_9606/", "environmental_hit"]
    assert unpaired_seqs == ["ACDE", "ACDE", "ACdDE"]


def test_split_monomer_msa_rejects_query_length_mismatch():
    with pytest.raises(ValueError, match="query row aligned length"):
        split_monomer_msa(
            query_sequence="ACDE",
            a3m=">query\nACD\n",
            source_name="bad.a3m",
        )


def test_need_msa_search_accepts_existing_monomer_msa_path(tmp_path):
    msa_path = tmp_path / "chain_a.a3m"
    msa_path.write_text(">query\nACDE\n", encoding="utf-8")

    json_data = {
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "ACDE",
                    "count": 1,
                    "monomerMsaPath": str(msa_path),
                }
            }
        ]
    }

    assert need_msa_search(json_data) is False


def test_need_msa_search_researches_missing_monomer_msa_path(tmp_path):
    json_data = {
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "ACDE",
                    "count": 1,
                    "monomerMsaPath": str(tmp_path / "missing.a3m"),
                }
            }
        ]
    }

    assert need_msa_search(json_data) is True


def test_need_msa_search_honors_split_path_precedence(tmp_path):
    paired_msa_path = tmp_path / "pairing.a3m"
    paired_msa_path.write_text(">query\nACDE\n", encoding="utf-8")

    json_data = {
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "ACDE",
                    "count": 1,
                    "pairedMsaPath": str(paired_msa_path),
                    "monomerMsaPath": str(tmp_path / "missing_monomer.a3m"),
                }
            }
        ]
    }

    assert need_msa_search(json_data) is False
