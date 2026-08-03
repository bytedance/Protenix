import logging

import numpy as np

from protenix.data.msa.msa_featurizer import InferenceMSAFeaturizer


class FakeAtomArray:
    def __init__(self, asym_id_int, chain_id, res_id):
        self.asym_id_int = np.asarray(asym_id_int)
        self.chain_id = np.asarray(chain_id)
        self.res_id = np.asarray(res_id)
        self.centre_atom_mask = np.ones(len(self.res_id), dtype=bool)

    def __getitem__(self, item):
        return FakeAtomArray(
            asym_id_int=self.asym_id_int[item],
            chain_id=self.chain_id[item],
            res_id=self.res_id[item],
        )


def _two_chain_atom_array(seq_a: str, seq_b: str) -> FakeAtomArray:
    return FakeAtomArray(
        asym_id_int=[0] * len(seq_a) + [1] * len(seq_b),
        chain_id=["A"] * len(seq_a) + ["B"] * len(seq_b),
        res_id=list(range(1, len(seq_a) + 1)) + list(range(1, len(seq_b) + 1)),
    )


def test_inference_monomer_msa_path_produces_paired_rows(tmp_path):
    seq_a, seq_b = "ACDEF", "FGHIK"
    msa_a = tmp_path / "chain_a.a3m"
    msa_b = tmp_path / "chain_b.a3m"
    msa_a.write_text(
        ">query\nACDEF\n>UniRef100_HITA_9606/\nACDEF\n>env_a\nAC-EF\n",
        encoding="utf-8",
    )
    msa_b.write_text(
        ">query\nFGHIK\n>UniRef100_HITB_9606/\nFGHIK\n>env_b\nF-HIK\n",
        encoding="utf-8",
    )
    bioassembly = [
        {"proteinChain": {"sequence": seq_a, "count": 1, "monomerMsaPath": str(msa_a)}},
        {"proteinChain": {"sequence": seq_b, "count": 1, "monomerMsaPath": str(msa_b)}},
    ]

    features = InferenceMSAFeaturizer.make_msa_feature(
        bioassembly=bioassembly,
        atom_array=_two_chain_atom_array(seq_a, seq_b),
        msa_pair_as_unpair=False,
    )

    assert int(features["prot_paired_num_alignments"]) == 2
    assert int(features["prot_unpaired_num_alignments"]) >= 1
    assert features["msa"].shape[1] == len(seq_a) + len(seq_b)


def test_inference_monomer_msa_path_without_taxonomy_is_unpaired_only(tmp_path, caplog):
    seq_a, seq_b = "ACDEF", "FGHIK"
    msa_a = tmp_path / "chain_a.a3m"
    msa_b = tmp_path / "chain_b.a3m"
    msa_a.write_text(">query\nACDEF\n>env_a\nAC-EF\n", encoding="utf-8")
    msa_b.write_text(">query\nFGHIK\n>env_b\nF-HIK\n", encoding="utf-8")
    bioassembly = [
        {"proteinChain": {"sequence": seq_a, "count": 1, "monomerMsaPath": str(msa_a)}},
        {"proteinChain": {"sequence": seq_b, "count": 1, "monomerMsaPath": str(msa_b)}},
    ]

    with caplog.at_level(logging.WARNING):
        features = InferenceMSAFeaturizer.make_msa_feature(
            bioassembly=bioassembly,
            atom_array=_two_chain_atom_array(seq_a, seq_b),
            msa_pair_as_unpair=False,
        )

    assert int(features["prot_paired_num_alignments"]) == 1
    assert int(features["prot_unpaired_num_alignments"]) == 1
    assert "using it as unpaired MSA only" in caplog.text
