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

"""Regression tests for ordered ColabFold MSA conversion, without live HTTP."""

import contextlib
import copy
import io
import re
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from protenix.web_service import colab_request_utils as writer


PAIRED_CHAINS = [
    [
        ("query_0", "ACDEFG"),
        ("UniRef100_A1", "ACDEFA"),
        ("DUMMY", "------"),
        ("UniRef100_A3", "ACDEFC"),
        ("DUMMY", "------"),
        ("UniRef100_A5", "ACDEFD"),
    ],
    [
        ("query_1", "HIKLMN"),
        ("UniRef100_B1", "HIKLMA"),
        ("UniRef100_B2", "HIKLMC"),
        ("UniRef100_B3", "HIKLMD"),
        ("UniRef100_B4", "HIKLME"),
        ("UniRef100_B5", "HIKLMF"),
    ],
    [
        ("query_2", "PQRSTV"),
        ("UniRef100_C1", "PQRSTA"),
        ("UniRef100_C2", "PQRSTC"),
        ("UniRef100_C3", "PQRSTD"),
        ("UniRef100_C4", "PQRSTE"),
        ("UniRef100_C5", "PQRSTF"),
    ],
]


def fasta(records):
    return "".join(f">{header}\n{sequence}\n" for header, sequence in records)


def paired_cases():
    """Include equal record-loss counts at different positions across chains."""
    for case in ["gap", "tabbed_gap", "ordinary", "different_gap_positions", "unique"]:
        chains = copy.deepcopy(PAIRED_CHAINS)
        if case == "tabbed_gap":
            for i in [2, 4]:
                chains[0][i] = ("DUMMY\t0\t0\t0\t1\t0\t0\t1\t0\t0", "------")
        elif case == "ordinary":
            chains[0][2] = ("UniRef100_REPEAT", "ACDEFE")
            chains[0][4] = ("UniRef100_REPEAT", "ACDEFF")
        elif case == "different_gap_positions":
            for i in [1, 3]:
                chains[1][i] = ("DUMMY", "------")
        elif case == "unique":
            for i in [2, 4]:
                chains[0][i] = (f"DUMMY_ROW{i}", "------")
        yield case, chains


def convert_response(files, query, prefix, use_pairing):
    """Run the real writer with a successful synthetic server response."""
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w:gz") as handle:
        for name, data in files.items():
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            handle.addfile(entry, io.BytesIO(data))
    response = Mock()
    response.json.return_value = {"status": "COMPLETE", "id": "synthetic"}
    response.content = archive.getvalue()
    with patch.object(writer.requests, "post", return_value=response), patch.object(
        writer.requests, "get", return_value=response
    ), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
        io.StringIO()
    ):
        return writer.run_mmseqs2_service(
            query,
            str(prefix),
            use_pairing=use_pairing,
            pairing_strategy="greedy",
            server_mode="colabfold",
            host_url="https://offline.invalid",
            user_agent="msa-conversion-test",
        )


def convert_paired(chains, prefix):
    raw = "".join(fasta(chain) + "\0" for chain in chains).encode()
    query = fasta([chain[0] for chain in chains]).rstrip("\n")
    convert_response({"pair.a3m": raw}, query, prefix, use_pairing=True)
    return [(prefix / str(i) / "pairing.a3m").read_text() for i in range(len(chains))]


class TestColabFoldConversion(unittest.TestCase):
    def test_paired_rows_keep_original_sequences_and_shared_ids(self):
        for case, chains in paired_cases():
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                outputs = convert_paired(chains, Path(directory) / "msa")
                for chain, output in zip(chains, outputs):
                    # The writer emits two lines per record, including gap rows.
                    lines = output.splitlines()
                    self.assertEqual(lines[1::2], [sequence for _, sequence in chain])
                    self.assertEqual(lines[0], ">query")
                    for row_id, header in enumerate(lines[0::2][1:], start=1):
                        match = re.search(r"_(\d+)/", header.split("\t", 1)[0])
                        self.assertIsNotNone(match)
                        self.assertEqual(int(match[1]), row_id)

    def test_unique_headers_keep_existing_output_format(self):
        chains = dict(paired_cases())["unique"]
        with tempfile.TemporaryDirectory() as directory:
            outputs = convert_paired(chains, Path(directory) / "msa")
        self.assertEqual(
            outputs[0],
            ">query\nACDEFG\n"
            ">UniRef100_A1_1/_1\nACDEFA\n"
            ">DUMMY_ROW2_2/_2\n------\n"
            ">UniRef100_A3_3/_3\nACDEFC\n"
            ">DUMMY_ROW4_4/_4\n------\n"
            ">UniRef100_A5_5/_5\nACDEFD\n",
        )

    def test_unpaired_conversion_keeps_distinct_hits_with_repeated_headers(self):
        files = {
            "bfd.mgnify30.metaeuk30.smag30.a3m": fasta(
                [("query_0", "ACDEFG"), ("env_hit", "ACDEFA"), ("env_hit", "ACDEFC")]
            ).encode(),
            "uniref.a3m": fasta(
                [
                    ("query_0", "ACDEFG"),
                    ("UniRef100_hit", "ACDEFD"),
                    ("UniRef100_hit", "ACDEFE"),
                ]
            ).encode(),
        }
        with tempfile.TemporaryDirectory() as directory:
            output_dir = convert_response(
                files, ">query_0\nACDEFG", Path(directory) / "msa", use_pairing=False
            )
            output = (Path(output_dir) / "non_pairing.a3m").read_text()
        self.assertEqual(
            output,
            ">query\nACDEFG\n"
            ">env_hit\nACDEFA\n>env_hit\nACDEFC\n"
            ">UniRef100_hit\nACDEFD\n>UniRef100_hit\nACDEFE\n",
        )

    def test_query_parser_retains_dictionary_interface(self):
        self.assertEqual(
            writer.parse_fasta_string(">query_0\nACD\nEFG\n>query_1\nHIKLMN\n"),
            {"query_0": "ACDEFG", "query_1": "HIKLMN"},
        )

    def test_query_mapping_preserves_duplicate_key_semantics(self):
        result = writer.parse_fasta_string(
            ">query_0\nACDEFG\n>query_1\nHIKLMN\n>query_0\nACDEFA\n"
        )
        self.assertEqual(list(result), ["query_0", "query_1"])
        self.assertEqual(result, {"query_0": "ACDEFA", "query_1": "HIKLMN"})


class TestFastaRecords(unittest.TestCase):
    def test_multiline_duplicate_headers_and_gap_records(self):
        self.assertEqual(
            writer.parse_fasta_records(
                ">query\nACD\nEFG\n>DUMMY\t0\n------\n"
                ">hit\nACdD\nEFG\n>DUMMY\t0\n------\n>hit\nACDEFA\n"
            ),
            [
                ("query", "ACDEFG"),
                ("DUMMY\t0", "------"),
                ("hit", "ACdDEFG"),
                ("DUMMY\t0", "------"),
                ("hit", "ACDEFA"),
            ],
        )

    def test_empty_input_and_empty_sequence_records(self):
        self.assertEqual(writer.parse_fasta_records("\n\n"), [])
        self.assertEqual(
            writer.parse_fasta_records(">hit\n>hit\nACD\n"),
            [("hit", ""), ("hit", "ACD")],
        )

    def test_header_metadata_whitespace_and_line_endings(self):
        self.assertEqual(
            writer.parse_fasta_records(
                "\r\n> hit\tfirst \r\n ACD \r\nEFG\r\n\r\n"
                ">hit\tsecond\r\nACDEFA\r\n>hit\tfirst\r\nACDEFC\r\n"
            ),
            [
                ("hit\tfirst", "ACDEFG"),
                ("hit\tsecond", "ACDEFA"),
                ("hit\tfirst", "ACDEFC"),
            ],
        )

    def test_fasta_comments_are_ignored(self):
        self.assertEqual(
            writer.parse_fasta_records(
                "# metadata\n>hit\nACD\n# comment\nEFG\n>hit\nACDEFA\n"
            ),
            [("hit", "ACDEFG"), ("hit", "ACDEFA")],
        )


class TestColabFoldInferenceFeatures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import numpy as np
            from biotite.structure import AtomArray
            from protenix.data.constants import PROTEIN_CHAIN
            from protenix.data.msa.msa_featurizer import InferenceMSAFeaturizer
            from protenix.data.msa.msa_utils import MSACore
        except ImportError as exc:
            raise unittest.SkipTest(
                f"Inference dependencies unavailable: {exc}"
            ) from exc
        cls.np = np
        cls.AtomArray = AtomArray
        cls.PROTEIN_CHAIN = PROTEIN_CHAIN
        cls.Featurizer = InferenceMSAFeaturizer
        cls.Core = MSACore

    def test_original_cross_chain_tuple_survives_final_feature_assembly(self):
        np = self.np
        correct, _ = self.Core.sequences_to_array(
            ["ACDEFDHIKLMFPQRSTF"], self.PROTEIN_CHAIN
        )
        wrong, _ = self.Core.sequences_to_array(
            ["ACDEFDHIKLMEPQRSTE"], self.PROTEIN_CHAIN
        )
        for case, chains in paired_cases():
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                outputs = convert_paired(chains, Path(directory) / "msa")
                atoms = self.AtomArray(18)
                atoms.set_annotation("asym_id_int", np.repeat(np.arange(3), 6))
                atoms.set_annotation("centre_atom_mask", np.ones(18, dtype=bool))
                atoms.chain_id = np.repeat(["A", "B", "C"], 6)
                atoms.res_id = np.tile(np.arange(1, 7), 3)
                entities = [
                    {
                        "proteinChain": dict(
                            sequence=chain[0][1],
                            count=1,
                            pairedMsa=output,
                            unpairedMsa=fasta([chain[0]]),
                        )
                    }
                    for chain, output in zip(chains, outputs)
                ]
                expected_rows, _ = self.Core.sequences_to_array(
                    [
                        "".join(chain[row][1] for chain in chains)
                        for row in range(len(chains[0]))
                    ],
                    self.PROTEIN_CHAIN,
                )
                original_tuples = {tuple(row) for row in expected_rows.tolist()}
                for flag in [False, True]:
                    with self.subTest(msa_pair_as_unpair=flag):
                        features = self.Featurizer.make_msa_feature(
                            entities, atoms, msa_pair_as_unpair=flag, use_rna_msa=False
                        )
                        self.assertTrue(
                            np.any(np.all(features["msa"] == correct[0], axis=1))
                        )
                        self.assertFalse(
                            np.any(np.all(features["msa"] == wrong[0], axis=1))
                        )
                        self.assertEqual(
                            {tuple(row) for row in features["msa"].tolist()},
                            original_tuples,
                        )


if __name__ == "__main__":
    unittest.main()
