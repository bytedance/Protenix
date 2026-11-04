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

"""Offline regressions for ColabFold search result directories and JSON paths."""

import contextlib
import io
import json
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from protenix.web_service import colab_request_utils as writer

SEQUENCES = ("ACDEFG", "HIKLMN", "PQRSTV")
PAIRED_HITS = dict(zip(SEQUENCES, ("ACDEFA", "HIKLMA", "PQRSTA")))
UNPAIRED_HITS = dict(zip(SEQUENCES, ("ACDEFC", "HIKLMC", "PQRSTC")))


def fasta(rows):
    return "".join(f">{header}\n{sequence}\n" for header, sequence in rows)


class ColabFoldServer:
    """Supply synthetic archives to the real writer through mocked HTTP."""

    def __init__(self, omit_paired_chain=None):
        self.omit_paired_chain = omit_paired_chain
        self.submissions = []
        self.archives = {}

    def post(self, url, **kwargs):
        rows = list(writer.parse_fasta_string(kwargs["data"]["q"]).items())
        paired = url.endswith("/ticket/pair")
        self.submissions.append(paired)
        if paired:
            data = "".join(
                fasta([row, (f"UniRef100_pair{i}", PAIRED_HITS[row[1]])]) + "\0"
                for i, row in enumerate(rows)
                if i != self.omit_paired_chain
            )
            files = {"pair.a3m": data.encode()}
        else:
            self.assert_single_query(rows)
            row = rows[0]
            data = fasta([row, ("UniRef100_unpaired", UNPAIRED_HITS[row[1]])])
            files = {
                "bfd.mgnify30.metaeuk30.smag30.a3m": data.encode(),
                "uniref.a3m": data.encode(),
            }
        archive = io.BytesIO()
        with tarfile.open(fileobj=archive, mode="w:gz") as handle:
            for name, data in files.items():
                entry = tarfile.TarInfo(name)
                entry.size = len(data)
                handle.addfile(entry, io.BytesIO(data))
        identifier = str(len(self.archives))
        self.archives[identifier] = archive.getvalue()
        response = Mock()
        response.json.return_value = {"status": "COMPLETE", "id": identifier}
        return response

    @staticmethod
    def assert_single_query(rows):
        if len(rows) != 1:
            raise AssertionError("The unpaired fixture expects one query per search")

    def get(self, url, **kwargs):
        if "/result/download/" not in url:
            raise AssertionError(f"Unexpected HTTP request: {url}")
        response = Mock()
        response.content = self.archives[url.rsplit("/", 1)[1]]
        return response

    @contextlib.contextmanager
    def activate(self):
        with patch.object(writer.requests, "post", side_effect=self.post), patch.object(
            writer.requests, "get", side_effect=self.get
        ), patch.object(
            writer.time, "sleep", side_effect=AssertionError("Unexpected retry")
        ), contextlib.redirect_stdout(
            io.StringIO()
        ), contextlib.redirect_stderr(
            io.StringIO()
        ):
            yield self


class TestColabFoldCachedDirectory(unittest.TestCase):
    def test_cached_unpaired_search_returns_the_original_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "msa"
            server = ColabFoldServer()
            with server.activate():
                first = writer.run_mmseqs2_service(
                    ">query_2\nPQRSTV",
                    str(prefix),
                    server_mode="colabfold",
                    host_url="https://offline.invalid",
                    user_agent="msa-handoff-test",
                )
                second = writer.run_mmseqs2_service(
                    ">query_2\nPQRSTV",
                    str(prefix),
                    server_mode="colabfold",
                    host_url="https://offline.invalid",
                    user_agent="msa-handoff-test",
                )
            self.assertEqual(second, first)
            self.assertEqual(Path(first), prefix / "2")
            self.assertEqual(server.submissions, [False])

    def test_cached_archive_without_unpaired_output_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory)
            (prefix / "out.tar.gz").write_bytes(b"cached archive")
            with ColabFoldServer().activate(), self.assertRaises(FileNotFoundError):
                writer.run_mmseqs2_service(
                    ">query_0\nACDEFG",
                    str(prefix),
                    server_mode="colabfold",
                    host_url="https://offline.invalid",
                    user_agent="msa-handoff-test",
                )

    def test_cached_unpaired_query_mismatch_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "msa"
            server = ColabFoldServer()
            with server.activate():
                writer.run_mmseqs2_service(
                    ">query_0\nACDEFG",
                    str(prefix),
                    server_mode="colabfold",
                    host_url="https://offline.invalid",
                    user_agent="msa-handoff-test",
                )
                with self.assertRaisesRegex(ValueError, "query"):
                    writer.run_mmseqs2_service(
                        ">query_0\nHIKLMN",
                        str(prefix),
                        server_mode="colabfold",
                        host_url="https://offline.invalid",
                        user_agent="msa-handoff-test",
                    )
            self.assertEqual(server.submissions, [False])


class TestColabFoldHandoff(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            from protenix.web_service import colab_request_parser
            from runner import msa_search
        except ImportError as exc:
            raise unittest.SkipTest(f"Search dependencies unavailable: {exc}") from exc
        cls.search = msa_search
        cls.parser = colab_request_parser

    def search_directories(self, count, root):
        return self.search.msa_search(SEQUENCES[:count], str(root), mode="colabfold")

    def assert_chain_files(self, directories, sequences):
        self.assertEqual(len(directories), len(sequences))
        for result_dir, sequence in zip(directories, sequences):
            directory = Path(result_dir)
            for name, hit in [
                ("pairing.a3m", PAIRED_HITS[sequence]),
                ("non_pairing.a3m", UNPAIRED_HITS[sequence]),
            ]:
                path = directory / name
                self.assertTrue(path.is_file(), str(path))
                records = writer.parse_fasta_string(path.read_text())
                self.assertEqual(records["query"], sequence)
                self.assertIn(hit, records.values())

    def test_returned_directories_contain_both_msas_for_each_chain(self):
        for count in [2, 3]:
            with self.subTest(chains=count), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "msa"
                with ColabFoldServer().activate():
                    results = self.search_directories(count, root)
                self.assert_chain_files(results, SEQUENCES[:count])
                for i, result_dir in enumerate(results):
                    self.assertEqual(
                        (Path(result_dir) / "pairing.a3m").read_bytes(),
                        (root / "complex" / str(i) / "pairing.a3m").read_bytes(),
                    )

    def test_single_chain_placeholder_is_in_the_returned_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            server = ColabFoldServer()
            with server.activate():
                results = self.search_directories(1, Path(directory) / "msa")
            self.assertEqual(server.submissions, [False])
            self.assertTrue((Path(results[0]) / "pairing.a3m").is_file())
            self.assertEqual(
                (Path(results[0]) / "pairing.a3m").read_text(), ">query\nACDEFG"
            )

    def test_repeated_search_reuses_cached_results_and_preserves_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "msa"
            server = ColabFoldServer()
            with server.activate():
                first = self.search_directories(3, root)
                for result_dir in first:
                    (Path(result_dir) / "pairing.a3m").unlink(missing_ok=True)
                second = self.search_directories(3, root)
            self.assertEqual(second, first)
            self.assertEqual(server.submissions, [False, False, False, True])
            self.assert_chain_files(second, SEQUENCES)

    def test_generated_pairs_replace_old_files_in_returned_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "msa"
            for i, sequence in enumerate(SEQUENCES):
                target = root / str(i) / str(i) / "pairing.a3m"
                target.parent.mkdir(parents=True)
                target.write_text(fasta([("query", sequence), ("stale", sequence)]))
            with ColabFoldServer().activate():
                results = self.search_directories(3, root)
            self.assert_chain_files(results, SEQUENCES)
            for result_dir in results:
                self.assertNotIn(
                    "stale", (Path(result_dir) / "pairing.a3m").read_text()
                )

    def test_search_failures_propagate_before_json_is_written(self):
        for paired_failure in [False, True]:
            with self.subTest(
                paired=paired_failure
            ), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                input_path = self.write_input(root, 3)
                server = ColabFoldServer()
                original_run = writer.run_mmseqs2_service

                def search(*args, **kwargs):
                    if kwargs.get("use_pairing", False) == paired_failure:
                        raise RuntimeError("synthetic search failure")
                    return original_run(*args, **kwargs)

                with server.activate(), patch.object(
                    self.parser, "run_mmseqs2_service", side_effect=search
                ), self.assertRaisesRegex(RuntimeError, "synthetic search failure"):
                    self.search.update_infer_json(
                        str(input_path), str(root / "output"), mode="colabfold"
                    )
                self.assertFalse((root / "input-update-msa.json").exists())

    def test_missing_chain_in_paired_response_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            with ColabFoldServer(omit_paired_chain=1).activate(), self.assertRaises(
                FileNotFoundError
            ):
                self.search_directories(3, Path(directory) / "msa")

    def test_retry_exhaustion_propagates_before_json_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = self.write_input(root, 3)
            server = ColabFoldServer()
            paired_attempts = 0

            def post(url, **kwargs):
                nonlocal paired_attempts
                if url.endswith("/ticket/pair"):
                    paired_attempts += 1
                    if paired_attempts <= 6:
                        raise writer.requests.exceptions.ConnectionError(
                            "synthetic paired connection error"
                        )
                return server.post(url, **kwargs)

            with server.activate(), patch.object(
                writer.requests, "post", side_effect=post
            ), patch.object(writer.time, "sleep"), self.assertRaisesRegex(
                writer.requests.exceptions.ConnectionError,
                "synthetic paired connection error",
            ):
                self.search.update_infer_json(
                    str(input_path), str(root / "output"), mode="colabfold"
                )
            self.assertEqual(paired_attempts, 6)
            self.assertFalse((root / "input-update-msa.json").exists())

    def test_missing_later_pair_does_not_replace_earlier_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "msa"
            old_outputs = {}
            for i, sequence in enumerate(SEQUENCES):
                path = root / str(i) / str(i) / "pairing.a3m"
                path.parent.mkdir(parents=True)
                path.write_text(fasta([("query", sequence), ("stale", sequence)]))
                old_outputs[path] = path.read_bytes()
            with ColabFoldServer(omit_paired_chain=2).activate(), self.assertRaises(
                FileNotFoundError
            ):
                self.search_directories(3, root)
            for path, content in old_outputs.items():
                self.assertEqual(path.read_bytes(), content)

    def test_cached_paired_query_mismatch_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "msa"
            server = ColabFoldServer()
            with server.activate():
                self.search_directories(2, root)
                for i in range(2):
                    shutil.rmtree(root / str(i))
                with self.assertRaisesRegex(ValueError, "query"):
                    self.search.msa_search(SEQUENCES[1:], str(root), mode="colabfold")
            self.assertEqual(server.submissions, [False, False, True, False, False])

    def test_cached_paired_query_set_change_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "msa"
            server = ColabFoldServer()
            with server.activate():
                results = self.search_directories(3, root)
                self.assert_chain_files(results, SEQUENCES)
                previous = {
                    Path(result)
                    / "pairing.a3m": (Path(result) / "pairing.a3m").read_bytes()
                    for result in results
                }
                with self.assertRaisesRegex(ValueError, "query"):
                    self.search_directories(2, root)
            self.assertEqual(server.submissions, [False, False, False, True])
            for path, content in previous.items():
                self.assertEqual(path.read_bytes(), content)

    def test_truncated_cached_paired_output_is_reported(self):
        for content in ["", ">query\n"]:
            with self.subTest(
                content=content
            ), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "msa"
                server = ColabFoldServer()
                with server.activate():
                    self.search_directories(2, root)
                    (root / "complex" / "1" / "pairing.a3m").write_text(content)
                    with self.assertRaisesRegex(ValueError, "query"):
                        self.search_directories(2, root)

    @staticmethod
    def write_input(root, count):
        path = root / "input.json"
        path.write_text(
            json.dumps(
                [
                    {
                        "name": "synthetic",
                        "sequences": [
                            {"proteinChain": {"sequence": sequence, "count": 1}}
                            for sequence in reversed(SEQUENCES[:count])
                        ],
                    }
                ]
            )
        )
        return path

    def test_json_paths_reach_real_inference_features(self):
        import numpy as np
        from biotite.structure import AtomArray
        from protenix.data.constants import PROTEIN_CHAIN
        from protenix.data.msa.msa_featurizer import InferenceMSAFeaturizer
        from protenix.data.msa.msa_utils import MSACore

        for count in [2, 3]:
            with self.subTest(chains=count), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                input_path = self.write_input(root, count)
                with ColabFoldServer().activate():
                    updated_path, searched = self.search.update_infer_json(
                        str(input_path), str(root / "output"), mode="colabfold"
                    )
                self.assertTrue(searched)
                task = json.loads(Path(updated_path).read_text())[0]
                chains = [entity["proteinChain"] for entity in task["sequences"]]
                for chain in chains:
                    self.assertIn("pairedMsaPath", chain)
                    self.assertIn("unpairedMsaPath", chain)
                    self.assert_chain_files(
                        [str(Path(chain["pairedMsaPath"]).parent)], [chain["sequence"]]
                    )
                self.assertFalse(self.search.need_msa_search(task))
                atoms = AtomArray(6 * count)
                atoms.set_annotation("asym_id_int", np.repeat(np.arange(count), 6))
                atoms.set_annotation("centre_atom_mask", np.ones(6 * count, dtype=bool))
                atoms.chain_id = np.repeat(list("ABC"[:count]), 6)
                atoms.res_id = np.tile(np.arange(1, 7), count)
                expected, _ = MSACore.sequences_to_array(
                    ["".join(PAIRED_HITS[chain["sequence"]] for chain in chains)],
                    PROTEIN_CHAIN,
                )
                for flag in [False, True]:
                    with self.subTest(msa_pair_as_unpair=flag):
                        features = InferenceMSAFeaturizer.make_msa_feature(
                            task["sequences"],
                            atoms,
                            msa_pair_as_unpair=flag,
                            use_rna_msa=False,
                        )
                        self.assertTrue(
                            np.any(np.all(features["msa"] == expected[0], axis=1))
                        )

    def test_valid_precomputed_paths_skip_search(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = self.write_input(root, 3)
            tasks = json.loads(input_path.read_text())
            for i, entity in enumerate(tasks[0]["sequences"]):
                for key in ["pairedMsaPath", "unpairedMsaPath"]:
                    path = root / f"provided_{i}_{key}.a3m"
                    path.write_text(
                        fasta([("query", entity["proteinChain"]["sequence"])])
                    )
                    entity["proteinChain"][key] = str(path)
            input_path.write_text(json.dumps(tasks))
            server = ColabFoldServer()
            with server.activate():
                result, searched = self.search.update_infer_json(
                    str(input_path), str(root / "output"), mode="colabfold"
                )
            self.assertFalse(searched)
            self.assertEqual(result, str(input_path))
            self.assertEqual(server.submissions, [])
            self.assertEqual(json.loads(Path(result).read_text()), tasks)

    def test_protenix_mode_keeps_its_postprocessing_path(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(
                self.parser.RequestParser, "msa_search"
            ) as search, patch.object(
                self.parser.RequestParser, "msa_postprocess", return_value=["processed"]
            ) as postprocess:
                results = self.search.msa_search(
                    [SEQUENCES[0]], directory, mode="protenix"
                )
            self.assertEqual(results, ["processed"])
            self.assertEqual(search.call_args.kwargs["mode"], "protenix")
            postprocess.assert_called_once_with(
                seqs_pending_msa=[SEQUENCES[0]], msa_res_dir=directory
            )


if __name__ == "__main__":
    unittest.main()
