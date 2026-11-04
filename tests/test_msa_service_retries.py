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

"""Offline transport retry checks for both MSA services and every request stage."""

import contextlib
import io
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from protenix.web_service import colab_request_utils as writer


MODES = ("colabfold", "protenix")
STAGES = ("submit", "status", "download")
TRANSPORT_ERRORS = (
    writer.requests.exceptions.ConnectionError,
    writer.requests.exceptions.Timeout,
)


class RetryServer:
    def __init__(self, mode, stage=None, errors=()):
        self.mode = mode
        self.stage = stage
        self.errors = list(errors)
        self.calls = dict.fromkeys(STAGES, 0)
        query = b">query_0\nACDEFG\n"
        files = (
            {"bfd.mgnify30.metaeuk30.smag30.a3m": query, "uniref.a3m": query}
            if mode == "colabfold"
            else {"0.a3m": query, "pdb70_220313_db.m8": b"", "uniref_tax.m8": b""}
        )
        stream = io.BytesIO()
        with tarfile.open(fileobj=stream, mode="w:gz") as handle:
            for name, data in files.items():
                entry = tarfile.TarInfo(name)
                entry.size = len(data)
                handle.addfile(entry, io.BytesIO(data))
        self.archive = stream.getvalue()

    def attempt(self, stage):
        self.calls[stage] += 1
        if stage == self.stage and self.errors:
            raise self.errors.pop(0)

    def post(self, url, **kwargs):
        self.attempt("submit")
        response = Mock()
        response.json.return_value = {
            "status": "RUNNING" if self.stage == "status" else "COMPLETE",
            "id": "synthetic",
        }
        return response

    def get(self, url, **kwargs):
        if url.endswith("/result/download/synthetic"):
            self.attempt("download")
            return Mock(content=self.archive)
        if not url.endswith("/ticket/synthetic"):
            raise AssertionError(f"Unexpected HTTP request: {url}")
        self.attempt("status")
        response = Mock()
        response.json.return_value = {"status": "COMPLETE", "id": "synthetic"}
        return response

    @contextlib.contextmanager
    def active(self):
        with patch.object(writer.requests, "post", side_effect=self.post), patch.object(
            writer.requests, "get", side_effect=self.get
        ), patch.object(writer.time, "sleep") as sleep, contextlib.redirect_stdout(
            io.StringIO()
        ), contextlib.redirect_stderr(
            io.StringIO()
        ):
            yield sleep

    def run(self, root):
        host = (
            "https://offline.invalid"
            if self.mode == "colabfold"
            else "https://protenix-server.com/api/msa"
        )
        return writer.run_mmseqs2_service(
            ">query_0\nACDEFG",
            str(root / "msa"),
            server_mode=self.mode,
            host_url=host,
            user_agent="msa-retry-test",
        )


class TestMSAServiceRetries(unittest.TestCase):
    @staticmethod
    def retry_sleeps(sleep):
        return [call.args[0] for call in sleep.call_args_list if call.args[0] == 5]

    def test_successful_requests_need_no_retries(self):
        for mode in MODES:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                server = RetryServer(mode)
                with server.active() as sleep:
                    result = server.run(Path(directory))
                self.assertEqual(server.calls, dict(submit=1, status=0, download=1))
                self.assertEqual(self.retry_sleeps(sleep), [])
                if mode == "colabfold":
                    self.assertTrue((Path(result) / "non_pairing.a3m").is_file())
                else:
                    self.assertTrue((Path(directory) / "msa/0.a3m").is_file())

    def test_transient_transport_errors_recover(self):
        self.assert_success_after_errors(2)

    def test_success_on_the_last_allowed_attempt(self):
        self.assert_success_after_errors(5)

    def assert_success_after_errors(self, count):
        for mode in MODES:
            for stage in STAGES:
                for error_type in TRANSPORT_ERRORS:
                    with self.subTest(
                        mode=mode, stage=stage, error=error_type.__name__
                    ), tempfile.TemporaryDirectory() as directory:
                        errors = [
                            error_type("synthetic transport error")
                            for _ in range(count)
                        ]
                        server = RetryServer(mode, stage, errors)
                        with server.active() as sleep:
                            server.run(Path(directory))
                        self.assertEqual(server.calls[stage], count + 1)
                        self.assertEqual(self.retry_sleeps(sleep), [5] * count)

    def test_transport_errors_exhaust_the_retry_budget(self):
        for mode in MODES:
            for stage in STAGES:
                for error_type in TRANSPORT_ERRORS:
                    with self.subTest(
                        mode=mode, stage=stage, error=error_type.__name__
                    ), tempfile.TemporaryDirectory() as directory:
                        errors = [
                            error_type("synthetic transport error") for _ in range(6)
                        ]
                        server = RetryServer(mode, stage, errors)
                        with server.active() as sleep, self.assertRaises(
                            error_type
                        ) as raised:
                            server.run(Path(directory))
                        self.assertIs(raised.exception, errors[-1])
                        self.assertEqual(server.calls[stage], 6)
                        self.assertEqual(self.retry_sleeps(sleep), [5] * 5)

    def test_mixed_transport_errors_share_one_budget(self):
        for mode in MODES:
            for stage in STAGES:
                with self.subTest(
                    mode=mode, stage=stage
                ), tempfile.TemporaryDirectory() as directory:
                    errors = [
                        error("synthetic transport error")
                        for _ in range(3)
                        for error in TRANSPORT_ERRORS
                    ]
                    server = RetryServer(mode, stage, errors)
                    with server.active() as sleep, self.assertRaises(
                        writer.requests.exceptions.RequestException
                    ) as raised:
                        server.run(Path(directory))
                    self.assertIs(raised.exception, errors[-1])
                    self.assertEqual(server.calls[stage], 6)
                    self.assertEqual(self.retry_sleeps(sleep), [5] * 5)

    def test_non_transport_errors_propagate_immediately(self):
        for mode in MODES:
            for stage in STAGES:
                with self.subTest(
                    mode=mode, stage=stage
                ), tempfile.TemporaryDirectory() as directory:
                    server = RetryServer(
                        mode, stage, [ValueError("synthetic invalid argument")]
                    )
                    with server.active() as sleep, self.assertRaisesRegex(
                        ValueError, "synthetic invalid argument"
                    ):
                        server.run(Path(directory))
                    self.assertEqual(server.calls[stage], 1)
                    self.assertEqual(self.retry_sleeps(sleep), [])


if __name__ == "__main__":
    unittest.main()
