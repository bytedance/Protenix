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

import os
import subprocess
import time
from contextlib import nullcontext
from typing import Sequence

from protenix.data.tools.common import parse_kalign_a3m, tmpdir_manager, to_a3m
from protenix.utils.logger import get_logger

logger = get_logger(__name__)


class Kalign:
    """
    Python wrapper of the Kalign binary.
    Adapted from openfold.data.tools.kalign

    Args:
      binary_path: The path to the Kalign binary.

    Raises:
      RuntimeError: If Kalign binary not found within the path.
    """

    def __init__(self, *, binary_path: str):
        # First check if kalign is available in the system PATH using 'which'
        found_path = None
        try:
            result = subprocess.run(["which", "kalign"], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                kalign_in_path = result.stdout.strip()
                if os.path.exists(kalign_in_path) and os.access(
                    kalign_in_path, os.X_OK
                ):
                    found_path = kalign_in_path
        except (subprocess.CalledProcessError, FileNotFoundError):
            # 'which' command not found or failed, continue to check other options
            pass

        # If not found in PATH, check the passed binary_path
        if found_path is None:
            if os.path.exists(binary_path) and os.access(binary_path, os.X_OK):
                found_path = binary_path
            else:
                # Generate error message with installation instructions
                raise RuntimeError(
                    f"Kalign binary not found. Neither kalign in system PATH "
                    f"nor provided path ({binary_path}) exist and are executable.\n"
                    f"To install kalign, you can use one of the following methods:\n"
                    f"1. Using conda: conda install -c bioconda kalign\n"
                    f"2. Using apt (Ubuntu/Debian): apt-get install kalign\n"
                    f"3. Download from: https://github.com/TimoLassmann/kalign\n"
                    f"After installation, make sure the binary is accessible either in PATH or at the specified location."
                )

        # Set the found path as the binary path to use
        self.binary_path = found_path

    def align(self, sequences: Sequence[str]) -> Sequence[str]:
        """Aligns the sequences and returns the alignment in A3M string.

        Args:
          sequences: A list of query sequence strings. The sequences have to be at
            least 6 residues long (Kalign requires this). Note that the order in
            which you give the sequences might alter the output slightly as
            different alignment tree might get constructed.

        Returns:
          A list of strings with the aligned sequences in a3m format.

        Raises:
          RuntimeError: If Kalign fails.
          ValueError: If any of the sequences is less than 6 residues long.
        """

        for s in sequences:
            if len(s) < 6:
                raise ValueError(
                    "Kalign requires all sequences to be at least 6 "
                    f"residues long. Got {s} ({len(s)} residues)."
                )

        with tmpdir_manager() as query_tmp_dir:
            input_fasta_path = os.path.join(query_tmp_dir, "input.fasta")
            output_a3m_path = os.path.join(query_tmp_dir, "output.a3m")

            logger.debug(
                "Kalign tmpdir=%s input_fasta=%s output_a3m=%s",
                query_tmp_dir,
                input_fasta_path,
                output_a3m_path,
            )

            write_st = time.perf_counter()
            with open(input_fasta_path, "w") as f:
                f.write(to_a3m(sequences))
            write_seconds = time.perf_counter() - write_st

            cmd = [
                self.binary_path,
                "-i",
                input_fasta_path,
                "-o",
                output_a3m_path,
                "-format",
                "fasta",
            ]

            popen_st = time.perf_counter()
            process = subprocess.Popen(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
            )

            with nullcontext():
                stdout, stderr = process.communicate()
                retcode = process.wait()
            popen_seconds = time.perf_counter() - popen_st

            if retcode:
                logger.debug(
                    "Kalign timings write=%.4fs run=%.4fs read=%.4fs",
                    write_seconds,
                    popen_seconds,
                    0.0,
                )
                raise RuntimeError(
                    f"Kalign failed\nstdout:\n{stdout.decode('utf-8')}\n\nstderr:\n{stderr.decode('utf-8')}\n"
                )

            read_st = time.perf_counter()
            with open(output_a3m_path) as f:
                a3m = f.read()
            read_seconds = time.perf_counter() - read_st

            logger.debug(
                "Kalign timings write=%.4fs run=%.4fs read=%.4fs",
                write_seconds,
                popen_seconds,
                read_seconds,
            )
            return parse_kalign_a3m(a3m)
