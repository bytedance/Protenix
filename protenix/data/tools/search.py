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

import dataclasses
import datetime
import os
import pathlib
import re
import subprocess
import tempfile
import time
from typing import Any, Final, List, Mapping, Optional, Protocol, Union

from protenix.data.constants import STANDARD_POLYMER_CHAIN_TYPES
from protenix.data.tools.common import (
    convert_a3m_to_stockholm,
    convert_stockholm_to_a3m,
    lazy_fasta_parse,
)
from protenix.utils.logger import get_logger

logger = get_logger(__name__)

SHORT_SEQ_LIMIT: Final[int] = 50


def run_shell(
    cmd: List[str],
    name: str,
    **kwargs: Any,
) -> subprocess.CompletedProcess:
    """
    Run a shell command and log its execution and results.

    Args:
        cmd: List of command arguments.
        name: Logical name for the command (for logging).
        **kwargs: Additional keyword arguments for subprocess.run.

    Returns:
        CompletedProcess object.

    Raises:
        RuntimeError: If the command fails.
    """
    logger.info(f"Executing {name}: {' '.join(cmd)}")
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True, **kwargs)
        logger.info(f"{name} finished in {time.time()-t0:.3f}s")
        return proc
    except subprocess.CalledProcessError as e:
        logger.error(f"{name} failed. Stderr: {e.stderr}")
        raise RuntimeError(
            f"{name} failed\nStdout: {e.stdout}\nStderr: {e.stderr}"
        ) from e


class BinaryWrapper:
    """
    Base class for wrapping external bioinformatics binaries.

    Args:
        binary_path: Path to the executable.
        name: Name of the tool.

    Raises:
        RuntimeError: If the binary path does not exist.
    """

    def __init__(self, binary_path: str, name: str) -> None:
        if not os.path.exists(binary_path):
            raise RuntimeError(f"{name} not found at {binary_path}")
        self.path = binary_path


class Hmmalign(BinaryWrapper):
    """
    Wrapper for HMMER hmmalign.

    Args:
        binary_path: Path to the hmmalign executable.
    """

    def __init__(self, binary_path: str) -> None:
        super().__init__(binary_path, "hmmalign")

    def align(
        self, a3m: str, profile: str, flags: Optional[Mapping[str, str]] = None
    ) -> str:
        """
        Align sequences in A3M format to a HMM profile.

        Args:
            a3m: Input sequences in A3M format.
            profile: HMM profile string.
            flags: Optional additional command line flags.

        Returns:
            Aligned sequences in A3M format.
        """
        with tempfile.TemporaryDirectory() as tmp:
            prof_p, in_p, out_p = f"{tmp}/p", f"{tmp}/i", f"{tmp}/o"
            pathlib.Path(prof_p).write_text(profile)
            pathlib.Path(in_p).write_text(a3m)
            cmd = [self.path, "-o", out_p, "--outformat", "A2M"]
            if flags:
                for k, v in flags.items():
                    cmd.extend([k, v])
            run_shell(cmd + [prof_p, in_p], "hmmalign")
            return pathlib.Path(out_p).read_text()

    def align_sequences_to_profile(self, profile: str, a3m: str) -> str:
        """
        Clean gaps and realign sequences to a profile.

        Args:
            profile: HMM profile string.
            a3m: Input sequences in A3M format.

        Returns:
            Realigned sequences.
        """
        clean_a3m = "\n".join(
            f">{d}\n{s.replace('-', '')}" for s, d in lazy_fasta_parse(a3m)
        )
        return self.align(clean_a3m, profile)


class Hmmbuild(BinaryWrapper):
    """
    Wrapper for HMMER hmmbuild.

    Args:
        binary_path: Path to the hmmbuild executable.
        singlemx: Whether to use --singlemx flag.
        alphabet: Optional alphabet specification (e.g., 'amino', 'rna', 'dna').
    """

    def __init__(
        self, binary_path: str, singlemx: bool = False, alphabet: Optional[str] = None
    ) -> None:
        super().__init__(binary_path, "hmmbuild")
        self.opts = (["--singlemx"] if singlemx else []) + (
            [f"--{alphabet}"] if alphabet else []
        )

    def build(self, msa: str, informat: str, model: str = "fast") -> str:
        """
        Build a HMM profile from a Multiple Sequence Alignment.

        Args:
            msa: Input MSA string.
            informat: Format of the input MSA (e.g., 'afa', 'stockholm').
            model: Model construction strategy ('fast' or 'hand').

        Returns:
            The HMM profile string.
        """
        with tempfile.TemporaryDirectory() as tmp:
            in_p, out_p = f"{tmp}/i", f"{tmp}/o"
            pathlib.Path(in_p).write_text(msa)
            cmd = [self.path, "--informat", informat] + self.opts
            if model == "hand":
                cmd.append("--hand")
            run_shell(cmd + [out_p, in_p], "hmmbuild")
            return pathlib.Path(out_p).read_text()

    def build_profile_from_a3m(self, a3m: str) -> str:
        """
        Build HMM from A3M format by removing insertions.

        Args:
            a3m: Input A3M string.

        Returns:
            HMM profile string.
        """
        clean = "".join(
            f">{d}\n{re.sub('[a-z]+', '', s)}\n" for s, d in lazy_fasta_parse(a3m)
        )
        return self.build(clean, "afa")

    def build_profile_from_sto(self, sto: str, model: str = "fast") -> str:
        """
        Build HMM from Stockholm format.

        Args:
            sto: Input Stockholm string.
            model: Model construction strategy.

        Returns:
            HMM profile string.
        """
        return self.build(sto, "stockholm", model)


class Hmmsearch(BinaryWrapper):
    """
    Wrapper for HMMER hmmsearch.

    Args:
        binary_path: Path to hmmsearch executable.
        hmmbuild_binary_path: Path to hmmbuild executable.
        database_path: Path to the sequence database to search.
        alphabet: Alphabet type.
        **filters: Search filters (e.g., e_value, inc_e, filter_f1, etc.).
    """

    def __init__(
        self,
        binary_path: str,
        hmmbuild_binary_path: str,
        database_path: str,
        alphabet: str = "amino",
        **filters: Any,
    ) -> None:
        super().__init__(binary_path, "hmmsearch")
        self.builder = Hmmbuild(binary_path=hmmbuild_binary_path, alphabet=alphabet)
        self.db = database_path
        self.flags = ["--max"] if filters.get("filter_max") else []
        if not self.flags:
            map_f = {
                "filter_f1": "--F1",
                "filter_f2": "--F2",
                "filter_f3": "--F3",
                "e_value": "-E",
                "inc_e": "--incE",
                "dom_e": "--domE",
                "incdom_e": "--incdomE",
            }
            for k, v in filters.items():
                if k in map_f and v is not None:
                    self.flags.extend([map_f[k], str(v)])

    def query_with_hmm(self, hmm: str) -> str:
        """
        Search a database using an HMM profile.

        Args:
            hmm: HMM profile string.

        Returns:
            Search result in A3M format.
        """
        with tempfile.TemporaryDirectory() as tmp:
            hmm_p, sto_p = f"{tmp}/q.hmm", f"{tmp}/o.sto"
            pathlib.Path(hmm_p).write_text(hmm)
            cmd = (
                [self.path, "--noali", "--cpu", "8"]
                + self.flags
                + ["-A", sto_p, hmm_p, self.db]
            )
            run_shell(cmd, "hmmsearch")
            with open(sto_p) as f:
                return convert_stockholm_to_a3m(f, remove_gaps=False, linewidth=60)

    def query_with_sto(self, sto: str, model: str = "fast") -> str:
        """
        Build an HMM from a Stockholm MSA and search the database.

        Args:
            sto: Input Stockholm MSA string.
            model: Model construction strategy.

        Returns:
            Search result in A3M format.
        """
        return self.query_with_hmm(self.builder.build_profile_from_sto(sto, model))


@dataclasses.dataclass(frozen=True, slots=True)
class MsaToolResult:
    """Data class for MSA tool results."""

    target_sequence: str
    e_value: Optional[float]
    a3m: str


class MsaTool(Protocol):
    """Protocol for MSA search tools."""

    def query(self, target_sequence: str) -> MsaToolResult:
        """Perform a sequence search query."""
        ...


class Jackhmmer(BinaryWrapper, MsaTool):
    """
    Wrapper for HMMER jackhmmer.

    Args:
        binary_path: Path to jackhmmer executable.
        database_path: Path to search database.
        n_cpu: Number of CPUs to use.
        n_iter: Number of iterations.
        e_value: E-value threshold.
        **kwargs: Additional parameters (e.g., z_value, filter_f1, etc.).
    """

    def __init__(
        self,
        binary_path: str,
        database_path: str,
        n_cpu: int = 8,
        n_iter: int = 3,
        e_value: Optional[float] = 1e-3,
        **kwargs: Any,
    ) -> None:
        super().__init__(binary_path, "jackhmmer")
        self.db, self.cpu, self.iter, self.e = database_path, n_cpu, n_iter, e_value
        self.extra = kwargs

    def query(self, seq: str) -> MsaToolResult:
        """
        Query a database using jackhmmer for a given sequence.

        Args:
            seq: Input sequence.

        Returns:
            MsaToolResult object.
        """
        with tempfile.TemporaryDirectory() as tmp:
            fa_p, sto_p = f"{tmp}/q.fa", f"{tmp}/o.sto"
            pathlib.Path(fa_p).write_text(f">query\n{seq}\n")
            cmd = [
                self.path,
                "-o",
                "/dev/null",
                "-A",
                sto_p,
                "--noali",
                "--cpu",
                str(self.cpu),
                "-N",
                str(self.iter),
            ]
            if self.e is not None:
                cmd.extend(["-E", str(self.e), "--incE", str(self.e)])
            for k, v in self.extra.items():
                flag = {
                    "filter_f1": "--F1",
                    "filter_f2": "--F2",
                    "filter_f3": "--F3",
                    "z_value": "-Z",
                }.get(k)
                if flag and v is not None:
                    cmd.extend([flag, str(v)])
            run_shell(cmd + [fa_p, self.db], "jackhmmer")
            with open(sto_p) as f:
                a3m = convert_stockholm_to_a3m(
                    f, max_seqs=self.extra.get("max_sequences", 5000)
                )
            return MsaToolResult(target_sequence=seq, e_value=self.e, a3m=a3m)


class Nhmmer(BinaryWrapper, MsaTool):
    """
    Wrapper for HMMER nhmmer.

    Args:
        binary_path: Path to nhmmer executable.
        hmmalign_binary_path: Path to hmmalign executable.
        hmmbuild_binary_path: Path to hmmbuild executable.
        database_path: Path to search database.
        n_cpu: Number of CPUs to use.
        e_value: E-value threshold.
        **kwargs: Additional parameters (e.g., alphabet, strand, filter_f3, etc.).
    """

    def __init__(
        self,
        binary_path: str,
        hmmalign_binary_path: str,
        hmmbuild_binary_path: str,
        database_path: str,
        n_cpu: int = 8,
        e_value: float = 1e-3,
        **kwargs: Any,
    ) -> None:
        super().__init__(binary_path, "nhmmer")
        self.aligner_p, self.builder_p, self.db, self.cpu, self.e = (
            hmmalign_binary_path,
            hmmbuild_binary_path,
            database_path,
            n_cpu,
            e_value,
        )
        self.opts = kwargs

    def query(self, seq: str) -> MsaToolResult:
        """
        Query a database using nhmmer for a given sequence.

        Args:
            seq: Input sequence.

        Returns:
            MsaToolResult object.
        """
        with tempfile.TemporaryDirectory() as tmp:
            fa_p, sto_p = f"{tmp}/q.fa", f"{tmp}/o.sto"
            pathlib.Path(fa_p).write_text(f">query\n{seq}\n")
            cmd = [
                self.path,
                "-o",
                "/dev/null",
                "--noali",
                "--cpu",
                str(self.cpu),
                "-E",
                str(self.e),
            ]
            if self.opts.get("alphabet"):
                cmd.append(f"--{self.opts['alphabet']}")
            if self.opts.get("strand"):
                cmd.append(f"--{self.opts['strand']}")
            f3 = (
                0.02
                if self.opts.get("alphabet") == "rna" and len(seq) < SHORT_SEQ_LIMIT
                else self.opts.get("filter_f3", 1e-5)
            )
            cmd.extend(["-A", sto_p, "--F3", str(f3), fa_p, self.db])
            run_shell(cmd, "nhmmer")

            if os.path.getsize(sto_p) > 0:
                with open(sto_p) as f:
                    a3m = convert_stockholm_to_a3m(
                        f, max_seqs=self.opts.get("max_sequences", 5000) - 1
                    )
                prof = Hmmbuild(
                    binary_path=self.builder_p, alphabet=self.opts.get("alphabet")
                ).build_profile_from_a3m(f">query\n{seq}\n")
                aligned = Hmmalign(
                    binary_path=self.aligner_p
                ).align_sequences_to_profile(prof, a3m)
                final_a3m = "\n".join(
                    f">{d}\n{s}"
                    for s, d in lazy_fasta_parse(f">query\n{seq}\n" + aligned)
                )
            else:
                final_a3m = f">query\n{seq}"
        return MsaToolResult(target_sequence=seq, e_value=self.e, a3m=final_a3m)


# Configuration Dataclasses
@dataclasses.dataclass(frozen=True, slots=True)
class DatabaseConfig:
    """Config for a biological sequence database."""

    name: str
    path: str


@dataclasses.dataclass(frozen=True, slots=True)
class JackhmmerConfig:
    """Config for a jackhmmer search run."""

    binary_path: str
    database_config: DatabaseConfig
    n_cpu: int
    n_iter: int
    e_value: float
    z_value: Optional[Union[float, int]]
    max_sequences: int


@dataclasses.dataclass(frozen=True, slots=True)
class NhmmerConfig:
    """Config for a nhmmer search run."""

    binary_path: str
    hmmalign_binary_path: str
    hmmbuild_binary_path: str
    database_config: DatabaseConfig
    n_cpu: int
    e_value: float
    max_sequences: int
    alphabet: Optional[str]


@dataclasses.dataclass(frozen=True, slots=True)
class RunConfig:
    """Generic config for an MSA search run."""

    config: Union[JackhmmerConfig, NhmmerConfig]
    chain_poly_type: str
    crop_size: Optional[int]

    def __post_init__(self) -> None:
        """Validate run configuration."""
        if self.crop_size is not None and self.crop_size < 2:
            raise ValueError("crop_size < 2")
        if self.chain_poly_type not in STANDARD_POLYMER_CHAIN_TYPES:
            raise ValueError("Invalid chain type")


@dataclasses.dataclass(frozen=True, slots=True)
class HmmsearchConfig:
    """Config for an hmmsearch run."""

    hmmsearch_binary_path: str
    hmmbuild_binary_path: str
    e_value: float = 100
    inc_e: float = 100
    dom_e: float = 100
    incdom_e: float = 100
    alphabet: str = "amino"
    filter_f1: Optional[float] = 0.1
    filter_f2: Optional[float] = 0.1
    filter_f3: Optional[float] = 0.1
    filter_max: bool = False


def run_hmmsearch_with_a3m(
    database_path: Union[str, os.PathLike],
    hmmsearch_config: HmmsearchConfig,
    max_a3m_query_sequences: Optional[int],
    a3m: Optional[str],
) -> str:
    """
    Search database using hmmsearch with an input A3M string.

    Args:
        database_path: Path to search database.
        hmmsearch_config: HmmsearchConfig object.
        max_a3m_query_sequences: Max sequences to use for building the HMM.
        a3m: Input A3M format string.

    Returns:
        Search results in A3M format.
    """
    cfg = dataclasses.asdict(hmmsearch_config)
    search_bin = cfg.pop("hmmsearch_binary_path")
    build_bin = cfg.pop("hmmbuild_binary_path")
    searcher = Hmmsearch(
        binary_path=search_bin,
        hmmbuild_binary_path=build_bin,
        database_path=database_path,
        **cfg,
    )
    sto = convert_a3m_to_stockholm(a3m, max_a3m_query_sequences)
    return searcher.query_with_sto(sto, model="hand")


# Legacy compatibility / Other classes if needed
@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class TemplateToolConfig:
    """Configuration for template search tool."""

    database_path: str
    chain_poly_type: str
    hmmsearch_config: HmmsearchConfig
    max_a3m_query_sequences: Optional[int] = 300

    def __post_init__(self) -> None:
        """Validate template tool configuration."""
        if self.chain_poly_type not in STANDARD_POLYMER_CHAIN_TYPES:
            raise ValueError("Invalid chain type")


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class TemplateFilterConfig:
    """Configuration for filtering template search results."""

    max_subsequence_ratio: Optional[float]
    min_align_ratio: Optional[float]
    min_hit_length: Optional[int]
    deduplicate_sequences: bool
    max_hits: Optional[int]
    max_template_date: datetime.date

    @classmethod
    def no_op_filter(cls) -> "TemplateFilterConfig":
        """Return a filter config that doesn't filter anything."""
        return cls(
            max_subsequence_ratio=None,
            min_align_ratio=None,
            min_hit_length=None,
            deduplicate_sequences=False,
            max_hits=None,
            max_template_date=datetime.date(3000, 1, 1),
        )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class TemplatesConfig:
    """Complete configuration for the template search pipeline."""

    template_tool_config: TemplateToolConfig
    filter_config: TemplateFilterConfig
