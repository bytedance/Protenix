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

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass
class LocalColabFoldConfig:
    """Configuration for ColabFold search."""

    colabsearch: str
    query_fpath: str
    db_dir: str
    results_dir: str
    mmseqs_path: Optional[str] = None
    db1: str = "uniref30_2302_db"
    db2: Optional[str] = None
    db3: Optional[str] = "colabfold_envdb_202108_db"
    use_env: int = 1
    filter: int = 1
    db_load_mode: int = 0


class A3MProcessor:
    """Processor for A3M file format."""

    def __init__(self, a3m_file: str, out_dir: str):
        self.out_dir = out_dir
        self.a3m_file = Path(a3m_file)
        self.a3m_content = self._read_a3m_file()
        self.chain_info = self._parse_header()

    def _read_a3m_file(self) -> str:
        """Read A3M file content."""
        return self.a3m_file.read_text()

    def _parse_header(self) -> Tuple[List[str], Dict[str, Tuple[int, int]]]:
        """Parse A3M header to get chain information."""
        first_line = self.a3m_content.split("\n")[0]
        if first_line[0] == "#":
            lengths, oligomeric_state = first_line.split("\t")

            chain_lengths = [int(x) for x in lengths[1:].split(",")]
            chain_names = [f"10{x+1}" for x in range(len(oligomeric_state.split(",")))]

            # Calculate sequence ranges for each chain
            seq_ranges = {}
            for i, name in enumerate(chain_names):
                start = sum(chain_lengths[:i])
                end = sum(chain_lengths[: i + 1])
                seq_ranges[name] = (start, end)

            return chain_names, seq_ranges
        else:
            non_pairing = ">query\n" + "\n".join(self.a3m_content.split("\n")[1:])
            query_seq = self.a3m_content.split("\n")[1]
            pairing = f">query\n{query_seq}"
            msa_path = Path(self.out_dir) / "msa"
            msa_path.mkdir(exist_ok=True)
            msa_path = msa_path / "0"
            msa_path.mkdir(exist_ok=True)
            with open(msa_path / "non_pairing.a3m", "w") as f:
                f.write(non_pairing)

            with open(msa_path / "pairing.a3m", "w") as f:
                f.write(pairing)

            return [None]

    def _extract_sequence(self, line: str, range_tuple: Tuple[int, int]) -> str:
        """Extract sequence for specific range."""
        seq = []
        no_insert_count = 0
        start, end = range_tuple

        for char in line:
            if char.isupper() or char == "-":
                no_insert_count += 1
            # we keep insertions
            if start < no_insert_count <= end:
                seq.append(char)
            elif no_insert_count > end:
                break

        return "".join(seq)

    def split_sequences(self) -> None:
        """Split A3M file into pairing and non-pairing sequences."""
        out_dir = Path(self.out_dir) / "msa"
        chain_names, seq_ranges = self.chain_info

        pairing_a3ms = {name: [] for name in chain_names}
        nonpairing_a3ms = {name: [] for name in chain_names}

        current_query = None
        for line in self.a3m_content.split("\n"):
            if line.startswith("#"):
                continue

            if line.startswith(">"):
                name = line[1:]
                if name in chain_names:
                    current_query = chain_names[chain_names.index(name)]
                elif name == "\t".join(chain_names):
                    current_query = None

                # Add header line to appropriate dictionary
                if current_query:
                    nonpairing_a3ms[current_query].append(line)
                else:
                    for name in chain_names:
                        pairing_a3ms[name].append(line)
                continue

            # Process sequence line
            if not line:
                continue

            if current_query:
                seq = self._extract_sequence(line, seq_ranges[current_query])
                nonpairing_a3ms[current_query].append(seq)
            else:
                for name in chain_names:
                    seq = self._extract_sequence(line, seq_ranges[name])
                    pairing_a3ms[name].append(seq)

        self._write_output_files(out_dir, nonpairing_a3ms, pairing_a3ms)

    def _write_output_files(
        self,
        out_dir: Path,
        nonpairing_a3ms: Dict[str, List[str]],
        pairing_a3ms: Dict[str, List[str]],
    ) -> None:
        """Write split sequences to output files."""
        out_dir.mkdir(exist_ok=True)

        # Write non-pairing sequences
        for i, (name, lines) in enumerate(nonpairing_a3ms.items()):
            chain_dir = out_dir / str(i)
            chain_dir.mkdir(exist_ok=True)

            with open(chain_dir / "non_pairing.a3m", "w") as f:
                query_seq = lines[1]
                f.write(">query\n")
                f.write(f"{query_seq}\n")
                f.write("\n".join(lines[2:]))

        # Write pairing sequences
        for i, (name, lines) in enumerate(pairing_a3ms.items()):
            chain_dir = out_dir / str(i)
            chain_dir.mkdir(exist_ok=True)

            with open(chain_dir / "pairing.a3m", "w") as f:
                query_seq = lines[1]
                f.write(">query\n")
                f.write(f"{query_seq}\n")

                # Process remaining sequences
                sequences = {}
                for j, line in enumerate(lines[2:]):
                    if line.startswith(">"):
                        current_name = f"UniRef100_{line[1:].split()[i]}_{j}"
                        sequences[current_name] = ""
                    elif line and "DUMMY" not in current_name:
                        sequences[current_name] = line

                # Write processed sequences
                for seq_name, seq in sequences.items():
                    if seq:  # Only write non-empty sequences
                        f.write(f">{seq_name}\n{seq}\n")


def run_colabfold_search(config: LocalColabFoldConfig) -> str:
    """Run ColabFold search with given configuration."""
    cmd = [config.colabsearch, config.query_fpath, config.db_dir, config.results_dir]

    # Add optional parameters
    if config.db1:
        cmd.extend(["--db1", config.db1])
    if config.db2:
        cmd.extend(["--db2", config.db2])
    if config.db3:
        cmd.extend(["--db3", config.db3])
    if config.mmseqs_path:
        cmd.extend(["--mmseqs", config.mmseqs_path])
    else:
        cmd.extend(["--mmseqs", "mmseqs"])
    if config.use_env:
        cmd.extend(["--use-env", str(config.use_env)])
    if config.filter:
        cmd.extend(["--filter", str(config.filter)])
    if config.db_load_mode:
        cmd.extend(["--db-load-mode", str(config.db_load_mode)])

    cmd = " ".join(cmd)
    os.system(cmd)

    # Return the first .a3m file found in results directory
    result_files = list(Path(config.results_dir).glob("*.a3m"))
    if not result_files:
        raise FileNotFoundError(f"No .a3m files found in {config.results_dir}")
    return str(result_files[0])


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="ColabFold search and A3M processing tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments
    parser.add_argument("query_fpath", help="Path to the query FASTA file")
    parser.add_argument("db_dir", help="Directory containing the databases")
    parser.add_argument("results_dir", help="Directory for storing results")

    # Optional arguments
    parser.add_argument(
        "--colabsearch", help="Path to colabfold_search", default="colabfold_search"
    )
    parser.add_argument(
        "--mmseqs_path", help="Path to MMseqs2 binary", default="mmseqs"
    )
    parser.add_argument("--db1", help="First database name", default="uniref30_2302_db")
    parser.add_argument("--db2", help="Templates database")
    parser.add_argument(
        "--db3", help="Environmental database (default: colabfold_envdb_202108_db)"
    )
    parser.add_argument(
        "--use_env", help="Use environment settings", type=int, default=1
    )
    parser.add_argument("--filter", help="Apply filtering", type=int, default=1)
    parser.add_argument(
        "--db_load_mode", help="Database load mode", type=int, default=0
    )
    parser.add_argument(
        "--output_split", help="Directory for split A3M files", default=None
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Create configuration from arguments
    config = LocalColabFoldConfig(
        colabsearch=args.colabsearch,
        query_fpath=args.query_fpath,
        db_dir=args.db_dir,
        results_dir=args.results_dir,
        mmseqs_path=args.mmseqs_path,
        db1=args.db1,
        db2=args.db2,
        db3=args.db3,
        use_env=args.use_env,
        filter=args.filter,
        db_load_mode=args.db_load_mode,
    )

    # Run search
    results_a3m = run_colabfold_search(config)

    processor = A3MProcessor(results_a3m, args.results_dir)
    if len(processor.chain_info) == 2:
        processor.split_sequences()
