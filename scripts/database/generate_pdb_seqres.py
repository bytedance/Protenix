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
import subprocess
from pathlib import Path

# https://www.rcsb.org/docs/programmatic-access/file-download-services#sequence-data
PDB_SEQRES_URL = "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt"


def filter_pdb_seqres(input_path, output_path):
    """
    Filters a PDB seqres file to keep only entries where 'mol:protein' is in the description.
    """
    if not os.path.exists(input_path):
        print(f"Error: Input file {input_path} does not exist.")
        print(f"You can download it manually from: {PDB_SEQRES_URL}")
        return

    print(f"Filtering {input_path} to {output_path}...")
    with open(input_path, "r") as f_in, open(output_path, "w") as f_out:
        header = None
        seq_lines = []
        count = 0
        total = 0
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    total += 1
                    if "mol:protein" in header:
                        f_out.write(header + "\n")
                        f_out.write("".join(seq_lines) + "\n")
                        count += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # Handle the last entry
        if header:
            total += 1
            if "mol:protein" in header:
                f_out.write(header + "\n")
                f_out.write("".join(seq_lines) + "\n")
                count += 1

    print(f"Processed {total} entries. Kept {count} protein entries.")


def download_pdb_seqres(target_path):
    """
    Downloads pdb_seqres.txt using wget.
    """
    print(f"Downloading {PDB_SEQRES_URL} to {target_path}...")
    try:
        subprocess.run(["wget", PDB_SEQRES_URL, "-O", target_path], check=True)
        print("Download completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading file: {e}")
        return False
    except FileNotFoundError:
        print(
            "Error: 'wget' command not found. Please install wget or download the file manually."
        )
        return False
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter pdb_seqres.txt to keep only proteins. "
        f"The input file can be downloaded from: {PDB_SEQRES_URL}"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Path to input pdb_seqres.txt. If not provided, will look for 'pdb_seqres.txt' locally or download it.",
    )
    parser.add_argument("-o", "--output", type=str, help="Path to output filtered file")

    args = parser.parse_args()

    output_file = args.output
    if output_file is None:
        PROTENIX_ROOT_DIR = os.environ.get("PROTENIX_ROOT_DIR", str(Path.home()))
        output_dir = os.path.join(PROTENIX_ROOT_DIR, "search_database")
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "pdb_seqres_protein.fasta")

    input_file = args.input
    if input_file is None:
        input_file = "pdb_seqres.txt"
        if not os.path.exists(input_file):
            print(f"Input file not specified and '{input_file}' not found locally.")
            if not download_pdb_seqres(input_file):
                exit(1)
    elif not os.path.exists(input_file):
        print(f"Specified input file '{input_file}' not found.")
        if (
            input_file == "pdb_seqres.txt"
            or os.path.basename(input_file) == "pdb_seqres.txt"
        ):
            print("Attempting to download it...")
            if not download_pdb_seqres(input_file):
                exit(1)
        else:
            exit(1)

    filter_pdb_seqres(input_file, output_file)
