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

import json
import os
from functools import partial
from os.path import join as opjoin
from typing import Callable, Tuple

with open("./scripts/msa/data/pdb_seqs/seq_to_pdb_id_entity_id.json", "r") as f:
    seq_to_pdbid = json.load(f)
first_pdbid_to_seq = {"_".join(v[0]): k for k, v in seq_to_pdbid.items()}

with open("./scripts/msa/data/pdb_seqs/seq_to_pdb_index.json", "r") as f:
    seq_to_pdb_index = json.load(f)


def rematch(pdb_line: str) -> Tuple[str, str]:
    pdb_id = pdb_line[1:-1]
    origin_query_seq = first_pdbid_to_seq[pdb_id]
    pdb_index = seq_to_pdb_index[origin_query_seq]
    return pdb_index, origin_query_seq


def write_log(
    msg: str,
    fname: str,
    log_root: str,
) -> None:
    basename = fname.split(".")[0]
    with open(opjoin(log_root, f"{basename}-{msg}"), "w") as f:
        pass


def process_one_file(
    fname: str, msa_root: str, save_root: str, logger: Callable
) -> None:
    with open(file_path := opjoin(msa_root, fname), "r") as f:
        for i, line in enumerate(f):
            if i == 0:
                pdb_line = line
            if i == 1:
                if len(line) == 1:
                    logger("empty_query_seq", fname)
                    return
                query_line = line
                break

    save_fname, origin_query_seq = rematch(pdb_line)

    os.makedirs(sub_dir_path := opjoin(save_root, f"{save_fname}"), exist_ok=True)
    uniref100_lines = [">query\n", f"{origin_query_seq}\n"]
    other_lines = [">query\n", f"{origin_query_seq}\n"]

    with open(file_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if i < 2:
            continue
        if i % 2 == 0:
            # header
            if not line.startswith(">"):
                logger(f"bad_header_{i}", fname)
                return
            seq = lines[i + 1]

            if line.startswith(">UniRef100"):
                uniref100_lines.extend([line, seq])
            else:
                other_lines.extend([line, seq])

    assert len(other_lines) + len(uniref100_lines) - 2 == len(lines)

    other_lines = other_lines[0:2] + other_lines[4:]
    for i, line in enumerate(other_lines):
        if i > 0 and i % 2 == 0:
            assert "\t" in line
    with open(opjoin(sub_dir_path, "uniref100_hits.a3m"), "w") as f:
        for line in uniref100_lines:
            f.write(line)
    with open(opjoin(sub_dir_path, "mmseqs_other_hits.a3m"), "w") as f:
        for line in other_lines:
            f.write(line)


if __name__ == "__main__":
    msa_root = "./scripts/msa/data/mmcif_msa_with_taxid"
    save_root = "./scripts/msa/data/mmcif_msa"
    log_root = "./scripts/msa/data/mmcif_msa_log"

    os.makedirs(log_root, exist_ok=True)
    os.makedirs(save_root, exist_ok=True)

    print("Loading file names...")

    logger = partial(write_log, log_root=log_root)
    for fname in os.listdir(msa_root):
        process_one_file(
            fname=fname, msa_root=msa_root, save_root=save_root, logger=logger
        )
