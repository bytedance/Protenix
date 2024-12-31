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
import uuid
from typing import Sequence

from protenix.utils.logger import get_logger
from protenix.web_service.colab_request_parser import RequestParser

logger = get_logger(__name__)


def contain_msa_res(json_file: str) -> bool:
    """
    check the json_path data has msa result or not.
    """
    if not os.path.exists(json_file):
        raise RuntimeError(f"`{json_file}` not exists.")
    with open(json_file, "r") as f:
        json_data = json.load(f)
    for seq in json_data:
        for sequence in seq["sequences"]:
            if "proteinChain" in sequence.keys():
                proteinChain = sequence["proteinChain"]
                if "msa" not in proteinChain.keys() or len(proteinChain["msa"]) == 0:
                    return False
    return True


def update_msa_res(seq: dict, protein_msa_res: dict) -> dict:
    for sequence in seq["sequences"]:
        if "proteinChain" in sequence.keys():
            sequence["proteinChain"]["msa"] = {
                "precomputed_msa_dir": protein_msa_res[
                    sequence["proteinChain"]["sequence"]
                ],
                "pairing_db": "uniref100",
            }
    return seq


def msa_search(seqs: Sequence[str], msa_res_dir: str) -> Sequence[str]:
    """
    do msa search with mmseqs and return result subdirs.
    """
    os.makedirs(msa_res_dir, exist_ok=True)
    tmp_fasta_fpath = os.path.join(msa_res_dir, f"tmp_{uuid.uuid4().hex}.fasta")
    RequestParser.msa_search(
        seqs_pending_msa=seqs,
        tmp_fasta_fpath=tmp_fasta_fpath,
        msa_res_dir=msa_res_dir,
    )
    msa_res_subdirs = RequestParser.msa_postprocess(
        seqs_pending_msa=seqs,
        msa_res_dir=msa_res_dir,
    )
    return msa_res_subdirs


def msa_search_update(json_file: str, out_dir: str) -> str:
    """
    do msa search with mmseqs from json input and update it.
    """
    assert os.path.exists(json_file), f"input file {json_file} not exists."
    if contain_msa_res(json_file):
        logger.warning(f"{json_file} has already msa result, skip.")
        return json_file
    with open(json_file, "r") as f:
        input_json_data = json.load(f)
    logger.info(f"starting to update msa result for {json_file}")
    for seq_idx, seq in enumerate(input_json_data):
        protein_seqs = []
        seq_name = seq.get("name", f"seq_{seq_idx}")
        for sequence in seq["sequences"]:
            if "proteinChain" in sequence.keys():
                protein_seqs.append(sequence["proteinChain"]["sequence"])
        if len(protein_seqs) > 0:
            protein_seqs = sorted(protein_seqs)
            msa_res_subdirs = msa_search(
                protein_seqs,
                os.path.join(out_dir, seq_name, "msa_res" f"msa_seq_{seq_idx}"),
            )
            assert len(msa_res_subdirs) == len(msa_res_subdirs), "msa search failed"
            update_msa_res(seq, dict(zip(protein_seqs, msa_res_subdirs)))
    msa_input_json = os.path.join(
        os.path.dirname(json_file),
        f"{os.path.splitext(os.path.basename(json_file))[0]}-add-msa.json",
    )
    with open(msa_input_json, "w") as f:
        json.dump(input_json_data, f, indent=4)
    logger.info(f"update msa result success and save to {msa_input_json}")
    return msa_input_json
