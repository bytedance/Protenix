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


def need_msa_search(json_data: dict) -> bool:
    need_msa = json_data.get("use_msa", True)
    # TODO: add esm check
    if not need_msa:
        return need_msa
    need_msa = False
    for sequence in json_data["sequences"]:
        if "proteinChain" in sequence.keys():
            proteinChain = sequence["proteinChain"]
            if "msa" not in proteinChain.keys() or len(proteinChain["msa"]) == 0:
                need_msa = True
    return need_msa


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


def update_seq_msa(infer_seq: dict, msa_res_dir: str) -> dict:
    protein_seqs = []
    for sequence in infer_seq["sequences"]:
        if "proteinChain" in sequence.keys():
            protein_seqs.append(sequence["proteinChain"]["sequence"])
    if len(protein_seqs) > 0:
        protein_seqs = sorted(protein_seqs)
        msa_res_subdirs = msa_search(protein_seqs, msa_res_dir)
        assert len(msa_res_subdirs) == len(msa_res_subdirs), "msa search failed"
        protein_msa_res = dict(zip(protein_seqs, msa_res_subdirs))
        for sequence in infer_seq["sequences"]:
            if "proteinChain" in sequence.keys():
                sequence["proteinChain"]["msa"] = {
                    "precomputed_msa_dir": protein_msa_res[
                        sequence["proteinChain"]["sequence"]
                    ],
                    "pairing_db": "uniref100",
                }
    return infer_seq


def update_infer_json(
    json_file: str, out_dir: str, use_msa_server: bool = False
) -> str:
    """
    update json file for inference.
    for every infer_data, if it needs to update msa result info,
    it will run msa searching if use_msa_server is True,
    else it will raise error.
    if it does not need to update msa result info, then pass.
    """
    if not os.path.exists(json_file):
        raise RuntimeError(f"`{json_file}` not exists.")
    with open(json_file, "r") as f:
        json_data = json.load(f)

    actual_updated = False
    for seq_idx, infer_data in enumerate(json_data):
        if need_msa_search(infer_data):
            actual_updated = True
            if use_msa_server:
                seq_name = infer_data.get("name", f"seq_{seq_idx}")
                logger.info(
                    f"starting to update msa result for seq {seq_idx} in {json_file}"
                )
                update_seq_msa(
                    infer_data,
                    os.path.join(out_dir, seq_name, "msa_res" f"msa_seq_{seq_idx}"),
                )
            else:
                raise RuntimeError(
                    f"infer seq {seq_idx} in `{json_file}` has no msa result, please add first."
                )
    if actual_updated:
        updated_json = os.path.join(
            os.path.dirname(os.path.abspath(json_file)),
            f"{os.path.splitext(os.path.basename(json_file))[0]}-add-msa.json",
        )
        with open(updated_json, "w") as f:
            json.dump(json_data, f, indent=4)
        logger.info(f"update msa result success and save to {updated_json}")
        return updated_json
    else:
        logger.info(f"do not need to update msa result, so return itself {json_file}")
        return json_file