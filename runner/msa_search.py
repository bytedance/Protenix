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
from typing import Any, Sequence, Tuple

from protenix.utils.logger import get_logger
from protenix.web_service.colab_request_parser import RequestParser

logger = get_logger(__name__)


def need_msa_search(json_data: dict) -> bool:
    """
    Check if the input JSON data needs an MSA search.

    Args:
        json_data (dict): The input JSON data for a task.

    Returns:
        bool: True if an MSA search is required, False otherwise.
    """
    need_msa = False
    # the new format of msa filed is `pairedMsaPath` and `unpairedMsaPath`
    # we need to check `pairedMsaPath` and `unpairedMsaPath`
    for sequence in json_data["sequences"]:
        if "proteinChain" in sequence:
            protein_chain = sequence["proteinChain"]
            paired_msa_path = protein_chain.get("pairedMsaPath")
            unpaired_msa_path = protein_chain.get("unpairedMsaPath")

            if paired_msa_path is None and unpaired_msa_path is None:
                need_msa = True
            else:
                if paired_msa_path is not None and not os.path.exists(paired_msa_path):
                    logger.warning(
                        f"pairedMsaPath {paired_msa_path} does not exist, will re-search MSA."
                    )
                    need_msa = True
                if unpaired_msa_path is not None and not os.path.exists(
                    unpaired_msa_path
                ):
                    logger.warning(
                        f"unpairedMsaPath {unpaired_msa_path} does not exist, will re-search MSA."
                    )
                    need_msa = True
    return need_msa


def convert_msa_to_new_format(
    data_list: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], bool]:
    """
    Convert MSA format from old format to new format in a list of task dictionaries.

    Args:
        data_list (list[dict[str, Any]]): List of task dictionaries.

    Returns:
        tuple[list[dict[str, Any]], bool]:
            - Updated list of task dictionaries.
            - True if any format conversion occurred, False otherwise.
    """
    result = []
    json_need_converted = False
    for item in data_list:
        # Process each dictionary item
        processed_item, format_converted = convert_one_json_dict(item)
        json_need_converted = json_need_converted or format_converted
        result.append(processed_item)

    return result, json_need_converted


def convert_one_json_dict(obj: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    """
    Process a single dictionary to convert 'msa' field to 'pairedMsaPath' and 'unpairedMsaPath'.

    Args:
        obj (dict[str, Any]): A single task dictionary.

    Returns:
        tuple[dict[str, Any], bool]: Updated dictionary and format conversion flag.
    """
    format_converted = False
    if "sequences" in obj and isinstance(obj["sequences"], list):
        for sequence in obj["sequences"]:
            if "proteinChain" in sequence and isinstance(
                sequence["proteinChain"], dict
            ):
                protein_chain = sequence["proteinChain"]
                if "msa" in protein_chain:
                    msa_info = protein_chain.pop("msa")  # Remove old msa field
                    logger.info(
                        f"Detecting old MSA format: {msa_info}, converting to new format."
                    )
                    precomputed_msa_dir = msa_info.get("precomputed_msa_dir")
                    format_converted = True
                    if precomputed_msa_dir:
                        # Build new paths
                        pairing_path = f"{precomputed_msa_dir}/pairing.a3m"
                        non_pairing_path = f"{precomputed_msa_dir}/non_pairing.a3m"

                        # Add new fields if the files exist
                        if os.path.exists(pairing_path):
                            protein_chain["pairedMsaPath"] = pairing_path
                        if os.path.exists(non_pairing_path):
                            protein_chain["unpairedMsaPath"] = non_pairing_path

    return obj, format_converted


def msa_search(
    seqs: Sequence[str], msa_res_dir: str, mode: str = "protenix"
) -> Sequence[str]:
    """
    Perform MSA search using MMseqs2 and return the resulting subdirectories.

    Args:
        seqs (Sequence[str]): List of protein sequences.
        msa_res_dir (str): Directory to save MSA results.
        mode (str): MSA search mode ('protenix' or 'colabfold').

    Returns:
        Sequence[str]: List of directories containing MSA results for each sequence.
    """
    os.makedirs(msa_res_dir, exist_ok=True)
    tmp_fasta_fpath = os.path.join(msa_res_dir, f"tmp_{uuid.uuid4().hex}.fasta")
    msa_res_subdirs = RequestParser.msa_search(
        seqs_pending_msa=seqs,
        tmp_fasta_fpath=tmp_fasta_fpath,
        msa_res_dir=msa_res_dir,
        mode=mode,
    )
    if mode == "protenix":
        msa_res_subdirs = RequestParser.msa_postprocess(
            seqs_pending_msa=seqs,
            msa_res_dir=msa_res_dir,
        )
    return msa_res_subdirs


def update_seq_msa(infer_seq: dict, msa_res_dir: str, mode: str) -> dict:
    """
    Update the sequences in the inference dictionary with their respective MSA paths.

    Args:
        infer_seq (dict): The task data containing sequences.
        msa_res_dir (str): Directory where MSA results are stored.
        mode (str): MSA search mode.

    Returns:
        dict: The updated task data.
    """
    protein_seqs = []
    for sequence in infer_seq["sequences"]:
        if "proteinChain" in sequence.keys():
            protein_seqs.append(sequence["proteinChain"]["sequence"])
    if len(protein_seqs) > 0:
        protein_seqs = sorted(protein_seqs)
        msa_res_subdirs = msa_search(protein_seqs, msa_res_dir, mode=mode)

        assert len(msa_res_subdirs) == len(protein_seqs), "msa search failed"
        protein_msa_res = dict(zip(protein_seqs, msa_res_subdirs))
        for sequence in infer_seq["sequences"]:
            if "proteinChain" in sequence.keys():
                precomputed_msa_dir = protein_msa_res[
                    sequence["proteinChain"]["sequence"]
                ]
                if os.path.exists(f"{precomputed_msa_dir}/pairing.a3m"):
                    sequence["proteinChain"][
                        "pairedMsaPath"
                    ] = f"{precomputed_msa_dir}/pairing.a3m"
                if os.path.exists(f"{precomputed_msa_dir}/non_pairing.a3m"):
                    sequence["proteinChain"][
                        "unpairedMsaPath"
                    ] = f"{precomputed_msa_dir}/non_pairing.a3m"

    return infer_seq


def update_infer_json(
    json_file: str, out_dir: str, use_msa: bool = True, mode: str = "protenix"
) -> Tuple[str, bool]:
    """
    Update the inference JSON file with MSA information.
    Iterates through tasks and runs MSA search if required and enabled.

    Args:
        json_file (str): Path to the input JSON file.
        out_dir (str): Directory to save MSA results.
        use_msa (bool): Whether to perform MSA search if missing.
        mode (str): MSA search mode ('protenix' or 'colabfold').

    Returns:
        Tuple[str, bool]:
            - Path to the updated (or original) JSON file.
            - Boolean indicating if any MSA search was actually performed.
    """
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"Input file `{json_file}` does not exist.")
    with open(json_file, "r") as f:
        json_data = json.load(f)

    # Change the old format of msa filed to new format
    json_data, json_need_converted = convert_msa_to_new_format(json_data)

    actual_updated = False
    for task_idx, infer_data in enumerate(json_data):
        if use_msa and need_msa_search(infer_data):
            actual_updated = True
            task_name = infer_data.get("name", f"task_{task_idx}")
            logger.info(
                f"starting to update msa result for task {task_idx} in {json_file}"
            )
            update_seq_msa(
                infer_data,
                os.path.join(out_dir, task_name, "msa"),
                mode,
            )
    if actual_updated or json_need_converted:
        updated_json = os.path.join(
            os.path.dirname(os.path.abspath(json_file)),
            f"{os.path.splitext(os.path.basename(json_file))[0]}-update-msa.json",
        )
        with open(updated_json, "w") as f:
            json.dump(json_data, f, indent=4)
        logger.info(f"update msa result success and save to {updated_json}")
        return updated_json, actual_updated
    elif not use_msa:
        logger.warning(
            f"the inference json file {json_file} \n"
            "do not contain msa and will not be updated,\n"
            "and you set not using msa, in this mode, \n"
            "if you do not use esm feature, model performance "
            "might degrade significantly"
        )
        return json_file, actual_updated
    else:
        logger.info(f"do not need to update msa result, so return itself {json_file}")
        return json_file, actual_updated
