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

import logging
import os
import tarfile
import time
from typing import List, Tuple, Dict

import requests
from requests.auth import HTTPBasicAuth
from tqdm import tqdm

TQDM_BAR_FORMAT = "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} estimate remaining: {remaining}]"
logger = logging.getLogger(__name__)

username = "example_user"
password = "example_password"


def parse_fasta_string(fasta_string: str) -> Dict:
    fasta_dict = {}
    lines = fasta_string.strip().split("\n")
    for line in lines:
        if line.startswith(">"):
            header = line[1:].strip()
            fasta_dict[header] = ""
        else:
            fasta_dict[header] += line.strip()
    return fasta_dict


def run_mmseqs2_service(
    x,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    pairing_strategy="complete",
    host_url="https://api.colabfold.com",
    user_agent: str = "",
    email: str = "",
    server_mode: str = "protenix",
) -> Tuple[List[str], List[str]]:
    if server_mode == "protenix":
        assert host_url == "https://protenix-server.com/api/msa"
    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"
    headers = {}
    if user_agent != "":
        headers["User-Agent"] = user_agent
    else:
        logger.warning(
            "No user agent specified. Please set a user agent (e.g., 'toolname/version contact@email') to help us debug in case of problems. This warning will become an error in the future."
        )

    def submit(seqs, mode, N=101):
        n, query = N, ""
        for seq in seqs:
            query += f"{seq}\n"
            n += 1

        while True:
            error_count = 0
            try:
                # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
                # "good practice to set connect timeouts to slightly larger than a multiple of 3"
                res = requests.post(
                    f"{host_url}/{submission_endpoint}",
                    data={"q": query, "mode": mode, "email": email},
                    timeout=6.02,
                    headers=headers,
                    auth=HTTPBasicAuth(username, password),
                )
            except requests.exceptions.Timeout:
                logger.warning("Timeout while submitting to MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break

        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        while True:
            error_count = 0
            try:
                res = requests.get(
                    f"{host_url}/ticket/{ID}",
                    timeout=6.02,
                    headers=headers,
                    auth=HTTPBasicAuth(username, password),
                )
            except requests.exceptions.Timeout:
                logger.warning(
                    "Timeout while fetching status from MSA server. Retrying..."
                )
                continue
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        error_count = 0
        while True:
            try:
                res = requests.get(
                    f"{host_url}/result/download/{ID}",
                    timeout=6.02,
                    headers=headers,
                    auth=HTTPBasicAuth(username, password),
                )
            except requests.exceptions.Timeout:
                logger.warning(
                    "Timeout while fetching result from MSA server. Retrying..."
                )
                continue
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        with open(path, "wb") as out:
            out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # compatibility to old option
    if filter is not None:
        use_filter = filter

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_pairing:
        use_templates = False
        use_env = False
        mode = ""
        # greedy is default, complete was the previous behavior
        if pairing_strategy == "greedy":
            mode = "pairgreedy"
        elif pairing_strategy == "complete":
            mode = "paircomplete"

    # define path
    path = prefix
    if not os.path.isdir(path):
        os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f"{path}/out.tar.gz"
    N, REDO = 101, True

    # deduplicate and keep track of order
    seqs_unique = []
    # TODO this might be slow for large sets
    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # lets do it!
    logger.info("Msa server is running.")
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 100
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 60
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    raise Exception(
                        f"MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later."
                    )

                if out["status"] == "MAINTENANCE":
                    raise Exception(
                        f"MMseqs2 API is undergoing maintenance. Please try again in a few minutes."
                    )

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 60
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                    pbar.n = min(99, int(100 * TIME / (30.0 * 60)))
                    pbar.refresh()
                if out["status"] == "COMPLETE":
                    pbar.n = 100
                    pbar.refresh()
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    raise Exception(
                        f"MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later."
                    )

            # Download results
            download(ID, tar_gz_file)
            with tarfile.open(tar_gz_file) as tar_gz:
                tar_gz.extractall(os.path.dirname(tar_gz_file))
            files = os.listdir(os.path.dirname(tar_gz_file))

            if server_mode == "protenix":
                if (
                    "0.a3m" not in files
                    or "pdb70_220313_db.m8" not in files
                    or "uniref_tax.m8" not in files
                ):
                    raise FileNotFoundError(
                        f"Files 0.a3m, pdb70_220313_db.m8, and uniref_tax.m8 not found in the directory."
                    )
                else:
                    print("Files downloaded and extracted successfully.")
            elif server_mode == "colabfold":
                if not use_pairing:
                    env_a3m_fpath = os.path.join(
                        prefix, "bfd.mgnify30.metaeuk30.smag30.a3m"
                    )
                    env_a3m_dict = parse_fasta_string(
                        open(env_a3m_fpath, "r").read().replace("\x00", "")
                    )
                    uniref_a3m_fpath = os.path.join(prefix, "uniref.a3m")
                    uniref_a3m_dict = parse_fasta_string(
                        open(uniref_a3m_fpath, "r").read().replace("\x00", "")
                    )
                    query_id = str(int(x.split("\n")[0].split("_")[-1]))
                    query_seq = x.split("\n")[1]
                    real_non_pairing_fpath = os.path.join(
                        prefix.split("msa_resmsa")[0],
                        "msa",
                        query_id,
                        "non_pairing.a3m",
                    )
                    output_dir = os.path.dirname(real_non_pairing_fpath)
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    with open(real_non_pairing_fpath, "w") as f:
                        f.write(f">query\n{query_seq}\n")
                        for k, v in env_a3m_dict.items():
                            if k.startswith("query_"):
                                continue
                            else:
                                f.write(f">{k}\n{v}\n")
                        for k, v in uniref_a3m_dict.items():
                            if k.startswith("query_"):
                                continue
                            else:
                                f.write(f">{k}\n{v}\n")
                    return os.path.abspath(os.path.dirname(real_non_pairing_fpath))
                else:
                    # pairing mode
                    pair_a3m = os.path.join(prefix, "pair.a3m")
                    pair_a3m_chunks = open(pair_a3m, "r").read().split("\x00")
                    for chunk in pair_a3m_chunks[:-1]:
                        real_pairing_fpath = os.path.join(
                            prefix.split("msa_resmsa")[0],
                            "msa",
                            str(int(chunk.split("\n")[0].split("_")[-1])),
                            "pairing.a3m",
                        )
                        output_dir = os.path.dirname(real_pairing_fpath)
                        if not os.path.exists(output_dir):
                            os.makedirs(output_dir)
                        chunk_fasta = parse_fasta_string(chunk)
                        with open(real_pairing_fpath, "w") as f:
                            for i, (k, v) in enumerate(chunk_fasta.items()):
                                if k.startswith("query_"):
                                    f.write(f">query\n{v}\n")
                                else:
                                    ks = k.split("\t")
                                    ks[0] = f"{ks[0]}_{i}/"
                                    k = "\t".join(ks)
                                    f.write(f">{k}_{i}\n{v}\n")
