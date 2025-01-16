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
from os.path import join as opjoin
from typing import Dict, List

from tqdm import tqdm


def read_a3m(a3m_file: str) -> tuple[List[str], List[str]]:
    """read a3m file from output of mmseqs

    Args:
        a3m_file (str): the a3m file searched by mmseqs(colabfold search)

    Returns:
        tuple[List[str], List[str]]: the header and seqs of a3m files
    """
    heads = []
    seqs = []
    # Record the row index. The index before this index is the MSA of Uniref30 DB,
    # and the index after this index is the MSA of ColabfoldDB.
    uniref_index = 0
    with open(a3m_file, "r") as infile:
        for idx, line in enumerate(infile):
            if line.startswith(">"):
                heads.append(line)
                if idx == 0:
                    query_name = line
                elif idx > 0 and line == query_name:
                    uniref_index = idx
            else:
                seqs.append(line)
    return heads, seqs, uniref_index


def read_m8(m8_file: str) -> Dict[str, str]:
    """the uniref_tax.m8 from output of mmseqs

    Args:
        m8_file (str): the uniref_tax.m8 from output of mmseqs(colabfold search)

    Returns:
        Dict[str, str]: the dict mapping uniref hit_name to NCBI TaxID
    """
    uniref_to_ncbi_taxid = {}
    with open(m8_file, "r") as infile:
        for line in infile:
            line_list = line.replace("\n", "").split("\t")
            hit_name = line_list[1]
            ncbi_taxid = line_list[2]
            uniref_to_ncbi_taxid[hit_name] = ncbi_taxid
    return uniref_to_ncbi_taxid


def update_a3m(
    a3m_path: str,
    uniref_to_ncbi_taxid: Dict,
    save_root: str,
) -> None:
    """add NCBI TaxID to header if "UniRef" in header

    Args:
        a3m_path (str): the original a3m path returned by mmseqs(colabfold search)
        uniref_to_ncbi_taxid (Dict): the dict mapping uniref hit_name to NCBI TaxID
        save_root (str): the updated a3m
    """
    heads, seqs, uniref_index = read_a3m(a3m_path)
    print(uniref_index)
    fname = a3m_path.split("/")[-1]
    out_a3m_path = opjoin(save_root, fname)
    with open(out_a3m_path, "w") as ofile:
        for idx, (head, seq) in enumerate(zip(heads, seqs)):
            uniref_id = head.split("\t")[0][1:]
            ncbi_taxid = uniref_to_ncbi_taxid.get(uniref_id, None)
            if (ncbi_taxid is not None) and (idx < (uniref_index // 2)):
                if not uniref_id.startswith("UniRef100_"):
                    head = head.replace(
                        uniref_id, f"UniRef100_{uniref_id}_{ncbi_taxid}/"
                    )
                else:
                    head = head.replace(uniref_id, f"{uniref_id}_{ncbi_taxid}/")
            ofile.write(f"{head}{seq}")


if __name__ == "__main__":
    input_msa_dir = "./scripts/msa/data/mmcif_msa_initial"

    output_msa_dir = "./scripts/msa/data/mmcif_msa_with_taxid"
    os.makedirs(output_msa_dir, exist_ok=True)

    a3m_paths = os.listdir(input_msa_dir)
    a3m_paths = [opjoin(input_msa_dir, x) for x in a3m_paths if x.endswith(".a3m")]
    m8_file = f"{input_msa_dir}/uniref_tax.m8"
    uniref_to_ncbi_taxid = read_m8(m8_file)
    for a3m_path in tqdm(a3m_paths):
        update_a3m(
            a3m_path=a3m_path,
            uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
            save_root=output_msa_dir,
        )
