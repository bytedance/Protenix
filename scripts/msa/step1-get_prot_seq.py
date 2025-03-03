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
from datetime import datetime
from pathlib import Path

import biotite.structure as struc
import pandas as pd
from joblib import Parallel, delayed

from protenix.data.ccd import get_one_letter_code
from protenix.data.parser import MMCIFParser

pd.options.mode.copy_on_write = True


def get_seqs(mmcif_file):
    mmcif_parser = MMCIFParser(mmcif_file)

    entity_poly = mmcif_parser.get_category_table("entity_poly")
    if entity_poly is None:
        pdb_id = mmcif_file.name.split(".")[0]
        return pdb_id, None
    entity_poly["mmcif_seq_old"] = entity_poly.pdbx_seq_one_letter_code_can.str.replace(
        "\n", ""
    )
    entity_poly["pdbx_type"] = entity_poly.type
    mol_type = []
    # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity_poly.type.html
    for i in entity_poly.type:
        if "ribonucleotide" in i:
            mol_type.append("na")
        elif "polypeptide" in i:
            mol_type.append("protein")
        else:
            mol_type.append("other")
    entity_poly["mol_type"] = mol_type

    entity = mmcif_parser.get_category_table("entity")
    info_df = pd.merge(
        entity, entity_poly, left_on="id", right_on="entity_id", how="inner"
    )
    pdb_id = mmcif_file.name.split(".")[0]
    info_df["pdb_id"] = pdb_id

    if "pdbx_audit_revision_history" in mmcif_parser.cif.block:
        history = mmcif_parser.cif.block["pdbx_audit_revision_history"]
        info_df["release_date"] = history["revision_date"].as_array()[0]
    else:
        # Handle non-official mmcif file which transform from pdb file
        info_df["release_date"] = datetime.now().strftime("%Y-%m-%d")

    if mmcif_parser.release_date:
        info_df["release_date_retrace_obsolete"] = mmcif_parser.release_date
    else:   
        # Handle non-official mmcif file which transform from pdb file
        info_df["release_date_retrace_obsolete"] = datetime.now().strftime("%Y-%m-%d")

    entity_poly_seq = mmcif_parser.get_category_table("entity_poly_seq")

    seq_from_resname = []
    diff_seq_mmcif_vs_atom_site = []
    has_alt_res = []
    diff_alt_res_seq_vs_atom_site = []
    for entity_id, mmcif_seq_old in zip(info_df.entity_id, info_df.mmcif_seq_old):
        chain_mask = entity_poly_seq.entity_id == entity_id
        res_names = entity_poly_seq.mon_id[chain_mask].to_numpy(dtype=str)
        res_ids = entity_poly_seq.num[chain_mask].to_numpy(dtype=int)

        seq = ""
        pre_res_id = 0
        id_2_name = {}
        for res_id, res_name in zip(res_ids, res_names):
            if res_id == pre_res_id:
                continue
            id_2_name[res_id] = res_name
            one = get_one_letter_code(res_name)
            if one is None:
                one = "X"
            if len(one) > 1:
                one = "X"
            seq += one
            pre_res_id = res_id
        assert len(seq) == max(res_ids)

        diff_seq_mmcif_vs_atom_site.append(seq != mmcif_seq_old)
        has_alt = False
        mismatch_res_name = False
        if len(seq) < len(res_ids):  # has altloc residue in same res_id
            has_alt = True
            # get_structure() return atom array only keep first altloc residue
            atom_array = mmcif_parser.get_structure()
            res_starts = struc.get_residue_starts(atom_array)
            for start in res_starts:
                if atom_array.label_entity_id[start] == entity_id:
                    first_res_in_seq = id_2_name[atom_array.res_id[start]]
                    first_res_in_atom = atom_array.res_name[start]
                    if first_res_in_seq != first_res_in_atom:
                        mismatch_res_name = True
                        break

        has_alt_res.append(has_alt)
        diff_alt_res_seq_vs_atom_site.append(mismatch_res_name)

        seq_from_resname.append(seq)
    info_df["seq"] = seq_from_resname
    info_df["length"] = [len(s) for s in info_df.seq]
    info_df["diff_seq_mmcif_vs_atom_site"] = diff_seq_mmcif_vs_atom_site
    info_df["has_alt_res"] = has_alt_res
    info_df["diff_alt_res_seq_vs_atom_site"] = diff_alt_res_seq_vs_atom_site
    info_df["auth_asym_id"] = info_df["pdbx_strand_id"]

    columns = [
        "pdb_id",
        "entity_id",
        "mol_type",
        "pdbx_type",
        "length",
        "mmcif_seq_old",
        "seq",
        "diff_seq_mmcif_vs_atom_site",
        "has_alt_res",
        "diff_alt_res_seq_vs_atom_site",
        "pdbx_description",
        "auth_asym_id",
        "release_date",
        "release_date_retrace_obsolete",
    ]
    info_df = info_df[columns]
    return pdb_id, info_df


def try_get_seqs(cif_file):
    pdb_id = cif_file.name.split(".")[0]
    try:
        return get_seqs(cif_file)
    except Exception as e:
        print("skip", pdb_id, e)
        return pdb_id, "Error:" + str(e)


def export_to_fasta(df, filename):
    df_protein = df[df["mol_type"] == "protein"]
    # drop duplicates sequence for avoiding duplicate msa search
    df_protein = df_protein.drop_duplicates(subset=["seq"])
    with open(filename, "w") as fasta_file:
        for _, row in df_protein.iterrows():
            header = f">{row['pdb_id']}_{row['entity_id']}\n"
            sequence = f"{row['seq']}\n"
            fasta_file.write(header)
            fasta_file.write(sequence)


def mapping_seqs_to_pdb_entity_id(df, output_json_file):
    df_protein = df[df["mol_type"] == "protein"]
    sequence_mapping = {}

    for _, row in df_protein.iterrows():
        seq = row["seq"]
        key = row["pdb_id"]
        value = row["entity_id"]

        if seq not in sequence_mapping:
            sequence_mapping[seq] = []
        sequence_mapping[seq].append([key, value])

    with open(output_json_file, "w") as json_file:
        json.dump(sequence_mapping, json_file, indent=4)
    return sequence_mapping


def mapping_seqs_to_integer_identifiers(
    sequence_mapping,
    pdb_index_to_seq_path,
    seq_to_pdb_index_path,
):
    seq_to_pdb_index = {}
    for idx, seq in enumerate(sorted(sequence_mapping.keys())):
        seq_to_pdb_index[seq] = idx
    pdb_index_to_seq = {v: k for k, v in seq_to_pdb_index.items()}
    with open(pdb_index_to_seq_path, "w") as f:
        json.dump(pdb_index_to_seq, f, indent=4)
    with open(seq_to_pdb_index_path, "w") as f:
        json.dump(seq_to_pdb_index, f, indent=4)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif_dir", type=str, default="./scripts/msa/data/mmcif")
    parser.add_argument("--out_dir", type=str, default="./scripts/msa/data/pdb_seqs")
    args = parser.parse_args()
    
    cif_dir = Path(args.cif_dir)
    cif_files = [x for x in cif_dir.iterdir() if x.is_file()]

    info_dfs = []
    none_list = []
    error_list = []
    with Parallel(n_jobs=-2, verbose=10) as parallel:
        for pdb_id, info_df in parallel([delayed(try_get_seqs)(f) for f in cif_files]):
            if info_df is None:
                none_list.append(pdb_id)
            elif isinstance(info_df, str) and info_df.startswith("Error:"):
                error_list.append((pdb_id, info_df))
            else:
                info_dfs.append(info_df)

    out_df = pd.concat(info_dfs)
    out_df = out_df.sort_values(["pdb_id", "entity_id"])

    out_dir = Path(args.out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    # 1. extract pdb sequence info
    seq_file = out_dir / "pdb_seq.csv"
    seq_df = out_df[out_df.mol_type != "other"]
    seq_df.to_csv(seq_file, index=False)
    # 2. generate protein fasta file as MSA input
    fasta_file = out_dir / "pdb_seq.fasta"
    export_to_fasta(seq_df, fasta_file)

    # 3. get seq_to_pdb_id_entity_id mapping
    seq_to_pdb_id_entity_id_json = out_dir / "seq_to_pdb_id_entity_id.json"
    sequence_mapping = mapping_seqs_to_pdb_entity_id(
        seq_df, seq_to_pdb_id_entity_id_json
    )

    # 4. mapping sequence with integers identifiers for saving MSA.
    # When we actually store MSA, we need to use simpler integers as
    # identifiers, It's much better than directly use the sequence as identifiers,
    # if there exists long sequences.
    pdb_index_to_seq_path = out_dir / "pdb_index_to_seq.json"
    seq_to_pdb_index_path = out_dir / "seq_to_pdb_index.json"
    mapping_seqs_to_integer_identifiers(
        sequence_mapping, pdb_index_to_seq_path, seq_to_pdb_index_path
    )
