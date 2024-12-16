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
import logging
import os
import tempfile
import time
import uuid
from pathlib import Path
from typing import List, Optional, Union

import click
import tqdm
from Bio import SeqIO
from configs.configs_base import configs as configs_base
from configs.configs_data import data_configs
from configs.configs_inference import inference_configs
from rdkit import Chem
from runner.inference import InferenceRunner, download_infercence_cache, infer_predict
from runner.msa_search import contain_msa_res, msa_search, msa_search_update

from protenix.config import parse_configs
from protenix.data.json_maker import cif_to_input_json
from protenix.data.json_parser import lig_file_to_atom_info
from protenix.data.utils import pdb_to_cif
from protenix.utils.logger import get_logger

logger = get_logger(__name__)


def init_logging():
    LOG_FORMAT = "%(asctime)s,%(msecs)-3d %(levelname)-8s [%(filename)s:%(lineno)s %(funcName)s] %(message)s"
    logging.basicConfig(
        format=LOG_FORMAT,
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        filemode="w",
    )


def generate_infer_jsons(
    protein_msa_res: dict, ligand_file: str, seeds: List[int] = [101]
) -> List[str]:
    protein_chains = []
    if len(protein_msa_res) <= 0:
        raise RuntimeError(f"invalid `protein_msa_res` data in {protein_msa_res}")
    for key, value in protein_msa_res.items():
        protein_chain = {}
        protein_chain["proteinChain"] = {}
        protein_chain["proteinChain"]["sequence"] = key
        protein_chain["proteinChain"]["count"] = value.get("count", 1)
        protein_chain["proteinChain"]["msa"] = value
        protein_chains.append(protein_chain)
    if os.path.isdir(ligand_file):
        ligand_files = [
            str(file) for file in Path(ligand_file).rglob("*") if file.is_file()
        ]
        if len(ligand_files) == 0:
            raise RuntimeError(
                f"can not read a valid `sdf` or `smi` ligand_file in {ligand_file}"
            )
    elif os.path.isfile(ligand_file):
        ligand_files = [ligand_file]
    else:
        raise RuntimeError(f"can not read a special ligand_file: {ligand_file}")

    invalid_ligand_files = []
    sdf_ligand_files = []
    smi_ligand_files = []
    tmp_json_name = uuid.uuid4().hex
    current_local_dir = (
        f"/tmp/{time.strftime('%Y-%m-%d', time.localtime())}/{tmp_json_name}"
    )
    current_local_json_dir = (
        f"/tmp/{time.strftime('%Y-%m-%d', time.localtime())}/{tmp_json_name}_jsons"
    )
    os.makedirs(current_local_dir, exist_ok=True)
    os.makedirs(current_local_json_dir, exist_ok=True)
    for li_file in ligand_files:
        try:
            if li_file.endswith(".smi"):
                smi_ligand_files.append(li_file)
            elif li_file.endswith(".sdf"):
                suppl = Chem.SDMolSupplier(li_file)
                if len(suppl) <= 1:
                    lig_file_to_atom_info(li_file)
                    sdf_ligand_files.append([li_file])
                else:
                    sdf_basename = os.path.join(
                        current_local_dir, os.path.basename(li_file).split(".")[0]
                    )
                    li_files = []
                    for idx, mol in enumerate(suppl):
                        p_sdf_path = f"{sdf_basename}_part_{idx}.sdf"
                        writer = Chem.SDWriter(p_sdf_path)
                        writer.write(mol)
                        writer.close()
                        li_files.append(p_sdf_path)
                        lig_file_to_atom_info(p_sdf_path)
                    sdf_ligand_files.append(li_files)
            else:
                lig_file_to_atom_info(li_file)
                sdf_ligand_files.append(li_file)
        except Exception as exc:
            logging.info(f" lig_file_to_atom_info failed with error info: {exc}")
            invalid_ligand_files.append(li_file)
    logger.info(f"the json to infer will be save to {current_local_json_dir}")
    infer_json_files = []
    for li_files in sdf_ligand_files:
        one_infer_seq = protein_chains[:]
        for li_file in li_files:
            ligand_name = os.path.basename(li_file).split(".")[0]
            ligand_chain = {}
            ligand_chain["ligand"] = {}
            ligand_chain["ligand"]["ligand"] = f"FILE_{li_file}"
            ligand_chain["ligand"]["count"] = 1
            one_infer_seq.append(ligand_chain)
        one_infer_json = [{"sequences": one_infer_seq, "name": ligand_name}]
        json_file_name = os.path.join(
            current_local_json_dir, f"{ligand_name}_sdf_{uuid.uuid4().hex}.json"
        )
        with open(json_file_name, "w") as f:
            json.dump(one_infer_json, f, indent=4)
        infer_json_files.append(json_file_name)

    for smi_ligand_file in smi_ligand_files:
        one_infer_seq = protein_chains[:]
        with open(smi_ligand_file, "r") as f:
            smile_list = f.readlines()
        one_infer_seq = protein_chains[:]
        ligand_name = os.path.basename(smi_ligand_file).split(".")[0]
        for smile in smile_list:
            normalize_smile = smile.replace("\n", "")
            ligand_chain = {}
            ligand_chain["ligand"] = {}
            ligand_chain["ligand"]["ligand"] = normalize_smile
            ligand_chain["ligand"]["count"] = 1
            one_infer_seq.append(ligand_chain)
        one_infer_json = [{"sequences": one_infer_seq, "name": ligand_name}]
        json_file_name = os.path.join(
            current_local_json_dir, f"{ligand_name}_smi_{uuid.uuid4().hex}.json"
        )
        with open(json_file_name, "w") as f:
            json.dump(one_infer_json, f, indent=4)
        infer_json_files.append(json_file_name)
    if len(invalid_ligand_files) > 0:
        logger.warning(
            f"{len(invalid_ligand_files)} sdf file is invaild, one of them is {invalid_ligand_files[0]}"
        )
    return infer_json_files


def get_default_runner(seeds: Optional[list] = None) -> InferenceRunner:
    inference_configs["load_checkpoint_path"] = "/af3-dev/release_model/model_v0.2.0.pt"
    configs_base["use_deepspeed_evo_attention"] = (
        os.environ.get("USE_DEEPSPEED_EVO_ATTTENTION", False) == "true"
    )
    configs_base["model"]["N_cycle"] = 10
    configs_base["sample_diffusion"]["N_sample"] = 5
    configs_base["sample_diffusion"]["N_step"] = 200
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}
    configs = parse_configs(
        configs=configs,
        fill_required_with_null=True,
    )
    if seeds is not None:
        configs.seeds = seeds
    download_infercence_cache(configs, model_version="v0.2.0")
    return InferenceRunner(configs)


def inference_jsons(
    json_file: str,
    out_dir: str = "./output",
    use_msa_server: bool = False,
    seeds: list = [101],
) -> None:
    """
    infer_json: json file or directory, will run infer with these jsons

    """
    infer_jsons = []
    if os.path.isdir(json_file):
        infer_jsons = [
            str(file) for file in Path(json_file).rglob("*") if file.is_file()
        ]
        if len(infer_jsons) == 0:
            raise RuntimeError(
                f"can not read a valid `sdf` or `smi` ligand_file in {json_file}"
            )
    elif os.path.isfile(json_file):
        infer_jsons = [json_file]
    else:
        raise RuntimeError(f"can not read a special ligand_file: {json_file}")
    infer_jsons = [file for file in infer_jsons if file.endswith(".json")]
    logger.info(f"will infer with {len(infer_jsons)} jsons")
    if len(infer_jsons) == 0:
        return

    infer_errors = {}
    inference_configs["dump_dir"] = out_dir
    inference_configs["input_json_path"] = infer_jsons[0]
    runner = get_default_runner(seeds)
    configs = runner.configs
    for idx, infer_json in enumerate(tqdm.tqdm(infer_jsons)):
        try:
            if use_msa_server:
                infer_json = msa_search_update(
                    infer_json, os.path.join(out_dir, f"msa_res_{idx}")
                )
            elif not contain_msa_res(infer_json):
                raise RuntimeError(f"can not find msa for {infer_json}")
            configs["input_json_path"] = infer_json
            if not contain_msa_res(infer_json):
                raise RuntimeError(
                    f"`{infer_json}` has no msa result for `proteinChain`, please add first."
                )
            infer_predict(runner, configs)
        except Exception as exc:
            infer_errors[infer_json] = str(exc)
    if len(infer_errors) > 0:
        logger.warning(f"run inference failed: {infer_errors}")


def batch_inference(
    protein_msa_res: dict,
    ligand_file: str,
    out_dir: str = "./output",
    seeds: List[int] = [101],
) -> None:
    """
    ligand_file: ligand file or directory, should be in sdf format or smi with smlies list;
    protein_msa_res: the msa result for `protein`, like:
        {  "MGHHHHHHHHHHSSGH": {
                "precomputed_msa_dir": "/path/to/msa_pairing/result/msa/1",
                "pairing_db": "uniref100"
            },
            "MAEVIRSSAFWRSFPIFEEFDSE": {
                "precomputed_msa_dir": "/path/to/msa_pairing/result/msa/2",
                "pairing_db": "uniref100"
            }
        }
    out_dir: the infer outout dir, default is `./output`
    """

    infer_jsons = generate_infer_jsons(protein_msa_res, ligand_file, seeds)
    logger.info(f"will infer with {len(infer_jsons)} jsons")
    if len(infer_jsons) == 0:
        return

    infer_errors = {}
    inference_configs["dump_dir"] = out_dir
    inference_configs["input_json_path"] = infer_jsons[0]
    runner = get_default_runner(seeds=seeds)
    configs = runner.configs
    for infer_json in tqdm.tqdm(infer_jsons):
        try:
            configs["input_json_path"] = infer_json
            if not contain_msa_res(infer_json):
                raise RuntimeError(
                    f"`{infer_json}` has no msa result for `proteinChain`, please add first."
                )
            infer_predict(runner, configs)
        except Exception as exc:
            infer_errors[infer_json] = str(exc)
    if len(infer_errors) > 0:
        logger.warning(f"run inference failed: {infer_errors}")


@click.group()
def protenix_cli():
    return


@click.command()
@click.option("--input", type=str, help="json files or dir for inference")
@click.option("--out_dir", default="./output", type=str, help="infer result dir")
@click.option(
    "--seeds", type=str, default="101", help="the inference seed, split by comma"
)
@click.option("--use_msa_server", is_flag=True, help="do msa search or not")
def predict(input, out_dir, seeds, use_msa_server):
    """
    predict: Run predictions with protenix.
    :param input, out_dir, use_msa_server
    :return:
    """
    init_logging()
    logger.info(
        f"run infer with input={input}, out_dir={out_dir}, use_msa_server={use_msa_server}"
    )
    seeds = list(map(int, seeds.split(",")))
    inference_jsons(input, out_dir, use_msa_server, seeds=seeds)


@click.command()
@click.option(
    "--input", type=str, help="pdb or cif files to generate jsons for inference"
)
@click.option("--out_dir", type=str, default="./output", help="dir to save json files")
@click.option(
    "--altloc",
    default="first",
    type=str,
    help=" Select the first altloc conformation of each residue in the input file, \
        or specify the altloc letter for selection. For example, 'first', 'A', 'B', etc.",
)
@click.option(
    "--assembly_id",
    default=None,
    type=str,
    help="Extends the structure based on the Assembly ID in \
                        the input file. The default is no extension",
)
def tojson(input, out_dir="./output", altloc="first", assembly_id=None):
    """
    tojson: convert pdb/cif files or dir to json files for predict.
    :param input, out_dir, altloc, assembly_id
    :return:
    """
    init_logging()
    logger.info(
        f"run tojson with input={input}, out_dir={out_dir}, altloc={altloc}, assembly_id={assembly_id}"
    )
    input_files = []
    if not os.path.exists(input):
        raise RuntimeError(f"input file {input} not exists.")
    if os.path.isdir(input):
        input_files.extend(
            [str(file) for file in Path(input).rglob("*") if file.is_file()]
        )
    elif os.path.isfile(input):
        input_files.append(input)
    else:
        raise RuntimeError(f"can not read a special file: {input}")

    input_files = [
        file for file in input_files if file.endswith(".pdb") or file.endswith(".cif")
    ]
    if len(input_files) == 0:
        raise RuntimeError(f"can not read a valid `pdb` or `cif` file from {input}")
    logger.info(
        f"will tojson jsons for {len(input_files)} input files with `pdb` or `cif` format."
    )
    output_jsons = []
    os.makedirs(out_dir, exist_ok=True)
    for input_file in input_files:
        stem, _ = os.path.splitext(os.path.basename(input_file))
        pdb_name = stem[:20]
        output_json = os.path.join(out_dir, f"{pdb_name}-{uuid.uuid4().hex}.json")
        if input_file.endswith(".pdb"):
            with tempfile.NamedTemporaryFile(suffix=".cif") as tmp:
                tmp_cif_file = tmp.name
                pdb_to_cif(input_file, tmp_cif_file)
                cif_to_input_json(
                    tmp_cif_file,
                    assembly_id=assembly_id,
                    altloc=altloc,
                    sample_name=pdb_name,
                    output_json=output_json,
                )
        elif input_file.endswith(".cif"):
            cif_to_input_json(
                input_file,
                assembly_id=assembly_id,
                altloc=altloc,
                output_json=output_json,
            )
        else:
            raise RuntimeError(f"can not read a special ligand_file: {input_file}")
        output_jsons.append(output_json)
    logger.info(f"{len(output_jsons)} generated jsons have been save to {out_dir}.")
    return output_jsons


@click.command()
@click.option(
    "--input", type=str, help="file to do msa search, support `json` or `fasta` format"
)
@click.option("--out_dir", type=str, default="./output", help="dir to save msa results")
def msa(input, out_dir) -> Union[str, dict]:
    """
    msa: do msa search by mmseqs. If input is in `fasta`, it should all be proteinChain.
    :param input, out_dir
    :return:
    """
    init_logging()
    out_dir = os.path.join(out_dir, uuid.uuid4().hex)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f"run msa with input={input}, out_dir={out_dir}")
    if input.endswith(".json"):
        msa_input_json = msa_search_update(input, out_dir)
        logger.info(f"msa results have been update to {msa_input_json}")
        return msa_input_json
    elif input.endswith(".fasta"):
        records = list(SeqIO.parse(input, "fasta"))
        protein_seqs = []
        for seq in records:
            protein_seqs.append(str(seq.seq))
        protein_seqs = sorted(protein_seqs)
        msa_res_subdirs = msa_search(protein_seqs, out_dir)
        assert len(msa_res_subdirs) == len(msa_res_subdirs), "msa search failed"
        fasta_msa_res = dict(zip(protein_seqs, msa_res_subdirs))
        logger.info(
            f"msa result is: {fasta_msa_res}, and it has been save to {out_dir}"
        )
        return fasta_msa_res
    else:
        raise RuntimeError(f"only support `json` or `fasta` format, but got : {input}")


protenix_cli.add_command(predict)
protenix_cli.add_command(tojson)
protenix_cli.add_command(msa)


def test_batch_inference():
    ligands_dir = "../examples/ligands"
    protein_msa_res = {
        "MASWSHPQFEKGGTHVAETSAPTRSEPDTRVLTLPGTASAPEFRLIDIDGLLNNRATTDVRDLGSGRLNAWGNSFPAAELPAPGSLITVAGIPFTWANAHARGDNIRCEGQVVDIPPGQYDWIYLLAASERRSEDTIWAHYDDGHADPLRVGISDFLDGTPAFGELSAFRTSRMHYPHHVQEGLPTTMWLTRVGMPRHGVARSLRLPRSVAMHVFALTLRTAAAVRLAEGATT": {
            "precomputed_msa_dir": "../examples/7wux/msa/1",
            "pairing_db": "uniref100",
        },
        "MGSSHHHHHHSQDPNSTTTAPPVELWTRDLGSCLHGTLATALIRDGHDPVTVLGAPWEFRRRPGAWSSEEYFFFAEPDSLAGRLALYHPFESTWHRSDGDGVDDLREALAAGVLPIAAVDNFHLPFRPAFHDVHAAHLLVVYRITETEVYVSDAQPPAFQGAIPLADFLASWGSLNPPDDADVFFSASPSGRRWLRTRMTGPVPEPDRHWVGRVIRENVARYRQEPPADTQTGLPGLRRYLDELCALTPGTNAASEALSELYVISWNIQAQSGLHAEFLRAHSVKWRIPELAEAAAGVDAVAHGWTGVRMTGAHSRVWQRHRPAELRGHATALVRRLEAALDLLELAADAVS": {
            "precomputed_msa_dir": "../examples/7wux/msa/2",
            "pairing_db": "uniref100",
        },
    }
    out_dir = "./infer_output"
    batch_inference(protein_msa_res, ligands_dir, out_dir=out_dir)


if __name__ == "__main__":
    init_logging()
    test_batch_inference()
