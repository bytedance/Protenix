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

import difflib
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
from configs.configs_model_type import model_configs
from ml_collections.config_dict import ConfigDict
from protenix.config.config import parse_configs
from protenix.data.inference.json_maker import cif_to_input_json
from protenix.data.inference.json_parser import lig_file_to_atom_info
from protenix.data.utils import pdb_to_cif
from protenix.utils.logger import get_logger
from rdkit import Chem

from runner.inference import (
    download_inference_cache,
    infer_predict,
    InferenceRunner,
    update_gpu_compatible_configs,
)
from runner.msa_search import msa_search, update_infer_json
from runner.rna_msa_search import update_rna_msa_info
from runner.template_search import update_template_info

logger = get_logger(__name__)


def init_logging() -> None:
    """Initialize logging configuration."""
    log_format = (
        "%(asctime)s,%(msecs)-3d %(levelname)-8s "
        "[%(filename)s:%(lineno)s %(funcName)s] %(message)s"
    )
    logging.basicConfig(
        format=log_format,
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        filemode="w",
    )


def preprocess_input(
    input_json: str,
    out_dir: str,
    use_msa: bool = True,
    use_template: bool = False,
    use_rna_msa: bool = False,
    msa_server_mode: str = "protenix",
    hmmsearch_binary_path: Optional[str] = None,
    hmmbuild_binary_path: Optional[str] = None,
    seqres_database_path: Optional[str] = None,
    nhmmer_binary_path: Optional[str] = None,
    hmmalign_binary_path: Optional[str] = None,
    hmmbuild_rna_binary_path: Optional[str] = None,
    ntrna_database_path: Optional[str] = None,
    rfam_database_path: Optional[str] = None,
    rna_central_database_path: Optional[str] = None,
    nhmmer_n_cpu: Optional[int] = None,
) -> str:
    """
    Preprocess the input JSON file by performing MSA, template, and RNA MSA searches as needed.

    Args:
        input_json (str): Path to the input JSON file.
        out_dir (str): Directory to save search results.
        use_msa (bool): Whether to use protein MSA.
        use_template (bool): Whether to use templates.
        use_rna_msa (bool): Whether to use RNA MSA.
        msa_server_mode (str): MSA search mode ('protenix' or 'colabfold').
        hmmsearch_binary_path (Optional[str]): Path to hmmsearch binary.
        hmmbuild_binary_path (Optional[str]): Path to hmmbuild binary.
        seqres_database_path (Optional[str]): Path to sequence database.
        nhmmer_binary_path (Optional[str]): Path to nhmmer binary.
        hmmalign_binary_path (Optional[str]): Path to hmmalign binary.
        hmmbuild_rna_binary_path (Optional[str]): Path to RNA hmmbuild binary.
        ntrna_database_path (Optional[str]): NT-RNA database path.
        rfam_database_path (Optional[str]): Rfam database path.
        rna_central_database_path (Optional[str]): RNAcentral database path.
        nhmmer_n_cpu (Optional[int]): Number of CPUs for nhmmer.

    Returns:
        str: Path to the updated JSON file.
    """
    # 1. Protein MSA search
    msa_updated_json, _ = update_infer_json(
        input_json, out_dir, use_msa=use_msa, mode=msa_server_mode
    )

    # Read the data (either original or updated)
    with open(msa_updated_json, "r") as f:
        json_data = json.load(f)

    actual_updated = False

    # 2. Template search
    if use_template:
        template_updated = update_template_info(
            json_data,
            hmmsearch_binary_path=hmmsearch_binary_path,
            hmmbuild_binary_path=hmmbuild_binary_path,
            seqres_database_path=seqres_database_path,
        )
        actual_updated = actual_updated or template_updated

    # 3. RNA MSA search
    if use_rna_msa:
        rna_updated = update_rna_msa_info(
            json_data,
            out_dir=out_dir,
            nhmmer_binary_path=nhmmer_binary_path,
            hmmalign_binary_path=hmmalign_binary_path,
            hmmbuild_binary_path=hmmbuild_rna_binary_path or hmmbuild_binary_path,
            ntrna_database_path=ntrna_database_path,
            rfam_database_path=rfam_database_path,
            rna_central_database_path=rna_central_database_path,
            nhmmer_n_cpu=nhmmer_n_cpu,
        )
        actual_updated = actual_updated or rna_updated

    if actual_updated:
        base, ext = os.path.splitext(os.path.basename(msa_updated_json))
        if "-update-msa" in base:
            output_json_name = base.replace("-update-msa", "-final-updated") + ext
        else:
            output_json_name = f"{base}-final-updated{ext}"

        output_json = os.path.join(
            os.path.dirname(os.path.abspath(msa_updated_json)), output_json_name
        )

        with open(output_json, "w") as f:
            json.dump(json_data, f, indent=4)
        logger.info(f"Input preprocessing completed, results saved to {output_json}")
        return output_json
    else:
        return msa_updated_json


def generate_infer_jsons(protein_msa_res: dict, ligand_file: str) -> List[str]:
    """
    Generate inference JSON files from protein MSA results and ligand files.

    Args:
        protein_msa_res (dict): Dictionary mapping protein sequences to their MSA results.
        ligand_file (str): Path to a ligand file (SDF or SMI) or directory containing ligand files.

    Returns:
        List[str]: List of paths to the generated inference JSON files.
    """
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
            f"Found {len(invalid_ligand_files)} invalid ligand files. "
            f"Example: {invalid_ligand_files[0]}"
        )
    return infer_json_files


def get_default_runner(
    seeds: Optional[list] = None,
    n_cycle: int = 10,
    n_step: int = 200,
    n_sample: int = 5,
    dtype: str = "bf16",
    model_name: str = "protenix_base_default_v1.0.0",
    use_msa: bool = True,
    trimul_kernel="cuequivariance",
    triatt_kernel="cuequivariance",
    enable_cache=True,
    enable_fusion=True,
    enable_tf32=True,
    use_template: bool = False,
    use_rna_msa: bool = False,
    use_seeds_in_json: bool = False,
) -> InferenceRunner:
    """
    Get a default InferenceRunner with the specified configurations.

    Args:
        seeds (Optional[list]): List of inference seeds.
        n_cycle (int): Number of Pairformer cycles.
        n_step (int): Number of diffusion steps.
        n_sample (int): Number of samples.
        dtype (str): Inference data type (e.g., 'bf16').
        model_name (str): Name of the model checkpoint.
        use_msa (bool): Whether to use MSA.
        trimul_kernel (str): Kernel for triangle multiplicative update.
        triatt_kernel (str): Kernel for triangle attention.
        enable_cache (bool): Whether to enable diffusion shared variables cache.
        enable_fusion (bool): Whether to enable diffusion transformer fusion.
        enable_tf32 (bool): Whether to enable TF32.
        use_template (bool): Whether to use templates.
        use_rna_msa (bool): Whether to use RNA MSA.
        use_seeds_in_json (bool): Whether to use seeds defined in the JSON file.

    Returns:
        InferenceRunner: An instance of InferenceRunner.
    """
    inference_configs["model_name"] = model_name
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}
    configs = parse_configs(
        configs=configs,
        fill_required_with_null=True,
    )
    if seeds is not None:
        configs.seeds = seeds
    model_name = configs.model_name
    _, model_size, model_feature, model_version = model_name.split("_")
    model_specfics_configs = ConfigDict(model_configs[model_name])
    # update model specific configs
    configs.update(model_specfics_configs)
    # the user input configs has the highest priority
    configs.model.N_cycle = n_cycle
    configs.sample_diffusion.N_sample = n_sample
    configs.sample_diffusion.N_step = n_step
    configs.dtype = dtype
    configs.use_msa = use_msa
    configs.triangle_multiplicative = trimul_kernel
    configs.triangle_attention = triatt_kernel
    configs.enable_diffusion_shared_vars_cache = enable_cache
    configs.enable_efficient_fusion = enable_fusion
    configs.enable_tf32 = enable_tf32
    configs.use_template = use_template
    configs.use_rna_msa = use_rna_msa
    configs.use_seeds_in_json = use_seeds_in_json

    configs = update_gpu_compatible_configs(configs)
    logger.info(
        f"Inference by Protenix: model_size: {model_size}, "
        f"with_feature: {model_feature.replace('-', ',')}, "
        f"model_version: {model_version}, dtype: {configs.dtype}"
    )
    logger.info(
        f"Triangle_multiplicative kernel: {trimul_kernel}, "
        f"Triangle_attention kernel: {triatt_kernel}"
    )
    logger.info(
        f"enable_diffusion_shared_vars_cache: {configs.enable_diffusion_shared_vars_cache}, "
        f"enable_efficient_fusion: {configs.enable_efficient_fusion}, "
        f"enable_tf32: {configs.enable_tf32}"
    )
    download_inference_cache(configs)
    return InferenceRunner(configs)


def inference_jsons(
    json_file: str,
    out_dir: str = "./output",
    use_msa: bool = True,
    seeds: list = [101],
    n_cycle: int = 10,
    n_step: int = 200,
    n_sample: int = 5,
    dtype: str = "bf16",
    model_name: str = "protenix_base_default_v1.0.0",
    trimul_kernel: str = "cuequivariance",
    triatt_kernel: str = "cuequivariance",
    enable_cache: bool = True,
    enable_fusion: bool = True,
    enable_tf32: bool = True,
    msa_server_mode: str = "protenix",
    use_template: bool = False,
    use_rna_msa: bool = False,
    use_seeds_in_json: bool = False,
    hmmsearch_binary_path: Optional[str] = None,
    hmmbuild_binary_path: Optional[str] = None,
    seqres_database_path: Optional[str] = None,
    nhmmer_binary_path: Optional[str] = None,
    hmmalign_binary_path: Optional[str] = None,
    hmmbuild_rna_binary_path: Optional[str] = None,
    ntrna_database_path: Optional[str] = None,
    rfam_database_path: Optional[str] = None,
    rna_central_database_path: Optional[str] = None,
    nhmmer_n_cpu: Optional[int] = None,
) -> None:
    """
    Run inference on a single JSON file or a directory of JSON files.

    Args:
        json_file (str): Path to a JSON file or directory containing JSON files.
        out_dir (str): Directory to save inference results.
        use_msa (bool): Whether to use MSA.
        seeds (list): List of inference seeds.
        n_cycle (int): Number of cycles.
        n_step (int): Number of diffusion steps.
        n_sample (int): Number of samples.
        dtype (str): Data type.
        model_name (str): Model name.
        trimul_kernel (str): Kernel for triangle multiplicative.
        triatt_kernel (str): Kernel for triangle attention.
        enable_cache (bool): Enable shared variables cache.
        enable_fusion (bool): Enable efficient fusion.
        enable_tf32 (bool): Enable TF32.
        msa_server_mode (str): MSA server mode.
        use_template (bool): Whether to use templates.
        use_rna_msa (bool): Whether to use RNA MSA.
        use_seeds_in_json (bool): Whether to use seeds from JSON.
        hmmsearch_binary_path (Optional[str]): Path to hmmsearch binary.
        hmmbuild_binary_path (Optional[str]): Path to hmmbuild binary.
        seqres_database_path (Optional[str]): Path to sequence database.
        nhmmer_binary_path (Optional[str]): Path to nhmmer binary.
        hmmalign_binary_path (Optional[str]): Path to hmmalign binary.
        hmmbuild_rna_binary_path (Optional[str]): Path to RNA hmmbuild binary.
        ntrna_database_path (Optional[str]): NT-RNA database path.
        rfam_database_path (Optional[str]): Rfam database path.
        rna_central_database_path (Optional[str]): RNAcentral database path.
        nhmmer_n_cpu (Optional[int]): Number of CPUs for nhmmer.
    """
    infer_jsons = []
    if os.path.isdir(json_file):
        infer_jsons = [
            str(file) for file in Path(json_file).rglob("*") if file.is_file()
        ]
        if len(infer_jsons) == 0:
            raise RuntimeError(f"Can not read a valid json file in {json_file}")
    elif os.path.isfile(json_file):
        infer_jsons = [json_file]
    else:
        raise RuntimeError(f"Can not read a special file: {json_file}")
    infer_jsons = [file for file in infer_jsons if file.endswith(".json")]
    logger.info(f"Will infer with {len(infer_jsons)} jsons")
    if len(infer_jsons) == 0:
        return

    infer_errors = {}
    inference_configs["dump_dir"] = out_dir
    runner = get_default_runner(
        seeds=seeds,
        n_cycle=n_cycle,
        n_step=n_step,
        n_sample=n_sample,
        dtype=dtype,
        model_name=model_name,
        use_msa=use_msa,
        trimul_kernel=trimul_kernel,
        triatt_kernel=triatt_kernel,
        enable_cache=enable_cache,
        enable_fusion=enable_fusion,
        enable_tf32=enable_tf32,
        use_template=use_template,
        use_rna_msa=use_rna_msa,
        use_seeds_in_json=use_seeds_in_json,
    )
    configs = runner.configs
    for _, infer_json in enumerate(tqdm.tqdm(infer_jsons)):
        try:
            configs["input_json_path"] = preprocess_input(
                infer_json,
                out_dir=out_dir,
                use_msa=use_msa,
                use_template=use_template,
                use_rna_msa=use_rna_msa,
                msa_server_mode=msa_server_mode,
                hmmsearch_binary_path=hmmsearch_binary_path,
                hmmbuild_binary_path=hmmbuild_binary_path,
                seqres_database_path=seqres_database_path,
                nhmmer_binary_path=nhmmer_binary_path,
                hmmalign_binary_path=hmmalign_binary_path,
                hmmbuild_rna_binary_path=hmmbuild_rna_binary_path,
                ntrna_database_path=ntrna_database_path,
                rfam_database_path=rfam_database_path,
                rna_central_database_path=rna_central_database_path,
                nhmmer_n_cpu=nhmmer_n_cpu,
            )
            infer_predict(runner, configs)
        except Exception as exc:
            infer_errors[infer_json] = str(exc)
    if len(infer_errors) > 0:
        logger.warning(f"Run inference failed: {infer_errors}")


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"], show_default=True)


class SuggestGroup(click.Group):
    """A Click group that suggests similar commands on error."""

    def resolve_command(self, ctx, args):
        """Try to resolve the command, and suggest matches if it fails."""
        try:
            return super().resolve_command(ctx, args)
        except click.UsageError as e:
            if len(args) > 0:
                cmd_name = args[0]
                all_commands = self.list_commands(ctx)
                matches = difflib.get_close_matches(cmd_name, all_commands)
                if matches:
                    e.message += (
                        f"\n\nDid you mean one of these?\n    {', '.join(matches)}"
                    )
            raise e


@click.group(cls=SuggestGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(version="1.0.0")
def protenix_cli() -> None:
    """
    Protenix: A trainable reproduction of AlphaFold 3.

    This CLI provides tools for structure prediction, data conversion,
    and MSA/template searching.
    """
    pass


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i", "--input", type=str, required=True, help="Input JSON file or directory."
)
@click.option("-o", "--out_dir", default="./output", type=str, help="Output directory.")
@click.option("-s", "--seeds", type=str, default="101", help="Seeds (comma-separated).")
@click.option("-c", "--cycle", type=int, default=10, help="Pairformer cycle number.")
@click.option("-p", "--step", type=int, default=200, help="Diffusion steps.")
@click.option("-e", "--sample", type=int, default=5, help="Number of samples.")
@click.option("-d", "--dtype", type=str, default="bf16", help="Inference dtype.")
@click.option(
    "-n",
    "--model_name",
    type=str,
    default="protenix_base_default_v1.0.0",
    help="Model checkpoint name.",
)
@click.option(
    "--use_msa",
    type=bool,
    default=True,
    help="Whether to use MSA for inference.",
)
@click.option(
    "--use_default_params",
    type=bool,
    default=False,
    help="Use recommended default parameters for the selected model.",
)
@click.option(
    "--trimul_kernel",
    type=str,
    default="cuequivariance",
    help="Triangle multiplicative update kernel ('cuequivariance' or 'torch').",
)
@click.option(
    "--triatt_kernel",
    type=str,
    default="cuequivariance",
    help=(
        "Triangle attention kernel ('triattention', 'cuequivariance', "
        "'deepspeed', or 'torch')."
    ),
)
@click.option(
    "--enable_cache",
    type=bool,
    default=True,
    help="Cache shareable variables in the diffusion module.",
)
@click.option(
    "--enable_fusion",
    type=bool,
    default=True,
    help="Enable efficient kernel fusion in the diffusion transformer.",
)
@click.option(
    "--enable_tf32",
    type=bool,
    default=True,
    help="Enable TF32 for FP32 matrix multiplications.",
)
@click.option(
    "--msa_server_mode",
    type=str,
    default="protenix",
    help="MSA search mode ('protenix' or 'colabfold').",
)
@click.option(
    "--use_template",
    type=bool,
    default=False,
    help="Use templates (requires templatesPath in input JSON).",
)
@click.option(
    "--use_rna_msa",
    type=bool,
    default=False,
    help="Use RNA MSA (requires rna_msa_path in input JSON).",
)
@click.option(
    "--use_seeds_in_json",
    type=bool,
    default=False,
    help="Priority to seeds defined in input JSON.",
)
@click.option(
    "--hmmsearch_binary_path",
    type=str,
    default=None,
    help="Path to hmmsearch (searches in PATH if not provided).",
)
@click.option(
    "--hmmbuild_binary_path",
    type=str,
    default=None,
    help="Path to hmmbuild (searches in PATH if not provided).",
)
@click.option(
    "--seqres_database_path",
    type=str,
    default=None,
    help="Path to the sequence database for template search.",
)
@click.option(
    "--nhmmer_binary_path",
    type=str,
    default=None,
    help="Path to nhmmer for RNA MSA search.",
)
@click.option(
    "--hmmalign_binary_path",
    type=str,
    default=None,
    help="Path to hmmalign for RNA MSA search.",
)
@click.option(
    "--hmmbuild_rna_binary_path",
    type=str,
    default=None,
    help="Path to RNA-specific hmmbuild.",
)
@click.option(
    "--ntrna_database_path",
    type=str,
    default=None,
    help="Path to the NT-RNA database.",
)
@click.option(
    "--rfam_database_path",
    type=str,
    default=None,
    help="Path to the Rfam database.",
)
@click.option(
    "--rna_central_database_path",
    type=str,
    default=None,
    help="Path to the RNAcentral database.",
)
@click.option(
    "--nhmmer_n_cpu",
    type=int,
    default=None,
    help="Number of CPUs for nhmmer.",
)
def predict(
    input: str,
    out_dir: str,
    seeds: str,
    cycle: int,
    step: int,
    sample: int,
    dtype: str,
    model_name: str,
    use_msa: bool,
    use_default_params: bool,
    trimul_kernel: str,
    triatt_kernel: str,
    enable_cache: bool,
    enable_fusion: bool,
    enable_tf32: bool,
    msa_server_mode: str,
    use_template: bool,
    use_rna_msa: bool,
    use_seeds_in_json: bool,
    hmmsearch_binary_path: Optional[str] = None,
    hmmbuild_binary_path: Optional[str] = None,
    seqres_database_path: Optional[str] = None,
    nhmmer_binary_path: Optional[str] = None,
    hmmalign_binary_path: Optional[str] = None,
    hmmbuild_rna_binary_path: Optional[str] = None,
    ntrna_database_path: Optional[str] = None,
    rfam_database_path: Optional[str] = None,
    rna_central_database_path: Optional[str] = None,
    nhmmer_n_cpu: Optional[int] = None,
) -> None:
    """
    Run predictions with Protenix using various input formats.

    Args:
        input (str): Input JSON file or directory.
        out_dir (str): Output directory for results.
        seeds (str): Comma-separated seeds.
        cycle (int): Number of cycles.
        step (int): Number of diffusion steps.
        sample (int): Number of samples.
        dtype (str): Data type.
        model_name (str): Model name.
        use_msa (bool): Use MSA.
        use_default_params (bool): Use default parameters for the model.
        trimul_kernel (str): Kernel for triangle multiplicative.
        triatt_kernel (str): Kernel for triangle attention.
        enable_cache (bool): Enable shared variables cache.
        enable_fusion (bool): Enable efficient fusion.
        enable_tf32 (bool): Enable TF32.
        msa_server_mode (str): MSA server mode.
        use_template (bool): Use templates.
        use_rna_msa (bool): Use RNA MSA.
        use_seeds_in_json (bool): Use seeds from JSON.
        hmmsearch_binary_path (Optional[str]): Path to hmmsearch binary.
        hmmbuild_binary_path (Optional[str]): Path to hmmbuild binary.
        seqres_database_path (Optional[str]): Path to sequence database.
        nhmmer_binary_path (Optional[str]): Path to nhmmer binary.
        hmmalign_binary_path (Optional[str]): Path to hmmalign binary.
        hmmbuild_rna_binary_path (Optional[str]): Path to RNA hmmbuild binary.
        ntrna_database_path (Optional[str]): NT-RNA database path.
        rfam_database_path (Optional[str]): Rfam database path.
        rna_central_database_path (Optional[str]): RNAcentral database path.
        nhmmer_n_cpu (Optional[int]): Number of CPUs for nhmmer.
    """
    init_logging()
    logger.info(f"Run infer with input={input}, out_dir={out_dir}, sample={sample}")
    if use_default_params:
        if model_name in [
            "protenix_base_default_v0.5.0",
            "protenix_base_constraint_v0.5.0",
            "protenix_base_default_v1.0.0",
            "protenix_base_20250630_v1.0.0",
        ]:
            cycle = 10
            step = 200
        elif model_name in [
            "protenix_mini_esm_v0.5.0",
            "protenix_mini_ism_v0.5.0",
            "protenix_mini_default_v0.5.0",
            "protenix_tiny_default_v0.5.0",
        ]:
            cycle = 4
            step = 5
            if model_name in [
                "protenix_mini_esm_v0.5.0",
                "protenix_mini_ism_v0.5.0",
            ]:
                use_msa = False
        else:
            raise RuntimeError(
                f"{model_name} is not supported for inference in our model list"
            )
    logger.info(
        f"Using default params for model {model_name}: "
        f"cycle={cycle}, step={step}, use_msa={use_msa}"
    )
    assert trimul_kernel in [
        "cuequivariance",
        "torch",
    ], "Invalid trimul_kernel. Options: 'cuequivariance', 'torch'."
    assert triatt_kernel in ["triattention", "cuequivariance", "deepspeed", "torch",], (
        "Invalid triatt_kernel. Options: 'triattention', "
        "'cuequivariance', 'deepspeed', 'torch'."
    )
    seeds = list(map(int, seeds.split(",")))

    if use_template:
        assert model_name in [
            "protenix_base_default_v1.0.0",
            "protenix_base_20250630_v1.0.0",
        ], "Only protenix_base_default_v1.0.0 and protenix_base_20250630_v1.0.0 supports template inference."
        logger.info("=" * 50)
        logger.info(
            "Using templates for inference. Template files should have "
            ".hrr or .a3m extensions and be specified in the JSON file.\n"
            "Example: /path/to/template.hrr or /path/to/template.a3m\n"
            "Note: Inference will proceed with automatic template search "
            "if none are provided and use_template is True."
        )
        logger.info("=" * 50)

    if use_rna_msa:
        assert model_name in [
            "protenix_base_default_v1.0.0",
            "protenix_base_20250630_v1.0.0",
        ], "Only protenix_base_default_v1.0.0 and protenix_base_20250630_v1.0.0 supports RNA MSA inference."
        logger.info("=" * 50)
        logger.info(
            "Using RNA MSA for inference. RNA MSA files should have .a3m "
            "extension and be specified in the JSON file.\n"
            "Example: /path/to/rna_msa.a3m\n"
            "Note: Inference will proceed with automatic RNA MSA search "
            "if none are provided and use_rna_msa is True."
        )
        logger.info("=" * 50)

    if use_seeds_in_json:
        logger.info("=" * 50)
        logger.info(
            "Using seeds defined in JSON file for inference.\n"
            "Note: This will override any seeds passed via command line, "
            "using seeds from modelSeeds defined in the JSON."
        )
        logger.info("=" * 50)

    inference_jsons(
        input,
        out_dir,
        use_msa,
        seeds=seeds,
        n_cycle=cycle,
        n_step=step,
        n_sample=sample,
        dtype=dtype,
        model_name=model_name,
        trimul_kernel=trimul_kernel,
        triatt_kernel=triatt_kernel,
        enable_cache=enable_cache,
        enable_fusion=enable_fusion,
        enable_tf32=enable_tf32,
        msa_server_mode=msa_server_mode,
        use_template=use_template,
        use_rna_msa=use_rna_msa,
        use_seeds_in_json=use_seeds_in_json,
        hmmsearch_binary_path=hmmsearch_binary_path,
        hmmbuild_binary_path=hmmbuild_binary_path,
        seqres_database_path=seqres_database_path,
        nhmmer_binary_path=nhmmer_binary_path,
        hmmalign_binary_path=hmmalign_binary_path,
        hmmbuild_rna_binary_path=hmmbuild_rna_binary_path,
        ntrna_database_path=ntrna_database_path,
        rfam_database_path=rfam_database_path,
        rna_central_database_path=rna_central_database_path,
        nhmmer_n_cpu=nhmmer_n_cpu,
    )


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input",
    type=str,
    required=True,
    help="PDB/CIF files or directory to generate inference JSONs.",
)
@click.option("-o", "--out_dir", type=str, default="./output", help="Output directory.")
@click.option(
    "--altloc",
    default="first",
    type=str,
    help=(
        "Select the first altloc conformation of each residue, "
        "or specify the altloc letter ('A', 'B', etc.)."
    ),
)
@click.option(
    "--assembly_id",
    default=None,
    type=str,
    help="Assembly ID for structure extension (default: no extension).",
)
def tojson(
    input: str,
    out_dir: str = "./output",
    altloc: str = "first",
    assembly_id: Optional[str] = None,
) -> List[str]:
    """
    Convert PDB or CIF files to JSON files for Protenix inference.

    Args:
        input (str): Input PDB/CIF file or directory.
        out_dir (str): Output directory for JSON files.
        altloc (str): Alternate location conformation selection.
        assembly_id (Optional[str]): Assembly ID for structure extension.

    Returns:
        List[str]: List of generated JSON file paths.
    """
    init_logging()
    logger.info(
        f"Run tojson with input={input}, out_dir={out_dir}, "
        f"altloc={altloc}, assembly_id={assembly_id}"
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
        output_json = os.path.join(out_dir, f"{pdb_name}.json")
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
    logger.info(f"{len(output_jsons)} generated jsons have been saved to {out_dir}.")
    return output_jsons


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input",
    type=str,
    required=True,
    help="JSON or FASTA file for MSA search.",
)
@click.option("-o", "--out_dir", type=str, default="./output", help="Output directory.")
@click.option(
    "-m",
    "--msa_server_mode",
    type=str,
    default="protenix",
    help="MSA search mode ('protenix' or 'colabfold').",
)
def msa(input: str, out_dir: str, msa_server_mode: str) -> Union[str, dict]:
    """
    Perform MSA search using MMseqs2.
    If input is a FASTA file, it should contain protein sequences.

    Args:
        input (str): Path to a JSON or FASTA file.
        out_dir (str): Directory to save MSA results.
        msa_server_mode (str): MSA search mode ('protenix' or 'colabfold').

    Returns:
        Union[str, dict]: Updated JSON path or dictionary of MSA results.
    """
    init_logging()
    logger.info(f"Run msa with input={input}, out_dir={out_dir}")
    if input.endswith(".json"):
        msa_input_json, _ = update_infer_json(
            input, out_dir, use_msa=True, mode=msa_server_mode
        )
        logger.info(f"msa results have been update to {msa_input_json}")
        return msa_input_json
    elif input.endswith(".fasta"):
        records = list(SeqIO.parse(input, "fasta"))
        protein_seqs = []
        for seq in records:
            protein_seqs.append(str(seq.seq))
        protein_seqs = sorted(protein_seqs)
        msa_res_subdirs = msa_search(protein_seqs, out_dir, msa_server_mode)
        assert len(msa_res_subdirs) == len(msa_res_subdirs), "msa search failed"
        fasta_msa_res = dict(zip(protein_seqs, msa_res_subdirs))
        logger.info(
            f"msa result is: {fasta_msa_res}, and it has been save to {out_dir}"
        )
        return fasta_msa_res
    else:
        raise RuntimeError(f"only support `json` or `fasta` format, but got : {input}")


# The new msatemplate command first performs MSA search, then performs template search
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input",
    type=str,
    required=True,
    help="JSON file for MSA and template search.",
)
@click.option(
    "-o",
    "--out_dir",
    type=str,
    default="./output",
    help="Output directory.",
)
@click.option(
    "-m",
    "--msa_server_mode",
    type=str,
    default="protenix",
    help="MSA search mode ('protenix' or 'colabfold').",
)
@click.option(
    "--hmmsearch_binary_path",
    type=str,
    default=None,
    help="Path to hmmsearch (searches in PATH if not provided).",
)
@click.option(
    "--hmmbuild_binary_path",
    type=str,
    default=None,
    help="Path to hmmbuild (searches in PATH if not provided).",
)
@click.option(
    "--seqres_database_path",
    type=str,
    default=None,
    help="Path to the sequence database for template search.",
)
def msatemplate(
    input: str,
    out_dir: str,
    msa_server_mode: str,
    hmmsearch_binary_path: Optional[str],
    hmmbuild_binary_path: Optional[str],
    seqres_database_path: Optional[str],
) -> str:
    """
    Perform MSA search followed by template search.

    Args:
        input (str): Path to the input JSON file.
        out_dir (str): Directory to save MSA and template results.
        msa_server_mode (str): MSA search mode ('protenix' or 'colabfold').
        hmmsearch_binary_path (Optional[str]): Path to hmmsearch binary.
        hmmbuild_binary_path (Optional[str]): Path to hmmbuild binary.
        seqres_database_path (Optional[str]): Path to sequence database.

    Returns:
        str: Updated JSON file path with template information.
    """
    logger.info(f"Run msa_template with input={input}, out_dir={out_dir}")

    if not input.endswith(".json"):
        raise RuntimeError(
            f"msa_template only supports `json` format, but got: {input}"
        )

    if not os.path.exists(input):
        raise RuntimeError(f"input file {input} does not exist")

    return preprocess_input(
        input_json=input,
        out_dir=out_dir,
        use_msa=True,
        use_template=True,
        use_rna_msa=False,
        msa_server_mode=msa_server_mode,
        hmmsearch_binary_path=hmmsearch_binary_path,
        hmmbuild_binary_path=hmmbuild_binary_path,
        seqres_database_path=seqres_database_path,
    )


# The new inputprep command calls the RNA MSA process after the MSA template process finishes
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input",
    type=str,
    required=True,
    help="JSON file to update with RNA MSA (supports 'json' format only).",
)
@click.option(
    "-o",
    "--out_dir",
    type=str,
    default="./output",
    help="Output directory.",
)
@click.option(
    "-m",
    "--msa_server_mode",
    type=str,
    default="protenix",
    help="MSA search mode ('protenix' or 'colabfold').",
)
@click.option(
    "--hmmsearch_binary_path",
    type=str,
    default=None,
    help="Path to hmmsearch.",
)
@click.option(
    "--hmmbuild_binary_path",
    type=str,
    default=None,
    help="Path to hmmbuild.",
)
@click.option(
    "--seqres_database_path",
    type=str,
    default=None,
    help="Path to the sequence database for template search.",
)
@click.option(
    "--nhmmer_binary_path",
    type=str,
    default=None,
    help="Path to nhmmer for RNA MSA search.",
)
@click.option(
    "--hmmalign_binary_path",
    type=str,
    default=None,
    help="Path to hmmalign for RNA MSA search.",
)
@click.option(
    "--hmmbuild_rna_binary_path",
    type=str,
    default=None,
    help="Path to RNA-specific hmmbuild.",
)
@click.option(
    "--ntrna_database_path",
    type=str,
    default=None,
    help="Path to the NT-RNA database.",
)
@click.option(
    "--rfam_database_path",
    type=str,
    default=None,
    help="Path to the Rfam database.",
)
@click.option(
    "--rna_central_database_path",
    type=str,
    default=None,
    help="Path to the RNAcentral database.",
)
@click.option(
    "--nhmmer_n_cpu",
    type=int,
    default=None,
    help="Number of CPUs for nhmmer.",
)
def inputprep(
    input: str,
    out_dir: str,
    msa_server_mode: str,
    hmmsearch_binary_path: Optional[str],
    hmmbuild_binary_path: Optional[str],
    seqres_database_path: Optional[str],
    nhmmer_binary_path: Optional[str],
    hmmalign_binary_path: Optional[str],
    hmmbuild_rna_binary_path: Optional[str],
    ntrna_database_path: Optional[str],
    rfam_database_path: Optional[str],
    rna_central_database_path: Optional[str],
    nhmmer_n_cpu: Optional[int],
) -> str:
    """
    Perform MSA search, template search, and RNA MSA search sequentially.

    Args:
        input (str): Path to the input JSON file.
        out_dir (str): Directory to save all search results.
        msa_server_mode (str): MSA search mode ('protenix' or 'colabfold').
        hmmsearch_binary_path (Optional[str]): Path to hmmsearch binary.
        hmmbuild_binary_path (Optional[str]): Path to hmmbuild binary.
        seqres_database_path (Optional[str]): Path to sequence database.
        nhmmer_binary_path (Optional[str]): Path to nhmmer binary.
        hmmalign_binary_path (Optional[str]): Path to hmmalign binary.
        hmmbuild_rna_binary_path (Optional[str]): Path to RNA hmmbuild binary.
        ntrna_database_path (Optional[str]): Path to NT-RNA database.
        rfam_database_path (Optional[str]): Path to Rfam database.
        rna_central_database_path (Optional[str]): Path to RNAcentral database.
        nhmmer_n_cpu (Optional[int]): Number of CPUs for nhmmer.

    Returns:
        str: Final updated JSON file path with all search information.
    """
    logger.info(f"Run inputprep with input={input}, out_dir={out_dir}")

    if not input.endswith(".json"):
        raise RuntimeError(f"inputprep only supports `json` format, but got: {input}")

    if not os.path.exists(input):
        raise RuntimeError(f"input file {input} does not exist")

    return preprocess_input(
        input_json=input,
        out_dir=out_dir,
        use_msa=True,
        use_template=True,
        use_rna_msa=True,
        msa_server_mode=msa_server_mode,
        hmmsearch_binary_path=hmmsearch_binary_path,
        hmmbuild_binary_path=hmmbuild_binary_path,
        seqres_database_path=seqres_database_path,
        nhmmer_binary_path=nhmmer_binary_path,
        hmmalign_binary_path=hmmalign_binary_path,
        hmmbuild_rna_binary_path=hmmbuild_rna_binary_path,
        ntrna_database_path=ntrna_database_path,
        rfam_database_path=rfam_database_path,
        rna_central_database_path=rna_central_database_path,
        nhmmer_n_cpu=nhmmer_n_cpu,
    )


protenix_cli.add_command(predict, name="pred")
protenix_cli.add_command(tojson, name="json")
protenix_cli.add_command(msa, name="msa")
protenix_cli.add_command(msatemplate, name="mt")
protenix_cli.add_command(inputprep, name="prep")

if __name__ == "__main__":
    predict()
