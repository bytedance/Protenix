import os, json, logging, uuid, time, tqdm, argparse
from pathlib import Path
from rdkit import Chem
from typing import Any, List
from protenix.data.json_parser import lig_file_to_atom_info
from runner.inference import download_infercence_cache, InferenceRunner, infer_predict
from configs.configs_base import configs as configs_base
from configs.configs_data import data_configs
from configs.configs_inference import inference_configs
from protenix.config import parse_configs, parse_sys_args
from protenix.utils.logger import get_logger

logger = get_logger(__name__)


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
        one_infer_json = [
            {"sequences": one_infer_seq, "modelSeeds": seeds, "name": ligand_name}
        ]
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
        one_infer_json = [
            {"sequences": one_infer_seq, "modelSeeds": seeds, "name": ligand_name}
        ]
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


def get_default_runner() -> InferenceRunner:
    inference_configs["load_checkpoint_path"] = "/af3-dev/release_model/model_v1.pt"
    configs_base["use_deepspeed_evo_attention"] = os.environ.get(
        "use_deepspeed_evo_attention", False
    )
    configs_base["model"]["N_cycle"] = 10
    configs_base["sample_diffusion"]["N_sample"] = 5
    configs_base["sample_diffusion"]["N_step"] = 200
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}
    configs = parse_configs(
        configs=configs,
        fill_required_with_null=True,
    )
    download_infercence_cache(configs)
    return InferenceRunner(configs)


def inference_jsons(json_file: str, out_dir: str = "./output") -> None:
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
    runner = get_default_runner()
    configs = runner.configs
    for infer_json in tqdm.tqdm(infer_jsons):
        try:
            configs["input_json_path"] = infer_json
            infer_predict(runner, configs)
        except Exception as exc:
            infer_errors[infer_json] = str(exc)
    logger.info(f"run inference failed jsons: {infer_errors}")


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
    runner = get_default_runner()
    configs = runner.configs
    for infer_json in tqdm.tqdm(infer_jsons):
        try:
            configs["input_json_path"] = infer_json
            infer_predict(runner, configs)
        except Exception as exc:
            infer_errors[infer_json] = str(exc)
    logger.info(f"run inference failed jsons: {infer_errors}")


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


def main():
    LOG_FORMAT = "%(asctime)s,%(msecs)-3d %(levelname)-8s [%(filename)s:%(lineno)s %(funcName)s] %(message)s"
    logging.basicConfig(
        format=LOG_FORMAT,
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        filemode="w",
    )
    parser = argparse.ArgumentParser(description="infer with jsons of argparse")
    parser.add_argument("--input_json_path", required=True, type=str)
    parser.add_argument("--dump_dir", default="./output", type=str)
    args = parser.parse_args()
    logger.info(
        f"run infer with input_json_path={args.input_json_path}, dump_dir={args.dump_dir}"
    )
    inference_jsons(args.input_json_path, args.dump_dir)


if __name__ == "__main__":
    LOG_FORMAT = "%(asctime)s,%(msecs)-3d %(levelname)-8s [%(filename)s:%(lineno)s %(funcName)s] %(message)s"
    logging.basicConfig(
        format=LOG_FORMAT,
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        filemode="w",
    )
    test_batch_inference()
