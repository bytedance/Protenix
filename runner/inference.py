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
import traceback
import urllib.request
from contextlib import nullcontext
from os.path import exists as opexists
from os.path import join as opjoin
from typing import Any, Mapping

import torch
import torch.distributed as dist
from configs.configs_base import configs as configs_base
from configs.configs_data import data_configs
from configs.configs_inference import inference_configs
from runner.dumper import DataDumper

from protenix.config import parse_configs, parse_sys_args
from protenix.data.infer_data_pipeline import get_inference_dataloader
from protenix.model.protenix import Protenix
from protenix.utils.distributed import DIST_WRAPPER
from protenix.utils.seed import seed_everything
from protenix.utils.torch_utils import to_device
from protenix.web_service.dependency_url import URL

logger = logging.getLogger(__name__)


class InferenceRunner(object):
    def __init__(self, configs: Any) -> None:
        self.configs = configs
        self.init_env()
        self.init_basics()
        self.init_model()
        self.load_checkpoint()
        self.init_dumper(
            need_atom_confidence=configs.need_atom_confidence,
            sorted_by_ranking_score=configs.sorted_by_ranking_score,
        )

    def init_env(self) -> None:
        self.print(
            f"Distributed environment: world size: {DIST_WRAPPER.world_size}, "
            + f"global rank: {DIST_WRAPPER.rank}, local rank: {DIST_WRAPPER.local_rank}"
        )
        self.use_cuda = torch.cuda.device_count() > 0
        if self.use_cuda:
            self.device = torch.device("cuda:{}".format(DIST_WRAPPER.local_rank))
            os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
            all_gpu_ids = ",".join(str(x) for x in range(torch.cuda.device_count()))
            devices = os.getenv("CUDA_VISIBLE_DEVICES", all_gpu_ids)
            logging.info(
                f"LOCAL_RANK: {DIST_WRAPPER.local_rank} - CUDA_VISIBLE_DEVICES: [{devices}]"
            )
            torch.cuda.set_device(self.device)
        else:
            self.device = torch.device("cpu")
        if DIST_WRAPPER.world_size > 1:
            dist.init_process_group(backend="nccl")
        if self.configs.use_deepspeed_evo_attention:
            env = os.getenv("CUTLASS_PATH", None)
            self.print(f"env: {env}")
            assert (
                env is not None
            ), "if use ds4sci, set `CUTLASS_PATH` env as https://www.deepspeed.ai/tutorials/ds4sci_evoformerattention/"
            if env is not None:
                logging.info(
                    "The kernels will be compiled when DS4Sci_EvoformerAttention is called for the first time."
                )
        use_fastlayernorm = os.getenv("LAYERNORM_TYPE", None)
        if use_fastlayernorm == "fast_layernorm":
            logging.info(
                "The kernels will be compiled when fast_layernorm is called for the first time."
            )

        logging.info("Finished init ENV.")

    def init_basics(self) -> None:
        self.dump_dir = self.configs.dump_dir
        self.error_dir = opjoin(self.dump_dir, "ERR")
        os.makedirs(self.dump_dir, exist_ok=True)
        os.makedirs(self.error_dir, exist_ok=True)

    def init_model(self) -> None:
        self.model = Protenix(self.configs).to(self.device)

    def load_checkpoint(self) -> None:
        checkpoint_path = self.configs.load_checkpoint_path
        if not os.path.exists(checkpoint_path):
            raise Exception(f"Given checkpoint path not exist [{checkpoint_path}]")
        self.print(
            f"Loading from {checkpoint_path}, strict: {self.configs.load_strict}"
        )
        checkpoint = torch.load(checkpoint_path, self.device)

        sample_key = [k for k in checkpoint["model"].keys()][0]
        self.print(f"Sampled key: {sample_key}")
        if sample_key.startswith("module."):  # DDP checkpoint has module. prefix
            checkpoint["model"] = {
                k[len("module.") :]: v for k, v in checkpoint["model"].items()
            }
        self.model.load_state_dict(
            state_dict=checkpoint["model"],
            strict=self.configs.load_strict,
        )
        self.model.eval()
        self.print(f"Finish loading checkpoint.")

    def init_dumper(
        self, need_atom_confidence: bool = False, sorted_by_ranking_score: bool = True
    ):
        self.dumper = DataDumper(
            base_dir=self.dump_dir,
            need_atom_confidence=need_atom_confidence,
            sorted_by_ranking_score=sorted_by_ranking_score,
        )

    # Adapted from runner.train.Trainer.evaluate
    @torch.no_grad()
    def predict(self, data: Mapping[str, Mapping[str, Any]]) -> dict[str, torch.Tensor]:
        eval_precision = {
            "fp32": torch.float32,
            "bf16": torch.bfloat16,
            "fp16": torch.float16,
        }[self.configs.dtype]

        enable_amp = (
            torch.autocast(device_type="cuda", dtype=eval_precision)
            if torch.cuda.is_available()
            else nullcontext()
        )

        data = to_device(data, self.device)
        with enable_amp:
            prediction, _, _ = self.model(
                input_feature_dict=data["input_feature_dict"],
                label_full_dict=None,
                label_dict=None,
                mode="inference",
            )

        return prediction

    def print(self, msg: str):
        if DIST_WRAPPER.rank == 0:
            logger.info(msg)

    def update_model_configs(self, new_configs: Any) -> None:
        self.model.configs = new_configs


def download_infercence_cache(configs: Any, model_version: str = "v0.2.0") -> None:
    current_file_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_file_path)
    code_directory = os.path.dirname(current_directory)

    data_cache_dir = os.path.join(code_directory, "release_data/ccd_cache")
    os.makedirs(data_cache_dir, exist_ok=True)
    for cache_name, fname in [
        ("ccd_components_file", "components.v20240608.cif"),
        ("ccd_components_rdkit_mol_file", "components.v20240608.cif.rdkit_mol.pkl"),
    ]:
        if not opexists(cache_path := os.path.abspath(opjoin(data_cache_dir, fname))):
            tos_url = URL[cache_name]
            logger.info(f"Downloading data cache from\n {tos_url}...")
            urllib.request.urlretrieve(tos_url, cache_path)

    checkpoint_path = configs.load_checkpoint_path

    if not opexists(checkpoint_path):
        checkpoint_path = os.path.join(
            code_directory, f"release_data/checkpoint/model_{model_version}.pt"
        )
        os.makedirs(os.path.dirname(checkpoint_path), exist_ok=True)
        tos_url = URL[f"model_{model_version}"]
        logger.info(f"Downloading model checkpoint from\n {tos_url}...")
        urllib.request.urlretrieve(tos_url, checkpoint_path)
        try:
            ckpt = torch.load(checkpoint_path)
            del ckpt
        except:
            os.remove(checkpoint_path)
            raise RuntimeError(
                "Download model checkpoint failed, please download by yourself with "
                f"wget {tos_url} -O {checkpoint_path}"
            )
        configs.load_checkpoint_path = checkpoint_path


def update_inference_configs(configs: Any, N_token: int):
    # Setting the default inference configs for different N_token and N_atom
    # when N_token is larger than 3000, the default config might OOM even on a
    # A100 80G GPUS,
    if N_token > 3840:
        configs.skip_amp.confidence_head = False
        configs.skip_amp.sample_diffusion = False
    elif N_token > 2560:
        configs.skip_amp.confidence_head = False
        configs.skip_amp.sample_diffusion = True
    else:
        configs.skip_amp.confidence_head = True
        configs.skip_amp.sample_diffusion = True
    return configs


def infer_predict(runner: InferenceRunner, configs: Any) -> None:
    # Data
    logger.info(f"Loading data from\n{configs.input_json_path}")
    dataloader = get_inference_dataloader(configs=configs)

    num_data = len(dataloader.dataset)
    for seed in configs.seeds:
        seed_everything(seed=seed, deterministic=configs.deterministic)
        for batch in dataloader:
            try:
                data, atom_array, data_error_message = batch[0]

                if len(data_error_message) > 0:
                    logger.info(data_error_message)
                    with open(
                        opjoin(runner.error_dir, f"{data['sample_name']}.txt"),
                        "w",
                    ) as f:
                        f.write(data_error_message)
                    continue

                sample_name = data["sample_name"]
                logger.info(
                    (
                        f"[Rank {DIST_WRAPPER.rank} ({data['sample_index'] + 1}/{num_data})] {sample_name}: "
                        f"N_asym {data['N_asym'].item()}, N_token {data['N_token'].item()}, "
                        f"N_atom {data['N_atom'].item()}, N_msa {data['N_msa'].item()}"
                    )
                )
                new_configs = update_inference_configs(configs, data["N_token"].item())
                runner.update_model_configs(new_configs)
                prediction = runner.predict(data)
                runner.dumper.dump(
                    dataset_name="",
                    pdb_id=sample_name,
                    seed=seed,
                    pred_dict=prediction,
                    atom_array=atom_array,
                    entity_poly_type=data["entity_poly_type"],
                )

                logger.info(
                    f"[Rank {DIST_WRAPPER.rank}] {data['sample_name']} succeeded.\n"
                    f"Results saved to {configs.dump_dir}"
                )
                torch.cuda.empty_cache()
            except Exception as e:
                error_message = f"[Rank {DIST_WRAPPER.rank}]{data['sample_name']} {e}:\n{traceback.format_exc()}"
                logger.info(error_message)
                # Save error info
                if opexists(
                    error_path := opjoin(runner.error_dir, f"{sample_name}.txt")
                ):
                    os.remove(error_path)
                with open(error_path, "w") as f:
                    f.write(error_message)
                if hasattr(torch.cuda, "empty_cache"):
                    torch.cuda.empty_cache()
                raise RuntimeError(f"run infer failed: {str(e)}")


def main(configs: Any) -> None:
    # Runner
    runner = InferenceRunner(configs)
    infer_predict(runner, configs)


def run() -> None:
    LOG_FORMAT = "%(asctime)s,%(msecs)-3d %(levelname)-8s [%(filename)s:%(lineno)s %(funcName)s] %(message)s"
    logging.basicConfig(
        format=LOG_FORMAT,
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        filemode="w",
    )
    configs_base["use_deepspeed_evo_attention"] = (
        os.environ.get("USE_DEEPSPEED_EVO_ATTTENTION", False) == "true"
    )
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}
    configs = parse_configs(
        configs=configs,
        arg_str=parse_sys_args(),
        fill_required_with_null=True,
    )
    download_infercence_cache(configs, model_version="v0.2.0")
    main(configs)


if __name__ == "__main__":
    run()
