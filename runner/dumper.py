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
from pathlib import Path

import numpy as np
import torch
from biotite.structure import AtomArray

from protenix.data.utils import save_structure_cif
from protenix.utils.file_io import save_json
from protenix.utils.torch_utils import round_values


def get_clean_full_confidence(full_confidence_dict: dict) -> dict:
    """
    Clean and format the full confidence dictionary by removing unnecessary keys and rounding values.

    Args:
        full_confidence_dict (dict): The dictionary containing full confidence data.

    Returns:
        dict: The cleaned and formatted dictionary.
    """
    # Remove atom_coordinate
    full_confidence_dict.pop("atom_coordinate")
    # Remove atom_is_polymer
    full_confidence_dict.pop("atom_is_polymer")
    # Keep two decimal places
    full_confidence_dict = round_values(full_confidence_dict)
    return full_confidence_dict


class DataDumper:
    def __init__(self, base_dir, need_atom_confidence: bool = False):
        self.base_dir = base_dir
        self.need_atom_confidence = need_atom_confidence

    def dump(
        self,
        dataset_name: str,
        pdb_id: str,
        seed: int,
        pred_dict: dict,
        atom_array: AtomArray,
        entity_poly_type: dict[str, str],
    ):
        """
        Dump the predictions and related data to the specified directory.

        Args:
            dataset_name (str): The name of the dataset.
            pdb_id (str): The PDB ID of the sample.
            seed (int): The seed used for randomization.
            pred_dict (dict): The dictionary containing the predictions.
            atom_array (AtomArray): The AtomArray object containing the structure data.
            entity_poly_type (dict[str, str]): The entity poly type information.
        """
        dump_dir = self._get_dump_dir(dataset_name, pdb_id, seed)
        Path(dump_dir).mkdir(parents=True, exist_ok=True)

        self.dump_predictions(
            pred_dict=pred_dict,
            dump_dir=dump_dir,
            pdb_id=pdb_id,
            atom_array=atom_array,
            entity_poly_type=entity_poly_type,
            seed=seed,
        )

    def _get_dump_dir(self, dataset_name: str, sample_name: str, seed: int) -> str:
        """
        Generate the directory path for dumping data based on the dataset name, sample name, and seed.
        """
        dump_dir = os.path.join(
            self.base_dir, dataset_name, sample_name, f"seed_{seed}"
        )
        return dump_dir

    def dump_predictions(
        self,
        pred_dict: dict,
        dump_dir: str,
        pdb_id: str,
        atom_array: AtomArray,
        entity_poly_type: dict[str, str],
        seed: int,
    ):
        """
        Dump raw predictions from the model:
            structure: Save the predicted coordinates as CIF files.
            confidence: Save the confidence data as JSON files.
        """
        prediction_save_dir = os.path.join(dump_dir, "predictions")
        os.makedirs(prediction_save_dir, exist_ok=True)

        # Dump structure
        b_factor = None
        if "full_data" in pred_dict:
            all_atom_plddt = []
            # len(pred_dict["full_data"]) == N_sample
            for each_sample_dict in pred_dict["full_data"]:
                if "atom_plddt" in each_sample_dict:
                    # atom_plddt.shape == [N_atom]
                    atom_plddt = each_sample_dict["atom_plddt"]
                    if atom_plddt.dtype == torch.bfloat16:
                        atom_plddt = atom_plddt.to(torch.float32)
                    all_atom_plddt.append(atom_plddt.cpu().numpy() * 100.0)

            if len(all_atom_plddt) == len(pred_dict["full_data"]):
                b_factor = all_atom_plddt

        self._save_structure(
            pred_coordinates=pred_dict["coordinate"],
            prediction_save_dir=prediction_save_dir,
            sample_name=pdb_id,
            atom_array=atom_array,
            entity_poly_type=entity_poly_type,
            seed=seed,
            b_factor=b_factor,
        )
        # Dump confidence
        self._save_confidence(
            data=pred_dict,
            prediction_save_dir=prediction_save_dir,
            sample_name=pdb_id,
            seed=seed,
        )

    def _save_structure(
        self,
        pred_coordinates: torch.Tensor,
        prediction_save_dir: str,
        sample_name: str,
        atom_array: AtomArray,
        entity_poly_type: dict[str, str],
        seed: int,
        b_factor: torch.Tensor = None,
    ):
        assert atom_array is not None
        N_sample = pred_coordinates.shape[0]
        for sample_idx in range(N_sample):
            output_fpath = os.path.join(
                prediction_save_dir,
                f"{sample_name}_seed_{seed}_sample_{sample_idx}.cif",
            )
            if b_factor is not None:
                # b_factor.shape == [N_sample, N_atom]
                atom_array.set_annotation("b_factor", np.round(b_factor[sample_idx], 2))

            save_structure_cif(
                atom_array=atom_array,
                pred_coordinate=pred_coordinates[sample_idx],
                output_fpath=output_fpath,
                entity_poly_type=entity_poly_type,
                pdb_id=sample_name,
            )

    def _save_confidence(
        self,
        data: dict,
        prediction_save_dir: str,
        sample_name: str,
        seed: int,
        sorted_by_ranking_score: bool = True,
    ):
        N_sample = len(data["summary_confidence"])
        for idx in range(N_sample):
            if self.need_atom_confidence:
                data["full_data"][idx] = get_clean_full_confidence(
                    data["full_data"][idx]
                )
        sorted_indices = range(N_sample)
        if sorted_by_ranking_score:
            sorted_indices = sorted(
                range(N_sample),
                key=lambda i: data["summary_confidence"][i]["ranking_score"],
                reverse=True,
            )

        for rank, idx in enumerate(sorted_indices):
            output_fpath = os.path.join(
                prediction_save_dir,
                f"{sample_name}_seed_{seed}_summary_confidence_sample_{rank}.json",
            )
            save_json(data["summary_confidence"][idx], output_fpath, indent=4)
            if self.need_atom_confidence:
                output_fpath = os.path.join(
                    prediction_save_dir,
                    f"{sample_name}_full_data_sample_{idx}.json",
                )
                save_json(data["full_data"][idx], output_fpath, indent=None)
