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

import traceback
from typing import Union

import torch

from protenix.utils.logger import get_logger
from protenix.utils.permutation.utils import save_permutation_error

from .heuristic import correct_symmetric_chains
from .pocket_based_permutation import permute_pred_to_optimize_pocket_aligned_rmsd

logger = get_logger(__name__)


def run(
    pred_coord: torch.Tensor,
    input_feature_dict: dict[str, Union[torch.Tensor, int, float, dict]],
    label_full_dict: dict[str, Union[torch.Tensor, int, float, dict]],
    max_num_chains: int = -1,
    permute_label: bool = True,
    permute_by_pocket: bool = False,
    error_dir: str = None,
    **kwargs,
) -> tuple[dict]:
    """
    Run chain permutation.


    Args:
        pred_coord (torch.Tensor): The predicted coordinates. Shape: [N_atoms, 3].
        input_feature_dict (dict[str, Union[torch.Tensor, int, float, dict]]): A dictionary containing input features.
        label_full_dict (dict[str, Union[torch.Tensor, int, float, dict]]): A dictionary containing full label information.
        max_num_chains (int, optional): The maximum number of chains to consider. Defaults to -1 (no limit).
        permute_label (bool, optional): Whether to permute the label. Defaults to True.
        permute_by_pocket (bool, optional): Whether to permute by pocket (for PoseBusters dataset). Defaults to False.
        error_dir (str, optional): Directory to save error data. Defaults to None.
        **kwargs: Additional keyword arguments.

    Returns:
        tuple[dict]: A tuple containing the output dictionary, log dictionary, permuted prediction indices, and permuted label indices.
    """

    if pred_coord.dim() > 2:
        assert (
            permute_label is False
        ), "Only supports prediction permutations in batch mode."

    try:

        if permute_by_pocket:
            """Optimize the chain assignment on pocket-ligand interface"""
            assert not permute_label

            if label_full_dict["pocket_mask"].dim() == 2:
                # first pocket is the `main` pocket
                pocket_mask = label_full_dict["pocket_mask"][0]
                ligand_mask = label_full_dict["interested_ligand_mask"][0]
            else:
                pocket_mask = label_full_dict["pocket_mask"]
                ligand_mask = label_full_dict["interested_ligand_mask"]

            permute_pred_indices, permuted_aligned_pred_coord, log_dict = (
                permute_pred_to_optimize_pocket_aligned_rmsd(
                    pred_coord=pred_coord,
                    true_coord=label_full_dict["coordinate"],
                    true_coord_mask=label_full_dict["coordinate_mask"],
                    true_pocket_mask=pocket_mask,
                    true_ligand_mask=ligand_mask,
                    atom_entity_id=input_feature_dict["entity_mol_id"],
                    atom_asym_id=input_feature_dict["mol_id"],
                    mol_atom_index=input_feature_dict["mol_atom_index"],
                    use_center_rmsd=kwargs.get("use_center_rmsd", False),
                )
            )
            output_dict = {"coordinate": permuted_aligned_pred_coord}
            permute_label_indices = []

        else:
            """Optimize the chain assignment on all chains"""
            output_dict, log_dict, permute_pred_indices, permute_label_indices = (
                correct_symmetric_chains(
                    pred_dict={**input_feature_dict, "coordinate": pred_coord},
                    label_full_dict=label_full_dict,
                    max_num_chains=max_num_chains,
                    permute_label=permute_label,
                    **kwargs,
                )
            )

    except Exception as e:
        error_message = f"{e}:\n{traceback.format_exc()}"
        logger.warning(error_message)
        save_permutation_error(
            data={
                "error_message": error_message,
                "pred_dict": {**input_feature_dict, "coordinate": pred_coord},
                "label_full_dict": label_full_dict,
                "max_num_chains": max_num_chains,
                "permute_label": permute_label,
                "dataset_name": input_feature_dict.get("dataset_name", None),
                "pdb_id": input_feature_dict.get("pdb_id", None),
            },
            error_dir=error_dir,
        )
        output_dict, log_dict, permute_pred_indices, permute_label_indices = (
            {},
            {},
            [],
            [],
        )

    return output_dict, log_dict, permute_pred_indices, permute_label_indices
