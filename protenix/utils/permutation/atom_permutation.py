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

import torch

from protenix.metrics.rmsd import rmsd, self_aligned_rmsd
from protenix.model.utils import expand_at_dim, pad_at_dim
from protenix.utils.logger import get_logger
from protenix.utils.permutation.utils import Checker, save_permutation_error

logger = get_logger(__name__)


def run(
    pred_coord: torch.Tensor,
    true_coord: torch.Tensor,
    true_coord_mask: torch.Tensor,
    ref_space_uid: torch.Tensor,
    atom_perm_list: torch.Tensor,
    permute_label: bool = True,
    alignment_mask: torch.Tensor = None,
    error_dir: str = None,
    dataset_name: str = None,
    pdb_id: str = None,
    global_align_wo_symmetric_atom: bool = False,
):
    """apply a permutation to correct symmetric atoms in residues.

    Args:
        Please refer to the args of `correct_symmetric_atoms`.

    Returns:
        if permute_label = True,
            output_dict: a dictionary in the following form, recording the permuted label.
                {
                    "coordinate": permuted_coord,
                    "coordinate_mask": permuted_mask,
                }
            info_dict: a dictionary of logging.
        if permute_label = False,
            output_dict: a dictionary in the following form, recording the permuted prediction.
                {
                    "coordinate": permuted_coord,
                }
            info_dict: a dictionary of logging.
    """
    try:
        permuted_coord, permuted_mask, info_dict, indices_permutation = (
            correct_symmetric_atoms(
                pred_coord=pred_coord,
                true_coord=true_coord,
                true_coord_mask=true_coord_mask,
                ref_space_uid=ref_space_uid,
                atom_perm_list=atom_perm_list,
                permute_label=permute_label,
                alignment_mask=alignment_mask,
                global_align_wo_symmetric_atom=global_align_wo_symmetric_atom,
            )
        )
        if permute_label:
            return (
                {
                    "coordinate": permuted_coord,
                    "coordinate_mask": permuted_mask,
                },
                info_dict,
                indices_permutation,
            )
        else:
            return {"coordinate": permuted_coord}, info_dict, indices_permutation

    except Exception as e:
        error_message = str(e)
        if dataset_name:
            logger.warning(f"dataset: {dataset_name}")
        if pdb_id:
            logger.warning(f"pdb id: {pdb_id}")
        logger.warning(error_message)
        save_permutation_error(
            data={
                "error_message": error_message,
                "pred_coord": pred_coord,
                "true_coord": true_coord,
                "true_coord_mask": true_coord_mask,
                "ref_space_uid": ref_space_uid,
                "atom_perm_list": atom_perm_list,
                "permute_label": permute_label,
                "alignment_mask": alignment_mask,
                "dataset_name": dataset_name,
                "pdb_id": pdb_id,
            },
            error_dir=error_dir,
        )
        return {}, {}, None


def collect_residues_with_symmetric_atoms(
    coord: torch.Tensor,
    coord_mask: torch.Tensor,
    ref_space_uid: torch.Tensor,
    atom_perm_list: list[list],
    run_checker: bool = False,
) -> tuple[list]:
    """Convert atom-level permutation attributes to residue-level attributes.
    Only residues that require symmetric corrections are returned.

    Args:
        coord (torch.Tensor): Coordinates of atoms.
            [N_atom, 3]
        coord_mask (torch.Tensor): The mask indicating whether the atom is resolved in GT.
            [N_atom]
        ref_space_uid (torch.Tensor): Each (chain id, residue index) tuple has a unique ID.
            [N_atom]
        atom_perm_list (list[list]): The atom permutation list, where each sublist contains
                                   the permutation information of the corresponding residue.
            len(atom_perm_list) = N_atom.
            len(atom_perm_list[i]) = N_perm for the residue of atom i.

    """

    device = coord_mask.device

    # Find start & end positions of each residue
    diff = torch.tensor([True] + (ref_space_uid[1:] != ref_space_uid[:-1]).tolist())
    start_positions = torch.cat(
        (torch.nonzero(diff, as_tuple=True)[0], torch.tensor([len(ref_space_uid)]))
    )
    res_start_end = list(
        zip(start_positions[:-1].tolist(), start_positions[1:].tolist())
    )  # [N_res, 2]
    N_res = len(res_start_end)
    assert N_res == len(torch.unique(ref_space_uid))

    position_list = []
    perm_list = []
    coord_list = []
    coord_mask_list = []

    # Traverse residues and store the corresponding data
    for start, end in res_start_end:

        assert len(torch.unique(ref_space_uid[start:end])) == 1

        # Skip if this residue contains < 3 resolved atoms.
        # Alignment requires at least 3 atoms to obtain a reasonable result.
        res_coord_mask = coord_mask[start:end].bool()  # [N_res_atom]
        if res_coord_mask.sum() < 3:
            continue

        # Drop duplicated permutations
        perm = torch.tensor(atom_perm_list[start:end], device=device, dtype=torch.long)
        perm = torch.unique(perm, dim=-1)  # [N_res_atom, N_perm]
        N_res_atom, N_perm = perm.size()

        # Basic checks
        assert perm.min().item() == 0
        assert perm.max().item() == N_res_atom - 1

        # Development checks
        if run_checker:
            Checker.are_permutations(perm, dim=0)
            Checker.contains_identity(perm, dim=0)

        # If all symmetric atoms are unresolved, drop the permutation
        identity_perm = torch.arange(len(perm), device=device).unsqueeze(
            dim=-1
        )  # [N_res_atom, 1]
        is_sym_atom = perm != identity_perm  # [N_res_atom, N_perm]
        is_sym_atom_resolved = is_sym_atom * res_coord_mask.unsqueeze(dim=-1)
        is_valid_perm = is_sym_atom_resolved.any(dim=0)
        if not is_valid_perm.any():
            # Skip if no valid permutation (other than identity) exists
            continue
        perm = perm[..., is_valid_perm]
        perm = torch.cat([identity_perm, perm], dim=-1)  # Put identity to the first
        perm = perm.transpose(-1, -2)  # [N_perm, N_res_atom]

        position_list.append((start, end))
        perm_list.append(perm)
        coord_mask_list.append(res_coord_mask)
        coord_list.append(coord[start:end, :])

    return position_list, coord_list, coord_mask_list, perm_list, N_res


def collect_permuted_coords(
    coord_list: list[torch.Tensor],
    coord_mask_list: list[torch.Tensor],
    perm_list: list[torch.Tensor],
    run_checker: bool = False,
) -> tuple[torch.Tensor]:
    """Apply permutations to coordinates and coordinate masks

    Args:
        coord_list (list[torch.Tensor]): A list of coordinates.
            Each element is a tensor of shape [N_res_atom, 3]. The value N_res_atom can
            vary across different residues.
        coord_mask_list (list[torch.Tensor]): A list of coordinate masks.
            Each element is a tensor of shape [N_res_atom].
        perm_list (list[torch.Tensor]): list of permutations.
            Each element is a long tensor of shape [N_perm, N_res_atom]. The value N_perm
            can vary across different residues.

    Returns:
        torch.Tensor:
            [N_total_perm, MAX_N_res_atom, 3]
            [N_total_perm, MAX_N_res_atom]
    """

    MAX_N_res_atom = max(perm.size(-1) for perm in perm_list)
    perm_coord = []  # [N_total_perm, N_res_atom, 3]
    perm_coord_mask = []  # [N_total_perm, N_res_atom]

    N_total_perm = 0
    for perm, res_coord, res_coord_mask in zip(perm_list, coord_list, coord_mask_list):

        # Basic shape checks
        N_perm, N_res_atom = perm.size()
        assert res_coord.size(-1) == 3
        assert res_coord.size(0) == res_coord_mask.size(0) == perm.size(-1)

        # Permute coordinates & masks
        res_coord_permuted = res_coord[perm]  # [N_perm, N_res_atom, 3]
        res_coord_mask_permuted = res_coord_mask[perm]  # [N_perm, N_res_atom]
        assert res_coord_permuted.size() == (N_perm, N_res_atom, 3)
        assert res_coord_mask_permuted.size() == (N_perm, N_res_atom)

        if run_checker:
            Checker.are_permutations(perm, dim=-1)
            Checker.batch_permute(perm, res_coord, res_coord_permuted)
            Checker.batch_permute(perm, res_coord_mask, res_coord_mask_permuted)

        # Pad to MAX_N_res_atom
        N_res_atom = perm.size(dim=-1)
        if N_res_atom < MAX_N_res_atom:
            pad_length = (0, MAX_N_res_atom - N_res_atom)
            res_coord_permuted = pad_at_dim(
                res_coord_permuted, dim=-2, pad_length=pad_length
            )  # [N_perm, MAX_N_res_atom, 3]
            res_coord_mask_permuted = pad_at_dim(
                res_coord_mask_permuted, dim=-1, pad_length=pad_length
            )

        N_total_perm += N_perm
        perm_coord.append(res_coord_permuted)
        perm_coord_mask.append(res_coord_mask_permuted)

    perm_coord = torch.cat(perm_coord, dim=0)
    perm_coord_mask = torch.cat(perm_coord_mask, dim=0)

    # Shape check
    assert perm_coord.size() == (N_total_perm, MAX_N_res_atom, 3)
    assert perm_coord_mask.size() == (N_total_perm, MAX_N_res_atom)

    return perm_coord, perm_coord_mask


class AtomPermutation(object):
    def __init__(
        self,
        eps: float = 1e-8,
        run_checker: bool = False,
        global_align_wo_symmetric_atom: bool = False,
    ):
        """Class for assigning the optimal permutations of true coordinates/pred coordinates and coordinate masks.

        Args:
            eps (float): A small number used in alignment.
            run_checker (bool): If true, it applies more checkers to ensure the correctness.
            global_align_wo_symmetric_atom (bool):  If true, the global alignment before AtomPermutation will not consider atoms with permutation.
        """

        self.eps = eps
        self.run_checker = run_checker
        self.global_align_wo_symmetric_atom = global_align_wo_symmetric_atom

    @staticmethod
    def check_input_shape(
        pred_coord: torch.Tensor,
        true_coord: torch.Tensor,
        true_coord_mask: torch.Tensor,
        ref_space_uid: torch.Tensor,
        atom_perm_list: list[list],
    ):

        N_atom = len(true_coord)
        assert true_coord.dim() == 2
        assert true_coord_mask.dim() == 1
        assert ref_space_uid.dim() == 1

        assert true_coord.size(-1) == 3
        assert true_coord.size(-2) == N_atom
        assert ref_space_uid.size(-1) == N_atom
        assert len(atom_perm_list) == N_atom

        assert pred_coord.dim() in [2, 3]  # for simplicity
        assert pred_coord.shape[-2:] == (N_atom, 3)

    @staticmethod
    def global_align_pred_to_true(
        pred_coord: torch.Tensor,
        true_coord: torch.Tensor,
        true_coord_mask: torch.Tensor,
        eps: float = 1e-8,
    ) -> tuple[torch.Tensor]:
        """Align the predicted coordinates to true coordinates

        Args:
            pred_coord (torch.Tensor):
                [Batch, N_atom, 3] or [N_atom, 3]
            true_coord (torch.Tensor):
                [N_atom, 3]
            true_coord_mask (torch.Tensor):
                [N_atom, 3]

        Returns:
            aligned_rmsd (torch.Tensor):
                [Batch] or []
            transformed_pred_coord (torch.Tensor): having the same shape as pred_coord.
                [Batch, N_atom, 3] or [N_atom, 3]
        """

        if true_coord.dim() < pred_coord.dim():
            assert pred_coord.dim() == 3  # [Batch, N_atom, 3]
            Batch = pred_coord.size(0)
            expand_func = lambda x: expand_at_dim(x, dim=0, n=Batch)
        else:
            expand_func = lambda x: x

        with torch.cuda.amp.autocast(enabled=False):
            aligned_rmsd, transformed_pred_coord, _, _ = self_aligned_rmsd(
                pred_pose=pred_coord.to(torch.float32),
                true_pose=expand_func(true_coord.to(torch.float32)),
                atom_mask=expand_func(true_coord_mask),
                allowing_reflection=False,
                reduce=False,
                eps=eps,
            )  # [Batch], [Batch, N_atom, 3]

        return aligned_rmsd, transformed_pred_coord

    @staticmethod
    def get_identity_permutation(batch_shape, N_atom, device):
        """Return identity permutation indices if no multiple-permutation exists for every residue

        Returns:
            torch.Tensor: identity permutation of indices
                [N_atom] or [Batch, N_atom]
        """
        identity = torch.arange(N_atom, device=device)
        if len(batch_shape) == 0:
            return identity
        else:
            assert len(batch_shape) == 1
            return torch.stack([identity for _ in range(batch_shape[0])], dim=0)

    @staticmethod
    def _optimize_per_residue_permutation_by_rmsd(
        per_residue_pred_coord_list: list[torch.Tensor],
        per_residue_coord_list: list[torch.Tensor],
        per_residue_coord_mask_list: list[torch.Tensor],
        per_residue_perm_list: list[torch.Tensor],
        eps: float = 1e-8,
        run_checker: bool = False,
    ) -> tuple[list[torch.Tensor]]:
        """Find the optimal permutations of true coordinates and coordinate masks to minimize the
        RMSD between true coordinates and predicted coordinates.

        Args:
            per_residue_pred_coord_list (torch.Tensor): List of residues. Each element records
                the predicted atom coordinates of one residue. Each element has shape
                [N_res_atom, 3] or [Batch, N_res_atom, 3]
            per_residue_coord_list (list[torch.Tensor]): List of residues. Each element records
                the atom coordinates of one residue. Each element has shape [N_res_atom, 3].
            per_residue_coord_mask_list (list[torch.Tensor]): List of residues. Each element records
                the atom coordinate masks of one residue. Each element has shape [N_res_atom].
            per_residue_perm_list (list[torch.Tensor]): List of residues. Each element records
                the atom permutations of one residue. Each element has shape [N_perm, N_res_atom].
            eps (float, optional): A small number, used in alignment. Defaults to 1e-8.
            run_checker (bool, optional): If True, run extensive checks.

        Returns:
            best_permutation_list (list[torch.Tensor]): List of residues. Each element records the
                optimal permutation that should apply to true coordinates for one residue. Each element
                has shape
                [N_res_atom] or [Batch, N_res_atom]
            is_permuted_list (list[torch.Tensor]): List of residues. Each element records whether the
                atoms in this residue is permuted. Each element has shape
                [] or [Batch]
            optimized_rmsd_list (list[torch.Tensor]): List of residues. Each element records the optimized
                rmsd of the residue. Each element has shape
                [] or [Batch]
            original_rmsd_list (list[torch.Tensor]): List of residues. Each element records the original
                rmsd of the residue. Each element has shape
                [] or [Batch]
        """

        # Find max number of per-residue atoms
        per_residue_N_perm = [perm.size(0) for perm in per_residue_perm_list]
        per_residue_N_atom = [perm.size(1) for perm in per_residue_perm_list]
        N_max_atom = max(per_residue_N_atom)

        # Permute true coordinates & masks according to the permutations in per_residue_perm_list
        permuted_coord, permuted_coord_mask = collect_permuted_coords(
            coord_list=per_residue_coord_list,
            coord_mask_list=per_residue_coord_mask_list,
            perm_list=per_residue_perm_list,
            run_checker=run_checker,
        )  # [N_total_perm, N_max_atom, 3], [N_total_perm, N_max_atom]
        assert permuted_coord.size(-2) == permuted_coord_mask.size(-1) == N_max_atom
        N_total_perm = permuted_coord.size(0)

        # Pad 'pred_coord' to the same shape as 'permuted_coord'
        per_residue_pred_coord_list = [
            pad_at_dim(
                p_coord, dim=-2, pad_length=(0, N_max_atom - p_coord.size(dim=-2))
            )
            for p_coord in per_residue_pred_coord_list
        ]
        # Repeat N_perm times for each residue
        pred_coord = torch.stack(
            sum(
                [
                    [p_coord] * N_perm
                    for N_perm, p_coord in zip(
                        per_residue_N_perm, per_residue_pred_coord_list
                    )
                ],
                [],
            ),
            dim=-3,
        )  # [N_total_perm, N_max_atom, 3] or [Batch, N_total_perm, N_max_atom, 3]
        assert pred_coord.shape[-3:] == (N_total_perm, N_max_atom, 3)

        batch_shape = pred_coord.shape[:-3]
        assert len(batch_shape) in [0, 1]
        if len(batch_shape) == 1:
            # expand true coord & mask to have the same batch size as pred coord
            Batch = pred_coord.size(0)
            permuted_coord = expand_at_dim(permuted_coord, dim=0, n=Batch)
            permuted_coord_mask = expand_at_dim(permuted_coord_mask, dim=0, n=Batch)

        # Compute per-residue rmsd
        with torch.cuda.amp.autocast(enabled=False):
            per_res_rmsd = rmsd(
                pred_pose=pred_coord.to(torch.float32),
                true_pose=permuted_coord.to(torch.float32),
                mask=permuted_coord_mask,
                eps=eps,
                reduce=False,
            )  # [N_total_perm] or [Batch, N_total_perm]
            assert per_res_rmsd.size() == batch_shape + (N_total_perm,)

        # Find the best permutation
        best_permutation_list = []
        is_permuted_list = []
        original_rmsd_list = []
        optimized_rmsd_list = []
        i = 0

        # Enumerate over all residues (could be improved by scatter)
        for N_perm, N_res_atom, perm in zip(
            per_residue_N_perm, per_residue_N_atom, per_residue_perm_list
        ):
            cur_res_rmsd = per_res_rmsd[..., i : i + N_perm]  # [batch_shape, N_perm]
            best_rmsd, best_j = torch.min(cur_res_rmsd, dim=-1)  # [batch_shape]
            best_perm = perm[best_j]  # [batch_shape, N_res_atom]
            best_permutation_list.append(best_perm)

            is_permuted_list.append(
                best_j > 0
            )  # The first of the perm lists is the identity

            optimized_rmsd_list.append(best_rmsd)
            original_rmsd_list.append(cur_res_rmsd[..., 0])

            i += N_perm

            if run_checker:
                assert perm.size() == (N_perm, N_res_atom)
                assert cur_res_rmsd.size() == batch_shape + (N_perm,)
                assert best_rmsd.size() == batch_shape
                assert best_j.size() == batch_shape
                assert best_perm.size() == batch_shape + (N_res_atom,)
                Checker.are_permutations(best_perm, dim=-1)

                def _check_identity(j_value, perm):
                    if j_value > 0:
                        Checker.not_contain_identity(perm)
                    else:
                        Checker.contains_identity(perm)

                if best_j.dim() == 0:
                    _check_identity(best_j, best_perm)
                else:
                    for j_value, perm_j in zip(best_j, best_perm):
                        _check_identity(j_value, perm_j)

        return (
            best_permutation_list,
            is_permuted_list,
            optimized_rmsd_list,
            original_rmsd_list,
        )

    def __call__(
        self,
        pred_coord: torch.Tensor,
        true_coord: torch.Tensor,
        true_coord_mask: torch.Tensor,
        ref_space_uid: torch.Tensor,
        atom_perm_list: list[list],
        alignment_mask: torch.Tensor,
        verbose: bool = False,
        run_checker: bool = False,
    ):
        """

        Args:
            pred_coord (torch.Tensor): Predicted coordinates of atoms.
                [N_atom, 3] or [Batch, atom, 3]
            true_coord (torch.Tensor): true coordinates of atoms.
                [N_atom, 3]
            true_coord_mask (torch.Tensor): The mask indicating whether the atom is resolved.
                [N_atom]
            ref_space_uid (torch.Tensor): Each (chain id, residue index) tuple has a unique ID.
                [N_atom]
            atom_perm_list (list[list]): The atom permutation list, where each sublist contains
                                          the permutation information of the corresponding residue.
                len(atom_perm_list) = N_atom.
                len(atom_perm_list[i]) = N_perm for the residue of atom i.
            permute_label (bool, optional): If true, return indices permutations of the true coordinate.
                Otherwise, return indices permutations for the predicted coordinate. Defaults to True.
            alignment_mask (torch.Tensor, optional): Defaults to None. A mask indicating which atoms to
                consider while performing the alignment.
            verbose (bool, optional): Defaults to False.
            run_checker (bool, optional): Whether running more checks for debug. Defaults to False.

        Returns:
            permutation (torch.Tensor): the optimized permutation of atoms.
                [N_atom] or [Batch, N_atom]
            log_dict (Dict): a dictionary recording the permutation stats.
        """

        # Basic Info & Shape checker
        device = pred_coord.device
        batch_shape = pred_coord.shape[:-2]
        N_atom = pred_coord.size(-2)
        self.check_input_shape(
            pred_coord, true_coord, true_coord_mask, ref_space_uid, atom_perm_list
        )

        # Initialize log dict
        log_dict = {}

        # Initialize the permutation as identity
        permutation = self.get_identity_permutation(
            batch_shape, N_atom=N_atom, device=device
        )

        # Collect residues that require permutations
        (
            per_residue_position_list,
            per_residue_coord_list,
            per_residue_coord_mask_list,
            per_residue_perm_list,
            N_res,
        ) = collect_residues_with_symmetric_atoms(
            coord=true_coord,
            coord_mask=true_coord_mask,
            ref_space_uid=ref_space_uid,
            atom_perm_list=atom_perm_list,
            run_checker=run_checker,
        )
        log_dict["N_res"] = N_res
        log_dict["N_res_with_symmetry"] = len(per_residue_coord_list)
        log_dict["N_res_permuted"] = 0.0
        log_dict["has_res_permuted"] = 0

        # If no residues contain symmetry, return now.
        if not per_residue_perm_list:
            print("No atom permutation is needed. Return the identity permutation.")
            return (permutation, log_dict)

        # no_permute_atom_mask: 1 represent this atom can not be permuted
        no_permute_atom_mask = torch.ones_like(true_coord_mask)
        for (start, end), per_residue_perm in zip(
            per_residue_position_list, per_residue_perm_list
        ):
            no_permute_atom_mask[start:end] = 1 - (
                (per_residue_perm != per_residue_perm[0]).sum(dim=0) > 0
            ).to(torch.int32)

        # Perform a global alignment of predictions to true coordinates
        if alignment_mask is None:
            alignment_mask = true_coord_mask
        else:
            alignment_mask = true_coord_mask * alignment_mask.bool()
        if self.global_align_wo_symmetric_atom:
            alignment_mask = no_permute_atom_mask * alignment_mask

        if alignment_mask.sum().item() < 3:
            print("No atom permutation is needed. Return the identity permutation.")
            return (permutation, log_dict)

        # This is for atom permutation, use mask with different strategies
        _, transformed_pred_coord = self.global_align_pred_to_true(
            pred_coord,
            true_coord,
            alignment_mask,
            eps=self.eps,
        )
        # This is for unpermuted all-atom baseline calculation
        aligned_rmsd, _ = self.global_align_pred_to_true(
            pred_coord,
            true_coord,
            true_coord_mask,
            eps=self.eps,
        )
        log_dict["unpermuted_rmsd"] = aligned_rmsd.mean().item()  # [Batch]

        """ 
        To efficiently optimize the residues parallely, group the residues
        according to the number of atoms in each residue.
        """
        per_residue_N_atom = [coord.size(0) for coord in per_residue_coord_list]
        res_atom_cutoff = [15, 30, 50, 100, 100000]
        grouped_indices = {}
        for i, n in enumerate(per_residue_N_atom):
            for atom_cutoff in res_atom_cutoff:
                if n <= atom_cutoff:
                    break
            grouped_indices.setdefault(atom_cutoff, []).append(i)

        assert len(sum(list(grouped_indices.values()), [])) == len(
            per_residue_perm_list
        )

        residue_position_list = []
        residue_best_permutation_list = []
        residue_is_permuted_list = []
        residue_optimized_rmsd_list = []
        residue_original_rmsd_list = []
        for atom_cutoff, residue_group in grouped_indices.items():

            if verbose:
                print(f"{len(residue_group)} residues have <={atom_cutoff} atoms.")

            # Enumerte permutations within each residue to minimize per-residue RMSD
            per_res_pos_list = [per_residue_position_list[i] for i in residue_group]
            (
                per_res_best_permutation,
                per_res_is_permuted,
                per_res_optimized_rmsd,
                per_res_ori_rmsd,
            ) = self._optimize_per_residue_permutation_by_rmsd(
                per_residue_pred_coord_list=[
                    transformed_pred_coord[..., pos[0] : pos[1], :]
                    for pos in per_res_pos_list
                ],
                per_residue_coord_list=[
                    per_residue_coord_list[i] for i in residue_group
                ],
                per_residue_coord_mask_list=[
                    per_residue_coord_mask_list[i] for i in residue_group
                ],
                per_residue_perm_list=[per_residue_perm_list[i] for i in residue_group],
                eps=self.eps,
                run_checker=self.run_checker,
            )
            residue_position_list.extend(per_res_pos_list)
            residue_best_permutation_list.extend(per_res_best_permutation)
            residue_is_permuted_list.extend(per_res_is_permuted)
            residue_optimized_rmsd_list.extend(per_res_optimized_rmsd)
            residue_original_rmsd_list.extend(per_res_ori_rmsd)

        # Aggregate per_residue results
        # 1. Best permutation
        indices_list = [
            torch.arange(pos[0], pos[1], device=device) for pos in residue_position_list
        ]
        residue_atom_indices = torch.cat(indices_list, dim=-1)  # [N_perm_atom]
        residue_best_permutation = torch.cat(
            [
                ind[perm]
                for ind, perm in zip(indices_list, residue_best_permutation_list)
            ],
            dim=-1,
        )  # [Batch, N_perm_atom] or [N_perm_atom]
        permutation[..., residue_atom_indices] = residue_best_permutation

        # 2. Other statistics
        is_res_permuted = torch.stack(residue_is_permuted_list, dim=-1).float()
        log_dict["N_res_permuted"] = is_res_permuted.sum(dim=-1).mean().item()
        log_dict["has_res_permuted"] = (
            (is_res_permuted.sum(dim=-1) > 0).float().mean().item()
        )

        return permutation, log_dict


def correct_symmetric_atoms(
    pred_coord: torch.Tensor,
    true_coord: torch.Tensor,
    true_coord_mask: torch.Tensor,
    ref_space_uid: torch.Tensor,
    atom_perm_list: list[list],
    permute_label: bool = True,
    alignment_mask: torch.Tensor = None,
    verbose: bool = False,
    run_checker: bool = False,
    eps: float = 1e-8,
    global_align_wo_symmetric_atom: bool = False,
):
    """
    Return optimally permuted true coordinates and masks according to the predicted coordinates
    Or, return optimalled permuted predicted coordinates if permute_label is False.

    Args:
        pred_coord (torch.Tensor): predicted atom positions
            [Batch, N_atom, 3] or [N_atom, 3]
        true_coord (torch.Tensor): true atom positions
        true_coord_mask (torch.Tensor): a mask indicating whether the atom is resolved.
        ref_space_uid (torch.Tensor): unique residue ID for each atom.
            [N_atom]
        atom_perm_list (list[list]): The atom permutation list, where each sublist contains
                                     the permutation information of the corresponding residue.
            len(atom_perm_list) = N_atom.
            len(atom_perm_list[i]) = N_perm for the residue of atom i.
        permute_label (bool): indicates whether permuted true coordinates are returned or
            predicted coordinates are returned.
        alignment_mask (torch.Tensor, optional): a mask indicating which atoms are considered while
            performing the alignment.
            [N_atom]
        eps (float, optional): A small number used in alignment. Defaults to 1e-8.
        global_align_wo_symmetric_atom (bool):  If true, the global alignment before AtomPermutation will not consider atoms has permutation.

    Returns:
        If permute_label is True, it returns
            coordinate (torch.Tensor): permuted true coordinates.
                [Batch, N_atom, 3] or [N_atom, 3]
            coordinate_mask (torch.Tensor): permuted true coordinate masks.
                [Batch, N_atom] or [N_atom]
        If permuted_label is False, it returns the permuted prediction.
            [Batch, N_atom, 3] or [N_atom, 3]

        log_dict: logging info for the permutation
            percent_res_permuted (torch.Tensor): percentage of residues (excluding those with less than 3 atoms or identity perm only) that have been permuted
            best_aligned_rmsd_improved: rmsd improved after permutation, using self_aligned_rmsd
    """

    assert pred_coord.dim() in [2, 3]
    assert pred_coord.size(-1) == 3

    if alignment_mask is not None:
        alignment_mask = (true_coord_mask * alignment_mask).bool()
    else:
        alignment_mask = true_coord_mask.bool()

    with torch.no_grad():
        # Do not compute gradient while optimizing the permutation
        atom_perm = AtomPermutation(
            run_checker=run_checker,
            eps=eps,
            global_align_wo_symmetric_atom=global_align_wo_symmetric_atom,
        )
        indices_permutation, log_dict = atom_perm(
            pred_coord,
            true_coord,
            true_coord_mask,
            ref_space_uid,
            atom_perm_list,
            alignment_mask=alignment_mask,
            verbose=verbose,
        )

    # Log aligned rmsd after permutation
    if "unpermuted_rmsd" in log_dict:
        # This is the final permuted all-atom rmsd
        permuted_rmsd, _ = AtomPermutation.global_align_pred_to_true(
            pred_coord,
            true_coord[indices_permutation],
            true_coord_mask[indices_permutation],
            eps=eps,
        )
        log_dict["permuted_rmsd"] = permuted_rmsd.mean().item()
        log_dict["improved_rmsd"] = (
            log_dict["unpermuted_rmsd"] - log_dict["permuted_rmsd"]
        )

    if permute_label:
        return (
            true_coord[indices_permutation],
            true_coord_mask[indices_permutation],
            log_dict,
            indices_permutation,
        )
    else:
        # Find the permutation of the prediction
        if pred_coord.dim() == 2:
            # Inverse permutation for 1D case
            indices_permutation = torch.argsort(indices_permutation)
            pred_coord_permuted = pred_coord[indices_permutation]
        else:
            # Inverse permutation for 2D case (batch mode)
            indices_permutation = torch.argsort(indices_permutation, dim=1)
            indices_permutation_expanded = expand_at_dim(
                indices_permutation, dim=-1, n=3
            )  # [Batch, N_atom, 3]
            pred_coord_permuted = pred_coord.gather(1, indices_permutation_expanded)

        return pred_coord_permuted, None, log_dict, indices_permutation
