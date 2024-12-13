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

from typing import Optional

import torch
import torch.nn as nn

from protenix.model import sample_confidence


def get_complex_level_rankers(scores, keys):
    assert all([k in ["plddt", "gpde", "ranking_score"] for k in keys])
    rankers = {}
    for key in keys:
        if key == "gpde":
            descending = False
        else:
            descending = True
        ranking = scores[key].argsort(dim=0, descending=descending)
        rankers[f"{key}.rank1"] = lambda x, rank1_idx=ranking[0].item(): x[
            ..., rank1_idx
        ]
    return rankers


def add_diff_metrics(scores, ranker_keys):
    diff_metrics = {
        "diff/best_worst": scores["best"] - scores["worst"],
        "diff/best_random": scores["best"] - scores["random"],
        "diff/best_median": scores["best"] - scores["median"],
    }

    for key in ranker_keys:
        diff_metrics.update(
            {
                f"diff/best_{key}": scores["best"] - scores[f"{key}.rank1"],
                f"diff/{key}_median": scores[f"{key}.rank1"] - scores["median"],
            }
        )
    scores.update(diff_metrics)
    return scores


class LDDTMetrics(nn.Module):
    """LDDT: evaluated on chains and interfaces"""

    def __init__(self, configs):
        super(LDDTMetrics, self).__init__()
        self.eps = configs.metrics.lddt.eps
        self.configs = configs
        self.chunk_size = self.configs.infer_setting.lddt_metrics_chunk_size
        self.lddt_base = LDDT(eps=self.eps)

        self.complex_ranker_keys = configs.metrics.get(
            "complex_ranker_keys", ["plddt", "gpde", "ranking_score"]
        )

    def compute_lddt(self, pred_dict: dict, label_dict: dict):
        """compute complex-level and chain/interface-level lddt

        Args:
            pred_dict (Dict): a dictionary containing
                coordinate: [N_sample, N_atom, 3]
            label_dict (Dict): a dictionary containing
                coordinate: [N_sample, N_atom, 3]
                lddt_mask: [N_atom, N_atom]
        """

        out = {}

        # Complex-level
        lddt = self.lddt_base.forward(
            pred_coordinate=pred_dict["coordinate"],
            true_coordinate=label_dict["coordinate"],
            lddt_mask=label_dict["lddt_mask"],
            chunk_size=self.chunk_size,
        )  # [N_sample]
        out["complex"] = lddt

        return out

    def aggregate(
        self,
        vals,
        dim: int = -1,
        aggregators: dict = {},
    ):
        N_sample = vals.size(dim)
        median_index = N_sample // 2
        basic_sample_aggregators = {
            "best": lambda x: x.max(dim=dim)[0],
            "worst": lambda x: x.min(dim=dim)[0],
            "random": lambda x: x.select(dim=dim, index=0),
            "mean": lambda x: x.mean(dim=dim),
            "median": lambda x: x.sort(dim=dim, descending=True)[0].select(
                dim=dim, index=median_index
            ),
        }
        sample_aggregators = {**basic_sample_aggregators, **aggregators}

        return {
            agg_name: agg_func(vals)
            for agg_name, agg_func in sample_aggregators.items()
        }

    def aggregate_lddt(self, lddt_dict, per_sample_summary_confidence):

        # Merge summary_confidence results
        confidence_scores = sample_confidence.merge_per_sample_confidence_scores(
            per_sample_summary_confidence
        )

        # Complex-level LDDT
        complex_level_ranker = get_complex_level_rankers(
            confidence_scores, self.complex_ranker_keys
        )

        complex_lddt = self.aggregate(
            lddt_dict["complex"], aggregators=complex_level_ranker
        )
        complex_lddt = add_diff_metrics(complex_lddt, self.complex_ranker_keys)
        # Log metrics
        complex_lddt = {
            f"lddt/complex/{name}": value for name, value in complex_lddt.items()
        }
        return complex_lddt, {}


class LDDT(nn.Module):
    """LDDT base metrics"""

    def __init__(self, eps: float = 1e-10):
        super(LDDT, self).__init__()
        self.eps = eps

    def _chunk_base_forward(self, pred_distance, true_distance) -> torch.Tensor:
        distance_error_l1 = torch.abs(
            pred_distance - true_distance
        )  # [N_sample, N_pair_sparse]
        thresholds = [0.5, 1, 2, 4]
        sparse_pair_lddt = (
            torch.stack([distance_error_l1 < t for t in thresholds], dim=-1)
            .to(dtype=distance_error_l1.dtype)
            .mean(dim=-1)
        )  # [N_sample, N_pair_sparse]
        del distance_error_l1
        # Compute mean
        if sparse_pair_lddt.numel() == 0:  # corespand to all zero in dense mask
            sparse_pair_lddt = torch.zeros_like(sparse_pair_lddt)
        lddt = torch.mean(sparse_pair_lddt, dim=-1)
        return lddt

    def _chunk_forward(
        self, pred_distance, true_distance, chunk_size: Optional[int] = None
    ) -> torch.Tensor:
        if chunk_size is None:
            return self._chunk_base_forward(pred_distance, true_distance)
        else:
            lddt = []
            N_sample = pred_distance.shape[-2]
            no_chunks = N_sample // chunk_size + (N_sample % chunk_size != 0)
            for i in range(no_chunks):
                lddt_i = self._chunk_base_forward(
                    pred_distance[
                        ...,
                        i * chunk_size : (i + 1) * chunk_size,
                        :,
                    ],
                    true_distance,
                )
                lddt.append(lddt_i)
            lddt = torch.cat(lddt, dim=-1)  # [N_sample]
            return lddt

    def _calc_sparse_dist(self, pred_coordinate, true_coordinate, l_index, m_index):
        pred_coords_l = pred_coordinate.index_select(
            -2, l_index
        )  # [N_sample, N_atom_sparse_l, 3]
        pred_coords_m = pred_coordinate.index_select(
            -2, m_index
        )  # [N_sample, N_atom_sparse_m, 3]
        true_coords_l = true_coordinate.index_select(
            -2, l_index
        )  # [N_atom_sparse_l, 3]
        true_coords_m = true_coordinate.index_select(
            -2, m_index
        )  # [N_atom_sparse_m, 3]

        pred_distance_sparse_lm = torch.norm(
            pred_coords_l - pred_coords_m, p=2, dim=-1
        )  # [N_sample, N_pair_sparse]
        true_distance_sparse_lm = torch.norm(
            true_coords_l - true_coords_m, p=2, dim=-1
        )  # [N_sample, N_pair_sparse]
        return pred_distance_sparse_lm, true_distance_sparse_lm

    def forward(
        self,
        pred_coordinate: torch.Tensor,
        true_coordinate: torch.Tensor,
        lddt_mask: torch.Tensor,
        chunk_size: Optional[int] = None,
    ) -> dict[str, torch.Tensor]:
        """LDDT: evaluated on complex, chains and interfaces
        sparse implementation, which largely reduce cuda memory when atom num reaches 10^4 +

        Args:
            pred_coordinate (torch.Tensor): the pred coordinates
                [N_sample, N_atom, 3]
            true_coordinate (torch.Tensor): the ground truth atom coordinates
                [N_atom, 3]
            lddt_mask (torch.Tensor):
                sparse version of [N_atom, N_atom] atompair mask based on bespoke radius of true distance
                [N_nonzero_mask, 2]

        Returns:
            Dict[str, torch.Tensor]:
                "best": [N_eval]
                "worst": [N_eval]
        """
        lddt_indices = torch.nonzero(lddt_mask, as_tuple=True)
        l_index = lddt_indices[0]
        m_index = lddt_indices[1]
        pred_distance_sparse_lm, true_distance_sparse_lm = self._calc_sparse_dist(
            pred_coordinate, true_coordinate, l_index, m_index
        )
        group_lddt = self._chunk_forward(
            pred_distance_sparse_lm, true_distance_sparse_lm, chunk_size=chunk_size
        )  # [N_sample]
        return group_lddt

    @staticmethod
    def compute_lddt_mask(
        true_coordinate: torch.Tensor,
        true_coordinate_mask: torch.Tensor,
        is_nucleotide: torch.Tensor = None,
        is_nucleotide_threshold: float = 30.0,
        threshold: float = 15.0,
    ):
        # Distance mask
        distance_mask = (
            true_coordinate_mask[..., None] * true_coordinate_mask[..., None, :]
        )
        # Distances for all atom pairs
        # Note: we convert to bf16 for saving cuda memory, if performance drops, do not convert it
        distance = (torch.cdist(true_coordinate, true_coordinate) * distance_mask).to(
            true_coordinate.dtype
        )  # [..., N_atom, N_atom]

        # Local mask
        c_lm = distance < threshold  # [..., N_atom, N_atom]
        if is_nucleotide is not None:
            # Use a different radius for nucleotide
            is_nucleotide_mask = is_nucleotide.bool()[..., None]
            c_lm = (distance < is_nucleotide_threshold) * is_nucleotide_mask + c_lm * (
                ~is_nucleotide_mask
            )

        # Zero-out diagonals of c_lm and cast to float
        c_lm = c_lm * (
            1 - torch.eye(n=c_lm.size(-1), device=c_lm.device, dtype=distance.dtype)
        )
        # Zero-out atom pairs without true coordinates
        c_lm = c_lm * distance_mask  # [..., N_atom, N_atom]
        return c_lm
