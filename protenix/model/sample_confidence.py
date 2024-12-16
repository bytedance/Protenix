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

from typing import Optional, Union

import torch
from ml_collections.config_dict import ConfigDict

from protenix.metrics.clash import Clash
from protenix.utils.distributed import traverse_and_aggregate


def merge_per_sample_confidence_scores(summary_confidence_list: list[dict]) -> dict:
    """
    Merge confidence scores from multiple samples into a single dictionary.

    Args:
        summary_confidence_list (list[dict]): List of dictionaries containing confidence scores for each sample.

    Returns:
        dict: Merged dictionary of confidence scores.
    """

    def stack_score(tensor_list: list):
        if tensor_list[0].dim() == 0:
            tensor_list = [x.unsqueeze(0) for x in tensor_list]
        score = torch.stack(tensor_list, dim=0)
        return score

    return traverse_and_aggregate(summary_confidence_list, aggregation_func=stack_score)


def _compute_full_data_and_summary(
    configs: ConfigDict,
    pae_logits: torch.Tensor,
    plddt_logits: torch.Tensor,
    pde_logits: torch.Tensor,
    contact_probs: torch.Tensor,
    token_asym_id: torch.Tensor,
    token_has_frame: torch.Tensor,
    atom_coordinate: torch.Tensor,
    atom_to_token_idx: torch.Tensor,
    atom_is_polymer: torch.Tensor,
    N_recycle: int,
    interested_atom_mask: Optional[torch.Tensor] = None,
    elements_one_hot: Optional[torch.Tensor] = None,
    mol_id: Optional[torch.Tensor] = None,
    return_full_data: bool = False,
) -> tuple[list[dict], list[dict]]:
    """
    Compute full data and summary confidence scores for the given inputs.

    Args:
        configs: Configuration object.
        pae_logits (torch.Tensor): Logits for PAE (Predicted Aligned Error).
        plddt_logits (torch.Tensor): Logits for pLDDT (Predicted Local Distance Difference Test).
        pde_logits (torch.Tensor): Logits for PDE (Predicted Distance Error).
        contact_probs (torch.Tensor): Contact probabilities.
        token_asym_id (torch.Tensor): Asymmetric ID for tokens.
        token_has_frame (torch.Tensor): Indicator for tokens having a frame.
        atom_coordinate (torch.Tensor): Atom coordinates.
        atom_to_token_idx (torch.Tensor): Mapping from atoms to tokens.
        atom_is_polymer (torch.Tensor): Indicator for atoms being part of a polymer.
        N_recycle (int): Number of recycles.
        interested_atom_mask (Optional[torch.Tensor]): Mask for interested atoms. Defaults to None.
        elements_one_hot (Optional[torch.Tensor]): One-hot encoding for elements. Defaults to None.
        mol_id (Optional[torch.Tensor]): Molecular ID. Defaults to None.
        return_full_data (bool): Whether to return full data. Defaults to False.

    Returns:
        tuple[list[dict], list[dict]]:
            - summary_confidence: List of dictionaries containing summary confidence scores.
            - full_data: List of dictionaries containing full data if `return_full_data` is True.
    """
    atom_is_ligand = (1 - atom_is_polymer).long()
    token_is_ligand = torch.zeros_like(token_asym_id).scatter_add(
        0, atom_to_token_idx, atom_is_ligand
    )
    token_is_ligand = token_is_ligand > 0

    full_data = {}
    full_data["atom_plddt"] = logits_to_score(
        plddt_logits, **get_bin_params(configs.loss.plddt)
    )  # [N_s, N_atom]
    # Cpu offload for saving cuda memory
    pde_logits = pde_logits.to(plddt_logits.device)
    full_data["token_pair_pde"] = logits_to_score(
        pde_logits, **get_bin_params(configs.loss.pde)
    )  # [N_s, N_token, N_token]
    del pde_logits
    full_data["contact_probs"] = contact_probs.clone()  # [N_token, N_token]
    pae_logits = pae_logits.to(plddt_logits.device)
    full_data["token_pair_pae"], pae_prob = logits_to_score(
        pae_logits, **get_bin_params(configs.loss.pae), return_prob=True
    )  # [N_s, N_token, N_token]
    del pae_logits

    summary_confidence = {}
    summary_confidence["plddt"] = full_data["atom_plddt"].mean(dim=-1) * 100  # [N_s, ]
    summary_confidence["gpde"] = (
        full_data["token_pair_pde"] * full_data["contact_probs"]
    ).sum(dim=[-1, -2]) / full_data["contact_probs"].sum(dim=[-1, -2])

    summary_confidence["ptm"] = calculate_ptm(
        pae_prob, has_frame=token_has_frame, **get_bin_params(configs.loss.pae)
    )  # [N_s, ]
    summary_confidence["iptm"] = calculate_iptm(
        pae_prob,
        has_frame=token_has_frame,
        asym_id=token_asym_id,
        **get_bin_params(configs.loss.pae)
    )  # [N_s, ]

    # Add: 'chain_pair_iptm', 'chain_pair_iptm_global' 'chain_iptm', 'chain_ptm'
    summary_confidence.update(
        calculate_chain_based_ptm(
            pae_prob,
            has_frame=token_has_frame,
            asym_id=token_asym_id,
            token_is_ligand=token_is_ligand,
            **get_bin_params(configs.loss.pae)
        )
    )
    # Add: 'chain_plddt', 'chain_pair_plddt'
    summary_confidence.update(
        calculate_chain_based_plddt(
            full_data["atom_plddt"], token_asym_id, atom_to_token_idx
        )
    )
    del pae_prob
    summary_confidence["has_clash"] = calculate_clash(
        atom_coordinate,
        token_asym_id,
        atom_to_token_idx,
        atom_is_polymer,
        configs.metrics.clash.af3_clash_threshold,
    )
    summary_confidence["num_recycles"] = torch.tensor(
        N_recycle, device=atom_coordinate.device
    )
    # TODO: disorder
    summary_confidence["disorder"] = torch.zeros_like(summary_confidence["ptm"])
    summary_confidence["ranking_score"] = (
        0.8 * summary_confidence["iptm"]
        + 0.2 * summary_confidence["ptm"]
        + 0.5 * summary_confidence["disorder"]
        - 100 * summary_confidence["has_clash"]
    )
    if interested_atom_mask is not None:
        token_idx = atom_to_token_idx[interested_atom_mask[0].bool()].long()
        asym_ids = token_asym_id[token_idx]
        assert len(torch.unique(asym_ids)) == 1
        interested_asym_id = asym_ids[0].item()
        N_chains = token_asym_id.max().long().item() + 1
        pb_ranking_score = summary_confidence["chain_pair_iptm_global"][
            :, interested_asym_id, torch.arange(N_chains) != interested_asym_id
        ]  # [N_s, N_chain - 1]
        summary_confidence["pb_ranking_score"] = pb_ranking_score[:, 0]
        if elements_one_hot is not None and mol_id is not None:
            vdw_clash = calculate_vdw_clash(
                pred_coordinate=atom_coordinate,
                asym_id=token_asym_id,
                mol_id=mol_id,
                is_polymer=atom_is_polymer,
                atom_token_idx=atom_to_token_idx,
                elements_one_hot=elements_one_hot,
                threshold=configs.metrics.clash.vdw_clash_threshold,
            )
            N_sample = atom_coordinate.shape[0]
            vdw_clash_per_sample_flag = (
                vdw_clash[:, interested_asym_id, :].reshape(N_sample, -1).max(dim=-1)[0]
            )
            summary_confidence["has_vdw_pl_clash"] = vdw_clash_per_sample_flag
            summary_confidence["pb_ranking_score_vdw_penalized"] = (
                summary_confidence["pb_ranking_score"] - 100 * vdw_clash_per_sample_flag
            )

    summary_confidence = break_down_to_per_sample_dict(
        summary_confidence, shared_keys=["num_recycles"]
    )
    torch.cuda.empty_cache()
    if return_full_data:
        # save extra inputs that are used for computing summary_confidence
        full_data["token_has_frame"] = token_has_frame.clone()
        full_data["token_asym_id"] = token_asym_id.clone()
        full_data["atom_to_token_idx"] = atom_to_token_idx.clone()
        full_data["atom_is_polymer"] = atom_is_polymer.clone()
        full_data["atom_coordinate"] = atom_coordinate.clone()

        full_data = break_down_to_per_sample_dict(
            full_data,
            shared_keys=[
                "contact_probs",
                "token_has_frame",
                "token_asym_id",
                "atom_to_token_idx",
                "atom_is_polymer",
            ],
        )
        return summary_confidence, full_data
    else:
        return summary_confidence, [{}]


def get_bin_params(cfg: ConfigDict) -> dict:
    """
    Extract bin parameters from the configuration object.
    """
    return {"min_bin": cfg.min_bin, "max_bin": cfg.max_bin, "no_bins": cfg.no_bins}


def compute_contact_prob(
    distogram_logits: torch.Tensor,
    min_bin: float,
    max_bin: float,
    no_bins: int,
    thres=8.0,
) -> torch.Tensor:
    """
    Compute the contact probability from distogram logits.

    Args:
        distogram_logits (torch.Tensor): Logits for the distogram.
            Shape: [N_token, N_token, N_bins]
        min_bin (float): Minimum bin value.
        max_bin (float): Maximum bin value.
        no_bins (int): Number of bins.
        thres (float): Threshold distance for contact probability. Defaults to 8.0.

    Returns:
        torch.Tensor: Contact probability.
            Shape: [N_token, N_token]
    """
    distogram_prob = torch.nn.functional.softmax(
        distogram_logits, dim=-1
    )  # [N_token, N_token, N_bins]
    distogram_bins = get_bin_centers(min_bin, max_bin, no_bins)
    thres_idx = (distogram_bins < thres).sum()
    contact_prob = distogram_prob[..., :thres_idx].sum(-1)
    return contact_prob


def get_bin_centers(min_bin: float, max_bin: float, no_bins: int) -> torch.Tensor:
    """
    Calculate the centers of the bins for a given range and number of bins.

    Args:
        min_bin (float): The minimum value of the bin range.
        max_bin (float): The maximum value of the bin range.
        no_bins (int): The number of bins.

    Returns:
        torch.Tensor: The centers of the bins.
            Shape: [no_bins]
    """
    bin_width = (max_bin - min_bin) / no_bins
    boundaries = torch.linspace(
        start=min_bin,
        end=max_bin - bin_width,
        steps=no_bins,
    )
    bin_centers = boundaries + 0.5 * bin_width
    return bin_centers


def logits_to_prob(logits: torch.Tensor, dim=-1) -> torch.Tensor:
    return torch.nn.functional.softmax(logits, dim=dim)


def logits_to_score(
    logits: torch.Tensor,
    min_bin: float,
    max_bin: float,
    no_bins: int,
    return_prob=False,
) -> Union[torch.Tensor, tuple[torch.Tensor, torch.Tensor]]:
    """
    Convert logits to a score using bin centers.

    Args:
        logits (torch.Tensor): Logits tensor.
            Shape: [..., no_bins]
        min_bin (float): Minimum bin value.
        max_bin (float): Maximum bin value.
        no_bins (int): Number of bins.
        return_prob (bool): Whether to return the probability distribution. Defaults to False.

    Returns:
        score (torch.Tensor): Converted score.
            Shape: [...]
        prob (torch.Tensor, optional): Probability distribution if `return_prob` is True.
            Shape: [..., no_bins]
    """
    prob = logits_to_prob(logits, dim=-1)
    bin_centers = get_bin_centers(min_bin, max_bin, no_bins).to(logits.device)
    score = prob @ bin_centers
    if return_prob:
        return score, prob
    else:
        return score


def calculate_normalization(N):
    # TM-score normalization constant
    return 1.24 * (max(N, 19) - 15) ** (1 / 3) - 1.8


def calculate_vdw_clash(
    pred_coordinate: torch.Tensor,
    asym_id: torch.LongTensor,
    mol_id: torch.LongTensor,
    atom_token_idx: torch.LongTensor,
    is_polymer: torch.BoolTensor,
    elements_one_hot: torch.Tensor,
    threshold: float,
) -> torch.Tensor:
    """
    Calculate Van der Waals (VDW) clash for predicted coordinates.

    Args:
        pred_coordinate (torch.Tensor): Predicted coordinates of atoms.
            Shape: [N_sample, N_atom, 3]
        asym_id (torch.LongTensor): Asymmetric ID for tokens.
            Shape: [N_token]
        mol_id (torch.LongTensor): Molecular ID.
            Shape: [N_atom]
        atom_token_idx (torch.LongTensor): Mapping from atoms to tokens.
            Shape: [N_atom]
        is_polymer (torch.BoolTensor): Indicator for atoms being part of a polymer.
            Shape: [N_atom]
        elements_one_hot (torch.Tensor): One-hot encoding for elements.
            Shape: [N_atom, N_elements]
        threshold (float): Threshold for VDW clash detection.

    Returns:
        torch.Tensor: VDW clash summary.
            Shape: [N_sample]
    """
    clash_calculator = Clash(vdw_clash_threshold=threshold, compute_af3_clash=False)
    # Check ligand-polymer VDW clash
    N_sample = pred_coordinate.shape[0]
    dummy_is_dna = torch.zeros_like(is_polymer)
    dummy_is_rna = torch.zeros_like(is_polymer)
    clash_dict = clash_calculator(
        pred_coordinate=pred_coordinate,
        asym_id=asym_id,
        atom_to_token_idx=atom_token_idx,
        mol_id=mol_id,
        is_ligand=1 - is_polymer,
        is_protein=is_polymer,
        is_dna=dummy_is_dna,
        is_rna=dummy_is_rna,
        elements_one_hot=elements_one_hot,
    )
    return clash_dict["summary"]["vdw_clash"]


def calculate_clash(
    pred_coordinate: torch.Tensor,
    asym_id: torch.LongTensor,
    atom_to_token_idx: torch.LongTensor,
    is_polymer: torch.BoolTensor,
    threshold: float,
) -> torch.Tensor:
    """Check complex clash

    Args:
        pred_coordinate (torch.Tensor): [N_sample, N_atom, 3]
        asym_id (torch.LongTensor): [N_token, ]
        atom_to_token_idx (torch.LongTensor): [N_atom, ]
        is_polymer (torch.BoolTensor): [N_atom, ]
        threshold: (float)

    Returns:
        torch.Tensor: [N_sample] whether there is a clash in the complex
    """
    N_sample = pred_coordinate.shape[0]
    dummy_is_dna = torch.zeros_like(is_polymer)
    dummy_is_rna = torch.zeros_like(is_polymer)
    clash_calculator = Clash(vdw_clash_threshold=threshold, compute_vdw_clash=False)
    clash_dict = clash_calculator(
        pred_coordinate,
        asym_id,
        atom_to_token_idx,
        1 - is_polymer,
        is_polymer,
        dummy_is_dna,
        dummy_is_rna,
    )
    return clash_dict["summary"]["af3_clash"].reshape(N_sample, -1).max(dim=-1)[0]


def calculate_ptm(
    pae_prob: torch.Tensor,
    has_frame: torch.BoolTensor,
    min_bin: float,
    max_bin: float,
    no_bins: int,
    token_mask: Optional[torch.BoolTensor] = None,
) -> torch.Tensor:
    """Compute pTM score

    Args:
        pae_prob (torch.Tensor): Predicted probability from PAE loss head.
            Shape: [..., N_token, N_token, N_bins]
        has_frame (torch.BoolTensor): Indicator for tokens having a frame.
            Shape: [N_token, ]
        min_bin (float): Minimum bin value.
        max_bin (float): Maximum bin value.
        no_bins (int): Number of bins.
        token_mask (Optional[torch.BoolTensor]): Mask for tokens.
            Shape: [N_token, ] or None

    Returns:
        torch.Tensor: pTM score. Higher values indicate better ranking.
            Shape: [...]
    """
    has_frame = has_frame.bool()

    if token_mask is not None:
        token_mask = token_mask.bool()
        pae_prob = pae_prob[..., token_mask, :, :][
            ..., :, token_mask, :
        ]  # [..., N_d, N_d, N_bins]
        has_frame = has_frame[token_mask]  # [N_d, ]

    if has_frame.sum() == 0:
        return torch.zeros(size=pae_prob.shape[:-3], device=pae_prob.device)

    N_d = has_frame.shape[-1]
    ptm_norm = calculate_normalization(N_d)

    bin_center = get_bin_centers(min_bin, max_bin, no_bins)
    per_bin_weight = (1 / (1 + (bin_center / ptm_norm) ** 2)).to(
        pae_prob.device
    )  # [N_bins]

    token_token_ptm = (pae_prob * per_bin_weight).sum(dim=-1)  # [..., N_d, N_d]

    ptm = token_token_ptm.mean(dim=-1)[..., has_frame].max(dim=-1).values
    return ptm


def calculate_chain_based_ptm(
    pae_prob: torch.Tensor,
    has_frame: torch.BoolTensor,
    asym_id: torch.LongTensor,
    token_is_ligand: torch.BoolTensor,
    min_bin: float,
    max_bin: float,
    no_bins: int,
) -> dict[str, torch.Tensor]:
    """
    Compute chain-based pTM scores.

    Args:
        pae_prob (torch.Tensor): Predicted probability from PAE loss head.
            Shape: [..., N_token, N_token, N_bins]
        has_frame (torch.BoolTensor): Indicator for tokens having a frame.
            Shape: [N_token, ]
        asym_id (torch.LongTensor): Asymmetric ID for tokens.
            Shape: [N_token, ]
        token_is_ligand (torch.BoolTensor): Indicator for tokens being ligands.
            Shape: [N_token, ]
        min_bin (float): Minimum bin value.
        max_bin (float): Maximum bin value.
        no_bins (int): Number of bins.

    Returns:
        dict: Dictionary containing chain-based pTM scores.
            - chain_ptm (torch.Tensor): pTM scores for each chain.
            - chain_iptm (torch.Tensor): ipTM scores for chain interface.
            - chain_pair_iptm (torch.Tensor): Pairwise ipTM scores between chains.
            - chain_pair_iptm_global (torch.Tensor): Global pairwise ipTM scores between chains.
    """

    has_frame = has_frame.bool()
    asym_id = asym_id.long()
    asym_id_to_asym_mask = {aid.item(): asym_id == aid for aid in torch.unique(asym_id)}
    chain_is_ligand = {
        aid.item(): token_is_ligand[asym_id == aid].sum() >= (asym_id == aid).sum() // 2
        for aid in torch.unique(asym_id)
    }

    batch_shape = pae_prob.shape[:-3]

    # Chain_pair_iptm
    # Change to dense tensor, otherwise it's troublesome in break_down_to_per_sample_dict and traverse_and_aggregate across different devices
    N_chain = len(asym_id_to_asym_mask)
    chain_pair_iptm = torch.zeros(size=batch_shape + (N_chain, N_chain)).to(
        pae_prob.device
    )
    for aid_1 in range(N_chain):
        for aid_2 in range(N_chain):
            if aid_1 == aid_2:
                continue
            if aid_1 > aid_2:
                chain_pair_iptm[:, aid_1, aid_2] = chain_pair_iptm[:, aid_2, aid_1]
                continue
            pair_mask = asym_id_to_asym_mask[aid_1] + asym_id_to_asym_mask[aid_2]
            chain_pair_iptm[:, aid_1, aid_2] = calculate_iptm(
                pae_prob,
                has_frame,
                asym_id,
                min_bin,
                max_bin,
                no_bins,
                token_mask=pair_mask,
            )

    # chain_ptm
    chain_ptm = torch.zeros(size=batch_shape + (N_chain,)).to(pae_prob.device)
    for aid, asym_mask in asym_id_to_asym_mask.items():
        chain_ptm[:, aid] = calculate_ptm(
            pae_prob,
            has_frame,
            min_bin,
            max_bin,
            no_bins,
            token_mask=asym_mask,
        )

    # Chain iptm
    chain_has_frame = [
        (asym_id_to_asym_mask[i] * has_frame).any() for i in range(N_chain)
    ]

    chain_iptm = torch.zeros(size=batch_shape + (N_chain,)).to(pae_prob.device)
    for aid, asym_mask in asym_id_to_asym_mask.items():
        pairs = [
            (i, j)
            for i in range(N_chain)
            for j in range(N_chain)
            if (i == aid or j == aid) and (i != j) and chain_has_frame[i]
        ]
        vals = [chain_pair_iptm[:, i, j] for (i, j) in pairs]
        if len(vals) > 0:
            chain_iptm[:, aid] = torch.stack(vals, dim=-1).mean(dim=-1)

    # Chain_pair_iptm_global
    chain_pair_iptm_global = torch.zeros(size=batch_shape + (N_chain, N_chain)).to(
        pae_prob.device
    )
    for aid_1 in range(N_chain):
        for aid_2 in range(N_chain):
            if aid_1 == aid_2:
                continue
            if chain_is_ligand[aid_1]:
                chain_pair_iptm_global[:, aid_1, aid_2] = chain_iptm[:, aid_1]
            elif chain_is_ligand[aid_2]:
                chain_pair_iptm_global[:, aid_1, aid_2] = chain_iptm[:, aid_2]
            else:
                chain_pair_iptm_global[:, aid_1, aid_2] = (
                    chain_iptm[:, aid_1] + chain_iptm[:, aid_2]
                ) * 0.5

    return {
        "chain_ptm": chain_ptm,
        "chain_iptm": chain_iptm,
        "chain_pair_iptm": chain_pair_iptm,
        "chain_pair_iptm_global": chain_pair_iptm_global,
    }


def calculate_chain_based_plddt(
    atom_plddt: torch.Tensor,
    asym_id: torch.LongTensor,
    atom_to_token_idx: torch.LongTensor,
) -> dict[str, torch.Tensor]:
    """
    Calculate chain-based pLDDT scores.

    Args:
        atom_plddt (torch.Tensor): Predicted pLDDT scores for atoms.
            Shape: [N_sample, N_atom]
        asym_id (torch.LongTensor): Asymmetric ID for tokens.
            Shape: [N_token]
        atom_to_token_idx (torch.LongTensor): Mapping from atoms to tokens.
            Shape: [N_atom]

    Returns:
        dict: Dictionary containing chain-based pLDDT scores.
            - chain_plddt (torch.Tensor): pLDDT scores for each chain.
            - chain_pair_plddt (torch.Tensor): Pairwise pLDDT scores between chains.
    """

    asym_id = asym_id.long()
    asym_id_to_asym_mask = {aid.item(): asym_id == aid for aid in torch.unique(asym_id)}
    N_chain = len(asym_id_to_asym_mask)
    assert N_chain == asym_id.max() + 1  # make sure it is from 0 to N_chain-1

    def _calculate_lddt_with_token_mask(token_mask):
        atom_mask = token_mask[atom_to_token_idx]
        sub_plddt = atom_plddt[:, atom_mask].mean(-1)
        return sub_plddt

    batch_shape = atom_plddt.shape[:-1]
    # Chain_plddt
    chain_plddt = torch.zeros(size=batch_shape + (N_chain,)).to(atom_plddt.device)
    for aid, asym_mask in asym_id_to_asym_mask.items():
        chain_plddt[:, aid] = _calculate_lddt_with_token_mask(token_mask=asym_mask)

    # Chain_pair_plddt
    chain_pair_plddt = torch.zeros(size=batch_shape + (N_chain, N_chain)).to(
        atom_plddt.device
    )
    for aid_1 in asym_id_to_asym_mask:
        for aid_2 in asym_id_to_asym_mask:
            if aid_1 == aid_2:
                continue
            pair_mask = asym_id_to_asym_mask[aid_1] + asym_id_to_asym_mask[aid_2]
            chain_pair_plddt[:, aid_1, aid_2] = _calculate_lddt_with_token_mask(
                token_mask=pair_mask
            )

    return {"chain_plddt": chain_plddt, "chain_pair_plddt": chain_pair_plddt}


def calculate_iptm(
    pae_prob: torch.Tensor,
    has_frame: torch.BoolTensor,
    asym_id: torch.LongTensor,
    min_bin: float,
    max_bin: float,
    no_bins: int,
    token_mask: Optional[torch.BoolTensor] = None,
    eps: float = 1e-8,
):
    """
    Compute ipTM score.

    Args:
        pae_prob (torch.Tensor): Predicted probability from PAE loss head.
            Shape: [..., N_token, N_token, N_bins]
        has_frame (torch.BoolTensor): Indicator for tokens having a frame.
            Shape: [N_token, ]
        asym_id (torch.LongTensor): Asymmetric ID for tokens.
            Shape: [N_token, ]
        min_bin (float): Minimum bin value.
        max_bin (float): Maximum bin value.
        no_bins (int): Number of bins.
        token_mask (Optional[torch.BoolTensor]): Mask for tokens.
            Shape: [N_token, ] or None
        eps (float): Small value to avoid division by zero. Defaults to 1e-8.

    Returns:
        torch.Tensor: ipTM score. Higher values indicate better ranking.
            Shape: [...]
    """
    has_frame = has_frame.bool()
    if token_mask is not None:
        token_mask = token_mask.bool()
        pae_prob = pae_prob[..., token_mask, :, :][
            ..., :, token_mask, :
        ]  # [..., N_d, N_d, N_bins]
        has_frame = has_frame[token_mask]  # [N_d, ]
        asym_id = asym_id[token_mask]  # [N_d, ]

    if has_frame.sum() == 0:
        return torch.zeros(size=pae_prob.shape[:-3], device=pae_prob.device)

    N_d = has_frame.shape[-1]
    ptm_norm = calculate_normalization(N_d)

    bin_center = get_bin_centers(min_bin, max_bin, no_bins)
    per_bin_weight = (1 / (1 + (bin_center / ptm_norm) ** 2)).to(
        pae_prob.device
    )  # [N_bins]

    token_token_ptm = (pae_prob * per_bin_weight).sum(dim=-1)  # [..., N_d, N_d]

    is_diff_chain = asym_id[None, :] != asym_id[:, None]  # [N_d, N_d]

    iptm = (token_token_ptm * is_diff_chain).sum(dim=-1) / (
        eps + is_diff_chain.sum(dim=-1)
    )  # [..., N_d]
    iptm = iptm[..., has_frame].max(dim=-1).values

    return iptm


def break_down_to_per_sample_dict(input_dict: dict, shared_keys=[]) -> list[dict]:
    """
    Break down a dictionary containing tensors into a list of dictionaries, each corresponding to a sample.

    Args:
        input_dict (dict): Dictionary containing tensors.
        shared_keys (list): List of keys that are shared across all samples. Defaults to an empty list.

    Returns:
        list[dict]: List of dictionaries, each containing data for a single sample.
    """
    per_sample_keys = [key for key in input_dict if key not in shared_keys]
    assert len(per_sample_keys) > 0
    N_sample = input_dict[per_sample_keys[0]].size(0)
    for key in per_sample_keys:
        assert input_dict[key].size(0) == N_sample

    per_sample_dict_list = []
    for i in range(N_sample):
        sample_dict = {key: input_dict[key][i] for key in per_sample_keys}
        sample_dict.update({key: input_dict[key] for key in shared_keys})
        per_sample_dict_list.append(sample_dict)

    return per_sample_dict_list


@torch.no_grad()
def compute_full_data_and_summary(
    configs,
    pae_logits,
    plddt_logits,
    pde_logits,
    contact_probs,
    token_asym_id,
    token_has_frame,
    atom_coordinate,
    atom_to_token_idx,
    atom_is_polymer,
    N_recycle,
    return_full_data: bool = False,
    interested_atom_mask=None,
    mol_id=None,
    elements_one_hot=None,
):
    """Wrapper of `_compute_full_data_and_summary` by enumerating over N samples"""

    N_sample = pae_logits.size(0)
    if contact_probs.dim() == 2:
        # Convert to [N_sample, N_token, N_token]
        contact_probs = contact_probs.unsqueeze(dim=0).expand(N_sample, -1, -1)
    else:
        assert contact_probs.dim() == 3
    assert (
        contact_probs.size(0) == plddt_logits.size(0) == pde_logits.size(0) == N_sample
    )

    summary_confidence = []
    full_data = []
    for i in range(N_sample):
        summary_confidence_i, full_data_i = _compute_full_data_and_summary(
            configs=configs,
            pae_logits=pae_logits[i : i + 1],
            plddt_logits=plddt_logits[i : i + 1],
            pde_logits=pde_logits[i : i + 1],
            contact_probs=contact_probs[i],
            token_asym_id=token_asym_id,
            token_has_frame=token_has_frame,
            atom_coordinate=atom_coordinate[i : i + 1],
            atom_to_token_idx=atom_to_token_idx,
            atom_is_polymer=atom_is_polymer,
            N_recycle=N_recycle,
            interested_atom_mask=interested_atom_mask,
            return_full_data=return_full_data,
            mol_id=mol_id,
            elements_one_hot=elements_one_hot,
        )
        summary_confidence.extend(summary_confidence_i)
        full_data.extend(full_data_i)
    return summary_confidence, full_data
