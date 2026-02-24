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

import dataclasses
import json
import logging
import os
import pickle
import random
import time
from datetime import datetime, timezone
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

import requests

import numpy as np
from typing_extensions import Final, TypeAlias

from protenix.data.constants import (
    ATOM37_NUM,
    ATOM37_ORDER,
    PROTEIN_AATYPE_DENSE_ATOM_TO_ATOM37,
    PROTEIN_CHAIN,
    RESTYPE_PSEUDOBETA_INDEX,
    RESTYPE_RIGIDGROUP_DENSE_ATOM_IDX,
    STD_RESIDUES_WITH_GAP,
)
from protenix.data.template.template_parser import (
    AlignRatioError,
    CaDistanceError,
    DateError,
    DuplicateError,
    encode_template_restype,
    get_pdb_id_and_chain,
    LengthError,
    MmcifObject,
    MultipleChainsError,
    NoAtomDataInTemplateError,
    NoChainsError,
    PrefilterError,
    PrefilterResult,
    QueryToTemplateAlignError,
    SingleHitResult,
    TemplateAtomMaskAllZerosError,
    TemplateHit,
    TemplateParser,
    TemplateSearchResult,
)
from protenix.data.tools.kalign import Kalign
from protenix.utils.logger import get_logger

logger = get_logger(__name__)
logger.setLevel(logging.INFO)

# Global caches for release dates and obsolete PDBs to avoid redundant loading.
_RELEASE_DATES_CACHE: Dict[str, Mapping[str, datetime]] = {}
_OBSOLETE_PDBS_CACHE: Dict[str, Dict[str, str]] = {}

FeatureDict: TypeAlias = Mapping[str, np.ndarray]

DAYS_BEFORE_QUERY_DATE: Final[int] = 60

TEMPLATE_FEATURES: Final[Tuple[str, ...]] = (
    "template_aatype",
    "template_atom_positions",
    "template_atom_mask",
)

_POLYMER_FEATURES: Final[Mapping[str, Union[np.float64, np.int32, object]]] = {
    "template_aatype": np.int32,
    "template_all_atom_masks": np.float64,
    "template_all_atom_positions": np.float64,
    "template_domain_names": object,
    "template_release_date": object,
    "template_sequence": object,
}


@dataclasses.dataclass(frozen=True)
class DistogramFeaturesConfig:
    # The left edge of the first bin.
    min_bin: float = 3.25
    # The left edge of the final bin. The final bin catches everything larger than
    # `max_bin`.
    max_bin: float = 50.75
    # The number of bins in the distogram.
    num_bins: int = 39


class TemplateFeatures:
    """Utility functions for AlphaFold 3 template processing."""

    @staticmethod
    def get_timestamp(date_str: str) -> float:
        """Converts an ISO format date string to a UTC timestamp."""
        dt = datetime.fromisoformat(date_str)
        dt = dt.replace(tzinfo=timezone.utc)
        return dt.timestamp()

    @staticmethod
    def package_template_features(
        *, hit_features: Sequence[Mapping[str, Any]]
    ) -> Mapping[str, Any]:
        """Stacks polymer features, adds empty and keeps ligand features unstacked."""
        features_to_include = set(_POLYMER_FEATURES)
        features = {
            feat: [single_hit_features[feat] for single_hit_features in hit_features]
            for feat in features_to_include
        }

        stacked_features = {}
        for k, v in features.items():
            if k in _POLYMER_FEATURES:
                v = (
                    np.stack(v, axis=0)
                    if v
                    else np.array([], dtype=_POLYMER_FEATURES[k])
                )
            stacked_features[k] = v

        return stacked_features

    @staticmethod
    def fix_template_features(
        template_features: FeatureDict, num_res: int
    ) -> FeatureDict:
        """Convert template features to AlphaFold 3 format."""
        if not template_features["template_aatype"].shape[0]:
            template_features = TemplateFeatures.empty_template_features(num_res)
        else:
            template_release_timestamp = [
                TemplateFeatures.get_timestamp(x.decode("utf-8"))
                for x in template_features["template_release_date"]
            ]

            # Convert from atom37 to dense atom
            dense_atom_indices = np.take(
                PROTEIN_AATYPE_DENSE_ATOM_TO_ATOM37,
                template_features["template_aatype"],
                axis=0,
            )

            atom_mask = np.take_along_axis(
                template_features["template_all_atom_masks"], dense_atom_indices, axis=2
            )
            atom_positions = np.take_along_axis(
                template_features["template_all_atom_positions"],
                dense_atom_indices[..., None],
                axis=2,
            )
            atom_positions *= atom_mask[..., None]

            template_features = {
                "template_aatype": template_features["template_aatype"],
                "template_atom_mask": atom_mask.astype(np.int32),
                "template_atom_positions": atom_positions.astype(np.float32),
                "template_domain_names": np.array(
                    template_features["template_domain_names"], dtype=object
                ),
                "template_release_timestamp": np.array(
                    template_release_timestamp, dtype=np.float32
                ),
            }
        return template_features

    @staticmethod
    def empty_template_features(num_res: int, num_dense: int = 24) -> FeatureDict:
        """Creates fully masked out template features."""
        template_features = {
            "template_aatype": np.array(
                [STD_RESIDUES_WITH_GAP["-"]] * num_res, dtype=np.int32
            )[None, ...],
            "template_atom_mask": np.zeros((num_res, num_dense), dtype=np.int32)[
                None, ...
            ],
            "template_atom_positions": np.zeros(
                (num_res, num_dense, 3), dtype=np.float32
            )[None, ...],
            "template_domain_names": np.array([b""], dtype=object),
            "template_release_timestamp": np.array([0.0], dtype=np.float32),
        }
        return template_features

    @staticmethod
    def pseudo_beta_fn(
        aatype: np.ndarray,
        dense_atom_positions: np.ndarray,
        dense_atom_masks: np.ndarray,
        is_ligand: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Create pseudo beta atom positions and optionally mask."""
        if is_ligand is None:
            is_ligand = np.zeros_like(aatype)

        pseudobeta_index_polymer = np.take(
            RESTYPE_PSEUDOBETA_INDEX, aatype, axis=0
        ).astype(np.int32)

        pseudobeta_index = np.where(
            is_ligand,
            np.zeros_like(pseudobeta_index_polymer),
            pseudobeta_index_polymer,
        )

        pseudo_beta = np.take_along_axis(
            dense_atom_positions, pseudobeta_index[..., None, None], axis=-2
        )
        pseudo_beta = np.squeeze(pseudo_beta, axis=-2)

        pseudo_beta_mask = np.take_along_axis(
            dense_atom_masks, pseudobeta_index[..., None], axis=-1
        ).astype(np.float32)
        pseudo_beta_mask = np.squeeze(pseudo_beta_mask, axis=-1)

        return pseudo_beta, pseudo_beta_mask

    @staticmethod
    def dgram_from_positions(
        positions: np.ndarray, config: DistogramFeaturesConfig
    ) -> np.ndarray:
        """Compute distogram from amino acid positions."""
        lower_breaks = np.linspace(config.min_bin, config.max_bin, config.num_bins)
        lower_breaks = np.square(lower_breaks)
        upper_breaks = np.concatenate(
            [lower_breaks[1:], np.array([1e8], dtype=np.float32)], axis=-1
        )
        dist2 = np.sum(
            np.square(
                np.expand_dims(positions, axis=-2) - np.expand_dims(positions, axis=-3)
            ),
            axis=-1,
            keepdims=True,
        )

        dgram = (dist2 > lower_breaks).astype(np.float32) * (
            dist2 < upper_breaks
        ).astype(np.float32)
        return dgram

    @staticmethod
    def compute_template_unit_vector(
        aatype: np.ndarray,
        atom_positions: np.ndarray,
        atom_mask: np.ndarray,
        epsilon: float = 1e-6,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Simplified calculation of template unit vector."""
        # Get backbone indices (C, CA, N) for each residue from protein_data_processing
        # Group 0: [C, CA, N]
        backbone_indices = RESTYPE_RIGIDGROUP_DENSE_ATOM_IDX[aatype, 0]  # [num_res, 3]

        # Indices according to protein_data_processing.py: C is 0, CA is 1, N is 2
        c_idx = backbone_indices[:, 0]
        ca_idx = backbone_indices[:, 1]
        n_idx = backbone_indices[:, 2]

        num_res = aatype.shape[0]
        res_indices = np.arange(num_res)

        c_pos = atom_positions[res_indices, c_idx]
        ca_pos = atom_positions[res_indices, ca_idx]
        n_pos = atom_positions[res_indices, n_idx]

        c_mask = atom_mask[res_indices, c_idx]
        ca_mask = atom_mask[res_indices, ca_idx]
        n_mask = atom_mask[res_indices, n_idx]

        mask = (c_mask * ca_mask * n_mask).astype(np.float32)

        # Local frame: origin at CA
        # x-axis along C-CA (following original code convention)
        v1 = c_pos - ca_pos
        v2 = n_pos - ca_pos

        e1 = v1 / (np.linalg.norm(v1, axis=-1, keepdims=True) + epsilon)
        # Orthogonalize v2 against e1
        e2 = v2 - np.sum(v2 * e1, axis=-1, keepdims=True) * e1
        e2 = e2 / (np.linalg.norm(e2, axis=-1, keepdims=True) + epsilon)
        # e3 = e1 x e2
        e3 = np.cross(e1, e2)

        # Relative positions of all CA atoms to all local frames
        # diff[i, j] = CA[j] - CA[i]
        diff = ca_pos[None, :, :] - ca_pos[:, None, :]  # [num_res, num_res, 3]

        # Transform to local frame: P' = R^T @ diff
        # R = [e1 | e2 | e3] -> x' = e1 . diff, y' = e2 . diff, z' = e3 . diff
        ux = np.sum(e1[:, None, :] * diff, axis=-1)
        uy = np.sum(e2[:, None, :] * diff, axis=-1)
        uz = np.sum(e3[:, None, :] * diff, axis=-1)

        unit_vector = np.stack([ux, uy, uz], axis=-1)  # [num_res, num_res, 3]

        # Normalize to unit vector
        uv_norm = np.linalg.norm(unit_vector, axis=-1, keepdims=True)
        unit_vector = unit_vector / (uv_norm + epsilon)

        # 2D mask
        mask_2d = mask[:, None] * mask[None, :]

        return unit_vector, mask_2d


# --- Filtering and Processing ---


class TemplateHitFilter:
    """Filters template hits based on date, alignment ratio, and other criteria."""

    def __init__(
        self,
        release_dates: Mapping[str, datetime],
        obsolete_pdbs: Mapping[str, str],
        strict: bool = False,
    ):
        self.release_dates = release_dates
        self.obsolete_pdbs = obsolete_pdbs
        self.strict = strict

    def _is_after_cutoff(self, pdb_id: str, cutoff: Optional[datetime]) -> bool:
        """Checks if a PDB entry was released after the cutoff date."""
        if cutoff is None:
            return False
        date = self.release_dates.get(pdb_id)
        return date > cutoff if date else False

    def _assess_hit(
        self,
        hit: TemplateHit,
        pdb_code: str,
        query_seq: str,
        cutoff: datetime,
        max_subseq_ratio: float = 0.95,
        min_align_ratio: float = 0.1,
    ) -> bool:
        """Performs quick pre-filtering on a hit."""
        align_ratio = hit.aligned_cols / len(query_seq)
        t_seq = hit.hit_sequence.replace("-", "")
        len_ratio = len(t_seq) / len(query_seq)

        if self._is_after_cutoff(pdb_code, cutoff):
            raise DateError(f"Release date for {pdb_code} is after cutoff.")
        if align_ratio <= min_align_ratio:
            raise AlignRatioError(f"Align ratio {align_ratio:.2f} <= {min_align_ratio}")
        if t_seq in query_seq and len_ratio > max_subseq_ratio:
            raise DuplicateError("Hit is a large duplicate of the query.")
        if len(t_seq) < 10:
            raise LengthError("Template sequence too short.")
        return True

    def prefilter(
        self,
        query_seq: str,
        hit: TemplateHit,
        max_date: datetime,
    ) -> PrefilterResult:
        """Prefilters a hit and handles obsolete PDBs."""
        try:
            pdb_id, chain_id = get_pdb_id_and_chain(hit)
        except ValueError as e:
            return PrefilterResult(False, error=str(e))

        if pdb_id not in self.release_dates:
            pdb_id = self.obsolete_pdbs.get(pdb_id, pdb_id)
            if pdb_id not in self.release_dates:
                return PrefilterResult(False, error=f"{pdb_id} date unknown.")

        try:
            self._assess_hit(hit, pdb_id, query_seq, max_date)
        except PrefilterError as e:
            msg = f"Hit {pdb_id}_{chain_id} failed prefilter: {e}"
            if self.strict and isinstance(e, (DateError, DuplicateError)):
                return PrefilterResult(False, error=msg)
            return PrefilterResult(False, warning=msg)
        return PrefilterResult(True)


class TemplateHitProcessor:
    """Processes a template hit to generate features."""

    _PDBE_CIF_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"

    def __init__(
        self,
        mmcif_dir: str,
        template_cache_dir: Optional[str] = None,
        kalign_binary_path: Optional[str] = None,
        _zero_center_positions: bool = True,
        fetch_remote: bool = False,
    ):
        self._mmcif_dir = mmcif_dir
        self._template_cache_dir = template_cache_dir
        self._kalign_binary_path = kalign_binary_path
        self._zero_center_positions = _zero_center_positions
        self._fetch_remote = fetch_remote

    def _read_file(self, path: str) -> str:
        """Reads file content."""
        with open(path, "r") as f:
            return f.read()

    def _fetch_or_read_cif(self, pdb_id: str) -> str:
        """Read mmCIF from local dir, falling back to PDBe API if fetch_remote is enabled."""
        cif_path = os.path.join(self._mmcif_dir, f"{pdb_id}.cif")
        if os.path.exists(cif_path):
            return self._read_file(cif_path)

        if not self._fetch_remote:
            raise FileNotFoundError(f"CIF not found: {cif_path}")

        url = self._PDBE_CIF_URL.format(pdb_id=pdb_id.lower())
        logger.info(f"Fetching mmCIF for {pdb_id} from {url}")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        cif_str = response.text

        os.makedirs(self._mmcif_dir, exist_ok=True)
        with open(cif_path, "w") as f:
            f.write(cif_str)
        logger.info(f"Cached mmCIF for {pdb_id} at {cif_path}")

        return cif_str

    def _get_atom_coords(
        self,
        mmcif_object: MmcifObject,
        chain_id: str,
        _zero_center_positions: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Extracts all-atom coordinates and mask for a specific chain."""
        chains = list(mmcif_object.structure.get_chains())
        relevant_chains = [c for c in chains if c.id == chain_id]
        if len(relevant_chains) != 1:
            raise MultipleChainsError(
                f"Expected 1 chain with id {chain_id}, found {len(relevant_chains)}"
            )
        chain = relevant_chains[0]

        num_res = len(mmcif_object.chain_to_seqres[chain_id])
        all_pos = np.zeros((num_res, ATOM37_NUM, 3), dtype=np.float32)
        all_mask = np.zeros((num_res, ATOM37_NUM), dtype=np.float32)

        for i in range(num_res):
            res_info = mmcif_object.seqres_to_structure[chain_id][i]
            if not res_info.is_missing:
                try:
                    res = chain[
                        (
                            res_info.hetflag,
                            res_info.position.residue_number,
                            res_info.position.insertion_code,
                        )
                    ]
                    for atom in res.get_atoms():
                        name = atom.get_name()
                        coord = atom.get_coord()
                        if name in ATOM37_ORDER:
                            idx = ATOM37_ORDER[name]
                            all_pos[i, idx] = coord
                            all_mask[i, idx] = 1.0
                        elif name.upper() == "SE" and res.get_resname() == "MSE":
                            idx = ATOM37_ORDER["SD"]
                            all_pos[i, idx] = coord
                            all_mask[i, idx] = 1.0

                    # Correct Arginine NH1/NH2 if swapped based on distance to CD.
                    cd, nh1, nh2 = (
                        ATOM37_ORDER["CD"],
                        ATOM37_ORDER["NH1"],
                        ATOM37_ORDER["NH2"],
                    )
                    if res.get_resname() == "ARG" and all(all_mask[i, [cd, nh1, nh2]]):
                        if np.linalg.norm(
                            all_pos[i, nh1] - all_pos[i, cd]
                        ) > np.linalg.norm(all_pos[i, nh2] - all_pos[i, cd]):
                            all_pos[i, [nh1, nh2]] = all_pos[i, [nh2, nh1]]
                            all_mask[i, [nh1, nh2]] = all_mask[i, [nh2, nh1]]
                except KeyError:
                    continue

        if _zero_center_positions:
            mask_bool = all_mask.astype(bool)
            if np.any(mask_bool):
                center = all_pos[mask_bool].mean(axis=0)
                all_pos[mask_bool] -= center

        return all_pos, all_mask

    def _check_residue_distances(
        self, pos: np.ndarray, mask: np.ndarray, max_dist: float
    ):
        """Verifies that distance between consecutive CA atoms is within limits."""
        ca_idx = ATOM37_ORDER["CA"]
        prev_pos = None
        for i, (p, m) in enumerate(zip(pos, mask)):
            if m[ca_idx]:
                curr_pos = p[ca_idx]
                if prev_pos is not None:
                    dist = np.linalg.norm(curr_pos - prev_pos)
                    if dist > max_dist:
                        raise CaDistanceError(
                            f"Distance between residues {i} and previous is {dist:.2f} > {max_dist}"
                        )
                prev_pos = curr_pos

    def _get_atom_positions(
        self,
        mmcif_obj: MmcifObject,
        auth_chain_id: str,
        max_ca_dist: float,
        _zero_center: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        pos, mask = self._get_atom_coords(mmcif_obj, auth_chain_id, _zero_center)
        self._check_residue_distances(pos, mask, max_ca_dist)
        return pos, mask

    def _align_query_to_hit_index_mapping(
        self,
        query_seq: str,
        chain_id: str,
        mmcif_obj: MmcifObject,
        min_ratio: float,
    ) -> Tuple[str, Dict[int, int], str]:
        """
        Aligns the query sequence directly to the chain sequence in mmCIF and establishes index mapping.

        This method ensures that each residue in the query sequence finds its corresponding
        position in the template structure sequence via alignment.

        Args:
            query_seq: The query sequence string.
            chain_id: The ID of the target chain.
            mmcif_obj: The parsed mmCIF object.
            min_ratio: Minimum allowed sequence identity ratio.

        Returns:
            target_seq: The original chain sequence from mmCIF.
            mapping: A dictionary establishing the mapping from query sequence index to
                    template sequence index (Query Index -> Template Index).
                    key: Residue index in the query sequence (0-based).
                    value: Residue index in the template sequence (0-based).
            actual_chain_id: The actual chain ID used from the mmCIF.
        """
        aligner = Kalign(binary_path=self._kalign_binary_path)
        target_seq = mmcif_obj.chain_to_seqres.get(chain_id, "")
        actual_chain_id = chain_id
        if not target_seq:
            if len(mmcif_obj.chain_to_seqres) == 1:
                actual_chain_id = list(mmcif_obj.chain_to_seqres.keys())[0]
                target_seq = mmcif_obj.chain_to_seqres[actual_chain_id]
            else:
                raise QueryToTemplateAlignError(f"Chain {chain_id} not found.")

        try:
            q_aln, t_aln = aligner.align([query_seq, target_seq])
        except Exception as e:
            raise QueryToTemplateAlignError(f"Query alignment failed: {e}")

        mapping = {}
        q_idx, t_idx, same, count = -1, -1, 0, 0
        for qa, ta in zip(q_aln, t_aln):
            if qa != "-":
                q_idx += 1
            if ta != "-":
                t_idx += 1
            if qa != "-" and ta != "-":
                # Establish mapping from Query index to Template (mmCIF) index
                mapping[q_idx] = t_idx
                count += 1
                if qa == ta:
                    same += 1

        if count > 0 and (same / count < min_ratio):
            raise QueryToTemplateAlignError("Insufficient similarity to query.")

        return target_seq, mapping, actual_chain_id

    def _extract_template_features(
        self,
        mmcif_obj: MmcifObject,
        pdb_id: str,
        mapping: Mapping[int, int],
        template_seq: str,
        query_seq: str,
        chain_id: str,
        _zero_center: bool = True,
    ) -> Tuple[Dict[str, Any], Optional[str]]:
        """Generates features for a template hit."""
        if not mmcif_obj or not mmcif_obj.chain_to_seqres:
            raise NoChainsError(f"No chains in {pdb_id}")

        warning = None
        try:
            all_pos, all_mask = self._get_atom_positions(
                mmcif_obj, chain_id, 150.0, _zero_center
            )
        except (CaDistanceError, KeyError) as e:
            raise NoAtomDataInTemplateError(f"Failed to get atom data: {e}")

        num_query = len(query_seq)
        out_pos = np.zeros((num_query, ATOM37_NUM, 3), dtype=np.float32)
        out_mask = np.zeros((num_query, ATOM37_NUM), dtype=np.float32)
        out_seq = ["-"] * num_query

        for q_idx, t_idx in mapping.items():
            if t_idx != -1:
                out_pos[q_idx] = all_pos[t_idx]
                out_mask[q_idx] = all_mask[t_idx]
                out_seq[q_idx] = template_seq[t_idx]

        if np.sum(out_mask) < 5:
            raise TemplateAtomMaskAllZerosError(
                f"Empty atom mask for {pdb_id}_{chain_id}"
            )

        out_seq_str = "".join(out_seq)
        aatype = encode_template_restype(PROTEIN_CHAIN, out_seq_str)

        features = {
            "template_all_atom_positions": out_pos,
            "template_all_atom_masks": out_mask,
            "template_sequence": out_seq_str.encode(),
            "template_aatype": np.array(aatype, dtype=np.int32),
            "template_domain_names": np.array(
                f"{pdb_id.lower()}_{chain_id}".encode(), dtype=object
            ),
        }
        return features, warning

    def _update_realigned_hit(
        self, hit: TemplateHit, seq: str, mapping: Mapping[int, int]
    ) -> TemplateHit:
        """
        Updates the TemplateHit object with new alignment information.

        This method synchronizes the hit information (alignment columns, hit sequence,
        and index mappings) after a re-alignment between the query and the template structure.

        Args:
            hit: The original TemplateHit object.
            seq: The new hit sequence. Note: This sequence MUST NOT contain GAPs ('-').
            mapping: A dictionary mapping Query Index -> Template Index.

        Returns:
            A new TemplateHit object with updated alignment details.
        """
        q_indices = list(range(len(hit.query)))
        h_indices = [-1] * len(hit.query)
        for i, j in mapping.items():
            h_indices[i] = j
        return TemplateHit(
            index=hit.index,
            name=hit.name,
            aligned_cols=len(mapping),
            sum_probs=hit.sum_probs,
            query=hit.query,
            hit_sequence=seq,  # GAP-free sequence
            indices_query=q_indices,
            indices_hit=h_indices,
        )

    def process(
        self,
        query_seq: str,
        hit: TemplateHit,
        max_date: datetime,
        release_dates: Mapping[str, datetime],
        obsolete_pdbs: Mapping[str, str],
        strict: bool = False,
    ) -> Tuple[SingleHitResult, Dict[str, float]]:
        """Processes a single hit to generate features."""
        track = {}
        t_start = time.time()
        pdb_id, chain_id = get_pdb_id_and_chain(hit)
        pdb_id = obsolete_pdbs.get(pdb_id, pdb_id)

        # Load mmCIF parsing result.
        if self._template_cache_dir:
            cache_path = os.path.join(
                self._template_cache_dir, f"{pdb_id}_{chain_id}.pkl"
            )
            if os.path.exists(cache_path):
                with open(cache_path, "rb") as f:
                    res = pickle.load(f)
                track["load"] = time.time() - t_start
            else:
                return SingleHitResult(None, None, None, None), track
        else:
            try:
                cif_str = self._fetch_or_read_cif(pdb_id)
                res = TemplateParser.parse(
                    file_id=pdb_id, mmcif_string=cif_str, auth_chain_id=chain_id
                )
                track["load"] = time.time() - t_start
            except (FileNotFoundError, requests.RequestException) as e:
                return (
                    SingleHitResult(None, None, f"CIF load failed for {pdb_id}: {e}", None),
                    track,
                )

        if not res.mmcif_object:
            return (
                SingleHitResult(None, None, f"Parsing failed for {pdb_id}", None),
                track,
            )

        date_str = res.mmcif_object.header.get("release_date", "9999-12-31")
        hit_date = datetime.strptime(date_str, "%Y-%m-%d")
        if hit_date > max_date:
            err = f"Hit {pdb_id} date {hit_date} > {max_date}"
            return (
                SingleHitResult(None, None, err if strict else None, None),
                track,
            )

        # Re-align and extract features.
        t_aln = time.time()
        try:
            (
                target_seq,
                mapping,
                actual_chain_id,
            ) = self._align_query_to_hit_index_mapping(
                query_seq, chain_id, res.mmcif_object, 0.0
            )
            hit = self._update_realigned_hit(hit, target_seq, mapping)
            track["align"] = time.time() - t_aln

            t_feat = time.time()
            feats, warn = self._extract_template_features(
                res.mmcif_object,
                pdb_id,
                mapping,
                target_seq,
                query_seq,
                actual_chain_id,
                self._zero_center_positions,
            )
            feats["template_sum_probs"] = [
                hit.sum_probs if hit.sum_probs is not None else 0.0
            ]
            feats["template_release_date"] = np.array(date_str.encode(), dtype=object)
            track["feature"] = time.time() - t_feat
            return SingleHitResult(hit, feats, None, warn), track
        except Exception as e:
            return (
                SingleHitResult(None, None, f"Error processing hit: {e}", None),
                track,
            )


class TemplateHitFeaturizer:
    """
    Featurizes template hits (e.g., from hmmsearch or HHSearch) into features.

    Args:
        mmcif_dir: Directory containing mmCIF files.
        template_cache_dir: Optional directory to cache parsed mmCIF objects.
        max_hits: Maximum number of hits to process.
        kalign_binary_path: Path to the kalign binary.
        max_template_date: Optional cutoff date for templates.
        release_dates_path: Path to the JSON file containing PDB release dates.
        obsolete_pdbs_path: Path to the JSON file containing obsolete PDB mappings.
        strict_error_check: If True, raise errors for filtering failures.
        _shuffle_top_k_prefiltered: If set, shuffle the top K hits after pre-filtering.
        _zero_center_positions: If True, center atom positions at zero.
        _max_template_candidates_num: Maximum number of candidate hits to consider.
    """

    def __init__(
        self,
        mmcif_dir: str,
        template_cache_dir: Optional[str] = None,
        max_hits: Optional[int] = None,
        kalign_binary_path: Optional[str] = None,
        max_template_date: Optional[Union[str, datetime]] = None,
        release_dates_path: Optional[str] = None,
        obsolete_pdbs_path: Optional[str] = None,
        strict_error_check: bool = False,
        _shuffle_top_k_prefiltered: Optional[int] = None,
        _zero_center_positions: bool = True,
        _max_template_candidates_num: Optional[int] = None,
        fetch_remote: bool = False,
    ):
        self._mmcif_dir = mmcif_dir
        self._template_cache_dir = template_cache_dir
        self._max_hits = max_hits
        self._kalign_binary_path = kalign_binary_path
        self._strict_error_check = strict_error_check
        self._shuffle_top_k_prefiltered = _shuffle_top_k_prefiltered
        self._zero_center_positions = _zero_center_positions
        self._max_template_candidates_num = _max_template_candidates_num
        self._fetch_remote = fetch_remote

        if max_template_date:
            if isinstance(max_template_date, str):
                self._max_template_date = datetime.strptime(
                    max_template_date, "%Y-%m-%d"
                )
            else:
                self._max_template_date = max_template_date
        else:
            self._max_template_date = None

        # Load or reuse release dates.
        if release_dates_path:
            if release_dates_path not in _RELEASE_DATES_CACHE:
                _RELEASE_DATES_CACHE[release_dates_path] = self._parse_release_dates(
                    release_dates_path
                )
            self._release_dates = _RELEASE_DATES_CACHE[release_dates_path]
        else:
            self._release_dates = {}

        # Load or reuse obsolete PDB mapping.
        if obsolete_pdbs_path and os.path.exists(obsolete_pdbs_path):
            if obsolete_pdbs_path not in _OBSOLETE_PDBS_CACHE:
                with open(obsolete_pdbs_path, "r") as f:
                    _OBSOLETE_PDBS_CACHE[obsolete_pdbs_path] = json.load(f)
            self._obsolete_pdbs = _OBSOLETE_PDBS_CACHE[obsolete_pdbs_path]
        else:
            self._obsolete_pdbs = {}

        self._hit_filter = TemplateHitFilter(
            release_dates=self._release_dates,
            obsolete_pdbs=self._obsolete_pdbs,
            strict=self._strict_error_check,
        )

        self._hit_processor = TemplateHitProcessor(
            mmcif_dir=self._mmcif_dir,
            template_cache_dir=self._template_cache_dir,
            kalign_binary_path=self._kalign_binary_path,
            _zero_center_positions=self._zero_center_positions,
            fetch_remote=self._fetch_remote,
        )

    def _parse_release_dates(self, path: str) -> Dict[str, datetime]:
        if not os.path.exists(path):
            return {}
        with open(path, "r") as f:
            data = json.load(f)
        return {
            pdb.lower(): datetime.strptime(v["release_date"], "%Y-%m-%d")
            for pdb, v in data.items()
            if "release_date" in v
        }

    def get_templates(
        self,
        sequence_uid: str,
        query_sequence: str,
        hits: Sequence[TemplateHit],
        max_template_date: Optional[Union[str, datetime]] = None,
    ) -> Tuple[TemplateSearchResult, Dict[str, float]]:
        """
        Processes hits to generate template features.

        Args:
            sequence_uid: Unique identifier for the sequence.
            query_sequence: The query sequence.
            hits: A sequence of TemplateHit objects.
            max_template_date: Optional cutoff date to override the default.

        Returns:
            A tuple of (TemplateSearchResult, timing_dict).
        """
        cutoff = self._max_template_date
        if max_template_date:
            if isinstance(max_template_date, str):
                cutoff = datetime.strptime(max_template_date, "%Y-%m-%d")
            else:
                cutoff = max_template_date

        # Prefilter hits.
        valid_hits = []
        errors, warnings = [], []
        for hit in hits:
            res = self._hit_filter.prefilter(
                query_sequence,
                hit,
                cutoff,
            )
            if res.valid:
                valid_hits.append(hit)
            if res.error:
                errors.append(res.error)
            if res.warning:
                warnings.append(res.warning)

        # Sort hits by sum_probs.
        valid_hits.sort(
            key=lambda x: x.sum_probs if x.sum_probs is not None else 0.0, reverse=True
        )

        # De-duplicate by hit sequence.
        deduped, seen_seq = [], set()
        for hit in valid_hits:
            seq = hit.hit_sequence.replace("-", "")
            if seq not in seen_seq:
                deduped.append(hit)
                seen_seq.add(seq)

        indices = list(range(len(deduped)))
        if self._shuffle_top_k_prefiltered:
            k = self._shuffle_top_k_prefiltered
            indices[:k] = np.random.permutation(indices[:k])

        if self._max_template_candidates_num:
            indices = indices[: self._max_template_candidates_num]

        # Determine number of templates to collect.
        max_to_collect = (
            self._max_hits
            if self._max_hits is not None
            else min(4, random.randint(0, min(len(indices), 20)))
        )

        features, final_hits = [], []
        already_seen_seqs = set()
        last_track = {}

        for i in indices:
            if len(features) >= max_to_collect:
                break
            hit = deduped[i]
            res, track = self._hit_processor.process(
                query_sequence,
                hit,
                cutoff,
                self._release_dates,
                self._obsolete_pdbs,
                self._strict_error_check,
            )
            last_track = track
            if res.error:
                errors.append(res.error)
            if res.warning:
                warnings.append(res.warning)
            if res.features:
                seq_key = res.features["template_sequence"]
                if seq_key not in already_seen_seqs:
                    already_seen_seqs.add(seq_key)
                    features.append(res.features)
                    final_hits.append(res.hit)

        return TemplateSearchResult(features, final_hits, errors, warnings), last_track
