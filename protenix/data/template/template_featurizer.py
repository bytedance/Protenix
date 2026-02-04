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
import time
from datetime import datetime, timedelta
from os.path import exists as opexists, join as opjoin
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
from biotite.structure import AtomArray
from typing_extensions import Self, TypeAlias

from protenix.data.constants import (
    DNA_CHAIN,
    LIGAND_CHAIN_TYPES,
    PROTEIN_CHAIN,
    RNA_CHAIN,
)
from protenix.data.msa.msa_utils import map_to_standard
from protenix.data.template.template_parser import HHRParser, HmmsearchA3MParser
from protenix.data.template.template_utils import (
    DAYS_BEFORE_QUERY_DATE,
    DistogramFeaturesConfig,
    TEMPLATE_FEATURES,
    TemplateFeatures,
    TemplateHitFeaturizer,
)
from protenix.data.utils import pad_to
from protenix.utils.file_io import load_json_cached
from protenix.utils.logger import get_logger

logger = get_logger(__name__)

BatchDict: TypeAlias = dict[str, np.ndarray]
FeatureDict: TypeAlias = Mapping[str, np.ndarray]


class TemplateSourceManager:
    """
    Manages template data retrieval and loading from multiple sources.

    Args:
        raw_paths: List of base paths for template storage.
        indexing_methods: List of indexing methods (e.g., 'sequence' or 'pdb_id').
        mappings: Dictionary mapping source index to its respective lookup table.
        enabled: Whether template loading is enabled.
    """

    def __init__(
        self,
        raw_paths: Sequence[str],
        indexing_methods: Sequence[str],
        mappings: Dict[int, Dict[str, Any]],
        enabled: bool = True,
    ) -> None:
        self.raw_paths = raw_paths
        self.indexing_methods = indexing_methods
        self.mappings = mappings
        self.enabled = enabled

    def fetch_template_paths(
        self, pdb_id: str, query_sequence: str, chain_entity_type: str
    ) -> List[str]:
        """
        Fetches template file paths from the configured sources.

        Args:
            pdb_id: PDB identifier of the query.
            query_sequence: Query sequence.
            chain_entity_type: Type of the chain (e.g., PROTEIN_CHAIN).

        Returns:
            A list of paths to found template files (.a3m or .hhr).
        """
        if not self.enabled or chain_entity_type != PROTEIN_CHAIN:
            return []

        template_paths = []
        for i, (path, method) in enumerate(zip(self.raw_paths, self.indexing_methods)):
            mapping = self.mappings.get(i, {})
            key = pdb_id if method == "pdb_id" else query_sequence

            if key not in mapping:
                continue

            dir_path = opjoin(path, str(mapping[key]))
            # Check for multiple possible template filenames
            possible_subpaths = [
                "hmmsearch.a3m",
                "concat.hhr",
            ]
            for subpath in possible_subpaths:
                full_path = opjoin(dir_path, subpath)
                if opexists(full_path):
                    template_paths.append(full_path)
        return template_paths


class TemplateFeatureAssemblyLine:
    """
    Orchestrates the conversion of raw templates into finalized Protenix features.

    Args:
        max_templates: Maximum number of templates to include in the features.
    """

    def __init__(self, max_templates: int = 4) -> None:
        self.max_templates = max_templates

    def assemble(
        self,
        bioassembly: Mapping[int, Mapping[str, Any]],
        standard_token_idxs: np.ndarray,
    ) -> "Templates":
        """
        Executes the complete feature assembly pipeline.

        Args:
            bioassembly: Mapping of asymmetric IDs to chain information.
            standard_token_idxs: Array of standardized residue indices.

        Returns:
            An assembled Templates object.
        """
        np_chains_list = []
        polymer_entity_features = {True: {}, False: {}}
        # Identify entities where template features can be safely copied (same sequence)
        safe_entity_ids = get_safe_entity_id_for_template_copy(bioassembly)

        for asym_id, info in bioassembly.items():
            chain_id = info["chain_id"]
            entity_id = info["entity_id"]
            chain_type = info["chain_entity_type"]
            num_tokens = len(info["sequence"])

            # Templates are currently only supported for protein chains with sufficient length
            skip_chain = chain_type != PROTEIN_CHAIN or num_tokens <= 4

            if (entity_id not in polymer_entity_features[skip_chain]) or (
                entity_id not in safe_entity_ids
            ):
                templates = info["templates"]
                if skip_chain or not templates:
                    template_features = TemplateFeatures.empty_template_features(
                        num_tokens
                    )
                else:
                    # Package and fix template features
                    template_features = TemplateFeatures.package_template_features(
                        hit_features=templates
                    )
                    template_features = TemplateFeatures.fix_template_features(
                        template_features=template_features,
                        num_res=num_tokens,
                    )
                # Reduce to requested maximum number of templates
                template_features = _reduce_template_features(
                    template_features, self.max_templates
                )
                if entity_id in safe_entity_ids:
                    polymer_entity_features[skip_chain][entity_id] = template_features

            if entity_id in safe_entity_ids:
                feats = polymer_entity_features[skip_chain][entity_id].copy()
            else:
                feats = template_features

            feats["chain_id"] = chain_id
            np_chains_list.append(feats)

        # Pad the number of templates to max_templates for each chain to allow concatenation
        for chain in np_chains_list:
            chain["template_aatype"] = pad_to(
                chain["template_aatype"], (self.max_templates, None)
            )
            chain["template_atom_positions"] = pad_to(
                chain["template_atom_positions"],
                (self.max_templates, None, None, None),
            )
            chain["template_atom_mask"] = pad_to(
                chain["template_atom_mask"], (self.max_templates, None, None)
            )

        # Concatenate features along the residue dimension
        merged_example = {
            ft: np.concatenate([c[ft] for c in np_chains_list], axis=1)
            for ft in np_chains_list[0]
            if ft in TEMPLATE_FEATURES
        }

        # Crop/index merged features using standard token indices
        for feature_name, v in merged_example.items():
            merged_example[feature_name] = v[
                : self.max_templates, standard_token_idxs, ...
            ]

        return Templates(
            aatype=merged_example["template_aatype"],
            atom_positions=merged_example["template_atom_positions"],
            atom_mask=merged_example["template_atom_mask"].astype(bool),
        )


class TemplateFeaturizer:
    """
    Main class for template featurization in training and evaluation.

    Args:
        dataset_name: Name of the dataset.
        prot_template_mmcif_dir: Directory containing template mmCIF files.
        prot_template_cache_dir: Directory for cached template data.
        prot_template_raw_paths: Base paths for protein template storage.
        prot_seq_or_filename_to_templatedir_jsons: JSON files mapping sequences/PDB IDs to template directories.
        prot_indexing_methods: Indexing methods for each source.
        enable_prot_template: Whether to enable protein templates.
        template_dropout_rate: Dropout rate for templates during training.
        stage: 'train' or 'eval'.
        kalign_binary_path: Path to kalign binary.
        release_dates_path: Path to PDB release dates JSON.
        obsolete_pdbs_path: Path to obsolete PDBs mapping JSON.
        **kwargs: Additional arguments.
    """

    def __init__(
        self,
        dataset_name: str = "",
        prot_template_mmcif_dir: str = "",
        prot_template_cache_dir: str = None,
        prot_template_raw_paths: Sequence[str] = [""],
        prot_seq_or_filename_to_templatedir_jsons: Sequence[str] = [""],
        prot_indexing_methods: Sequence[str] = ["sequence"],
        enable_prot_template: bool = True,
        template_dropout_rate: float = 0.0,
        stage: str = "train",
        kalign_binary_path: str = "",
        release_dates_path: str = "",
        obsolete_pdbs_path: str = "",
        **kwargs,
    ) -> None:
        # Validate inputs
        assert len(prot_template_raw_paths) == len(
            prot_seq_or_filename_to_templatedir_jsons
        ), "Mismatch between raw paths and JSON mapping files."
        assert len(kalign_binary_path) > 0, "kalign_binary_path must be provided."
        assert len(release_dates_path) > 0, "release_dates_path must be provided."
        assert len(obsolete_pdbs_path) > 0, "obsolete_pdbs_path must be provided."
        assert len(prot_template_raw_paths) == len(
            prot_indexing_methods
        ), "Mismatch between raw paths and indexing methods."

        # Currently, multiple source training strategy is not determined
        assert (
            len(prot_template_raw_paths) == 1
        ), "Only single source template supported."

        if stage != "train":
            max_template_date = "2021-09-30"
        else:
            if ("distillation" in dataset_name) or ("openprotein" in dataset_name):
                max_template_date = "2018-04-30"
            else:
                max_template_date = "3000-01-01"

        self.stage = stage
        self.max_template_date = datetime.strptime(max_template_date, "%Y-%m-%d")
        self.dataset_name = dataset_name
        self.enable_prot_template = enable_prot_template
        self.template_dropout_rate = template_dropout_rate

        # Initialize source manager
        mappings = {
            i: load_json_cached(p)
            for i, p in enumerate(prot_seq_or_filename_to_templatedir_jsons)
        }
        self.source_mgr = TemplateSourceManager(
            raw_paths=prot_template_raw_paths,
            indexing_methods=prot_indexing_methods,
            mappings=mappings,
            enabled=enable_prot_template,
        )

        # Initialize online featurizer
        self.online_featurizer = TemplateHitFeaturizer(
            mmcif_dir=prot_template_mmcif_dir,
            template_cache_dir=prot_template_cache_dir,
            max_hits=None if self.stage == "train" else 4,
            kalign_binary_path=kalign_binary_path,
            max_template_date=max_template_date,
            release_dates_path=release_dates_path,
            obsolete_pdbs_path=obsolete_pdbs_path,
            _shuffle_top_k_prefiltered=20 if self.stage == "train" else None,
        )

        self._last_profile: Dict[str, Any] = {
            "elapsed_seconds": 0.0,
            "raw_templates": 0,
            "selected_templates": 0,
            "chain_profiles": [],
        }

    def get_template(
        self,
        pdb_id: str,
        query_sequence: str,
        sequence_uid: str,
        query_release_date: Optional[datetime],
        chain_entity_type: str,
    ) -> Tuple[List[Dict[str, Any]], int, Dict[str, float]]:
        """
        Fetches and processes templates for a single sequence.

        Args:
            pdb_id: PDB identifier.
            query_sequence: Query sequence string.
            sequence_uid: Unique identifier for the sequence.
            query_release_date: Release date of the query structure.
            chain_entity_type: Type of the chain.

        Returns:
            A tuple of (template_features_list, raw_hit_count, timing_track).
        """
        if not self.enable_prot_template or chain_entity_type != PROTEIN_CHAIN:
            return [], 0, {}

        # Fetch template paths from sources
        template_paths = self.source_mgr.fetch_template_paths(
            pdb_id, query_sequence, chain_entity_type
        )

        hmmsearched_a3m, hmmsearched_hhr = "", ""
        for path in template_paths:
            with open(path, "r") as f:
                content = f.read()
            if not content:
                continue
            if path.endswith(".a3m"):
                hmmsearched_a3m += content
            elif path.endswith(".hhr"):
                hmmsearched_hhr += content

        # Determine cutoff date for filtering templates
        if query_release_date is None:
            cutoff_date = self.max_template_date
        else:
            cutoff_date = min(
                self.max_template_date,
                query_release_date - timedelta(days=DAYS_BEFORE_QUERY_DATE),
            )

        # Parse hits from a3m or hhr content
        if hmmsearched_hhr:
            hits = HHRParser.parse(hhr_string=hmmsearched_hhr)
        else:
            hits = HmmsearchA3MParser.parse(
                query_seq=query_sequence, a3m_str=hmmsearched_a3m, skip_first=False
            )

        raw_hit_count = len(hits)
        if raw_hit_count == 0:
            return [], 0, {}

        # Process hits into features using the online featurizer
        result, time_track = self.online_featurizer.get_templates(
            sequence_uid=sequence_uid,
            query_sequence=query_sequence,
            hits=hits,
            max_template_date=cutoff_date,
        )
        return result.features, raw_hit_count, time_track

    def make_template_features(
        self,
        bioassembly_dict: Dict[str, Any],
        selected_token_indices: Optional[np.ndarray],
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
    ) -> Dict[str, np.ndarray]:
        """
        Main method to generate template features for a bioassembly.

        Args:
            bioassembly_dict: Dictionary containing bioassembly data.
            selected_token_indices: Indices of selected tokens (for cropping).
            entity_to_asym_id_int: Mapping from entity ID to asymmetric IDs.

        Returns:
            Dictionary of processed template features.
        """
        profile_start = time.time()
        atom_array = bioassembly_dict["atom_array"]
        token_array = bioassembly_dict["token_array"]

        # Identify selected asymmetric IDs
        centre_atom_indices = (
            token_array[selected_token_indices]
            if selected_token_indices is not None
            else token_array
        ).get_annotation("centre_atom_index")
        selected_asym_ids = set(atom_array[centre_atom_indices].asym_id_int)

        # Build mapping from asym_id to entity_id and sequence
        asym_to_entity_id: Dict[int, str] = {}
        for eid, asyms in entity_to_asym_id_int.items():
            for aid in asyms:
                if aid in selected_asym_ids:
                    asym_to_entity_id[aid] = eid

        entity_id_to_sequence = {
            eid: bioassembly_dict["sequences"].get(eid)
            for eid in asym_to_entity_id.values()
        }

        # Patch sequences for ligands if missing
        asym_id_to_sequence = {}
        for aid, eid in asym_to_entity_id.items():
            seq = entity_id_to_sequence.get(eid)
            if seq is None:
                ligand_mask = atom_array.asym_id_int == aid
                assert atom_array.is_ligand[ligand_mask].all(), "Unknown molecule type"
                seq = "X" * ligand_mask.sum()
            asym_id_to_sequence[aid] = seq

        poly_types_mapping = {
            "polypeptide(L)": PROTEIN_CHAIN,
            "polyribonucleotide": RNA_CHAIN,
            "polydeoxyribonucleotide": DNA_CHAIN,
        }
        entity_poly_type = bioassembly_dict.get("entity_poly_type", {})

        query_release_date = None
        if self.stage == "train" and "release_date" in bioassembly_dict:
            query_release_date = datetime.strptime(
                bioassembly_dict["release_date"], "%Y-%m-%d"
            )

        # Global template dropout
        drop_template = False
        if self.stage == "train" and self.template_dropout_rate > 0.0:
            drop_template = np.random.uniform(0, 1) < self.template_dropout_rate

        template_meta_infos = {}
        seq_to_templates_cache: Dict[str, Tuple[List[Dict], int, Dict]] = {}
        profile = {
            "elapsed_seconds": 0.0,
            "raw_templates": 0,
            "selected_templates": 0,
            "chain_profiles": [],
        }

        pdb_id = bioassembly_dict.get("pdb_id", "unknown")

        for aid, sequence in asym_id_to_sequence.items():
            chain_start = time.time()
            eid = asym_to_entity_id[aid]
            chain_type = poly_types_mapping.get(
                entity_poly_type.get(eid, "non-polymer"), LIGAND_CHAIN_TYPES
            )

            if drop_template or chain_type != PROTEIN_CHAIN:
                templates, raw_hit_count, time_track = [], 0, {}
            else:
                if sequence in seq_to_templates_cache:
                    templates, raw_hit_count, time_track = seq_to_templates_cache[
                        sequence
                    ]
                else:
                    templates, raw_hit_count, time_track = self.get_template(
                        pdb_id=pdb_id,
                        query_sequence=sequence,
                        sequence_uid=f"{pdb_id}_{aid}",
                        query_release_date=query_release_date,
                        chain_entity_type=chain_type,
                    )
                    seq_to_templates_cache[sequence] = (
                        templates,
                        raw_hit_count,
                        time_track,
                    )

            chain_id = atom_array.chain_id[atom_array.asym_id_int == aid][0]
            if chain_type == PROTEIN_CHAIN:
                num_sel = len(templates) if isinstance(templates, list) else 0
                profile["raw_templates"] += raw_hit_count
                profile["selected_templates"] += num_sel
                profile["chain_profiles"].append(
                    {
                        "asym_id": aid,
                        "chain_id": chain_id,
                        "elapsed_seconds": time.time() - chain_start,
                        "raw_templates": raw_hit_count,
                        "selected_templates": num_sel,
                        "chain_entity_type": chain_type,
                        "time_track": time_track,
                    }
                )

            template_meta_infos[aid] = {
                "entity_id": eid,
                "chain_id": chain_id,
                "sequence": sequence,
                "chain_entity_type": chain_type,
                "templates": templates,
            }

        if not get_safe_entity_id_for_template_copy(template_meta_infos):
            logger.warning(
                f"PDB {pdb_id} has inconsistent sequences within same entity."
            )

        # Coordinate mapping to standardized indices
        ca = atom_array[centre_atom_indices]
        std_idxs = map_to_standard(ca.asym_id_int, ca.res_id, template_meta_infos)

        # Assemble final features
        max_t = 4 if self.stage != "train" else 4  # Default to 4
        assembly_line = TemplateFeatureAssemblyLine(max_templates=max_t)
        template_features = assembly_line.assemble(
            template_meta_infos, std_idxs
        ).as_protenix_dict()

        profile["elapsed_seconds"] = time.time() - profile_start
        self._last_profile = profile
        return template_features

    def __call__(
        self,
        bioassembly_dict: Dict[str, Any],
        selected_indices: np.ndarray,
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
    ) -> Dict[str, np.ndarray]:
        """Call method proxying to make_template_features."""
        return self.make_template_features(
            bioassembly_dict, selected_indices, entity_to_asym_id_int
        )

    def get_last_profile(self) -> Dict[str, Union[float, int]]:
        """Returns profiling information for the last call."""
        return self._last_profile


@dataclasses.dataclass(frozen=True)
class Templates:
    """Dataclass containing template features."""

    # aatype: [num_templates, num_res]
    aatype: np.ndarray
    # atom_positions: [num_templates, num_res, 24, 3]
    atom_positions: np.ndarray
    # atom_mask: [num_templates, num_res, 24]
    atom_mask: np.ndarray

    @classmethod
    def from_data_dict(cls, batch: BatchDict) -> Self:
        """Construct instance from a data dictionary."""
        return cls(
            aatype=batch["template_aatype"],
            atom_positions=batch["template_atom_positions"],
            atom_mask=batch["template_atom_mask"],
        )

    def as_data_dict(self) -> BatchDict:
        """Convert to a standard data dictionary."""
        return {
            "template_aatype": self.aatype,
            "template_atom_positions": self.atom_positions,
            "template_atom_mask": self.atom_mask,
        }

    def as_protenix_dict(self) -> BatchDict:
        """Compute additional features and return as Protenix dictionary."""
        features = self.as_data_dict()
        dgrams, pb_masks = [], []
        unit_vectors, bb_masks = [], []

        num_templates = self.aatype.shape[0]
        for i in range(num_templates):
            aatype = self.aatype[i]
            mask = self.atom_mask[i]
            pos = self.atom_positions[i] * mask[..., None]

            # Compute pseudo-beta positions and mask
            pb_pos, pb_mask = TemplateFeatures.pseudo_beta_fn(aatype, pos, mask)
            pb_mask_2d = pb_mask[:, None] * pb_mask[None, :]

            # Compute distogram
            dgram = TemplateFeatures.dgram_from_positions(
                pb_pos,
                config=DistogramFeaturesConfig(
                    min_bin=3.25, max_bin=50.75, num_bins=39
                ),
            )
            dgrams.append(dgram * pb_mask_2d[..., None])
            pb_masks.append(pb_mask_2d)

            # Compute normalized unit vectors between residues
            uv, bb_mask_2d = TemplateFeatures.compute_template_unit_vector(
                aatype, pos, mask
            )
            unit_vectors.append(uv * bb_mask_2d[..., None])
            bb_masks.append(bb_mask_2d)

        features.update(
            {
                "template_pseudo_beta_mask": np.stack(pb_masks),
                "template_distogram": np.stack(dgrams),
                "template_unit_vector": np.stack(unit_vectors),
                "template_backbone_frame_mask": np.stack(bb_masks),
            }
        )
        return features


def _reduce_template_features(
    template_features: FeatureDict, max_templates: int
) -> FeatureDict:
    """Reduces templates to the requested maximum number."""
    num_t = template_features["template_aatype"].shape[0]
    keep_mask = np.arange(num_t) < max_templates
    fields = TEMPLATE_FEATURES + ("template_release_timestamp",)
    return {k: v[keep_mask] for k, v in template_features.items() if k in fields}


def get_safe_entity_id_for_template_copy(
    bioassembly: Mapping[int, Mapping[str, Any]],
) -> List[str]:
    """Identifies entity IDs that have consistent sequences across all chains."""
    eid_to_seqs = {}
    for aid, info in bioassembly.items():
        eid = info["entity_id"]
        seq = info["sequence"]
        eid_to_seqs.setdefault(eid, set()).add(seq)
    return [eid for eid, seqs in eid_to_seqs.items() if len(seqs) == 1]


class InferenceTemplateFeaturizer:
    """Simplified featurizer for inference, leveraging the same assembly logic."""

    @staticmethod
    def make_template_feature(
        bioassembly: Sequence[Mapping[str, Any]],
        atom_array: AtomArray,
        use_template: bool = True,
        online_template_featurizer: Optional[TemplateHitFeaturizer] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Generates template features during inference.

        Args:
            bioassembly: List of entity information from the input JSON.
            atom_array: Parsed atom structure.
            use_template: Whether to use templates.
            online_template_featurizer: Featurizer for processing template hits.

        Returns:
            Dictionary of template features.
        """
        logger.info("Calling InferenceTemplateFeaturizer.make_template_feature")
        template_meta_infos = {}
        curr_asym_id = 0

        for eid, info in enumerate(bioassembly):
            seq, count, ctype, t_path = "", 0, LIGAND_CHAIN_TYPES, ""

            if "proteinChain" in info:
                c = info["proteinChain"]
                seq, count, ctype, t_path = (
                    c["sequence"],
                    c["count"],
                    PROTEIN_CHAIN,
                    c.get("templatesPath", ""),
                )
            elif "rnaSequence" in info:
                c = info["rnaSequence"]
                seq, count, ctype = c["sequence"], c["count"], RNA_CHAIN
            elif "dnaSequence" in info:
                c = info["dnaSequence"]
                seq, count, ctype = c["sequence"], c["count"], DNA_CHAIN
            elif "ligand" in info:
                count, ctype = info["ligand"]["count"], LIGAND_CHAIN_TYPES
                seq = "X" * (atom_array.asym_id_int == curr_asym_id).sum()

            templates = []
            if t_path and use_template and online_template_featurizer:
                assert ctype == PROTEIN_CHAIN, "Only protein templates are supported."
                with open(t_path, "r") as f:
                    content = f.read()

                if t_path.endswith(".hhr"):
                    hits = HHRParser.parse(hhr_string=content)
                elif t_path.endswith(".a3m"):
                    hits = HmmsearchA3MParser.parse(
                        query_seq=seq, a3m_str=content, skip_first=False
                    )
                else:
                    raise ValueError(f"Unsupported template format: {t_path}")

                result, _ = online_template_featurizer.get_templates(
                    sequence_uid=seq,
                    query_sequence=seq,
                    hits=hits,
                    max_template_date=None,
                )
                templates = result.features
                logger.info(f"Found {len(templates)} templates for sequence {seq}")

            for i in range(count):
                aid = curr_asym_id + i
                template_meta_infos[aid] = {
                    "entity_id": eid,
                    "chain_id": atom_array.chain_id[atom_array.asym_id_int == aid][0],
                    "sequence": seq,
                    "chain_entity_type": ctype,
                    "templates": templates,
                }
            curr_asym_id += count

        # Coordinate mapping
        ca = atom_array[atom_array.centre_atom_mask.astype(bool)]
        std_idxs = map_to_standard(ca.asym_id_int, ca.res_id, template_meta_infos)

        # Assemble features
        return (
            TemplateFeatureAssemblyLine(max_templates=4)
            .assemble(template_meta_infos, std_idxs)
            .as_protenix_dict()
        )
