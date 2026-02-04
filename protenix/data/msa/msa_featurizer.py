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
from os.path import exists as opexists, join as opjoin
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from biotite.structure import AtomArray

from protenix.data.constants import (
    DNA_CHAIN,
    LIGAND_CHAIN_TYPES,
    PROTEIN_CHAIN,
    RNA_CHAIN,
    STANDARD_POLYMER_CHAIN_TYPES,
    STD_RESIDUES_WITH_GAP,
)
from protenix.data.msa.msa_utils import (
    map_to_standard,
    MSA_GAP_IDX,
    MSAPairingEngine,
    NUM_SEQ_NUM_RES_MSA_FEATURES,
    RawMsa,
)
from protenix.utils.file_io import load_json_cached
from protenix.utils.logger import get_logger

logger = get_logger(__name__)


class MSASourceManager:
    """
    Manages MSA data retrieval and loading from multiple sources.

    Args:
        raw_paths: List of base paths for MSA storage.
        methods: List of indexing methods (e.g., 'sequence' or 'pdb_id').
        mappings: Dictionary mapping source index to its respective lookup table.
        enabled: Whether MSA loading is enabled.
    """

    def __init__(
        self,
        raw_paths: Sequence[str],
        methods: Sequence[str],
        mappings: Dict[int, Dict[Any, Any]],
        enabled: bool,
    ) -> None:
        self.raw_paths = raw_paths
        self.methods = methods
        self.mappings = mappings
        self.enabled = enabled

    def fetch_msas(
        self,
        sequence: str,
        pdb_id: str,
        chain_type: str,
        p_dbs: Optional[Sequence[str]] = None,
        np_dbs: Optional[Sequence[str]] = None,
    ) -> Tuple[List[RawMsa], List[RawMsa]]:
        """
        Fetches MSAs from the configured paths based on chain type and indexing methods.

        Args:
            sequence: Query sequence.
            pdb_id: PDB identifier.
            chain_type: Type of the chain (e.g., PROTEIN_CHAIN, RNA_CHAIN).
            p_dbs: Databases for paired MSA search.
            np_dbs: Databases for unpaired MSA search.

        Returns:
            A tuple of (unpaired_msas, paired_msas).
        """
        if not self.enabled:
            return [], []
        unpaired, paired = [], []

        # RNA-specific loading logic
        if chain_type == RNA_CHAIN:
            for path, method, m_key in zip(self.raw_paths, self.methods, self.mappings):
                if method != "sequence":
                    continue
                mapping = self.mappings[m_key]
                if sequence not in mapping:
                    continue
                eid = str(mapping[sequence][0])
                fpath = opjoin(path, eid, f"{eid}_all.a3m")
                if opexists(fpath):
                    with open(fpath, "r") as f:
                        content = f.read()
                    if content:
                        unpaired.append(
                            RawMsa.from_a3m(
                                sequence,
                                RNA_CHAIN,
                                content,
                                depth_limit=30000,
                                dedup=False,
                            )
                        )
            return unpaired, []

        # Protein-specific loading logic
        if chain_type == PROTEIN_CHAIN:
            p_dbs = p_dbs or []
            np_dbs = np_dbs or []
            for p_db, np_db, path, method, m_key in zip(
                p_dbs, np_dbs, self.raw_paths, self.methods, self.mappings
            ):
                key = sequence if method == "sequence" else str(pdb_id)
                mapping = self.mappings[m_key]
                if key not in mapping:
                    continue
                dir_path = opjoin(path, str(mapping[key]))

                if p_db:
                    for pat in [f"{p_db}.a3m"]:
                        fpath = opjoin(dir_path, pat)
                        if opexists(fpath):
                            with open(fpath, "r") as f:
                                content = f.read()
                            if content:
                                paired.append(
                                    RawMsa.from_a3m(
                                        sequence, PROTEIN_CHAIN, content, dedup=False
                                    )
                                )
                            break
                if np_db:
                    for db in np_db.split("-"):
                        for pat in [f"{db}.a3m"]:
                            fpath = opjoin(dir_path, pat)
                            if opexists(fpath):
                                with open(fpath, "r") as f:
                                    content = f.read()
                                if content:
                                    unpaired.append(
                                        RawMsa.from_a3m(
                                            sequence,
                                            PROTEIN_CHAIN,
                                            content,
                                            dedup=False,
                                        )
                                    )
                                break
        return unpaired, paired


class FeatureAssemblyLine:
    """
    Orchestrates the conversion of Raw MSAs into finalized Protenix features.

    Args:
        max_msa_size: Maximum number of sequences allowed in the final MSA.
        max_paired_per_species: Maximum number of paired sequences per species.
    """

    def __init__(
        self, max_msa_size: int = 16384, max_paired_per_species: int = 600
    ) -> None:
        self.max_size = max_msa_size
        self.max_paired_per_sp = max_paired_per_species

    def assemble(
        self, bioassembly: Mapping[int, Mapping[str, Any]], std_idxs: np.ndarray
    ) -> "MSAFeat":
        """
        Executes the complete feature assembly pipeline.

        Args:
            bioassembly: Mapping of asymmetric IDs to chain information.
            std_idxs: Array of standardized residue indices.

        Returns:
            An assembled MSAFeat object.
        """
        # 1. Base featurization
        unique_prot_seqs = {
            v["sequence"]
            for v in bioassembly.values()
            if v["chain_entity_type"] == PROTEIN_CHAIN
        }
        need_pairing = len(unique_prot_seqs) > 1
        active_chain_ids = {v["chain_id"] for v in bioassembly.values()}

        raw_chains = []
        for aid, info in bioassembly.items():
            ctype, seq = info["chain_entity_type"], info["sequence"]
            skip = ctype not in STANDARD_POLYMER_CHAIN_TYPES or len(seq) <= 4

            if ctype in STANDARD_POLYMER_CHAIN_TYPES:
                up_msa = RawMsa.from_a3m(
                    seq,
                    ctype,
                    (
                        info["unpaired_msa"]
                        if not skip and ctype in [PROTEIN_CHAIN, RNA_CHAIN]
                        else ""
                    ),
                    dedup=True,
                )
                p_msa = RawMsa.from_a3m(
                    seq,
                    ctype,
                    (
                        info["paired_msa"]
                        if not skip and need_pairing and ctype == PROTEIN_CHAIN
                        else ""
                    ),
                    dedup=False,
                )
            else:
                up_msa = p_msa = RawMsa(
                    seq, PROTEIN_CHAIN, [], [], deduplicate=False
                )  # Ligand placeholders

            u_f, p_f = up_msa.featurize(), p_msa.featurize()
            chain_feat = dict(u_f)
            chain_feat.update({f"{k}_all_seq": v for k, v in p_f.items()})
            chain_feat.update(
                {
                    "asym_id": np.full(len(seq), aid),
                    "chain_id": info["chain_id"],
                    "entity_id": info["entity_id"],
                }
            )
            # Compute Profile
            msa = chain_feat["msa"]
            prof = (msa[..., None] == np.arange(len(STD_RESIDUES_WITH_GAP))).sum(
                axis=0
            ) / msa.shape[0]
            chain_feat.update(
                {
                    "profile": prof.astype(np.float32),
                    "deletion_mean": np.mean(chain_feat["deletion_matrix"], axis=0),
                }
            )
            raw_chains.append(chain_feat)

        # 2. Pairing and cleanup
        max_p = self.max_size // 2
        if need_pairing:
            raw_chains = MSAPairingEngine.pair_chains_by_species(
                raw_chains, max_p, active_chain_ids, self.max_paired_per_sp
            )
            raw_chains = MSAPairingEngine.cleanup_unpaired_features(raw_chains)

        # 3. Filter all-gap rows
        nonempty_asyms = [
            c["asym_id"][0] for c in raw_chains if c["chain_id"] in active_chain_ids
        ]
        if "msa_all_seq" in raw_chains[0]:
            raw_chains = MSAPairingEngine.filter_all_gapped_rows(
                raw_chains, nonempty_asyms
            )

        # 4. Cropping and merging
        cropped = []
        for c in raw_chains:
            p_msa = c.get("msa_all_seq")
            ps = min(p_msa.shape[0], max_p) if p_msa is not None else 0
            us = max(0, min(c["msa"].shape[0], self.max_size - ps))

            cr = {
                "asym_id": c["asym_id"],
                "chain_id": c["chain_id"],
                "profile": c["profile"],
                "deletion_mean": c["deletion_mean"],
            }
            for k in NUM_SEQ_NUM_RES_MSA_FEATURES:
                if k in c:
                    cr[k] = c[k][:us]
                if f"{k}_all_seq" in c:
                    cr[f"{k}_all_seq"] = c[f"{k}_all_seq"][:ps]
            cropped.append(cr)

        merged = {"asym_id": np.concatenate([c["asym_id"] for c in cropped])}
        for base in NUM_SEQ_NUM_RES_MSA_FEATURES:
            for f in [base, f"{base}_all_seq"]:
                if f in cropped[0]:
                    merged[f] = MSAPairingEngine.merge_chain_features(cropped, f)
        for f in ["profile", "deletion_mean"]:
            merged[f] = np.concatenate([c[f] for c in cropped])

        # 5. Depth tracking
        active_set = set(nonempty_asyms)
        max_u = max([len(c["msa"]) for c in cropped if c["asym_id"][0] in active_set])
        rna_u = max(
            [1]
            + [
                len(c["msa"])
                for c in cropped
                if bioassembly[c["asym_id"][0]]["chain_entity_type"] == RNA_CHAIN
            ]
        )
        prot_u = max(
            [1]
            + [
                len(c["msa"])
                for c in cropped
                if bioassembly[c["asym_id"][0]]["chain_entity_type"] == PROTEIN_CHAIN
            ]
        )

        merged["msa"] = merged["msa"][:max_u]
        prot_p = 1
        if "msa_all_seq" in merged:
            max_p_actual = max(
                [
                    len(c["msa_all_seq"])
                    for c in cropped
                    if c["asym_id"][0] in active_set
                ]
            )
            merged["msa_all_seq"] = merged["msa_all_seq"][:max_p_actual]
            prot_p = max_p_actual

        # 6. Final integration and coordinate mapping
        for k in NUM_SEQ_NUM_RES_MSA_FEATURES:
            if k in merged and f"{k}_all_seq" in merged:
                merged[k] = np.concatenate([merged[f"{k}_all_seq"], merged[k]], axis=0)

        # Forward compatibility patch for non-protein entities
        for aid in [
            aid
            for aid, info in bioassembly.items()
            if info["chain_entity_type"] != PROTEIN_CHAIN
        ]:
            cols = np.where(merged["asym_id"] == aid)[0]
            if cols.size > 0:
                gap_mask = np.all(merged["msa"][:, cols] == MSA_GAP_IDX, axis=1)
                merged["msa"][np.ix_(np.where(gap_mask)[0], cols)] = merged["msa"][
                    0, cols
                ]

        for f in NUM_SEQ_NUM_RES_MSA_FEATURES:
            if f in merged:
                merged[f] = merged[f][:, std_idxs].copy()
        for f in ["profile", "deletion_mean"]:
            merged[f] = merged[f][std_idxs]

        def to_i8(x: np.ndarray) -> np.ndarray:
            return np.clip(x, -128, 127).astype(np.int8)

        return MSAFeat(
            rows=to_i8(merged["msa"]),
            mask=np.ones_like(merged["msa"], dtype=bool),
            deletion_matrix=to_i8(merged["deletion_matrix"]),
            profile=merged["profile"],
            deletion_mean=merged["deletion_mean"],
            prot_unpaired_num_alignments=np.array(prot_u, dtype=np.int32),
            prot_paired_num_alignments=np.array(prot_p, dtype=np.int32),
            rna_unpaired_num_alignments=np.array(rna_u, dtype=np.int32),
        )


class MSAFeaturizer:
    """
    Main entry point for MSA featurization, coordinating source management and assembly.

    Args:
        dataset_name: Name of the dataset.
        prot_seq_or_filename_to_msadir_jsons: JSON maps for protein MSA lookups.
        prot_msadir_raw_paths: Base paths for protein MSAs.
        rna_seq_or_filename_to_msadir_jsons: JSON maps for RNA MSA lookups.
        rna_msadir_raw_paths: Base paths for RNA MSAs.
        prot_pairing_dbs: List of databases for protein pairing.
        prot_non_pairing_dbs: List of databases for protein non-pairing.
        prot_indexing_methods: Methods for protein MSA indexing.
        rna_indexing_methods: Methods for RNA MSA indexing.
        enable_prot_msa: Whether to enable protein MSA processing.
        enable_rna_msa: Whether to enable RNA MSA processing.
    """

    def __init__(
        self,
        dataset_name: str = "",
        prot_seq_or_filename_to_msadir_jsons: Sequence[str] = [""],
        prot_msadir_raw_paths: Sequence[str] = [""],
        rna_seq_or_filename_to_msadir_jsons: Sequence[str] = [""],
        rna_msadir_raw_paths: Sequence[str] = [""],
        prot_pairing_dbs: Sequence[str] = [""],
        prot_non_pairing_dbs: Sequence[str] = [""],
        prot_indexing_methods: Sequence[str] = ["sequence"],
        rna_indexing_methods: Sequence[str] = ["sequence"],
        enable_prot_msa: bool = True,
        enable_rna_msa: bool = True,
    ) -> None:
        self.dataset_name = dataset_name
        super().__init__()
        # Initialize source managers for protein and RNA
        self.prot_mgr = MSASourceManager(
            prot_msadir_raw_paths,
            prot_indexing_methods,
            {
                i: load_json_cached(p)
                for i, p in enumerate(prot_seq_or_filename_to_msadir_jsons)
            },
            enable_prot_msa,
        )
        self.rna_mgr = MSASourceManager(
            rna_msadir_raw_paths,
            rna_indexing_methods,
            {
                i: load_json_cached(p)
                for i, p in enumerate(rna_seq_or_filename_to_msadir_jsons)
            },
            enable_rna_msa,
        )
        self.prot_p_dbs = prot_pairing_dbs
        self.prot_np_dbs = prot_non_pairing_dbs
        self._profile: Dict[str, Any] = {}
        logger.info(f"MSAFeaturizer for {dataset_name} initialized.")

    def set_last_profile(self, p: Dict[str, Any]) -> None:
        """Sets the internal profile state."""
        self._profile = p

    def get_last_profile(self) -> Dict[str, Any]:
        """Returns the internal profile state."""
        return self._profile

    def make_msa_features(
        self,
        bioassembly_dict: Dict[str, Any],
        selected_indices: Optional[np.ndarray],
        entity_to_asym_id_int: Mapping[str, Sequence[int]],
    ) -> Dict[str, Any]:
        """
        Processes bioassembly information into a dictionary of MSA features.

        Args:
            bioassembly_dict: Dictionary containing biological assembly data.
            selected_indices: Optional array of indices to select from the token array.
            entity_to_asym_id_int: Mapping from entity ID to asymmetric IDs.

        Returns:
            A dictionary containing processed MSA features.
        """
        atom_array, token_array = (
            bioassembly_dict["atom_array"],
            bioassembly_dict["token_array"],
        )
        sel_tokens = (
            token_array[selected_indices]
            if selected_indices is not None
            else token_array
        )
        sel_asyms = set(
            atom_array[sel_tokens.get_annotation("centre_atom_index")].asym_id_int
        )

        # 1. Resolve metadata and fetch MSAs
        meta = {}
        poly_map = {
            "polypeptide(L)": PROTEIN_CHAIN,
            "polyribonucleotide": RNA_CHAIN,
            "polydeoxyribonucleotide": DNA_CHAIN,
        }
        for eid, asyms in entity_to_asym_id_int.items():
            for aid in [a for a in asyms if a in sel_asyms]:
                seq = bioassembly_dict["sequences"].get(eid) or (
                    "X" * (atom_array.asym_id_int == aid).sum()
                )
                ctype = poly_map.get(
                    bioassembly_dict["entity_poly_type"].get(eid, "non-polymer"),
                    LIGAND_CHAIN_TYPES,
                )

                up_msas, p_msas = [], []
                if ctype == RNA_CHAIN:
                    up_msas, _ = self.rna_mgr.fetch_msas(seq, "", ctype)
                elif ctype == PROTEIN_CHAIN:
                    up_msas, p_msas = self.prot_mgr.fetch_msas(
                        seq,
                        bioassembly_dict["pdb_id"],
                        ctype,
                        self.prot_p_dbs,
                        self.prot_np_dbs,
                    )

                meta[aid] = {
                    "entity_id": eid,
                    "chain_id": atom_array.chain_id[atom_array.asym_id_int == aid][0],
                    "sequence": seq,
                    "chain_entity_type": ctype,
                    "paired_msa": RawMsa.merge(p_msas).to_a3m() if p_msas else "",
                    "unpaired_msa": RawMsa.merge(up_msas).to_a3m() if up_msas else "",
                }

        # 2. Map coordinates and assemble features
        ca = atom_array[sel_tokens.get_annotation("centre_atom_index")]
        std_idxs = map_to_standard(ca.asym_id_int, ca.res_id, meta)

        res = FeatureAssemblyLine().assemble(meta, std_idxs).to_dict()
        keep = {
            "msa",
            "has_deletion",
            "deletion_value",
            "deletion_mean",
            "profile",
            "prot_pair_num_alignments",
            "prot_unpair_num_alignments",
            "rna_pair_num_alignments",
            "rna_unpair_num_alignments",
        }
        return {k: v for k, v in res.items() if k in keep}

    def __call__(self, *args: Any, **kwargs: Any) -> Dict[str, Any]:
        """Convenience method to call make_msa_features."""
        return self.make_msa_features(*args, **kwargs)


@dataclasses.dataclass(frozen=True)
class MSAFeat:
    """Container for finalized numerical MSA features."""

    rows: np.ndarray
    mask: np.ndarray
    deletion_matrix: np.ndarray
    profile: np.ndarray
    deletion_mean: np.ndarray
    prot_unpaired_num_alignments: np.ndarray
    prot_paired_num_alignments: np.ndarray
    rna_unpaired_num_alignments: np.ndarray

    def to_dict(self) -> Dict[str, Any]:
        """Converts the MSA object into a standard Protenix data dictionary."""
        return {
            "msa": self.rows,
            "msa_mask": self.mask,
            "deletion_matrix": self.deletion_matrix,
            "deletion_value": (np.arctan(self.deletion_matrix / 3.0) * (2.0 / np.pi)),
            "has_deletion": np.clip(self.deletion_matrix, 0.0, 1.0),
            "profile": self.profile,
            "deletion_mean": self.deletion_mean,
            "prot_unpaired_num_alignments": self.prot_unpaired_num_alignments,
            "prot_paired_num_alignments": self.prot_paired_num_alignments,
            "rna_unpaired_num_alignments": self.rna_unpaired_num_alignments,
            "rna_pair_num_alignments": np.asarray(1, dtype=np.int32),
            "prot_pair_num_alignments": self.prot_paired_num_alignments,
            "prot_unpair_num_alignments": self.prot_unpaired_num_alignments,
            "rna_unpair_num_alignments": self.rna_unpaired_num_alignments,
        }


def ensure_ends_with_newline(s: Optional[str]) -> Optional[str]:
    """
    Ensure the given string ends with a newline character.

    If the string is non-empty and does not already end with '\\n',
    append '\\n'. Empty strings are returned unchanged.

    Args:
        s (str): Input string.

    Returns:
        str: The input string guaranteed to end with '\\n' when non-empty.
    """
    if not s:
        return s
    if not s.endswith("\n"):
        s += "\n"
    return s


class InferenceMSAFeaturizer:
    """Specialized featurizer for inference scenarios, leveraging the unified assembly line."""

    @staticmethod
    def make_msa_feature(
        bioassembly: Sequence[Dict[str, Any]],
        atom_array: AtomArray,
        msa_pair_as_unpair: bool = False,
        use_rna_msa: bool = True,
    ) -> Dict[str, Any]:
        """
        Prepares MSA features during inference from bioassembly structure.

        Args:
            bioassembly: List of entities in the biological assembly.
            atom_array: Structural data array.
            msa_pair_as_unpair: Whether to treat paired MSA as unpaired.
            use_rna_msa: Whether to use MSA for RNA chains.

        Returns:
            Dictionary of processed MSA features.
        """
        meta, curr_aid = {}, 0
        for eid, info in enumerate(bioassembly):
            seq, count, ctype, u_a3m, p_a3m = "", 0, LIGAND_CHAIN_TYPES, None, None
            if "proteinChain" in info:
                c = info["proteinChain"]
                seq, count, ctype, u_a3m, p_a3m = (
                    c["sequence"],
                    c["count"],
                    PROTEIN_CHAIN,
                    c.get("unpairedMsa"),
                    c.get("pairedMsa"),
                )
                if u_a3m is None and c.get("unpairedMsaPath"):
                    with open(c["unpairedMsaPath"]) as f:
                        u_a3m = f.read()
                if p_a3m is None and c.get("pairedMsaPath"):
                    with open(c["pairedMsaPath"]) as f:
                        p_a3m = f.read()
                if u_a3m is None and (p_a3m is None):
                    if c.get("msa"):
                        msa_dir = c["msa"].get("precomputed_msa_dir")
                        if msa_dir and opexists(msa_dir):
                            logger.warning(
                                "Use the old msa json format, change to pairedMsaPath/unpairedMsaPath field for future use."
                            )
                            if opexists(opjoin(msa_dir, "pairing.a3m")):
                                with open(opjoin(msa_dir, "pairing.a3m")) as f:
                                    p_a3m = f.read()
                            if opexists(opjoin(msa_dir, "non_pairing.a3m")):
                                with open(opjoin(msa_dir, "non_pairing.a3m")) as f:
                                    u_a3m = f.read()

            elif "rnaSequence" in info:
                c = info["rnaSequence"]
                seq, count, ctype = c["sequence"], c["count"], RNA_CHAIN
                if use_rna_msa:
                    u_a3m = c.get("unpairedMsa")
                    if u_a3m is None and c.get("unpairedMsaPath"):
                        with open(c["unpairedMsaPath"]) as f:
                            u_a3m = f.read()
            elif "dnaSequence" in info:
                c = info["dnaSequence"]
                seq, count, ctype = c["sequence"], c["count"], DNA_CHAIN
            elif "ligand" in info:
                count, ctype, seq = (
                    info["ligand"]["count"],
                    LIGAND_CHAIN_TYPES,
                    "X" * (atom_array.asym_id_int == curr_aid).sum(),
                )

            p_a3m = ensure_ends_with_newline(p_a3m)
            u_a3m = ensure_ends_with_newline(u_a3m)

            if msa_pair_as_unpair and p_a3m:
                u_a3m = RawMsa.from_a3m(
                    seq, ctype, p_a3m + (u_a3m or ""), dedup=True
                ).to_a3m()

            for c_idx in range(count):
                aid = curr_aid + c_idx
                meta[aid] = {
                    "entity_id": eid,
                    "chain_id": atom_array.chain_id[atom_array.asym_id_int == aid][0],
                    "sequence": seq,
                    "paired_msa": p_a3m or "",
                    "unpaired_msa": u_a3m or "",
                    "chain_entity_type": ctype,
                }
            curr_aid += count

        ca = atom_array[atom_array.centre_atom_mask.astype(bool)]
        std_idxs = map_to_standard(ca.asym_id_int, ca.res_id, meta)
        return FeatureAssemblyLine().assemble(meta, std_idxs).to_dict()
