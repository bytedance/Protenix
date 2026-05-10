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

import copy
import csv
import json
import logging
import tempfile
import traceback
from pathlib import Path
from typing import Any, Optional

import numpy as np
import torch

from protenix.data.inference.infer_dataloader import InferenceDataset
from protenix.data.inference.structure_scoring import (
    build_chain_id_map,
    build_source_coord_map,
    collect_structure_files,
    map_source_coords_to_internal,
    sanitize_sample_name,
    structure_to_score_input_json,
)
from protenix.data.utils import save_structure_cif
from protenix.utils.file_io import save_json
from protenix.utils.torch_utils import round_values
from runner.dumper import get_clean_full_confidence
from runner.inference import update_inference_configs

logger = logging.getLogger(__name__)


def score_structures(
    input_path: str,
    out_dir: str = "./output",
    recursive: bool = False,
    glob: str = "*.pdb,*.cif",
    sample_name: Optional[str] = None,
    missing_atom_policy: str = "error",
    write_full_confidence: bool = False,
    write_scored_cif: bool = False,
    assembly_id: Optional[str] = None,
    altloc: str = "first",
    include_discont_poly_poly_bonds: bool = True,
    use_msa: bool = True,
    use_template: bool = False,
    use_rna_msa: bool = False,
    msa_server_mode: str = "protenix",
    hmmsearch_binary_path: Optional[str] = None,
    hmmbuild_binary_path: Optional[str] = None,
    seqres_database_path: Optional[str] = None,
    nhmmer_binary_path: Optional[str] = None,
    hmmalign_binary_path: Optional[str] = None,
    hmmbuild_rna_binary_path: Optional[str] = None,
    ntrna_database_path: Optional[str] = None,
    rfam_database_path: Optional[str] = None,
    rna_central_database_path: Optional[str] = None,
    nhmmer_n_cpu: Optional[int] = None,
    cycle: int = 10,
    dtype: str = "bf16",
    model_name: str = "protenix_base_default_v1.0.0",
    trimul_kernel: str = "cuequivariance",
    triatt_kernel: str = "cuequivariance",
    enable_cache: bool = True,
    enable_fusion: bool = True,
    enable_tf32: bool = True,
    kalign_binary_path: Optional[str] = None,
) -> dict[str, Any]:
    """
    Score PDB/mmCIF structures without running diffusion sampling.
    """
    from runner.batch_inference import get_default_runner, preprocess_input

    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)

    globs = [pattern.strip() for pattern in glob.split(",") if pattern.strip()]
    structure_files = collect_structure_files(
        input_path=input_path,
        recursive=recursive,
        globs=globs,
    )
    if not structure_files:
        raise RuntimeError(
            f"can not read a valid `pdb` or `cif` file from {input_path}"
        )
    if sample_name is not None and len(structure_files) > 1:
        raise RuntimeError("--sample_name can only be used when scoring one structure")

    runner = get_default_runner(
        seeds=[101],
        n_cycle=cycle,
        n_step=1,
        n_sample=1,
        dtype=dtype,
        model_name=model_name,
        use_msa=use_msa,
        trimul_kernel=trimul_kernel,
        triatt_kernel=triatt_kernel,
        enable_cache=enable_cache,
        enable_fusion=enable_fusion,
        enable_tf32=enable_tf32,
        use_template=use_template,
        use_rna_msa=use_rna_msa,
        need_atom_confidence=write_full_confidence,
        kalign_binary_path=kalign_binary_path,
    )

    rows: list[dict[str, Any]] = []
    failed_records: list[str] = []
    seen_names: set[str] = set()
    with tempfile.TemporaryDirectory(prefix="protenix-score-") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for structure_file in structure_files:
            base_name = sample_name or sanitize_sample_name(structure_file.stem)
            current_sample_name = _unique_sample_name(base_name, seen_names)
            sample_out_dir = out_dir_path / current_sample_name
            sample_out_dir.mkdir(parents=True, exist_ok=True)
            try:
                logger.info("Scoring %s as %s", structure_file, current_sample_name)
                input_json, source_atom_array, chain_ids_preserved = (
                    structure_to_score_input_json(
                        input_file=structure_file,
                        sample_name=current_sample_name,
                        assembly_id=assembly_id,
                        altloc=altloc,
                        include_discont_poly_poly_bonds=include_discont_poly_poly_bonds,
                    )
                )
                intermediate_json = tmp_dir_path / f"{current_sample_name}.json"
                with open(intermediate_json, "w") as f:
                    json.dump([input_json], f, indent=4)

                processed_json = preprocess_input(
                    input_json=str(intermediate_json),
                    out_dir=str(sample_out_dir),
                    use_msa=use_msa,
                    use_template=use_template,
                    use_rna_msa=use_rna_msa,
                    msa_server_mode=msa_server_mode,
                    hmmsearch_binary_path=hmmsearch_binary_path,
                    hmmbuild_binary_path=hmmbuild_binary_path,
                    seqres_database_path=seqres_database_path,
                    nhmmer_binary_path=nhmmer_binary_path,
                    hmmalign_binary_path=hmmalign_binary_path,
                    hmmbuild_rna_binary_path=hmmbuild_rna_binary_path,
                    ntrna_database_path=ntrna_database_path,
                    rfam_database_path=rfam_database_path,
                    rna_central_database_path=rna_central_database_path,
                    nhmmer_n_cpu=nhmmer_n_cpu,
                )

                runner.configs.input_json_path = processed_json
                runner.configs.dump_dir = str(out_dir_path)
                dataset = InferenceDataset(runner.configs)
                data, internal_atom_array, data_error_message = dataset[0]
                if data_error_message:
                    raise RuntimeError(data_error_message)

                new_configs = update_inference_configs(
                    runner.configs,
                    data["N_token"].item(),
                )
                runner.update_model_configs(new_configs)

                chain_id_map = build_chain_id_map(
                    source_atom_array=source_atom_array,
                    internal_atom_array=internal_atom_array,
                    prefer_identity=chain_ids_preserved,
                )
                chain_id_map["source_chain_ids_preserved_in_json"] = chain_ids_preserved
                source_coord_map = build_source_coord_map(source_atom_array)
                score_coordinates, missing_atoms = map_source_coords_to_internal(
                    source_coord_map=source_coord_map,
                    internal_atom_array=internal_atom_array,
                    chain_id_map=chain_id_map,
                    missing_atom_policy=missing_atom_policy,
                )

                prediction = runner.score(
                    data=data,
                    coordinates=torch.from_numpy(score_coordinates),
                )
                _save_score_outputs(
                    prediction=prediction,
                    sample_out_dir=sample_out_dir,
                    sample_name=current_sample_name,
                    atom_array=internal_atom_array,
                    entity_poly_type=data["entity_poly_type"],
                    chain_id_map=chain_id_map,
                    missing_atoms=missing_atoms,
                    write_full_confidence=write_full_confidence,
                    write_scored_cif=write_scored_cif,
                )
                rows.append(
                    _summary_row(
                        sample_name=current_sample_name,
                        input_path=structure_file,
                        summary_confidence=prediction["summary_confidence"][0],
                    )
                )
                logger.info("Score-only run succeeded for %s", current_sample_name)
            except Exception as exc:
                error_message = f"{structure_file}\t{type(exc).__name__}: {exc}"
                failed_records.append(error_message)
                with open(sample_out_dir / "error.txt", "w", encoding="utf-8") as f:
                    f.write(f"{error_message}\n{traceback.format_exc()}")
                logger.error(
                    "Score-only run failed for %s: %s",
                    structure_file,
                    traceback.format_exc(),
                )
            finally:
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

    _write_summary_csv(rows=rows, output_path=out_dir_path / "summary.csv")
    if failed_records:
        with open(out_dir_path / "failed_records.txt", "w", encoding="utf-8") as f:
            f.write("\n".join(failed_records) + "\n")
    if not rows and failed_records:
        raise RuntimeError(
            f"all {len(failed_records)} score-only jobs failed; "
            f"see {out_dir_path / 'failed_records.txt'}"
        )
    return {
        "succeeded": len(rows),
        "failed": len(failed_records),
        "out_dir": str(out_dir_path),
    }


def _save_score_outputs(
    prediction: dict[str, Any],
    sample_out_dir: Path,
    sample_name: str,
    atom_array: Any,
    entity_poly_type: dict[str, str],
    chain_id_map: dict[str, Any],
    missing_atoms: list[dict[str, Any]],
    write_full_confidence: bool,
    write_scored_cif: bool,
) -> None:
    save_json(
        _jsonable(round_values(copy.deepcopy(prediction["summary_confidence"][0]))),
        sample_out_dir / "summary_confidence.json",
        indent=4,
    )
    save_json(
        chain_id_map,
        sample_out_dir / "chain_id_map.json",
        indent=4,
    )
    if missing_atoms:
        save_json(
            {"missing_atoms": missing_atoms},
            sample_out_dir / "missing_atoms.json",
            indent=4,
        )
    if write_full_confidence:
        save_json(
            get_clean_full_confidence(copy.deepcopy(prediction["full_data"][0])),
            sample_out_dir / "full_confidence.json",
            indent=None,
        )
    if write_scored_cif:
        b_factor = _atom_plddt_b_factor(prediction["full_data"][0])
        if b_factor is not None:
            atom_array.set_annotation("b_factor", b_factor)
        save_structure_cif(
            atom_array=atom_array,
            pred_coordinate=prediction["coordinate"][0],
            output_fpath=str(sample_out_dir / "scored.cif"),
            entity_poly_type={
                k: v for k, v in entity_poly_type.items() if v != "non-polymer"
            },
            pdb_id=sample_name,
        )


def _atom_plddt_b_factor(full_data: dict[str, Any]) -> Optional[np.ndarray]:
    atom_plddt = full_data.get("atom_plddt")
    if atom_plddt is None:
        return None
    if isinstance(atom_plddt, torch.Tensor):
        if atom_plddt.dtype == torch.bfloat16:
            atom_plddt = atom_plddt.float()
        atom_plddt = atom_plddt.detach().cpu().numpy()
    return np.round(np.asarray(atom_plddt, dtype=np.float32) * 100.0, 2)


def _summary_row(
    sample_name: str,
    input_path: Path,
    summary_confidence: dict[str, Any],
) -> dict[str, Any]:
    summary = _jsonable(round_values(copy.deepcopy(summary_confidence)))
    row = {
        "name": sample_name,
        "input_path": str(input_path),
    }
    row.update(
        {
            key: value if _is_csv_scalar(value) else json.dumps(value)
            for key, value in summary.items()
        }
    )
    return row


def _write_summary_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    if not rows:
        return
    fieldnames = ["name", "input_path"]
    fieldnames.extend(
        sorted({key for row in rows for key in row.keys()} - set(fieldnames))
    )
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _jsonable(value: Any) -> Any:
    if isinstance(value, torch.Tensor):
        if value.dtype == torch.bfloat16:
            value = value.float()
        value = value.detach().cpu().numpy()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(v) for key, v in value.items()}
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    return value


def _is_csv_scalar(value: Any) -> bool:
    return value is None or isinstance(value, (str, int, float, bool))


def _unique_sample_name(name: str, seen_names: set[str]) -> str:
    candidate = name
    suffix = 2
    while candidate in seen_names:
        candidate = f"{name}_{suffix}"
        suffix += 1
    seen_names.add(candidate)
    return candidate
