#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)"
if [[ -z "${PYTHON_BIN:-}" ]]; then
  if [[ -x ".venv/bin/python" ]]; then
    PYTHON_BIN=".venv/bin/python"
  else
    PYTHON_BIN="$(command -v python3.12 || command -v python3.11 || command -v python3 || command -v python)"
  fi
fi

missing="$("${PYTHON_BIN}" - <<'PY'
missing = []
for module in ("pandas", "numpy", "biotite", "tqdm"):
    try:
        __import__(module)
    except Exception:
        missing.append(module)
print(",".join(missing))
PY
)"
if [[ -n "${missing}" ]]; then
  echo "Skipping nanobody smoke: ${PYTHON_BIN} is missing required smoke dependencies: ${missing}."
  echo "Install Protenix runtime requirements, then rerun this script."
  exit 0
fi

tmp_root="${TMPDIR:-/tmp}/protenix_nanobody_smoke_$$"
data_root="${tmp_root}/protenix_data"
mkdir -p "${tmp_root}/raw/structures" "${tmp_root}/preds/2LWU_single/seed_101/predictions" "${data_root}"
trap 'if [[ "${KEEP_NANOBODY_SMOKE:-false}" != "true" ]]; then rm -rf "${tmp_root}"; fi' EXIT

cache_root="${PROTENIX_SMOKE_CACHE_ROOT:-}"
if [[ -z "${cache_root}" && -f "output/smoke_uv/data/common/components.cif" ]]; then
  cache_root="$(pwd)/output/smoke_uv/data"
fi
prepare_extra_args=(--skip-inference-jsons)
if [[ -n "${cache_root}" && -f "${cache_root}/common/components.cif" ]]; then
  for cache_dir in common rna_msa mmcif_msa_template; do
    if [[ -e "${cache_root}/${cache_dir}" ]]; then
      ln -s "${cache_root}/${cache_dir}" "${data_root}/${cache_dir}"
    fi
  done
  export PROTENIX_ROOT_DIR="${data_root}"
  echo "Using Protenix smoke cache: ${cache_root}/common"
else
  export PROTENIX_ROOT_DIR="${data_root}"
  prepare_extra_args+=(--skip-protenix-conversion)
  echo "No Protenix common cache found; smoke will skip DataPipeline conversion."
fi

cp examples/2lwu.cif "${tmp_root}/raw/structures/2lwu.cif"
cat > "${tmp_root}/raw/sabdab_metadata.csv" <<EOF
pdb,release_date,Hchain,antigen_chain,structure_path,cdrh3_start,cdrh3_end
2LWU,2021-01-01,A,,${tmp_root}/raw/structures/2lwu.cif,3,8
2LWU,2022-02-01,A,,${tmp_root}/raw/structures/2lwu.cif,3,8
EOF

"${PYTHON_BIN}" scripts/nanobody/prepare_nanobody_data.py \
  --metadata "${tmp_root}/raw/sabdab_metadata.csv" \
  --structure-dir "${tmp_root}/raw/structures" \
  --out-root "${PROTENIX_ROOT_DIR}/nanobody" \
  "${prepare_extra_args[@]}"

cp "${PROTENIX_ROOT_DIR}/nanobody/mmcif/single/2lwu.cif" \
  "${tmp_root}/preds/2LWU_single/seed_101/predictions/2LWU_single_sample_0.cif"
echo "Skipping real inference in smoke; using the reference structure as a fake prediction because released checkpoints may be unavailable."
cat > "${tmp_root}/preds/2LWU_single/seed_101/predictions/2LWU_single_summary_confidence_sample_0.json" <<EOF
{"ranking_score": 1.0, "plddt": 100.0, "gpde": 0.0, "ptm": 1.0}
EOF

"${PYTHON_BIN}" scripts/nanobody/evaluate_nanobody_predictions.py \
  --manifest "${PROTENIX_ROOT_DIR}/nanobody/manifests/nanobody_single_val_2022h1.csv" \
  --pred-dir "${tmp_root}/preds" \
  --out-csv "${tmp_root}/metrics.csv" \
  --out-json "${tmp_root}/summary.json" \
  --sample-mode rank0

run_training="${PROTENIX_SMOKE_RUN_TRAINING:-auto}"
if [[ "${run_training}" == "auto" ]]; then
  if [[ -f "${PROTENIX_ROOT_DIR}/nanobody/bioassembly/single/2lwu.pkl.gz" ]] \
    && { command -v torchrun >/dev/null 2>&1 || [[ -x ".venv/bin/torchrun" ]]; }; then
    run_training="true"
  else
    run_training="false"
  fi
fi

if [[ "${run_training}" == "true" ]]; then
  MODEL_NAME="${MODEL_NAME:-protenix_tiny_default_v0.5.0}" \
  FROM_SCRATCH_SMOKE=true \
  LAYERNORM_TYPE=torch \
  MAX_STEPS=1 \
  N_STEP=1 \
  N_CYCLE=1 \
  TRAIN_CROP_SIZE=64 \
  DIFFUSION_BATCH_SIZE=1 \
  TRIANGLE_ATTENTION=torch \
  TRIANGLE_MULTIPLICATIVE=torch \
  TRAIN_SETS=nanobody_single_train_pre2022 \
  TEST_SETS=nanobody_single_val_2022h1 \
  scripts/nanobody/pretrain_protenix_v2_openfold_subset.sh
else
  echo "Skipping one-step training smoke because DataPipeline conversion or torchrun is unavailable. Set PROTENIX_SMOKE_RUN_TRAINING=true to force it."
fi

echo "Nanobody smoke outputs: ${tmp_root}"
