#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}${PWD}"

MODEL_NAME="${MODEL_NAME:-protenix-v2}"
PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-${PWD}/output/protenix_v2_ft_smoke/data}"
BASE_DIR="${BASE_DIR:-${PWD}/output/protenix_v2_ft_smoke/runs}"
SMOKE_ROOT="${SMOKE_ROOT:-${PWD}/output/protenix_v2_ft_smoke}"
NPROC_PER_NODE="${NPROC_PER_NODE:-4}"
MAX_STEPS="${MAX_STEPS:-10}"
DTYPE="${DTYPE:-fp16}"
TRAIN_CROP_SIZE="${TRAIN_CROP_SIZE:-128}"
DIFFUSION_BATCH_SIZE="${DIFFUSION_BATCH_SIZE:-1}"
DIFFUSION_CHUNK_SIZE="${DIFFUSION_CHUNK_SIZE:-1}"
N_CYCLE="${N_CYCLE:-1}"
N_STEP="${N_STEP:-1}"
N_SAMPLE="${N_SAMPLE:-1}"
USE_WANDB="${USE_WANDB:-false}"
TRIANGLE_ATTENTION="${TRIANGLE_ATTENTION:-torch}"
TRIANGLE_MULTIPLICATIVE="${TRIANGLE_MULTIPLICATIVE:-torch}"
LAYERNORM_TYPE="${LAYERNORM_TYPE:-torch}"
RUN_NAME="${RUN_NAME:-protenix_v2_finetune_smoke}"
EPOCH_SIZE="${EPOCH_SIZE:-$((MAX_STEPS * NPROC_PER_NODE))}"

export PROTENIX_ROOT_DIR
export TRIANGLE_ATTENTION
export TRIANGLE_MULTIPLICATIVE
export LAYERNORM_TYPE

CHECKPOINT_PATH="${CHECKPOINT_PATH:-${PROTENIX_ROOT_DIR}/checkpoint/${MODEL_NAME}.pt}"
MIN_CHECKPOINT_BYTES="${MIN_CHECKPOINT_BYTES:-1000000}"
DATA_ROOT="${DATA_ROOT:-${PROTENIX_ROOT_DIR}/tiny_finetune}"
DATA_SUMMARY="${DATA_ROOT}/prepare_summary.json"
SMOKE_SUMMARY="${SMOKE_SUMMARY:-${SMOKE_ROOT}/smoke_summary.json}"
TRAIN_LOG="${TRAIN_LOG:-${SMOKE_ROOT}/logs/smoke_train.log}"

if [[ -z "${PYTHON:-}" ]]; then
  if [[ -x ".venv/bin/python" ]]; then
    PYTHON=".venv/bin/python"
  elif [[ -x "../.venv/bin/python" ]]; then
    PYTHON="../.venv/bin/python"
  else
    PYTHON="python"
  fi
elif ! command -v "${PYTHON}" >/dev/null 2>&1; then
  echo "Python executable not found: ${PYTHON}" >&2
  exit 2
fi

if [[ -z "${TORCHRUN:-}" ]]; then
  if [[ -x ".venv/bin/torchrun" ]]; then
    TORCHRUN=".venv/bin/torchrun"
  elif [[ -x "../.venv/bin/torchrun" ]]; then
    TORCHRUN="../.venv/bin/torchrun"
  else
    TORCHRUN="torchrun"
  fi
elif ! command -v "${TORCHRUN}" >/dev/null 2>&1; then
  echo "torchrun executable not found: ${TORCHRUN}" >&2
  exit 2
fi

mkdir -p "$(dirname "${CHECKPOINT_PATH}")" "${BASE_DIR}" "$(dirname "${TRAIN_LOG}")"

checkpoint_needs_download=false
if [[ ! -f "${CHECKPOINT_PATH}" ]]; then
  checkpoint_needs_download=true
elif [[ "$(stat -Lc%s "${CHECKPOINT_PATH}")" -lt "${MIN_CHECKPOINT_BYTES}" ]]; then
  echo "Checkpoint exists but is too small to be valid: ${CHECKPOINT_PATH}" >&2
  checkpoint_needs_download=true
fi

if [[ "${checkpoint_needs_download}" == "true" ]]; then
  if [[ "${AUTO_DOWNLOAD_CHECKPOINT:-false}" == "true" ]]; then
    rm -f "${CHECKPOINT_PATH}"
    "${PYTHON}" - "${MODEL_NAME}" "${CHECKPOINT_PATH}" <<'PY'
import sys
import urllib.request
from pathlib import Path
from protenix.web_service.dependency_url import URL

model_name, checkpoint_path = sys.argv[1:3]
url = URL[model_name]
Path(checkpoint_path).parent.mkdir(parents=True, exist_ok=True)
print(f"Downloading {url} -> {checkpoint_path}")
urllib.request.urlretrieve(url, checkpoint_path)
PY
    if [[ ! -f "${CHECKPOINT_PATH}" || "$(stat -Lc%s "${CHECKPOINT_PATH}")" -lt "${MIN_CHECKPOINT_BYTES}" ]]; then
      echo "Downloaded checkpoint is missing or too small: ${CHECKPOINT_PATH}" >&2
      exit 2
    fi
  else
    checkpoint_url="$("${PYTHON}" - "${MODEL_NAME}" <<'PY'
import sys
from protenix.web_service.dependency_url import URL
print(URL[sys.argv[1]])
PY
)"
    echo "Missing checkpoint: ${CHECKPOINT_PATH}" >&2
    echo "Download it with:" >&2
    echo "  curl -L -o '${CHECKPOINT_PATH}' '${checkpoint_url}'" >&2
    echo "Or rerun with AUTO_DOWNLOAD_CHECKPOINT=true." >&2
    exit 2
  fi
fi

prepare_args=(
  scripts/smoke/prepare_tiny_finetune_data.py
  --out-root "${DATA_ROOT}"
  --max-structures "${MAX_STRUCTURES:-6}"
  --n-train "${N_TRAIN:-4}"
  --n-val "${N_VAL:-2}"
  --max-n-token "${MAX_N_TOKEN:-256}"
  --n-cpu "${N_CPU:-1}"
)
if [[ "${PREPARE_FORCE:-false}" == "true" ]]; then
  prepare_args+=(--force)
fi
"${PYTHON}" "${prepare_args[@]}"

torchrun_args=(
  --nproc_per_node "${NPROC_PER_NODE}"
  --nnodes "${NNODES:-1}"
  --node_rank "${NODE_RANK:-0}"
  --master_addr "${MASTER_ADDR:-127.0.0.1}"
  --master_port "${MASTER_PORT:-29511}"
  runner/train.py
  --model_name "${MODEL_NAME}"
  --run_name "${RUN_NAME}"
  --seed "${SEED:-42}"
  --base_dir "${BASE_DIR}"
  --project "${WANDB_PROJECT:-protenix_v2_finetune_smoke}"
  --use_wandb "${USE_WANDB}"
  --dtype "${DTYPE}"
  --diffusion_batch_size "${DIFFUSION_BATCH_SIZE}"
  --diffusion_chunk_size "${DIFFUSION_CHUNK_SIZE}"
  --eval_interval "${EVAL_INTERVAL:-10}"
  --log_interval "${LOG_INTERVAL:-1}"
  --checkpoint_interval "${CHECKPOINT_INTERVAL:-10}"
  --ema_decay "${EMA_DECAY:--1}"
  --train_crop_size "${TRAIN_CROP_SIZE}"
  --test_max_n_token "${TEST_MAX_N_TOKEN:-256}"
  --max_steps "${MAX_STEPS}"
  --warmup_steps "${WARMUP_STEPS:-2}"
  --lr "${LR:-1e-5}"
  --iters_to_accumulate "${ITERS_TO_ACCUMULATE:-1}"
  --data.num_dl_workers "${NUM_DL_WORKERS:-0}"
  --data.epoch_size "${EPOCH_SIZE}"
  --model.N_cycle "${N_CYCLE}"
  --sample_diffusion.N_step "${N_STEP}"
  --sample_diffusion.N_sample "${N_SAMPLE}"
  --sample_diffusion.N_step_mini_rollout "${N_STEP_MINI_ROLLOUT:-1}"
  --sample_diffusion.N_sample_mini_rollout "${N_SAMPLE_MINI_ROLLOUT:-1}"
  --triangle_attention "${TRIANGLE_ATTENTION}"
  --triangle_multiplicative "${TRIANGLE_MULTIPLICATIVE}"
  --load_checkpoint_path "${CHECKPOINT_PATH}"
  --load_params_only true
  --skip_load_optimizer true
  --skip_load_scheduler true
  --skip_load_step true
  --data.train_sets tiny_protenix_finetune_train
  --data.train_sampler.train_sample_weights 1.0
  --data.test_sets tiny_protenix_finetune_val
)

"${TORCHRUN}" "${torchrun_args[@]}" 2>&1 | tee "${TRAIN_LOG}"

"${PYTHON}" - \
  "${SMOKE_SUMMARY}" \
  "${SMOKE_ROOT}" \
  "${BASE_DIR}" \
  "${RUN_NAME}" \
  "${CHECKPOINT_PATH}" \
  "${DATA_SUMMARY}" \
  "${TRAIN_LOG}" \
  "${MODEL_NAME}" \
  "${NPROC_PER_NODE}" \
  "${MAX_STEPS}" \
  "${DTYPE}" \
  "${TRAIN_CROP_SIZE}" \
  "${DIFFUSION_BATCH_SIZE}" \
  "${DIFFUSION_CHUNK_SIZE}" \
  "${N_CYCLE}" \
  "${N_STEP}" \
  "${N_SAMPLE}" \
  "${USE_WANDB}" <<'PY'
import json
import sys
from pathlib import Path

(
    summary_path,
    smoke_root,
    base_dir,
    run_name,
    checkpoint_path,
    data_summary,
    train_log,
    model_name,
    nproc_per_node,
    max_steps,
    dtype,
    train_crop_size,
    diffusion_batch_size,
    diffusion_chunk_size,
    n_cycle,
    n_step,
    n_sample,
    use_wandb,
) = sys.argv[1:]

base = Path(base_dir)
run_dirs = sorted(
    [p for p in base.glob(f"{run_name}_*") if p.is_dir()],
    key=lambda p: p.stat().st_mtime,
    reverse=True,
)
run_dir = run_dirs[0] if run_dirs else None
checkpoint_files = (
    sorted(str(p) for p in (run_dir / "checkpoints").glob("*.pt"))
    if run_dir is not None
    else []
)
log_text = Path(train_log).read_text(errors="replace") if Path(train_log).exists() else ""
finished_marker = f"Finished training after {max_steps} steps"
summary = {
    "model_name": model_name,
    "settings": {
        "nproc_per_node": int(nproc_per_node),
        "max_steps": int(max_steps),
        "dtype": dtype,
        "train_crop_size": int(train_crop_size),
        "diffusion_batch_size": int(diffusion_batch_size),
        "diffusion_chunk_size": int(diffusion_chunk_size),
        "n_cycle": int(n_cycle),
        "n_step": int(n_step),
        "n_sample": int(n_sample),
        "use_wandb": use_wandb,
    },
    "smoke_root": str(Path(smoke_root).resolve()),
    "data_summary": str(Path(data_summary).resolve()),
    "run_dir": str(run_dir.resolve()) if run_dir is not None else None,
    "checkpoint_path": str(Path(checkpoint_path).absolute()),
    "checkpoint_realpath": str(Path(checkpoint_path).resolve()),
    "saved_checkpoints": checkpoint_files,
    "final_checkpoint_exists": bool(checkpoint_files),
    "train_log": str(Path(train_log).resolve()),
    "finished_marker": finished_marker,
    "finished_marker_found": finished_marker in log_text,
}
Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
Path(summary_path).write_text(json.dumps(summary, indent=2) + "\n")
print(f"Wrote smoke summary to {summary_path}")
if not summary["final_checkpoint_exists"]:
    raise SystemExit("No checkpoint was saved by the smoke run.")
if not summary["finished_marker_found"]:
    raise SystemExit(f"Did not find log marker: {finished_marker}")
PY
