#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)"
export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-$(pwd)}"

MODEL_NAME="${MODEL_NAME:-protenix-v2}"
CHECKPOINT_PATH="${CHECKPOINT_PATH:-${PROTENIX_ROOT_DIR}/checkpoint/${MODEL_NAME}.pt}"
EMA_CHECKPOINT_PATH="${EMA_CHECKPOINT_PATH:-${CHECKPOINT_PATH}}"
TORCHRUN="${TORCHRUN:-torchrun}"
if ! command -v "${TORCHRUN}" >/dev/null 2>&1 && [[ -x ".venv/bin/torchrun" ]]; then
  TORCHRUN=".venv/bin/torchrun"
fi

if [[ ! -f "${CHECKPOINT_PATH}" ]]; then
  mkdir -p "$(dirname "${CHECKPOINT_PATH}")"
  if [[ "${AUTO_DOWNLOAD_CHECKPOINT:-false}" == "true" ]]; then
    python - "${MODEL_NAME}" "${CHECKPOINT_PATH}" <<'PY'
import sys
import urllib.request
from protenix.web_service.dependency_url import URL

model_name, checkpoint_path = sys.argv[1:3]
url = URL[model_name]
print(f"Downloading {url} -> {checkpoint_path}")
urllib.request.urlretrieve(url, checkpoint_path)
PY
  else
    url="$(python - "${MODEL_NAME}" <<'PY'
import sys
from protenix.web_service.dependency_url import URL
print(URL[sys.argv[1]])
PY
)"
    echo "Missing checkpoint: ${CHECKPOINT_PATH}" >&2
    echo "Download it with: wget -O '${CHECKPOINT_PATH}' '${url}'" >&2
    echo "Or rerun with AUTO_DOWNLOAD_CHECKPOINT=true." >&2
    exit 2
  fi
fi

"${TORCHRUN}" \
  --nproc_per_node "${NPROC_PER_NODE:-1}" \
  --nnodes "${NNODES:-1}" \
  --node_rank "${NODE_RANK:-0}" \
  --master_addr "${MASTER_ADDR:-127.0.0.1}" \
  --master_port "${MASTER_PORT:-29500}" \
  runner/train.py \
  --model_name "${MODEL_NAME}" \
  --run_name "${RUN_NAME:-protenix_v2_nanobody_ft}" \
  --seed "${SEED:-42}" \
  --base_dir "${BASE_DIR:-./outputs/nanobody}" \
  --project "${WANDB_PROJECT:-protenix_nanobody}" \
  --use_wandb "${USE_WANDB:-false}" \
  --dtype "${DTYPE:-bf16}" \
  --diffusion_batch_size "${DIFFUSION_BATCH_SIZE:-16}" \
  --eval_interval "${EVAL_INTERVAL:-400}" \
  --log_interval "${LOG_INTERVAL:-50}" \
  --checkpoint_interval "${CHECKPOINT_INTERVAL:-400}" \
  --ema_decay "${EMA_DECAY:-0.999}" \
  --train_crop_size "${TRAIN_CROP_SIZE:-384}" \
  --max_steps "${MAX_STEPS:-20000}" \
  --warmup_steps "${WARMUP_STEPS:-1000}" \
  --lr "${LR:-0.0003}" \
  --model.N_cycle "${N_CYCLE:-4}" \
  --sample_diffusion.N_step "${N_STEP:-20}" \
  --triangle_attention "${TRIANGLE_ATTENTION:-cuequivariance}" \
  --triangle_multiplicative "${TRIANGLE_MULTIPLICATIVE:-cuequivariance}" \
  --load_checkpoint_path "${CHECKPOINT_PATH}" \
  --load_ema_checkpoint_path "${EMA_CHECKPOINT_PATH}" \
  --load_params_only true \
  --skip_load_optimizer true \
  --skip_load_scheduler true \
  --skip_load_step true \
  --data.train_sets "${TRAIN_SETS:-nanobody_single_train_pre2022,nanobody_complex_train_pre2022}" \
  --data.train_sampler.train_sample_weights "${TRAIN_SAMPLE_WEIGHTS:-0.35,0.65}" \
  --data.test_sets "${TEST_SETS:-nanobody_single_val_2022h1,nanobody_complex_val_2022h1}"
