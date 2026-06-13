#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)"
export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-$(pwd)}"

MODEL_NAME="${MODEL_NAME:-protenix-v2}"
CHECKPOINT_PATH="${CHECKPOINT_PATH:-${PROTENIX_ROOT_DIR}/checkpoint/${MODEL_NAME}.pt}"
FROM_SCRATCH_SMOKE="${FROM_SCRATCH_SMOKE:-false}"
TORCHRUN="${TORCHRUN:-torchrun}"
if ! command -v "${TORCHRUN}" >/dev/null 2>&1 && [[ -x ".venv/bin/torchrun" ]]; then
  TORCHRUN=".venv/bin/torchrun"
fi

checkpoint_args=()
if [[ "${FROM_SCRATCH_SMOKE}" != "true" ]]; then
  if [[ ! -f "${CHECKPOINT_PATH}" ]]; then
    echo "Missing checkpoint for continued pre-training: ${CHECKPOINT_PATH}" >&2
    echo "Set FROM_SCRATCH_SMOKE=true for a no-checkpoint debug run, or download the checkpoint first." >&2
    exit 2
  fi
  checkpoint_args=(
    --load_checkpoint_path "${CHECKPOINT_PATH}"
    --load_ema_checkpoint_path "${EMA_CHECKPOINT_PATH:-${CHECKPOINT_PATH}}"
    --load_params_only true
    --skip_load_optimizer true
    --skip_load_scheduler true
    --skip_load_step true
  )
fi

"${TORCHRUN}" \
  --nproc_per_node "${NPROC_PER_NODE:-1}" \
  --nnodes "${NNODES:-1}" \
  --node_rank "${NODE_RANK:-0}" \
  --master_addr "${MASTER_ADDR:-127.0.0.1}" \
  --master_port "${MASTER_PORT:-29501}" \
  runner/train.py \
  --model_name "${MODEL_NAME}" \
  --run_name "${RUN_NAME:-protenix_v2_openfold_subset_pretrain}" \
  --seed "${SEED:-42}" \
  --base_dir "${BASE_DIR:-./outputs/nanobody}" \
  --project "${WANDB_PROJECT:-protenix_nanobody_pretrain}" \
  --use_wandb "${USE_WANDB:-false}" \
  --dtype "${DTYPE:-bf16}" \
  --diffusion_batch_size "${DIFFUSION_BATCH_SIZE:-16}" \
  --eval_interval "${EVAL_INTERVAL:-1000}" \
  --log_interval "${LOG_INTERVAL:-50}" \
  --checkpoint_interval "${CHECKPOINT_INTERVAL:-1000}" \
  --ema_decay "${EMA_DECAY:-0.999}" \
  --train_crop_size "${TRAIN_CROP_SIZE:-384}" \
  --max_steps "${MAX_STEPS:-50000}" \
  --warmup_steps "${WARMUP_STEPS:-2000}" \
  --lr "${LR:-0.0003}" \
  --model.N_cycle "${N_CYCLE:-4}" \
  --sample_diffusion.N_step "${N_STEP:-20}" \
  --triangle_attention "${TRIANGLE_ATTENTION:-cuequivariance}" \
  --triangle_multiplicative "${TRIANGLE_MULTIPLICATIVE:-cuequivariance}" \
  "${checkpoint_args[@]}" \
  --data.train_sets "${TRAIN_SETS:-openfold_subset_train}" \
  --data.train_sampler.train_sample_weights "${TRAIN_SAMPLE_WEIGHTS:-1.0}" \
  --data.test_sets "${TEST_SETS:-nanobody_single_val_2022h1,nanobody_complex_val_2022h1}"

cat <<'EOF'
To evaluate a saved pre-training checkpoint on nanobody validation splits, run for example:
python scripts/nanobody/run_inference_checkpoint.py \
  --checkpoint outputs/nanobody/<run>/checkpoints/<step>_ema_0.999.pt \
  --model-name protenix-v2 \
  --input-json ${PROTENIX_ROOT_DIR}/nanobody/inference_jsons/complex_val_2022h1.json \
  --out-dir outputs/nanobody/<run>/val_complex_preds
EOF
