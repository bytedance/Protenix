# Protenix-v2 Fine-tuning Smoke Test

This workflow verifies that a pretrained `protenix-v2` checkpoint can be loaded and fine-tuned end-to-end on a tiny local dataset. It is a debugging smoke test only; it is not intended to produce a useful scientific checkpoint or metric.

## Hardware

The acceptance target is a single node with 4 NVIDIA V100 32GB GPUs. Defaults use FP16, per-GPU dataloader batch size 1, diffusion batch size 1, one recycling cycle, one diffusion step, and 10 optimizer steps.

## Setup

Use the repository's normal environment setup. For an editable install:

```bash
pip install -e .
```

## Checkpoint

You can download the checkpoint manually:

```bash
mkdir -p "$PWD/output/protenix_v2_ft_smoke/data/checkpoint"
curl -L -o "$PWD/output/protenix_v2_ft_smoke/data/checkpoint/protenix-v2.pt" \
  "https://protenix.tos-cn-beijing.volces.com/checkpoint/protenix-v2.pt"
```

Or let the smoke script download it:

```bash
AUTO_DOWNLOAD_CHECKPOINT=true bash scripts/smoke_finetune_protenix_v2.sh
```

## 4-GPU Smoke Run

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3 \
PROTENIX_ROOT_DIR="$PWD/output/protenix_v2_ft_smoke/data" \
AUTO_DOWNLOAD_CHECKPOINT=true \
NPROC_PER_NODE=4 \
MAX_STEPS=10 \
DTYPE=fp16 \
TRAIN_CROP_SIZE=128 \
DIFFUSION_BATCH_SIZE=1 \
DIFFUSION_CHUNK_SIZE=1 \
N_CYCLE=1 \
N_STEP=1 \
N_SAMPLE=1 \
USE_WANDB=false \
bash scripts/smoke_finetune_protenix_v2.sh
```

The script prepares `output/protenix_v2_ft_smoke/data/tiny_finetune`, launches `torchrun`, writes `output/protenix_v2_ft_smoke/logs/smoke_train.log`, and writes `output/protenix_v2_ft_smoke/smoke_summary.json`.

After it finishes:

```bash
find output/protenix_v2_ft_smoke/runs -path '*/checkpoints/*.pt' -type f -print
cat output/protenix_v2_ft_smoke/smoke_summary.json
```

## Local Debug

This one-GPU command is only a developer convenience. It does not replace the 4xV100 acceptance test.

```bash
NPROC_PER_NODE=1 MAX_STEPS=1 TRAIN_CROP_SIZE=64 \
bash scripts/smoke_finetune_protenix_v2.sh
```

## Troubleshooting

- V100 BF16 errors: keep `DTYPE=fp16`.
- Custom kernel errors: keep `TRIANGLE_ATTENTION=torch`, `TRIANGLE_MULTIPLICATIVE=torch`, and `LAYERNORM_TYPE=torch`.
- OOM: first try `TRAIN_CROP_SIZE=64`; keep `DIFFUSION_BATCH_SIZE=1`, `N_CYCLE=1`, and `N_STEP=1`.
- Checkpoint `module.` key errors: verify the prefix adaptation in `runner/train.py`.
