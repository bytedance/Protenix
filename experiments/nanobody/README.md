# Protenix Nanobody Workflow

This workflow lets three teams run comparable Protenix-v2 experiments for nanobody folding and nanobody-antigen complex prediction without any ESMFold, Biohub ESM, ESMC, or external repository context.

## Install And Paths

Install Protenix in the normal local development environment:

```bash
pip install -e .
pip install -r requirements.txt
```

Set the data root used by `configs/configs_data.py`:

```bash
export PROTENIX_ROOT_DIR=/path/to/protenix_data
mkdir -p ${PROTENIX_ROOT_DIR}/nanobody
```

Large inputs and outputs should stay below `${PROTENIX_ROOT_DIR}/nanobody` or `outputs/nanobody`. Do not commit checkpoints, SAbDab downloads, OpenFold structures, MSA/template databases, or generated predictions.

## Prepare SAbDab Data

Place local SAbDab-like metadata and structures here:

```text
${PROTENIX_ROOT_DIR}/nanobody/raw/sabdab_metadata.csv
${PROTENIX_ROOT_DIR}/nanobody/raw/structures/
```

Run:

```bash
python scripts/nanobody/prepare_nanobody_data.py \
  --metadata ${PROTENIX_ROOT_DIR}/nanobody/raw/sabdab_metadata.csv \
  --structure-dir ${PROTENIX_ROOT_DIR}/nanobody/raw/structures \
  --out-root ${PROTENIX_ROOT_DIR}/nanobody \
  --n-cpu 16
```

The script accepts common metadata aliases for PDB ID, release date, nanobody/heavy chain, antigen chains, structure path, and optional CDR-H3 bounds. It writes canonical manifests, Protenix index CSVs, bioassembly pickles, and inference JSONs under `${PROTENIX_ROOT_DIR}/nanobody`.

Use `--strict-protenix-cutoff` for a 2021-09-30 train cutoff.

## Download Structure Datasets

Additional structure-dataset download workflows live under `scripts/nanobody/datasets/`. They default to writing below `${PROTENIX_ROOT_DIR}/nanobody`, or `./nanobody` if `PROTENIX_ROOT_DIR` is unset. Override with `OUT_ROOT=/path/to/output`.

Download and build the temporal SAbDab-nano bundle:

```bash
OUT_ROOT=${PROTENIX_ROOT_DIR}/nanobody \
bash scripts/nanobody/datasets/download_sabdab_nano_temporal.sh
```

Download and verify the tFold SAbDab-22H2 nanobody benchmarks:

```bash
OUT_ROOT=${PROTENIX_ROOT_DIR}/nanobody \
bash scripts/nanobody/datasets/download_tfold_sabdab22h2.sh
```

These scripts write source URLs, SHA256 files, canonical manifests, and verification reports under `data/sabdab_nano_temporal/` or `data/tfold_sabdab22h2/` inside the selected output root. Keep those downloaded datasets outside git.

## Baseline Validation

```bash
python scripts/nanobody/run_inference_checkpoint.py \
  --checkpoint ${PROTENIX_ROOT_DIR}/checkpoint/protenix-v2.pt \
  --model-name protenix-v2 \
  --input-json ${PROTENIX_ROOT_DIR}/nanobody/inference_jsons/complex_val_2022h1.json \
  --out-dir outputs/nanobody/baseline/complex_val_preds \
  --seeds 101,102,103 \
  --n-sample 5 \
  --n-step 200 \
  --n-cycle 10 \
  --use-msa true \
  --use-template true
```

Evaluate complex predictions:

```bash
python scripts/nanobody/evaluate_nanobody_predictions.py \
  --manifest ${PROTENIX_ROOT_DIR}/nanobody/manifests/nanobody_complex_val_2022h1.csv \
  --pred-dir outputs/nanobody/baseline/complex_val_preds \
  --out-csv outputs/nanobody/baseline/complex_val_metrics.csv \
  --out-json outputs/nanobody/baseline/complex_val_summary.json \
  --sample-mode rank0
```

Evaluate single-chain folding:

```bash
python scripts/nanobody/evaluate_nanobody_predictions.py \
  --manifest ${PROTENIX_ROOT_DIR}/nanobody/manifests/nanobody_single_val_2022h1.csv \
  --pred-dir outputs/nanobody/baseline/single_val_preds \
  --out-csv outputs/nanobody/baseline/single_val_metrics.csv \
  --out-json outputs/nanobody/baseline/single_val_summary.json \
  --sample-mode rank0
```

DockQ is optional. If it is not installed, the evaluator still writes CSV/JSON with confidence fields and records `dockq_unavailable`.

## Post-Training Team

Smoke fine-tuning:

```bash
MAX_STEPS=1 N_STEP=1 N_CYCLE=1 DIFFUSION_BATCH_SIZE=1 TRAIN_CROP_SIZE=64 \
TRIANGLE_ATTENTION=torch TRIANGLE_MULTIPLICATIVE=torch \
scripts/nanobody/finetune_protenix_v2_nanobody.sh
```

Real fine-tuning:

```bash
RUN_NAME=protenix_v2_nanobody_ft \
BASE_DIR=outputs/nanobody/posttrain \
MAX_STEPS=20000 \
LR=0.0003 \
TRAIN_SETS=nanobody_single_train_pre2022,nanobody_complex_train_pre2022 \
TRAIN_SAMPLE_WEIGHTS=0.35,0.65 \
TEST_SETS=nanobody_single_val_2022h1,nanobody_complex_val_2022h1 \
scripts/nanobody/finetune_protenix_v2_nanobody.sh
```

The wrapper defaults to `MODEL_NAME=protenix-v2`, loads `${PROTENIX_ROOT_DIR}/checkpoint/protenix-v2.pt`, skips optimizer/scheduler/step loading, and uses the nanobody dataset config entries.

## Pre-Training Team

Prepare an OpenFold-style subset:

```bash
python scripts/nanobody/prepare_openfold_subset.py \
  --input-dir /path/to/openfold/mmcif \
  --out-root ${PROTENIX_ROOT_DIR}/nanobody \
  --max-tokens 1536 \
  --max-resolution 4.5 \
  --max-release-date 2022-01-01
```

Continued pre-training:

```bash
RUN_NAME=protenix_v2_openfold_subset_pretrain \
BASE_DIR=outputs/nanobody/pretrain \
TRAIN_SETS=openfold_subset_train \
TEST_SETS=nanobody_single_val_2022h1,nanobody_complex_val_2022h1 \
scripts/nanobody/pretrain_protenix_v2_openfold_subset.sh
```

For from-scratch debug:

```bash
FROM_SCRATCH_SMOKE=true MODEL_NAME=protenix_tiny_default_v0.5.0 MAX_STEPS=1 \
N_STEP=1 N_CYCLE=1 DIFFUSION_BATCH_SIZE=1 \
scripts/nanobody/pretrain_protenix_v2_openfold_subset.sh
```

Evaluate the resulting checkpoint with `scripts/nanobody/run_inference_checkpoint.py` and `scripts/nanobody/evaluate_nanobody_predictions.py` as shown above.

## Template Team

Run template and MSA sweeps without changing weights:

```bash
python scripts/nanobody/run_template_sweep.py \
  --manifest ${PROTENIX_ROOT_DIR}/nanobody/manifests/nanobody_complex_val_2022h1.csv \
  --base-input-json ${PROTENIX_ROOT_DIR}/nanobody/inference_jsons/complex_val_2022h1.json \
  --checkpoint ${PROTENIX_ROOT_DIR}/checkpoint/protenix-v2.pt \
  --out-root outputs/nanobody/template_sweeps \
  --modes no_template,default_template,top1,top5 \
  --seeds 101,102,103
```

`no_template` removes template fields and runs with `--use_template false`. `default_template` leaves the input JSON unchanged and runs with `--use_template true`. `top1` and `top5` are recorded as best effort because Protenix can consume JSON templates, A3M, and HHR-style inputs whose filtering rules are format-specific; this wrapper records the limitation in the sweep summary instead of silently modifying unsupported files.

For validation/test sweeps, keep `--template-max-release-date 2022-01-01` unless there is a documented reason to change it. Do not use validation or test targets themselves as templates.

## Experiment Outputs

Each team should write comparable artifacts under:

```text
outputs/nanobody/experiments/<experiment_id>/
  command.sh
  config.yaml
  checkpoint_paths.txt
  validation_metrics.csv
  validation_summary.json
  test_metrics.csv
  experiment_summary.json
```

`experiment_summary.json` should contain:

```json
{
  "experiment_id": "posttrain_nanobody_v2_example",
  "team": "posttrain",
  "model_name": "protenix-v2",
  "base_checkpoint": "${PROTENIX_ROOT_DIR}/checkpoint/protenix-v2.pt",
  "output_checkpoint": "outputs/nanobody/posttrain/<run>/checkpoints/<step>_ema_0.999.pt",
  "train_sets": ["nanobody_single_train_pre2022", "nanobody_complex_train_pre2022"],
  "validation_sets": ["nanobody_single_val_2022h1", "nanobody_complex_val_2022h1"],
  "template_mode": "default_template",
  "mean_dockq": null,
  "dockq_success_rate_0_23": null,
  "median_cdrh3_rmsd_ca": null,
  "mean_cdrh3_rmsd_ca": null,
  "notes": ""
}
```

Use `scripts/nanobody/summarize_nanobody_results.py` to collect summaries.

## Optional Dependencies

DockQ enables DockQ, iRMSD, LRMSD, and fnat metrics for complexes. ANARCI/anarcii can be installed for future robust CDR-H3 numbering; explicit metadata CDR-H3 bounds remain the preferred reproducible source. Neither package is imported at module import time or required for the pipeline to run.

## Smoke Test

```bash
bash scripts/nanobody/smoke_nanobody_pipeline.sh
```

The default smoke creates a tiny local manifest from `examples/2lwu.cif`, prepares canonical files, copies the reference structure as a fake prediction, and verifies evaluator CSV/JSON output. Set `PROTENIX_SMOKE_RUN_TRAINING=true` only in an environment with Protenix runtime dependencies and training data caches installed.
