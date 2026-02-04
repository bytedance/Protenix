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

# ==============================================================================
# Protenix Model Inference Test Script
#
# Purpose:
#   This script provides usage examples for running inference with various
#   Protenix model versions and configurations.
#
# Arguments Summary (for 'protenix pred' or 'runner/inference.py'):
#   -i, --input (str):       [Required] Input JSON file or directory.
#   -o, --out_dir (str):     [Default: ./output] Output directory for results.
#   -s, --seeds (str):       [Default: 101] Inference seeds (e.g., "101,102").
#   -c, --cycle (int):       [Default: 10] Number of Pairformer cycles.
#   -p, --step (int):        [Default: 200] Number of diffusion steps.
#   -e, --sample (int):      [Default: 5] Samples per seed.
#   -d, --dtype (str):       [Default: bf16] Inference data type (bf16, fp32).
#   -n, --model_name (str):  [Default: protenix_base_default_v1.0.0] Model name.
#                            NOTE: protenix_base_default_v1.0.0 is the RECOMMENDED default.
#   --use_msa (bool):        Whether to use protein MSA features.
#   --use_default_params:    Auto-load recommended defaults for the model.
#   --trimul_kernel (str):   Triangle multiplicative kernel ('cuequivariance', 'torch').
#   --triatt_kernel (str):   Triangle attention kernel ('triattention', 'cuequivariance', etc.).
#   --use_template (bool):   Enable template features (v1.0.0+ only).
#   --use_rna_msa (bool):    Enable RNA MSA features (v1.0.0+ only).
#   --use_seeds_in_json:     Prioritize seeds defined in the input JSON.
#
# Available Models (Ref: configs/configs_model_type.py, docs/supported_models.md):
#   * protenix_base_default_v1.0.0:    [DEFAULT] Advanced model supporting Template & RNA MSA (Training Data Cutoff: 2021-09-30).
#   1. protenix_base_20250630_v1.0.0:  Latest model for practical scenarios (Training Data Cutoff: 2025-06-30).
#   2. protenix_base_default_v0.5.0:   Standard base model (Training Data Cutoff: 2021-09-30).
#   3. protenix_base_constraint_v0.5.0: Base model with constraint support (Training Data Cutoff: 2021-09-30).
#   4. protenix_mini_esm_v0.5.0:       Lightweight ESM-only model (no MSA) (Training Data Cutoff: 2021-09-30).
#   5. protenix_mini_ism_v0.5.0:       Lightweight ISM-only model (no MSA) (Training Data Cutoff: 2021-09-30).
#   6. protenix_mini_default_v0.5.0:   Standard lightweight model (Training Data Cutoff: 2021-09-30).
#   7. protenix_tiny_default_v0.5.0:   Ultra-lightweight model (Training Data Cutoff: 2021-09-30).
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 1: Running via Protenix CLI (protenix pred)
# ------------------------------------------------------------------------------

# IMPORTANT: Ensure environment variables (e.g., PROTENIX_ROOT_DIR, CUTLASS_PATH)
# are correctly set. Missing or incorrect settings may lead to execution errors.
# Specify your data root directory by uncommenting the following line.
# export PROTENIX_ROOT_DIR="/modify/to/your/data_root_dir"
# Environment Setup for deepspeed
# export CUTLASS_PATH=/opt/cutlass/

echo "Starting Section 1: CLI-based inference tests..."

# Example 1.1: Standard inference with Template support (v1.0.0)
protenix pred \
    -i examples/input.json \
    -o ./test_outputs/cmd/output_base_v1 \
    -s 101 \
    -n protenix_base_default_v1.0.0 \
    --use_template true \
    --use_default_params true


# Example 1.2: Inference using seeds defined in JSON
protenix pred \
    -i examples/examples_with_template/example_mgyp004658859411.json \
    -o ./test_outputs/cmd/output_base_v1 \
    -s 101 \
    -n protenix_base_default_v1.0.0 \
    --use_template true \
    --use_seeds_in_json true \
    --use_default_params true

# Example 1.3: RNA MSA support (v1.0.0 exclusive)
protenix pred \
    -i examples/examples_with_rna_msa/example_9gmw_2.json \
    -o ./test_outputs/cmd/output_base_v1 \
    -n protenix_base_default_v1.0.0 \
    --use_rna_msa true  \
    --use_default_params true

# Example 1.4: Latest model v1.0.0 with 2025-06-30 cutoff
protenix pred \
    -i examples/input.json \
    -o ./test_outputs/cmd/output_base_v1_20250630 \
    -s 101 \
    -n protenix_base_20250630_v1.0.0 \
    -c 4 \
    -p 20 \
    --use_template true

# Example 1.5: Base model v0.5.0 with precomputed MSA
protenix pred \
    -i examples/example.json \
    -o ./test_outputs/cmd/output_base \
    -s 101 \
    -c 4 \
    -p 20 \
    -n "protenix_base_default_v0.5.0" \
    --use_default_params true

# Example 1.6: Mini model with ESM features only
protenix pred \
    -i examples/example.json \
    -o ./test_outputs/cmd/output_mini_esm \
    -s 102 \
    -n "protenix_mini_esm_v0.5.0" \
    --use_default_params true

# Example 1.7: Mini model with ISM features only
protenix pred \
    -i examples/example.json \
    -o ./test_outputs/cmd/output_mini_ism \
    -s 103 \
    -n "protenix_mini_ism_v0.5.0" \
    --use_default_params true

# Example 1.8: Base constraint model
protenix pred \
    -i examples/example_constraint_msa.json \
    -o ./test_outputs/cmd/output_constraint \
    -s 104 \
    -n "protenix_base_constraint_v0.5.0" \
    --use_default_params true

# Example 1.9: Tiny default model
protenix pred \
    -i examples/example.json \
    -o ./test_outputs/cmd/output_tiny \
    -s 106 \
    -n "protenix_tiny_default_v0.5.0" \
    --use_default_params true


# ------------------------------------------------------------------------------
# Section 2: Running via Runner Script (runner/inference.py)
#
# IMPORTANT:
#   Direct script execution requires features (MSA, templates, RNA MSA, etc.)
#   to be pre-prepared in the input JSON. This mode is optimized for GPU-only
#   computation.
#   If features are NOT ready, please use the preprocessing command first:
#   Example: protenix prep --input examples/input.json --out_dir ./output
# ------------------------------------------------------------------------------

echo "Starting Section 2: Script-based inference tests..."
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Test 2.1: Base v1.0.0 with Template support
# Features: Template enabled, cuequivariance attention
N_sample=5
N_step=200
N_cycle=10
seed=103
input_json_path="./examples/examples_with_template/example_9fm7.json"
dump_dir="./test_outputs/sh/output_m_9fm7"
model_name="protenix_base_default_v1.0.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "cuequivariance" \
    --use_seeds_in_json true \
    --triangle_multiplicative "cuequivariance" \
    --use_template true

# Test 2.2: Latest model v1.0.0 with 2025-06-30 cutoff
N_sample=1
N_step=200
N_cycle=10
seed=101
input_json_path="./examples/input.json"
dump_dir="./test_outputs/sh/output_base_20250630"
model_name="protenix_base_20250630_v1.0.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "cuequivariance" \
    --triangle_multiplicative "cuequivariance" \
    --use_template true

# Test 2.3: Base v0.5.0 with triattention
N_sample=1
N_step=200
N_cycle=10
seed=101
input_json_path="./examples/example.json"
dump_dir="./test_outputs/sh/output_base"
model_name="protenix_base_default_v0.5.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "triattention" \
    --triangle_multiplicative "cuequivariance"

# Test 2.4: Mini ESM v0.5.0 with cuequivariance
N_sample=1
N_step=5
N_cycle=4
seed=101
input_json_path="./examples/example.json"
dump_dir="./test_outputs/sh/output_mini_esm"
model_name="protenix_mini_esm_v0.5.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "cuequivariance" \
    --triangle_multiplicative "cuequivariance"

# Test 2.5: Mini ISM v0.5.0 with deepspeed
N_sample=1
N_step=5
N_cycle=4
seed=101
input_json_path="./examples/example.json"
dump_dir="./test_outputs/sh/output_mini_ism"
model_name="protenix_mini_ism_v0.5.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "deepspeed" \
    --triangle_multiplicative "cuequivariance"

# Test 2.6: Base Constraint v0.5.0 with torch attention
N_sample=1
N_step=200
N_cycle=10
seed=101
input_json_path="./examples/example_constraint_msa.json"
dump_dir="./test_outputs/sh/output_constraint"
model_name="protenix_base_constraint_v0.5.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "torch" \
    --triangle_multiplicative "cuequivariance"

# Test 2.7: Mini Default v0.5.0 with torch attention/multiplicative
N_sample=1
N_step=5
N_cycle=4
seed=101
input_json_path="./examples/example.json"
dump_dir="./test_outputs/sh/output_mini"
model_name="protenix_mini_default_v0.5.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "torch" \
    --triangle_multiplicative "torch"

# Test 2.8: Tiny Default v0.5.0 with torch attention/multiplicative
N_sample=1
N_step=5
N_cycle=4
seed=101
input_json_path="./examples/example.json"
dump_dir="./test_outputs/sh/output_tiny"
model_name="protenix_tiny_default_v0.5.0"

python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "torch" \
    --triangle_multiplicative "torch"

echo "All inference tests completed."


# The following is a demo to use DDP for inference
# torchrun \
#     --nproc_per_node $NPROC \
#     --master_addr $WORKER_0_HOST \
#     --master_port $WORKER_0_PORT \
#     --node_rank=$ID \
#     --nnodes=$WORKER_NUM \
#     runner/inference.py \
#     --seeds ${seed} \
#     --dump_dir ${dump_dir} \
#     --input_json_path ${input_json_path} \
#     --model.N_cycle ${N_cycle} \
#     --sample_diffusion.N_sample ${N_sample} \
#     --sample_diffusion.N_step ${N_step}
