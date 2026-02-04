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
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
# fast_layernorm is used by default, no need to set explicitly. Set LAYERNORM_TYPE=torch to disable.
# export LAYERNORM_TYPE=fast_layernorm
# Kernel options:
# - triangle_attention: supports 'triattention', 'cuequivariance', 'deepspeed', 'torch'
# - triangle_multiplicative: supports 'cuequivariance', 'torch'

# Specify your data root directory by uncommenting the following line.
# export PROTENIX_ROOT_DIR="/modify/to/your/data_root_dir"

python3 ./runner/train.py \
--run_name protenix_train \
--model_name "protenix_base_default_v1.0.0" \
--seed 42 \
--base_dir ./output \
--dtype bf16 \
--project protenix \
--use_wandb false \
--diffusion_batch_size 48 \
--eval_interval 400 \
--log_interval 50 \
--checkpoint_interval 400 \
--ema_decay 0.999 \
--train_crop_size 384 \
--max_steps 100000 \
--warmup_steps 2000 \
--lr 0.001 \
--model.N_cycle 4 \
--sample_diffusion.N_step 20 \
--triangle_attention "cuequivariance" \
--triangle_multiplicative "cuequivariance" \
--data.train_sets weightedPDB_before2109_wopb_nometalc_0925 \
--data.test_sets recentPDB_1536_sample384_0925,posebusters_0925 \
--data.posebusters_0925.base_info.max_n_token 768
