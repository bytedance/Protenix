# Training and Inference Instructions

This document provides detailed instructions for installing Protenix, performing inference, and training or fine-tuning the model.

## 🛠 Installation

### From PyPI (Stable Version)
```bash
pip3 install protenix
```

### From GitHub (Latest Development Version)
To install the absolute latest code directly from the repository:
```bash
pip3 install git+https://github.com/bytedance/Protenix.git
```

### From Local Source (For Developers)
```bash
git clone https://github.com/bytedance/Protenix.git
cd Protenix
pip3 install -e .
```

### Docker (Recommended for Training)
Check the detailed guide: [<u> Docker Installation</u>](./docker_installation.md).


## 🚀 Inference & CLI Usage

Protenix provides a unified CLI for structure prediction, data preprocessing, and database searching. If installed via `pip`, use the `protenix` command.

### CLI Commands Overview
| Command | Alias | Description |
|---------|-------|-------------|
| `predict` | `pred` | Perform model inference on JSON input(s). |
| `tojson` | `json` | Convert PDB or CIF files to Protenix-compatible JSON. |
| `msa` | `msa` | Generate Multiple Sequence Alignments (MSA) for proteins. |
| `msatemplate` | `mt` | Run sequential MSA and template search. |
| `inputprep` | `prep` | Full preprocessing: MSA, Template, and RNA MSA search. |

### 1. Data Conversion (`tojson`)
Convert structural files into the required JSON format.
```bash
# Convert PDB/CIF to JSON
protenix json --input ./examples/7pzb.pdb --out_dir ./output --altloc first

# Advanced: Specify assembly ID for biological assemblies
wget -P ./examples/ https://files.rcsb.org/download/7pzb.cif 
protenix json --input ./examples/7pzb.cif --out_dir ./output --altloc first
```

### 2. Input Preprocessing (`prep`, `mt`, `msa`)
Protenix requires MSA and template information for optimal accuracy.
```bash
# Full preprocessing (Protein MSA + Template + RNA MSA)
protenix prep --input examples/input.json --out_dir ./output

# Sequential Protein MSA and Template search
protenix mt --input examples/input.json --out_dir ./output

# Independent MSA search (supports JSON or Protein FASTA)
protenix msa --input examples/prot.fasta --out_dir ./output --msa_server_mode protenix
```

> **Note**: For `prep` and `mt`, you may need to specify paths to external databases (e.g., `--seqres_database_path`) and HMMER binaries if they are not in your system PATH.

### 3. Model Inference (`predict`)
Run the prediction engine with customizable configurations.
```bash
# Standard inference (using default model and parameters)
protenix pred -i examples/input.json -o ./output -s 101 -n  protenix_base_default_v1.0.0 --use_template true

# Standard inference with seeds from JSON
protenix pred -i  examples/examples_with_template/example_mgyp004658859411.json --use_seeds_in_json true

# Standard inference with RNA MSA, only support protenix_base_default_v1.0.0
protenix pred -i  examples/examples_with_rna_msa/example_9gmw_2.json --use_rna_msa true -n protenix_base_default_v1.0.0

# Inference using Protenix-Mini (faster, lightweight)
protenix pred --input examples/input.json --model_name protenix_mini_default_v0.5.0

# Customized inference: Disable MSA and use shared variable caching
protenix pred --input examples/input.json --use_msa false --enable_cache true
```

#### Key Inference Flags
- `--seeds`: Comma-separated list of random seeds (e.g., `101,102`).
- `--model_name`: Model variant selection (e.g., `protenix_base_default_v1.0.0`, `protenix_mini_default_v0.5.0`).
- `--use_default_params`: (Default: `true`) Automatically configures cycles and steps based on the selected model. Set to `false` to manually override `--cycle` and `--step`.
- `--use_msa` / `--use_template` / `--use_rna_msa`: (Default: `true`/`false`/`false`) Toggle specific features for inference.
- `--dtype`: Set data type to `bf16` (default) or `fp32`.
- `--trimul_kernel` / `--triatt_kernel`: Choose specialized kernels (e.g., `cuequivariance`, `triattention`) for hardware acceleration.
- `--enable_cache` / `--enable_fusion`: Enable memory/speed optimizations (recommended for GPU).

### Inference via Bash Script
Alternatively, use the provided demo script for automated runs:
```bash
bash inference_demo.sh <model_name> <input_json> <output_dir> <dtype> <use_msa>
```

Key arguments in `inference_demo.sh`:
* `model_name`: Name of the model to use for inference.
* `input_json_path`: Path to a JSON file that fully specifies the input structure.
* `dump_dir`: Directory where inference results will be saved.
* `dtype`: Data type used during inference. Supported options: `bf16` and `fp32`.
* `use_msa`: Whether to enable MSA features (default: true).
* `sample_diffusion.N_sample`: Number of samples to generate for each structure.
* `sample_diffusion.N_step`: Number of steps for the diffusion process (e.g., 200).
* `model.N_cycle`: Number of recycling steps.
* `use_template`: Whether to use structural templates (requires `templatesPath` in the input JSON).

> **Performance Tip**: By default, specialized CUDA kernels are enabled. For significant speedups on NVIDIA GPUs, follow the [**Kernels Setup Guide**](./kernels.md).


## 🧬 Training

### Preparing the datasets
To download the [wwPDB dataset](https://www.wwpdb.org/) and preprocessed training data, you need at least 1.5T disk space.

We recommend setting an environment variable `PROTENIX_ROOT_DIR` to specify the data directory. You can then download and extract the preprocessed databases using the following commands:

```bash
# Set your data root directory
export PROTENIX_ROOT_DIR=/path/to/your/data_root
mkdir -p $PROTENIX_ROOT_DIR

# Download and extract data components using the provided script
# Use --inference_only (default) or --full for training/finetuning
bash scripts/database/download_protenix_data.sh --full
```

Data hierarchy after extraction should be as follows:

  ```bash
  ├── common
  │   ├── clusters-by-entity-40.txt # cluster file for identity antibody
  │   ├── components.cif # ccd source file
  │   ├── components.cif.rdkit_mol.pkl # rdkit Mol object
  │   ├── obsolete_release_date.csv # release date of obsolete pdb
  │   ├── obsolete_to_successor.json # mapping from obsolete pdb to successor pdb
  │   ├── release_date_cache.json # cache of release date
  │   └── seq_to_pdb_index.json # mapping from sequence to directory ID
  ├── indices
  │   ├── posebusters_indices_mainchain_interface.csv # indices for posebusters dataset
  │   ├── recentPDB_low_homology_maxtoken1024_sample384_pdb_id.txt # indices for recentPDB evaluation
  │   ├── recentPDB_low_homology_maxtoken1536.csv # indices for recentPDB evaluation
  │   └── weightedPDB_indices_before_2021-09-30_wo_posebusters_resolution_below_9.csv.gz # indices for wwpdb training dataset
  ├── mmcif # raw mmcif data / template search database
  ├── mmcif_bioassembly # preprocessed wwPDB structural data/ traning cache
  ├── mmcif_msa_template # msa / template files for training
  ├── posebusters_bioassembly # preprocessed posebusters structural data
  ├── posebusters_mmcif # raw mmcif data
  ├── recentPDB_bioassembly # preprocessed recentPDB structural data
  └── search_database
      ├── nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta # NT-RNA database for RNA MSA search
      ├── pdb_seqres_2022_09_28.fasta # PDB seqres database for template search
      ├── rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta # Rfam database for RNA MSA search
      └── rnacentral_active_seq_id_90_cov_80_linclust.fasta # RNAcentral database for RNA MSA search
  ```

Data processing scripts have also been released. you can refer to [prepare_training_data.md](./prepare_training_data.md) for generating `{dataset}_bioassembly` and `indices`. And you can refer to [msa_template_pipeline.md](./msa_template_pipeline.md) for pipelines to get `mmcif_msa_template` and `seq_to_pdb_index.json`.

### Training demo
After installation and data preparation, you can run the following command to train the model from scratch:

```bash
bash train_demo.sh 
```

The script sets the `PYTHONPATH` and provides an option to set `PROTENIX_ROOT_DIR`. Key arguments in `train_demo.sh` are explained as follows:

* `model_name`: The configuration name for the model (e.g., `"protenix_base_default_v1.0.0"`).
* `dtype`: Data type used in training. Valid options include `"bf16"` (default) and `"fp32"`. 
  * `--dtype fp32`: Full FP32 precision.
  * `--dtype bf16`: BF16 Mixed precision. By default, `SampleDiffusion`, `ConfidenceHead`, `Mini-rollout`, and `Loss` parts remain in FP32 for stability.
* `ema_decay`: The decay rate for the Exponential Moving Average (EMA) of model weights, default is 0.999.
* `sample_diffusion.N_step`: Number of steps for the diffusion process during evaluation. Set to 20 in the demo for efficiency.
* `data.train_sets` / `data.test_sets`: Datasets used for training and evaluation (e.g., `weightedPDB_before2109_wopb_nometalc_0925`).
* `triangle_attention` / `triangle_multiplicative`: Optimization kernels for triangular modules. We recommend using `"cuequivariance"` for better performance.
* `use_wandb`: Whether to use Weights & Biases for logging.

The model also supports distributed training with PyTorch's `torchrun`:
```bash
torchrun --nproc_per_node=8 runner/train.py --model_name "protenix_base_default_v1.0.0" [OTHER_ARGS]
```


If you want to speed up training, see [<u> setting up kernels documentation </u>](./kernels.md).

### Finetune demo

If you want to fine-tune the model on a specific subset (e.g., a specific set of PDBs), you can use `finetune_demo.sh`. This script loads pretrained weights and restricts the training to a provided list of PDB IDs.

```bash
bash finetune_demo.sh
```

Key arguments for fine-tuning:
* `load_checkpoint_path` / `load_ema_checkpoint_path`: Paths to the pretrained model and its EMA weights.
* `data.<dataset_name>.base_info.pdb_list`: Path to a text file containing the specific PDB IDs to be used from the dataset.

For example, to fine-tune on a subset defined in `examples/finetune_subset.txt`:
```bash
--data.weightedPDB_before2109_wopb_nometalc_0925.base_info.pdb_list examples/finetune_subset.txt
```
The `subset.txt` should contain PDB IDs (one per line):
```text
6hvq
5mqc
5zin
```

## 📊 Training and Inference Cost

### Training Cost

The training configurations largely adhere to the specifications described in the [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w) paper. The table below summarizes the hyperparameters and performance metrics across different training stages:

| Hyperparameters | Initial Training | Fine-tuning Stage 1 | Fine-tuning Stage 2 | Fine-tuning Stage 3 |
| :--- | :---: | :---: | :---: | :---: |
| `train_crop_size` | 384 | 640 | 768 | 768 |
| `diffusion_batch_size` | 48 | 32 | 32 | 32 |
| `loss.weight.alpha_pae` | 0 | 0 | 0 | 1.0 |
| `loss.weight.alpha_bond` | 0 | 1.0 | 1.0 | 0 |
| `loss.weight.smooth_lddt` | 1.0 | 0 | 0 | 0 |
| `loss.weight.alpha_confidence` | 1e-4 | 1e-4 | 1e-4 | 1e-4 |
| `loss.weight.alpha_diffusion` | 4.0 | 4.0 | 4.0 | 0 |
| `loss.weight.alpha_distogram` | 0.03 | 0.03 | 0.03 | 0 |
| `train_confidence_only` | False | False | False | True |
| Throughput (A100-80G, s/step) | ~12 | ~30 | ~44 | ~13 |
| Peak GPU Memory (GB) | ~34 | ~35 | ~48 | ~24 |

We recommend training on high-performance GPUs such as NVIDIA A100 (80GB), H20, or H100. When using BF16 mixed-precision, the initial training stage can also be conducted on NVIDIA A800 (40GB). For GPUs with smaller memory capacities (e.g., NVIDIA A30), it is necessary to adjust the model configuration—such as reducing `model.pairformer.nblocks` or `diffusion_batch_size`—to ensure compatibility.

### Inference Cost

By default, the model performs inference in mixed-precision (BF16). However, the `SampleDiffusion` and `ConfidenceHead` modules are executed in full-precision (FP32) to maintain numerical stability and prediction accuracy.

The table below provides benchmark data for GPU memory utilization and inference latency across various input sizes.

| `N_token` | `N_atom` | Peak Memory (GB) | Latency (s) |
| :--- | :--- | :---: | :---: |
| 500 | 5000 | 6.1 | 17 |
| 1000 | 10000 | 18.2 | 59 |
| 2000 | 20000 | 66.6 | 226 |
| 3000 | 30000 | 60.8 | 935 |
| 4000 | 40000 | 78.1 | 1424 |

To mitigate potential Out-of-Memory (OOM) issues during large-scale inference, the inference script ([runner/inference.py](../runner/inference.py)) dynamically adjusts the precision for `SampleDiffusion` and `ConfidenceHead` based on the token count (`N_token`):
```python
def update_inference_configs(configs: Any, n_token: int) -> Any:
    """
    Adjust inference configurations based on the number of tokens to manage memory usage and prevent OOM.
    
    Args:
        configs (Any): Original configurations.
        n_token (int): Number of tokens in the sample.

    Returns:
        Any: Updated configurations.
    """
    if n_token > 3840:
        # Enable AMP for both modules to save memory for extremely large sequences
        configs.skip_amp.confidence_head = False
        configs.skip_amp.sample_diffusion = False
    elif n_token > 2560:
        # Enable AMP only for ConfidenceHead
        configs.skip_amp.confidence_head = False
        configs.skip_amp.sample_diffusion = True
    else:
        # Default: Disable AMP for both (run in FP32) to prioritize accuracy
        configs.skip_amp.confidence_head = True
        configs.skip_amp.sample_diffusion = True

    return configs
```
