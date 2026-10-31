# Protenix on Ascend NPU

This document describes how to install and run Protenix training and inference on Ascend NPU.

## Repository layout

```
Protenix/
├── npu_installation/
│   ├── npu_docker/
│   │   ├── Dockerfile              # NPU (950) base image
│   │   ├── build.sh                # wrapper around docker build
│   │   ├── check_npu_env.sh        # in-container environment self check
│   │   └── docker-entrypoint.sh    # container entrypoint (env setup + --selftest)
│   └── patches/
│       ├── fag_dpse_1nss.patch     # dpse support for ops-transformer (external repo)
│       └── Protenix.patch          # NPU adaptation for this repo (unmodified upstream only)
└── migrate_to_npu/
    ├── README.md                   # short notes for the NPU assets
    ├── env_npu.sh                  # one-shot NPU environment (CANN / toggles)
    ├── requirements_npu.txt        # NPU python dependencies
    └── test/
        ├── inference_demo.sh       # NPU inference example
        └── train_demo.sh           # NPU training example
```

`migrate_to_npu/` is created by `npu_installation/patches/Protenix.patch`, together with the
NPU-related changes to the model code, so it appears in the tree once the patch is applied.

## Component versions

```shell
Ascend NPU: 910B / 910C / 950
hdk: 26.0.RC1
cann: 9.3.0
python: 3.11
torch: 2.7.1
torch-npu: 2.7.1.post6
```

## Environment setup

The environment can be installed either manually or with the provided Dockerfile.

### Manual installation

1. Activate the CANN environment.

    ```shell
    # The deault path is for reference, please modify it to your actual path.
    source /usr/local/Ascend/ascend-toolkit/set_env.sh
    ```

2. Create a conda environment and install the Python dependencies.

    ```shell
    conda create --name protenix python=3.11
    conda activate protenix

    pip install torch==2.7.1
    pip install torch_npu==2.7.1.post6

    # Optional third-party acceleration library: 910B/910C -> mx-driving, 950 -> mx-driving-a5
    pip install mx-driving
    pip install mx-driving-a5
    ```

3. Apply the NPU adaptation patch and install the NPU requirements.
    Run this from the repository root; the patch adds `migrate_to_npu/` and the
    NPU-related changes to the model code.

    ```shell
    git apply --ignore-whitespace npu_installation/patches/Protenix.patch
    pip install -r migrate_to_npu/requirements_npu.txt
    ```

4. Install the CANN patch package for dpse (required for training).
    Pick the patch package based on NPU type and CANN version.

    ```shell
    1) NPU: 910B/910C  CANN: 9.3.0
    git clone https://gitcode.com/cann/ops-transformer.git && cd ops-transformer
    git fetch origin refs/merge-requests/2576/head:pr-2576
    git checkout pr-2576
    cd ${git_clone_path}/experimental/attention
    mv flash_attention_score_grad_dpse flash_attention_score_grad
    cd ../..
    # ascend910b is an example; change it to your actual SoC
    bash build.sh --pkg --experimental --soc=ascend910b --ops=flash_attention_score_grad
    ./build_out/cann-ops-transformer-*linux*.run

    2) NPU: 950  CANN: 9.3.0
    git clone https://gitcode.com/cann/ops-transformer.git && cd ops-transformer
    git checkout dc004cb8b
    git apply ${repo_root}/npu_installation/patches/fag_dpse_1nss.patch
    bash build.sh --pkg --experimental --soc=ascend950 --ops=flash_attention_score_grad --vendor_name=fagtest -j64
    ./build_out/cann-ops-transformer-fagtest_linux-aarch64.run
    ```

### Installation with the Dockerfile

The Dockerfile only supports 950. Training requires CANN 9.3.0; inference requires
8.3.RC1 or later. Download the base image from the Ascend website:
https://www.hiascend.com/developer/ascendhub

The patch under `npu_installation/patches/` is applied during the build, so there is no
need to apply it manually.

1. Build the image.

    ```shell
    BASE=<your-Base-image> IMAGE=protenix_env bash npu_installation/npu_docker/build.sh
    ```

2. If the base image ships an older CANN, let the build install the 9.3.0 packages.
    The installers are consumed from `canns/` in the build context.

    ```shell
    mkdir canns
    ln {install_path}/Ascend-cann-toolkit_*.run  canns/
    ln {install_path}/Ascend-cann-950-ops_*.run  canns/
    BASE=<your-Base-image> IMAGE=protenix_env bash npu_installation/npu_docker/build.sh
    ```

3. Start the container. Change the data path below to your own path.

    ```shell
    docker run -itd --ipc=host \
      --name protenix_docker \
      --network=host \
      --privileged -u=root \
      --device=/dev/davinci0 \
      --device=/dev/davinci1 \
      --device=/dev/davinci2 \
      --device=/dev/davinci3 \
      --device=/dev/davinci4 \
      --device=/dev/davinci5 \
      --device=/dev/davinci6 \
      --device=/dev/davinci7 \
      --device=/dev/davinci_manager \
      --device=/dev/devmm_svm \
      --device=/dev/hisi_hdc \
      -v /usr/local/sbin/npu-smi:/usr/local/sbin/npu-smi \
      -v /usr/local/Ascend/driver:/usr/local/Ascend/driver \
      -v /usr/bin/hccn_tool:/usr/bin/hccn_tool \
      -v /usr/lib64/liburma.so.0:/usr/lib/aarch64-linux-gnu/liburma.so.0:ro \
      -v /usr/lib64/liburma_common.so.0:/usr/lib/aarch64-linux-gnu/liburma_common.so.0:ro \
      -v /your/data:/data \
      protenix_env \
      /bin/bash
    ```

    The container uses `/data/Protenix` as the data root, so `-v /your/data:/data`
    exposes `/your/data/Protenix` at that location.

4. Check the environment inside the container.

    ```shell
    docker exec -it protenix_docker check_npu_env
    ```

## Dataset preparation

```shell
a. Go to the Protenix directory
cd Protenix/scripts/database
b. Set the data root directory and download the data
export PROTENIX_ROOT_DIR=/data/Protenix
# Download the inference dataset only
bash download_protenix_data.sh
# Download the full dataset (training + inference)
bash download_protenix_data.sh --full
```

## Weights preparation

```shell
The checkpoint path is detected automatically when running inference; download the
weights yourself and place them under $PROTENIX_ROOT_DIR/checkpoint/.
```

## Inference

```shell
cd Protenix
bash ./migrate_to_npu/test/inference_demo.sh
```

## Training

```shell
cd Protenix
bash ./migrate_to_npu/test/train_demo.sh
```

## Training benchmark

|     Chip      | Cards | Average step time (s) |
| :-----------: | :---: | :-------------------: |
|     A100      |  1p   |          22           |
|     910C      |  1p   |          25           |
|     950DT     |  1p   |          12           |

*Measured with sequence length 593 aa and fixed N_cycle=4; the reported value is the step-14 time.

## Inference benchmark

|     Chip      | Sequence length |  Configuration    | Forward time (s) |
| :-----------: | :-------------: | :---------------: | :--------------: |
|     A100      |       693       | cuequivariance on |      25.65       |
|     950DT     |       693       |      default      |      17.21       |

*Model weights: base_v1.0.0; test input: example input.json.

## Q&A

```shell
1. During training you may hit errors related to unfold_grad or cdist:
   This is caused by a CANN version mismatch. Switch to CANN 8.3.RC1 or 9.3.0.

2. With the Dockerfile setup, running the model may fail with "507033" or "TsdOpen failed":
   This is a base-image problem, not a device problem: changing ASCEND_RT_VISIBLE_DEVICES does not
   help and npu-smi still reports the cards as healthy. It is also unrelated to the CANN version
   (switching CANN or torch_npu inside the same image does not help either).
   Workaround: use the manual installation path on the bare metal, or rebuild on a different base image.
```
