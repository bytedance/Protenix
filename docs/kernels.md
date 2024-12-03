### Setting up kernels

- **Custom CUDA layernorm kernels** modified from [FastFold](https://github.com/hpcaitech/FastFold) and [Oneflow](https://github.com/Oneflow-Inc/oneflow) accelerate about 30%-50% during different training stages. To use this feature, run the following command:
  ```bash
  export LAYERNORM_TYPE=fast_layernorm
  ```
  If the environment variable `LAYERNORM_TYPE` is set to `fast_layernorm`, the model will employ the layernorm we have developed; otherwise, the naive PyTorch layernorm will be adopted. The kernels will be compiled when `fast_layernorm` is called for the first time.
- **[DeepSpeed DS4Sci_EvoformerAttention kernel](https://www.deepspeed.ai/tutorials/ds4sci_evoformerattention/)** is a memory-efficient attention kernel developed as part of a collaboration between OpenFold and the DeepSpeed4Science initiative. To use this feature, simply pass: 
  ```bash
  --use_deepspeed_evo_attention true
  ```
  into the command line. DS4Sci_EvoformerAttention is implemented based on [CUTLASS](https://github.com/NVIDIA/cutlass). You need to clone the CUTLASS repository and specify the path to it in the environment variable CUTLASS_PATH. The [Dockerfile](Dockerfile) has already include this setting:
  ```bash
  RUN git clone -b v3.5.1 https://github.com/NVIDIA/cutlass.git  /opt/cutlass
  ENV CUTLASS_PATH=/opt/cutlass
  ```
  If you set up `Protenix` by `pip`, you can set environment variable `CUTLASS_PATH` as follows:

  ```bash
  git clone -b v3.5.1 https://github.com/NVIDIA/cutlass.git  /path/to/cutlass
  export CUTLASS_PATH=/path/to/cutlass
  ```

  The kernels will be compiled when DS4Sci_EvoformerAttention is called for the first time.