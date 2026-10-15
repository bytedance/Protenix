### Setting up kernels

- **Custom CUDA layernorm kernels** modified from [FastFold](https://github.com/hpcaitech/FastFold) and [Oneflow](https://github.com/Oneflow-Inc/oneflow) accelerate about 30%-50% during different training stages.
  - **`fast_layernorm` is used by default**, and no explicit setting is required.
  - If you wish to disable it and use the native PyTorch layernorm, you can set the environment variable to `torch`:
    ```bash
    export LAYERNORM_TYPE=torch
    ```
  - The kernels will be JIT-compiled when `fast_layernorm` is called for the first time.
  - On non-CUDA runtimes (for example Intel XPU), Protenix automatically falls back to the native PyTorch layernorm implementation.
  - CUDA keeps its existing fast-layernorm default; set `LAYERNORM_TYPE=torch` to use native PyTorch layernorm if the optional CUDA extension is unavailable.
- **Triangle_attention Kernel Options**
  The model supports four implementations for triangle attention, configurable in [configs_base.py](../configs/configs_base.py):
  ```python
  triangle_attention = "cuequivariance"  # or "triattention"/"deepspeed"/"torch"
  ```

  1. **TriAttention kernel**

      Custom kernel implementation from [protenix/model/tri_attention/](../protenix/model/tri_attention/)

  2. **cuEquivariance Kernel**
      Optimized implementation using NVIDIA's [cuEquivariance](https://github.com/NVIDIA/cuEquivariance) library.
  
  3. **[DeepSpeed DS4Sci_EvoformerAttention kernel](https://www.deepspeed.ai/tutorials/ds4sci_evoformerattention/)** is a memory-efficient attention kernel developed as part of a collaboration between OpenFold and the DeepSpeed4Science initiative.

      DS4Sci_EvoformerAttention is implemented based on [CUTLASS](https://github.com/NVIDIA/cutlass). If you use this feature, you need to clone the CUTLASS repository and specify the path to it in the environment variable `CUTLASS_PATH`. The [Dockerfile](Dockerfile) already includes this setting:
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

      This backend depends on NVIDIA CUTLASS kernels and is not available on Intel XPU.

- **Triangle_multiplicative Kernel Options**
    
    The Triangle Multiplicative operation supports two implementations, configurable in [configs_base.py](../configs/configs_base.py):
    ```python
    triangle_multiplicative = "cuequivariance"  # or "torch"
    ```

    1. cuEquivariance Kernel (Default)
    Optimized implementation using NVIDIA's [cuEquivariance](https://github.com/NVIDIA/cuEquivariance) library.

    2. Torch Native
    Standard PyTorch implementation.

- **Intel XPU note**
  - Install the portable dependency set from [requirements-xpu.txt](../requirements-xpu.txt) together with the PyTorch XPU wheel.
  - For single-device XPU inference, omitted triangle-kernel options default to `torch`. Explicit non-`torch` selections are rejected; CUDA retains its existing `cuequivariance` defaults.
