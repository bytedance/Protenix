### Setting up kernels

- **Custom CUDA layernorm kernels** modified from [FastFold](https://github.com/hpcaitech/FastFold) and [Oneflow](https://github.com/Oneflow-Inc/oneflow) accelerate about 30%-50% during different training stages.
  - **`fast_layernorm` is used by default**, and no explicit setting is required.
  - If you wish to disable it and use the native PyTorch layernorm, you can set the environment variable to `torch`:
    ```bash
    export LAYERNORM_TYPE=torch
    ```
  - The kernels will be JIT-compiled when `fast_layernorm` is called for the first time.
- **Triangle_attention Kernel Options**
  The model supports four implementations for triangle attention, configurable in [configs_base.py](../configs/configs_base.py):
  ```python
  triangle_attention = "cuequivariance"  # or "triattention"/"deepspeed"/"torch"
  ```

  1. **TriAttention kernel(Default)**

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

- **Triangle_multiplicative Kernel Options**
    
    The Triangle Multiplicative operation supports three implementations, configurable in [configs_base.py](../configs/configs_base.py):
    ```python
    triangle_multiplicative = "cuequivariance"  # or "torch" / "fast_trimul"
    ```

    1. cuEquivariance Kernel (Default)
    Optimized implementation using NVIDIA's [cuEquivariance](https://github.com/NVIDIA/cuEquivariance) library.

    2. Torch Native
    Standard PyTorch implementation.

    3. fast_trimul (optional, third-party)
    Fused fp16 [CuTe DSL](https://github.com/tiagomonteiro0715/fast_trimul) kernel. **Inference only** (its
    backward is a correct but slow torch recompute), and it automatically falls back to the
    torch path when unavailable or during training. Install with:
    ```bash
    pip install "protenix[fast_trimul]"   # brings in fast_trimul + cuda-python<13
    ```
    The fast path targets NVIDIA Ampere (sm80) with `N` divisible by 8; other
    hardware/shapes/dtypes use fast_trimul's built-in pure-torch fallback.

