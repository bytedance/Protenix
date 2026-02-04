### Run with Docker

1. Install Docker (with GPU Support)

    Ensure that Docker is installed and configured with GPU support. Follow these steps:
    *   Install [Docker](https://www.docker.com/) if not already installed.
    *   Install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) to enable GPU support.
    *   Verify the setup with (using a version close to our environment):
        ```bash
        docker run --rm --gpus all nvidia/cuda:12.6.3-base-ubuntu22.04 nvidia-smi
        ```

2. Pull the Docker image
    The image contains all necessary dependencies (PyTorch, HMMER, Kalign, CUTLASS, etc.), but does not include the Protenix source code by default.
    ```bash
    docker pull ai4s-share-public-cn-beijing.cr.volces.com/release/protenix:1.0.0.4
    ```

3. Clone this repository
    ```bash
    git clone https://github.com/bytedance/protenix.git 
    cd ./protenix
    ```

4. Run Docker with an interactive shell
    Mount the current directory to `/app` inside the container. If you have external data or weights (e.g., in `/root/protenix`), consider mounting them as well.
    ```bash
    docker run --gpus all -it \
        -v "$(pwd)":/app \
        -v /dev/shm:/dev/shm \
        ai4s-share-public-cn-beijing.cr.volces.com/release/protenix:1.0.0.4 \
        /bin/bash
    ```

5. Install Protenix and Verify
    Once inside the container, install Protenix in editable mode and verify the installation:
    ```bash
    cd /app
    pip install -e .
    
    # Verify the installation by checking the help message
    protenix --help
    ```

After completing these steps, you can proceed with inference or training. See [Inference Guide](infer_json_format.md) for more details.
