### Run with Docker

1. Install Docker (with GPU Support)

    Ensure that Docker is installed and configured with GPU support. Follow these steps:
    *  Install [Docker](https://www.docker.com/) if not already installed.
    *  Install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) to enable GPU support.
    *  Verify the setup with:
        ```bash
        docker run --rm --gpus all nvidia/cuda:12.1.0-base-ubuntu22.04 nvidia-smi
        ```
        
2. Pull the Docker image, which was built based on this [Dockerfile](../Dockerfile)
    ```bash
    docker pull ai4s-cn-beijing.cr.volces.com/infra/protenix:v0.0.1
    ```

3. Clone this repository and `cd` into it
    ```bash
    git clone https://github.com/bytedance/protenix.git 
    cd ./protenix
    pip install -e .
    ```

4. Run Docker with an interactive shell
    ```bash
    docker run --gpus all -it -v $(pwd):/workspace -v /dev/shm:/dev/shm ai4s-cn-beijing.cr.volces.com/infra/protenix:v0.0.1 /bin/bash
    ```
  
  After running above commands, you’ll be inside the container’s environment and can execute commands as you would on a normal Linux terminal.