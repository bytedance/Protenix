FROM ai4s-cn-beijing.cr.volces.com/pytorch-mirror/pytorch:2.3.1-cuda12.1-cudnn8-devel

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    g++ \
    gcc \
    libc6-dev \
    make zlib1g zlib1g-dev \
    git git-lfs expect zsh vim wget curl unzip zip cmake cmake-curses-gui libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

RUN apt update && apt -y install postgresql

RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        hmmer cmake cmake-curses-gui   \
    && git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
    && mkdir /tmp/hh-suite/build \
    && cd /tmp/hh-suite/build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
    && make -j 32 && make install \
    && ln -s /opt/hhsuite/bin/* /usr/bin \
    && cd - \
    && rm -rf /tmp/hh-suite

RUN apt-get install -yq --no-install-recommends iproute2 curl
# Add PIP Package
RUN pip3 --no-cache-dir install \
    scipy \
    ml_collections \
    tqdm \
    pandas \
    optree \
    rdkit

# Add openfold dependency
RUN pip3 --no-cache-dir install \
    biopython==1.83 \
    modelcif==0.7 
  
# Add datapipeline dependency
RUN pip3 --no-cache-dir install \
    biotite==1.0.1 \
    scikit-learn \
    scikit-learn-extra \
    deepspeed>=0.15.1 \
    protobuf==3.20.2 tos icecream ipdb wandb numpy==1.26.3 matplotlib==3.9.2 ipywidgets py3Dmol

# For H20 compatibility
RUN pip3 install --no-cache-dir nvidia-cublas-cu12==12.4.5.8 --no-deps
RUN git clone  -b v3.5.1 https://github.com/NVIDIA/cutlass.git  /opt/cutlass
ENV CUTLASS_PATH=/opt/cutlass
