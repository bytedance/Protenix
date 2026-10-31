# syntax=docker/dockerfile:1.7
# =============================================================================
# Protenix on Ascend 950 (A5)
#
# Mirrors the manual installation described in docs/npu_installation.md
#
# The build context is the repository root (the directory that contains
# npu_installation/), so every COPY below is relative to it.
#
# Optional context directory:
#   canns/    CANN installers (*.run). Empty means "use the CANN that ships with
#             the base image".
#
# Build arguments:
#   BASE               base image
#   PY_ENV             conda env inside the base image (python 3.11 + torch + torch_npu)
#   ASCEND_SET_ENV     path of the CANN set_env.sh (point it to a newer CANN if the
#                      base image already contains several CANN versions)
#   TORCH_NPU_SPEC     override the torch_npu version (must match the CANN version)
#   KEEP_BASE_TORCH    1 = keep the torch/torch_npu from the base image (avoids the
#                      downgrade performed by the requirements file)
#   INSTALL_MX_DRIVING 1 = install mx_driving (A5: mx-driving-a5, A2/A3: mx-driving)
# =============================================================================

ARG BASE
ARG PY_ENV=torch2.7.1
ARG ASCEND_SET_ENV=/usr/local/Ascend/ascend-toolkit/set_env.sh

# =============================================================================
# Stage 0 (optional): install a newer CANN into the image
#
# Put both installers into canns/ of the build context (same version):
#   Ascend-cann-toolkit_<ver>_linux-aarch64.run
#   Ascend-cann-<soc>-ops_<ver>_linux-aarch64.run      e.g. Ascend-cann-950-ops_9.3.0_...
# When canns/ is empty this stage is skipped and the CANN of the base image is used.
# After the installation ascend-toolkit/latest points to the new version, so the
# default value of ASCEND_SET_ENV keeps working.
# =============================================================================
FROM ${BASE} AS cannbase
SHELL ["/bin/bash", "-lc"]

# build.sh creates canns/ even when it is empty, so this COPY always succeeds
COPY canns/ /tmp/cann/
RUN if ls /tmp/cann/*.run >/dev/null 2>&1; then \
      set -e; \
      chmod +x /tmp/cann/*.run; \
      for f in /tmp/cann/Ascend-cann-toolkit_*.run; do \
        echo "[info] installing $(basename "$f")"; \
        "$f" --install --quiet --install-path=/usr/local/Ascend --force; \
      done; \
      for f in /tmp/cann/Ascend-cann-*-ops_*.run; do \
        echo "[info] installing $(basename "$f")"; \
        "$f" --install --quiet --install-path=/usr/local/Ascend --force; \
      done; \
      rm -rf /tmp/cann; \
      echo "[ OK ] CANN installed, latest -> $(readlink -f /usr/local/Ascend/ascend-toolkit/latest)"; \
    else \
      echo "[skip] no .run installer in canns/, using the CANN of the base image"; \
      rm -rf /tmp/cann; \
    fi

# =============================================================================
# Stage 1: builder -- needs the CANN toolchain, only produces the op package
# =============================================================================
FROM cannbase AS builder
SHELL ["/bin/bash", "-lc"]

ARG PY_ENV=torch2.7.1
ARG ASCEND_SET_ENV
ARG OPS_TRANSFORMER_COMMIT=dc004cb8b
ENV PY=/opt/conda/envs/${PY_ENV}/bin/python

# Build prerequisites. The cmake logic of ops-transformer expects a GNU patch
# binary, which some base images do not ship, so install it when missing.
RUN source ${ASCEND_SET_ENV} \
 && ${PY} -m pip install --no-cache-dir ninja wheel

RUN if ! command -v patch >/dev/null 2>&1; then \
      echo "[info] GNU patch not found, installing it"; \
      if command -v apt-get >/dev/null 2>&1; then \
          apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends patch; \
      elif command -v dnf >/dev/null 2>&1; then dnf install -y patch; \
      elif command -v yum >/dev/null 2>&1; then yum install -y patch; \
      else echo "[ERROR] no package manager available, install patch manually"; exit 1; fi; \
    fi \
 && command -v patch && command -v cmake && command -v g++

# ---------- dpse patch package (1NSS broadcast pse support) ----------
COPY npu_installation/patches/fag_dpse_1nss.patch /work/patches/
RUN source ${ASCEND_SET_ENV} \
 && git clone https://gitcode.com/cann/ops-transformer.git /work/ops-transformer \
 && cd /work/ops-transformer \
 && git checkout ${OPS_TRANSFORMER_COMMIT} \
 && git apply /work/patches/fag_dpse_1nss.patch \
 && bash build.sh --pkg --experimental --soc=ascend950 \
      --ops=flash_attention_score_grad --vendor_name=fagtest -j"$(nproc)" \
 && ./build_out/cann-ops-transformer-fagtest_linux-aarch64.run \
      --quiet --install-path=/opt/ops/fag

# =============================================================================
# Stage 2: runtime -- only what is needed at run time
# =============================================================================
FROM cannbase AS runtime
SHELL ["/bin/bash", "-lc"]

ARG PY_ENV=torch2.7.1
ARG ASCEND_SET_ENV
ARG PROTENIX_COMMIT=HEAD
# Source of the Protenix code. github.com is not reachable from inside the build
# container in some networks, so the gitcode mirror is the default.
#   alternatives: https://gitee.com/mirrors/Protenix.git
#                 https://ghfast.top/https://github.com/bytedance/Protenix.git
ARG PROTENIX_GIT_URL=https://gitcode.com/gh_mirrors/pr/Protenix.git
ARG KEEP_BASE_TORCH=1
ARG TORCH_NPU_SPEC=""
ARG INSTALL_MX_DRIVING=1
ARG MX_DRIVING_SPEC=mx-driving-a5
ARG PIP_INDEX_URL=""
ENV PY=/opt/conda/envs/${PY_ENV}/bin/python \
    PROTENIX_ROOT_DIR=/data/Protenix \
    PROTENIX_CONDA_ENV=${PY_ENV} \
    PIP_NO_CACHE_DIR=1 \
    LANG=C.UTF-8 \
    TZ=Asia/Shanghai \
    LAYERNORM_TYPE=torch

# 1) Op package produced by the builder stage (libraries / kernels / headers only)
COPY --from=builder /opt/ops /opt/ops

# 2) Model code + NPU adaptation patch + python dependencies
#    - Dependencies come from migrate_to_npu/requirements_npu.txt, which is the full
#      list: upstream requirements.txt minus the CUDA-only packages, plus torch_npu
#      and the NPU specific additions.
#    - KEEP_BASE_TORCH=1 skips torch/torch_npu: the base image already provides a
#      pair that matches its CANN, and reinstalling would downgrade it.
COPY npu_installation/patches/Protenix.patch /work/patches/
# The base image may already contain /work/Protenix (rebuild on top of this image),
# so remove it first. cd out of it: the WORKDIR of such a base image is /work/Protenix.
RUN cd / \
 && rm -rf /work/Protenix \
 && git -c advice.detachedHead=false clone ${PROTENIX_GIT_URL} /work/Protenix \
 && cd /work/Protenix \
 && git -c advice.detachedHead=false checkout ${PROTENIX_COMMIT} \
 && git apply --ignore-whitespace /work/patches/Protenix.patch \
 && if [ "${KEEP_BASE_TORCH}" = "1" ]; then \
        grep -vE '^(torch|torch_npu|torchvision|torchaudio)([=<>!~ ]|$)' \
             ./migrate_to_npu/requirements_npu.txt > /tmp/req.txt; \
    else cp ./migrate_to_npu/requirements_npu.txt /tmp/req.txt; fi \
 && echo "=== dependencies to install ===" && cat /tmp/req.txt \
 && if [ -n "${PIP_INDEX_URL}" ]; then ${PY} -m pip install -r /tmp/req.txt -i "${PIP_INDEX_URL}"; \
    else ${PY} -m pip install -r /tmp/req.txt; fi \
 && rm -f /tmp/req.txt

# 3) torch_npu override (must match the CANN version installed above)
RUN if [ -n "${TORCH_NPU_SPEC}" ]; then ${PY} -m pip install --no-cache-dir "${TORCH_NPU_SPEC}"; fi

# 3.5) Optional acceleration library mx_driving (scatter)
#      A5 -> mx-driving-a5, A2/A3 -> mx-driving. The code falls back to torch when
#      the module is missing, so a failure here is only a warning.
RUN if [ "${INSTALL_MX_DRIVING}" = "1" ]; then \
      echo "[info] installing ${MX_DRIVING_SPEC} (A5: mx-driving-a5, A2/A3: mx-driving)"; \
      if [ -n "${PIP_INDEX_URL}" ]; then \
          ${PY} -m pip install --no-cache-dir -i "${PIP_INDEX_URL}" "${MX_DRIVING_SPEC}" \
          && ${PY} -c "import mx_driving; print('[ OK ] mx_driving', getattr(mx_driving,'__version__','?'))"; \
      else \
          ${PY} -m pip install --no-cache-dir "${MX_DRIVING_SPEC}" \
          && ${PY} -c "import mx_driving; print('[ OK ] mx_driving', getattr(mx_driving,'__version__','?'))"; \
      fi || echo "[WARN] mx_driving not installed, scatter falls back to torch (no functional impact)"; \
    else \
      echo "[skip] INSTALL_MX_DRIVING=0"; \
    fi

# 4) Unified environment script + entrypoint + self check
#    /etc/profile.d/ptx-npu.sh gives `docker exec ... bash -l` the same environment:
#    docker exec does not go through the ENTRYPOINT, so without it the shell would
#    land in the base conda env, which has no torch.
RUN { \
      printf '%s\n' "[ -f ${ASCEND_SET_ENV} ] && source ${ASCEND_SET_ENV}"; \
      printf '%s\n' 'export ASCEND_CUSTOM_OPP_PATH="${ASCEND_CUSTOM_OPP_PATH:-}"'; \
      printf '%s\n' 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"'; \
      printf '%s\n' 'for _f in /opt/ops/*/vendors/*/bin/set_env.bash; do [ -f "$_f" ] && source "$_f"; done'; \
      printf '%s\n' "[ -f /opt/conda/etc/profile.d/conda.sh ] && source /opt/conda/etc/profile.d/conda.sh && conda activate ${PY_ENV} 2>/dev/null || true"; \
    } > /etc/profile.d/ptx-npu.sh

COPY npu_installation/npu_docker/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
COPY npu_installation/npu_docker/check_npu_env.sh     /usr/local/bin/check_npu_env
RUN chmod +x /usr/local/bin/docker-entrypoint.sh /usr/local/bin/check_npu_env \
 && chmod +x /etc/profile.d/ptx-npu.sh

WORKDIR /work/Protenix
LABEL org.opencontainers.image.title="protenix-npu" \
      org.opencontainers.image.description="Protenix on Ascend 950 (A5)" \
      org.opencontainers.image.base.name="${BASE}" \
      protenix.commit="${PROTENIX_COMMIT}" \
      ops.transformer.commit="dc004cb8b" \
      ops.dpse.patch="fag_dpse_1nss.patch (1NSS dpse support)"

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["--selftest"]
