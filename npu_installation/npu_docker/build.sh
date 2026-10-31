#!/bin/bash
# Build the Protenix NPU image.
#
# Usage (run from the repository root):
#   bash npu_installation/npu_docker/build.sh
#
# Environment variables:
#   BASE=<base image>       required, e.g. BASE=<repository:tag>
#   IMAGE=<image name>      default: protenix_env:1.0.0
#   PY_ENV=torch2.7.1       conda env inside the base image (python 3.11 + torch + torch_npu)
#   ASCEND_SET_ENV=...      path of the CANN set_env.sh (default: ascend-toolkit/set_env.sh)
#   TORCH_NPU_SPEC=...      override the torch_npu version (must match the CANN version)
#   KEEP_BASE_TORCH=1|0     keep the torch/torch_npu of the base image (default 1)
#   INSTALL_MX_DRIVING=1|0  install mx_driving (default 1)
#   MX_DRIVING_SPEC=...     package name: 950 -> mx-driving-a5, 910B/910C -> mx-driving
#   PIP_INDEX_URL=<url>     pip mirror for restricted networks
#   TARGET=builder|runtime  build stage (default runtime)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
# The build context is the repository root: <repo>/npu_installation/npu_docker -> <repo>
REPO="$(cd "$HERE/../.." && pwd)"

BASE="${BASE:?BASE must be set to an Ascend base image, e.g. BASE=<repository:tag>}"
PY_ENV="${PY_ENV:-torch2.7.1}"
IMAGE="${IMAGE:-protenix_env:1.0.0}"
TARGET="${TARGET:-runtime}"
PIP_INDEX_URL="${PIP_INDEX_URL:-}"
KEEP_BASE_TORCH="${KEEP_BASE_TORCH:-1}"
TORCH_NPU_SPEC="${TORCH_NPU_SPEC:-}"
ASCEND_SET_ENV="${ASCEND_SET_ENV:-/usr/local/Ascend/ascend-toolkit/set_env.sh}"
# Source of the Protenix code. github.com may be unreachable during the build, so
# the gitcode mirror is used by default.
PROTENIX_GIT_URL="${PROTENIX_GIT_URL:-https://gitcode.com/gh_mirrors/pr/Protenix.git}"
# Optional scatter acceleration library: A5 -> mx-driving-a5, A2/A3 -> mx-driving
INSTALL_MX_DRIVING="${INSTALL_MX_DRIVING:-1}"
MX_DRIVING_SPEC="${MX_DRIVING_SPEC:-mx-driving-a5}"

# --- BASE must be "repository:tag" and must exist locally; otherwise docker tries
#     to pull it from Docker Hub, which times out in offline environments.
case "$BASE" in
  *:*) ;;
  *) echo "[ERROR] BASE is missing the repository prefix, expected something like drivingsdk:<tag>; got: $BASE"; exit 1;;
esac
if ! docker image inspect "$BASE" >/dev/null 2>&1; then
    echo "[ERROR] image $BASE does not exist locally"
    docker images --format "  {{.Repository}}:{{.Tag}}" | grep -iE "drivingsdk|cann" | head -10
    exit 1
fi

CTX="$REPO"
DOCKERFILE="$HERE/Dockerfile"

# canns/ must exist in the build context because the Dockerfile copies it:
#   contains Ascend-cann-toolkit_*.run + Ascend-cann-<soc>-ops_*.run to install a
#   newer CANN into the image; an empty directory means "use the CANN of the base
#   image". The installers are large, so link them instead of copying if possible.
mkdir -p "$REPO/canns"
[ -z "$(ls -A "$REPO/canns")" ] && touch "$REPO/canns/.gitkeep"

BUILDFLAGS=(
  --build-arg BASE="$BASE"
  --build-arg PY_ENV="$PY_ENV"
  --build-arg ASCEND_SET_ENV="$ASCEND_SET_ENV"
  --build-arg TORCH_NPU_SPEC="$TORCH_NPU_SPEC"
  --build-arg KEEP_BASE_TORCH="$KEEP_BASE_TORCH"
  --build-arg PROTENIX_GIT_URL="$PROTENIX_GIT_URL"
  --build-arg INSTALL_MX_DRIVING="$INSTALL_MX_DRIVING"
  --build-arg MX_DRIVING_SPEC="$MX_DRIVING_SPEC"
  --network host
)

echo "=============================================================="
echo " BASE             = $BASE"
echo " IMAGE            = $IMAGE"
echo " PY_ENV           = $PY_ENV"
echo " ASCEND_SET_ENV   = $ASCEND_SET_ENV"
echo " TORCH_NPU_SPEC   = ${TORCH_NPU_SPEC:-<from base image>}"
echo " KEEP_BASE_TORCH  = $KEEP_BASE_TORCH"
echo " INSTALL_MX_DRIV  = $INSTALL_MX_DRIVING ($MX_DRIVING_SPEC)"
echo " PROTENIX_GIT_URL = $PROTENIX_GIT_URL"
echo " TARGET           = $TARGET"
echo " CTX              = $CTX"
echo " FILE             = $DOCKERFILE"
echo "=============================================================="

docker build \
  -f "$DOCKERFILE" \
  --target "$TARGET" \
  "${BUILDFLAGS[@]}" \
  -t "$IMAGE" \
  "$CTX"

echo
echo "Build finished: $IMAGE"
echo "Self check: docker run --rm -it --ipc=host --privileged -u=root \\"
echo "              -v /usr/local/Ascend/driver:/usr/local/Ascend/driver $IMAGE --selftest"
