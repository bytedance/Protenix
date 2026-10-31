#!/bin/bash
# One-shot environment self check: CANN / op vendor / torch_npu / dpse patch / code
# Usage: docker exec <container> check_npu_env      (does not rely on the caller shell)
#        or: docker run ... <image> --selftest
set -u

# ---- Prepare the environment here: docker exec does not run the ENTRYPOINT ----
# The vendor set_env.bash uses ${ASCEND_CUSTOM_OPP_PATH}:... syntax, which triggers
# "unbound variable" under set -u, so seed the variables and relax -u while sourcing.
export ASCEND_CUSTOM_OPP_PATH="${ASCEND_CUSTOM_OPP_PATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
[ -f /usr/local/Ascend/ascend-toolkit/set_env.sh ] && source /usr/local/Ascend/ascend-toolkit/set_env.sh
set +u
for f in /opt/ops/*/vendors/*/bin/set_env.bash; do [ -f "$f" ] && source "$f"; done
set -u

# Use an absolute python path: the base conda env has no torch
PYBIN="${PY:-/opt/conda/envs/${PROTENIX_CONDA_ENV:-torch2.7.1}/bin/python}"
[ -x "$PYBIN" ] || PYBIN="$(command -v python3 || echo /usr/bin/python3)"

FAIL=0
ok()  { echo "  [ OK ] $*"; }
bad() { echo "  [FAIL] $*"; FAIL=1; }

echo "=============================================================="
echo " Protenix NPU environment self check"
echo "=============================================================="

echo "== 1. CANN =="
if [ -n "${ASCEND_HOME_PATH:-}" ] && [ -d "${ASCEND_HOME_PATH}" ]; then
    ok "ASCEND_HOME_PATH=${ASCEND_HOME_PATH}"
else
    bad "ASCEND_HOME_PATH is unset or does not exist"
fi
echo "  ASCEND_OPP_PATH=${ASCEND_OPP_PATH:-<unset>}"
echo "  ASCEND_CUSTOM_OPP_PATH=${ASCEND_CUSTOM_OPP_PATH:-<unset>}"
echo "  python=${PYBIN}"

echo "== 2. Custom op vendor =="
ls -d /opt/ops/*/vendors/* 2>/dev/null | sed 's/^/  /'
if compgen -G "/opt/ops/*/vendors/*" > /dev/null; then
    ok "at least one vendor is installed"
else
    bad "no /opt/ops/*/vendors/* found"
fi

echo "== 3. torch / torch_npu / NPU visibility =="
"$PYBIN" - <<'PY'
import sys
try:
    import torch, torch_npu
    print("  torch      :", torch.__version__)
    print("  torch_npu  :", torch_npu.__version__)
    n = torch.npu.device_count()
    print("  device_cnt :", n)
    if torch.npu.is_available() and n > 0:
        print("  [ OK ] NPU is_available")
    else:
        print("  [FAIL] NPU is not visible (check --device=/dev/davinci*, the driver mount and liburma.so)")
        sys.exit(1)
except Exception as e:
    print("  [FAIL] importing torch/torch_npu failed:", e)
    sys.exit(1)
PY
[ $? -ne 0 ] && FAIL=1

echo "== 4. dpse patch (1NSS broadcast pse, B=2) =="
"$PYBIN" - <<'PY'
import sys, torch, torch_npu
try:
    torch.npu.set_device(0)
    torch.manual_seed(0)
    B, N, S, D = 2, 4, 256, 32
    q = torch.randn(B, N, S, D, dtype=torch.bfloat16, device="npu:0").requires_grad_()
    k = torch.randn_like(q).requires_grad_()
    v = torch.randn_like(q).requires_grad_()
    pse = torch.zeros(1, N, S, S, dtype=torch.bfloat16, device="npu:0").requires_grad_()
    o = torch_npu.npu_fusion_attention(q, k, v, N, "BNSD", pse=pse)[0]
    o.sum().backward()
    nz = int((pse.grad != 0).sum()); total = pse.grad.numel()
    print(f"  dpse nonzero = {nz} / {total}")
    print("  [ OK ] dpse gets gradients (patch is effective)" if nz > 0 else "  [FAIL] dpse is all zeros (patch is not effective)")
    sys.exit(0 if nz > 0 else 1)
except Exception as e:
    print("  [FAIL] dpse self check raised:", e); sys.exit(1)
PY
[ $? -ne 0 ] && FAIL=1

echo "== 5. Protenix code =="
if [ -f /work/Protenix/runner/train.py ] && [ -d /work/Protenix/migrate_to_npu ]; then
    ok "Protenix code and migrate_to_npu are in place"
else
    bad "Protenix code is missing (expected at /work/Protenix)"
fi

echo "=============================================================="
[ "$FAIL" -eq 0 ] && echo " Result: PASS" || echo " Result: FAIL (see the [FAIL] lines above)"
echo "=============================================================="
exit "$FAIL"
