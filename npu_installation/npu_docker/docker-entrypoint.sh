#!/bin/bash
# Container entrypoint: sources the full environment so that users do not have to
# remember any environment variable.
set -e

# 1) Unified environment (CANN + custom ops + conda), provided by the image
if [ -f /etc/profile.d/ptx-npu.sh ]; then
    source /etc/profile.d/ptx-npu.sh
fi

# 2) Self check sub command: docker run ... --selftest [--] [command]
if [ "$1" = "--selftest" ] || [ "$1" = "selftest" ]; then
    shift
    [ "${1:-}" = "--" ] && shift
    check_npu_env || echo "[WARN] self check failed, dropping into a shell for debugging"
    if [ $# -eq 0 ]; then
        exec bash -l
    fi
fi

# 3) Otherwise: run the user command from the Protenix directory
cd /work/Protenix 2>/dev/null || true
exec "$@"
