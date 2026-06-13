#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 0 ]]; then
  echo "Usage: OUT_ROOT=/path/to/output bash scripts/nanobody/datasets/download_sabdab_nano_temporal.sh" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if command -v python3 >/dev/null 2>&1; then
  PYTHON="python3"
elif command -v python >/dev/null 2>&1; then
  PYTHON="python"
else
  echo "ERROR: neither python3 nor python is available on PATH." >&2
  exit 1
fi

OUT_ROOT="${OUT_ROOT:-${PROTENIX_ROOT_DIR:-$(pwd)}/nanobody}"
mkdir -p "$OUT_ROOT"

"$PYTHON" "$SCRIPT_DIR/download_sabdab_nano_temporal.py" --output-root "$OUT_ROOT"
"$PYTHON" "$SCRIPT_DIR/build_sabdab_nano_temporal_splits.py" --output-root "$OUT_ROOT"
"$PYTHON" "$SCRIPT_DIR/verify_sabdab_nano_temporal.py" --output-root "$OUT_ROOT"

echo "Split manifests:"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/single_nano_train_pre2022_manifest.csv"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/single_nano_val_2022h1_manifest.csv"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/single_nano_test_2022h2_manifest.csv"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/complex_nanoag_train_pre2022_manifest.csv"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/complex_nanoag_val_2022h1_manifest.csv"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/complex_nanoag_test_2022h2_manifest.csv"
echo "Verification report:"
echo "$OUT_ROOT/data/sabdab_nano_temporal/manifests/verify_report.json"
