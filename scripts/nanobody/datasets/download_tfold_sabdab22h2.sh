#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 0 ]]; then
  echo "Usage: OUT_ROOT=/path/to/output bash scripts/nanobody/datasets/download_tfold_sabdab22h2.sh" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="${OUT_ROOT:-${PROTENIX_ROOT_DIR:-$(pwd)}/nanobody}"
mkdir -p "$OUT_ROOT"

if command -v python3 >/dev/null 2>&1; then
  PYTHON="python3"
elif command -v python >/dev/null 2>&1; then
  PYTHON="python"
else
  echo "ERROR: neither python3 nor python is available on PATH." >&2
  exit 1
fi

DATASET_DIR="$OUT_ROOT/data/tfold_sabdab22h2"
ARCHIVE="$DATASET_DIR/raw/tFold_test_set.tar.gz"
NANO_SUBSET="$DATASET_DIR/subsets/SAbDab-22H2-Nano"
NANOAG_SUBSET="$DATASET_DIR/subsets/SAbDab-22H2-NanoAg"
NANO_MANIFEST="$DATASET_DIR/manifests/SAbDab-22H2-Nano_manifest.csv"
NANOAG_MANIFEST="$DATASET_DIR/manifests/SAbDab-22H2-NanoAg_manifest.csv"
VERIFY_REPORT="$DATASET_DIR/manifests/verify_report.json"

echo "Output root: $OUT_ROOT"
echo "Official source URLs:"
echo "  https://github.com/TencentAI4S/tfold"
echo "  https://drive.google.com/file/d/1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ/view?usp=drive_link"
echo "  https://drive.google.com/uc?export=download&id=1szSr5bjP3Y6XbhUpbfZEb9ZL9UMPXtvZ"
echo "  https://share.weiyun.com/zycZDrfA"
echo "Archive path: $ARCHIVE"

"$PYTHON" "$SCRIPT_DIR/download_tfold_sabdab22h2.py" --output-root "$OUT_ROOT"
"$PYTHON" "$SCRIPT_DIR/verify_tfold_sabdab22h2.py" --output-root "$OUT_ROOT"

echo "Dataset directory: $DATASET_DIR"
echo "Subset: $NANO_SUBSET"
echo "Subset: $NANOAG_SUBSET"
echo "Manifest: $NANO_MANIFEST"
echo "Manifest: $NANOAG_MANIFEST"
echo "Verification report: $VERIFY_REPORT"
