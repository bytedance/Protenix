#!/bin/bash

# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# --- Script Configuration & Safety Settings ---
set -euo pipefail

# Log functions for formatted output
info() { echo -e "\033[1;34m[INFO]\033[0m $1"; }
err() { echo -e "\033[1;31m[ERROR]\033[0m $1"; exit 1; }

# --- Help Message ---
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --inference_only    Download only the data needed for inference (common, mmcif, search_database) (default)"
    echo "  --full              Download all data components for training and finetuning"
    echo "  -h, --help          Show this help message"
    echo ""
    exit 0
}

# --- Parse Arguments ---
DOWNLOAD_MODE="inference"
while [[ $# -gt 0 ]]; do
    case $1 in
        --inference_only)
            DOWNLOAD_MODE="inference"
            shift
            ;;
        --full)
            DOWNLOAD_MODE="full"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            err "Unknown option: $1"
            ;;
    esac
done

# --- Check Environment Variable ---
if [ -z "${PROTENIX_ROOT_DIR:-}" ]; then
    err "PROTENIX_ROOT_DIR is not set. Please set it before running this script.\nExample: export PROTENIX_ROOT_DIR=/path/to/your/data_root"
fi

info "Using PROTENIX_ROOT_DIR: $PROTENIX_ROOT_DIR"
info "Download mode: $DOWNLOAD_MODE"
mkdir -p "$PROTENIX_ROOT_DIR"
# --- Download & Extract Data Components ---
if [ "$DOWNLOAD_MODE" == "inference" ]; then
    DATA_FILES=(
        "common.tar.gz"
#        "mmcif.tar.gz"
        "search_database.tar.gz"
    )
else
    DATA_FILES=(
        "rna_msa.tar.gz"
        "common.tar.gz"
        "indices.tar.gz"
        "mmcif.tar.gz"
        "mmcif_bioassembly.tar.gz"
        "mmcif_msa_template.tar"
        "posebusters_bioassembly.tar.gz"
        "posebusters_mmcif.tar.gz"
        "recentPDB_bioassembly.tar.gz"
        "search_database.tar.gz"
    )
fi

BASE_URL="https://protenix.tos-cn-beijing.volces.com"

for file in "${DATA_FILES[@]}"; do
    info "Downloading $file from $BASE_URL/$file ..."
    wget -c -P "$PROTENIX_ROOT_DIR" "$BASE_URL/$file"

    info "Extracting $file to $PROTENIX_ROOT_DIR ..."
    if [[ "$file" == *.tar.gz ]]; then
        tar -xzvf "$PROTENIX_ROOT_DIR/$file" -C "$PROTENIX_ROOT_DIR"
    else
        tar -xvf "$PROTENIX_ROOT_DIR/$file" -C "$PROTENIX_ROOT_DIR"
    fi

    info "Removing $file ..."
    rm "$PROTENIX_ROOT_DIR/$file"
done

info "Selected data components ($DOWNLOAD_MODE) have been downloaded and extracted successfully to $PROTENIX_ROOT_DIR"
