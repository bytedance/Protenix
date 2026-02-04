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
warn() { echo -e "\033[1;33m[WARN]\033[0m $1"; }
err() { echo -e "\033[1;31m[ERROR]\033[0m $1"; exit 1; }

# --- Environment Setup ---
LOCAL_BIN="$HOME/.local/bin"
mkdir -p "$LOCAL_BIN"
export PATH="$LOCAL_BIN:$PATH"

# --- 1. Check & Install Base Tools ---
info "Checking for base dependencies..."
for tool in curl wget zcat; do
    if ! command -v "$tool" &> /dev/null; then
        warn "Tool '$tool' not found in PATH."
        if [ "$tool" == "zcat" ]; then
            info "Checking for 'gzip' as a fallback..."
            command -v gzip &> /dev/null || err "gzip is missing."
        else
            err "Required tool '$tool' is missing."
        fi
    fi
done

# --- 2. Local Installation of MMseqs2 ---
if ! command -v mmseqs &> /dev/null; then
    info "MMseqs2 not found. Installing to $LOCAL_BIN ..."
    wget -q https://mmseqs.com/latest/mmseqs-linux64.tar.gz
    tar -xzf mmseqs-linux64.tar.gz
    cp mmseqs/bin/mmseqs "$LOCAL_BIN/"
    rm -rf mmseqs mmseqs-linux64.tar.gz
    info "MMseqs2 installed successfully."
else
    info "MMseqs2 is already installed. Skipping installation."
fi

# --- 3. Download & Extract RNAcentral Data ---
info "Downloading RNAcentral active sequences..."
# Using pget for faster download if lftp is available
if command -v lftp &> /dev/null; then
    lftp -c "pget -n 10 https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz"
else
    wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
fi

info "Extracting RNAcentral data..."
if command -v pigz &> /dev/null; then
    pigz -d rnacentral_active.fasta.gz
else
    gunzip rnacentral_active.fasta.gz
fi

# --- 4. MMseqs2 Linclust Pipeline (Multiple Thresholds) ---
info "Starting MMseqs2 linclust with multiple thresholds..."

# Step: Define parameter sets: ID COV SUFFIX
declare -a configs=(
    "0.9 0.8 strict"
    "0.8 0.7 medium"
    "0.7 0.6 loose"
)

input_fasta="rnacentral_active.fasta"

for config in "${configs[@]}"; do
    read min_id cov suffix <<< "$config"
    info "Running easy-linclust: identity=${min_id}, coverage=${cov} (${suffix})"

    # Output prefix and final filename
    out_prefix="rnacentral_clust_${suffix}"
    tmp_dir="tmp_linclust_${suffix}"

    # Convert identity/coverage to integer-style for filename (e.g., 0.9 -> 90)
    id_int=$(awk "BEGIN {printf \"%.0f\", ${min_id}*100}")
    cov_int=$(awk "BEGIN {printf \"%.0f\", ${cov}*100}")
    final_fasta="rnacentral_active_seq_id_${id_int}_cov_${cov_int}_linclust.fasta"

    # Run easy-linclust
    mmseqs easy-linclust "$input_fasta" "$out_prefix" "$tmp_dir" \
        --min-seq-id "$min_id" \
        -c "$cov" \
        --cov-mode 1

    # Rename representative sequences to desired format
    mv "${out_prefix}_rep_seq.fasta" "$final_fasta"

    # Optional: Print stats
    rep_count=$(grep -c "^>" "$final_fasta")
    echo "  -> Representative sequences (${suffix}): $rep_count"

    # Clean up intermediate files
    rm -rf "${out_prefix}_cluster.tsv" "${out_prefix}_all_seqs.fasta" "$tmp_dir"
done

echo "----------------------------------------------------"
echo "All easy-linclust runs completed!"
for f in rnacentral_active_seq_id_*_linclust.fasta; do
    if [ -f "$f" ]; then
        count=$(grep -c "^>" "$f")
        echo "File: $f -> $count representatives"
    fi
done
echo "----------------------------------------------------"
