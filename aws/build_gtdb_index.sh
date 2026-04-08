#!/usr/bin/env bash
# Build Dragon Index on AWS
# =========================
# Downloads GTDB genomes and builds a pre-built Dragon index.
#
# Run this AFTER setup_index_builder.sh completes.
#
# Recommended: run inside tmux/screen so it survives SSH disconnects:
#   tmux new -s dragon
#   bash build_gtdb_index.sh
#
# The index will be saved to ~/dragon_index/ and can be compressed
# + uploaded to Zenodo for distribution.

set -euo pipefail

DRAGON="${DRAGON:-$HOME/Dragon/target/release/dragon}"
INDEX_DIR="${1:-$HOME/dragon_index}"
WORKSPACE="$(dirname "$INDEX_DIR")/gtdb_workspace"
GENOME_DIR="$WORKSPACE/genomes"

# System info
echo "====================================="
echo "Dragon GTDB Index Build"
echo "====================================="
echo "CPUs: $(nproc)"
echo "RAM: $(free -h | awk '/Mem:/{print $2}')"
echo "Disk: $(df -h $HOME | awk 'NR==2{print $4}') free"
echo ""

RAM_GB=$(free -g | awk '/Mem:/{print $2}')
# Reserve 8 GB for OS, use rest for Dragon
MAX_RAM=$((RAM_GB - 8))
if [ "$MAX_RAM" -lt 4 ]; then
    MAX_RAM=4
fi
echo "Using --max-ram $MAX_RAM GB for index construction"

# Step 1: Download genomes
echo ""
echo "--- Step 1: Downloading GTDB genomes ---"
echo "This may take 1-2 hours depending on bandwidth..."
mkdir -p "$GENOME_DIR"

$DRAGON download --database gtdb-r220 --output "$GENOME_DIR" 2>&1 | tee "$WORKSPACE/download.log"

GENOME_COUNT=$(find "$GENOME_DIR" -name "*.fa" -o -name "*.fna" -o -name "*.fasta" | wc -l)
echo "Downloaded $GENOME_COUNT genome files"

# Step 2: Build index
echo ""
echo "--- Step 2: Building Dragon index ---"
echo "This may take 2-4 hours..."
START=$(date +%s)

$DRAGON index \
    -i "$GENOME_DIR" \
    -o "$INDEX_DIR" \
    -k 31 \
    --low-memory \
    --max-ram "$MAX_RAM" \
    --threads "$(nproc)" \
    2>&1 | tee "$WORKSPACE/index_build.log"

END=$(date +%s)
ELAPSED=$(( (END - START) / 60 ))
echo "Index build completed in ${ELAPSED} minutes"

# Step 3: Report
echo ""
echo "--- Index Summary ---"
INDEX_SIZE=$(du -sh "$INDEX_DIR" | awk '{print $1}')
echo "Index directory: $INDEX_DIR"
echo "Index size: $INDEX_SIZE"
echo "Genome count: $GENOME_COUNT"
ls -lh "$INDEX_DIR"/

# Step 4: Package for distribution
echo ""
echo "--- Packaging index for distribution ---"
TARBALL="$WORKSPACE/dragon_gtdb_index_$(date +%Y%m%d).tar.gz"
tar -czf "$TARBALL" -C "$(dirname "$INDEX_DIR")" "$(basename "$INDEX_DIR")"
echo "Tarball: $TARBALL ($(du -sh "$TARBALL" | awk '{print $1}'))"
echo ""
echo "====================================="
echo "Done! Upload $TARBALL to Zenodo or S3 for distribution."
echo "Users can download with: dragon download --database gtdb-prebuilt --output ./my_index"
echo "====================================="
