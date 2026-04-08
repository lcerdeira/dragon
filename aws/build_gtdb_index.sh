#!/usr/bin/env bash
# Build Dragon Index from GTDB r220
# ===================================
# Downloads GTDB representative genomes and builds a Dragon index.
#
# Disk requirements:
#   - Compressed download: ~85 GB
#   - Extracted genomes: ~200 GB
#   - Dragon index: ~50 GB
#   - Total: ~400 GB recommended
#
# Run inside tmux/screen so it survives SSH disconnects:
#   screen -S gtdb
#   bash build_gtdb_index.sh --output /data/gtdb_index

set -euo pipefail

# ---- Configuration ----
GTDB_VERSION="220"
GTDB_URL="https://data.gtdb.ecogenomic.org/releases/release${GTDB_VERSION}/${GTDB_VERSION}.0/genomic_files_reps/gtdb_genomes_reps_r${GTDB_VERSION}.tar.gz"
DRAGON="${DRAGON:-$HOME/Dragon/target/release/dragon}"
INDEX_DIR=""
WORKSPACE=""
GENOME_DIR=""
DOWNLOAD_DIR=""
KMER_SIZE="${KMER_SIZE:-31}"
LOW_MEMORY=false
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output) INDEX_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --kmer-size) KMER_SIZE="$2"; shift 2 ;;
        --low-memory) LOW_MEMORY=true; shift ;;
        --version) GTDB_VERSION="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Default paths — derive workspace from output directory
INDEX_DIR="${INDEX_DIR:-$HOME/dragon_index_gtdb}"
WORKSPACE="$(dirname "$INDEX_DIR")/gtdb_workspace"
GENOME_DIR="$WORKSPACE/genomes"
DOWNLOAD_DIR="$WORKSPACE/downloads"

echo "====================================="
echo "Dragon GTDB r${GTDB_VERSION} Index Build"
echo "====================================="
echo "GTDB version: r${GTDB_VERSION}"
echo "K-mer size: $KMER_SIZE"
echo "Threads: $THREADS"
echo "Low-memory mode: $LOW_MEMORY"
echo "Output: $INDEX_DIR"
echo "Workspace: $WORKSPACE"
if command -v nproc &>/dev/null; then
    echo "CPUs: $(nproc)"
fi
if command -v free &>/dev/null; then
    echo "RAM: $(free -h | awk '/Mem:/{print $2}')"
fi
echo "Disk: $(df -h "$(dirname "$INDEX_DIR")" | awk 'NR==2{print $4}') free"
echo ""

# Check Dragon binary
if [ ! -x "$DRAGON" ]; then
    echo "ERROR: Dragon binary not found at $DRAGON"
    echo "Set DRAGON= environment variable to the path of the dragon binary."
    exit 1
fi

# RAM calculation for --low-memory
if command -v free &>/dev/null; then
    RAM_GB=$(free -g | awk '/Mem:/{print $2}')
    MAX_RAM=$((RAM_GB - 8))
    if [ "$MAX_RAM" -lt 4 ]; then MAX_RAM=4; fi
else
    MAX_RAM=8
fi

# ===========================================================================
# Step 1: Download GTDB representative genomes
# ===========================================================================
echo "--- Step 1: Downloading GTDB r${GTDB_VERSION} representative genomes ---"
echo "Source: $GTDB_URL"
echo "This may take 1-2 hours depending on bandwidth..."
mkdir -p "$DOWNLOAD_DIR" "$GENOME_DIR"

TARBALL="$DOWNLOAD_DIR/gtdb_genomes_reps_r${GTDB_VERSION}.tar.gz"

if [ -f "$TARBALL" ]; then
    echo "  Tarball already exists, skipping download"
else
    if command -v curl &>/dev/null; then
        curl -L --retry 3 --retry-delay 10 -o "$TARBALL" "$GTDB_URL"
    elif command -v wget &>/dev/null; then
        wget --retry-connrefused --tries=3 -O "$TARBALL" "$GTDB_URL"
    else
        echo "ERROR: Neither curl nor wget found"
        exit 1
    fi
fi

echo "Download complete: $(du -sh "$TARBALL" | awk '{print $1}')"

# ===========================================================================
# Step 2: Extract genomes
# ===========================================================================
echo ""
echo "--- Step 2: Extracting genomes ---"
echo "This may take 30-60 minutes..."

tar -xzf "$TARBALL" -C "$GENOME_DIR"

# GTDB tarball has nested structure: GCA/GCF dirs with .fna.gz files
# Decompress and flatten
echo "Decompressing genome files..."
find "$GENOME_DIR" -name "*.fna.gz" -exec gunzip -f {} \; 2>/dev/null || true
find "$GENOME_DIR" -name "*.fa.gz" -exec gunzip -f {} \; 2>/dev/null || true
find "$GENOME_DIR" -name "*.fasta.gz" -exec gunzip -f {} \; 2>/dev/null || true

GENOME_COUNT=$(find "$GENOME_DIR" -name "*.fna" -o -name "*.fa" -o -name "*.fasta" | wc -l | tr -d ' ')
echo "Extracted $GENOME_COUNT genome files"
echo "Genome directory size: $(du -sh "$GENOME_DIR" | awk '{print $1}')"

if [ "$GENOME_COUNT" -eq 0 ]; then
    echo "ERROR: No genome files found after extraction."
    echo "Check the tarball structure:"
    tar -tzf "$TARBALL" | head -20
    exit 1
fi

# Delete tarball to save disk space
echo "Removing tarball to free disk space..."
rm -f "$TARBALL"

# ===========================================================================
# Step 3: Build Dragon index
# ===========================================================================
echo ""
echo "--- Step 3: Building Dragon index ---"
echo "This may take 2-6 hours depending on genome count and hardware..."
START=$(date +%s)

INDEX_CMD="$DRAGON index -i $GENOME_DIR -o $INDEX_DIR -k $KMER_SIZE --threads $THREADS"
if [ "$LOW_MEMORY" = true ]; then
    INDEX_CMD="$INDEX_CMD --low-memory --max-ram $MAX_RAM"
fi

echo "Running: $INDEX_CMD"
eval "$INDEX_CMD" 2>&1 | tee "$WORKSPACE/index_build.log"

END=$(date +%s)
ELAPSED=$(( (END - START) / 60 ))
echo "Index build completed in ${ELAPSED} minutes"

# ===========================================================================
# Step 4: Report and package
# ===========================================================================
echo ""
echo "--- Index Summary ---"
INDEX_SIZE=$(du -sh "$INDEX_DIR" | awk '{print $1}')
echo "Index directory: $INDEX_DIR"
echo "Index size: $INDEX_SIZE"
echo "Genome count: $GENOME_COUNT"
ls -lh "$INDEX_DIR"/

echo ""
echo "--- Packaging index for distribution ---"
OUTPUT_TARBALL="$(dirname "$INDEX_DIR")/dragon_gtdb_r${GTDB_VERSION}_k${KMER_SIZE}_$(date +%Y%m%d).tar.gz"
tar -czf "$OUTPUT_TARBALL" -C "$(dirname "$INDEX_DIR")" "$(basename "$INDEX_DIR")"
echo "Tarball: $OUTPUT_TARBALL ($(du -sh "$OUTPUT_TARBALL" | awk '{print $1}'))"

# Write metadata
cat > "$INDEX_DIR/metadata.json" <<METAEOF
{
    "database": "gtdb",
    "version": "r${GTDB_VERSION}",
    "genome_count": $GENOME_COUNT,
    "kmer_size": $KMER_SIZE,
    "build_date": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "build_time_minutes": $ELAPSED,
    "index_size": "$INDEX_SIZE",
    "dragon_version": "$($DRAGON --version 2>/dev/null || echo unknown)",
    "source": "$GTDB_URL"
}
METAEOF

echo ""
echo "====================================="
echo "Done! Upload to S3 or Zenodo for distribution:"
echo "  aws s3 cp $OUTPUT_TARBALL s3://your-bucket/dragon-indices/"
echo ""
echo "Users can download with:"
echo "  dragon download --database gtdb-r220 --output ./my_index"
echo "====================================="
