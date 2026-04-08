#!/usr/bin/env bash
# Build Dragon Index from AllTheBacteria
# ========================================
# Downloads ~1.9M bacterial assemblies from AllTheBacteria (EBI),
# builds a Dragon index, and packages it as a distributable tarball.
#
# AllTheBacteria: https://github.com/AllTheBacteria/AllTheBacteria
# Data hosted at: ftp.ebi.ac.uk/pub/databases/AllTheBacteria/
#
# Recommended instance: r6i.4xlarge (16 vCPU, 128 GB RAM)
#   With --low-memory: r6i.2xlarge (8 vCPU, 64 GB RAM)
#
# Disk requirements:
#   Compressed downloads:   ~100 GB
#   Extracted genomes:      ~700 GB
#   Dragon index:           ~50-100 GB
#   Total:                  ~1 TB EBS gp3 recommended
#
# Estimated time: 4-8 hours (download) + 4-8 hours (index build)
#
# Usage:
#   tmux new -s dragon
#   bash build_allthebacteria_index.sh [--low-memory] [--release VERSION]
#
# Run AFTER setup_index_builder.sh completes.

set -euo pipefail

# ---- Configuration ----
ATB_RELEASE="${ATB_RELEASE:-0.2}"
ATB_FTP="https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases"
DRAGON="${DRAGON:-$HOME/Dragon/target/release/dragon}"
WORKSPACE="$HOME/dragon_workspace/allthebacteria"
GENOME_DIR="$WORKSPACE/genomes"
INDEX_DIR="$HOME/dragon_index_allthebacteria"
DOWNLOAD_DIR="$WORKSPACE/downloads"
KMER_SIZE="${KMER_SIZE:-31}"
LOW_MEMORY=false
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --low-memory) LOW_MEMORY=true; shift ;;
        --release) ATB_RELEASE="$2"; shift 2 ;;
        --kmer-size) KMER_SIZE="$2"; shift 2 ;;
        --output) INDEX_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "====================================="
echo "Dragon AllTheBacteria Index Build"
echo "====================================="
echo "AllTheBacteria release: $ATB_RELEASE"
echo "K-mer size: $KMER_SIZE"
echo "Threads: $THREADS"
echo "Low-memory mode: $LOW_MEMORY"
echo "CPUs: $THREADS"
if command -v free &>/dev/null; then
    echo "RAM: $(free -h | awk '/Mem:/{print $2}')"
fi
echo "Disk: $(df -h "$HOME" | awk 'NR==2{print $4}') free"
echo ""

# Compute RAM budget
if command -v free &>/dev/null; then
    RAM_GB=$(free -g | awk '/Mem:/{print $2}')
elif [[ "$(uname)" == "Darwin" ]]; then
    RAM_GB=$(( $(sysctl -n hw.memsize) / 1073741824 ))
else
    RAM_GB=16
fi
MAX_RAM=$((RAM_GB - 8))
if [ "$MAX_RAM" -lt 4 ]; then MAX_RAM=4; fi

# Verify dragon binary
if [ ! -x "$DRAGON" ]; then
    echo "ERROR: Dragon binary not found at $DRAGON"
    echo "Run setup_index_builder.sh first, or set DRAGON= environment variable."
    exit 1
fi

mkdir -p "$WORKSPACE" "$GENOME_DIR" "$DOWNLOAD_DIR" "$INDEX_DIR"

# ===========================================================================
# Step 1: Download AllTheBacteria assemblies
# ===========================================================================
echo ""
echo "--- Step 1: Downloading AllTheBacteria release $ATB_RELEASE ---"
echo "Source: $ATB_FTP/$ATB_RELEASE/assembly/"
echo "This may take several hours depending on bandwidth..."

# AllTheBacteria distributes assemblies as per-species-cluster tarballs.
# The file listing is at the assembly/ directory.
MANIFEST_URL="$ATB_FTP/$ATB_RELEASE/assembly/"
MANIFEST_FILE="$DOWNLOAD_DIR/file_listing.html"

echo "Fetching file listing from $MANIFEST_URL ..."
curl -sL "$MANIFEST_URL" -o "$MANIFEST_FILE"

# Extract tarball filenames (*.tar.xz or *.tar.gz) from the HTML directory listing
TARBALL_LIST="$DOWNLOAD_DIR/tarballs.txt"
grep -oP 'href="\K[^"]+\.(tar\.xz|tar\.gz|tar\.zst)' "$MANIFEST_FILE" \
    | sort -u > "$TARBALL_LIST" 2>/dev/null || true

# Fallback: try simple <a href> parsing if grep -P not available
if [ ! -s "$TARBALL_LIST" ]; then
    grep -o 'href="[^"]*\.tar\.[gxz]\{2,3\}"' "$MANIFEST_FILE" \
        | sed 's/href="//;s/"//' \
        | sort -u > "$TARBALL_LIST" 2>/dev/null || true
fi

# If we still have no tarballs, try the known AllTheBacteria structure:
# Releases/0.2/assembly/ contains files like assembly_*.tar.xz
if [ ! -s "$TARBALL_LIST" ]; then
    echo "Could not parse directory listing. Trying alternative download strategy..."
    # AllTheBacteria v0.2 uses a single large archive or per-chunk archives
    # Try the known patterns
    for CHUNK_URL in \
        "$ATB_FTP/$ATB_RELEASE/assembly/allthebacteria_assemblies.tar.xz" \
        "$ATB_FTP/$ATB_RELEASE/assembly/all_assemblies.tar.xz" \
        "$ATB_FTP/$ATB_RELEASE/assembly/"; do
        echo "  Trying: $CHUNK_URL"
        if curl -sIf "$CHUNK_URL" > /dev/null 2>&1; then
            echo "$CHUNK_URL" > "$TARBALL_LIST"
            break
        fi
    done
fi

TOTAL_FILES=$(wc -l < "$TARBALL_LIST" | tr -d ' ')
echo "Found $TOTAL_FILES archive(s) to download"

if [ "$TOTAL_FILES" -eq 0 ]; then
    echo "ERROR: No assembly archives found at $MANIFEST_URL"
    echo "Check the AllTheBacteria release page and update ATB_RELEASE."
    echo "You can also manually download genomes into $GENOME_DIR and re-run with --skip-download"
    exit 1
fi

# Download and extract each archive
DOWNLOAD_COUNT=0
for TARBALL_NAME in $(cat "$TARBALL_LIST"); do
    DOWNLOAD_COUNT=$((DOWNLOAD_COUNT + 1))

    # Handle relative vs absolute URLs
    if [[ "$TARBALL_NAME" == http* ]]; then
        TARBALL_URL="$TARBALL_NAME"
    else
        TARBALL_URL="$ATB_FTP/$ATB_RELEASE/assembly/$TARBALL_NAME"
    fi

    LOCAL_FILE="$DOWNLOAD_DIR/$(basename "$TARBALL_NAME")"

    echo "  [$DOWNLOAD_COUNT/$TOTAL_FILES] Downloading $(basename "$TARBALL_NAME")..."

    # Skip if already downloaded
    if [ -f "$LOCAL_FILE" ]; then
        echo "    (already exists, skipping download)"
    else
        curl -L --retry 3 --retry-delay 5 -o "$LOCAL_FILE" "$TARBALL_URL" --progress-bar
    fi

    # Extract FASTA files into GENOME_DIR
    echo "    Extracting..."
    case "$LOCAL_FILE" in
        *.tar.xz)  tar -xJf "$LOCAL_FILE" -C "$GENOME_DIR" --wildcards '*.fna' '*.fa' '*.fasta' '*.fna.gz' '*.fa.gz' 2>/dev/null || tar -xJf "$LOCAL_FILE" -C "$GENOME_DIR" ;;
        *.tar.gz)   tar -xzf "$LOCAL_FILE" -C "$GENOME_DIR" --wildcards '*.fna' '*.fa' '*.fasta' '*.fna.gz' '*.fa.gz' 2>/dev/null || tar -xzf "$LOCAL_FILE" -C "$GENOME_DIR" ;;
        *.tar.zst)  zstd -d "$LOCAL_FILE" --stdout | tar -xf - -C "$GENOME_DIR" --wildcards '*.fna' '*.fa' '*.fasta' 2>/dev/null || { zstd -d "$LOCAL_FILE" -o "${LOCAL_FILE%.zst}" && tar -xf "${LOCAL_FILE%.zst}" -C "$GENOME_DIR"; } ;;
    esac

    # Remove archive after extraction to save disk
    rm -f "$LOCAL_FILE"
done

# Decompress any gzipped FASTAs
echo "Decompressing any .gz FASTA files..."
find "$GENOME_DIR" -name "*.fna.gz" -exec gunzip -f {} \; 2>/dev/null || true
find "$GENOME_DIR" -name "*.fa.gz" -exec gunzip -f {} \; 2>/dev/null || true
find "$GENOME_DIR" -name "*.fasta.gz" -exec gunzip -f {} \; 2>/dev/null || true

# Flatten nested directories: move all FASTA files to GENOME_DIR root
find "$GENOME_DIR" -mindepth 2 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) \
    -exec mv -n {} "$GENOME_DIR/" \; 2>/dev/null || true
# Clean empty subdirectories
find "$GENOME_DIR" -mindepth 1 -type d -empty -delete 2>/dev/null || true

GENOME_COUNT=$(find "$GENOME_DIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | wc -l | tr -d ' ')
echo ""
echo "Downloaded and extracted $GENOME_COUNT genome files"
echo "Genome directory size: $(du -sh "$GENOME_DIR" | awk '{print $1}')"

if [ "$GENOME_COUNT" -eq 0 ]; then
    echo "ERROR: No genome FASTA files found after extraction."
    echo "Check download logs in $DOWNLOAD_DIR/"
    exit 1
fi

# ===========================================================================
# Step 2: Build Dragon index
# ===========================================================================
echo ""
echo "--- Step 2: Building Dragon index ---"
echo "Genomes: $GENOME_COUNT"
echo "K-mer size: $KMER_SIZE"
echo "Output: $INDEX_DIR"

BUILD_START=$(date +%s)

INDEX_CMD="$DRAGON index -i $GENOME_DIR -o $INDEX_DIR -k $KMER_SIZE --threads $THREADS"
if [ "$LOW_MEMORY" = true ]; then
    INDEX_CMD="$INDEX_CMD --low-memory --max-ram $MAX_RAM"
    echo "Low-memory mode: --max-ram $MAX_RAM GB"
fi

echo "Running: $INDEX_CMD"
echo "This may take several hours..."
eval "$INDEX_CMD" 2>&1 | tee "$WORKSPACE/index_build.log"

BUILD_END=$(date +%s)
BUILD_ELAPSED=$(( (BUILD_END - BUILD_START) / 60 ))
echo ""
echo "Index build completed in ${BUILD_ELAPSED} minutes"

# ===========================================================================
# Step 3: Summary
# ===========================================================================
echo ""
echo "--- Index Summary ---"
INDEX_SIZE=$(du -sh "$INDEX_DIR" | awk '{print $1}')
echo "Index directory: $INDEX_DIR"
echo "Index size: $INDEX_SIZE"
echo "Genome count: $GENOME_COUNT"
echo "K-mer size: $KMER_SIZE"
ls -lh "$INDEX_DIR"/

# ===========================================================================
# Step 4: Package for distribution
# ===========================================================================
echo ""
echo "--- Packaging index for distribution ---"
TARBALL="$WORKSPACE/dragon_allthebacteria_v${ATB_RELEASE}_k${KMER_SIZE}_$(date +%Y%m%d).tar.gz"
tar -czf "$TARBALL" -C "$(dirname "$INDEX_DIR")" "$(basename "$INDEX_DIR")"
TARBALL_SIZE=$(du -sh "$TARBALL" | awk '{print $1}')
echo "Tarball: $TARBALL ($TARBALL_SIZE)"

# Write metadata
cat > "$WORKSPACE/index_metadata.json" <<METAEOF
{
    "database": "AllTheBacteria",
    "release": "$ATB_RELEASE",
    "genome_count": $GENOME_COUNT,
    "kmer_size": $KMER_SIZE,
    "index_size": "$INDEX_SIZE",
    "tarball_size": "$TARBALL_SIZE",
    "build_time_minutes": $BUILD_ELAPSED,
    "build_date": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "dragon_version": "$($DRAGON --version 2>/dev/null || echo unknown)",
    "low_memory": $LOW_MEMORY
}
METAEOF
echo "Metadata: $WORKSPACE/index_metadata.json"

echo ""
echo "====================================="
echo "Done! AllTheBacteria index built successfully."
echo ""
echo "Upload to Zenodo/S3:"
echo "  $TARBALL"
echo ""
echo "Users can download with:"
echo "  dragon download --database allthebacteria-v2 --output ./my_index"
echo ""
echo "Or search directly:"
echo "  dragon search --index $INDEX_DIR --query <reads.fasta>"
echo "====================================="
