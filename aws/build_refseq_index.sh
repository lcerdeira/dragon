#!/usr/bin/env bash
# Build Dragon Index from NCBI RefSeq Bacteria
# ==============================================
# Downloads complete/representative bacterial genomes from NCBI RefSeq,
# builds a Dragon index, and packages it as a distributable tarball.
#
# RefSeq: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
# Assembly summary: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
#
# Recommended instance: r6i.4xlarge (16 vCPU, 128 GB RAM)
#   With --low-memory: r6i.2xlarge (8 vCPU, 64 GB RAM)
#
# Disk requirements:
#   Compressed downloads:   ~200 GB
#   Extracted genomes:      ~800 GB  (depends on --representative-only)
#   Dragon index:           ~50-100 GB
#   Total:                  ~1.2 TB EBS gp3 recommended
#                           ~300 GB with --representative-only
#
# Estimated time:
#   Representative only (~100K genomes): 2-4h download + 2-4h index
#   All complete (~300K genomes):        6-12h download + 6-12h index
#
# Usage:
#   tmux new -s dragon
#   bash build_refseq_index.sh [--representative-only] [--low-memory]
#
# Run AFTER setup_index_builder.sh completes.

set -euo pipefail

# ---- Configuration ----
REFSEQ_FTP="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria"
ASSEMBLY_SUMMARY_URL="$REFSEQ_FTP/assembly_summary.txt"
DRAGON="${DRAGON:-$HOME/Dragon/target/release/dragon}"
WORKSPACE="$HOME/dragon_workspace/refseq"
GENOME_DIR="$WORKSPACE/genomes"
INDEX_DIR="$HOME/dragon_index_refseq"
DOWNLOAD_DIR="$WORKSPACE/downloads"
KMER_SIZE="${KMER_SIZE:-31}"
LOW_MEMORY=false
REPRESENTATIVE_ONLY=false
PARALLEL_DOWNLOADS=8
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --low-memory)           LOW_MEMORY=true; shift ;;
        --representative-only)  REPRESENTATIVE_ONLY=true; shift ;;
        --kmer-size)            KMER_SIZE="$2"; shift 2 ;;
        --output)               INDEX_DIR="$2"; shift 2 ;;
        --threads)              THREADS="$2"; shift 2 ;;
        --parallel-downloads)   PARALLEL_DOWNLOADS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "====================================="
echo "Dragon RefSeq Bacteria Index Build"
echo "====================================="
echo "K-mer size: $KMER_SIZE"
echo "Threads: $THREADS"
echo "Low-memory mode: $LOW_MEMORY"
echo "Representative only: $REPRESENTATIVE_ONLY"
echo "Parallel downloads: $PARALLEL_DOWNLOADS"
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

mkdir -p "$WORKSPACE" "$GENOME_DIR" "$DOWNLOAD_DIR"

# ===========================================================================
# Step 1: Download and parse assembly summary
# ===========================================================================
echo ""
echo "--- Step 1: Downloading RefSeq assembly summary ---"

SUMMARY_FILE="$DOWNLOAD_DIR/assembly_summary.txt"
echo "Fetching $ASSEMBLY_SUMMARY_URL ..."
curl -L --retry 3 -o "$SUMMARY_FILE" "$ASSEMBLY_SUMMARY_URL" --progress-bar

# Parse assembly summary to get FTP paths
# Columns (tab-separated, # header lines):
#   1:  assembly_accession
#   5:  refseq_category (representative genome, reference genome, na)
#   12: assembly_level (Complete Genome, Chromosome, Scaffold, Contig)
#   20: ftp_path
URL_LIST="$DOWNLOAD_DIR/ftp_urls.txt"

echo "Parsing assembly summary..."
if [ "$REPRESENTATIVE_ONLY" = true ]; then
    # Only representative and reference genomes
    awk -F'\t' '
        !/^#/ && ($5 == "representative genome" || $5 == "reference genome") && $20 != "na" {
            # Build the URL: ftp_path + "/" + basename + "_genomic.fna.gz"
            n = split($20, parts, "/")
            basename = parts[n]
            print $20 "/" basename "_genomic.fna.gz"
        }
    ' "$SUMMARY_FILE" > "$URL_LIST"
    echo "Filter: representative + reference genomes only"
else
    # All complete + chromosome-level genomes
    awk -F'\t' '
        !/^#/ && ($12 == "Complete Genome" || $12 == "Chromosome") && $20 != "na" {
            n = split($20, parts, "/")
            basename = parts[n]
            print $20 "/" basename "_genomic.fna.gz"
        }
    ' "$SUMMARY_FILE" > "$URL_LIST"
    echo "Filter: Complete Genome + Chromosome-level assemblies"
fi

# Convert FTP to HTTPS (NCBI supports both)
sed -i 's|^ftp://|https://|' "$URL_LIST"

TOTAL_GENOMES=$(wc -l < "$URL_LIST" | tr -d ' ')
echo "Found $TOTAL_GENOMES genomes to download"

if [ "$TOTAL_GENOMES" -eq 0 ]; then
    echo "ERROR: No genomes found in assembly summary."
    echo "Check $SUMMARY_FILE for issues."
    exit 1
fi

# ===========================================================================
# Step 2: Download genome FASTA files
# ===========================================================================
echo ""
echo "--- Step 2: Downloading $TOTAL_GENOMES genome files ---"
echo "Using $PARALLEL_DOWNLOADS parallel downloads..."
echo "This may take several hours..."

DOWNLOAD_START=$(date +%s)

# Download function for a single genome
download_genome() {
    local URL="$1"
    local BASENAME
    BASENAME=$(basename "$URL")
    local OUTFILE="$GENOME_DIR/$BASENAME"
    local FASTA="${OUTFILE%.gz}"

    # Skip if already extracted
    if [ -f "$FASTA" ]; then
        return 0
    fi

    # Skip if already downloaded (compressed)
    if [ -f "$OUTFILE" ]; then
        gunzip -f "$OUTFILE" 2>/dev/null || true
        return 0
    fi

    # Download + decompress in one step
    curl -sL --retry 3 --retry-delay 2 "$URL" 2>/dev/null | gunzip -c > "$FASTA" 2>/dev/null
    if [ ! -s "$FASTA" ]; then
        # Retry with file-based download
        curl -sL --retry 3 -o "$OUTFILE" "$URL" 2>/dev/null && gunzip -f "$OUTFILE" 2>/dev/null
    fi

    # Remove if empty (failed download)
    if [ -f "$FASTA" ] && [ ! -s "$FASTA" ]; then
        rm -f "$FASTA"
    fi
}
export -f download_genome
export GENOME_DIR

# Use xargs for parallel downloads (portable across Linux/macOS)
PROGRESS_FILE="$DOWNLOAD_DIR/progress.log"
> "$PROGRESS_FILE"

if command -v parallel &>/dev/null; then
    # GNU parallel available — use it for best performance
    parallel -j "$PARALLEL_DOWNLOADS" --bar --joblog "$DOWNLOAD_DIR/parallel.log" \
        download_genome {} < "$URL_LIST"
else
    # Fallback: xargs-based parallelism with progress counter
    DONE=0
    while IFS= read -r URL; do
        # Limit concurrent background jobs
        while [ "$(jobs -r | wc -l)" -ge "$PARALLEL_DOWNLOADS" ]; do
            sleep 0.5
        done

        (
            download_genome "$URL"
            DONE=$((DONE + 1))
            if [ $((DONE % 1000)) -eq 0 ]; then
                echo "  Downloaded $DONE / $TOTAL_GENOMES ..."
            fi
        ) &
    done < "$URL_LIST"

    # Wait for remaining jobs
    wait
fi

DOWNLOAD_END=$(date +%s)
DOWNLOAD_ELAPSED=$(( (DOWNLOAD_END - DOWNLOAD_START) / 60 ))

GENOME_COUNT=$(find "$GENOME_DIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | wc -l | tr -d ' ')
echo ""
echo "Downloaded $GENOME_COUNT genome files in ${DOWNLOAD_ELAPSED} minutes"
echo "Genome directory size: $(du -sh "$GENOME_DIR" | awk '{print $1}')"

FAILED=$((TOTAL_GENOMES - GENOME_COUNT))
if [ "$FAILED" -gt 0 ]; then
    echo "WARNING: $FAILED genomes failed to download"
    echo "Re-running may recover some. Continuing with $GENOME_COUNT genomes."
fi

if [ "$GENOME_COUNT" -eq 0 ]; then
    echo "ERROR: No genome FASTA files found after download."
    exit 1
fi

# ===========================================================================
# Step 3: Build Dragon index
# ===========================================================================
echo ""
echo "--- Step 3: Building Dragon index ---"
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
# Step 4: Summary
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
# Step 5: Package for distribution
# ===========================================================================
echo ""
echo "--- Packaging index for distribution ---"
MODE="complete"
if [ "$REPRESENTATIVE_ONLY" = true ]; then MODE="representative"; fi

TARBALL="$WORKSPACE/dragon_refseq_${MODE}_k${KMER_SIZE}_$(date +%Y%m%d).tar.gz"
tar -czf "$TARBALL" -C "$(dirname "$INDEX_DIR")" "$(basename "$INDEX_DIR")"
TARBALL_SIZE=$(du -sh "$TARBALL" | awk '{print $1}')
echo "Tarball: $TARBALL ($TARBALL_SIZE)"

# Write metadata
cat > "$WORKSPACE/index_metadata.json" <<METAEOF
{
    "database": "RefSeq bacteria",
    "mode": "$MODE",
    "genome_count": $GENOME_COUNT,
    "total_available": $TOTAL_GENOMES,
    "failed_downloads": $FAILED,
    "kmer_size": $KMER_SIZE,
    "index_size": "$INDEX_SIZE",
    "tarball_size": "$TARBALL_SIZE",
    "download_time_minutes": $DOWNLOAD_ELAPSED,
    "build_time_minutes": $BUILD_ELAPSED,
    "build_date": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "dragon_version": "$($DRAGON --version 2>/dev/null || echo unknown)",
    "assembly_summary_url": "$ASSEMBLY_SUMMARY_URL",
    "low_memory": $LOW_MEMORY
}
METAEOF
echo "Metadata: $WORKSPACE/index_metadata.json"

echo ""
echo "====================================="
echo "Done! RefSeq bacteria index built successfully."
echo ""
echo "Upload to Zenodo/S3:"
echo "  $TARBALL"
echo ""
echo "Users can download with:"
echo "  dragon download --database refseq-bacteria --output ./my_index"
echo ""
echo "Or search directly:"
echo "  dragon search --index $INDEX_DIR --query <reads.fasta>"
echo "====================================="
