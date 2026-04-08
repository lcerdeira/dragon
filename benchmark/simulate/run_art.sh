#!/bin/bash
# Simulate Illumina short reads using ART.
# Usage: bash run_art.sh <genome_dir> <output_prefix>

set -euo pipefail

GENOME_DIR=$1
OUTPUT_PREFIX=$2
NUM_GENOMES=50
READ_LENGTH=150
COVERAGE=10
FRAGMENT_LEN=350
FRAGMENT_SD=50

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

# Concatenate selected genomes
TEMP_REF=$(mktemp /tmp/art_ref_XXXXXX.fa)
GENOME_FILES=($(find "$GENOME_DIR" -name "*.fa" -o -name "*.fasta" -o -name "*.fna" | sort | head -n "$NUM_GENOMES"))

echo "Simulating short reads from ${#GENOME_FILES[@]} genomes..."

cat "${GENOME_FILES[@]}" > "$TEMP_REF"

# Run ART
art_illumina \
    -ss HS25 \
    -i "$TEMP_REF" \
    -p \
    -l "$READ_LENGTH" \
    -f "$COVERAGE" \
    -m "$FRAGMENT_LEN" \
    -s "$FRAGMENT_SD" \
    -o "${OUTPUT_PREFIX}" \
    -rs 42

rm -f "$TEMP_REF"

echo "Short reads written to ${OUTPUT_PREFIX}1.fq and ${OUTPUT_PREFIX}2.fq"
