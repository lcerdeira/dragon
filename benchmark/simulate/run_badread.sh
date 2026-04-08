#!/bin/bash
# Simulate ONT long reads using Badread.
# Usage: bash run_badread.sh <genome_dir> <output.fa>

set -euo pipefail

GENOME_DIR=$1
OUTPUT=$2
NUM_GENOMES=50
MEAN_LENGTH=5000
IDENTITY="85,99"
CHIMERA_RATE=1

mkdir -p "$(dirname "$OUTPUT")"

# Concatenate selected genomes
TEMP_REF=$(mktemp /tmp/badread_ref_XXXXXX.fa)
GENOME_FILES=($(find "$GENOME_DIR" -name "*.fa" -o -name "*.fasta" -o -name "*.fna" | sort | head -n "$NUM_GENOMES"))

echo "Simulating long reads from ${#GENOME_FILES[@]} genomes..."

cat "${GENOME_FILES[@]}" > "$TEMP_REF"

# Run Badread
badread simulate \
    --reference "$TEMP_REF" \
    --quantity 10x \
    --length "$MEAN_LENGTH,2000" \
    --identity "$IDENTITY,4" \
    --chimeras "$CHIMERA_RATE" \
    --junk_reads 1 \
    --random_reads 1 \
    --seed 42 \
    > "$OUTPUT"

rm -f "$TEMP_REF"

echo "Long reads written to $OUTPUT"
echo "Reads: $(grep -c '^>' "$OUTPUT")"
