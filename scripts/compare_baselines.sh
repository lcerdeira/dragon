#!/usr/bin/env bash
# compare_baselines.sh — Reproducible baseline comparison for Dragon paper
#
# Runs Dragon, LexicMap, Minimap2, and BLAST on the same queries and index,
# collects timing + accuracy, and produces a comparison table.
#
# Prerequisites:
#   - dragon (built: cargo build --release)
#   - lexicmap (https://github.com/shenwei356/lexicmap)
#   - minimap2 (https://github.com/lh3/minimap2)
#   - blastn   (NCBI BLAST+)
#
# Usage:
#   bash scripts/compare_baselines.sh \
#     --index /path/to/dragon_index \
#     --genomes /path/to/genomes/ \
#     --query /path/to/queries.fa \
#     --output /path/to/results/

set -euo pipefail

# ---- Defaults ----
DRAGON="${DRAGON:-dragon}"
INDEX=""
GENOMES=""
QUERY=""
OUTPUT="baseline_results"
THREADS=4

# ---- Parse args ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --index) INDEX="$2"; shift 2 ;;
        --genomes) GENOMES="$2"; shift 2 ;;
        --query) QUERY="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --dragon) DRAGON="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ -z "$INDEX" ] || [ -z "$QUERY" ]; then
    echo "Usage: $0 --index <dragon_index> --query <queries.fa> [--genomes <genomes_dir>] [--output <dir>]"
    exit 1
fi

mkdir -p "$OUTPUT"

echo "================================================="
echo " Dragon Baseline Comparison"
echo "================================================="
echo "Index:   $INDEX"
echo "Query:   $QUERY"
echo "Output:  $OUTPUT"
echo "Threads: $THREADS"
echo ""

# ---- Helper: time a command and capture peak RAM ----
time_cmd() {
    local name="$1"
    shift
    local time_file="$OUTPUT/${name}_time.txt"

    echo "--- Running $name ---"
    # GNU time for peak RSS; fall back to bash time
    if command -v gtime &>/dev/null; then
        gtime -v "$@" 2>"$time_file" || true
        wall=$(grep "wall clock" "$time_file" | awk '{print $NF}')
        ram=$(grep "Maximum resident" "$time_file" | awk '{print $NF}')
        echo "  Wall time: $wall"
        echo "  Peak RAM:  $ram KB"
    elif command -v /usr/bin/time &>/dev/null; then
        /usr/bin/time -l "$@" 2>"$time_file" || true
        ram=$(grep "maximum resident" "$time_file" | awk '{print $1}')
        echo "  Peak RAM: $ram bytes"
    else
        time "$@" 2>"$time_file" || true
    fi
    echo ""
}

# ====================================================================
# 1. DRAGON
# ====================================================================
echo "=== DRAGON ==="
time_cmd "dragon" "$DRAGON" search \
    --index "$INDEX" \
    --query "$QUERY" \
    --output "$OUTPUT/dragon.paf" \
    --threads "$THREADS" \
    --min-identity 0.0 \
    --min-query-coverage 0.0

# Generate surveillance summary
"$DRAGON" summarize \
    --input "$OUTPUT/dragon.paf" \
    --index "$INDEX" \
    --format tsv \
    --output "$OUTPUT/dragon_summary.tsv" 2>/dev/null || true

"$DRAGON" summarize \
    --input "$OUTPUT/dragon.paf" \
    --index "$INDEX" \
    --format json \
    --output "$OUTPUT/dragon_summary.json" 2>/dev/null || true

# Generate GFA context
"$DRAGON" search \
    --index "$INDEX" \
    --query "$QUERY" \
    --output "$OUTPUT/dragon_context.gfa" \
    --format gfa \
    --threads "$THREADS" \
    --min-identity 0.0 \
    --min-query-coverage 0.0 2>/dev/null || true

DRAGON_HITS=$(grep -c "^[^#@]" "$OUTPUT/dragon.paf" 2>/dev/null || echo 0)
echo "Dragon hits: $DRAGON_HITS"

# ====================================================================
# 2. LEXICMAP (if available)
# ====================================================================
if command -v lexicmap &>/dev/null; then
    echo "=== LEXICMAP ==="

    LEXICMAP_INDEX="$OUTPUT/lexicmap_index"
    if [ -n "$GENOMES" ] && [ ! -d "$LEXICMAP_INDEX" ]; then
        echo "Building LexicMap index..."
        time_cmd "lexicmap_index" lexicmap index \
            -I "$GENOMES" \
            -O "$LEXICMAP_INDEX" \
            -j "$THREADS"
    fi

    if [ -d "$LEXICMAP_INDEX" ]; then
        time_cmd "lexicmap_search" lexicmap search \
            -d "$LEXICMAP_INDEX" \
            -q "$QUERY" \
            -o "$OUTPUT/lexicmap.tsv" \
            -j "$THREADS"
        LEXICMAP_HITS=$(wc -l < "$OUTPUT/lexicmap.tsv" 2>/dev/null | tr -d ' ')
        echo "LexicMap hits: $LEXICMAP_HITS"
    fi
else
    echo "=== LEXICMAP: not found, skipping ==="
fi

# ====================================================================
# 3. MINIMAP2 (if available)
# ====================================================================
if command -v minimap2 &>/dev/null && [ -n "$GENOMES" ]; then
    echo "=== MINIMAP2 ==="

    # Concatenate all genomes for minimap2
    CONCAT_REF="$OUTPUT/concat_ref.fa"
    if [ ! -f "$CONCAT_REF" ]; then
        echo "Concatenating genomes for minimap2..."
        cat "$GENOMES"/*.{fa,fasta,fna} > "$CONCAT_REF" 2>/dev/null || true
    fi

    time_cmd "minimap2" minimap2 \
        -t "$THREADS" \
        -c \
        "$CONCAT_REF" "$QUERY" \
        > "$OUTPUT/minimap2.paf" 2>/dev/null

    MINIMAP2_HITS=$(grep -c "^[^#@]" "$OUTPUT/minimap2.paf" 2>/dev/null || echo 0)
    echo "Minimap2 hits: $MINIMAP2_HITS"
else
    echo "=== MINIMAP2: not found or no --genomes, skipping ==="
fi

# ====================================================================
# 4. BLAST (if available)
# ====================================================================
if command -v blastn &>/dev/null && [ -n "$GENOMES" ]; then
    echo "=== BLAST ==="

    BLAST_DB="$OUTPUT/blast_db"
    CONCAT_REF="$OUTPUT/concat_ref.fa"
    if [ ! -f "${BLAST_DB}.ndb" ] && [ -f "$CONCAT_REF" ]; then
        echo "Building BLAST database..."
        time_cmd "blast_index" makeblastdb -in "$CONCAT_REF" -dbtype nucl -out "$BLAST_DB"
    fi

    if [ -f "${BLAST_DB}.ndb" ]; then
        time_cmd "blastn" blastn \
            -db "$BLAST_DB" \
            -query "$QUERY" \
            -outfmt 6 \
            -num_threads "$THREADS" \
            -evalue 1e-5 \
            -out "$OUTPUT/blast.tsv"
        BLAST_HITS=$(wc -l < "$OUTPUT/blast.tsv" 2>/dev/null | tr -d ' ')
        echo "BLAST hits: $BLAST_HITS"
    fi
else
    echo "=== BLAST: not found or no --genomes, skipping ==="
fi

# ====================================================================
# 5. COMPARISON TABLE
# ====================================================================
echo ""
echo "================================================="
echo " Results Summary"
echo "================================================="
echo ""

# Collect metrics
{
    echo -e "tool\thits\tindex_size\twall_time\tpeak_ram_kb"

    # Dragon
    idx_size=$(du -sk "$INDEX" 2>/dev/null | awk '{print $1}' || echo "N/A")
    dragon_wall=$(grep "wall clock\|real" "$OUTPUT/dragon_time.txt" 2>/dev/null | head -1 | awk '{print $NF}' || echo "N/A")
    dragon_ram=$(grep -i "maximum resident\|Maximum resident" "$OUTPUT/dragon_time.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
    echo -e "dragon\t$DRAGON_HITS\t${idx_size}K\t$dragon_wall\t$dragon_ram"

    # LexicMap
    if [ -f "$OUTPUT/lexicmap.tsv" ]; then
        lm_size=$(du -sk "$LEXICMAP_INDEX" 2>/dev/null | awk '{print $1}' || echo "N/A")
        lm_wall=$(grep "wall clock\|real" "$OUTPUT/lexicmap_search_time.txt" 2>/dev/null | head -1 | awk '{print $NF}' || echo "N/A")
        lm_ram=$(grep -i "maximum resident" "$OUTPUT/lexicmap_search_time.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
        echo -e "lexicmap\t${LEXICMAP_HITS:-0}\t${lm_size}K\t$lm_wall\t$lm_ram"
    fi

    # Minimap2
    if [ -f "$OUTPUT/minimap2.paf" ]; then
        mm_wall=$(grep "wall clock\|real" "$OUTPUT/minimap2_time.txt" 2>/dev/null | head -1 | awk '{print $NF}' || echo "N/A")
        mm_ram=$(grep -i "maximum resident" "$OUTPUT/minimap2_time.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
        echo -e "minimap2\t${MINIMAP2_HITS:-0}\tN/A\t$mm_wall\t$mm_ram"
    fi

    # BLAST
    if [ -f "$OUTPUT/blast.tsv" ]; then
        bl_wall=$(grep "wall clock\|real" "$OUTPUT/blastn_time.txt" 2>/dev/null | head -1 | awk '{print $NF}' || echo "N/A")
        bl_ram=$(grep -i "maximum resident" "$OUTPUT/blastn_time.txt" 2>/dev/null | awk '{print $NF}' || echo "N/A")
        echo -e "blastn\t${BLAST_HITS:-0}\tN/A\t$bl_wall\t$bl_ram"
    fi
} | tee "$OUTPUT/comparison.tsv"

echo ""
echo "Full results in: $OUTPUT/"
echo "  dragon.paf           — Dragon alignments"
echo "  dragon_summary.tsv   — Surveillance summary (TSV)"
echo "  dragon_summary.json  — Surveillance summary (JSON)"
echo "  dragon_context.gfa   — GFA subgraph context"
echo "  comparison.tsv       — Tool comparison table"
echo ""
echo "Done."
