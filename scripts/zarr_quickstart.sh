#!/usr/bin/env bash
# Dragon Zarr quickstart — build a tiny demo index, export it to a Zarr store,
# and query it both locally and (optionally) from public S3 over HTTPS.
#
#   bash scripts/zarr_quickstart.sh
#
# Requires: the `dragon` binary (cargo build --release) and Python 3 (stdlib only).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DRAGON="${DRAGON:-$ROOT/target/release/dragon}"
WORK="${WORK:-$ROOT/demo_run}"

[ -x "$DRAGON" ] || { echo "dragon binary not found at $DRAGON — run 'cargo build --release' or set \$DRAGON"; exit 1; }

echo "==> 1/4  generating demo data (6 small genomes + 2 queries)"
python3 "$ROOT/scripts/make_demo_data.py" -o "$WORK"

echo "==> 2/4  building index"
"$DRAGON" index -i "$WORK/genomes" -o "$WORK/index" -k 31 -j 4 >/dev/null 2>&1

echo "==> 3/4  exporting to Zarr (chunked + Zstd)"
"$DRAGON" export-zarr -i "$WORK/index" -o "$WORK/index.zarr" >/dev/null 2>&1

echo "==> 4/4  querying the Zarr store (local)"
echo "--- dragon search-zarr (cloud-native seed discovery) ---"
"$DRAGON" search-zarr -z "$WORK/index.zarr" -q "$WORK/query.fa" 2>/dev/null
echo
echo "--- dragon search (full base-level alignment) ---"
"$DRAGON" search -i "$WORK/index" -q "$WORK/query.fa" --min-identity 0.0 --min-query-coverage 0.0 -o - 2>/dev/null \
  | grep -v '^#' | awk -F'\t' '{printf "  %-16s -> %-10s id=%.3f\n",$1,$6,$10/$11}'

cat <<'EOF'

Expected: 'core_fragment' hits ALL 6 genomes (shared core), 'resistance_gene'
hits ONLY genome_3 and genome_5 (the two carriers).

To test the SAME store from public S3 over HTTPS (no credentials, no download):

  dragon search-zarr \
    --zarr https://dragon-zarr.s3.eu-west-2.amazonaws.com/demo/index.zarr \
    -q demo_run/query.fa

To publish your own store:  dragon export-zarr ... && aws s3 cp --recursive my.zarr s3://your-bucket/my.zarr
EOF
