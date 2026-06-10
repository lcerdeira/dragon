#!/usr/bin/env bash
# build_vfdb_index.sh — build a Dragon index of VFDB LOCALLY.
#
# VFDB (mgc.ac.cn/VFs) is "free for academic use" but lacks an explicit open-redistribution
# license, so (pending confirmation from the VFDB authors) we do not redistribute a prebuilt
# index. Build it locally here. If VFDB authors confirm redistribution-with-attribution is OK,
# this index can instead be published under CC-BY-4.0 alongside the others.
set -euo pipefail

DRAGON="${DRAGON:-$HOME/dragon/target/release/dragon}"
OUT="${1:-./vfdb_index}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
# setB = full dataset (nucleotide). Use VFDB_setA_nt.fas.gz for the core set.
VFDB_URL="${VFDB_URL:-http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz}"   # «verificar URL»

echo ">> downloading VFDB…"
curl -fSL "$VFDB_URL" -o "$WORK/vfdb.fas.gz"
gunzip -f "$WORK/vfdb.fas.gz"

echo ">> building Dragon index -> $OUT"
"$DRAGON" index -i "$WORK/vfdb.fas" -o "$OUT" -k 31 --threads "${THREADS:-8}"

echo ">> done. Index at: $OUT"
echo ">> cite: Chen L, et al. VFDB. http://www.mgc.ac.cn/VFs/"
