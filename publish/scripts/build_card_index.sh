#!/usr/bin/env bash
# build_card_index.sh — build a Dragon index of CARD LOCALLY.
#
# WHY THIS SCRIPT EXISTS:
#   CARD (card.mcmaster.ca) is distributed under a Data Usage Agreement that PROHIBITS
#   redistribution and restricts commercial use. We therefore do NOT ship a prebuilt CARD
#   index. You download CARD yourself (accepting its terms) and build the index locally.
#
# Usage:  CARD_OK=1 ./build_card_index.sh [outdir]
set -euo pipefail

DRAGON="${DRAGON:-$HOME/dragon/target/release/dragon}"
OUT="${1:-./card_index}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
CARD_URL="${CARD_URL:-https://card.mcmaster.ca/latest/data}"   # «verificar URL atual»

cat <<'EOF'
-------------------------------------------------------------------------------
CARD Data Usage Agreement: by running this you confirm you have read and accept
the CARD license at https://card.mcmaster.ca/about . The resulting index is for
YOUR use and must NOT be redistributed. Commercial use requires a CARD license.
-------------------------------------------------------------------------------
EOF
if [ "${CARD_OK:-0}" != "1" ]; then
  echo "Refusing to proceed. Re-run with CARD_OK=1 to confirm you accept CARD's terms." >&2
  exit 1
fi

echo ">> downloading CARD…"
curl -fSL "$CARD_URL" -o "$WORK/card-data.tar.bz2"
tar -xjf "$WORK/card-data.tar.bz2" -C "$WORK"

# Homolog-model nucleotide sequences (adjust if you want protein/variant models):
FASTA="$WORK/nucleotide_fasta_protein_homolog_model.fasta"
[ -f "$FASTA" ] || { echo "expected $FASTA not found; check CARD archive layout" >&2; ls "$WORK" >&2; exit 1; }

echo ">> building Dragon index -> $OUT"
# «confirmar sintaxe do Dragon para indexar um multi-FASTA de genes»
"$DRAGON" index -i "$FASTA" -o "$OUT" -k 31 --threads "${THREADS:-8}"

echo ">> done. Index at: $OUT  (LOCAL ONLY — do not redistribute)"
echo ">> optional: $DRAGON export-zarr -i $OUT -o ${OUT}.zarr"
