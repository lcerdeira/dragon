#!/usr/bin/env bash
# gen_manifest.sh — regenerate publish/manifest.csv from the live S3 bucket.
# Usage:  ./gen_manifest.sh [bucket] [region]   (defaults: dragon-zarr eu-west-2)
# Requires: awscli v2. Run on a host with bucket read access (e.g. loginhpc).
set -euo pipefail

BUCKET="${1:-dragon-zarr}"
REGION="${2:-eu-west-2}"
AWS="${AWS_BIN:-aws}"
OUT="${OUT:-manifest.csv}"

echo "index,type,s3_uri,region,objects,size_bytes,size_human,status" > "$OUT"

human() { # bytes -> human
  awk -v b="$1" 'BEGIN{u="B KB MB GB TB PB";split(u,a," ");i=1;while(b>=1024&&i<6){b/=1024;i++}printf "%.1f%s",b,a[i]}'
}

# top-level prefixes in the bucket
mapfile -t PREFIXES < <("$AWS" s3 ls "s3://$BUCKET/" --region "$REGION" | awk '/PRE /{print $2}' | sed 's#/##')

for p in "${PREFIXES[@]}"; do
  summ=$("$AWS" s3 ls "s3://$BUCKET/$p/" --region "$REGION" --recursive --summarize 2>/dev/null \
         | grep -E 'Total (Objects|Size)' || true)
  objs=$(awk -F': *' '/Total Objects/{print $2}' <<<"$summ" | tr -d ' ')
  size=$(awk -F': *' '/Total Size/{print $2}' <<<"$summ" | tr -d ' ')
  objs=${objs:-0}; size=${size:-0}
  status="published"; [ "$size" = "0" ] && status="empty"
  # crude type guess
  case "$p" in
    *index*|atb|gtdb*|kpneumo|saureus*|ecoli|shigella|salmonella*) type="genome_collection";;
    amrfinderplus*|card|vfdb|plasmidfinder) type="reference_db";;
    demo|smoke) type="test";;
    *) type="unknown";;
  esac
  printf '%s,%s,s3://%s/%s/,%s,%s,%s,%s,%s\n' \
    "$p" "$type" "$BUCKET" "$p" "$REGION" "$objs" "$size" "$(human "$size")" "$status" >> "$OUT"
  echo "  $p: $objs objects, $(human "$size")" >&2
done

echo "wrote $OUT" >&2

# ---- optional: generate + upload SHA256SUMS per published prefix ----
# for p in "${PREFIXES[@]}"; do
#   "$AWS" s3 cp --recursive "s3://$BUCKET/$p/" - --region "$REGION" 2>/dev/null \
#     | sha256sum  # (placeholder — real impl should checksum each object/chunk index)
# done
