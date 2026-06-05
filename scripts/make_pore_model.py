#!/usr/bin/env python3
"""
Fetch ONT's real R10.4.1 9-mer pore model and convert it to Dragon's JSON format
for `dragon signal-index --pore-model`.

ONT publishes k-mer level models at github.com/nanoporetech/kmer_models. The
R10.4.1 (e8.2, 400 bps) model is a 9-mer table of 4^9 = 262,144 entries, one
normalized current level per 9-mer, in lexicographic (= base-4: A<C<G<T) order —
which is exactly the index order Dragon uses, so the level column maps directly
to Dragon's flat `levels` array.

Usage:
    python scripts/make_pore_model.py -o r10.4.1_9mer.json
    dragon signal-index -i genomes/ -o sig_index/ --pore-model r10.4.1_9mer.json

This replaces Dragon's built-in synthetic pore model with the real ONT physics,
which is required for searching REAL nanopore reads (the synthetic model only
round-trips against simulated signal generated from itself).
"""
import argparse, json, os, sys, urllib.request

URL = ("https://raw.githubusercontent.com/nanoporetech/kmer_models/master/"
       "dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt")
ENC = {"A": 0, "C": 1, "G": 2, "T": 3}

def kmer_index(k):
    n = 0
    for c in k:
        n = n * 4 + ENC[c]
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="r10.4.1_9mer.json")
    ap.add_argument("--url", default=URL)
    args = ap.parse_args()

    print(f"downloading {args.url}", file=sys.stderr)
    txt = urllib.request.urlopen(args.url).read().decode()
    rows = [ln.split("\t") for ln in txt.splitlines() if ln.strip()]
    k = len(rows[0][0])
    expected = 4 ** k
    if len(rows) != expected:
        sys.exit(f"expected {expected} rows for {k}-mer, got {len(rows)}")

    # Verify file order == base-4 index order, then take the level column directly.
    levels = [0.0] * expected
    for kmer, lvl in rows:
        levels[kmer_index(kmer)] = float(lvl)
    # sanity: first/last rows should already be in order
    if kmer_index(rows[0][0]) != 0 or kmer_index(rows[-1][0]) != expected - 1:
        print("warning: file not in base-4 order; reindexed by k-mer", file=sys.stderr)

    model = {"kmer_size": k, "levels": levels,
             "name": "dna_r10.4.1_e8.2_400bps_9mer_v1"}
    json.dump(model, open(args.out, "w"))
    print(f"wrote {args.out} ({os.path.getsize(args.out)/1e6:.1f} MB, k={k}, "
          f"{len(levels)} levels)", file=sys.stderr)

if __name__ == "__main__":
    main()
