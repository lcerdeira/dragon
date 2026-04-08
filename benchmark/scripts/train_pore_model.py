#!/usr/bin/env python3
"""
Train a data-driven pore model from aligned ONT data.

Learns a 5-mer lookup table (1024 entries) mapping each DNA 5-mer to its
median observed pA current level, based on aligned nanopore reads.

Usage:
    1. Align ONT reads to reference:
       minimap2 -a ref.fa reads.fastq | samtools sort > aligned.bam

    2. Extract signal-to-reference alignment (requires f5c or nanopolish):
       f5c eventalign --rna -b aligned.bam -g ref.fa -r reads.fastq > eventalign.tsv

    3. Train the pore model:
       python3 train_pore_model.py --eventalign eventalign.tsv --output pore_model.json

    4. Use with Dragon:
       dragon signal-index --input genomes/ --output idx/ --pore-model pore_model.json
"""

import argparse
import json
import sys
from collections import defaultdict

import numpy as np


BASES = "ACGT"
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def kmer_to_index(kmer):
    """Encode a DNA k-mer as a base-4 integer."""
    idx = 0
    for base in kmer.upper():
        if base not in BASE_TO_IDX:
            return None
        idx = idx * 4 + BASE_TO_IDX[base]
    return idx


def index_to_kmer(idx, k=5):
    """Decode a base-4 integer back to a DNA k-mer."""
    kmer = []
    for _ in range(k):
        kmer.append(BASES[idx % 4])
        idx //= 4
    return "".join(reversed(kmer))


def load_eventalign(path, max_records=10_000_000):
    """Load nanopolish/f5c eventalign TSV.

    Expected columns: contig, position, reference_kmer, model_kmer, event_level_mean, ...
    We need: reference_kmer (or model_kmer) and event_level_mean.
    """
    kmer_signals = defaultdict(list)
    count = 0

    with open(path) as f:
        header = f.readline().strip().split("\t")

        # Find column indices
        kmer_col = None
        signal_col = None
        for i, name in enumerate(header):
            if name in ("reference_kmer", "model_kmer"):
                kmer_col = i
            if name == "event_level_mean":
                signal_col = i

        if kmer_col is None or signal_col is None:
            print(f"ERROR: Could not find required columns in {path}")
            print(f"  Found columns: {header}")
            print("  Need: reference_kmer/model_kmer and event_level_mean")
            sys.exit(1)

        for line in f:
            fields = line.strip().split("\t")
            if len(fields) <= max(kmer_col, signal_col):
                continue

            kmer = fields[kmer_col].upper()
            try:
                signal = float(fields[signal_col])
            except ValueError:
                continue

            if len(kmer) != 5:
                continue

            idx = kmer_to_index(kmer)
            if idx is None:
                continue

            kmer_signals[kmer].append(signal)
            count += 1

            if count >= max_records:
                break

    return kmer_signals


def load_simple_tsv(path):
    """Load a simple TSV with columns: kmer, pA_level."""
    kmer_signals = defaultdict(list)
    with open(path) as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            kmer = fields[0].upper()
            try:
                signal = float(fields[1])
            except ValueError:
                continue
            kmer_signals[kmer].append(signal)
    return kmer_signals


def build_model_table(kmer_signals, k=5):
    """Compute median pA per k-mer from collected signal observations."""
    num_kmers = 4 ** k
    levels = np.zeros(num_kmers, dtype=np.float32)
    coverage = 0

    for idx in range(num_kmers):
        kmer = index_to_kmer(idx, k)
        if kmer in kmer_signals and len(kmer_signals[kmer]) > 0:
            levels[idx] = np.median(kmer_signals[kmer])
            coverage += 1
        else:
            # Fill missing k-mers with a default (will be interpolated)
            levels[idx] = 90.0  # approximate baseline

    return levels, coverage


def interpolate_missing(levels, kmer_signals, k=5):
    """Fill in gaps using nearest-neighbor averaging for k-mers with no data."""
    num_kmers = 4 ** k
    filled = levels.copy()

    for idx in range(num_kmers):
        kmer = index_to_kmer(idx, k)
        if kmer in kmer_signals and len(kmer_signals[kmer]) > 0:
            continue

        # Average over all single-base substitution neighbors
        neighbor_levels = []
        for pos in range(k):
            for base in BASES:
                neighbor = kmer[:pos] + base + kmer[pos + 1:]
                if neighbor != kmer and neighbor in kmer_signals and len(kmer_signals[neighbor]) > 0:
                    nidx = kmer_to_index(neighbor)
                    neighbor_levels.append(levels[nidx])

        if neighbor_levels:
            filled[idx] = np.mean(neighbor_levels)

    return filled


def main():
    parser = argparse.ArgumentParser(description="Train Dragon pore model from ONT data")
    parser.add_argument("--eventalign", help="nanopolish/f5c eventalign TSV")
    parser.add_argument("--simple-tsv", help="Simple TSV with columns: kmer, pA_level")
    parser.add_argument("--output", default="pore_model.json", help="Output model JSON")
    parser.add_argument("--kmer-size", type=int, default=5, help="K-mer size (default: 5)")
    parser.add_argument("--max-records", type=int, default=10_000_000,
                        help="Max records to read from eventalign")
    args = parser.parse_args()

    if not args.eventalign and not args.simple_tsv:
        print("ERROR: Provide --eventalign or --simple-tsv")
        sys.exit(1)

    k = args.kmer_size
    num_kmers = 4 ** k

    if args.eventalign:
        print(f"Loading eventalign data from {args.eventalign}...")
        kmer_signals = load_eventalign(args.eventalign, args.max_records)
    else:
        print(f"Loading simple TSV from {args.simple_tsv}...")
        kmer_signals = load_simple_tsv(args.simple_tsv)

    total_obs = sum(len(v) for v in kmer_signals.values())
    print(f"  {len(kmer_signals)} unique k-mers, {total_obs} total observations")

    print("Computing median pA per k-mer...")
    levels, coverage = build_model_table(kmer_signals, k)
    print(f"  {coverage}/{num_kmers} k-mers have data ({100*coverage/num_kmers:.1f}%)")

    if coverage < num_kmers:
        print("Interpolating missing k-mers from neighbors...")
        levels = interpolate_missing(levels, kmer_signals, k)

    # Stats
    print(f"\nModel statistics:")
    print(f"  Range: [{levels.min():.1f}, {levels.max():.1f}] pA")
    print(f"  Mean:  {levels.mean():.1f} pA")
    print(f"  Std:   {levels.std():.1f} pA")

    # Export as full PoreModel JSON
    model = {
        "kmer_size": k,
        "levels": [float(l) for l in levels],
        "name": f"learned-{k}mer-{coverage}obs",
    }

    with open(args.output, "w") as f:
        json.dump(model, f, indent=2)
    print(f"\nPore model saved to {args.output}")
    print(f"Use with: dragon signal-index --pore-model {args.output} ...")


if __name__ == "__main__":
    main()
