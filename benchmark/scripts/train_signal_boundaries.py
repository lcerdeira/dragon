#!/usr/bin/env python3
"""
Learn discretization boundaries for Dragon signal-level indexing.

Reads normalized signal values from TSV files (one value per line or multi-column
with a 'signal' column), fits k-means (or quantile-based) boundaries, and writes
a JSON array of boundary thresholds.

Usage:
    python3 train_signal_boundaries.py --signals signal_values.tsv --num-levels 16 --output boundaries.json

    # Then build signal index with learned boundaries:
    dragon signal-index --input genomes/ --output signal_idx/ --signal-boundaries boundaries.json
"""

import argparse
import json
import sys

import numpy as np


def load_signal_values(path):
    """Load normalized signal values from a TSV/CSV file."""
    values = []
    with open(path) as f:
        header = f.readline().strip()
        # Check if it's a header or a number
        try:
            values.append(float(header))
        except ValueError:
            pass  # skip header

        for line in f:
            line = line.strip()
            if not line:
                continue
            # Take first column if multi-column
            parts = line.split("\t")
            if not parts:
                parts = line.split(",")
            try:
                values.append(float(parts[0]))
            except ValueError:
                continue

    return np.array(values)


def quantile_boundaries(values, num_levels):
    """Compute boundaries using quantiles (equal-frequency bins)."""
    quantiles = np.linspace(0, 1, num_levels + 1)[1:-1]  # num_levels - 1 boundaries
    boundaries = np.quantile(values, quantiles)
    # Remove duplicates (can happen with discrete or concentrated data)
    boundaries = np.unique(boundaries)
    return boundaries


def kmeans_boundaries(values, num_levels):
    """Compute boundaries using k-means clustering."""
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        print("ERROR: scikit-learn required for k-means. Install: pip install scikit-learn")
        print("Falling back to quantile boundaries.")
        return quantile_boundaries(values, num_levels)

    # Subsample for speed if too many values
    if len(values) > 100_000:
        rng = np.random.RandomState(42)
        values = rng.choice(values, 100_000, replace=False)

    km = KMeans(n_clusters=num_levels, n_init=10, random_state=42)
    km.fit(values.reshape(-1, 1))

    centers = np.sort(km.cluster_centers_.ravel())
    # Boundaries are midpoints between consecutive cluster centers
    boundaries = (centers[:-1] + centers[1:]) / 2.0
    return boundaries


def main():
    parser = argparse.ArgumentParser(description="Learn signal discretization boundaries")
    parser.add_argument("--signals", required=True, help="File with normalized signal values")
    parser.add_argument("--num-levels", type=int, default=16, help="Number of discrete levels")
    parser.add_argument("--method", choices=["quantile", "kmeans"], default="quantile",
                        help="Boundary computation method")
    parser.add_argument("--output", default="boundaries.json", help="Output JSON file")
    args = parser.parse_args()

    print(f"Loading signal values from {args.signals}...")
    values = load_signal_values(args.signals)
    print(f"  {len(values)} values, range [{values.min():.3f}, {values.max():.3f}]")
    print(f"  mean={values.mean():.3f}, std={values.std():.3f}")

    if args.method == "kmeans":
        print(f"Computing k-means boundaries (k={args.num_levels})...")
        boundaries = kmeans_boundaries(values, args.num_levels)
    else:
        print(f"Computing quantile boundaries ({args.num_levels} levels)...")
        boundaries = quantile_boundaries(values, args.num_levels)

    print(f"  {len(boundaries)} boundaries:")
    for i, b in enumerate(boundaries):
        print(f"    [{i}] {b:+.4f}")

    with open(args.output, "w") as f:
        json.dump([float(b) for b in boundaries], f, indent=2)
    print(f"\nBoundaries saved to {args.output}")
    print(f"Use with: dragon signal-index --signal-boundaries {args.output} ...")


if __name__ == "__main__":
    main()
