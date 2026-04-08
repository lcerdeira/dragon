#!/usr/bin/env python3
"""
Train a signal scorer from Dragon signal search results.

The signal scorer replaces the hand-tuned compute_genome_score() formula
with a learned logistic regression over 4 features:
  num_hits, coverage, density, avg_match_len

Usage:
    1. Run Dragon signal search with ground truth to collect (genome_id, score, label) pairs:
       dragon signal-search --index <idx> --query <reads.tsv> --output results.tsv

    2. Prepare a truth TSV mapping read_id -> true_genome_name

    3. Train:
       python3 train_signal_scorer.py --results results.tsv --truth truth.tsv --output signal_scorer.json

    4. (Optional) Place signal_scorer.json in the signal index directory.
       Dragon will load it automatically during signal search.
"""

import argparse
import json
import sys

import numpy as np


def load_results(path):
    """Load signal search results TSV."""
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue
            records.append({
                "read_id": fields[0],
                "signal_length": int(fields[1]),
                "genome_name": fields[2],
                "genome_id": int(fields[3]),
                "position": int(fields[4]),
                "score": float(fields[5]),
                "strand": fields[6],
                "match_len": int(fields[7]),
            })
    return records


def load_truth(path):
    """Load truth TSV: read_id -> genome_name."""
    truth = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 2:
                truth[fields[0]] = fields[1]
    return truth


def main():
    parser = argparse.ArgumentParser(description="Train Dragon signal scorer")
    parser.add_argument("--results", required=True, help="Signal search results TSV")
    parser.add_argument("--truth", required=True, help="Ground truth TSV (read_id -> genome)")
    parser.add_argument("--output", default="signal_scorer.json", help="Output weights JSON")
    args = parser.parse_args()

    try:
        from sklearn.linear_model import LogisticRegression
    except ImportError:
        print("ERROR: scikit-learn required. Install: pip install scikit-learn")
        sys.exit(1)

    print("Loading results...")
    results = load_results(args.results)
    print(f"  {len(results)} result records")

    print("Loading truth...")
    truth = load_truth(args.truth)
    print(f"  {len(truth)} read-genome mappings")

    # Label each result
    features = []
    labels = []
    for r in results:
        true_genome = truth.get(r["read_id"])
        if true_genome is None:
            continue

        label = 1 if true_genome in r["genome_name"] or r["genome_name"] in true_genome else 0

        # Reconstruct features from result fields
        num_hits = r["score"]  # approximate: score ~ num_hits in default scorer
        coverage = 0.5  # placeholder (not in output)
        density = 0.5   # placeholder
        avg_match_len = r["match_len"] / max(num_hits, 1)

        features.append([num_hits, coverage, density, avg_match_len])
        labels.append(label)

    X = np.array(features)
    y = np.array(labels)

    print(f"\nLabeled data: {len(y)} records, {sum(y)} TP, {len(y) - sum(y)} FP")

    if len(y) == 0 or sum(y) == 0:
        print("ERROR: No positive labels found.")
        sys.exit(1)

    # Standardize
    mean = X.mean(axis=0)
    std = X.std(axis=0) + 1e-8
    X_norm = (X - mean) / std

    model = LogisticRegression(max_iter=1000, class_weight="balanced")
    model.fit(X_norm, y)

    # Convert to raw weights (absorb standardization)
    raw_weights = model.coef_[0] / std
    raw_bias = model.intercept_[0] - np.dot(model.coef_[0], mean / std)

    weights = [float(raw_bias)] + [float(w) for w in raw_weights]

    print("\nTrained weights:")
    names = ["bias", "num_hits", "coverage", "density", "avg_match_len"]
    for name, w in zip(names, weights):
        print(f"  {name:20s}: {w:+.6f}")

    y_pred = model.predict(X_norm)
    acc = np.mean(y_pred == y)
    print(f"\nTraining accuracy: {acc:.3f}")

    with open(args.output, "w") as f:
        json.dump(weights, f, indent=2)
    print(f"\nWeights saved to {args.output}")


if __name__ == "__main__":
    main()
