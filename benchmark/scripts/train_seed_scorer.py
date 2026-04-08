#!/usr/bin/env python3
"""
Train a logistic regression seed scorer from Dragon --dump-seeds output.

Usage:
    1. Run Dragon with seed dump:
       dragon search --index <idx> --query <queries.fa> --dump-seeds seeds.tsv --no-ml

    2. Train the scorer:
       python3 train_seed_scorer.py --seeds seeds.tsv --truth truth.tsv --output seed_scorer.json

    3. Use trained weights:
       dragon search --index <idx> --query <queries.fa> --ml-weights seed_scorer.json
"""

import argparse
import json
import sys
from collections import defaultdict

import numpy as np


def load_seeds(path):
    """Load seed dump TSV."""
    records = []
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) != len(header):
                continue
            record = dict(zip(header, fields))
            records.append(record)
    return records


def load_truth(path):
    """Load truth TSV. Returns {query_id: genome_name}."""
    truth = {}
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            row = dict(zip(header, fields))
            query_id = row.get("query_id", fields[0])
            genome = row.get("genome", fields[1])
            truth[query_id] = genome
    return truth


def label_seeds(seeds, truth):
    """Label each seed as TP (1) or FP (0) based on truth genome mapping."""
    labeled = []
    for seed in seeds:
        query = seed["query_name"]
        genome_id = int(seed["genome_id"])

        # Check if this seed's genome matches the truth genome
        true_genome = truth.get(query)
        if true_genome is None:
            continue

        # genome_id is an integer index; truth genome is a name.
        # Try matching by checking if the genome name contains the ID or vice versa.
        # In practice, the benchmark truth has the genome file stem.
        label = 1 if str(genome_id) in true_genome or true_genome in str(genome_id) else 0

        labeled.append({
            "match_len": float(seed["match_len"]),
            "log_sa_count": np.log2(max(float(seed["sa_count"]), 1)),
            "query_pos_frac": 0.5,  # placeholder (not in dump currently)
            "match_frac": 0.0,  # placeholder
            "color_card": float(seed["color_card"]),
            "gc_content": float(seed["gc_content"]),
            "label": label,
        })

    return labeled


def train(labeled):
    """Train logistic regression on labeled seed features."""
    try:
        from sklearn.linear_model import LogisticRegression  # noqa
    except ImportError:
        print("ERROR: scikit-learn is required. Install with: pip install scikit-learn")
        sys.exit(1)

    if not labeled:
        print("ERROR: No labeled seeds to train on.")
        sys.exit(1)

    feature_names = ["match_len", "log_sa_count", "query_pos_frac", "match_frac", "color_card", "gc_content"]
    X = np.array([[s[f] for f in feature_names] for s in labeled])
    y = np.array([s["label"] for s in labeled])

    print(f"Training data: {len(labeled)} seeds, {sum(y)} TP, {len(y) - sum(y)} FP")
    print(f"TP rate: {sum(y) / len(y):.3f}")

    # Standardize features for better training
    mean = X.mean(axis=0)
    std = X.std(axis=0) + 1e-8
    X_norm = (X - mean) / std

    model = LogisticRegression(max_iter=1000, class_weight="balanced")
    model.fit(X_norm, y)

    # Convert to raw-feature weights by absorbing standardization:
    #   z = w_norm . (X - mean) / std + b_norm
    #     = (w_norm / std) . X + (b_norm - w_norm . mean / std)
    raw_weights = model.coef_[0] / std
    raw_bias = model.intercept_[0] - np.dot(model.coef_[0], mean / std)

    # Format as [bias, w1, w2, ..., w6]
    weights = [float(raw_bias)] + [float(w) for w in raw_weights]

    print(f"\nTrained weights:")
    names = ["bias"] + feature_names
    for name, w in zip(names, weights):
        print(f"  {name:20s}: {w:+.6f}")

    # Training accuracy
    y_pred = model.predict(X_norm)
    acc = np.mean(y_pred == y)
    print(f"\nTraining accuracy: {acc:.3f}")

    return weights


def main():
    parser = argparse.ArgumentParser(description="Train Dragon seed scorer")
    parser.add_argument("--seeds", required=True, help="Seed dump TSV from --dump-seeds")
    parser.add_argument("--truth", required=True, help="Ground truth TSV (query_id -> genome)")
    parser.add_argument("--output", default="seed_scorer.json", help="Output weights JSON")
    args = parser.parse_args()

    print("Loading seeds...")
    seeds = load_seeds(args.seeds)
    print(f"  {len(seeds)} seed records")

    print("Loading truth...")
    truth = load_truth(args.truth)
    print(f"  {len(truth)} query-genome mappings")

    print("Labeling seeds...")
    labeled = label_seeds(seeds, truth)
    print(f"  {len(labeled)} labeled seeds")

    if not labeled:
        print("ERROR: No seeds could be labeled. Check that query names match between --seeds and --truth.")
        sys.exit(1)

    print("\nTraining logistic regression...")
    weights = train(labeled)

    with open(args.output, "w") as f:
        json.dump(weights, f, indent=2)
    print(f"\nWeights saved to {args.output}")
    print(f"Use with: dragon search --ml-weights {args.output} ...")


if __name__ == "__main__":
    main()
