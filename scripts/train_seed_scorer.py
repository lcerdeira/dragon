#!/usr/bin/env python3
"""Train logistic regression seed scorer weights from Dragon seed dump data.

Usage:
    1. Generate training data:
       dragon search --index <index> --query <query.fa> \
           --dump-seeds seeds.tsv --ground-truth <genome_name> \
           --no-ml --min-identity 0.0 --min-query-coverage 0.0

    2. Train weights:
       python scripts/train_seed_scorer.py seeds.tsv -o seed_scorer.json

    3. Use trained weights:
       dragon search --index <index> --query <query.fa> --ml-weights seed_scorer.json

    Multiple TSV files can be concatenated for training across diverse queries.

Output: JSON array of 11 floats [bias, w1, ..., w10] ready for seed_scorer.json.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def load_data(tsv_paths: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Load and merge training data from one or more TSV files."""
    all_X = []
    all_y = []
    feature_names = None

    for path in tsv_paths:
        with open(path) as f:
            header = f.readline().strip().split("\t")

        # Find feature columns (between genome_name and label)
        # Header: query_name, genome_id, genome_name, <features...>, label
        label_idx = header.index("label")
        feat_start = 3  # after query_name, genome_id, genome_name
        feat_end = label_idx
        current_features = header[feat_start:feat_end]

        if feature_names is None:
            feature_names = current_features
        else:
            assert current_features == feature_names, (
                f"Feature mismatch: {current_features} vs {feature_names}"
            )

        # Load with numpy, skip header
        data = np.genfromtxt(
            path,
            delimiter="\t",
            skip_header=1,
            usecols=list(range(feat_start, feat_end)) + [label_idx],
            dtype=float,
        )

        if data.ndim == 1:
            data = data.reshape(1, -1)

        # Filter out unlabeled rows (label == -1)
        labeled = data[:, -1] >= 0
        data = data[labeled]

        if len(data) == 0:
            print(f"  Warning: {path} has no labeled data, skipping")
            continue

        X = data[:, :-1]
        y = data[:, -1]

        # Remove rows with NaN/inf
        valid = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
        X = X[valid]
        y = y[valid]

        all_X.append(X)
        all_y.append(y)
        print(f"  {path}: {len(X)} seeds ({int(y.sum())} positive, {int(len(y) - y.sum())} negative)")

    if not all_X:
        print("Error: no valid training data found")
        sys.exit(1)

    X = np.vstack(all_X)
    y = np.concatenate(all_y)
    print(f"\nTotal: {len(X)} seeds, {int(y.sum())} positive ({y.mean():.1%}), "
          f"{feature_names}")

    return X, y


def train_logistic_regression(
    X: np.ndarray,
    y: np.ndarray,
    lr: float = 0.01,
    epochs: int = 1000,
    l2_reg: float = 0.01,
) -> np.ndarray:
    """Train logistic regression via gradient descent (no sklearn dependency).

    Returns weights array of shape [1 + n_features] = [bias, w1, ..., wN].
    """
    n_samples, n_features = X.shape

    # Feature normalization (z-score) for stable training
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    sigma[sigma == 0] = 1.0  # avoid division by zero
    X_norm = (X - mu) / sigma

    # Add bias column
    X_b = np.column_stack([np.ones(n_samples), X_norm])

    # Initialize weights
    w = np.zeros(n_features + 1)

    # Class weights to handle imbalance (positive seeds are rare)
    n_pos = max(y.sum(), 1)
    n_neg = max(len(y) - n_pos, 1)
    pos_weight = len(y) / (2 * n_pos)
    neg_weight = len(y) / (2 * n_neg)
    sample_weights = np.where(y == 1, pos_weight, neg_weight)

    best_w = w.copy()
    best_loss = float("inf")

    for epoch in range(epochs):
        # Forward pass
        z = X_b @ w
        z = np.clip(z, -500, 500)  # prevent overflow
        p = 1.0 / (1.0 + np.exp(-z))

        # Weighted cross-entropy loss + L2 regularization
        eps = 1e-15
        loss = -np.mean(
            sample_weights * (y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps))
        ) + l2_reg * np.sum(w[1:] ** 2)

        if loss < best_loss:
            best_loss = loss
            best_w = w.copy()

        # Gradient
        grad = X_b.T @ (sample_weights * (p - y)) / n_samples
        grad[1:] += 2 * l2_reg * w[1:]  # L2 on non-bias

        w -= lr * grad

        if epoch % 200 == 0:
            acc = ((p > 0.5) == y).mean()
            print(f"  Epoch {epoch:4d}: loss={loss:.4f}, acc={acc:.3f}")

    # Convert back to original scale: w_orig = w_norm / sigma, bias adjusted
    w_orig = np.zeros(n_features + 1)
    w_orig[0] = best_w[0] - np.sum(best_w[1:] * mu / sigma)
    w_orig[1:] = best_w[1:] / sigma

    return w_orig


def evaluate(X: np.ndarray, y: np.ndarray, weights: np.ndarray) -> dict:
    """Evaluate trained weights on data."""
    X_b = np.column_stack([np.ones(len(X)), X])
    z = X_b @ weights
    p = 1.0 / (1.0 + np.exp(-np.clip(z, -500, 500)))
    pred = p > 0.5

    tp = ((pred == 1) & (y == 1)).sum()
    fp = ((pred == 1) & (y == 0)).sum()
    fn = ((pred == 0) & (y == 1)).sum()
    tn = ((pred == 0) & (y == 0)).sum()

    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-15)
    accuracy = (tp + tn) / len(y)

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "tn": int(tn),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Train Dragon seed scorer weights from labeled seed dumps"
    )
    parser.add_argument(
        "input",
        nargs="+",
        help="TSV file(s) from dragon search --dump-seeds --ground-truth",
    )
    parser.add_argument(
        "-o", "--output",
        default="seed_scorer.json",
        help="Output JSON weights file (default: seed_scorer.json)",
    )
    parser.add_argument("--lr", type=float, default=0.01, help="Learning rate")
    parser.add_argument("--epochs", type=int, default=2000, help="Training epochs")
    parser.add_argument("--l2", type=float, default=0.01, help="L2 regularization")

    args = parser.parse_args()

    print("Loading training data...")
    X, y = load_data(args.input)

    if y.sum() == 0:
        print("Error: no positive labels found. Check --ground-truth matches genome names.")
        sys.exit(1)

    print(f"\nTraining logistic regression ({X.shape[1]} features, {args.epochs} epochs)...")
    weights = train_logistic_regression(X, y, lr=args.lr, epochs=args.epochs, l2_reg=args.l2)

    # Evaluate
    print("\nEvaluation on training data:")
    metrics = evaluate(X, y, weights)
    for k, v in metrics.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
        else:
            print(f"  {k}: {v}")

    # Print weight interpretation
    print("\nLearned weights:")
    print(f"  bias: {weights[0]:.6f}")
    for i, name in enumerate(FEATURE_NAMES):
        print(f"  {name}: {weights[i+1]:+.6f}")

    # Save
    weights_list = [round(float(w), 8) for w in weights]
    with open(args.output, "w") as f:
        json.dump(weights_list, f, indent=2)
    print(f"\nSaved {len(weights_list)} weights to {args.output}")
    print(f"Copy to index directory or use: dragon search --ml-weights {args.output}")


# Feature names must match ml_score.rs FEATURE_NAMES exactly
FEATURE_NAMES = [
    "match_len",
    "log_sa_count",
    "query_pos_frac",
    "match_frac",
    "color_cardinality",
    "gc_content",
    "information_content",
    "inverse_sa_count",
    "local_seed_density",
    "seed_uniqueness",
]

if __name__ == "__main__":
    main()
