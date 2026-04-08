#!/usr/bin/env python3
"""
Evaluate Dragon seed scorer by running search with and without ML scoring
on benchmark data and comparing precision/recall.

Usage:
    python3 evaluate_seed_scorer.py --index <idx> --queries <queries.fa> --truth <truth.tsv> \
        [--ml-weights seed_scorer.json]
"""

import argparse
import json
import subprocess
import sys
import tempfile
from collections import defaultdict


def run_dragon_search(index_dir, query_file, extra_args=None, output_format="paf"):
    """Run Dragon search and return PAF records."""
    cmd = [
        "cargo", "run", "--release", "--", "search",
        "--index", index_dir,
        "--query", query_file,
        "--format", output_format,
        "--output", "-",
    ]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        print(f"Dragon search failed: {result.stderr}", file=sys.stderr)
        return []

    records = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) >= 12:
            records.append({
                "query": fields[0],
                "target": fields[5],
                "identity": int(fields[9]) / max(int(fields[10]), 1),
                "mapq": int(fields[11]),
            })
    return records


def load_truth(path):
    """Load truth TSV: query_id -> genome_name."""
    truth = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                truth[parts[0]] = parts[1]
    return truth


def evaluate(records, truth):
    """Compute precision and recall from search results vs truth."""
    tp = 0
    fp = 0
    queries_found = set()

    for r in records:
        query = r["query"]
        target = r["target"]
        true_genome = truth.get(query)

        if true_genome is None:
            continue

        queries_found.add(query)

        if true_genome in target or target in true_genome:
            tp += 1
        else:
            fp += 1

    total_queries = len(truth)
    fn = total_queries - len(queries_found)

    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-9)

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "total_hits": len(records),
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate Dragon seed scorer")
    parser.add_argument("--index", required=True, help="Path to Dragon index")
    parser.add_argument("--queries", required=True, help="Query FASTA file")
    parser.add_argument("--truth", required=True, help="Truth TSV (query_id -> genome)")
    parser.add_argument("--ml-weights", help="Custom ML weights file")
    args = parser.parse_args()

    truth = load_truth(args.truth)
    print(f"Truth: {len(truth)} query-genome mappings")

    # Run WITHOUT ML scoring
    print("\n--- Without ML scoring (--no-ml) ---")
    records_noml = run_dragon_search(args.index, args.queries, ["--no-ml"])
    metrics_noml = evaluate(records_noml, truth)
    print(f"  Hits: {metrics_noml['total_hits']}")
    print(f"  TP: {metrics_noml['tp']}, FP: {metrics_noml['fp']}, FN: {metrics_noml['fn']}")
    print(f"  Precision: {metrics_noml['precision']:.3f}")
    print(f"  Recall:    {metrics_noml['recall']:.3f}")
    print(f"  F1:        {metrics_noml['f1']:.3f}")

    # Run WITH ML scoring
    extra = []
    if args.ml_weights:
        extra = ["--ml-weights", args.ml_weights]
    print("\n--- With ML scoring ---")
    records_ml = run_dragon_search(args.index, args.queries, extra)
    metrics_ml = evaluate(records_ml, truth)
    print(f"  Hits: {metrics_ml['total_hits']}")
    print(f"  TP: {metrics_ml['tp']}, FP: {metrics_ml['fp']}, FN: {metrics_ml['fn']}")
    print(f"  Precision: {metrics_ml['precision']:.3f}")
    print(f"  Recall:    {metrics_ml['recall']:.3f}")
    print(f"  F1:        {metrics_ml['f1']:.3f}")

    # Summary
    print("\n--- Comparison ---")
    dp = metrics_ml["precision"] - metrics_noml["precision"]
    dr = metrics_ml["recall"] - metrics_noml["recall"]
    df = metrics_ml["f1"] - metrics_noml["f1"]
    print(f"  Precision change: {dp:+.3f}")
    print(f"  Recall change:    {dr:+.3f}")
    print(f"  F1 change:        {df:+.3f}")


if __name__ == "__main__":
    main()
