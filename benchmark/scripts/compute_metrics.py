#!/usr/bin/env python3
"""Compute accuracy metrics (sensitivity, precision, F1) from search results."""

import argparse
import csv
import os


def parse_paf(path):
    """Parse PAF format results."""
    hits = {}
    if not os.path.exists(path):
        return hits
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            query = parts[0]
            target = parts[5]
            identity = int(parts[9]) / max(int(parts[10]), 1)
            if query not in hits:
                hits[query] = []
            hits[query].append({
                "target": target,
                "identity": identity,
                "query_start": int(parts[2]),
                "query_end": int(parts[3]),
                "target_start": int(parts[7]),
                "target_end": int(parts[8]),
            })
    return hits


def parse_blast6(path):
    """Parse BLAST tabular format results."""
    hits = {}
    if not os.path.exists(path):
        return hits
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            query = parts[0]
            target = parts[1]
            identity = float(parts[2]) / 100.0
            if query not in hits:
                hits[query] = []
            hits[query].append({
                "target": target,
                "identity": identity,
            })
    return hits


def parse_truth(path):
    """Parse ground truth TSV."""
    truth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            query_id = row["query_id"]
            truth[query_id] = {
                "genome": row["genome"],
                "contig": row.get("contig", ""),
                "start": int(row.get("start", 0)),
                "end": int(row.get("end", 0)),
            }
    return truth


def compute_metrics(hits, truth):
    """Compute sensitivity, precision, and F1."""
    tp, fp, fn = 0, 0, 0

    for query_id, true_info in truth.items():
        true_genome = true_info["genome"]

        if query_id in hits and len(hits[query_id]) > 0:
            # Check if any hit matches the true genome
            found = any(
                true_genome in h["target"] for h in hits[query_id]
            )
            if found:
                tp += 1
            else:
                fp += len(hits[query_id])
                fn += 1
        else:
            fn += 1

    # Count extra hits for queries not in truth
    for query_id in hits:
        if query_id not in truth:
            fp += len(hits[query_id])

    sensitivity = tp / max(tp + fn, 1)
    precision = tp / max(tp + fp, 1)
    f1 = 2 * precision * sensitivity / max(precision + sensitivity, 1e-10)

    return {
        "sensitivity": sensitivity,
        "precision": precision,
        "f1_score": f1,
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", nargs="+", required=True)
    parser.add_argument("--truths", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, "w") as out:
        out.write("tool\tdataset\tdivergence\tsensitivity\tprecision\tf1_score\t"
                  "tp\tfp\tfn\n")

        for result_path, truth_path in zip(args.results, args.truths):
            # Extract divergence from filename
            div = "0.0"
            for part in os.path.basename(result_path).split("_"):
                if part.startswith("div"):
                    div = part.replace("div", "").replace(".out", "")

            # Auto-detect format
            hits = parse_paf(result_path)
            if not hits:
                hits = parse_blast6(result_path)

            truth = parse_truth(truth_path)
            metrics = compute_metrics(hits, truth)

            out.write(
                f"{args.tool}\t{args.dataset}\t{div}\t"
                f"{metrics['sensitivity']:.6f}\t{metrics['precision']:.6f}\t"
                f"{metrics['f1_score']:.6f}\t"
                f"{metrics['true_positives']}\t{metrics['false_positives']}\t"
                f"{metrics['false_negatives']}\n"
            )

    print(f"Metrics written to {args.output}")


if __name__ == "__main__":
    main()
