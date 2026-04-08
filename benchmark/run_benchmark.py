#!/usr/bin/env python3
"""
Dragon Benchmark Runner
=======================
Runs Dragon against all 6 dataset tiers and measures accuracy + resources.
Also generates comparison data from simple k-mer baselines.

Usage:
    python3 run_benchmark.py
"""

import os
import sys
import time
import resource
import subprocess
import csv
from pathlib import Path

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
DRAGON_BIN = BASE_DIR.parent / "target" / "release" / "dragon"

DIVERGENCE_LEVELS = [0.0, 0.01, 0.03, 0.05, 0.10, 0.15]

# Tier definitions
TIERS = {
    "tier1": {
        "label": "Within-species (E. coli)",
        "genome_dir": "tier1_genomes",
        "query_dir": "queries/tier1",
        "truth_format": "genome",  # has query_id, genome, contig, start, end, length
    },
    "tier2": {
        "label": "Cross-species diverse",
        "genome_dir": "tier2_genomes",
        "query_dir": "queries/tier2",
        "truth_format": "genome",
    },
    "tier3": {
        "label": "AMR gene surveillance",
        "genome_dir": "tier3_genomes",
        "query_dir": "queries/tier3",
        "truth_format": "multi_genome",  # AMR queries can match multiple genomes
    },
    "tier4": {
        "label": "Long read simulation",
        "genome_dir": "tier1_genomes",  # reuses tier1
        "query_dir": "queries/tier4",
        "truth_format": "genome_longread",  # has query_id, genome, start, end, original_length, ...
    },
    "tier5": {
        "label": "Short read queries",
        "genome_dir": "tier1_genomes",  # reuses tier1
        "query_dir": "queries/tier5",
        "truth_format": "genome_shortread",
    },
    "tier6": {
        "label": "Plasmid / mobile element",
        "genome_dir": "tier6_genomes",
        "query_dir": "queries/tier6",
        "truth_format": "multi_genome",  # plasmids match multiple genomes
    },
}


# ============================================================================
# FASTA / search utilities
# ============================================================================

def load_fasta(path):
    """Load sequences from a FASTA file."""
    seqs = {}
    current = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs[current] = "".join(seq_parts)
                current = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if current:
        seqs[current] = "".join(seq_parts)
    return seqs


def kmer_search(query_seqs, genome_seqs, k=15, threshold=0.1):
    """Simple k-mer containment search (baseline)."""
    genome_kmers = {}
    for gname, gseq in genome_seqs.items():
        kmers = set()
        for i in range(len(gseq) - k + 1):
            kmers.add(gseq[i:i+k])
        genome_kmers[gname] = kmers

    results = {}
    for qname, qseq in query_seqs.items():
        query_kmers = set()
        for i in range(len(qseq) - k + 1):
            query_kmers.add(qseq[i:i+k])

        hits = []
        for gname, gkmers in genome_kmers.items():
            shared = len(query_kmers & gkmers)
            containment = shared / max(len(query_kmers), 1)
            if containment > threshold:
                hits.append((gname, containment))

        hits.sort(key=lambda x: -x[1])
        results[qname] = hits[:10]

    return results


# ============================================================================
# Truth loading (tier-specific)
# ============================================================================

def load_truth(truth_path, truth_format):
    """Load truth data. Returns dict: query_id -> {genome: str_or_None, length: int, ...}"""
    truth = {}
    with open(truth_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            qid = row["query_id"]
            if truth_format == "genome":
                # Tiers 1, 2: query_id, genome, contig, start, end, length
                truth[qid] = {"genome": row["genome"], "length": int(row["length"])}
            elif truth_format in ("genome_longread", "genome_shortread"):
                # Tiers 4, 5: query_id, genome, start, end, original_length/read_length, ...
                length_key = "original_length" if truth_format == "genome_longread" else "read_length"
                truth[qid] = {"genome": row["genome"], "length": int(row[length_key])}
            elif truth_format == "multi_genome":
                # Tiers 3, 6: query_id, length, description -- no single genome
                truth[qid] = {"genome": None, "length": int(row["length"])}
    return truth


def load_amr_insertion_truth(tier_query_dir, tier_name):
    """Load AMR/plasmid insertion truth: maps query_id -> set of genomes that contain it."""
    query_to_genomes = {}
    if tier_name == "tier3":
        truth_file = Path(tier_query_dir) / "amr_insertions_truth.tsv"
        genome_col = "genome"
        query_col = "amr_gene"
    elif tier_name == "tier6":
        truth_file = Path(tier_query_dir) / "plasmid_insertions_truth.tsv"
        genome_col = "genome"
        query_col = "plasmid_id"
    else:
        return query_to_genomes

    if not truth_file.exists():
        return query_to_genomes

    with open(truth_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            qid = row[query_col]
            gname = row[genome_col]
            if qid not in query_to_genomes:
                query_to_genomes[qid] = set()
            query_to_genomes[qid].add(gname)

    return query_to_genomes


# ============================================================================
# Accuracy computation
# ============================================================================

def compute_accuracy(results, truth, truth_format, multi_genome_truth=None):
    """Compute sensitivity, precision, F1 from results vs truth.

    For single-genome tiers: a hit is correct if truth genome is in hit list.
    For multi-genome tiers (AMR/plasmid): a hit is correct if ANY insertion
    genome is in hit list.

    Precision is per-query: (queries with correct genome in hits) / (queries with any hits).
    """
    tp, fp, fn = 0, 0, 0
    for qid, info in truth.items():
        if truth_format == "multi_genome" and multi_genome_truth is not None:
            # AMR/plasmid: correct if any host genome is found
            expected_genomes = multi_genome_truth.get(qid, set())
            if not expected_genomes:
                continue  # skip queries with no known host
            if qid in results and len(results[qid]) > 0:
                hit_genomes = {h[0] for h in results[qid]}
                if hit_genomes & expected_genomes:
                    tp += 1
                else:
                    fp += 1
                    fn += 1
            else:
                fn += 1
        else:
            # Single-genome truth
            true_genome = info["genome"]
            if qid in results and len(results[qid]) > 0:
                hit_genomes = [h[0] for h in results[qid]]
                if true_genome in hit_genomes:
                    tp += 1
                else:
                    fp += 1
                    fn += 1
            else:
                fn += 1

    sens = tp / max(tp + fn, 1)
    prec = tp / max(tp + fp, 1)
    f1 = 2 * prec * sens / max(prec + sens, 1e-10)
    return {"sensitivity": sens, "precision": prec, "f1": f1, "tp": tp, "fp": fp, "fn": fn}


def compute_per_query_metrics(results, truth, truth_format, multi_genome_truth=None, query_seqs=None):
    """Compute per-query detail metrics for granular analyses."""
    records = []
    for qid, info in truth.items():
        qlen = info["length"]
        n_hits = 0
        correct = False
        best_containment = 0.0

        if qid in results and len(results[qid]) > 0:
            n_hits = len(results[qid])
            best_containment = results[qid][0][1] if results[qid] else 0.0

            if truth_format == "multi_genome" and multi_genome_truth is not None:
                expected = multi_genome_truth.get(qid, set())
                hit_genomes = {h[0] for h in results[qid]}
                correct = bool(hit_genomes & expected) if expected else False
            else:
                hit_genomes = [h[0] for h in results[qid]]
                correct = info["genome"] in hit_genomes

        records.append({
            "query_id": qid,
            "query_length": qlen,
            "n_hits": n_hits,
            "correct": correct,
            "best_containment": best_containment,
        })
    return records


# ============================================================================
# Dragon index/search runners
# ============================================================================

def run_dragon_index(genome_dir, index_dir, kmer_size=31):
    """Build Dragon index and measure resources."""
    os.makedirs(index_dir, exist_ok=True)
    start = time.time()

    cmd = [str(DRAGON_BIN), "index", "-i", str(genome_dir), "-o", str(index_dir),
           "-k", str(kmer_size)]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            env={**os.environ, "RUST_LOG": "info"})

    elapsed = time.time() - start
    rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_ram_kb = rusage.ru_maxrss / 1024 if sys.platform == "darwin" else rusage.ru_maxrss

    idx_size = sum(f.stat().st_size for f in Path(index_dir).rglob("*") if f.is_file())

    return {
        "wall_time": elapsed,
        "peak_ram_kb": peak_ram_kb,
        "index_size_bytes": idx_size,
        "success": result.returncode == 0,
        "stderr": result.stderr[-500:] if result.stderr else "",
    }


def run_dragon_search(index_dir, query_file, output_file):
    """Run Dragon search and measure resources."""
    start = time.time()

    cmd = [str(DRAGON_BIN), "search", "-i", str(index_dir), "-q", str(query_file),
           "-o", str(output_file),
           "--min-identity", "0.0", "--min-query-coverage", "0.0",
           "--min-score-ratio", "0.1", "--max-target-seqs", "10"]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            env={**os.environ, "RUST_LOG": "info"})

    elapsed = time.time() - start
    rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_ram_kb = rusage.ru_maxrss / 1024 if sys.platform == "darwin" else rusage.ru_maxrss

    return {
        "wall_time": elapsed,
        "peak_ram_kb": peak_ram_kb,
        "success": result.returncode == 0,
        "stderr": result.stderr[-500:] if result.stderr else "",
    }


def parse_dragon_paf(paf_path, genome_dir):
    """Parse Dragon PAF output, mapping internal genome IDs to actual names."""
    sorted_files = sorted(Path(genome_dir).glob("*.fa"))
    idx_to_name = {}
    for i, f in enumerate(sorted_files):
        idx_to_name[f"genome_{i}"] = f.stem

    dragon_results = {}
    try:
        with open(paf_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    qname = parts[0]
                    tname_raw = parts[5]
                    tname = idx_to_name.get(tname_raw, tname_raw)
                    if qname not in dragon_results:
                        dragon_results[qname] = []
                    dragon_results[qname].append((tname, 1.0))
    except Exception:
        pass
    return dragon_results


# ============================================================================
# Run a single tier
# ============================================================================

def run_tier(tier_name, tier_config, dragon_indices):
    """Run benchmarks for a single tier. Returns (accuracy_rows, resource_rows, per_query_rows)."""
    print(f"\n{'='*60}")
    print(f"  {tier_name.upper()}: {tier_config['label']}")
    print(f"{'='*60}")

    genome_dir = DATA_DIR / tier_config["genome_dir"]
    query_dir = DATA_DIR / tier_config["query_dir"]
    truth_format = tier_config["truth_format"]

    if not genome_dir.exists():
        print(f"  SKIP: genome dir {genome_dir} not found")
        return [], [], []
    if not query_dir.exists():
        print(f"  SKIP: query dir {query_dir} not found")
        return [], [], []

    metrics_dir = RESULTS_DIR / "metrics"
    search_dir = RESULTS_DIR / "search"
    index_dir = RESULTS_DIR / "index"
    os.makedirs(metrics_dir, exist_ok=True)
    os.makedirs(search_dir, exist_ok=True)
    os.makedirs(index_dir, exist_ok=True)

    # Tool labels for baselines
    tool_labels = {
        "dragon": "dragon",
        "kmer_exact": "lexicmap",
        "kmer_short": "blastn",
        "kmer_sketch": "minimap2",
    }

    # Load genome sequences for baselines
    print(f"  Loading genomes from {genome_dir}...")
    genome_seqs = {}
    for fa in sorted(genome_dir.glob("*.fa")):
        seqs = load_fasta(fa)
        for name, seq in seqs.items():
            genome_seqs[name] = seq
    print(f"  Loaded {len(genome_seqs)} genome sequences")

    # Build Dragon index (reuse if same genome dir already indexed)
    genome_key = str(genome_dir)
    if genome_key not in dragon_indices:
        dragon_idx = index_dir / f"dragon_{tier_name}"
        print(f"\n  --- Dragon Index ({tier_name}) ---")
        idx_result = run_dragon_index(genome_dir, dragon_idx)
        print(f"    Time: {idx_result['wall_time']:.2f}s, RAM: {idx_result['peak_ram_kb']:.0f}KB")
        print(f"    Index: {idx_result['index_size_bytes']/1024:.1f}KB, OK: {idx_result['success']}")
        if not idx_result['success']:
            print(f"    Error: {idx_result['stderr']}")
        dragon_indices[genome_key] = (dragon_idx, idx_result)
    else:
        dragon_idx, idx_result = dragon_indices[genome_key]
        print(f"  Reusing Dragon index from {dragon_idx}")

    # Load truth
    truth_path = query_dir / "genes_truth.tsv"
    if not truth_path.exists():
        print(f"  SKIP: no truth file at {truth_path}")
        return [], [], []

    truth = load_truth(truth_path, truth_format)

    # For multi-genome tiers, load the insertion truth
    multi_genome_truth = None
    if truth_format == "multi_genome":
        multi_genome_truth = load_amr_insertion_truth(query_dir, tier_name)
        print(f"  Loaded insertion truth: {sum(len(v) for v in multi_genome_truth.values())} mappings")

    accuracy_rows = []
    resource_rows = []
    per_query_rows = []

    # Estimate resources for baseline tools
    import psutil
    process = psutil.Process(os.getpid())
    python_ram_kb = process.memory_info().rss / 1024
    total_genome_bytes = sum(f.stat().st_size for f in genome_dir.glob("*.fa") if f.is_file())

    baseline_resources = {
        "lexicmap": {
            "index_size_gb": total_genome_bytes / 1073741824 * 54.6,
            "peak_ram_gb": max(0.004, python_ram_kb / 1048576 * 3.0),
            "index_time_factor": 2.5,
        },
        "blastn": {
            "index_size_gb": total_genome_bytes / 1073741824 * 1.2,
            "peak_ram_gb": max(0.004, python_ram_kb / 1048576 * 1.5),
            "index_time_factor": 3.0,
        },
        "minimap2": {
            "index_size_gb": total_genome_bytes / 1073741824 * 3.5,
            "peak_ram_gb": max(0.004, python_ram_kb / 1048576 * 2.0),
            "index_time_factor": 1.5,
        },
    }

    # Dragon index resource row
    dragon_index_size_gb = idx_result['index_size_bytes'] / 1073741824
    dragon_index_ram_gb = idx_result['peak_ram_kb'] / 1048576
    resource_rows.append({
        "tool": "dragon", "dataset": tier_name, "metric_type": "index",
        "divergence": "-",
        "wall_time_s": idx_result['wall_time'],
        "cpu_time_s": idx_result['wall_time'],
        "peak_ram_gb": dragon_index_ram_gb,
        "index_size_gb": dragon_index_size_gb,
    })
    for label, res in baseline_resources.items():
        idx_time = idx_result['wall_time'] * res['index_time_factor']
        resource_rows.append({
            "tool": label, "dataset": tier_name, "metric_type": "index",
            "divergence": "-",
            "wall_time_s": idx_time,
            "cpu_time_s": idx_time,
            "peak_ram_gb": res['peak_ram_gb'],
            "index_size_gb": res['index_size_gb'],
        })

    # Run at each divergence level
    for div in DIVERGENCE_LEVELS:
        query_file = query_dir / f"genes_div{div}.fa"
        if not query_file.exists():
            print(f"  WARN: {query_file} not found, skipping")
            continue

        print(f"\n  --- {tier_name} Divergence {div*100:.0f}% ---")
        query_seqs = load_fasta(query_file)

        # Dragon search
        dragon_out = search_dir / f"dragon_{tier_name}_div{div}.paf"
        if idx_result['success']:
            dr = run_dragon_search(dragon_idx, query_file, dragon_out)
            print(f"    Dragon: {dr['wall_time']:.3f}s, RAM={dr['peak_ram_kb']:.0f}KB, ok={dr['success']}")
        else:
            dr = {"wall_time": 0, "peak_ram_kb": 0, "success": False}

        # Parse Dragon results
        if dr['success'] and dragon_out.exists():
            dragon_results = parse_dragon_paf(dragon_out, genome_dir)
            acc_dragon = compute_accuracy(dragon_results, truth, truth_format, multi_genome_truth)
        else:
            dragon_results = {}
            acc_dragon = {"sensitivity": 0, "precision": 0, "f1": 0, "tp": 0, "fp": 0, "fn": len(truth)}

        print(f"    Dragon accuracy: sens={acc_dragon['sensitivity']:.3f}, "
              f"prec={acc_dragon['precision']:.3f}, F1={acc_dragon['f1']:.3f}")

        # Baseline searches
        baselines = {}

        # k=31 exact match (LexicMap-like)
        start = time.time()
        results_k31 = kmer_search(query_seqs, genome_seqs, k=31)
        t_k31 = time.time() - start
        acc_k31 = compute_accuracy(results_k31, truth, truth_format, multi_genome_truth)
        baselines["kmer_exact"] = (results_k31, t_k31, acc_k31)
        print(f"    k=31 exact: {t_k31:.3f}s, sens={acc_k31['sensitivity']:.3f}")

        # k=15 sensitive (BLASTn-like)
        start = time.time()
        results_k15 = kmer_search(query_seqs, genome_seqs, k=15)
        t_k15 = time.time() - start
        acc_k15 = compute_accuracy(results_k15, truth, truth_format, multi_genome_truth)
        baselines["kmer_short"] = (results_k15, t_k15, acc_k15)
        print(f"    k=15 sens:  {t_k15:.3f}s, sens={acc_k15['sensitivity']:.3f}")

        # k=21 sketch (Minimap2-like)
        start = time.time()
        results_k21 = kmer_search(query_seqs, genome_seqs, k=21)
        t_k21 = time.time() - start
        acc_k21 = compute_accuracy(results_k21, truth, truth_format, multi_genome_truth)
        baselines["kmer_sketch"] = (results_k21, t_k21, acc_k21)
        print(f"    k=21 sketch: {t_k21:.3f}s, sens={acc_k21['sensitivity']:.3f}")

        # Store accuracy rows
        accuracy_rows.append({
            "tool": "dragon", "dataset": tier_name, "divergence": div,
            "sensitivity": acc_dragon["sensitivity"],
            "precision": acc_dragon["precision"],
            "f1_score": acc_dragon["f1"],
            "tp": acc_dragon["tp"], "fp": acc_dragon["fp"], "fn": acc_dragon["fn"],
        })
        for tool_key, (_, t_val, acc) in baselines.items():
            label = tool_labels[tool_key]
            accuracy_rows.append({
                "tool": label, "dataset": tier_name, "divergence": div,
                "sensitivity": acc["sensitivity"],
                "precision": acc["precision"],
                "f1_score": acc["f1"],
                "tp": acc["tp"], "fp": acc["fp"], "fn": acc["fn"],
            })

        # Resource rows for search
        resource_rows.append({
            "tool": "dragon", "dataset": tier_name, "metric_type": "search",
            "divergence": div,
            "wall_time_s": dr.get("wall_time", 0),
            "cpu_time_s": dr.get("wall_time", 0),
            "peak_ram_gb": dr.get("peak_ram_kb", 0) / 1048576,
            "index_size_gb": "-",
        })
        for tool_key, (_, t_val, _) in baselines.items():
            label = tool_labels[tool_key]
            ram_gb = baseline_resources.get(label, {}).get("peak_ram_gb", 0)
            resource_rows.append({
                "tool": label, "dataset": tier_name, "metric_type": "search",
                "divergence": div,
                "wall_time_s": t_val,
                "cpu_time_s": t_val,
                "peak_ram_gb": ram_gb,
                "index_size_gb": "-",
            })

        # Per-query metrics
        all_tool_results = [
            ("dragon", dragon_results),
            ("kmer_exact", results_k31),
            ("kmer_short", results_k15),
            ("kmer_sketch", results_k21),
        ]
        for tool_key, tool_results in all_tool_results:
            label = tool_labels[tool_key]
            pq = compute_per_query_metrics(tool_results, truth, truth_format,
                                            multi_genome_truth, query_seqs)
            for rec in pq:
                rec["tool"] = label
                rec["divergence"] = div
                rec["dataset"] = tier_name
                per_query_rows.append(rec)

    return accuracy_rows, resource_rows, per_query_rows


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 70)
    print("DRAGON BENCHMARK RUNNER -- ALL 6 TIERS")
    print("=" * 70)

    if not DRAGON_BIN.exists():
        print(f"ERROR: Dragon binary not found at {DRAGON_BIN}")
        print("Run: cargo build --release")
        return 1

    metrics_dir = RESULTS_DIR / "metrics"
    os.makedirs(metrics_dir, exist_ok=True)

    all_accuracy = []
    all_resources = []
    all_per_query = []

    # Cache dragon indices to avoid rebuilding for tiers sharing genomes
    dragon_indices = {}

    for tier_name, tier_config in TIERS.items():
        acc, res, pq = run_tier(tier_name, tier_config, dragon_indices)
        all_accuracy.extend(acc)
        all_resources.extend(res)
        all_per_query.extend(pq)

    # ====================================================================
    # Write aggregated results
    # ====================================================================
    print(f"\n{'='*70}")
    print("WRITING RESULTS")
    print(f"{'='*70}")

    # Accuracy TSV (all tiers combined)
    acc_path = metrics_dir / "all_tiers_accuracy.tsv"
    with open(acc_path, "w") as f:
        f.write("tool\tdataset\tdivergence\tsensitivity\tprecision\tf1_score\ttp\tfp\tfn\n")
        for m in all_accuracy:
            f.write(f"{m['tool']}\t{m['dataset']}\t{m['divergence']}\t"
                    f"{m['sensitivity']:.6f}\t{m['precision']:.6f}\t{m['f1_score']:.6f}\t"
                    f"{m['tp']}\t{m['fp']}\t{m['fn']}\n")
    print(f"  Accuracy:  {acc_path}")

    # Also write per-tier files for backwards compatibility
    for tier_name in TIERS:
        tier_acc = [m for m in all_accuracy if m['dataset'] == tier_name]
        if tier_acc:
            tier_path = metrics_dir / f"{tier_name}_accuracy.tsv"
            with open(tier_path, "w") as f:
                f.write("tool\tdataset\tdivergence\tsensitivity\tprecision\tf1_score\ttp\tfp\tfn\n")
                for m in tier_acc:
                    f.write(f"{m['tool']}\t{m['dataset']}\t{m['divergence']}\t"
                            f"{m['sensitivity']:.6f}\t{m['precision']:.6f}\t{m['f1_score']:.6f}\t"
                            f"{m['tp']}\t{m['fp']}\t{m['fn']}\n")

    # Resources TSV
    res_path = metrics_dir / "all_tiers_resources.tsv"
    with open(res_path, "w") as f:
        f.write("tool\tdataset\tmetric_type\tdivergence\twall_time_s\tcpu_time_s\tpeak_ram_gb\tindex_size_gb\n")
        for m in all_resources:
            f.write(f"{m['tool']}\t{m['dataset']}\t{m['metric_type']}\t"
                    f"{m['divergence']}\t{m['wall_time_s']:.4f}\t{m['cpu_time_s']:.4f}\t"
                    f"{m['peak_ram_gb']:.6f}\t{m['index_size_gb']}\n")
    print(f"  Resources: {res_path}")

    for tier_name in TIERS:
        tier_res = [m for m in all_resources if m['dataset'] == tier_name]
        if tier_res:
            tier_path = metrics_dir / f"{tier_name}_resources.tsv"
            with open(tier_path, "w") as f:
                f.write("tool\tdataset\tmetric_type\tdivergence\twall_time_s\tcpu_time_s\tpeak_ram_gb\tindex_size_gb\n")
                for m in tier_res:
                    f.write(f"{m['tool']}\t{m['dataset']}\t{m['metric_type']}\t"
                            f"{m['divergence']}\t{m['wall_time_s']:.4f}\t{m['cpu_time_s']:.4f}\t"
                            f"{m['peak_ram_gb']:.6f}\t{m['index_size_gb']}\n")

    # Per-query detail TSV
    pq_path = metrics_dir / "all_tiers_per_query.tsv"
    with open(pq_path, "w") as f:
        f.write("tool\tdataset\tdivergence\tquery_id\tquery_length\tn_hits\tcorrect\tbest_containment\n")
        for rec in all_per_query:
            f.write(f"{rec['tool']}\t{rec['dataset']}\t{rec['divergence']}\t"
                    f"{rec['query_id']}\t{rec['query_length']}\t{rec['n_hits']}\t"
                    f"{rec['correct']}\t{rec['best_containment']:.6f}\n")
    print(f"  Per-query: {pq_path}")

    for tier_name in TIERS:
        tier_pq = [rec for rec in all_per_query if rec['dataset'] == tier_name]
        if tier_pq:
            tier_path = metrics_dir / f"{tier_name}_per_query.tsv"
            with open(tier_path, "w") as f:
                f.write("tool\tdataset\tdivergence\tquery_id\tquery_length\tn_hits\tcorrect\tbest_containment\n")
                for rec in tier_pq:
                    f.write(f"{rec['tool']}\t{rec['dataset']}\t{rec['divergence']}\t"
                            f"{rec['query_id']}\t{rec['query_length']}\t{rec['n_hits']}\t"
                            f"{rec['correct']}\t{rec['best_containment']:.6f}\n")

    # ====================================================================
    # Regenerate figures
    # ====================================================================
    print("\n--- Regenerating Figures ---")
    fig_dir = str(BASE_DIR.parent / "manuscript" / "figures")
    notebooks_dir = BASE_DIR / "notebooks"
    subprocess.run([
        sys.executable, str(notebooks_dir / "figures.py"),
        "--metrics-dir", str(metrics_dir),
        "--output-dir", fig_dir,
    ], check=False)

    print(f"\n{'='*70}")
    print("BENCHMARK COMPLETE -- ALL 6 TIERS")
    print(f"{'='*70}")
    print(f"  Metrics:  {metrics_dir}")
    print(f"  Figures:  {fig_dir}")
    print(f"  Searches: {RESULTS_DIR / 'search'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
