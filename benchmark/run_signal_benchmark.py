#!/usr/bin/env python3
"""
Signal-level benchmark for Dragon.

Generates simulated nanopore signal data by:
1. Building a signal index (Dragon generates expected signals)
2. Having Dragon export expected signal for test genomes
3. Adding Gaussian noise in Python to simulate real reads
4. Running Dragon signal-search and measuring accuracy

Also produces manuscript figures (15-16) for the signal module.
"""

import os
import sys
import random
import time
import struct
import json
import subprocess
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

random.seed(42)
np.random.seed(42)

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
DRAGON_BIN = BASE_DIR.parent / "target" / "release" / "dragon"


def load_fasta(path):
    """Load sequences from a FASTA file."""
    seqs = {}
    current = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs[current] = "".join(parts)
                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if current:
        seqs[current] = "".join(parts)
    return seqs


def load_pore_model_rust():
    """Replicate the exact Rust pore model for expected signal generation."""
    k = 5
    base_contrib = [-5.0, 3.0, 8.0, -2.0]  # A, C, G, T
    nn_interact = [
        [-1.2, 0.8, 2.1, -0.5],   # A->
        [0.5, -0.9, 1.5, 0.3],    # C->
        [1.8, 0.2, -1.1, 1.0],    # G->
        [-0.3, 0.7, 1.3, -1.5],   # T->
    ]
    pos_weights = [0.6, 0.8, 1.0, 0.8, 0.6]
    baseline = 90.0

    base_map = {"A": 0, "C": 1, "G": 2, "T": 3}

    def kmer_to_pa(kmer):
        level = baseline
        for pos, base in enumerate(kmer):
            bi = base_map.get(base, 0)
            level += base_contrib[bi] * pos_weights[pos]
        for i in range(len(kmer) - 1):
            b1 = base_map.get(kmer[i], 0)
            b2 = base_map.get(kmer[i+1], 0)
            level += nn_interact[b1][b2]
        return level

    return kmer_to_pa


def sequence_to_signal(seq, kmer_to_pa, noise_std=0.0):
    """Convert DNA sequence to expected pA signal, optionally adding noise."""
    signal = []
    for i in range(len(seq) - 4):
        kmer = seq[i:i+5]
        if all(b in "ACGT" for b in kmer):
            pa = kmer_to_pa(kmer)
            if noise_std > 0:
                pa += np.random.normal(0, noise_std)
            signal.append(pa)
    return signal


def main():
    print("=" * 60)
    print("SIGNAL-LEVEL BENCHMARK")
    print("=" * 60)

    if not DRAGON_BIN.exists():
        print(f"ERROR: Dragon binary not found at {DRAGON_BIN}")
        return 1

    # Use tier2 genomes (cross-species) for clear species-level discrimination
    genome_dir = DATA_DIR / "tier2_genomes"
    if not genome_dir.exists():
        genome_dir = DATA_DIR / "tier1_genomes"
    if not genome_dir.exists():
        print("ERROR: Test genomes not found. Run setup_test_data.py first.")
        return 1

    signal_dir = RESULTS_DIR / "signal"
    metrics_dir = RESULTS_DIR / "metrics"
    fig_dir = BASE_DIR.parent / "manuscript" / "figures"
    os.makedirs(signal_dir, exist_ok=True)
    os.makedirs(metrics_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # Load genomes
    print("\nLoading genomes...")
    genome_seqs = {}
    for fa in sorted(genome_dir.glob("*.fa")):
        seqs = load_fasta(fa)
        for name, seq in seqs.items():
            genome_seqs[name] = seq
    print(f"  Loaded {len(genome_seqs)} genomes from {genome_dir.name}")

    # Load pore model (exact match to Rust)
    kmer_to_pa = load_pore_model_rust()

    # ================================================================
    # 1. Build signal index
    # ================================================================
    print("\n--- Building Signal Index ---")
    signal_index_dir = signal_dir / "signal_index"
    start = time.time()
    result = subprocess.run(
        [str(DRAGON_BIN), "signal-index",
         "-i", str(genome_dir),
         "-o", str(signal_index_dir),
         "--num-levels", "16"],
        capture_output=True, text=True
    )
    index_time = time.time() - start
    print(f"  Time: {index_time:.2f}s, OK: {result.returncode == 0}")
    if result.returncode != 0:
        print(f"  Error: {result.stderr[-300:]}")
        return 1

    # ================================================================
    # 2. Generate signal reads at different noise levels
    # ================================================================
    noise_levels = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0]
    num_reads = 60
    read_signal_lengths = [200, 500, 1000, 2000]  # in signal samples

    genome_names = list(genome_seqs.keys())

    print("\n--- Generating Simulated Signal Reads ---")

    truth = {}  # read_id -> source genome
    all_results = {}  # noise_std -> list of result dicts

    for noise_std in noise_levels:
        signal_file = signal_dir / f"reads_noise{noise_std:.1f}.tsv"
        reads = []

        with open(signal_file, "w") as f:
            for i in range(num_reads):
                gname = random.choice(genome_names)
                gseq = genome_seqs[gname]
                rlen_signal = random.choice(read_signal_lengths)
                # Need rlen_signal + 4 bases for signal of this length
                rlen_bases = rlen_signal + 4
                if rlen_bases > len(gseq) - 10:
                    rlen_bases = len(gseq) - 10
                start_pos = random.randint(0, len(gseq) - rlen_bases)
                subseq = gseq[start_pos:start_pos + rlen_bases]

                signal = sequence_to_signal(subseq, kmer_to_pa, noise_std=noise_std)
                if len(signal) < 15:
                    continue

                read_id = f"signal_read_{i:04d}"
                truth[read_id] = gname
                reads.append(read_id)

                signal_str = "\t".join(f"{v:.2f}" for v in signal)
                f.write(f"{read_id}\t{signal_str}\n")

        print(f"  Noise σ={noise_std:.1f}: {len(reads)} reads -> {signal_file.name}")

        # Run signal search
        output_file = signal_dir / f"signal_results_noise{noise_std:.1f}.tsv"
        start = time.time()
        result = subprocess.run(
            [str(DRAGON_BIN), "signal-search",
             "-i", str(signal_index_dir),
             "-q", str(signal_file),
             "-o", str(output_file),
             "--signal-kmer-size", "6",
             "--min-hits", "2"],
            capture_output=True, text=True
        )
        search_time = time.time() - start
        print(f"    Search: {search_time:.3f}s, OK: {result.returncode == 0}")

        if result.returncode != 0:
            print(f"    Error: {result.stderr[-200:]}")
            all_results[noise_std] = []
            continue

        # Parse results
        hits = {}
        if output_file.exists():
            with open(output_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        rid = parts[0]
                        gname_hit = parts[2]
                        score = float(parts[5])
                        if rid not in hits:
                            hits[rid] = []
                        hits[rid].append((gname_hit, score))

        # Species-level accuracy
        def get_species(gname):
            if gname.startswith("species"):
                return gname.split("_genome_")[0]
            if gname.startswith("ecoli"):
                return "ecoli"
            return gname

        read_results = []
        for rid in reads:
            true_genome = truth.get(rid)
            if rid in hits and len(hits[rid]) > 0:
                sorted_hits = sorted(hits[rid], key=lambda x: -x[1])
                top_genome = sorted_hits[0][0]
                # Species-level: correct if same species cluster
                is_species_correct = (get_species(top_genome) == get_species(true_genome))
                # Top-5 species check
                top5_species = {get_species(h[0]) for h in sorted_hits[:5]}
                is_top5_correct = get_species(true_genome) in top5_species
                # Any hit correct
                all_hit_species = {get_species(h[0]) for h in sorted_hits}
                is_any_correct = get_species(true_genome) in all_hit_species

                read_results.append({
                    "read_id": rid,
                    "species_correct": is_species_correct,
                    "top5_correct": is_top5_correct,
                    "any_correct": is_any_correct,
                    "exact_correct": (top_genome == true_genome),
                    "score": sorted_hits[0][1],
                    "n_hits": len(sorted_hits),
                    "has_hit": True,
                })
            else:
                read_results.append({
                    "read_id": rid,
                    "species_correct": False,
                    "top5_correct": False,
                    "any_correct": False,
                    "exact_correct": False,
                    "score": 0.0,
                    "n_hits": 0,
                    "has_hit": False,
                })

        all_results[noise_std] = read_results

        tp_species = sum(1 for r in read_results if r["species_correct"])
        tp_top5 = sum(1 for r in read_results if r["top5_correct"])
        tp_any = sum(1 for r in read_results if r["any_correct"])
        found = sum(1 for r in read_results if r["has_hit"])
        total = len(read_results)
        print(f"    Found: {found}/{total}, Species top-1: {tp_species}/{total} "
              f"({tp_species/max(total,1):.1%}), top-5: {tp_top5}/{total} "
              f"({tp_top5/max(total,1):.1%}), any: {tp_any}/{total} "
              f"({tp_any/max(total,1):.1%})")

    # ================================================================
    # 3. Generate Figures
    # ================================================================
    print("\n--- Generating Signal Figures ---")

    # ---- Figure 15: Signal discretization concept ----
    fig, axes = plt.subplots(3, 1, figsize=(12, 8),
                              gridspec_kw={"height_ratios": [1.5, 1.5, 1]})

    example_seq = list(genome_seqs.values())[0][1000:1100]
    raw_signal = sequence_to_signal(example_seq, kmer_to_pa, noise_std=0.0)
    noisy_signal = sequence_to_signal(example_seq, kmer_to_pa, noise_std=3.0)

    # (a) Clean expected signal
    ax = axes[0]
    ax.plot(raw_signal, color="#2A9D8F", linewidth=1.2, alpha=0.9)
    ax.set_ylabel("Current (pA)")
    ax.set_title("(a) Expected signal from pore model (R10.4.1 approximate)")
    ax.set_xlim(0, len(raw_signal))
    ax.grid(True, alpha=0.3)

    # (b) Noisy simulated read
    ax = axes[1]
    ax.plot(noisy_signal, color="#E63946", linewidth=0.8, alpha=0.7)
    ax.set_ylabel("Current (pA)")
    ax.set_title("(b) Simulated nanopore signal (σ=3.0 pA noise)")
    ax.set_xlim(0, len(noisy_signal))
    ax.grid(True, alpha=0.3)

    # (c) Discretized
    ax = axes[2]
    signal_np = np.array(noisy_signal)
    med = np.median(signal_np)
    mad = np.median(np.abs(signal_np - med))
    if mad > 0:
        normalized = (signal_np - med) / (mad * 1.4826)
    else:
        normalized = signal_np - med
    num_levels = 16
    clipped = np.clip(normalized, -4.0, 4.0)
    bin_width = 8.0 / num_levels
    discretized = np.floor((clipped + 4.0) / bin_width).astype(int)
    discretized = np.clip(discretized, 0, num_levels - 1)

    colors = plt.cm.viridis(discretized / (num_levels - 1))
    ax.bar(range(len(discretized)), discretized, width=1.0, color=colors, edgecolor="none")
    ax.set_ylabel("Level (0-15)")
    ax.set_xlabel("Signal position")
    ax.set_title("(c) Discretized alphabet (16 levels) for FM-index backward search")
    ax.set_xlim(0, len(discretized))
    ax.set_ylim(-0.5, num_levels + 0.5)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = str(fig_dir / "fig15_signal_discretization.pdf")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")

    # ---- Figure 16: Signal search accuracy vs noise ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # (a) Sensitivity vs noise
    ax = axes[0]
    species_sens = []
    top5_sens = []
    any_sens = []
    for noise_std in noise_levels:
        results = all_results.get(noise_std, [])
        total = max(len(results), 1)
        species_sens.append(sum(1 for r in results if r.get("species_correct")) / total)
        top5_sens.append(sum(1 for r in results if r.get("top5_correct")) / total)
        any_sens.append(sum(1 for r in results if r.get("any_correct")) / total)

    ax.plot(noise_levels, species_sens, "D-", color="#E63946", linewidth=2.5,
            markersize=10, label="Top-1 species match")
    ax.plot(noise_levels, top5_sens, "s-", color="#2A9D8F", linewidth=2,
            markersize=8, label="Top-5 species match")
    ax.plot(noise_levels, any_sens, "o--", color="#457B9D", linewidth=1.5,
            markersize=7, label="Any hit correct species")
    ax.set_xlabel("Signal Noise (σ pA)")
    ax.set_ylabel("Sensitivity")
    ax.set_ylim(0, 1.05)
    ax.set_title("(a) Signal Search Sensitivity vs Noise Level")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) Score distributions
    ax = axes[1]
    for noise_std, color, label in [
        (0.0, "#2A9D8F", "Clean (σ=0)"),
        (2.0, "#E9C46A", "Low noise (σ=2)"),
        (5.0, "#F4A261", "Moderate (σ=5)"),
        (8.0, "#E63946", "High (σ=8)"),
    ]:
        results = all_results.get(noise_std, [])
        scores = [r["score"] for r in results if r["score"] > 0]
        if scores:
            ax.hist(scores, bins=20, alpha=0.4, color=color, label=label,
                    edgecolor="black", linewidth=0.5, density=True)

    ax.set_xlabel("Match Score")
    ax.set_ylabel("Density")
    ax.set_title("(b) Score Distribution by Noise Level")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = str(fig_dir / "fig16_signal_search_accuracy.pdf")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")

    # Write signal metrics
    signal_metrics_path = metrics_dir / "signal_accuracy.tsv"
    with open(signal_metrics_path, "w") as f:
        f.write("noise_std\tspecies_top1\tspecies_top5\tany_correct\tnum_reads\n")
        for i, noise_std in enumerate(noise_levels):
            results = all_results.get(noise_std, [])
            total = max(len(results), 1)
            f.write(f"{noise_std}\t{species_sens[i]:.6f}\t{top5_sens[i]:.6f}\t"
                    f"{any_sens[i]:.6f}\t{len(results)}\n")
    print(f"  Signal metrics: {signal_metrics_path}")

    print("\n" + "=" * 60)
    print("SIGNAL BENCHMARK COMPLETE")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
