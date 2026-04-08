#!/usr/bin/env python3
"""
Generate benchmark data from REAL genomes downloaded by download_genomes.py.

Uses existing tools (simulate/extract_genes.py, simulate/mutate_sequences.py)
to extract gene-length subsequences from real genomes and mutate them at
controlled divergence levels.

Produces the same 6-tier structure as setup_test_data.py but with real biology.
"""

import os
import random
import subprocess
import sys
from pathlib import Path

random.seed(42)

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
REAL_GENOMES_DIR = DATA_DIR / "real_genomes"
SIMULATE_DIR = BASE_DIR / "simulate"

DIVERGENCE_LEVELS = [0.0, 0.01, 0.03, 0.05, 0.10, 0.15]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, **kw):
    """Run a command, printing it first."""
    print(f"  -> {' '.join(str(c) for c in cmd)}")
    subprocess.check_call([str(c) for c in cmd], **kw)


def read_fasta(path):
    """Simple FASTA reader returning list of (name, seq)."""
    records = []
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(seq)))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name:
        records.append((name, "".join(seq)))
    return records


def write_fasta(path, records):
    """Write list of (name, seq) to FASTA."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def extract_random_genes(genome_dir, output_fa, truth_tsv, num_genes, min_len, max_len, seed=42):
    """Extract random subsequences from genomes in a directory."""
    run([
        sys.executable, SIMULATE_DIR / "extract_genes.py",
        "--genome-dir", genome_dir,
        "--output", output_fa,
        "--truth", truth_tsv,
        "--num-genes", str(num_genes),
        "--min-length", str(min_len),
        "--max-length", str(max_len),
        "--seed", str(seed),
    ])


def mutate_queries(input_fa, output_fa, truth_tsv, divergence, seed=42):
    """Mutate sequences at a given divergence level."""
    run([
        sys.executable, SIMULATE_DIR / "mutate_sequences.py",
        "--input", input_fa,
        "--output", output_fa,
        "--truth", truth_tsv,
        "--divergence", str(divergence),
        "--seed", str(seed),
    ])


# ---------------------------------------------------------------------------
# Tier builders
# ---------------------------------------------------------------------------

def setup_tier1(ecoli_dir):
    """Tier 1: Within-species (real E. coli genomes)."""
    tier_dir = DATA_DIR / "tier1_within_species"
    genome_dir = tier_dir / "genomes"
    os.makedirs(genome_dir, exist_ok=True)

    # Symlink real genomes
    fasta_files = sorted(ecoli_dir.glob("*.fna")) + sorted(ecoli_dir.glob("*.fa")) + sorted(ecoli_dir.glob("*.fasta"))
    if not fasta_files:
        print(f"WARNING: No FASTA files in {ecoli_dir}. Run download_genomes.py first.")
        return
    for f in fasta_files[:50]:
        dest = genome_dir / f.name
        if not dest.exists():
            os.symlink(f.resolve(), dest)

    # Extract 100 gene queries (500-2000 bp)
    base_fa = tier_dir / "queries_base.fa"
    base_truth = tier_dir / "queries_base_truth.tsv"
    extract_random_genes(genome_dir, base_fa, base_truth, 100, 500, 2000)

    # Mutate at each divergence level
    for div in DIVERGENCE_LEVELS:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = tier_dir / f"queries_div{div_str}.fa"
        out_truth = tier_dir / f"queries_div{div_str}_truth.tsv"
        if div == 0.0:
            # Just copy the base (no mutation)
            import shutil
            shutil.copy2(base_fa, out_fa)
            shutil.copy2(base_truth, out_truth)
        else:
            mutate_queries(base_fa, out_fa, out_truth, div, seed=42 + int(div * 100))

    print(f"Tier 1: {len(fasta_files[:50])} E. coli genomes, 100 queries x {len(DIVERGENCE_LEVELS)} divergence levels")


def setup_tier2(cross_species_dir):
    """Tier 2: Cross-species (real mixed-species genomes)."""
    tier_dir = DATA_DIR / "tier2_cross_species"
    genome_dir = tier_dir / "genomes"
    os.makedirs(genome_dir, exist_ok=True)

    # Collect all FASTA files from species subdirectories
    fasta_files = []
    for species_dir in sorted(cross_species_dir.iterdir()):
        if species_dir.is_dir():
            for f in sorted(species_dir.glob("*.fna")) + sorted(species_dir.glob("*.fa")):
                fasta_files.append(f)

    if not fasta_files:
        # Try flat directory
        fasta_files = sorted(cross_species_dir.glob("*.fna")) + sorted(cross_species_dir.glob("*.fa"))

    if not fasta_files:
        print(f"WARNING: No FASTA files in {cross_species_dir}. Run download_genomes.py first.")
        return

    for f in fasta_files[:30]:
        dest = genome_dir / f.name
        if not dest.exists():
            os.symlink(f.resolve(), dest)

    base_fa = tier_dir / "queries_base.fa"
    base_truth = tier_dir / "queries_base_truth.tsv"
    extract_random_genes(genome_dir, base_fa, base_truth, 80, 500, 2000, seed=123)

    for div in DIVERGENCE_LEVELS:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = tier_dir / f"queries_div{div_str}.fa"
        out_truth = tier_dir / f"queries_div{div_str}_truth.tsv"
        if div == 0.0:
            import shutil
            shutil.copy2(base_fa, out_fa)
            shutil.copy2(base_truth, out_truth)
        else:
            mutate_queries(base_fa, out_fa, out_truth, div, seed=123 + int(div * 100))

    print(f"Tier 2: {len(fasta_files[:30])} cross-species genomes, 80 queries x {len(DIVERGENCE_LEVELS)} divergence levels")


def setup_tier3(ecoli_dir):
    """Tier 3: AMR gene-like queries from real genomes (longer, diverse regions)."""
    tier_dir = DATA_DIR / "tier3_amr_genes"
    genome_dir = tier_dir / "genomes"
    os.makedirs(genome_dir, exist_ok=True)

    fasta_files = sorted(ecoli_dir.glob("*.fna")) + sorted(ecoli_dir.glob("*.fa"))
    for f in fasta_files[:30]:
        dest = genome_dir / f.name
        if not dest.exists():
            os.symlink(f.resolve(), dest)

    base_fa = tier_dir / "queries_base.fa"
    base_truth = tier_dir / "queries_base_truth.tsv"
    extract_random_genes(genome_dir, base_fa, base_truth, 20, 600, 3000, seed=200)

    for div in DIVERGENCE_LEVELS:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = tier_dir / f"queries_div{div_str}.fa"
        out_truth = tier_dir / f"queries_div{div_str}_truth.tsv"
        if div == 0.0:
            import shutil
            shutil.copy2(base_fa, out_fa)
            shutil.copy2(base_truth, out_truth)
        else:
            mutate_queries(base_fa, out_fa, out_truth, div, seed=200 + int(div * 100))

    print(f"Tier 3: 20 AMR-like queries x {len(DIVERGENCE_LEVELS)} divergence levels")


def setup_tier4(ecoli_dir):
    """Tier 4: Long read simulation (2-10 Kbp) from real genomes."""
    tier_dir = DATA_DIR / "tier4_long_reads"
    genome_dir = tier_dir / "genomes"
    os.makedirs(genome_dir, exist_ok=True)

    fasta_files = sorted(ecoli_dir.glob("*.fna")) + sorted(ecoli_dir.glob("*.fa"))
    for f in fasta_files[:50]:
        dest = genome_dir / f.name
        if not dest.exists():
            os.symlink(f.resolve(), dest)

    base_fa = tier_dir / "queries_base.fa"
    base_truth = tier_dir / "queries_base_truth.tsv"
    extract_random_genes(genome_dir, base_fa, base_truth, 40, 2000, 10000, seed=300)

    for div in DIVERGENCE_LEVELS:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = tier_dir / f"queries_div{div_str}.fa"
        out_truth = tier_dir / f"queries_div{div_str}_truth.tsv"
        if div == 0.0:
            import shutil
            shutil.copy2(base_fa, out_fa)
            shutil.copy2(base_truth, out_truth)
        else:
            # Long reads: higher insertion rate to simulate ONT error profile
            run([
                sys.executable, SIMULATE_DIR / "mutate_sequences.py",
                "--input", base_fa,
                "--output", out_fa,
                "--truth", out_truth,
                "--divergence", str(div),
                "--sub-rate", "0.5",
                "--ins-rate", "0.35",
                "--del-rate", "0.15",
                "--seed", str(300 + int(div * 100)),
            ])

    print(f"Tier 4: 40 long-read queries x {len(DIVERGENCE_LEVELS)} divergence levels")


def setup_tier5(ecoli_dir):
    """Tier 5: Short read queries (150/250 bp) from real genomes."""
    tier_dir = DATA_DIR / "tier5_short_reads"
    genome_dir = tier_dir / "genomes"
    os.makedirs(genome_dir, exist_ok=True)

    fasta_files = sorted(ecoli_dir.glob("*.fna")) + sorted(ecoli_dir.glob("*.fa"))
    for f in fasta_files[:50]:
        dest = genome_dir / f.name
        if not dest.exists():
            os.symlink(f.resolve(), dest)

    base_fa = tier_dir / "queries_base.fa"
    base_truth = tier_dir / "queries_base_truth.tsv"
    extract_random_genes(genome_dir, base_fa, base_truth, 200, 150, 250, seed=400)

    for div in DIVERGENCE_LEVELS:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = tier_dir / f"queries_div{div_str}.fa"
        out_truth = tier_dir / f"queries_div{div_str}_truth.tsv"
        if div == 0.0:
            import shutil
            shutil.copy2(base_fa, out_fa)
            shutil.copy2(base_truth, out_truth)
        else:
            mutate_queries(base_fa, out_fa, out_truth, div, seed=400 + int(div * 100))

    print(f"Tier 5: 200 short-read queries x {len(DIVERGENCE_LEVELS)} divergence levels")


def setup_tier6(ecoli_dir):
    """Tier 6: Plasmid-like mobile element queries (2-5 Kbp)."""
    tier_dir = DATA_DIR / "tier6_plasmids"
    genome_dir = tier_dir / "genomes"
    os.makedirs(genome_dir, exist_ok=True)

    fasta_files = sorted(ecoli_dir.glob("*.fna")) + sorted(ecoli_dir.glob("*.fa"))
    for f in fasta_files[:20]:
        dest = genome_dir / f.name
        if not dest.exists():
            os.symlink(f.resolve(), dest)

    base_fa = tier_dir / "queries_base.fa"
    base_truth = tier_dir / "queries_base_truth.tsv"
    extract_random_genes(genome_dir, base_fa, base_truth, 30, 2000, 5000, seed=500)

    for div in DIVERGENCE_LEVELS:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = tier_dir / f"queries_div{div_str}.fa"
        out_truth = tier_dir / f"queries_div{div_str}_truth.tsv"
        if div == 0.0:
            import shutil
            shutil.copy2(base_fa, out_fa)
            shutil.copy2(base_truth, out_truth)
        else:
            mutate_queries(base_fa, out_fa, out_truth, div, seed=500 + int(div * 100))

    print(f"Tier 6: 30 plasmid-like queries x {len(DIVERGENCE_LEVELS)} divergence levels")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Dragon Benchmark: Real Genome Data Setup")
    print("=" * 60)

    ecoli_dir = REAL_GENOMES_DIR / "ecoli"
    cross_species_dir = REAL_GENOMES_DIR / "cross_species"

    # Check that genomes exist
    has_ecoli = ecoli_dir.exists() and any(ecoli_dir.glob("*.fn*"))
    has_cross = cross_species_dir.exists()

    if not has_ecoli:
        print(f"\nERROR: No real genomes found at {ecoli_dir}")
        print("Run: python3 benchmark/download_genomes.py")
        print("This downloads real E. coli and cross-species genomes from NCBI.\n")
        sys.exit(1)

    print(f"\nUsing real genomes from: {REAL_GENOMES_DIR}")
    print()

    setup_tier1(ecoli_dir)
    setup_tier2(cross_species_dir if has_cross else ecoli_dir)
    setup_tier3(ecoli_dir)
    setup_tier4(ecoli_dir)
    setup_tier5(ecoli_dir)
    setup_tier6(ecoli_dir)

    print()
    print("=" * 60)
    print("Setup complete! Real genome benchmark data ready.")
    print(f"Data directory: {DATA_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
