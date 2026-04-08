#!/usr/bin/env python3
"""
Generate synthetic test genomes and query data for the Dragon benchmark.

Creates 6 diverse dataset tiers for robust benchmarking:

  Tier 1 -- Within-species (E. coli-like)
             50 genomes (~50 Kbp), 95-99% ANI, 100 gene queries
  Tier 2 -- Cross-species diverse
             30 genomes from 6 "species clusters", 80 gene queries
  Tier 3 -- AMR gene surveillance
             30 genomes with AMR cassettes, 20 AMR query genes
  Tier 4 -- Long read simulation
             Reuses Tier 1 genomes, 2-10 Kbp queries with 5-15% error
  Tier 5 -- Short read queries
             Reuses Tier 1 genomes, 150/250 bp Illumina-like queries
  Tier 6 -- Plasmid / mobile element queries
             20 genomes with shared plasmid sequences (2-5 Kbp)

All query files are written at divergence levels
[0.0, 0.01, 0.03, 0.05, 0.10, 0.15].
"""

import os
import random
import shutil
import sys
import time

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
random.seed(42)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BASES = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
DIVERGENCE_LEVELS = [0.0, 0.01, 0.03, 0.05, 0.10, 0.15]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")


# ===================================================================
# Utility helpers
# ===================================================================

def random_seq(length):
    """Return a random DNA string of *length* bases."""
    return "".join(random.choice(BASES) for _ in range(length))


def reverse_complement(seq):
    """Return the reverse complement of *seq*."""
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq))


def mutate(seq, rate, indel_fraction=0.3):
    """Introduce point mutations (subs/ins/dels) at the given *rate*.

    *indel_fraction* controls how many of the mutations are indels
    (split roughly evenly between insertions and deletions); the rest
    are substitutions.

    Returns (mutated_seq, n_subs, n_ins, n_dels).
    """
    seq = list(seq)
    n_muts = max(int(len(seq) * rate), 0)
    if n_muts == 0:
        return "".join(seq), 0, 0, 0

    positions = random.sample(range(len(seq)), min(n_muts, len(seq)))
    subs, ins, dels = 0, 0, 0

    sub_threshold = 1.0 - indel_fraction
    ins_threshold = sub_threshold + indel_fraction / 2.0

    for pos in sorted(positions, reverse=True):
        r = random.random()
        if r < sub_threshold:
            orig = seq[pos]
            seq[pos] = random.choice([b for b in BASES if b != orig])
            subs += 1
        elif r < ins_threshold:
            seq.insert(pos, random.choice(BASES))
            ins += 1
        else:
            del seq[pos]
            dels += 1

    return "".join(seq), subs, ins, dels


def mutate_long_read(seq, error_rate):
    """Mutate *seq* with a higher indel fraction typical of long reads.

    Long-read sequencers (ONT / PacBio) have ~60 % indel errors.
    """
    return mutate(seq, error_rate, indel_fraction=0.6)


def mutate_short_read(seq, error_rate):
    """Mutate *seq* with a substitution-dominated profile (Illumina)."""
    return mutate(seq, error_rate, indel_fraction=0.05)


def write_fasta(path, name, seq, mode="w"):
    """Write a single FASTA record to *path* (creating dirs as needed)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode) as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def append_fasta(path, name, seq):
    """Append a single FASTA record to an existing file."""
    with open(path, "a") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def write_fasta_records(path, records):
    """Write many FASTA records to *path*.

    *records* is an iterable of (name, seq) tuples.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def emit_divergence_files(query_dir, tier_label, queries, divergence_levels,
                          mutate_fn=mutate):
    """Write per-divergence query FASTA + truth TSV files.

    *queries*  -- list of (query_id, original_sequence)
    *mutate_fn* -- callable(seq, rate) -> (mutated_seq, subs, ins, dels)
    """
    for div in divergence_levels:
        div_str = f"{div:.2f}".replace(".", "")
        out_fa = os.path.join(query_dir, f"genes_div{div}.fa")
        out_truth = os.path.join(query_dir, f"genes_div{div}_truth.tsv")

        with open(out_fa, "w") as fa, open(out_truth, "w") as tsv:
            tsv.write(
                "query_id\toriginal_length\tmutated_length\tdivergence\t"
                "num_subs\tnum_ins\tnum_dels\tactual_divergence\n"
            )
            for qid, qseq in queries:
                if div == 0.0:
                    mutated = qseq
                    s, ins_count, d = 0, 0, 0
                else:
                    mutated, s, ins_count, d = mutate_fn(qseq, div)

                fa.write(f">{qid}\n")
                for j in range(0, len(mutated), 80):
                    fa.write(mutated[j : j + 80] + "\n")

                total = s + ins_count + d
                actual_div = total / max(len(qseq), 1)
                tsv.write(
                    f"{qid}\t{len(qseq)}\t{len(mutated)}\t{div}\t"
                    f"{s}\t{ins_count}\t{d}\t{actual_div:.6f}\n"
                )

        print(f"    Divergence {div * 100:5.1f}%: {os.path.basename(out_fa)}")


def banner(text):
    width = 64
    print()
    print("=" * width)
    print(f"  {text}")
    print("=" * width)


# ===================================================================
# Tier 1 -- Within-species (E. coli-like)
# ===================================================================

def generate_tier1():
    """50 genomes (~50 Kbp), 95-99 % ANI, 100 gene queries."""
    banner("Tier 1: Within-species (E. coli-like)")

    NUM_GENOMES = 50
    GENOME_SIZE = 50_000       # 50 Kbp (scaled from real 5 Mbp)
    SHARED_CORE_FRAC = 0.90    # 90 % shared core
    CORE_SIZE = int(GENOME_SIZE * SHARED_CORE_FRAC)
    NUM_QUERIES = 100
    QUERY_MIN_LEN = 300
    QUERY_MAX_LEN = 2000

    genome_dir = os.path.join(DATA_DIR, "tier1_genomes")
    query_dir = os.path.join(DATA_DIR, "queries", "tier1")
    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(query_dir, exist_ok=True)

    # Shared core
    core = random_seq(CORE_SIZE)
    print(f"  Shared core: {CORE_SIZE:,} bp ({SHARED_CORE_FRAC * 100:.0f}% of genome)")

    genome_names = []
    genome_seqs = {}

    for i in range(NUM_GENOMES):
        name = f"ecoli_genome_{i:04d}"
        # Species-level divergence: 1-5 % from core  (ANI 95-99 %)
        div_rate = random.uniform(0.01, 0.05)
        mutated_core, _, _, _ = mutate(core, div_rate)
        # Unique flanking regions fill the rest of the genome
        flank_total = GENOME_SIZE - len(mutated_core)
        flank_5 = random_seq(flank_total // 2)
        flank_3 = random_seq(flank_total - flank_total // 2)
        seq = flank_5 + mutated_core + flank_3

        write_fasta(os.path.join(genome_dir, f"{name}.fa"), name, seq)
        genome_names.append(name)
        genome_seqs[name] = seq

    print(f"  Wrote {NUM_GENOMES} genomes to {genome_dir}")

    # ---- queries ----
    queries = []
    truth_records = []
    for i in range(NUM_QUERIES):
        gname = random.choice(genome_names)
        gseq = genome_seqs[gname]
        qlen = random.randint(QUERY_MIN_LEN, min(QUERY_MAX_LEN, len(gseq) - 100))
        start = random.randint(0, len(gseq) - qlen)
        qseq = gseq[start : start + qlen]
        qid = f"tier1_gene_{i:06d}"
        queries.append((qid, qseq))
        truth_records.append((qid, gname, f"{gname}_contig", start, start + qlen, qlen))

    # Raw queries FASTA
    write_fasta_records(os.path.join(query_dir, "genes.fa"), queries)

    # Truth TSV
    truth_path = os.path.join(query_dir, "genes_truth.tsv")
    with open(truth_path, "w") as f:
        f.write("query_id\tgenome\tcontig\tstart\tend\tlength\n")
        for rec in truth_records:
            f.write("\t".join(str(x) for x in rec) + "\n")

    print(f"  {NUM_QUERIES} gene queries ({QUERY_MIN_LEN}-{QUERY_MAX_LEN} bp)")

    # Per-divergence files
    emit_divergence_files(query_dir, "tier1", queries, DIVERGENCE_LEVELS)

    # Return for reuse by tiers 4 & 5
    return genome_names, genome_seqs


# ===================================================================
# Tier 2 -- Cross-species diverse
# ===================================================================

def generate_tier2():
    """30 genomes from 6 species clusters, 80 gene queries."""
    banner("Tier 2: Cross-species diverse")

    NUM_CLUSTERS = 6
    GENOMES_PER_CLUSTER = 5
    GENOME_SIZE = 50_000
    CORE_SIZE = int(GENOME_SIZE * 0.50)  # inter-species shared fraction
    NUM_QUERIES = 80
    QUERY_MIN_LEN = 300
    QUERY_MAX_LEN = 2000

    genome_dir = os.path.join(DATA_DIR, "tier2_genomes")
    query_dir = os.path.join(DATA_DIR, "queries", "tier2")
    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(query_dir, exist_ok=True)

    # Universal core shared among all species (low identity between clusters)
    universal_core = random_seq(CORE_SIZE)
    print(f"  Universal core: {CORE_SIZE:,} bp")
    print(f"  {NUM_CLUSTERS} species clusters x {GENOMES_PER_CLUSTER} genomes each")

    genome_names = []
    genome_seqs = {}

    for cluster_idx in range(NUM_CLUSTERS):
        # Between-cluster divergence: 15-30 % from universal core (ANI 70-85 %)
        cluster_div = random.uniform(0.15, 0.30)
        cluster_core, _, _, _ = mutate(universal_core, cluster_div)

        # Additional cluster-specific sequence
        cluster_unique = random_seq(GENOME_SIZE - len(cluster_core))

        for g in range(GENOMES_PER_CLUSTER):
            name = f"species{cluster_idx}_genome_{g:02d}"
            # Within-cluster divergence: 1-5 % (ANI 95-99 %)
            within_div = random.uniform(0.01, 0.05)
            this_core, _, _, _ = mutate(cluster_core, within_div)
            this_unique, _, _, _ = mutate(cluster_unique, within_div)
            seq = this_core + this_unique

            write_fasta(os.path.join(genome_dir, f"{name}.fa"), name, seq)
            genome_names.append(name)
            genome_seqs[name] = seq

    print(f"  Wrote {len(genome_names)} genomes to {genome_dir}")

    # ---- queries ----
    queries = []
    truth_records = []
    for i in range(NUM_QUERIES):
        gname = random.choice(genome_names)
        gseq = genome_seqs[gname]
        qlen = random.randint(QUERY_MIN_LEN, min(QUERY_MAX_LEN, len(gseq) - 100))
        start = random.randint(0, len(gseq) - qlen)
        qseq = gseq[start : start + qlen]
        qid = f"tier2_gene_{i:06d}"
        queries.append((qid, qseq))
        truth_records.append((qid, gname, f"{gname}_contig", start, start + qlen, qlen))

    write_fasta_records(os.path.join(query_dir, "genes.fa"), queries)

    truth_path = os.path.join(query_dir, "genes_truth.tsv")
    with open(truth_path, "w") as f:
        f.write("query_id\tgenome\tcontig\tstart\tend\tlength\n")
        for rec in truth_records:
            f.write("\t".join(str(x) for x in rec) + "\n")

    print(f"  {NUM_QUERIES} gene queries ({QUERY_MIN_LEN}-{QUERY_MAX_LEN} bp)")
    emit_divergence_files(query_dir, "tier2", queries, DIVERGENCE_LEVELS)


# ===================================================================
# Tier 3 -- AMR gene surveillance
# ===================================================================

def generate_tier3():
    """30 genomes with AMR cassettes; 20 AMR-like query genes."""
    banner("Tier 3: AMR gene surveillance")

    NUM_GENOMES = 30
    GENOME_SIZE = 50_000
    NUM_AMR_GENES = 20
    AMR_MIN_LEN = 500
    AMR_MAX_LEN = 1500
    MAX_CASSETTES_PER_GENOME = 8

    genome_dir = os.path.join(DATA_DIR, "tier3_genomes")
    query_dir = os.path.join(DATA_DIR, "queries", "tier3")
    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(query_dir, exist_ok=True)

    # Create a library of AMR gene cassettes
    amr_cassettes = []
    for k in range(NUM_AMR_GENES):
        length = random.randint(AMR_MIN_LEN, AMR_MAX_LEN)
        amr_cassettes.append((f"AMR_gene_{k:03d}", random_seq(length)))

    print(f"  AMR gene library: {NUM_AMR_GENES} cassettes ({AMR_MIN_LEN}-{AMR_MAX_LEN} bp)")

    # Build genomes with random AMR cassette insertions
    genome_names = []
    genome_seqs = {}
    amr_truth = []  # (genome, amr_gene_id, start, end, exact_or_mutated)

    for i in range(NUM_GENOMES):
        name = f"amr_genome_{i:04d}"
        backbone = random_seq(GENOME_SIZE)

        # Insert a random subset of AMR cassettes
        n_insert = random.randint(1, MAX_CASSETTES_PER_GENOME)
        chosen = random.sample(range(NUM_AMR_GENES), min(n_insert, NUM_AMR_GENES))

        # Collect insertion points (non-overlapping)
        insert_positions = sorted(random.sample(
            range(1000, GENOME_SIZE - 1000), len(chosen)
        ))

        seq_parts = []
        prev_end = 0
        for idx, (pos, amr_idx) in enumerate(zip(insert_positions, chosen)):
            amr_id, amr_seq = amr_cassettes[amr_idx]
            # Optionally mutate the inserted cassette (0-10 %)
            cassette_div = random.uniform(0.0, 0.10)
            if cassette_div > 0:
                inserted_seq, _, _, _ = mutate(amr_seq, cassette_div)
                variant = "mutated"
            else:
                inserted_seq = amr_seq
                variant = "exact"

            seq_parts.append(backbone[prev_end:pos])
            insert_start = sum(len(p) for p in seq_parts)
            seq_parts.append(inserted_seq)
            insert_end = insert_start + len(inserted_seq)
            prev_end = pos

            amr_truth.append((name, amr_id, insert_start, insert_end,
                              variant, f"{cassette_div:.4f}"))

        seq_parts.append(backbone[prev_end:])
        full_seq = "".join(seq_parts)

        write_fasta(os.path.join(genome_dir, f"{name}.fa"), name, full_seq)
        genome_names.append(name)
        genome_seqs[name] = full_seq

    print(f"  Wrote {NUM_GENOMES} genomes to {genome_dir}")
    print(f"  Total AMR insertions: {len(amr_truth)}")

    # Write AMR insertion truth
    amr_truth_path = os.path.join(query_dir, "amr_insertions_truth.tsv")
    os.makedirs(query_dir, exist_ok=True)
    with open(amr_truth_path, "w") as f:
        f.write("genome\tamr_gene\tinsert_start\tinsert_end\tvariant\tcassette_divergence\n")
        for rec in amr_truth:
            f.write("\t".join(str(x) for x in rec) + "\n")

    # ---- queries: the AMR cassettes themselves ----
    queries = [(aid, aseq) for aid, aseq in amr_cassettes]

    write_fasta_records(os.path.join(query_dir, "genes.fa"), queries)

    # Simple truth for queries
    truth_path = os.path.join(query_dir, "genes_truth.tsv")
    with open(truth_path, "w") as f:
        f.write("query_id\tlength\tdescription\n")
        for aid, aseq in amr_cassettes:
            f.write(f"{aid}\t{len(aseq)}\tAMR reference cassette\n")

    print(f"  {NUM_AMR_GENES} AMR query genes")
    emit_divergence_files(query_dir, "tier3", queries, DIVERGENCE_LEVELS)


# ===================================================================
# Tier 4 -- Long read simulation
# ===================================================================

def generate_tier4(tier1_genome_names, tier1_genome_seqs):
    """Long read-like queries (2-10 Kbp) with 5-15 % error (ONT/PacBio)."""
    banner("Tier 4: Long read simulation")

    NUM_QUERIES = 80
    QUERY_MIN_LEN = 2000
    QUERY_MAX_LEN = 10_000
    ERROR_MIN = 0.05
    ERROR_MAX = 0.15

    query_dir = os.path.join(DATA_DIR, "queries", "tier4")
    os.makedirs(query_dir, exist_ok=True)

    print(f"  Reusing {len(tier1_genome_names)} Tier 1 genomes")
    print(f"  {NUM_QUERIES} long-read queries ({QUERY_MIN_LEN}-{QUERY_MAX_LEN} bp)")
    print(f"  Base error rate: {ERROR_MIN * 100:.0f}-{ERROR_MAX * 100:.0f}%")

    # Symlink or note that tier 4 shares tier 1 genomes
    genome_link_note = os.path.join(DATA_DIR, "tier4_genomes_NOTE.txt")
    with open(genome_link_note, "w") as f:
        f.write("Tier 4 reuses Tier 1 genomes from tier1_genomes/\n")

    queries = []
    truth_records = []
    for i in range(NUM_QUERIES):
        gname = random.choice(tier1_genome_names)
        gseq = tier1_genome_seqs[gname]
        max_qlen = min(QUERY_MAX_LEN, len(gseq) - 100)
        if max_qlen < QUERY_MIN_LEN:
            max_qlen = QUERY_MIN_LEN
        qlen = random.randint(QUERY_MIN_LEN, max_qlen)
        start = random.randint(0, max(len(gseq) - qlen, 0))
        qseq = gseq[start : start + qlen]

        # Apply a per-read base error rate
        read_error = random.uniform(ERROR_MIN, ERROR_MAX)
        noisy_seq, s, ins_c, d = mutate_long_read(qseq, read_error)

        qid = f"tier4_longread_{i:06d}"
        queries.append((qid, noisy_seq))
        truth_records.append((qid, gname, start, start + qlen, qlen,
                              f"{read_error:.4f}", s, ins_c, d))

    # Write base queries (already noisy -- these are the "reads")
    write_fasta_records(os.path.join(query_dir, "genes.fa"), queries)

    truth_path = os.path.join(query_dir, "genes_truth.tsv")
    with open(truth_path, "w") as f:
        f.write("query_id\tgenome\tstart\tend\toriginal_length\t"
                "base_error_rate\tnum_subs\tnum_ins\tnum_dels\n")
        for rec in truth_records:
            f.write("\t".join(str(x) for x in rec) + "\n")

    print(f"  Wrote {NUM_QUERIES} long-read queries")

    # Additional divergence on top of the already-noisy reads
    emit_divergence_files(query_dir, "tier4", queries, DIVERGENCE_LEVELS,
                          mutate_fn=mutate_long_read)


# ===================================================================
# Tier 5 -- Short read queries
# ===================================================================

def generate_tier5(tier1_genome_names, tier1_genome_seqs):
    """150 bp and 250 bp Illumina-like queries with 0.1-1 % error."""
    banner("Tier 5: Short read queries")

    READ_LENGTHS = [150, 250]
    READS_PER_LENGTH = 100
    ERROR_MIN = 0.001
    ERROR_MAX = 0.01

    query_dir = os.path.join(DATA_DIR, "queries", "tier5")
    os.makedirs(query_dir, exist_ok=True)

    print(f"  Reusing {len(tier1_genome_names)} Tier 1 genomes")
    print(f"  Read lengths: {READ_LENGTHS}")
    print(f"  {READS_PER_LENGTH} reads per length = {len(READ_LENGTHS) * READS_PER_LENGTH} total")
    print(f"  Error rate: {ERROR_MIN * 100:.1f}-{ERROR_MAX * 100:.1f}%")

    genome_link_note = os.path.join(DATA_DIR, "tier5_genomes_NOTE.txt")
    with open(genome_link_note, "w") as f:
        f.write("Tier 5 reuses Tier 1 genomes from tier1_genomes/\n")

    queries = []
    truth_records = []
    idx = 0

    for rlen in READ_LENGTHS:
        for _ in range(READS_PER_LENGTH):
            gname = random.choice(tier1_genome_names)
            gseq = tier1_genome_seqs[gname]
            if len(gseq) < rlen + 10:
                continue
            start = random.randint(0, len(gseq) - rlen)
            qseq = gseq[start : start + rlen]

            read_error = random.uniform(ERROR_MIN, ERROR_MAX)
            noisy_seq, s, ins_c, d = mutate_short_read(qseq, read_error)

            qid = f"tier5_read{rlen}_{idx:06d}"
            queries.append((qid, noisy_seq))
            truth_records.append((qid, gname, start, start + rlen, rlen,
                                  f"{read_error:.6f}", s, ins_c, d))
            idx += 1

    write_fasta_records(os.path.join(query_dir, "genes.fa"), queries)

    truth_path = os.path.join(query_dir, "genes_truth.tsv")
    with open(truth_path, "w") as f:
        f.write("query_id\tgenome\tstart\tend\tread_length\t"
                "base_error_rate\tnum_subs\tnum_ins\tnum_dels\n")
        for rec in truth_records:
            f.write("\t".join(str(x) for x in rec) + "\n")

    print(f"  Wrote {len(queries)} short-read queries")

    emit_divergence_files(query_dir, "tier5", queries, DIVERGENCE_LEVELS,
                          mutate_fn=mutate_short_read)


# ===================================================================
# Tier 6 -- Plasmid / mobile element queries
# ===================================================================

def generate_tier6():
    """20 genomes with shared plasmid sequences (2-5 Kbp), query the plasmids."""
    banner("Tier 6: Plasmid / mobile element queries")

    NUM_GENOMES = 20
    GENOME_SIZE = 50_000
    NUM_PLASMIDS = 8
    PLASMID_MIN_LEN = 2000
    PLASMID_MAX_LEN = 5000
    MIN_PLASMIDS_PER_GENOME = 1
    MAX_PLASMIDS_PER_GENOME = 5

    genome_dir = os.path.join(DATA_DIR, "tier6_genomes")
    query_dir = os.path.join(DATA_DIR, "queries", "tier6")
    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(query_dir, exist_ok=True)

    # Create a library of plasmid sequences
    plasmids = []
    for k in range(NUM_PLASMIDS):
        length = random.randint(PLASMID_MIN_LEN, PLASMID_MAX_LEN)
        plasmids.append((f"plasmid_{k:03d}", random_seq(length)))

    print(f"  Plasmid library: {NUM_PLASMIDS} plasmids ({PLASMID_MIN_LEN}-{PLASMID_MAX_LEN} bp)")

    genome_names = []
    genome_seqs = {}
    plasmid_truth = []  # (genome, plasmid_id, contig, start, end)

    for i in range(NUM_GENOMES):
        name = f"plasmid_host_{i:04d}"
        chromosome = random_seq(GENOME_SIZE)

        # Pick a random subset of plasmids for this genome
        n_plas = random.randint(MIN_PLASMIDS_PER_GENOME, MAX_PLASMIDS_PER_GENOME)
        chosen_indices = random.sample(range(NUM_PLASMIDS), n_plas)

        # Insert plasmids at random positions in the chromosome
        insert_positions = sorted(random.sample(
            range(2000, GENOME_SIZE - 2000), len(chosen_indices)
        ))

        seq_parts = []
        prev_end = 0
        for pos, pidx in zip(insert_positions, chosen_indices):
            pid, pseq = plasmids[pidx]
            # Slight variation per genome (0-3 %)
            plas_div = random.uniform(0.0, 0.03)
            if plas_div > 0:
                inserted_seq, _, _, _ = mutate(pseq, plas_div)
            else:
                inserted_seq = pseq

            seq_parts.append(chromosome[prev_end:pos])
            insert_start = sum(len(p) for p in seq_parts)
            seq_parts.append(inserted_seq)
            insert_end = insert_start + len(inserted_seq)
            prev_end = pos

            plasmid_truth.append((name, pid, f"{name}_chr", insert_start, insert_end,
                                  len(inserted_seq), f"{plas_div:.4f}"))

        seq_parts.append(chromosome[prev_end:])
        full_seq = "".join(seq_parts)

        write_fasta(os.path.join(genome_dir, f"{name}.fa"), name, full_seq)
        genome_names.append(name)
        genome_seqs[name] = full_seq

    print(f"  Wrote {NUM_GENOMES} genomes to {genome_dir}")
    print(f"  Total plasmid insertions: {len(plasmid_truth)}")

    # Plasmid insertion truth
    plas_truth_path = os.path.join(query_dir, "plasmid_insertions_truth.tsv")
    with open(plas_truth_path, "w") as f:
        f.write("genome\tplasmid_id\tcontig\tinsert_start\tinsert_end\t"
                "inserted_length\tplasmid_divergence\n")
        for rec in plasmid_truth:
            f.write("\t".join(str(x) for x in rec) + "\n")

    # ---- queries: the reference plasmid sequences ----
    queries = [(pid, pseq) for pid, pseq in plasmids]

    write_fasta_records(os.path.join(query_dir, "genes.fa"), queries)

    truth_path = os.path.join(query_dir, "genes_truth.tsv")
    with open(truth_path, "w") as f:
        f.write("query_id\tlength\tdescription\n")
        for pid, pseq in plasmids:
            f.write(f"{pid}\t{len(pseq)}\tplasmid reference sequence\n")

    print(f"  {NUM_PLASMIDS} plasmid query sequences")
    emit_divergence_files(query_dir, "tier6", queries, DIVERGENCE_LEVELS)


# ===================================================================
# Main entry point
# ===================================================================

def main():
    t0 = time.time()

    banner("Dragon Benchmark -- Synthetic Data Generator")
    print(f"  Output directory : {DATA_DIR}")
    print(f"  Divergence levels: {DIVERGENCE_LEVELS}")
    print(f"  Random seed      : 42")

    # Clean previous data so stale files don't accumulate
    if os.path.isdir(DATA_DIR):
        print(f"\n  Cleaning previous data in {DATA_DIR} ...")
        shutil.rmtree(DATA_DIR)
    os.makedirs(DATA_DIR, exist_ok=True)

    # Tier 1 returns its genomes for reuse in tiers 4 & 5
    tier1_names, tier1_seqs = generate_tier1()

    generate_tier2()
    generate_tier3()
    generate_tier4(tier1_names, tier1_seqs)
    generate_tier5(tier1_names, tier1_seqs)
    generate_tier6()

    elapsed = time.time() - t0

    banner("Summary")
    tiers = [
        ("Tier 1", "Within-species (E. coli-like)",    "50 genomes, 100 queries"),
        ("Tier 2", "Cross-species diverse",             "30 genomes,  80 queries"),
        ("Tier 3", "AMR gene surveillance",             "30 genomes,  20 AMR queries"),
        ("Tier 4", "Long read simulation",              "reuses Tier 1, 80 long-read queries"),
        ("Tier 5", "Short read queries",                "reuses Tier 1, 200 short-read queries"),
        ("Tier 6", "Plasmid / mobile element",          "20 genomes,   8 plasmid queries"),
    ]
    for label, desc, stats in tiers:
        print(f"  {label:8s} | {desc:35s} | {stats}")

    print(f"\n  Divergence levels: {DIVERGENCE_LEVELS}")
    print(f"  Total elapsed time: {elapsed:.1f}s")
    print(f"\n  All data written to: {DATA_DIR}")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
