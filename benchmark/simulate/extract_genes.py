#!/usr/bin/env python3
"""Extract random genes (subsequences) from reference genomes as query sequences."""

import argparse
import os
import random
from pathlib import Path


def read_fasta(path):
    """Simple FASTA reader."""
    sequences = {}
    current_name = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name:
        sequences[current_name] = "".join(current_seq)
    return sequences


def main():
    parser = argparse.ArgumentParser(description="Extract random gene-like subsequences")
    parser.add_argument("--genome-dir", required=True, help="Directory of FASTA genomes")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--truth", required=True, help="Ground truth TSV file")
    parser.add_argument("--num-genes", type=int, default=1000)
    parser.add_argument("--min-length", type=int, default=500)
    parser.add_argument("--max-length", type=int, default=5000)
    parser.add_argument("--num-genomes", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)

    # List genome files
    genome_dir = Path(args.genome_dir)
    extensions = {".fa", ".fasta", ".fna", ".fsa"}
    genome_files = sorted(
        [f for f in genome_dir.iterdir() if f.suffix.lower() in extensions]
    )

    if len(genome_files) == 0:
        print(f"No FASTA files found in {genome_dir}")
        return

    # Sample genomes
    selected = random.sample(genome_files, min(args.num_genomes, len(genome_files)))

    # Load sequences
    genome_seqs = {}
    for gf in selected:
        seqs = read_fasta(gf)
        genome_name = gf.stem
        for seq_name, seq in seqs.items():
            genome_seqs[(genome_name, seq_name)] = seq

    if not genome_seqs:
        print("No sequences loaded!")
        return

    # Extract random subsequences
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    keys = list(genome_seqs.keys())

    with open(args.output, "w") as fa_out, open(args.truth, "w") as tsv_out:
        tsv_out.write("query_id\tgenome\tcontig\tstart\tend\tlength\n")

        for i in range(args.num_genes):
            genome_name, seq_name = random.choice(keys)
            seq = genome_seqs[(genome_name, seq_name)]

            length = random.randint(args.min_length, min(args.max_length, len(seq)))
            start = random.randint(0, max(0, len(seq) - length))
            end = start + length

            subseq = seq[start:end]
            query_id = f"gene_{i:06d}"

            fa_out.write(f">{query_id}\n")
            # Write in 80-char lines
            for j in range(0, len(subseq), 80):
                fa_out.write(subseq[j : j + 80] + "\n")

            tsv_out.write(
                f"{query_id}\t{genome_name}\t{seq_name}\t{start}\t{end}\t{length}\n"
            )

    print(f"Extracted {args.num_genes} gene queries from {len(selected)} genomes")
    print(f"Output: {args.output}")
    print(f"Truth: {args.truth}")


if __name__ == "__main__":
    main()
