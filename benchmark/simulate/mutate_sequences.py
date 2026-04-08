#!/usr/bin/env python3
"""Introduce controlled mutations into DNA sequences at specified divergence levels."""

import argparse
import random


BASES = "ACGT"


def mutate_sequence(seq, divergence, sub_rate=0.7, ins_rate=0.2, del_rate=0.1):
    """Apply random mutations to a sequence at a given divergence rate.

    Returns (mutated_seq, num_substitutions, num_insertions, num_deletions).
    """
    seq = list(seq.upper())
    n_mutations = int(len(seq) * divergence)

    subs, ins, dels = 0, 0, 0
    positions = random.sample(range(len(seq)), min(n_mutations, len(seq)))
    positions.sort(reverse=True)  # Process from end to preserve indices

    for pos in positions:
        r = random.random()
        if r < sub_rate:
            # Substitution
            original = seq[pos]
            alternatives = [b for b in BASES if b != original]
            seq[pos] = random.choice(alternatives)
            subs += 1
        elif r < sub_rate + ins_rate:
            # Insertion (1-3 bases)
            insert_len = random.choices([1, 2, 3], weights=[0.7, 0.2, 0.1])[0]
            insert_seq = [random.choice(BASES) for _ in range(insert_len)]
            seq[pos:pos] = insert_seq
            ins += 1
        else:
            # Deletion (1-3 bases)
            del_len = random.choices([1, 2, 3], weights=[0.7, 0.2, 0.1])[0]
            del seq[pos : pos + del_len]
            dels += 1

    return "".join(seq), subs, ins, dels


def read_fasta(path):
    """Read FASTA file, return list of (name, sequence)."""
    records = []
    current_name = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    records.append((current_name, "".join(current_seq)))
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name:
        records.append((current_name, "".join(current_seq)))
    return records


def main():
    parser = argparse.ArgumentParser(description="Mutate sequences at controlled divergence")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--truth", required=True, help="Mutation truth TSV file")
    parser.add_argument("--divergence", type=float, required=True, help="Target divergence (0-1)")
    parser.add_argument("--sub-rate", type=float, default=0.70, help="Fraction of mutations that are substitutions")
    parser.add_argument("--ins-rate", type=float, default=0.20, help="Fraction that are insertions")
    parser.add_argument("--del-rate", type=float, default=0.10, help="Fraction that are deletions")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)

    records = read_fasta(args.input)

    with open(args.output, "w") as fa_out, open(args.truth, "w") as tsv_out:
        tsv_out.write("query_id\toriginal_length\tmutated_length\tdivergence\t"
                      "num_subs\tnum_ins\tnum_dels\tactual_divergence\n")

        for name, seq in records:
            mutated, subs, ins, dels = mutate_sequence(
                seq, args.divergence,
                args.sub_rate, args.ins_rate, args.del_rate
            )

            total_ops = subs + ins + dels
            actual_div = total_ops / max(len(seq), 1)

            fa_out.write(f">{name}\n")
            for j in range(0, len(mutated), 80):
                fa_out.write(mutated[j : j + 80] + "\n")

            tsv_out.write(
                f"{name}\t{len(seq)}\t{len(mutated)}\t{args.divergence}\t"
                f"{subs}\t{ins}\t{dels}\t{actual_div:.6f}\n"
            )

    print(f"Mutated {len(records)} sequences at {args.divergence*100}% divergence")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
