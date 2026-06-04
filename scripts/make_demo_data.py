#!/usr/bin/env python3
"""
Generate a tiny, self-contained demo dataset for testing Dragon's Zarr backend.

Creates 6 small "genomes" (~30 kb each) that share a common core region (so the
coloured de Bruijn graph stores it once and colours it by all 6), plus a small
"resistance gene" carried by only 2 of them. Two query records exercise both
cases: the core (hits all 6 genomes) and the gene (hits the 2 carriers).

Output:
  demo/genomes/genome_{1..6}.fa
  demo/query.fa

Deterministic (fixed seed) so the demo is reproducible.
"""
import os, random, argparse

def rnd(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))

def mutate(s, d, rng):
    s = list(s)
    for i in range(len(s)):
        if rng.random() < d:
            s[i] = rng.choice("ACGT")
    return "".join(s)

def wrap(s, w=70):
    return "\n".join(s[i:i+w] for i in range(0, len(s), w))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", default="demo")
    args = ap.parse_args()
    rng = random.Random(20260604)

    gdir = os.path.join(args.outdir, "genomes")
    os.makedirs(gdir, exist_ok=True)

    # Shared CORE region present (near-identically) in every genome.
    core = rnd(15000, rng)
    # A "resistance gene" carried by genomes 3 and 5 only.
    gene = rnd(1200, rng)

    n_genomes = 6
    for g in range(1, n_genomes + 1):
        # each genome: unique left flank + core (lightly diverged) + unique right flank
        left = rnd(rng.randint(6000, 9000), rng)
        right = rnd(rng.randint(6000, 9000), rng)
        core_g = mutate(core, 0.01, rng)            # ~1% intra-species drift
        parts = [left, core_g]
        if g in (3, 5):
            parts.append(mutate(gene, 0.005, rng))    # the AMR-like gene
        parts.append(right)
        seq = "".join(parts)
        with open(os.path.join(gdir, f"genome_{g}.fa"), "w") as f:
            f.write(f">genome_{g} demo Staphylo-like assembly\n{wrap(seq)}\n")

    # Queries: (1) a 1 kb slice of the core -> should hit all 6 genomes;
    #          (2) the resistance gene -> should hit genomes 3 and 5 only.
    core_slice = mutate(core[4000:5000], 0.02, rng)   # 2% diverged query
    gene_query = mutate(gene, 0.01, rng)
    with open(os.path.join(args.outdir, "query.fa"), "w") as f:
        f.write(f">core_fragment expect: all 6 genomes\n{wrap(core_slice)}\n")
        f.write(f">resistance_gene expect: genome_3 and genome_5\n{wrap(gene_query)}\n")

    print(f"wrote {n_genomes} genomes to {gdir}/ and 2 queries to {args.outdir}/query.fa")

if __name__ == "__main__":
    main()
