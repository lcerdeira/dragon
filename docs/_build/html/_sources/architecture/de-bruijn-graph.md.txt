# Coloured compacted de Bruijn graph

## What is a de Bruijn graph?

A de Bruijn graph represents a DNA sequence collection where:
- Each **node** is a unique k-mer (substring of length k)
- Each **edge** connects k-mers that overlap by k-1 bases
- A **compacted** graph merges non-branching paths into **unitigs** (maximal linear chains)
- A **coloured** graph annotates each unitig with the set of genomes containing it

## Why use a de Bruijn graph?

For 2 million prokaryotic genomes totalling ~10 Tbp of sequence:
- **Without graph**: store 10 Tbp of sequence = ~2.5 TB at 2 bits/base
- **With graph**: store ~5 Gbp of unique unitig sequence = ~1.2 GB at 2 bits/base

The graph collapses shared content, achieving a **~2,000-fold reduction** in sequence storage.

## Construction

Dragon uses [GGCAT](https://github.com/algbio/ggcat) for graph construction:

1. **Input**: N genome FASTA files
2. **K-mer extraction**: extract all k-mers (default k=31) from all genomes
3. **Graph construction**: build the de Bruijn graph, merge non-branching paths into unitigs
4. **Colour annotation**: record which genomes contain each unitig
5. **Output**: unitig FASTA file + colour mapping

### GGCAT advantages

- Written in Rust (same as Dragon)
- 5-39x faster than Bifrost
- Low memory via minimiser-based partitioning
- Native colour support

### Fallback builder

Without GGCAT, Dragon uses a built-in builder that:
- Collects all k-mers with genome IDs via a BTreeMap
- Outputs each unique k-mer as a minimal unitig
- Works correctly but produces less compacted graphs
- Suitable for datasets up to ~10,000 genomes

## Colour storage

Each unitig's colour set (the genomes containing it) is stored as a **Roaring Bitmap** — a compressed integer set that efficiently handles:
- **Dense ranges**: consecutive genome IDs (e.g., all *E. coli* genomes)
- **Sparse sets**: scattered genome IDs
- **Fast operations**: intersection, union, cardinality in O(n/64) time

The colour index is stored as a memory-mapped file with an offset table for O(1) per-unitig access.
