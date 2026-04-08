# Indexing

## Overview

Dragon's index construction transforms a collection of genome FASTA files into a compact, queryable index through five stages:

1. **De Bruijn graph construction** &mdash; builds a coloured compacted de Bruijn graph (ccdBG) from all genomes
2. **Unitig encoding** &mdash; encodes unitig sequences in 2-bit packed format
3. **Colour index** &mdash; creates Roaring Bitmap mappings from unitigs to genomes
4. **FM-index** &mdash; builds a run-length FM-index over concatenated unitigs
5. **Genome path index** &mdash; records each genome's traversal through the graph

## Command

```bash
dragon index [OPTIONS] --input <DIR> --output <DIR>
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `--input`, `-i` | Directory containing genome FASTA files |
| `--output`, `-o` | Output directory for the index |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--kmer-size`, `-k` | 31 | K-mer size for the de Bruijn graph. Must be odd. Range: 15-31. |
| `--threads`, `-j` | 4 | Number of threads for parallel processing |

## Input format

- One FASTA file per genome
- Supported extensions: `.fa`, `.fasta`, `.fna`, `.fsa`
- Genomes may contain multiple contigs/chromosomes
- Ambiguous bases (N, R, Y, etc.) are treated as A

## Choosing k-mer size

| k | Sensitivity | Specificity | Index size | Use case |
|---|-------------|-------------|------------|----------|
| 15 | Highest | Lowest | Larger | Highly divergent queries (>10%) |
| 21 | High | Medium | Medium | General purpose |
| 31 | Medium | Highest | Smallest | Closely related genomes (<5% divergence) |

**Default (k=31)** is recommended for most prokaryotic applications.

## Index files

| File | Description | Typical size (2M genomes) |
|------|-------------|--------------------------|
| `fm_index.bin` | Run-length FM-index | 2-4 GB |
| `colors.drgn` | Roaring Bitmap colour index | 50-100 GB |
| `paths.bin` | Genome path index | 8-15 GB |
| `unitigs.fa` | Unitig FASTA sequences | 1-2 GB |
| `colors.tsv` | Unitig-to-genome mapping (text) | varies |
| `metadata.json` | Index metadata | <1 KB |

## Resource requirements

| Database scale | Build time | Build RAM | Index size |
|----------------|-----------|-----------|------------|
| 500 genomes | ~10 seconds | <1 GB | ~1.5 GB |
| 85K genomes | ~1 hour | ~8 GB | ~15 GB |
| 2.34M genomes | ~12 hours | ~64 GB | ~100 GB |

## GGCAT integration

For databases larger than ~10,000 genomes, Dragon automatically uses [GGCAT](https://github.com/algbio/ggcat) if available in `PATH`. GGCAT provides:

- 5-39x faster graph construction than alternatives
- Lower memory usage via minimiser-based partitioning
- Native Rust implementation

Without GGCAT, Dragon falls back to a built-in graph builder that constructs each unique k-mer as a minimal unitig. This works correctly but produces a less compacted graph.
