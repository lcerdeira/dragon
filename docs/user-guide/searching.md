# Searching

## Overview

Dragon's query pipeline finds sequences in the indexed database through four stages:

1. **Seed finding** &mdash; FM-index backward search with variable-length extension
2. **Candidate filtering** &mdash; colour-based voting to identify promising genomes
3. **Colinear chaining** &mdash; graph-aware dynamic programming to find the best alignment chain
4. **Alignment** &mdash; banded wavefront alignment for base-level accuracy

## Command

```bash
dragon search [OPTIONS] --index <DIR> --query <FILE>
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `--index`, `-i` | Path to Dragon index directory |
| `--query`, `-q` | Query FASTA or FASTQ file |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--output`, `-o` | `-` (stdout) | Output file path |
| `--format`, `-f` | `paf` | Output format: `paf` or `blast6` |
| `--threads`, `-j` | 4 | Number of parallel threads |
| `--max-ram` | 4.0 | Maximum RAM budget in GB |
| `--min-seed-len` | 15 | Minimum seed match length |
| `--max-seed-freq` | 10,000 | Skip seeds more frequent than this |
| `--min-chain-score` | 50 | Minimum chain score to report |
| `--max-target-seqs` | 100 | Maximum target genomes per query |

## Query types

Dragon handles various query types:

| Query type | Typical length | Recommended parameters |
|-----------|----------------|----------------------|
| Single gene | 500-5,000 bp | defaults |
| 16S rRNA | ~1,500 bp | `--max-seed-freq 50000` (highly conserved) |
| Plasmid | 2,000-200,000 bp | defaults |
| Long read (ONT/PacBio) | 1,000-50,000 bp | `--min-seed-len 12` |
| AMR gene panel | batch of 1,000+ | `--threads 16` |

## Seed finding details

Dragon uses **variable-length seed matching**: rather than searching for fixed k-mers, it extends each FM-index backward search character by character until the suffix array interval is empty. This naturally produces longer, more specific seeds in conserved regions and shorter seeds in variable regions.

Seeds with suffix array interval width exceeding `--max-seed-freq` are discarded as too repetitive (e.g., rRNA-derived seeds that match millions of genomes).

## Candidate filtering

For each seed hit, Dragon retrieves the unitig's Roaring Bitmap to identify which genomes contain it. Genomes accumulate "votes" proportional to the number of distinct unitig hits. Only genomes exceeding a vote threshold proceed to the chaining stage, dramatically reducing computation.

## Chaining

Dragon performs colinear chaining along each candidate genome's **path through the de Bruijn graph**, not in linear coordinate space. This provides two advantages:

1. Seeds are chained in graph topology order, naturally handling structural rearrangements
2. The Fenwick tree DP runs in O(h log h) time where h is the number of anchors

## Memory management

Dragon keeps query-time RAM below `--max-ram` via:

- **Memory-mapped index files** &mdash; only active pages consume RAM
- **On-demand colour loading** &mdash; Roaring Bitmaps deserialised per-unitig
- **Streaming query processing** &mdash; queries processed one at a time
- **Configurable batch size** &mdash; parallel queries share the same index
