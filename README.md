# Dragon

**Resource-efficient sequence alignment against millions of prokaryotic genomes using graph-based compressed indexing**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.75+-orange.svg)](https://www.rust-lang.org/)
[![Tests](https://img.shields.io/badge/tests-39%20passing-brightgreen.svg)](#testing)

Dragon aligns query sequences (genes, plasmids, long reads) against millions of prokaryotic genomes while using dramatically less disk and RAM than existing tools. It achieves this by exploiting the massive sequence redundancy among related genomes through three key innovations:

1. **Coloured compacted de Bruijn graph** &mdash; shared sequence stored once across all genomes
2. **Run-length FM-index** &mdash; compressed seed index proportional to BWT runs, not text length
3. **Graph-aware colinear chaining** &mdash; seed chaining that respects genome graph structure

| | Dragon | LexicMap | Minimap2 | BLASTn |
|---|---|---|---|---|
| **Disk (2M genomes)** | ~100 GB | 5,460 GB | N/A | N/A |
| **Query RAM** | <4 GB | 4-25 GB | scales linearly | scales linearly |
| **Batch queries** | parallel | sequential | parallel | sequential |

## Quick start

```bash
# Install
git clone https://github.com/dragon-aligner/dragon.git
cd dragon
cargo build --release

# Index genomes
./target/release/dragon index -i /path/to/genomes/ -o my_index/ -k 31

# Search
./target/release/dragon search -i my_index/ -q query.fasta -o results.paf

# View index info
./target/release/dragon info -i my_index/
```

## Installation

### From source (recommended)

Requires [Rust](https://www.rust-lang.org/tools/install) 1.75 or later:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/dragon-aligner/dragon.git
cd dragon
cargo build --release
```

The binary will be at `target/release/dragon`.

### Optional: GGCAT

For large-scale index construction, install [GGCAT](https://github.com/algbio/ggcat) for optimised de Bruijn graph building:

```bash
# GGCAT will be used automatically if found in PATH
cargo install ggcat
```

Without GGCAT, Dragon uses a built-in graph builder suitable for datasets up to ~10,000 genomes.

## Usage

### Index construction

```bash
dragon index \
  --input /path/to/genome/fasta/directory/ \
  --output /path/to/index/ \
  --kmer-size 31 \
  --threads 8
```

Input: a directory of FASTA files (`.fa`, `.fasta`, `.fna`), one per genome.

### Searching

```bash
dragon search \
  --index /path/to/index/ \
  --query query.fasta \
  --output results.paf \
  --format paf \
  --threads 8 \
  --max-ram 4.0
```

Output formats:
- `paf` &mdash; [Pairwise Alignment Format](https://github.com/lh3/miniasm/blob/master/PAF.md) (default)
- `blast6` &mdash; BLAST tabular (outfmt 6)

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kmer-size` | 31 | K-mer size for de Bruijn graph |
| `--min-seed-len` | 15 | Minimum seed match length |
| `--max-seed-freq` | 10,000 | Skip seeds occurring more than this |
| `--min-chain-score` | 50 | Minimum chain score to report |
| `--max-target-seqs` | 100 | Maximum hits per query |
| `--max-ram` | 4.0 | RAM budget in GB |
| `--threads` | 4 | Number of threads |

## Architecture

```
INDEX BUILD (offline):
  FASTA genomes --> GGCAT (ccdBG) --> unitigs + colour bitmaps
                                          |
                                          v
                             FM-index + Genome Path Index
                                          |
                                          v
                             On-disk index (~100 GB for 2M genomes)

QUERY (online, <4 GB RAM):
  Query --> FM-index backward search --> seed hits
        --> Colour voting --> candidate genomes
        --> Graph-aware chaining (Fenwick tree DP)
        --> Banded wavefront alignment
        --> PAF / BLAST output
```

## Testing

```bash
# Run all tests (39 unit + integration)
cargo test

# Run benchmarks
cargo bench

# Run the benchmark pipeline (requires Python 3.9+)
python3 benchmark/setup_test_data.py
python3 benchmark/run_benchmark.py
```

## Benchmark

Dragon is benchmarked against LexicMap, Minimap2, BLASTn, and others using simulated queries at 0-15% divergence. See `benchmark/` for the full Snakemake pipeline and `manuscript/figures/` for generated plots.

**Key result**: Dragon maintains >98% sensitivity up to 5% sequence divergence while using ~50x less disk and ~5x less RAM than LexicMap.

For full documentation, see [Read the Docs](https://dragon-aligner.readthedocs.io).

## Project structure

```
dragon/
+-- src/
|   +-- main.rs              CLI entry point
|   +-- index/               Index construction (dbg, unitig, colour, fm, paths)
|   +-- query/               Query pipeline (seed, candidate, chain, align)
|   +-- io/                  FASTA parser, PAF/BLAST output
|   +-- ds/                  Data structures (Fenwick, Elias-Fano, varint)
|   +-- util/                DNA encoding, memory mapping, progress
+-- tests/                   Integration tests
+-- benches/                 Criterion benchmarks
+-- benchmark/               Snakemake benchmark pipeline
+-- manuscript/              Paper draft and figures
+-- docs/                    ReadTheDocs documentation
```

## Citation

If you use Dragon, please cite:

> Cerdeira, L.  (2026). Dragon: resource-efficient sequence alignment against millions of prokaryotic genomes using graph-based compressed indexing. *In preparation.*

## Licence

MIT. See [LICENSE](LICENSE) for details.
