<p align="center">
  <img src="assets/dragon-logo.png" alt="Dragon Logo" width="300">
</p>

<p align="center"><strong>Dragon: a cloud-native, signal-aware aligner for surveillance-scale microbial genomics</strong></p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19478347.svg)](https://doi.org/10.5281/zenodo.19478347)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.75+-orange.svg)](https://www.rust-lang.org/)
[![Tests](https://img.shields.io/badge/tests-99%20passing-brightgreen.svg)](#testing)

Dragon aligns query sequences (genes, plasmids, long/short reads, raw nanopore signal) against millions of prokaryotic genomes while using dramatically less disk and RAM than existing tools.

It exploits the redundancy among related genomes through:

1. **Coloured compacted de Bruijn graph** — shared sequence stored once across all genomes (built via [GGCAT](https://github.com/algbio/ggcat) for >10K-genome scale).
2. **FM-index over concatenated unitigs** — variable-length seed extension via backward search.
3. **Graph-aware colinear chaining** — anchor chaining that respects the de Bruijn graph topology, with ML-weighted seed scoring.
4. **Roaring-bitmap colour index** — O(1) genome-membership lookups per unitig.
5. **Streaming, mmap-friendly on-disk format** (`paths.bin v2`) — O(1) cold-load via per-genome offset table; queries fault in only the chunks they touch.
6. **Cloud-native Zarr backend** (`dragon export-zarr`) — chunked + Zstd-compressed; readable from any Zarr-aware tool (zarr-python, xarray) and queryable directly from S3 / GCS.

| | Dragon | LexicMap | Minimap2 | BLASTn |
|---|---|---|---|---|
| **Disk (2M genomes)** | ~100 GB | 5,460 GB | scales linearly | scales linearly |
| **Query RAM** | <4 GB | 4–25 GB | scales linearly | scales linearly |
| **Multi-shard search** | Yes `--shard` | No | No | No |
| **Cloud-native (S3 random read)** | Yes / Zarr v3 | No | No | No |
| **Raw nanopore signal search** | Yes | No | No | No |
| **Per-species surveillance summary** | Yes | No | No | No |
| **Hardware profile (laptop mode)** | Yes | No | partial | partial |

A 16,000-genome demo index lives at `s3://dragon-zarr/saureus/b1/` (eu-west-2, public-read). No credentials needed:

```bash
pip install 'zarr>=3.0' s3fs numcodecs
python scripts/zarr_demo.py s3://dragon-zarr/saureus/b1
```

## Quick start

```bash
# Install
git clone https://github.com/lcerdeira/dragon.git
cd dragon
cargo build --release

# Index a directory of genomes
./target/release/dragon index -i /path/to/genomes/ -o my_index/ -k 31 -j 8

# Search
./target/release/dragon search -i my_index/ -q query.fa -o results.paf

# Search across multiple shards (for indices split by RAM/quota)
./target/release/dragon search -i shard_a/ --shard shard_b/ --shard shard_c/ \
    -q query.fa -o results.paf

# Export to Zarr for cloud deployment
./target/release/dragon export-zarr -i my_index/ -o my_index.zarr/

# Query a Zarr store (local or s3://)
./target/release/dragon search-zarr -z my_index.zarr/ -q query.fa
```

## Installation

Requires [Rust](https://www.rust-lang.org/tools/install) 1.75 or later:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/lcerdeira/dragon.git
cd dragon
cargo build --release
```

The binary is at `target/release/dragon`. Install system-wide with `cargo install --path .` or copy the binary into your `$PATH`.

### Optional: GGCAT

For databases >10K genomes, install [GGCAT](https://github.com/algbio/ggcat):

```bash
git clone https://github.com/algbio/ggcat
cd ggcat
cargo build --release
cp target/release/ggcat ~/.cargo/bin/   # or anywhere on PATH
```

Dragon detects GGCAT automatically. Without it, the built-in graph builder handles small datasets (~thousands of genomes).

## Subcommands

| Command | Purpose |
| --- | --- |
| `dragon index` | Build a Dragon index from a directory of FASTA genomes |
| `dragon search` | Align query sequences against an index (single or multi-shard) |
| `dragon info` | Print index metadata (genome count, k-mer size, on-disk size) |
| `dragon download` | Download genomes (RefSeq, AllTheBacteria) or pre-built indices |
| `dragon update` | Add new genomes as a lightweight overlay (no full rebuild) |
| `dragon compact` | Merge base + overlays back into a single optimised index |
| `dragon summarize` | Produce a per-species prevalence/identity report from PAF output |
| `dragon export-zarr` | Export an index as a Zarr v3 store (cloud-native, chunked) |
| `dragon search-zarr` | Pattern-search a Zarr-backed index (local path or `s3://` URI) |
| `dragon signal-index` | Build a signal-level index from FASTA via a pore model |
| `dragon signal-search` | Align raw nanopore current signals (TSV/CSV/SLOW5) directly |

Run `dragon <subcommand> --help` for the full option list.

### Key search options

| Option | Default | Description |
| --- | --- | --- |
| `--index` | required | Primary index directory |
| `--shard` (repeatable) | — | Additional shard directories for multi-index search |
| `--query` | required | Query FASTA/FASTQ file |
| `--format` | `paf` | Output: `paf`, `blast6`, `summary`, `gfa` |
| `--profile` | `workstation` | `laptop` (≤8 GB RAM, 4 threads) or `workstation` (full resources) |
| `--threads` | 4 | CPU threads |
| `--max-ram` | 4.0 | RAM budget in GB |
| `--min-seed-len` | 15 | Minimum seed match length |
| `--min-identity` | 0.7 | Minimum alignment identity to report |
| `--min-query-coverage` | 0.3 | Minimum query coverage to report |
| `--max-target-seqs` | 10 | Hits per query |
| `--no-ml` | off | Disable learned seed scoring (use raw match length) |

### Output formats

- **PAF** — minimap2-compatible pairwise alignment.
- **BLAST6** — BLAST-tabular `outfmt 6`.
- **summary** — per-species prevalence + identity distribution (surveillance-ready).
- **gfa** — graph-context unitigs around each hit (for mobile-element analysis).

## Architecture

```
INDEX BUILD (offline)
  FASTA genomes ──► GGCAT ccdBG ──► unitigs.fa + colormap.dat
                                          │
                                          ▼
                       fm_index.bin   colors.drgn (RoaringBitmaps)
                                          │
                                          ▼
                       paths.bin v2  (mmap'd, varint-encoded per-genome blobs)
                                          │
                                          ▼
                       specificity.drgn   metadata.json
                                          │
                                          ▼
                  ┌───────────────┴───────────────┐
                  ▼                               ▼
         on-disk Dragon index            dragon export-zarr
         (~100 GB / 2M genomes)          ──►  Zarr v3 store
                                              (chunked + Zstd, S3/GCS-ready)

QUERY (online, <4 GB RAM)
  Query FASTA ──► FM-index backward search ──► variable-length seeds
              ──► colour voting (RoaringBitmap) ──► candidate genomes
              ──► ML-weighted graph-aware chaining
              ──► banded WFA alignment + path-walking ref extraction
              ──► PAF / BLAST6 / summary / gfa output

SIGNAL SEARCH
  Raw nanopore pA ──► median-MAD normalise ──► 16-level discretise
                  ──► signal-FM-index backward search
                  ──► per-genome score ──► TSV
```

## Testing

```bash
cargo test --lib                # 99 unit tests
cargo test                      # + integration tests
cargo bench                     # criterion micro-benchmarks
```

## Documentation

Full documentation: <https://dragon-aligner.readthedocs.io>

Key references:

- [Quickstart](docs/getting-started/quickstart.md)
- [CLI reference](docs/api/cli.md)
- [Architecture overview](docs/architecture/overview.md)
- [Performance tuning](docs/user-guide/performance-tuning.md)

## Project structure

```
dragon/
├── src/
│   ├── main.rs              CLI entry point (11 subcommands)
│   ├── index/               Index construction
│   │   ├── dbg.rs           ccdBG via GGCAT (fallback: internal builder)
│   │   ├── unitig.rs        2-bit packed unitig encoding
│   │   ├── color.rs         RoaringBitmap colour index
│   │   ├── ggcat_colors.rs  GGCAT binary colormap → colors.drgn (no TSV)
│   │   ├── fm.rs            Suffix array + binary search FM-index
│   │   ├── paths.rs         Genome path index (legacy bincode loader)
│   │   ├── paths_v2.rs      Mmap-friendly v2 format (default for new builds)
│   │   ├── specificity.rs   Per-genome private-unitig sets
│   │   ├── auto_batch.rs    Auto-split large collections into overlay batches
│   │   ├── update.rs        Incremental overlay addition
│   │   └── zarr_backend.rs  Zarr v3 export + ZarrFmIndex / ZarrColorIndex
│   ├── query/               Query pipeline
│   │   ├── seed.rs          Variable-length backward search
│   │   ├── chain.rs         Graph-aware chaining + ML scoring + containment ranking
│   │   ├── align.rs         Banded WFA alignment
│   │   ├── containment.rs   K-mer containment ranking
│   │   ├── direct_align.rs  Direct alignment to candidate genome subsequences
│   │   └── mod.rs           Multi-shard orchestration
│   ├── signal/              Raw-current nanopore search (signal-index, signal-search)
│   ├── io/                  FASTA/FASTQ + PAF/BLAST6/GFA output
│   ├── ds/                  Fenwick tree, Elias-Fano, varint codecs
│   └── util/                DNA encoding, mmap, colorspace (SOLiD), progress
├── scripts/
│   ├── zarr_demo.py         Read a Zarr store from local or s3:// (paper §4.8 demo)
│   └── train_seed_scorer.py Train the logistic-regression ML seed weights
├── tests/                   Integration tests
├── benches/                 Criterion micro-benchmarks
└── docs/                    Sphinx + Read the Docs
```

Benchmarks, manuscript, and AWS build scripts live in the companion repo `lcerdeira/dragon-private` (private until publication).

## Citation

> Cerdeira, L. (2026). Dragon: a cloud-native, signal-aware aligner for surveillance-scale microbial genomics. *In preparation.*

## Licence

MIT. See [LICENSE](LICENSE) for details.
