<p align="center">
  <img src="assets/dragon-logo.png" alt="Dragon Logo" width="300">
</p>

<p align="center"><strong>Dragon: a cloud-native, signal-aware aligner for surveillance-scale microbial genomics</strong></p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19478347.svg)](https://doi.org/10.5281/zenodo.19478347)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.75+-orange.svg)](https://www.rust-lang.org/)
[![Tests](https://img.shields.io/badge/tests-163%20passing-brightgreen.svg)](#testing)

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

### Try it in 60 seconds

Query a live Dragon index hosted on public S3, **over HTTPS, with no credentials and nothing downloaded up front** — `search-zarr` fetches only the compressed chunks it needs:

```bash
cargo build --release
curl -sO https://dragon-zarr.s3.eu-west-2.amazonaws.com/demo/query.fa
./target/release/dragon search-zarr \
    --zarr https://dragon-zarr.s3.eu-west-2.amazonaws.com/demo/index.zarr \
    -q query.fa
# core_fragment   -> found in all 6 demo genomes (containment ~0.94)
# resistance_gene -> found only in the 2 carrier genomes (containment ~0.91)
```

Or run the whole `index → export-zarr → search-zarr` pipeline locally on bundled demo data (one command, no S3):

```bash
bash scripts/zarr_quickstart.sh
```

A Python/Zarr-native demo (zarr-python, s3fs, xarray) against a 16,000-genome *S. aureus* shard also works without credentials:

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
| `dragon migrate-paths` | Stream-convert a legacy `paths.bin` to the mmap-friendly format |

Run `dragon <subcommand> --help` for the full option list.

### Key search options

| Option | Default | Description |
| --- | --- | --- |
| `--index` | required | Primary index directory |
| `--shard` (repeatable) | — | Additional shard directories for multi-index search |
| `--query` | required | Query FASTA/FASTQ file |
| `--preset` | `default` | Tuning bundle: `default`, `cross-species`, `amr`, `fast` (see below) |
| `--format` | `paf` | Output: `paf`, `blast6`, `summary`, `gfa` |
| `--profile` | `workstation` | `laptop` (≤8 GB RAM, 4 threads) or `workstation` (full resources) |
| `--threads` | 4 | CPU threads |
| `--max-ram` | 4.0 | RAM budget in GB |
| `--min-seed-len` | 15 | Minimum seed match length |
| `--min-identity` | 0.7 | Minimum alignment identity to report |
| `--min-query-coverage` | 0.3 | Minimum query coverage to report |
| `--min-score-ratio` | 0.1 | Keep hits scoring ≥ ratio × best hit |
| `--max-target-seqs` | 10 | Hits per query (`0` = uncapped) |
| `--no-batch` | off | Disable batch-query KmerCache (per-query processing) |
| `--no-parallel-shards` | off | Load shards one at a time (lower peak RAM) |
| `--no-ml` | off | Disable learned seed scoring (use raw match length) |

### Search presets

`--preset` sets sensible defaults for common scenarios; explicit flags always override the preset.

| Preset | What it does | Use when |
| --- | --- | --- |
| `default` | Balanced: `--max-target-seqs 10`, `--min-identity 0.7`. Solid 31-mer seeds with a pigeonhole multi-anchor fallback that fires only when solid seeding is sparse (short reads / high divergence). | Within-species search, genes, plasmids, reads |
| `cross-species` | Short k=7–8 anchors (5×), `--min-identity 0.5` — recovers homologs at 70–85% ANI that solid 31-mers miss. | Cross-genus / distant-homolog detection |
| `amr` | `--max-target-seqs 1000`, `--min-identity 0.9`, `--min-query-coverage 0.8`, batch KmerCache + parallel shards. | AMR gene panels against a whole species |
| `fast` | Containment pre-filter, `--min-identity 0.5`, batch + parallel shards. | Quick triage / large query panels |

```bash
# AMR surveillance: 1,000+ hits per gene across a species, batch-accelerated
dragon search -i saureus/b1 --shard saureus/b2 ... -q card_amr.fa --preset amr -j 16

# Distant homolog search (cross-genus)
dragon search -i gtdb/b1 --shard gtdb/b2 ... -q gene.fa --preset cross-species
```

### Output formats

- **PAF** — minimap2-compatible pairwise alignment.
- **BLAST6** — BLAST-tabular `outfmt 6`.
- **summary** — per-species prevalence + identity distribution (surveillance-ready).
- **gfa** — graph-context unitigs around each hit (for mobile-element analysis).

## Cloud-native access (Zarr over S3)

`dragon export-zarr` rewrites an index as a [Zarr v3](https://zarr.dev) store — chunked and Zstd-compressed — so it can be **queried directly from object storage without downloading the whole index**. The FM-index, colour bitmaps, and unitig text become chunked arrays; a query fetches only the chunks it touches via HTTP range requests, decompresses them on the fly (Zstd magic-byte sniffing, so mixed-codec stores work), and caches them per-query.

```bash
# 1. Export a binary index to a Zarr store
dragon export-zarr -i saureus/b1 -o saureus_b1.zarr

# 2. Upload to your own bucket (public-read or credentialed)
aws s3 cp --recursive saureus_b1.zarr s3://dragon-zarr/saureus/b1.zarr

# 3a. Query straight from S3 over HTTPS — nothing downloaded up front
dragon search-zarr \
    --zarr https://dragon-zarr.s3.eu-west-2.amazonaws.com/saureus/b1.zarr \
    -q queries.fa -o hits.tsv

# 3b. ...or with an s3:// URI / local path
dragon search-zarr -z s3://dragon-zarr/saureus/b1.zarr -q queries.fa
```

The store is also readable by any Zarr-aware tool (zarr-python, xarray, s3fs) — see `scripts/zarr_demo.py`. This is what lets a laptop query a multi-terabyte database it could never hold locally: only the touched, compressed chunks cross the wire.

> `search-zarr` reports matching text positions + unitig IDs (cloud-native seed discovery). For full base-level alignment, run `dragon search` against a binary index (local or downloaded).

## Pre-built indices

We are publishing four single-/multi-species indices plus a full AllTheBacteria index (the same collection LexicMap ships) under the public-read bucket **`s3://dragon-zarr/`** (eu-west-2, `--no-sign-request`). Each is offered two ways: **(a)** a binary index to download for local/HPC `dragon search`, and **(b)** a Zarr store to query in place with `dragon search-zarr`.

| Database | Genomes | Scope | Shards | Index status | Location |
| --- | --- | --- | --- | --- | --- |
| **S. aureus** | 104,323 | single species | 7 × 16K | ✓ built | `s3://dragon-zarr/saureus/` |
| **K. pneumoniae** | 57,077 | single species | 2 | ✓ built | `s3://dragon-zarr/kpneumoniae/` |
| **E. coli** | 315,066 | single species | building | downloading/indexing | `s3://dragon-zarr/ecoli/` |
| **GTDB r220** | 113,106 | cross-species (≥1 representative / species) | 6 × 18.8K | ✓ built | `s3://dragon-zarr/gtdb/` |
| **AllTheBacteria v2** | ~2,290,000 | all bacteria (LexicMap parity) | sharded (in progress) | building | `s3://dragon-zarr/atb/` |

S3 publication is rolling out per database — a 16,000-genome S. aureus demo shard is already public at `s3://dragon-zarr/saureus/b1` (see the demo at the top of this README). Until a given store is uploaded, you can build any of these locally:

```bash
# Build any species set from a FASTA directory (.fa / .fa.gz both supported)
dragon index -i ./genomes/ -o ./my_index/ -k 31 -j 16 --auto   # --auto shards by RAM

# Or fetch a pre-built index / source genomes
dragon download -d gtdb-r220          -o ./gtdb/      # pre-built index
dragon download -d allthebacteria     -o ./atb_fasta/ # source genomes, then `dragon index`
```

Index size scales with sequence **diversity**, not raw genome count: single-species sets (high redundancy) stay compact and laptop-friendly (S. aureus FM-index ~400 MB for 8K genomes), whereas pan-kingdom sets like GTDB/ATB grow with the novel sequence each clade contributes and suit workstation/HPC RAM.

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
cargo test --lib                # 163 unit tests
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
│   ├── main.rs              CLI entry point (12 subcommands)
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
│   ├── make_demo_data.py    Generate the tiny self-contained demo dataset
│   ├── zarr_quickstart.sh   One-command index → export-zarr → search-zarr demo
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
