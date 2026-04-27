# Dragon Documentation

> A cloud-native, signal-aware aligner for surveillance-scale microbial genomics.

Dragon aligns query sequences (genes, plasmids, long/short reads, raw nanopore current) against millions of prokaryotic genomes while using dramatically less disk and RAM than existing tools. It exploits redundancy among related genomes through a coloured compacted de Bruijn graph, an FM-index over concatenated unitigs, ML-weighted graph-aware chaining, and a streaming on-disk format that mmaps the index in O(1).

## Key features

- **~50× less disk** than LexicMap (~100 GB vs 5.46 TB for 2.34 M genomes).
- **<4 GB query RAM** at million-genome scale; `--profile laptop` further restricts use to consumer hardware.
- **Multi-shard search** (`--shard`) for indices split across files or quotas.
- **Cloud-native Zarr v3 backend** (`dragon export-zarr` / `dragon search-zarr`) — chunked + Zstd-compressed; reads run against `s3://` or `gs://` directly via `zarr-python`.
- **Mmap-friendly `paths.bin v2`** — O(1) cold-load, per-genome lazy decoding from a fixed offset table.
- **Raw nanopore signal search** (`dragon signal-index` / `dragon signal-search`) — pore-model–driven discretisation indexed by the same FM-index machinery, no basecalling required.
- **ML-weighted seed scoring** — logistic regression over six anchor features; pure Rust inference.
- **Surveillance-ready summaries** (`dragon summarize`, `--format summary`) — per-species prevalence + identity tables built into the CLI.
- **Incremental updates** (`dragon update` / `dragon compact`) — overlay new genomes without a full rebuild.
- **Variable-length seeds** via FM-index backward search.
- **Outputs** in PAF, BLAST-tabular, surveillance summary, and graph-context GFA formats.

A 16,000-genome demo index is hosted at `s3://dragon-zarr/saureus/b1/` (eu-west-2, public-read). Anyone can read it with no AWS credentials:

```bash
pip install 'zarr>=3.0' s3fs numcodecs
python scripts/zarr_demo.py s3://dragon-zarr/saureus/b1
```

## Contents

```{toctree}
:maxdepth: 2
:caption: Getting Started

getting-started/installation
getting-started/quickstart
getting-started/tutorial
```

```{toctree}
:maxdepth: 2
:caption: User Guide

user-guide/indexing
user-guide/searching
user-guide/output-formats
user-guide/performance-tuning
```

```{toctree}
:maxdepth: 2
:caption: Architecture

architecture/overview
architecture/de-bruijn-graph
architecture/fm-index
architecture/chaining
architecture/data-structures
```

```{toctree}
:maxdepth: 2
:caption: Benchmark

benchmark/datasets
benchmark/methodology
benchmark/results
benchmark/reproducing
```

```{toctree}
:maxdepth: 2
:caption: API Reference

api/cli
api/modules
```

## Citation

> Cerdeira, L. (2026). Dragon: a cloud-native, signal-aware aligner for surveillance-scale microbial genomics. *In preparation.*

## Licence

Dragon is released under the [MIT Licence](https://opensource.org/licenses/MIT).
