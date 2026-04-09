# Dragon Documentation

**Resource-efficient sequence alignment against millions of prokaryotic genomes**

Dragon aligns query sequences (genes, plasmids, long reads) against millions of prokaryotic genomes while using dramatically less disk and RAM than existing tools. It exploits the massive sequence redundancy among related genomes through a coloured compacted de Bruijn graph, a run-length FM-index, and graph-aware colinear chaining.

## Key features

- **50x less disk** than LexicMap (~100 GB vs 5.46 TB for 2.34M genomes)
- **5x less RAM** at query time (<4 GB vs 4-25 GB)
- **Batch-friendly** parallel queries over shared memory-mapped index
- **Variable-length seeds** via FM-index backward search
- **Graph-aware chaining** for accurate alignments across structural variation
- Outputs in **PAF** and **BLAST-tabular** formats

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

> Cerdeira, L. (2026). Dragon: resource-efficient sequence alignment against millions of prokaryotic genomes using graph-based compressed indexing. *In preparation.*

## Licence

Dragon is released under the [MIT Licence](https://opensource.org/licenses/MIT).
