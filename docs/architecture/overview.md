# Architecture overview

## The redundancy problem

Millions of prokaryotic genomes share the vast majority of their sequence content. Within a single species, genomes typically share >95% average nucleotide identity (ANI). Existing tools like LexicMap, Minimap2, and BLASTn index each genome independently, duplicating shared content millions of times.

Dragon solves this by **storing shared sequence once** using a coloured compacted de Bruijn graph, then indexing the unique content with a compressed FM-index.

## Pipeline overview

Dragon has two phases: **offline index construction** (expensive, run once) and **online query** (fast, low RAM, run many times).

### Index construction

```text
FASTA genomes
    |
    v
GGCAT (coloured compacted de Bruijn graph)
    |
    +---> unitigs.fa (2-bit encoded unitigs)
    |         |
    |         v
    |     fm_index.bin (concatenated text + suffix array)
    |
    +---> colors.drgn      (Roaring Bitmaps: unitig -> genome set)
    +---> paths.bin v2     (genome -> unitig traversal, mmap'd offset table)
    +---> specificity.drgn (per-genome private-unitig sets)
    +---> metadata.json
    |
    v
On-disk Dragon index (memory-mapped at query time)
    |
    | dragon export-zarr
    v
Zarr v3 store (chunked + Zstd, served from local disk or s3:// / gs://)
```

### Query pipeline

```text
Query FASTA
    |
    v
Stage 1: FM-index backward search
    |     (variable-length seed matching)
    v
Stage 2: Candidate genome filtering
    |     (Roaring-bitmap colour voting)
    v
Stage 3: ML-weighted graph-aware colinear chaining
    |     (logistic-regression seed scoring + Fenwick / O(h^2) DP)
    v
Stage 4: Containment ranking (top-N candidates by total matched bases)
    |
    v
Stage 5: Banded wavefront alignment along genome paths
    |
    v
PAF / BLAST6 / surveillance summary / GFA output
```

Multi-shard search (`--shard`) drives this pipeline once per shard, then merges results with per-genome deduplication so quota- or RAM-split indices behave like a single logical database.

## Why this is efficient

| Component | What it compresses | Compression factor |
|-----------|-------------------|-------------------|
| De Bruijn graph | Shared sequence across genomes | ~2,000x (10 Tbp -> 5 Gbp) |
| Run-length FM-index | Repetitive BWT runs | 10-100x (r/n ~ 0.01-0.1) |
| Roaring Bitmaps | Clustered genome ID sets | ~10x vs raw bitvectors |
| Delta-coded paths | Similar genome traversals | 5-10x within species clusters |

**Net result**: ~50x total disk reduction vs LexicMap (100 GB vs 5.46 TB for 2.34M genomes).

## Module map

```text
src/
+-- index/                      Index construction
|   +-- dbg.rs                  GGCAT integration (with internal-builder fallback)
|   +-- unitig.rs               Unitig parsing, 2-bit encoding
|   +-- color.rs                Roaring-bitmap colour index
|   +-- ggcat_colors.rs         GGCAT binary colormap -> colors.drgn (no TSV)
|   +-- fm.rs                   FM-index construction and search
|   +-- paths.rs                Genome path index (legacy bincode loader)
|   +-- paths_v2.rs             Mmap-friendly v2 format (default for new builds)
|   +-- specificity.rs          Per-genome private-unitig sets
|   +-- auto_batch.rs           Auto-split large collections into overlay batches
|   +-- update.rs               Incremental overlay addition (dragon update / compact)
|   +-- zarr_backend.rs         Zarr v3 export + ZarrFmIndex / ZarrColorIndex readers
|   +-- metadata.rs             Index metadata (JSON)
|
+-- query/                      Query pipeline
|   +-- seed.rs                 FM-index backward search (variable-length)
|   +-- candidate.rs            Colour-based genome voting
|   +-- chain.rs                ML-weighted graph-aware chaining + containment ranking
|   +-- containment.rs          K-mer containment scoring
|   +-- direct_align.rs         Direct alignment to top candidate genomes
|   +-- align.rs                Banded wavefront alignment along genome paths
|   +-- mod.rs                  Multi-shard orchestration (search_multi_index)
|
+-- signal/                     Raw nanopore current search
|   +-- index.rs                Pore-model-driven discretisation + signal FM-index
|   +-- search.rs               Backward search over signal k-mers
|
+-- ds/                         Data structures
|   +-- elias_fano.rs           CumulativeLengthIndex (position -> unitig)
|   +-- fenwick.rs              Binary indexed tree for O(h log h) chaining
|   +-- varint.rs               LEB128 codecs (also used by paths_v2)
|
+-- io/                         Input / output
|   +-- fasta.rs                FASTA / FASTQ parser
|   +-- paf.rs                  PAF writer
|   +-- blast.rs                BLAST-tabular writer
|   +-- gfa.rs                  Graph-context GFA writer
|   +-- summary.rs              Surveillance summary writer
|
+-- util/                       Utilities
    +-- dna.rs                  2-bit DNA encoding, reverse complement
    +-- mmap.rs                 Bincode + memory-mapped helpers
    +-- colorspace.rs           SOLiD-style 2-base colour-space encoder/decoder
    +-- progress.rs             Progress bars
```
