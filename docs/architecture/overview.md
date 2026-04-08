# Architecture overview

## The redundancy problem

Millions of prokaryotic genomes share the vast majority of their sequence content. Within a single species, genomes typically share >95% average nucleotide identity (ANI). Existing tools like LexicMap, Minimap2, and BLASTn index each genome independently, duplicating shared content millions of times.

Dragon solves this by **storing shared sequence once** using a coloured compacted de Bruijn graph, then indexing the unique content with a compressed FM-index.

## Pipeline overview

Dragon has two phases: **offline index construction** (expensive, run once) and **online query** (fast, low RAM, run many times).

### Index construction

```
FASTA genomes
    |
    v
GGCAT (coloured compacted de Bruijn graph)
    |
    +---> Unitig sequences (2-bit encoded)
    |         |
    |         v
    |     Run-length FM-index (r-index)
    |
    +---> Colour index (Roaring Bitmaps: unitig -> genome set)
    |
    +---> Genome path index (genome -> unitig traversal order)
    |
    v
On-disk index (memory-mapped at query time)
```

### Query pipeline

```
Query FASTA
    |
    v
Stage 1: FM-index backward search
    |     (variable-length seed matching)
    v
Stage 2: Candidate genome filtering
    |     (colour-based voting)
    v
Stage 3: Graph-aware colinear chaining
    |     (Fenwick tree DP)
    v
Stage 4: Banded wavefront alignment
    |     (WFA along genome paths)
    v
Alignment results (PAF / BLAST-tabular)
```

## Why this is efficient

| Component | What it compresses | Compression factor |
|-----------|-------------------|-------------------|
| De Bruijn graph | Shared sequence across genomes | ~2,000x (10 Tbp -> 5 Gbp) |
| Run-length FM-index | Repetitive BWT runs | 10-100x (r/n ~ 0.01-0.1) |
| Roaring Bitmaps | Clustered genome ID sets | ~10x vs raw bitvectors |
| Delta-coded paths | Similar genome traversals | 5-10x within species clusters |

**Net result**: ~50x total disk reduction vs LexicMap (100 GB vs 5.46 TB for 2.34M genomes).

## Module map

```
src/
+-- index/          Index construction
|   +-- dbg.rs        GGCAT integration, de Bruijn graph
|   +-- unitig.rs     Unitig parsing, 2-bit encoding
|   +-- color.rs      Roaring Bitmap colour index
|   +-- fm.rs         FM-index construction and search
|   +-- paths.rs      Genome path index
|   +-- metadata.rs   Index metadata
|
+-- query/          Query pipeline
|   +-- seed.rs       FM-index backward search
|   +-- candidate.rs  Colour-based genome voting
|   +-- chain.rs      Colinear chaining (Fenwick tree DP)
|   +-- extract.rs    Reference path extraction
|   +-- align.rs      Banded wavefront alignment
|
+-- ds/             Data structures
|   +-- elias_fano.rs Position-to-unitig mapping
|   +-- fenwick.rs    Binary indexed tree for chaining
|   +-- varint.rs     Variable-length integer coding
|
+-- io/             Input/output
|   +-- fasta.rs      FASTA/FASTQ parser
|   +-- paf.rs        PAF output writer
|   +-- blast.rs      BLAST-tabular output writer
|
+-- util/           Utilities
    +-- dna.rs        2-bit DNA encoding, reverse complement
    +-- mmap.rs       Memory-mapped file helpers
    +-- progress.rs   Progress bars
```
