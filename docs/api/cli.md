# CLI reference

Dragon ships a single binary with eleven subcommands. Run `dragon <subcommand> --help` for the full option list at any time; this page describes each subcommand and its most-used flags.

| Command | Purpose |
| --- | --- |
| `dragon index` | Build a Dragon index from a directory of FASTA genomes |
| `dragon search` | Align query sequences (single or multi-shard) |
| `dragon info` | Print index metadata |
| `dragon download` | Download genomes (RefSeq, AllTheBacteria) or pre-built indices |
| `dragon update` | Add new genomes as a lightweight overlay |
| `dragon compact` | Merge base + overlays back into one optimised index |
| `dragon summarize` | Per-species prevalence/identity report from PAF output |
| `dragon export-zarr` | Export an index as a Zarr v3 store |
| `dragon search-zarr` | Pattern-search a Zarr-backed index (local or `s3://`) |
| `dragon signal-index` | Build a signal-level index from FASTA via a pore model |
| `dragon signal-search` | Align raw nanopore current signals (TSV/CSV/SLOW5) |

---

## `dragon index`

Build a Dragon index from a directory of genome FASTA files. Uses GGCAT for the colored compacted de Bruijn graph if available, falling back to an internal builder for small datasets.

```text
dragon index [OPTIONS] --input <DIR> --output <DIR>
```

| Option | Short | Default | Description |
| --- | --- | --- | --- |
| `--input` | `-i` | *required* | Directory of genome FASTA files (`.fa`, `.fasta`, `.fna`) |
| `--output` | `-o` | *required* | Output directory for the index |
| `--kmer-size` | `-k` | `31` | K-mer size for the de Bruijn graph |
| `--threads` | `-j` | `4` | Number of threads |
| `--low-memory` | | off | External-memory SA construction (≤8 GB RAM by default) |
| `--max-ram` | | `8.0` | RAM budget in GB for `--low-memory` mode |
| `--auto` | | off | Auto-batch large collections into overlays (transparent at query time) |

The index directory contains:

- `fm_index.bin` — concatenated unitig text + suffix array
- `colors.drgn` — Roaring-bitmap colour index per unitig
- `paths.bin` — per-genome unitig path (mmap-friendly v2 format on new builds)
- `specificity.drgn` — per-genome private-unitig sets
- `metadata.json` — version, k-mer size, genome count, total bases
- `unitigs.fa` (optional) — keep for resume/auto-batch; safe to delete after successful build

### Index examples

```bash
# Default: 31-mer, 4 threads, RAM-bounded only by the system
dragon index -i genomes/ -o my_index/

# Low-memory: external-memory SA construction with 8 GB cap
dragon index -i genomes/ -o my_index/ --low-memory --max-ram 8

# Use all cores
dragon index -i genomes/ -o my_index/ -j $(nproc)

# Auto-batch a million-genome collection into overlays
dragon index -i giant_dir/ -o giant_idx/ --auto --max-ram 64
```

### Resume

Index construction is resumable: if `fm_index.bin` and `colors.drgn` already exist in the output directory, Dragon skips Steps 1–4 (GGCAT + FM-index + colours) and resumes from Step 5 (path index). Useful when a job is killed during the long path-building step.

---

## `dragon search`

Search query sequences against a Dragon index. Supports multi-shard search via repeatable `--shard` arguments — each shard is searched independently and results are merged with per-genome deduplication.

```text
dragon search [OPTIONS] --index <DIR> --query <FILE>
```

### Core options

| Option | Short | Default | Description |
| --- | --- | --- | --- |
| `--index` | `-i` | *required* | Path to Dragon index directory |
| `--shard` | | — | Additional shard directory (repeatable) |
| `--query` | `-q` | *required* | Query FASTA/FASTQ file |
| `--output` | `-o` | `-` (stdout) | Output file |
| `--format` | `-f` | `paf` | `paf`, `blast6`, `summary`, or `gfa` |
| `--threads` | `-j` | `4` | Number of threads |
| `--max-ram` | | `4.0` | RAM budget in GB |
| `--profile` | | `workstation` | `laptop` (≤8 GB, 4 threads) or `workstation` |

### Filtering & scoring

| Option | Default | Description |
| --- | --- | --- |
| `--min-seed-len` | `15` | Minimum seed match length |
| `--max-seed-freq` | `10000` | Skip seeds occurring more than this many times |
| `--min-chain-score` | `50` | Minimum chain score to report |
| `--max-target-seqs` | `10` | Maximum hits per query |
| `--min-identity` | `0.7` | Minimum alignment identity (0.0–1.0) |
| `--min-query-coverage` | `0.3` | Minimum query coverage (0.0–1.0) |
| `--min-score-ratio` | `0.1` | Drop hits scoring below `ratio × best_score` |

### ML scoring & training

| Option | Default | Description |
| --- | --- | --- |
| `--no-ml` | off | Disable learned seed scoring (use raw match length) |
| `--ml-weights` | built-in | Path to a custom JSON of 7 scorer weights |
| `--dump-seeds` | — | Dump every seed + features to TSV for ML training |
| `--ground-truth` | — | Ground-truth genome name (with `--dump-seeds`) for labelled training data |
| `--gfa-radius` | `5` | Number of unitig steps around each hit (used with `--format gfa`) |

### Search examples

```bash
# Basic search
dragon search -i my_index/ -q query.fa -o results.paf

# Multi-shard against several species-level batches
dragon search -i saureus_b1/ \
    --shard saureus_b2/ --shard saureus_b3/ --shard kpneumo_b1/ \
    -q amr_genes.fa -o hits.paf

# Surveillance summary instead of PAF
dragon search -i my_index/ -q amr_genes.fa --format summary

# Laptop profile (clamps RAM and threads)
dragon search -i my_index/ -q query.fa --profile laptop

# Pipe through standard PAF tooling
dragon search -i my_index/ -q query.fa | awk '$12 >= 30' > filtered.paf
```

---

## `dragon info`

Display index metadata.

```text
dragon info --index <DIR>
```

### Example output

```text
Dragon Index Information
========================
Version:         0.1.0
K-mer size:      31
Genomes:         32000
Unitigs:         9137000
Total bases:     434072471
Index size:      783.86 GB
```

---

## `dragon download`

Download genomes or a pre-built Dragon index.

```text
dragon download [OPTIONS] --database <NAME> --output <DIR>
```

### Supported databases

| Name | Behaviour |
| --- | --- |
| `gtdb-r220` | Download a pre-built GTDB r220 representative-genomes index |
| `allthebacteria-v2` | Download a pre-built AllTheBacteria v2 index |
| `refseq-bacteria` | Download a pre-built RefSeq bacteria index |
| `allthebacteria` | Download genomes from EBI AllTheBacteria, then build |
| `refseq` | Download all RefSeq bacteria genomes, then build |
| `refseq-representative` | Download only RefSeq representatives, then build |
| `http(s)://...` | Custom URL to a pre-built index tarball |

### Example

```bash
dragon download -d gtdb-r220 -o gtdb_r220_index/
```

For the genome-download-and-build modes, the index construction step honours `--low-memory`, `--kmer-size`, `--threads`, and `--max-ram`.

---

## `dragon update`

Add new genomes as a lightweight overlay without re-running GGCAT / FM-index construction. Queries automatically search both the base index and all overlays.

```text
dragon update --index <DIR> --genomes <DIR> [--kmer-size 31]
```

When overlays exceed ~10 % of the base index, `dragon update` warns and recommends running `dragon compact`.

---

## `dragon compact`

Merge the base index and all overlays back into one optimised index. Run after `dragon update` has accumulated significant overlay growth.

```text
dragon compact --index <DIR> --genomes <DIR> [--kmer-size 31]
```

`--genomes` should point to **all** FASTA files (base + overlays), since `compact` rebuilds from scratch with the merged genome set.

---

## `dragon summarize`

Generate a per-species surveillance summary from PAF output produced by `dragon search`.

```text
dragon summarize --input <PAF> [--output <FILE>] [--format tsv|json]
                 [--index <DIR>] [--total-genomes <N>]
```

The summary contains, per species:

- Prevalence (fraction of database genomes carrying the query)
- Mean / min / max alignment identity
- Number of unique sequence variants

Designed for AMR-gene surveillance and other epidemiological queries.

---

## `dragon export-zarr`

Export a Dragon index as a Zarr v3 store (chunked, Zstd-compressed). The original index is not modified.

```text
dragon export-zarr --index <DIR> --output <DIR>
```

Store layout under `<output>/`:

```text
zarr.json                   root attrs (kmer_size, num_genomes, ...)
text/                       u8 unitig text, 1 MiB chunks, Zstd-3
suffix_array/               u64 SA, 131 072-entry chunks, Zstd-3
unitig_lengths/             u64 per-unitig lengths
colors/offsets/             u64 byte offsets per unitig
colors/bitmaps/             raw RoaringBitmap bytes, 1 MiB chunks, Zstd-3
```

Anyone with `zarr-python` and `s3fs` can open the store from a public bucket without AWS credentials — the on-disk format is the on-cloud format.

---

## `dragon search-zarr`

Pattern-search a Zarr-backed index. Reads only the chunks each query touches, making it cheap over remote object stores (S3, GCS) that expose HTTP range requests.

```text
dragon search-zarr --zarr <PATH | s3://...> --query <FILE> [--output <FILE>]
```

For full alignment use the binary backend (`dragon search`). `search-zarr` is intended as the cloud-native pattern-match demo: it returns matching text positions, the underlying unitig IDs, and the genomes carrying each unitig.

---

## `dragon signal-index`

Build a signal-level index from genome FASTA files by converting expected nanopore current via a pore model and discretising into a finite alphabet.

```text
dragon signal-index --input <DIR> --output <DIR> [OPTIONS]
```

| Option | Default | Description |
| --- | --- | --- |
| `--num-levels` | `16` | Discretisation alphabet size |
| `--threads` | `4` | Threads |
| `--signal-boundaries` | — | Path to learned discretisation boundaries (JSON) |
| `--pore-model` | built-in R10.4.1 | Path to a custom pore model JSON |

---

## `dragon signal-search`

Align raw nanopore current signals against a signal-level index. Inputs are TSV, CSV, or SLOW5 text format with auto-detection.

```text
dragon signal-search --index <DIR> --query <FILE> [OPTIONS]
```

| Option | Default | Description |
| --- | --- | --- |
| `--signal-kmer-size` | `10` | Signal k-mer size for backward search |
| `--min-hits` | `3` | Minimum k-mer hits to report a genome |
| `--max-seed-freq` | `10000` | Skip signal k-mers above this frequency |
| `--max-results` | `50` | Maximum results per query read |
| `--threads` | `4` | Threads |

---

## Environment variables

| Variable | Description |
| --- | --- |
| `RUST_LOG` | Logging level: `error`, `warn`, `info`, `debug`, `trace` |

```bash
RUST_LOG=debug dragon search -i my_index/ -q query.fa
RUST_LOG=warn  dragon search -i my_index/ -q query.fa
```
