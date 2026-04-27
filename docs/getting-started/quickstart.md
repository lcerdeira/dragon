# Quick start

This guide walks you through indexing a small set of genomes and searching for a query sequence.

## Step 1: Prepare genome files

Organise your reference genomes as individual FASTA files in a directory:

```
genomes/
  genome_001.fasta
  genome_002.fasta
  genome_003.fasta
  ...
```

Supported extensions: `.fa`, `.fasta`, `.fna`, `.fsa`

## Step 2: Build the index

```bash
dragon index \
  --input genomes/ \
  --output my_index/ \
  --kmer-size 31 \
  --threads 8
```

This creates the following files in `my_index/`:

| File | Description |
| --- | --- |
| `fm_index.bin` | FM-index over concatenated unitig sequences |
| `colors.drgn` | Roaring-bitmap colour index (unitig → genome mapping) |
| `paths.bin` | Genome path index (mmap-friendly v2 format) |
| `specificity.drgn` | Per-genome private-unitig sets |
| `unitigs.fa` | Unitig sequences from the de Bruijn graph (optional after build) |
| `metadata.json` | Index statistics (genome count, k-mer size, total bases) |

## Step 3: Search

```bash
dragon search \
  --index my_index/ \
  --query query_genes.fasta \
  --output results.paf \
  --threads 8
```

## Step 4: Inspect results

```bash
# View PAF output
head results.paf

# Count hits per query
cut -f1 results.paf | sort | uniq -c | sort -rn | head

# View index statistics
dragon info --index my_index/
```

## Example output

PAF format (tab-separated):

```
gene_001  1500  10  1490  +  genome_042  4800000  123456  124946  1450  1490  60  AS:i:2900
gene_001  1500  10  1490  +  genome_108  5100000  234567  236057  1430  1490  55  AS:i:2860
```

Columns: query name, query length, query start, query end, strand, target name, target length, target start, target end, matches, alignment length, mapping quality, tags.

## Step 5 (optional): Multi-shard search

If your collection is too large for a single index, build several shards and search them as one:

```bash
dragon search \
  --index shard_a/ \
  --shard shard_b/ \
  --shard shard_c/ \
  --query query_genes.fasta \
  --output results.paf
```

Each shard is loaded in turn (memory-bounded) and results are merged with per-genome deduplication.

## Step 6 (optional): Cloud-native deployment

Export an index as a Zarr v3 store for direct reading from S3 / GCS:

```bash
dragon export-zarr -i my_index/ -o my_index.zarr/
aws s3 sync my_index.zarr/ s3://your-bucket/my_index/

# Anywhere with internet (no AWS creds needed for public buckets):
pip install 'zarr>=3.0' s3fs numcodecs
python scripts/zarr_demo.py s3://your-bucket/my_index
```

A pre-built 16,000-genome demo lives at `s3://dragon-zarr/saureus/b1/` (eu-west-2, public-read). See [Architecture overview](../architecture/overview.md) for details.

## Step 7 (optional): Surveillance summary

For AMR-gene panels and similar epidemiological queries, ask for a per-species summary instead of raw PAF:

```bash
dragon search -i my_index/ -q amr_genes.fa --format summary > prevalence.tsv
```

Or post-process an existing PAF:

```bash
dragon summarize --input results.paf --format tsv > prevalence.tsv
```
