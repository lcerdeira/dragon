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
|------|-------------|
| `fm_index.bin` | Run-length FM-index over unitig sequences |
| `colors.drgn` | Roaring Bitmap colour index (unitig-to-genome mapping) |
| `paths.bin` | Genome path index (genome-to-unitig traversal) |
| `unitigs.fa` | Unitig sequences from the de Bruijn graph |
| `metadata.json` | Index statistics (genome count, k-mer size, etc.) |

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
