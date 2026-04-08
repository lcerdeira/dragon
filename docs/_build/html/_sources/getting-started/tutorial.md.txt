# Tutorial: AMR gene search

This tutorial demonstrates searching antimicrobial resistance (AMR) genes against a collection of bacterial genomes.

## Scenario

You have:
- A database of 500 *E. coli* genomes
- A set of AMR gene sequences from the CARD database
- Goal: identify which genomes carry which resistance genes

## Step 1: Download test data

```bash
# Create working directory
mkdir dragon_tutorial && cd dragon_tutorial

# Download example E. coli genomes (first 10 for this tutorial)
# In practice, download from NCBI/RefSeq
mkdir genomes
for i in $(seq 1 10); do
  # Replace with actual genome download commands
  echo ">genome_${i}" > genomes/genome_${i}.fa
  # ... download genome FASTA ...
done

# Download AMR genes from CARD
# https://card.mcmaster.ca/download
wget -O card_genes.fasta https://card.mcmaster.ca/latest/data
```

## Step 2: Build index

```bash
dragon index \
  --input genomes/ \
  --output ecoli_index/ \
  --kmer-size 31 \
  --threads 4

# Check the index
dragon info --index ecoli_index/
```

Expected output:
```
Dragon Index Information
========================
Version:         0.1.0
K-mer size:      31
Genomes:         10
Unitigs:         45230
Total bases:     48500000
Index size:      0.05 GB
```

## Step 3: Search AMR genes

```bash
dragon search \
  --index ecoli_index/ \
  --query card_genes.fasta \
  --output amr_hits.paf \
  --format paf \
  --threads 4 \
  --min-chain-score 100

# Count hits
echo "Total hits: $(wc -l < amr_hits.paf)"

# Top genomes with most AMR genes
cut -f6 amr_hits.paf | sort | uniq -c | sort -rn | head
```

## Step 4: BLAST-tabular output

```bash
dragon search \
  --index ecoli_index/ \
  --query card_genes.fasta \
  --output amr_hits.tsv \
  --format blast6

# Filter by identity > 90%
awk -F'\t' '$3 > 90' amr_hits.tsv | head
```

## Step 5: Batch analysis

Dragon handles batch queries efficiently via parallel processing:

```bash
# Search 1000 AMR genes with 8 threads
dragon search \
  --index ecoli_index/ \
  --query card_all_1003_genes.fasta \
  --output all_amr_hits.paf \
  --threads 8 \
  --max-ram 4.0
```

## Performance notes

- **Small datasets (<100 genomes)**: index builds in seconds, queries in milliseconds
- **Medium datasets (~85K genomes)**: index builds in ~1 hour, queries in seconds
- **Large datasets (>1M genomes)**: index builds in hours (use GGCAT), queries in minutes
- **RAM**: stays below 4 GB at query time regardless of database size
