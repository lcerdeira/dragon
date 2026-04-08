# Benchmark methodology

## Tools compared

| Tool | Type | Version | Why included |
|------|------|---------|--------------|
| **Dragon** | Graph + FM-index | 0.1.0 | Our tool |
| **LexicMap** | LexicHash probes | latest | Direct predecessor (Nature Biotech 2025) |
| **Minimap2** | Minimiser + chain | 2.28 | Gold standard for long-read alignment |
| **BLASTn** | Word + extend | 2.15.0 | Gold standard for sensitivity |
| **MMseqs2** | k-mer prefilter | 15 | Fast protein/nucleotide search |
| **COBS** | Bit-sliced signatures | latest | k-mer containment search |
| **sourmash** | FracMinHash sketches | 4.8 | Containment estimation |
| **skani** | Sparse chaining ANI | 0.2 | Fast ANI estimation |

## Accuracy metrics

All metrics computed against ground truth from simulation:

- **Sensitivity (recall)**: TP / (TP + FN) &mdash; fraction of true hits found
- **Precision**: TP / (TP + FP) &mdash; fraction of reported hits that are correct
- **F1 score**: harmonic mean of precision and recall
- **Alignment identity error**: |reported_identity - true_identity|

Metrics are stratified by:
- Divergence level (0-15%)
- Query length (short/medium/long)

## Resource metrics

- **Peak RAM**: maximum resident set size via `getrusage()` or `/usr/bin/time -v`
- **Wall-clock time**: elapsed time per query
- **CPU time**: user + system time
- **Index size on disk**: total bytes of all index files
- **Index construction time**: wall-clock time for `dragon index`

## Scalability metrics

- RAM and time vs number of indexed genomes (100, 1K, 10K, 100K, 1M, 2M)
- RAM and time vs batch size (1, 10, 100, 1000 queries)
- Performance on HDD vs SSD

## Read simulation

### Gene-level queries

```python
# Extract random genes
python3 benchmark/simulate/extract_genes.py \
  --genome-dir data/genomes/ \
  --output queries/genes.fa \
  --num-genes 1000 --min-length 500 --max-length 5000

# Introduce mutations
python3 benchmark/simulate/mutate_sequences.py \
  --input queries/genes.fa \
  --output queries/genes_div0.05.fa \
  --divergence 0.05
```

### Long reads (Badread)

```bash
bash benchmark/simulate/run_badread.sh data/genomes/ queries/long_reads.fa
```

## Statistical analysis

- Each measurement repeated 3 times; median reported
- Error bars show min/max across replicates
- Statistical significance tested with Wilcoxon signed-rank test where applicable
