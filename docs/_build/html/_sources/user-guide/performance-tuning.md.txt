# Performance tuning

## Hardware recommendations

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| RAM (query) | 4 GB | 8 GB | Dragon stays <4 GB; extra for OS |
| RAM (index build) | 8 GB | 32-64 GB | SA-IS construction |
| Storage | SSD | NVMe SSD | Memory-mapped access; HDD 5-10x slower |
| CPU | 2 cores | 8+ cores | Parallel queries via rayon |

## SSD vs HDD

Dragon's index is memory-mapped, so I/O speed directly affects performance:

| Storage | Single gene query | Batch (1000 genes) |
|---------|-------------------|-------------------|
| NVMe SSD | ~1 second | ~5 minutes |
| SATA SSD | ~3 seconds | ~15 minutes |
| HDD | ~15 seconds | ~2 hours |

**Recommendation**: use SSD whenever possible. If using HDD, increase `--max-ram` to keep more index pages resident.

## Tuning parameters

### For maximum sensitivity

```bash
dragon search \
  --min-seed-len 12 \
  --max-seed-freq 50000 \
  --min-chain-score 20 \
  --max-target-seqs 500
```

### For maximum speed

```bash
dragon search \
  --min-seed-len 21 \
  --max-seed-freq 1000 \
  --min-chain-score 100 \
  --max-target-seqs 10 \
  --threads 16
```

### For low-memory machines (4 GB total RAM)

```bash
dragon search \
  --max-ram 2.0 \
  --threads 2
```

## Scaling guidelines

### Number of genomes vs resources

| Genomes | Index disk | Query RAM | Build time |
|---------|-----------|-----------|------------|
| 100 | 200 MB | <500 MB | 5 seconds |
| 1,000 | 1 GB | <500 MB | 30 seconds |
| 10,000 | 5 GB | 1 GB | 10 minutes |
| 100,000 | 20 GB | 2 GB | 2 hours |
| 1,000,000 | 60 GB | 3 GB | 8 hours |
| 2,340,000 | 100 GB | 3.5 GB | 12 hours |

### Query length vs performance

| Query length | Seeds found | Chaining time | Total time |
|-------------|-------------|---------------|------------|
| 150 bp | ~5 | <1 ms | ~0.1 s |
| 1,000 bp | ~30 | ~5 ms | ~0.5 s |
| 10,000 bp | ~300 | ~50 ms | ~2 s |
| 100,000 bp | ~3,000 | ~500 ms | ~15 s |

## Distributing pre-built indices

Since index construction is expensive but only done once, consider:

1. **Build once on a server** with sufficient RAM
2. **Distribute the index** to query machines (e.g., via rsync, S3, or shared filesystem)
3. **Query on laptops** with just 4 GB RAM

```bash
# On build server
dragon index -i all_genomes/ -o dragon_index/ -j 32

# Transfer to query machine
rsync -avP dragon_index/ user@laptop:~/dragon_index/

# On laptop
dragon search -i ~/dragon_index/ -q my_query.fa -o results.paf
```
