# Benchmark datasets

## Tiered approach

Dragon is benchmarked at three scales to validate scalability from laptop to server:

### Tier 1: Small (development & CI)

| Property | Value |
|----------|-------|
| **Genomes** | 500 complete *E. coli* / *Shigella* from RefSeq |
| **Total sequence** | ~2.5 Gbp |
| **Redundancy** | High (~95% ANI within species) |
| **Index time** | <1 hour |
| **Index size** | ~1.5 GB |
| **Use case** | Unit testing, CI, rapid iteration |

### Tier 2: Medium (validation)

| Property | Value |
|----------|-------|
| **Genomes** | ~85,000 GTDB r220 representative genomes |
| **Total sequence** | ~250 Gbp |
| **Redundancy** | Medium (one genome per species) |
| **Index time** | ~1 hour |
| **Index size** | ~15 GB |
| **Use case** | Sensitivity/accuracy validation |

### Tier 3: Large (full benchmark)

| Property | Value |
|----------|-------|
| **Genomes** | ~2.34M GenBank + RefSeq prokaryotic assemblies |
| **Total sequence** | ~10 Tbp |
| **Redundancy** | Very high (many strains per species) |
| **Index time** | ~12 hours |
| **Index size** | ~100 GB |
| **Use case** | Full-scale comparison with LexicMap |

## Query types

### Gene-level queries (primary)

- 1,000 random subsequences (500-5,000 bp) extracted from 100 diverse genomes
- Controlled mutations at 6 divergence levels: 0%, 1%, 3%, 5%, 10%, 15%
- Mutation model: 70% substitutions, 20% insertions, 10% deletions

### Long reads (Badread)

- Simulated Oxford Nanopore reads from 50 genomes
- Mean length: 5 Kbp
- Identity range: 85-99%
- Chimera rate: 1%

### Challenging scenarios

- **16S rRNA**: highly conserved, extreme seed frequency
- **Plasmids**: 10-200 Kbp, multi-copy
- **AMR genes**: batch of 1,003 genes from CARD database
- **HGT events**: genes implanted from distant species
