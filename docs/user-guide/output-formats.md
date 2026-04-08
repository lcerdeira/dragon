# Output formats

## PAF (Pairwise Alignment Format)

Default output format. Tab-separated with 12 mandatory columns:

| Column | Description |
|--------|-------------|
| 1 | Query name |
| 2 | Query length |
| 3 | Query start (0-based) |
| 4 | Query end |
| 5 | Strand (`+` or `-`) |
| 6 | Target genome name |
| 7 | Target genome length |
| 8 | Target start (0-based) |
| 9 | Target end |
| 10 | Number of matching bases |
| 11 | Alignment block length |
| 12 | Mapping quality (0-60) |

### Optional tags

| Tag | Description |
|-----|-------------|
| `AS:i:<N>` | Chain alignment score |
| `cs:f:<F>` | Query coverage fraction |

### Example

```
gene_001  1500  10  1490  +  genome_042  4800000  123456  124946  1450  1490  60  AS:i:2900  cs:f:0.9867
```

## BLAST tabular (outfmt 6)

Use `--format blast6`. Tab-separated with 12 columns:

| Column | Description |
|--------|-------------|
| 1 | Query ID |
| 2 | Subject ID |
| 3 | % identity |
| 4 | Alignment length |
| 5 | Mismatches |
| 6 | Gap opens |
| 7 | Query start (1-based) |
| 8 | Query end |
| 9 | Subject start (1-based) |
| 10 | Subject end |
| 11 | E-value |
| 12 | Bit score |

### Example

```
gene_001  genome_042  96.67  1490  48  2  11  1490  123457  124946  0.00e+00  2900.0
```

## Parsing output

### Extract top hits per query

```bash
# PAF: best hit per query (highest mapping quality)
sort -k1,1 -k12,12rn results.paf | awk '!seen[$1]++' > best_hits.paf

# BLAST: best hit per query (highest bit score)
sort -k1,1 -k12,12rn results.tsv | awk '!seen[$1]++' > best_hits.tsv
```

### Filter by identity

```bash
# PAF: filter by >90% identity (matches/alignment_length)
awk -F'\t' '$10/$11 > 0.90' results.paf > filtered.paf

# BLAST: filter directly on column 3
awk -F'\t' '$3 > 90' results.tsv > filtered.tsv
```

### Count hits per genome

```bash
cut -f6 results.paf | sort | uniq -c | sort -rn | head -20
```
