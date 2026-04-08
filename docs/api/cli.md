# CLI reference

## `dragon index`

Build a Dragon index from a directory of genome FASTA files.

```
dragon index [OPTIONS] --input <DIR> --output <DIR>
```

### Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--input` | `-i` | *required* | Directory containing genome FASTA files |
| `--output` | `-o` | *required* | Output directory for the index |
| `--kmer-size` | `-k` | `31` | K-mer size for the de Bruijn graph |
| `--threads` | `-j` | `4` | Number of threads |

### Examples

```bash
# Basic usage
dragon index -i genomes/ -o my_index/

# Custom k-mer size for higher sensitivity
dragon index -i genomes/ -o my_index/ -k 21

# Use all available cores
dragon index -i genomes/ -o my_index/ -j $(nproc)
```

---

## `dragon search`

Search query sequences against a Dragon index.

```
dragon search [OPTIONS] --index <DIR> --query <FILE>
```

### Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--index` | `-i` | *required* | Path to Dragon index directory |
| `--query` | `-q` | *required* | Query FASTA/FASTQ file |
| `--output` | `-o` | `-` (stdout) | Output file |
| `--format` | `-f` | `paf` | Output format: `paf` or `blast6` |
| `--threads` | `-j` | `4` | Number of threads |
| `--max-ram` | | `4.0` | Maximum RAM in GB |
| `--min-seed-len` | | `15` | Minimum seed match length |
| `--max-seed-freq` | | `10000` | Skip seeds exceeding this frequency |
| `--min-chain-score` | | `50` | Minimum chain score to report |
| `--max-target-seqs` | | `100` | Maximum targets per query |

### Examples

```bash
# Basic search
dragon search -i my_index/ -q query.fa -o results.paf

# BLAST-tabular output
dragon search -i my_index/ -q query.fa -o results.tsv -f blast6

# High-sensitivity mode
dragon search -i my_index/ -q query.fa --min-seed-len 12 --max-seed-freq 50000

# Low-memory mode
dragon search -i my_index/ -q query.fa --max-ram 2.0 -j 2

# Pipe to downstream tools
dragon search -i my_index/ -q query.fa | awk '$12 >= 30' > filtered.paf
```

---

## `dragon info`

Display information about a Dragon index.

```
dragon info --index <DIR>
```

### Options

| Option | Short | Description |
|--------|-------|-------------|
| `--index` | `-i` | Path to Dragon index directory |

### Example output

```
Dragon Index Information
========================
Version:         0.1.0
K-mer size:      31
Genomes:         85205
Unitigs:         4523000
Total bases:     4850000000
Index size:      14.73 GB
```

---

## Environment variables

| Variable | Description |
|----------|-------------|
| `RUST_LOG` | Logging level: `error`, `warn`, `info`, `debug`, `trace` |

```bash
# Enable verbose logging
RUST_LOG=debug dragon search -i my_index/ -q query.fa

# Show only warnings and errors
RUST_LOG=warn dragon search -i my_index/ -q query.fa
```
