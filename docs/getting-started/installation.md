# Installation

## Requirements

- **Rust** 1.75 or later
- **Operating system**: Linux, macOS, or Windows (WSL)
- **RAM**: 4 GB minimum for querying; 8+ GB recommended for index construction
- **Disk**: depends on database size (see [Performance tuning](../user-guide/performance-tuning.md))

## From source (recommended)

### 1. Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### 2. Clone and build

```bash
git clone https://github.com/dragon-aligner/dragon.git
cd dragon
cargo build --release
```

The compiled binary is at `target/release/dragon`.

### 3. (Optional) Install system-wide

```bash
cargo install --path .
# or manually:
cp target/release/dragon /usr/local/bin/
```

### 4. Verify installation

```bash
dragon --version
dragon --help
```

## Optional dependencies

### GGCAT (recommended for large databases)

[GGCAT](https://github.com/algbio/ggcat) provides optimised coloured compacted de Bruijn graph construction. Dragon will use it automatically if found in `PATH`.

```bash
# Install via cargo
cargo install ggcat

# Or build from source
git clone https://github.com/algbio/ggcat.git
cd ggcat
cargo build --release
```

Without GGCAT, Dragon uses a built-in graph builder that works well for datasets up to ~10,000 genomes.

### Benchmark dependencies

To run the benchmark pipeline:

```bash
pip install snakemake matplotlib seaborn pandas numpy
```

## Troubleshooting

### `sux` crate build failure

On Rust 1.94+, the `sux` crate may fail due to `common_traits` ambiguity. Dragon does not depend on `sux` — it uses an internal Elias-Fano implementation. If you encounter this error, ensure your `Cargo.toml` does not list `sux` as a dependency.

### Memory issues during index construction

For very large databases (>100K genomes), index construction may require significant RAM for suffix array construction. Consider:

1. Using a machine with 32+ GB RAM for the one-time index build
2. Distributing the pre-built index to query machines
3. Using `--kmer-size 21` to reduce index memory (at some sensitivity cost)
