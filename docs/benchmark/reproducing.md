# Reproducing benchmarks

## Quick start (synthetic data)

Run the complete benchmark pipeline with synthetic test data in under 2 minutes:

```bash
# 1. Build Dragon
cargo build --release

# 2. Generate synthetic test data (20 genomes, 50 queries)
python3 benchmark/setup_test_data.py

# 3. Run benchmarks
python3 benchmark/run_benchmark.py
```

Results are saved to:
- `benchmark/results/metrics/` &mdash; accuracy and resource TSV files
- `manuscript/figures/` &mdash; regenerated PDF and PNG figures

## Full benchmark (Snakemake)

For the full benchmark with real data and all tools:

### Prerequisites

```bash
# Install benchmark dependencies
pip install snakemake matplotlib seaborn pandas numpy

# Install comparison tools
conda install -c bioconda lexicmap minimap2 blast mmseqs2 cobs sourmash skani

# Install read simulators
conda install -c bioconda badread art
```

### Running

```bash
cd benchmark

# Dry run (show what will be executed)
snakemake --cores 8 -n

# Full run
snakemake --cores 8 --use-conda

# Regenerate figures only
snakemake --cores 1 -R plot_figures
```

### Configuration

Edit `benchmark/config.yaml` to customise:

```yaml
# Dataset paths
datasets:
  tier1:
    genome_dir: "data/tier1_genomes"

# Simulation parameters
simulation:
  gene_queries:
    num_genes: 1000
    divergence_levels: [0.0, 0.01, 0.03, 0.05, 0.10, 0.15]

# Hardware
hardware:
  default_threads: 8
```

## Pipeline structure

```
benchmark/
+-- Snakefile                    Workflow orchestration
+-- config.yaml                  Configuration
+-- setup_test_data.py           Synthetic data generator
+-- run_benchmark.py             Standalone benchmark runner
+-- simulate/
|   +-- extract_genes.py         Extract random genes from genomes
|   +-- mutate_sequences.py      Introduce controlled mutations
|   +-- run_badread.sh           Long-read simulation
|   +-- run_art.sh               Short-read simulation
+-- scripts/
|   +-- compute_metrics.py       Sensitivity, precision, F1
|   +-- resource_usage.py        Parse /usr/bin/time output
+-- notebooks/
|   +-- figures.py               Generate all manuscript figures
+-- envs/
|   +-- simulation.yaml          Conda environment
+-- results/                     Output (auto-generated)
    +-- metrics/                 TSV metric files
    +-- figures/                 PDF/PNG plots
    +-- search/                  Raw search output files
```

## Adding a new tool

1. Add tool configuration to `config.yaml` under `tools:`
2. Create a conda environment file in `envs/`
3. The Snakemake pipeline will automatically include it in indexing, searching, and evaluation

```yaml
# Example: adding a new tool
tools:
  my_new_tool:
    binary: "mytool"
    version: "1.0"
    index_cmd: "mytool index {genome_dir} {index_dir}"
    search_cmd: "mytool search {index_dir} {query} > {output}"
```
