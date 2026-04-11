pub mod color;
pub mod dbg;
pub mod fm;
pub mod metadata;
pub mod paths;
pub mod specificity;
pub mod unitig;
pub mod update;

use anyhow::Result;
use std::path::Path;

/// Build a Dragon index from a directory of FASTA genome files.
pub fn build_index(genome_dir: &Path, output_dir: &Path, kmer_size: usize, threads: usize) -> Result<()> {
    build_index_with_options(genome_dir, output_dir, kmer_size, threads, None)
}

/// Build a Dragon index with optional RAM budget for low-memory construction.
///
/// If `max_ram_bytes` is `Some(n)`, the FM-index suffix array is built using
/// external-memory (disk-based) sorting, reducing peak RAM from ~8n bytes to ~n bytes.
/// This enables indexing 2M genomes on hardware with only ~8 GB RAM instead of 64 GB.
pub fn build_index_with_options(
    genome_dir: &Path,
    output_dir: &Path,
    kmer_size: usize,
    threads: usize,
    max_ram_bytes: Option<usize>,
) -> Result<()> {
    log::info!("Building Dragon index from {:?}", genome_dir);
    log::info!("Parameters: k={}, threads={}", kmer_size, threads);
    if let Some(ram) = max_ram_bytes {
        log::info!("Low-memory mode: RAM budget {} GB", ram as f64 / 1_073_741_824.0);
    }

    std::fs::create_dir_all(output_dir)?;

    // Check if Step 1 can be skipped (resume from existing unitigs + colors)
    let unitig_file = output_dir.join("unitigs.fa");
    let color_file = output_dir.join("colors.tsv");
    let dbg_result = if unitig_file.exists() && color_file.exists() {
        log::info!("Step 1/5: SKIPPED — existing unitigs.fa + colors.tsv found (resume mode)");
        let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
        let num_genomes = genome_files.len();
        // Count unitigs by counting '>' lines (fast, no RAM needed for huge files)
        let num_unitigs = {
            use std::io::BufRead;
            let f = std::io::BufReader::new(std::fs::File::open(&unitig_file)?);
            f.lines().filter(|l| l.as_ref().map(|s| s.starts_with('>')).unwrap_or(false)).count()
        };
        log::info!("  Resuming with {} unitigs, {} genomes", num_unitigs, num_genomes);
        dbg::DbgResult {
            unitig_file,
            color_file,
            num_genomes,
            num_unitigs,
            kmer_size,
        }
    } else {
        // Step 1: Build colored compacted de Bruijn graph
        log::info!("Step 1/5: Building colored compacted de Bruijn graph...");
        dbg::build_cdbg(genome_dir, output_dir, kmer_size, threads)?
    };

    // Step 2: Parse unitigs and encode in 2-bit format
    log::info!("Step 2/5: Encoding unitigs...");
    let unitigs = unitig::parse_and_encode_unitigs(&dbg_result.unitig_file)?;

    // Step 3: Build color index (Roaring Bitmaps)
    log::info!("Step 3/5: Building color index...");
    color::build_color_index(&dbg_result.color_file, output_dir, dbg_result.num_genomes)?;

    // Step 4: Build FM-index over concatenated unitigs
    log::info!("Step 4/5: Building FM-index...");
    fm::build_fm_index_with_ram(&unitigs, output_dir, max_ram_bytes)?;

    // Step 5: Build genome path index
    log::info!("Step 5/5: Building genome path index...");
    paths::build_path_index(genome_dir, &unitigs, output_dir, kmer_size)?;

    // Step 6: Build specificity index (private unitig sets per genome)
    log::info!("Step 6: Building specificity index...");
    let color_index = color::load_color_index(output_dir)?;
    let spec_index = specificity::SpecificityIndex::build(
        &color_index,
        specificity::DEFAULT_MAX_SHARING,
    )?;
    spec_index.save(output_dir)?;

    // Write metadata
    metadata::write_metadata(output_dir, &dbg_result, &unitigs)?;

    log::info!("Index construction complete. Output: {:?}", output_dir);
    Ok(())
}
