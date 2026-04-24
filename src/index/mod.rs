pub mod auto_batch;
pub mod color;
pub mod dbg;
pub mod fm;
pub mod ggcat_colors;
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

    // Check if Step 1 can be skipped (resume from existing outputs)
    let unitig_file = output_dir.join("unitigs.fa");
    let color_file = output_dir.join("colors.tsv");
    let colors_drgn_path = output_dir.join("colors.drgn");
    let fm_path = output_dir.join("fm_index.bin");
    let has_colors = color_file.exists() || colors_drgn_path.exists();

    // If fm_index.bin + colors.drgn already exist, we can skip GGCAT entirely.
    // unitigs.fa isn't needed (Step 2 reconstructs UnitigSet from fm_index.bin).
    let dbg_result = if fm_path.exists() && colors_drgn_path.exists() {
        log::info!("Step 1/5: SKIPPED — fm_index.bin + colors.drgn already exist (full resume mode)");
        let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
        let num_genomes = genome_files.len();
        log::info!("  Resuming with existing FM-index, {} genomes", num_genomes);
        dbg::DbgResult {
            unitig_file: unitig_file.clone(),  // may not exist, used only for Step 2 which also skips
            color_file: color_file.clone(),
            num_genomes,
            num_unitigs: 0,  // recomputed from fm_index in Step 2
            kmer_size,
        }
    } else if unitig_file.exists() && has_colors {
        log::info!("Step 1/5: SKIPPED — existing unitigs.fa + colors found (resume mode)");
        let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
        let num_genomes = genome_files.len();
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

    let fm_exists = output_dir.join("fm_index.bin").exists();
    let colors_drgn_exists = output_dir.join("colors.drgn").exists();

    // Step 2: Parse unitigs and encode in 2-bit format
    let unitigs = if fm_exists {
        log::info!("Step 2/5: SKIPPED — fm_index.bin exists, loading unitigs from it");
        let fm_index = fm::load_fm_index(output_dir)?;
        let unitig_lengths: Vec<u64> = fm_index.cumulative_lengths.lengths().to_vec();
        unitig::UnitigSet::from_fm_text(&fm_index.text, &unitig_lengths)
    } else {
        log::info!("Step 2/5: Encoding unitigs...");
        unitig::parse_and_encode_unitigs(&dbg_result.unitig_file)?
    };

    // Step 3: Build color index (Roaring Bitmaps)
    if colors_drgn_exists {
        log::info!("Step 3/5: SKIPPED — colors.drgn already exists");
    } else {
        log::info!("Step 3/5: Building color index...");
        color::build_color_index(&dbg_result.color_file, output_dir, dbg_result.num_genomes)?;
    }

    // Step 4: Build FM-index over concatenated unitigs
    if fm_exists {
        log::info!("Step 4/5: SKIPPED — fm_index.bin already exists");
    } else {
        log::info!("Step 4/5: Building FM-index...");
        fm::build_fm_index_with_ram(&unitigs, output_dir, max_ram_bytes)?;
    }

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
