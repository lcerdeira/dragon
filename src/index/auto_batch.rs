/// Auto-batching index builder for large genome collections.
///
/// When a user runs `dragon index` on a huge genome directory (e.g., 100K+
/// genomes), this module transparently splits the build into smaller batches
/// that fit available memory, then merges them via the overlay architecture.
///
/// The user interface is unchanged — they just run:
///
/// ```text
/// dragon index -i genomes/ -o my_index/ --auto
/// ```
///
/// Dragon:
/// 1. Counts input genomes
/// 2. Estimates per-batch RAM requirement (~30 MB per genome at k=31)
/// 3. If it fits in available RAM: single build (fast path)
/// 4. Otherwise: splits into N batches, builds each, merges as overlays
///
/// Queries automatically traverse the overlay structure (via
/// `search_with_overlays` from Phase 3B) so the end user sees one unified index.

use anyhow::{Context, Result};
use std::path::Path;

/// Estimate the number of batches needed based on genome count and available RAM.
///
/// Rule of thumb: each genome needs ~30 MB of peak RAM during index build at k=31.
/// We also reserve 20% of RAM for OS/overhead.
pub fn estimate_batches(num_genomes: usize, available_ram_bytes: usize) -> usize {
    const BYTES_PER_GENOME: usize = 30 * 1024 * 1024; // 30 MB
    const SAFETY_FACTOR: f64 = 0.8;

    let usable = (available_ram_bytes as f64 * SAFETY_FACTOR) as usize;
    let per_batch = usable / BYTES_PER_GENOME.max(1);
    if per_batch == 0 {
        return num_genomes; // one genome per batch as fallback
    }
    (num_genomes + per_batch - 1) / per_batch.max(1)
}

/// Build an index with automatic batching if the genome count exceeds what
/// fits in the given RAM budget.
pub fn build_index_auto(
    genome_dir: &Path,
    output_dir: &Path,
    kmer_size: usize,
    threads: usize,
    max_ram_bytes: Option<usize>,
) -> Result<()> {
    let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
    let num_genomes = genome_files.len();

    // Determine RAM budget
    let ram_budget = max_ram_bytes.unwrap_or_else(|| {
        // Default: 90% of total system memory
        let sys_ram = detect_system_ram().unwrap_or(64 * 1024 * 1024 * 1024);
        (sys_ram as f64 * 0.9) as usize
    });

    let num_batches = estimate_batches(num_genomes, ram_budget);

    if num_batches <= 1 {
        // Single build — fast path
        log::info!("Building single index ({} genomes, fits in {} GB RAM)",
            num_genomes, ram_budget / 1_073_741_824);
        return super::build_index_with_options(genome_dir, output_dir, kmer_size, threads, max_ram_bytes);
    }

    // Auto-batching path
    log::info!(
        "Large collection ({} genomes > ~{} per batch) — splitting into {} batches",
        num_genomes,
        num_genomes / num_batches,
        num_batches
    );

    std::fs::create_dir_all(output_dir)?;

    // Build each batch to a temp location, then overlay into the main index
    let batch_size = (num_genomes + num_batches - 1) / num_batches;
    let temp_root = output_dir.join(".auto_batch_staging");
    std::fs::create_dir_all(&temp_root)?;

    for batch_idx in 0..num_batches {
        let start = batch_idx * batch_size;
        let end = ((batch_idx + 1) * batch_size).min(num_genomes);
        log::info!("=== Batch {}/{}: genomes {}..{} ===",
            batch_idx + 1, num_batches, start, end);

        // Create a directory with symlinks to this batch's genomes
        let batch_input = temp_root.join(format!("batch_{:03}_input", batch_idx));
        std::fs::create_dir_all(&batch_input)?;
        for (i, src) in genome_files[start..end].iter().enumerate() {
            let dst = batch_input.join(format!("genome_{:06}.fa", i));
            // Use hard link when possible (same FS), fallback to symlink
            if std::fs::hard_link(src, &dst).is_err() {
                let _ = std::os::unix::fs::symlink(src, &dst);
            }
        }

        if batch_idx == 0 {
            // Batch 0 is the base index
            super::build_index_with_options(&batch_input, output_dir, kmer_size, threads, max_ram_bytes)
                .with_context(|| format!("batch {} (base index)", batch_idx))?;
        } else {
            // Subsequent batches become overlays
            log::info!("Adding as overlay...");
            super::update::add_genomes(output_dir, &batch_input, kmer_size, threads, max_ram_bytes)
                .with_context(|| format!("batch {} (overlay)", batch_idx))?;
        }

        // Clean up batch input to free disk
        std::fs::remove_dir_all(&batch_input)?;
    }

    // Clean up staging directory
    let _ = std::fs::remove_dir_all(&temp_root);

    log::info!(
        "Auto-batch index complete: {} genomes across {} overlays",
        num_genomes, num_batches - 1
    );
    Ok(())
}

/// Detect total system RAM in bytes.
fn detect_system_ram() -> Option<usize> {
    // Read /proc/meminfo on Linux
    if let Ok(contents) = std::fs::read_to_string("/proc/meminfo") {
        for line in contents.lines() {
            if let Some(rest) = line.strip_prefix("MemTotal:") {
                let trimmed = rest.trim();
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if let Some(first) = parts.first() {
                    if let Ok(kb) = first.parse::<usize>() {
                        return Some(kb * 1024);
                    }
                }
            }
        }
    }
    // Fallback for macOS/others: use sysctl if available
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_batches_small() {
        // 100 genomes in 64 GB RAM: easily fits in one batch
        assert_eq!(estimate_batches(100, 64 * 1024 * 1024 * 1024), 1);
    }

    #[test]
    fn test_estimate_batches_large() {
        // 100,000 genomes in 64 GB RAM: needs multiple batches
        let n = estimate_batches(100_000, 64 * 1024 * 1024 * 1024);
        assert!(n > 1, "Expected >1 batch, got {}", n);
        // At 30 MB/genome, 64 GB * 0.8 / 30 MB ≈ 1747 genomes/batch
        // 100000 / 1747 ≈ 58 batches
        assert!(n > 10 && n < 200, "Batch count {} outside expected range", n);
    }

    #[test]
    fn test_estimate_batches_tiny_ram() {
        // 100 genomes in 1 GB RAM: should still batch
        let n = estimate_batches(100, 1024 * 1024 * 1024);
        assert!(n >= 1);
    }

    #[test]
    fn test_estimate_batches_zero_ram() {
        // Edge case: 0 RAM budget
        let n = estimate_batches(50, 0);
        assert_eq!(n, 50);
    }
}
