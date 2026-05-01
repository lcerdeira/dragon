/// Colored compacted de Bruijn graph construction via GGCAT.

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

/// Result of de Bruijn graph construction.
pub struct DbgResult {
    pub unitig_file: PathBuf,
    pub color_file: PathBuf,
    pub num_genomes: usize,
    pub num_unitigs: usize,
    pub kmer_size: usize,
}

/// Build a colored compacted de Bruijn graph using GGCAT.
///
/// GGCAT is invoked as an external binary. If not available, falls back
/// to a simplified internal implementation for small datasets.
pub fn build_cdbg(
    genome_dir: &Path,
    output_dir: &Path,
    kmer_size: usize,
    threads: usize,
) -> Result<DbgResult> {
    let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
    let num_genomes = genome_files.len();

    log::info!(
        "Building ccdBG from {} genomes with k={}",
        num_genomes,
        kmer_size
    );

    // Check if GGCAT is available
    if is_ggcat_available() {
        build_with_ggcat(genome_dir, output_dir, kmer_size, threads, num_genomes)
    } else {
        log::warn!("GGCAT not found in PATH. Using internal de Bruijn graph builder.");
        log::warn!("For best performance, install GGCAT: https://github.com/algbio/ggcat");
        build_internal(genome_dir, output_dir, kmer_size, num_genomes)
    }
}

fn is_ggcat_available() -> bool {
    std::process::Command::new("ggcat")
        .arg("--version")
        .output()
        .is_ok()
}

fn build_with_ggcat(
    genome_dir: &Path,
    output_dir: &Path,
    kmer_size: usize,
    threads: usize,
    num_genomes: usize,
) -> Result<DbgResult> {
    let unitig_file = output_dir.join("unitigs.fa");

    // Create input file list with absolute paths (GGCAT resolves relative to its CWD)
    let file_list = output_dir.join("input_files.txt");
    let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
    let file_list_content: String = genome_files
        .iter()
        .map(|p| std::fs::canonicalize(p).unwrap_or_else(|_| p.clone()).to_string_lossy().to_string())
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(&file_list, file_list_content)?;

    // Use output_dir for GGCAT temp files to avoid filling the root disk
    let temp_dir = output_dir.join(".ggcat_temp");
    std::fs::create_dir_all(&temp_dir)?;

    let status = std::process::Command::new("ggcat")
        .args([
            "build",
            "-k",
            &kmer_size.to_string(),
            "-j",
            &threads.to_string(),
            "--colors",
            "--temp-dir",
            &temp_dir.to_string_lossy(),
            "-l",
            &file_list.to_string_lossy(),
            "-o",
            &unitig_file.to_string_lossy(),
            "-s",
            "1",
        ])
        .status()
        .context("Failed to run GGCAT")?;

    // Clean up GGCAT temp files
    let _ = std::fs::remove_dir_all(&temp_dir);

    if !status.success() {
        anyhow::bail!("GGCAT exited with non-zero status: {}", status);
    }

    // Build colors.drgn DIRECTLY from GGCAT binary colormap — skipping the
    // colors.tsv intermediate which can be 300+ GB for large indices.
    log::info!("Building color index directly from GGCAT colormap...");
    let colormap_path = output_dir.join("unitigs.colors.dat");

    if !colormap_path.exists() {
        anyhow::bail!("GGCAT colormap not found at {:?} (GGCAT may have failed)", colormap_path);
    }

    let header = crate::index::ggcat_colors::build_color_drgn_direct(
        &colormap_path,
        &unitig_file,
        output_dir,
    ).context("Failed to build color index from GGCAT colormap")?;
    log::info!(
        "  Built colors.drgn: {} color subsets across {} genomes",
        header.subsets_count, header.colors_count
    );

    // Empty colors.tsv placeholder (some downstream code expects its path)
    let color_file = output_dir.join("colors.tsv");
    if !color_file.exists() {
        std::fs::write(&color_file, b"")?;
    }

    // Count unitigs
    let seqs = crate::io::fasta::read_sequences(&unitig_file)?;
    let num_unitigs = seqs.len();

    Ok(DbgResult {
        unitig_file,
        color_file,
        num_genomes,
        num_unitigs,
        kmer_size,
    })
}

/// Internal simplified de Bruijn graph builder for small datasets.
/// Builds a proper compacted de Bruijn graph by:
/// 1. Collecting all k-mers with genome colors
/// 2. Building an adjacency graph (suffix -> prefix overlap)
/// 3. Compacting non-branching paths into unitigs
/// For large datasets, GGCAT should be used instead.
fn build_internal(
    genome_dir: &Path,
    output_dir: &Path,
    kmer_size: usize,
    num_genomes: usize,
) -> Result<DbgResult> {
    use std::collections::{BTreeSet, HashMap, HashSet};
    use std::io::Write;

    let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;

    // Step 1: Collect all canonical k-mers with their genome colors
    let mut kmer_colors: HashMap<Vec<u8>, BTreeSet<usize>> = HashMap::new();

    for (genome_id, path) in genome_files.iter().enumerate() {
        let sequences = crate::io::fasta::read_sequences(path)?;
        for seq in &sequences {
            if seq.seq.len() < kmer_size {
                continue;
            }
            let upper: Vec<u8> = seq.seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            for i in 0..=upper.len() - kmer_size {
                let kmer = upper[i..i + kmer_size].to_vec();
                // Only keep k-mers with valid bases
                if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    kmer_colors.entry(kmer).or_default().insert(genome_id);
                }
            }
        }
    }

    let total_kmers = kmer_colors.len();
    log::info!("Internal builder: {} unique k-mers, compacting...", total_kmers);

    // Step 2: Build adjacency graph using (k-1)-mer overlaps
    // For each k-mer, its suffix (last k-1 bases) can overlap with another k-mer's prefix (first k-1 bases)
    let k1 = kmer_size - 1;

    // Map from (k-1)-prefix to list of k-mers starting with that prefix
    // prefix_map: (k-1)-prefix -> k-mers starting with that prefix
    let mut prefix_map: HashMap<&[u8], Vec<&Vec<u8>>> = HashMap::new();
    // suffix_map: (k-1)-suffix -> k-mers ending with that suffix
    let mut suffix_map: HashMap<&[u8], Vec<&Vec<u8>>> = HashMap::new();

    for kmer in kmer_colors.keys() {
        prefix_map.entry(&kmer[..k1]).or_default().push(kmer);
        suffix_map.entry(&kmer[kmer_size - k1..]).or_default().push(kmer);
    }

    // Successors of kmer X: k-mers whose prefix = X's suffix (last k-1 bases)
    let get_successors = |kmer: &[u8]| -> Vec<&Vec<u8>> {
        let suf = &kmer[kmer_size - k1..]; // last k-1 bases
        prefix_map.get(suf).cloned().unwrap_or_default()
    };

    // Predecessors of kmer X: k-mers whose suffix = X's prefix (first k-1 bases)
    let get_predecessors = |kmer: &[u8]| -> Vec<&Vec<u8>> {
        let pre = &kmer[..k1]; // first k-1 bases
        suffix_map.get(pre).cloned().unwrap_or_default()
    };

    // Step 3: Compact non-branching paths
    // A non-branching node has exactly 1 successor AND 1 predecessor in the dBG.
    // We merge chains of such nodes into unitigs.
    // Unitig colors = intersection of all component k-mer colors.
    let mut used: HashSet<Vec<u8>> = HashSet::new();
    let mut unitigs: Vec<(Vec<u8>, BTreeSet<usize>)> = Vec::new();

    let kmer_list: Vec<Vec<u8>> = kmer_colors.keys().cloned().collect();

    for start_kmer in &kmer_list {
        if used.contains(start_kmer) {
            continue;
        }

        // Walk backward to find the true start of this chain
        let mut chain_start = start_kmer.clone();
        loop {
            let preds = get_predecessors(&chain_start);
            if preds.len() != 1 || used.contains(preds[0]) {
                break;
            }
            let pred = preds[0];
            // Check pred has exactly 1 successor (so the path is non-branching)
            let pred_succs = get_successors(pred);
            if pred_succs.len() != 1 {
                break;
            }
            if *pred == chain_start {
                break; // self-loop
            }
            chain_start = pred.clone();
        }

        if used.contains(&chain_start) {
            continue;
        }

        // Walk forward from chain_start, extending the unitig
        let mut unitig_seq = chain_start.clone();
        let unitig_colors = kmer_colors.get(&chain_start).cloned().unwrap_or_default();
        used.insert(chain_start.clone());

        let mut current = chain_start;
        loop {
            let succs = get_successors(&current);
            if succs.len() != 1 {
                break;
            }
            let next = succs[0];
            if used.contains(next) {
                break;
            }
            // Check that next has exactly 1 predecessor (non-branching)
            let next_preds = get_predecessors(next);
            if next_preds.len() != 1 {
                break;
            }
            // Only merge if color sets are identical (preserves color accuracy)
            if kmer_colors.get(next) != Some(&unitig_colors) {
                break;
            }
            // Extend unitig by appending the last base of the successor
            unitig_seq.push(next[kmer_size - 1]);
            used.insert(next.clone());
            current = next.clone();
        }

        unitigs.push((unitig_seq, unitig_colors));
    }

    // Add any remaining unused k-mers as single-kmer unitigs
    for kmer in &kmer_list {
        if !used.contains(kmer) {
            let colors = kmer_colors.get(kmer).cloned().unwrap_or_default();
            unitigs.push((kmer.clone(), colors));
        }
    }

    // Sort unitigs for deterministic output
    unitigs.sort_by(|a, b| a.0.cmp(&b.0));

    let avg_len = if unitigs.is_empty() { 0 } else {
        unitigs.iter().map(|(s, _)| s.len()).sum::<usize>() / unitigs.len()
    };
    let max_len = unitigs.iter().map(|(s, _)| s.len()).max().unwrap_or(0);
    log::info!(
        "Compacted {} k-mers into {} unitigs (avg len: {} bp, max: {} bp)",
        total_kmers, unitigs.len(), avg_len, max_len
    );

    // Step 4: Write output
    let unitig_file = output_dir.join("unitigs.fa");
    let color_file = output_dir.join("colors.tsv");

    let mut fa_writer = std::io::BufWriter::new(std::fs::File::create(&unitig_file)?);
    let mut color_writer = std::io::BufWriter::new(std::fs::File::create(&color_file)?);

    for (uid, (seq, colors)) in unitigs.iter().enumerate() {
        writeln!(fa_writer, ">unitig_{}", uid)?;
        fa_writer.write_all(seq)?;
        writeln!(fa_writer)?;

        let color_str: String = colors
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(color_writer, "{}\t{}", uid, color_str)?;
    }

    let num_unitigs = unitigs.len();

    Ok(DbgResult {
        unitig_file,
        color_file,
        num_genomes,
        num_unitigs,
        kmer_size,
    })
}
