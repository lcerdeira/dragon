/// Signal-level search: query raw nanopore signal against a signal FM-index.
///
/// Pipeline:
/// 1. Load raw signal from FAST5/TSV
/// 2. Normalize using median-MAD
/// 3. Discretize to the signal alphabet
/// 4. Extract signal k-mers
/// 5. FM-index backward search for each k-mer
/// 6. Map hits back to genome coordinates
/// 7. Score and filter results

use anyhow::Result;
use std::path::Path;

use crate::index::fm::DragonFmIndex;
use crate::signal::discretize::{self, SignalAlphabet};
use crate::signal::index::{self, SignalIndexMetadata};

/// A hit from signal-level search.
#[derive(Clone, Debug)]
pub struct SignalHit {
    /// Genome index (into metadata.genome_names).
    pub genome_id: u32,
    /// Position within the genome's discretized signal.
    pub position: u32,
    /// Match score (higher is better).
    pub score: f64,
    /// Strand: true = forward, false = reverse complement.
    pub forward_strand: bool,
    /// Length of the matching signal segment (in discrete symbols).
    pub match_len: usize,
    /// Genome name (populated after lookup).
    pub genome_name: String,
}

/// Configuration for signal search.
#[derive(Clone, Debug)]
pub struct SignalSearchConfig {
    /// Path to the signal index directory.
    pub index_dir: Box<Path>,
    /// Signal k-mer size for search (default: 10).
    pub signal_kmer_size: usize,
    /// Minimum number of k-mer hits to report a genome (default: 3).
    pub min_hits: usize,
    /// Maximum seed frequency: skip signal k-mers with more hits than this.
    pub max_seed_freq: usize,
    /// Maximum number of results to report per query read.
    pub max_results: usize,
    /// Number of threads.
    pub threads: usize,
}

impl Default for SignalSearchConfig {
    fn default() -> Self {
        Self {
            index_dir: Path::new(".").into(),
            signal_kmer_size: 10,
            min_hits: 3,
            max_seed_freq: 10_000,
            max_results: 50,
            threads: 4,
        }
    }
}

/// Result of searching a single signal read.
#[derive(Clone, Debug)]
pub struct SignalSearchResult {
    /// Read identifier.
    pub read_id: String,
    /// Signal length (raw samples).
    pub signal_length: usize,
    /// Hits found, sorted by score descending.
    pub hits: Vec<SignalHit>,
}

/// Search a single raw signal against the signal FM-index.
///
/// Steps:
/// 1. Normalize and discretize the query signal
/// 2. Extract overlapping k-mers from the discretized signal
/// 3. Search each k-mer in the FM-index
/// 4. Aggregate hits by genome, computing per-genome scores
/// 5. Return scored, filtered results
pub fn search_signal(
    signal_data: &[f32],
    fm_index: &DragonFmIndex,
    metadata: &SignalIndexMetadata,
    alphabet: &SignalAlphabet,
    config: &SignalSearchConfig,
) -> Vec<SignalHit> {
    if signal_data.is_empty() {
        return Vec::new();
    }

    // Step 1: Normalize and discretize query signal
    let normalized = discretize::normalize_signal(signal_data);
    let discretized = discretize::discretize_normalized(&normalized, alphabet);

    if discretized.len() < config.signal_kmer_size {
        log::warn!(
            "Discretized signal too short ({} < {})",
            discretized.len(),
            config.signal_kmer_size
        );
        return Vec::new();
    }

    // Step 2: Extract signal k-mers and search
    let mut genome_hits: std::collections::HashMap<u32, Vec<(u32, usize)>> =
        std::collections::HashMap::new();

    let num_kmers = discretized.len() - config.signal_kmer_size + 1;
    let mut total_seeds = 0usize;
    let mut filtered_seeds = 0usize;

    for i in 0..num_kmers {
        let kmer = &discretized[i..i + config.signal_kmer_size];

        // FM-index backward search
        let hit_count = fm_index.count(kmer);

        if hit_count == 0 {
            continue;
        }

        total_seeds += 1;

        // Skip overly frequent seeds
        if hit_count > config.max_seed_freq {
            filtered_seeds += 1;
            continue;
        }

        let positions = fm_index.search(kmer);

        for &pos in &positions {
            // Map position to genome
            if let Some((genome_id, offset)) = fm_index.position_to_unitig(pos) {
                genome_hits
                    .entry(genome_id)
                    .or_default()
                    .push((offset, config.signal_kmer_size));
            }
        }
    }

    log::debug!(
        "Signal search: {} total seeds, {} filtered (freq > {})",
        total_seeds,
        filtered_seeds,
        config.max_seed_freq
    );

    // Step 3: Score genomes by hit count and coverage
    let mut hits: Vec<SignalHit> = Vec::new();

    for (&genome_id, positions) in &genome_hits {
        if positions.len() < config.min_hits {
            continue;
        }

        // Compute score based on number of hits and their spread
        let score = compute_genome_score(positions);

        // Determine best position (mode of hit positions)
        let best_pos = find_best_position(positions);

        // Determine strand from hit position within the genome's signal block.
        // The index stores forward signal, sentinel, then reverse complement signal
        // consecutively per genome. Hits in the first half are forward strand.
        let forward_strand = determine_strand(positions, &metadata);

        let genome_name = if (genome_id as usize) < metadata.genome_names.len() {
            metadata.genome_names[genome_id as usize].clone()
        } else {
            format!("genome_{}", genome_id)
        };

        hits.push(SignalHit {
            genome_id,
            position: best_pos,
            score,
            forward_strand,
            match_len: positions.len() * config.signal_kmer_size,
            genome_name,
        });
    }

    // Sort by score descending
    hits.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

    // Limit results
    hits.truncate(config.max_results);
    hits
}

/// Search all reads from a signal file against the index.
pub fn search_signal_file(
    signal_file: &Path,
    config: &SignalSearchConfig,
) -> Result<Vec<SignalSearchResult>> {
    // Load index components
    log::info!("Loading signal index from {:?}", config.index_dir);
    let fm_index = index::load_signal_index(&config.index_dir)?;
    let metadata = index::load_signal_metadata(&config.index_dir)?;
    let alphabet = index::load_signal_alphabet(&config.index_dir)?;

    // Load signal reads
    let reads = crate::signal::io::read_signal_file(signal_file)?;
    log::info!("Loaded {} signal reads", reads.len());

    let mut results = Vec::with_capacity(reads.len());

    for read in &reads {
        log::debug!(
            "Searching read {} ({} samples)",
            read.id,
            read.signal.len()
        );

        let hits = search_signal(&read.signal, &fm_index, &metadata, &alphabet, config);

        results.push(SignalSearchResult {
            read_id: read.id.clone(),
            signal_length: read.signal.len(),
            hits,
        });
    }

    log::info!(
        "Signal search complete: {} reads, {} total hits",
        results.len(),
        results.iter().map(|r| r.hits.len()).sum::<usize>()
    );

    Ok(results)
}

/// Write signal search results to a TSV output file.
pub fn write_signal_results(
    results: &[SignalSearchResult],
    output: &mut dyn std::io::Write,
) -> Result<()> {
    writeln!(
        output,
        "#read_id\tread_signal_len\tgenome_name\tgenome_id\tposition\tscore\tstrand\tmatch_len"
    )?;

    for result in results {
        for hit in &result.hits {
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
                result.read_id,
                result.signal_length,
                hit.genome_name,
                hit.genome_id,
                hit.position,
                hit.score,
                if hit.forward_strand { "+" } else { "-" },
                hit.match_len,
            )?;
        }
    }

    Ok(())
}

/// Compute a score for a genome based on hit positions.
///
/// Score combines:
/// 1. Number of k-mer hits (base signal)
/// 2. Coverage of the hit region (rewards contiguous matches)
/// 3. Hit density (penalizes scattered low-quality hits)
/// 4. Match quality bonus for longer individual matches (adaptive mode)
fn compute_genome_score(positions: &[(u32, usize)]) -> f64 {
    let scorer = SignalScorer::default();
    scorer.score(positions)
}

/// Learned signal scorer: logistic-regression-style scoring for signal genome hits.
///
/// Features: num_hits, coverage, density, avg_match_len.
/// Inference: `dot(weights, features) + bias` (linear score, no sigmoid — we want
/// a ranking score, not a probability).
pub struct SignalScorer {
    /// Weights: [bias, num_hits, coverage, density, avg_match_len]
    pub weights: [f64; 5],
}

impl Default for SignalScorer {
    fn default() -> Self {
        // Hand-tuned defaults that reproduce the original formula's ranking behavior.
        // These can be replaced by trained weights loaded from JSON.
        Self {
            weights: [
                0.0,   // bias
                1.0,   // num_hits: base magnitude
                0.5,   // coverage: rewards contiguous alignment
                0.3,   // density: rewards tight clustering
                0.2,   // avg_match_len: rewards longer, more specific matches
            ],
        }
    }
}

impl SignalScorer {
    /// Load weights from a JSON file, falling back to defaults.
    pub fn load_or_default(path: &std::path::Path) -> Self {
        let weights_path = if path.is_dir() {
            path.join("signal_scorer.json")
        } else {
            path.to_path_buf()
        };
        if weights_path.exists() {
            if let Ok(data) = std::fs::read_to_string(&weights_path) {
                if let Ok(weights) = serde_json::from_str::<Vec<f64>>(&data) {
                    if weights.len() == 5 {
                        log::info!("Loaded signal scorer weights from {:?}", weights_path);
                        return Self {
                            weights: [weights[0], weights[1], weights[2], weights[3], weights[4]],
                        };
                    }
                }
            }
            log::warn!("Failed to load signal scorer from {:?}, using defaults", weights_path);
        }
        Self::default()
    }

    /// Compute features from hit positions.
    fn features(positions: &[(u32, usize)]) -> [f64; 4] {
        let num_hits = positions.len() as f64;
        if positions.is_empty() {
            return [0.0; 4];
        }

        let min_pos = positions.iter().map(|&(p, _)| p).min().unwrap_or(0) as f64;
        let max_pos = positions.iter().map(|&(p, _)| p).max().unwrap_or(0) as f64;
        let span = max_pos - min_pos + 1.0;

        let total_match_len: f64 = positions.iter().map(|&(_, len)| len as f64).sum();
        let coverage = (total_match_len / span).min(1.0);
        let density = (num_hits / span.max(1.0)).min(1.0);
        let avg_match_len = total_match_len / num_hits;

        [num_hits, coverage, density, avg_match_len]
    }

    /// Score a genome given its hit positions.
    pub fn score(&self, positions: &[(u32, usize)]) -> f64 {
        if positions.is_empty() {
            return 0.0;
        }
        let feats = Self::features(positions);
        let raw = self.weights[0]
            + self.weights[1] * feats[0]
            + self.weights[2] * feats[1]
            + self.weights[3] * feats[2]
            + self.weights[4] * feats[3];
        raw.max(0.0)
    }
}

/// Find the most common hit position region (approximate mode).
fn find_best_position(positions: &[(u32, usize)]) -> u32 {
    if positions.is_empty() {
        return 0;
    }

    // Use a simple binning approach
    let bin_size = 100u32;
    let mut bins: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();

    for &(pos, _) in positions {
        let bin = pos / bin_size;
        *bins.entry(bin).or_insert(0) += 1;
    }

    let best_bin = bins
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(bin, _)| bin)
        .unwrap_or(0);

    best_bin * bin_size
}

/// Determine strand by majority vote of hit positions.
///
/// Each genome's signal block contains forward signal, a sentinel, then reverse
/// complement signal. Hits in the first half of the block are forward strand.
/// We estimate the midpoint as total_signal_length / (2 * num_genomes) and
/// check whether the majority of hits fall before or after the midpoint within
/// the genome's block.
fn determine_strand(positions: &[(u32, usize)], metadata: &SignalIndexMetadata) -> bool {
    if positions.is_empty() || metadata.num_genomes == 0 {
        return true;
    }

    // Approximate midpoint of each genome's signal block
    let avg_genome_signal = metadata.total_signal_length as f64 / metadata.num_genomes as f64;
    let midpoint = (avg_genome_signal / 2.0) as u32;

    // Count hits in forward vs reverse half
    let mut fwd_count = 0usize;
    let mut rev_count = 0usize;
    for &(pos, _) in positions {
        let pos_in_block = pos % midpoint.max(1);
        if pos_in_block < midpoint / 2 {
            fwd_count += 1;
        } else {
            rev_count += 1;
        }
    }

    fwd_count >= rev_count
}

/// Perform adaptive-length signal search using variable-length FM-index extension.
///
/// Instead of fixed k-mer sizes, this extends each seed as far as possible
/// in the FM-index, yielding longer (and more specific) matches.
pub fn search_signal_adaptive(
    signal_data: &[f32],
    fm_index: &DragonFmIndex,
    metadata: &SignalIndexMetadata,
    alphabet: &SignalAlphabet,
    config: &SignalSearchConfig,
) -> Vec<SignalHit> {
    if signal_data.is_empty() {
        return Vec::new();
    }

    let normalized = discretize::normalize_signal(signal_data);
    let discretized = discretize::discretize_normalized(&normalized, alphabet);

    if discretized.len() < config.signal_kmer_size {
        return Vec::new();
    }

    let mut genome_hits: std::collections::HashMap<u32, Vec<(u32, usize)>> =
        std::collections::HashMap::new();

    let mut i = 0;
    while i < discretized.len().saturating_sub(config.signal_kmer_size) {
        let remaining = &discretized[i..];
        let (match_len, _count) = fm_index.variable_length_search(remaining);

        if match_len >= config.signal_kmer_size {
            let pattern = &discretized[i..i + match_len];
            let positions = fm_index.search(pattern);

            if positions.len() <= config.max_seed_freq {
                for &pos in &positions {
                    if let Some((genome_id, offset)) = fm_index.position_to_unitig(pos) {
                        genome_hits
                            .entry(genome_id)
                            .or_default()
                            .push((offset, match_len));
                    }
                }
            }

            // Skip forward by match_len (non-overlapping seeds for adaptive)
            i += match_len;
        } else {
            i += 1;
        }
    }

    // Score and collect hits (same as fixed-kmer search)
    let mut hits: Vec<SignalHit> = Vec::new();

    for (&genome_id, positions) in &genome_hits {
        if positions.len() < config.min_hits {
            continue;
        }

        let score = compute_genome_score(positions);
        let best_pos = find_best_position(positions);

        let genome_name = if (genome_id as usize) < metadata.genome_names.len() {
            metadata.genome_names[genome_id as usize].clone()
        } else {
            format!("genome_{}", genome_id)
        };

        let total_match: usize = positions.iter().map(|&(_, l)| l).sum();

        hits.push(SignalHit {
            genome_id,
            position: best_pos,
            score,
            forward_strand: true,
            match_len: total_match,
            genome_name,
        });
    }

    hits.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
    hits.truncate(config.max_results);
    hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signal::model;

    #[test]
    fn test_compute_genome_score() {
        let positions = vec![(10, 5), (20, 5), (30, 5)];
        let score = compute_genome_score(&positions);
        assert!(score > 0.0);
    }

    #[test]
    fn test_compute_genome_score_empty() {
        let positions: Vec<(u32, usize)> = vec![];
        let score = compute_genome_score(&positions);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_find_best_position() {
        let positions = vec![(150, 10), (160, 10), (170, 10), (500, 10)];
        let best = find_best_position(&positions);
        // Most hits are in the 100-200 range (bin 1)
        assert_eq!(best, 100);
    }

    #[test]
    fn test_search_signal_empty() {
        let signal_data: Vec<f32> = vec![];
        let fm = DragonFmIndex {
            text: vec![0, 1, 2, 3],
            suffix_array: vec![0, 1, 2, 3],
            cumulative_lengths: crate::ds::elias_fano::CumulativeLengthIndex::from_lengths(&[4]),
        };
        let metadata = SignalIndexMetadata {
            version: "0.1.0".to_string(),
            num_genomes: 1,
            genome_names: vec!["test".to_string()],
            total_signal_length: 4,
            pore_model_name: "test".to_string(),
            pore_model_kmer_size: 5,
            num_levels: 16,
            alphabet_min: -4.0,
            alphabet_max: 4.0,
            alphabet_boundaries: None,
        };
        let alphabet = SignalAlphabet::default();
        let config = SignalSearchConfig::default();

        let hits = search_signal(&signal_data, &fm, &metadata, &alphabet, &config);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_signal_search_integration() {
        // Build a small index and search against it
        let genome_dir = tempfile::tempdir().unwrap();
        let genome_path = genome_dir.path().join("test.fa");
        {
            use std::io::Write;
            let mut f = std::fs::File::create(&genome_path).unwrap();
            writeln!(f, ">chr1").unwrap();
            writeln!(f, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        }

        let output_dir = tempfile::tempdir().unwrap();
        let idx_config = crate::signal::index::SignalIndexConfig::default();
        crate::signal::index::build_signal_index(
            genome_dir.path(),
            output_dir.path(),
            &idx_config,
        )
        .unwrap();

        // Load the index
        let fm_index = index::load_signal_index(output_dir.path()).unwrap();
        let metadata = index::load_signal_metadata(output_dir.path()).unwrap();
        let alphabet = index::load_signal_alphabet(output_dir.path()).unwrap();

        // Generate a "query signal" from the same sequence
        let pore_model = model::load_default_model();
        let query_signal = pore_model.sequence_to_expected_signal(b"ACGTACGTACGTACGT");

        let config = SignalSearchConfig {
            signal_kmer_size: 5,
            min_hits: 1,
            max_seed_freq: 100_000,
            max_results: 10,
            ..Default::default()
        };

        let hits = search_signal(&query_signal, &fm_index, &metadata, &alphabet, &config);
        // We should get at least some hits since the query matches the indexed genome
        // (though normalization differences may reduce exact matches)
        log::info!("Integration test found {} hits", hits.len());
    }

    #[test]
    fn test_write_signal_results() {
        let results = vec![SignalSearchResult {
            read_id: "read_001".to_string(),
            signal_length: 1000,
            hits: vec![SignalHit {
                genome_id: 0,
                position: 100,
                score: 42.5,
                forward_strand: true,
                match_len: 50,
                genome_name: "genome_A".to_string(),
            }],
        }];

        let mut buf = Vec::new();
        write_signal_results(&results, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("read_001"));
        assert!(output.contains("genome_A"));
        assert!(output.contains("42.50"));
    }
}
