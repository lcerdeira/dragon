pub mod align;
pub mod candidate;
pub mod chain;
pub mod extract;
pub mod ml_score;
pub mod seed;

use anyhow::Result;
use std::collections::HashSet;
use std::path::Path;

use crate::io::paf::PafRecord;

/// Configuration for a search query.
pub struct SearchConfig {
    pub index_dir: Box<Path>,
    pub min_seed_len: usize,
    pub max_seed_freq: usize,
    pub min_chain_score: f64,
    pub max_target_seqs: usize,
    pub threads: usize,
    pub max_ram_gb: f64,
    /// Minimum alignment identity to report (0.0 - 1.0). Default: 0.7
    /// Uses real NW alignment identity when reference can be extracted,
    /// falls back to seed-density estimation for very large regions.
    pub min_identity: f64,
    /// Minimum query coverage to report (0.0 - 1.0). Default: 0.3
    pub min_query_coverage: f64,
    /// Minimum score ratio: only report hits with score >= ratio * best_score (0.0-1.0).
    pub min_score_ratio: f64,
    /// Disable ML seed scoring (use raw match_len instead).
    pub no_ml: bool,
    /// Path to custom ML weights file (overrides index-dir/seed_scorer.json).
    pub ml_weights_path: Option<String>,
    /// If set, dump all seeds with features to this TSV path (for ML training).
    pub dump_seeds_path: Option<String>,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            index_dir: Path::new(".").into(),
            min_seed_len: 15,
            max_seed_freq: 10_000,
            min_chain_score: 50.0,
            max_target_seqs: 10,
            threads: 4,
            max_ram_gb: 4.0,
            min_identity: 0.7,
            min_query_coverage: 0.3,
            min_score_ratio: 0.1,
            no_ml: false,
            ml_weights_path: None,
            dump_seeds_path: None,
        }
    }
}

/// Result of aligning a single query against the database.
pub struct QueryResult {
    pub query_name: String,
    pub query_len: usize,
    pub alignments: Vec<PafRecord>,
}

/// Search a set of query sequences against a Dragon index.
pub fn search(
    query_file: &Path,
    config: &SearchConfig,
) -> Result<Vec<QueryResult>> {
    log::info!("Loading index from {:?}", config.index_dir);

    // Load index components (memory-mapped)
    let fm_index = crate::index::fm::load_fm_index(&config.index_dir)?;
    let color_index = crate::index::color::load_color_index(&config.index_dir)?;
    let path_index = crate::index::paths::load_path_index(&config.index_dir)?;
    let _metadata = crate::index::metadata::load_metadata(&config.index_dir)?;

    // Reconstruct UnitigSet from FM-index text + cumulative lengths
    let unitig_lengths: Vec<u64> = fm_index.cumulative_lengths.lengths().to_vec();
    let unitigs = crate::index::unitig::UnitigSet::from_fm_text(&fm_index.text, &unitig_lengths);

    // Initialize ML seed scorer (unless disabled)
    let scorer = if config.no_ml {
        log::info!("ML seed scoring disabled (--no-ml)");
        None
    } else if let Some(ref weights_path) = config.ml_weights_path {
        let path = std::path::Path::new(weights_path);
        log::info!("Loading custom ML weights from {:?}", path);
        Some(ml_score::SeedScorer::load_or_default(path))
    } else {
        Some(ml_score::SeedScorer::load_or_default(&config.index_dir))
    };

    // Open seed dump file if requested
    let mut seed_dump: Option<Box<dyn std::io::Write>> = if let Some(ref path) = config.dump_seeds_path {
        log::info!("Dumping seed features to {:?}", path);
        let file = std::fs::File::create(path)?;
        let mut writer = Box::new(std::io::BufWriter::new(file)) as Box<dyn std::io::Write>;
        use std::io::Write;
        writeln!(writer, "query_name\tunitig_id\tquery_pos\tmatch_len\tsa_count\tgenome_id\tcolor_card\tgc_content")?;
        Some(writer)
    } else {
        None
    };

    // Parse queries
    let queries = crate::io::fasta::read_sequences(query_file)?;
    log::info!("Loaded {} queries", queries.len());

    let mut results = Vec::new();

    for query in &queries {
        log::debug!("Searching query: {} ({} bp)", query.name, query.seq.len());

        // Adaptive parameters for short reads
        let min_votes = if query.seq.len() < 300 { 2 } else { 10 };
        let effective_chain_score = if query.seq.len() < 300 {
            config.min_chain_score.min(20.0)
        } else {
            config.min_chain_score
        };

        // Stage 1: Find seeds via FM-index backward search
        let seeds = seed::find_seeds(
            &query.seq,
            &fm_index,
            config.min_seed_len,
            config.max_seed_freq,
        );

        // Stage 2: Identify candidate genomes via color voting
        let candidates = candidate::find_candidates(&seeds, &color_index, min_votes);

        // Dump seeds if requested (for ML training data generation)
        if let Some(ref mut writer) = seed_dump {
            use std::io::Write;
            for cand in &candidates {
                for seed in &cand.seeds {
                    let color_card = color_index
                        .get_colors(seed.unitig_id)
                        .map(|c| c.len())
                        .unwrap_or(0);
                    let gc = if seed.match_len > 0 && seed.query_pos + seed.match_len <= query.seq.len() {
                        let region = &query.seq[seed.query_pos..seed.query_pos + seed.match_len];
                        let gc = region.iter().filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c')).count();
                        gc as f64 / seed.match_len as f64
                    } else {
                        0.5
                    };
                    let _ = writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                        query.name, seed.unitig_id, seed.query_pos, seed.match_len,
                        seed.sa_count, cand.genome_id, color_card, gc
                    );
                }
            }
        }

        // Stage 3: Colinear chaining per candidate genome (with ML seed scoring)
        let chains = chain::chain_candidates_with_query_len(
            &seeds,
            &candidates,
            &path_index,
            effective_chain_score,
            query.seq.len(),
            scorer.as_ref(),
            &color_index,
            &query.seq,
        );

        // Stage 4: Extract reference and align
        let mut alignments = align::align_chains(
            &query.seq,
            &chains,
            &path_index,
            &unitigs,
        );

        // Stage 5: Post-alignment filtering — this is critical for precision
        alignments.retain(|record| {
            let identity = record.identity();
            identity >= config.min_identity
        });

        // Also filter by query coverage using the chain-computed coverage
        // (which is relative to full query length)
        // This must match against the chain; since we lose the chain reference
        // after alignment, use PAF query_start/query_end as the coverage measure.
        alignments.retain(|record| {
            let query_cov = if query.seq.len() > 0 {
                (record.query_end - record.query_start) as f64 / query.seq.len() as f64
            } else {
                0.0
            };
            query_cov >= config.min_query_coverage
        });

        // Score-ratio filter: discard hits with score < min_score_ratio * best_score
        if !alignments.is_empty() && config.min_score_ratio > 0.0 {
            let best_score = alignments[0]
                .tags
                .iter()
                .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            let threshold = best_score * config.min_score_ratio;
            alignments.retain(|record| {
                record
                    .tags
                    .iter()
                    .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                    .unwrap_or(0.0)
                    >= threshold
            });
        }

        // Per-genome deduplication: keep only the best chain per target genome
        // (alignments are sorted by score descending, so first occurrence is best)
        {
            let mut seen: HashSet<String> = HashSet::new();
            alignments.retain(|record| seen.insert(record.target_name.clone()));
        }

        // Enforce max_target_seqs: keep only the top N hits sorted by score
        // (alignments are already sorted by chain score descending from stage 3)
        alignments.truncate(config.max_target_seqs);

        results.push(QueryResult {
            query_name: query.name.clone(),
            query_len: query.seq.len(),
            alignments,
        });
    }

    Ok(results)
}
