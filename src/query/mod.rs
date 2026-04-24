pub mod align;
pub mod candidate;
pub mod chain;
pub mod containment;
pub mod direct_align;
pub mod extract;
pub mod ml_score;
pub mod seed;

use anyhow::Result;
use std::collections::HashSet;
use std::path::Path;

use crate::io::paf::PafRecord;

/// Configuration for a search query.
#[derive(Clone)]
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
    /// Ground truth genome name (for labeling training data).
    /// When set, seeds matching this genome get label=1, others get label=0.
    pub ground_truth_genome: Option<String>,
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
            ground_truth_genome: None,
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
    let _specificity_index = crate::index::specificity::SpecificityIndex::load_or_build(
        &config.index_dir,
        &color_index,
    ).ok();
    let metadata = crate::index::metadata::load_metadata(&config.index_dir)?;

    // Reconstruct UnitigSet from FM-index text + cumulative lengths
    let unitig_lengths: Vec<u64> = fm_index.cumulative_lengths.lengths().to_vec();
    let unitigs = crate::index::unitig::UnitigSet::from_fm_text(&fm_index.text, &unitig_lengths);

    // Initialize ML seed scorer (unless disabled)
    let _scorer = if config.no_ml {
        log::info!("ML seed scoring disabled (--no-ml)");
        None
    } else if let Some(ref weights_path) = config.ml_weights_path {
        let path = std::path::Path::new(weights_path);
        log::info!("Loading custom ML weights from {:?}", path);
        Some(ml_score::SeedScorer::load_or_default(path))
    } else {
        Some(ml_score::SeedScorer::load_or_default(&config.index_dir))
    };

    // Open seed dump file if requested (for ML training data)
    let mut seed_dump: Option<Box<dyn std::io::Write>> = if let Some(ref path) = config.dump_seeds_path {
        log::info!("Dumping seed features to {:?}", path);
        if let Some(ref gt) = config.ground_truth_genome {
            log::info!("Ground truth genome: {:?} (seeds will be labeled)", gt);
        }
        let file = std::fs::File::create(path)?;
        let mut writer = Box::new(std::io::BufWriter::new(file)) as Box<dyn std::io::Write>;
        use std::io::Write;
        writeln!(writer, "query_name\tgenome_id\tgenome_name\t{}\tlabel",
            ml_score::SeedFeatures::header())?;
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

        // Adaptive vote threshold: lower for short reads and small seed sets
        let base_min_votes = if query.seq.len() < 300 { 2u32 } else { 10u32 };
        let _effective_chain_score = if query.seq.len() < 300 {
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
        // Adaptive threshold: use at most 20% of seed count, but at least 2
        let min_votes = if seeds.is_empty() {
            base_min_votes
        } else {
            base_min_votes.min((seeds.len() as u32 / 5).max(2))
        };
        let candidates = candidate::find_candidates(&seeds, &color_index, min_votes);

        // Dump seeds with full ML features if requested (for training data generation)
        if let Some(ref mut writer) = seed_dump {
            use std::io::Write;
            // Collect all seed positions for local density computation
            let all_positions: Vec<usize> = candidates
                .iter()
                .flat_map(|c| c.seeds.iter().map(|s| s.query_pos))
                .collect();

            for cand in &candidates {
                let genome_name = metadata.genome_names
                    .get(cand.genome_id as usize)
                    .cloned()
                    .unwrap_or_else(|| format!("genome_{}", cand.genome_id));

                // Label: 1 if this genome matches ground truth, 0 otherwise, -1 if no ground truth
                let label = match &config.ground_truth_genome {
                    Some(gt) => if genome_name.contains(gt.as_str()) { 1i8 } else { 0i8 },
                    None => -1i8,
                };

                for seed in &cand.seeds {
                    let features = ml_score::SeedFeatures::from_seed_with_context(
                        seed,
                        query.seq.len(),
                        &query.seq,
                        &color_index,
                        &all_positions,
                    );
                    let _ = writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}",
                        query.name, cand.genome_id, genome_name,
                        features.as_tsv(), label
                    );
                }
            }
        }

        // Stage 3: Containment-based ranking (primary method)
        // Computes k-mer containment between query and each genome via color index.
        // This bypasses the unitig-boundary seeding problem by counting ALL k-mer
        // matches, not just those within single unitigs.
        let containment_hits = containment::containment_rank(
            &query.seq,
            &fm_index,
            &color_index,
            config.min_seed_len.min(31), // Use k for containment counting
            config.max_seed_freq,
        );

        // Stage 4: Direct alignment against top candidates
        // Extracts actual genome sequences and aligns directly, bypassing
        // the lossy unitig→path→coordinate mapping.
        let mut alignments = direct_align::direct_align_candidates(
            &query.seq,
            &query.name,
            &containment_hits,
            &path_index,
            &unitigs,
            config.max_target_seqs,
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

/// Search with overlay support: searches base index + all overlays, merges results.
///
/// If no overlays exist, behaves identically to `search()`.
pub fn search_with_overlays(
    query_file: &Path,
    config: &SearchConfig,
) -> Result<Vec<QueryResult>> {
    // Search the base index
    let mut results = search(query_file, config)?;

    // Check for overlays
    let manifest = crate::index::update::OverlayManifest::load_or_create(&config.index_dir);
    let manifest = match manifest {
        Ok(m) if m.overlay_count > 0 => m,
        _ => return Ok(results), // No overlays, return base results
    };

    log::info!(
        "Searching {} overlay(s) ({} additional genomes)",
        manifest.overlay_count,
        manifest.overlay_genomes
    );

    // Search each overlay
    for overlay_dir in manifest.overlay_dirs(&config.index_dir) {
        if !overlay_dir.exists() {
            log::warn!("Overlay directory {:?} not found, skipping", overlay_dir);
            continue;
        }

        let overlay_config = SearchConfig {
            index_dir: overlay_dir.into(),
            ..config.clone()
        };

        match search(query_file, &overlay_config) {
            Ok(overlay_results) => {
                // Merge overlay results into base results
                for (base_result, overlay_result) in results.iter_mut().zip(overlay_results) {
                    base_result.alignments.extend(overlay_result.alignments);
                }
            }
            Err(e) => {
                log::warn!("Overlay search failed: {}, skipping", e);
            }
        }
    }

    // Re-sort and deduplicate merged results
    for result in &mut results {
        // Sort by alignment score descending
        result.alignments.sort_by(|a, b| {
            let a_score = a.tags.iter()
                .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            let b_score = b.tags.iter()
                .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal)
        });

        // Per-genome deduplication
        let mut seen = HashSet::new();
        result.alignments.retain(|record| seen.insert(record.target_name.clone()));

        // Enforce max targets
        result.alignments.truncate(config.max_target_seqs);
    }

    Ok(results)
}

/// Search a query against multiple independent indices and merge results.
///
/// This enables distributed/sharded index architecture: when a genome
/// collection is too large to fit in a single index, it can be split into
/// N independent shards. Queries are searched against each shard sequentially
/// (memory-bounded; one index loaded at a time) and results are merged by
/// score with per-genome deduplication.
///
/// Each shard can itself have overlays (via `search_with_overlays`), so
/// hierarchical architectures are supported.
pub fn search_multi_index(
    query_file: &Path,
    index_dirs: &[std::path::PathBuf],
    base_config: &SearchConfig,
) -> Result<Vec<QueryResult>> {
    if index_dirs.is_empty() {
        anyhow::bail!("No index directories provided");
    }
    if index_dirs.len() == 1 {
        // Single index: use overlay search directly
        let mut config = base_config.clone();
        config.index_dir = index_dirs[0].clone().into();
        return search_with_overlays(query_file, &config);
    }

    log::info!("Searching across {} shard indices", index_dirs.len());

    // Merged results by query name
    let mut merged: std::collections::HashMap<String, QueryResult> =
        std::collections::HashMap::new();

    for (i, index_dir) in index_dirs.iter().enumerate() {
        log::info!("  Shard {}/{}: {:?}", i + 1, index_dirs.len(), index_dir);
        if !index_dir.exists() {
            log::warn!("  Shard {:?} not found, skipping", index_dir);
            continue;
        }

        let mut shard_config = base_config.clone();
        shard_config.index_dir = index_dir.clone().into();

        let shard_results = match search_with_overlays(query_file, &shard_config) {
            Ok(r) => r,
            Err(e) => {
                log::warn!("  Shard {:?} failed: {}, skipping", index_dir, e);
                continue;
            }
        };

        for qr in shard_results {
            merged
                .entry(qr.query_name.clone())
                .and_modify(|existing| existing.alignments.extend(qr.alignments.clone()))
                .or_insert(qr);
        }
    }

    // Post-process: sort by score + dedupe + truncate
    let mut results: Vec<QueryResult> = merged.into_values().collect();
    results.sort_by(|a, b| a.query_name.cmp(&b.query_name));

    for result in &mut results {
        result.alignments.sort_by(|a, b| {
            let a_score = a.tags.iter()
                .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            let b_score = b.tags.iter()
                .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut seen = HashSet::new();
        result.alignments.retain(|record| seen.insert(record.target_name.clone()));
        result.alignments.truncate(base_config.max_target_seqs);
    }

    log::info!(
        "Multi-index search complete: {} queries across {} shards",
        results.len(), index_dirs.len()
    );
    Ok(results)
}
