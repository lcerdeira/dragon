pub mod adaptive_kmer;
pub mod bayes;
pub mod candidate;
pub mod chain;
pub mod containment;
pub mod direct_align;
pub mod graph_align;
pub mod ml_score;
pub mod seed;
pub mod spaced_seed;
pub mod sprt;

use anyhow::Result;
use std::collections::HashSet;
use std::path::Path;

use crate::io::paf::PafRecord;

// ─── KmerCache ────────────────────────────────────────────────────────────────

/// Pre-computed FM-index search results AND unitig colour sets for every
/// unique k-mer (and unitig) sampled from a batch of queries.
///
/// Built once per shard before the parallel query loop, then shared read-only
/// across all query threads.
///
/// ## Two layers of caching
///
/// 1. **k-mer → SA positions** (`kmers`): avoids repeated FM-index backward
///    searches for the same k-mer across queries.  Keys with empty SA interval
///    or SA count > `max_seed_freq` are not stored.
///
/// 2. **unitig_id → RoaringBitmap** (`unitig_colors`): avoids repeated
///    `color_index.get_colors()` calls, each of which deserialises a bitmap
///    from the mmap'd file.  For AMR gene panels every unitig in a conserved
///    resistance domain is queried by tens to hundreds of different genes.
///    Pre-fetching all colours once (in parallel) and sharing the results
///    eliminates this O(queries × unitigs_per_query) mmap cost.
pub struct KmerCache {
    /// k-mer bytes → SA positions (filtered by max_seed_freq)
    kmers: std::collections::HashMap<Box<[u8]>, Vec<usize>>,
    /// unitig_id → pre-deserialised RoaringBitmap colour set
    pub unitig_colors: std::collections::HashMap<u32, roaring::RoaringBitmap>,
    /// IDF weight per unitig: log(N_genomes / cardinality). Higher = more discriminative.
    pub unitig_idf: std::collections::HashMap<u32, f32>,
    /// Combined centrality score per unitig: IDF(u) × ln(1 + degree(u)).
    /// Higher = rarer AND at more structural junctions = most discriminative.
    /// Empty when PathIndex is not v4 (degree not available).
    pub unitig_centrality: std::collections::HashMap<u32, f32>,
}

impl KmerCache {
    /// Build both caches for `queries` against `fm` and `colors`.
    ///
    /// Phase 1 — k-mer collection + FM search (parallel over unique k-mers).
    /// Phase 2 — unitig ID extraction from SA positions.
    /// Phase 3 — colour pre-fetch (parallel over unique unitig IDs).
    /// Phase 4 — IDF weight computation.
    /// Phase 5 — graph centrality scores (IDF × ln(1 + degree)).
    pub fn build(
        fm: &crate::index::fm::DragonFmIndex,
        colors: &crate::index::color::ColorIndex,
        path_index: &crate::index::paths::PathIndex,   // NEW parameter
        queries: &[crate::io::fasta::Sequence],
        kmer_size: usize,
        max_seed_freq: usize,
    ) -> Self {
        use rayon::prelude::*;

        const TARGET_SAMPLES: usize = 384;

        // ── Phase 1: collect unique k-mers ────────────────────────────────────
        let mut unique_kmers: std::collections::HashSet<Box<[u8]>> =
            std::collections::HashSet::new();

        for q in queries {
            if q.seq.len() < kmer_size {
                continue;
            }
            let total = q.seq.len() - kmer_size + 1;
            let stride = (total / TARGET_SAMPLES).max(1);
            let mut p = 0usize;
            while p + kmer_size <= q.seq.len() {
                let fwd = &q.seq[p..p + kmer_size];
                let has_ambig = fwd.iter().any(|&b| {
                    !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
                });
                if !has_ambig {
                    unique_kmers.insert(Box::from(fwd));
                    unique_kmers.insert(
                        containment::reverse_complement(fwd).into_boxed_slice(),
                    );
                }
                p += stride;
            }
        }

        let n_kmers_total = unique_kmers.len();

        // Parallel FM-index search → keep only (kmer, positions) pairs that
        // satisfy the frequency filter.
        let kmer_hits: Vec<(Box<[u8]>, Vec<usize>)> = unique_kmers
            .into_par_iter()
            .filter_map(|kmer| {
                let pos = fm.search(&kmer);
                if pos.is_empty() || pos.len() > max_seed_freq {
                    None
                } else {
                    Some((kmer, pos))
                }
            })
            .collect();

        log::info!(
            "KmerCache phase 1: {}/{} k-mers pass freq filter ({:.0}%)",
            kmer_hits.len(),
            n_kmers_total,
            100.0 * kmer_hits.len() as f64 / n_kmers_total.max(1) as f64,
        );

        // ── Phase 2: extract unique unitig IDs from SA positions ──────────────
        let mut unique_unitig_ids: std::collections::HashSet<u32> =
            std::collections::HashSet::new();

        for (_, positions) in &kmer_hits {
            for &pos in positions {
                if let Some((uid, _)) = fm.position_to_unitig(pos) {
                    unique_unitig_ids.insert(uid);
                }
            }
        }

        let n_unitigs = unique_unitig_ids.len();

        // ── Phase 3: parallel colour pre-fetch ────────────────────────────────
        // Each unitig's RoaringBitmap is deserialised exactly once here;
        // all queries that hit this unitig share the result.
        let unitig_colors: std::collections::HashMap<u32, roaring::RoaringBitmap> =
            unique_unitig_ids
                .into_par_iter()
                .filter_map(|uid| {
                    colors
                        .get_colors(uid)
                        .ok()
                        .map(|bm| (uid, bm))
                })
                .collect();

        log::info!(
            "KmerCache phase 3: {}/{} unitig colours pre-fetched",
            unitig_colors.len(),
            n_unitigs,
        );

        // ── Phase 4: compute IDF weights ─────────────────────────────────────
        // IDF(u) = ln(N_genomes / cardinality(u)) — higher means rarer = more discriminative.
        // N_genomes from the color index (already loaded).
        let n_genomes = colors.num_genomes() as f64;
        let unitig_idf: std::collections::HashMap<u32, f32> = unitig_colors
            .iter()
            .map(|(&uid, bm)| {
                let card = bm.len().max(1) as f64;
                let idf = (n_genomes / card).ln() as f32;
                (uid, idf.max(0.0))
            })
            .collect();
        log::info!("KmerCache phase 4: IDF computed for {} unitigs", unitig_idf.len());

        // ── Phase 5: graph centrality scores ─────────────────────────────────────
        // centrality(u) = IDF(u) × ln(1 + degree(u))
        // degree(u) = number of distinct successors in cDBG CSR table (paths.bin v4)
        // Unitigs at structural junctions AND rare → most discriminative for search.
        let unitig_centrality: std::collections::HashMap<u32, f32> = unitig_idf
            .iter()
            .map(|(&uid, &idf)| {
                let degree = path_index.unitig_successor_degree(uid);
                let centrality = idf * (1.0_f64 + degree as f64).ln() as f32;
                (uid, centrality.max(0.0))
            })
            .collect();
        let n_with_centrality = unitig_centrality.values().filter(|&&c| c > 0.0).count();
        log::info!(
            "KmerCache phase 5: centrality computed for {} unitigs ({} non-zero)",
            unitig_centrality.len(), n_with_centrality
        );

        Self {
            kmers: kmer_hits.into_iter().collect(),
            unitig_colors,
            unitig_idf,
            unitig_centrality,
        }
    }

    /// Returns cached SA positions for `kmer`, or `None` if absent/high-freq.
    #[inline]
    pub fn get(&self, kmer: &[u8]) -> Option<&Vec<usize>> {
        self.kmers.get(kmer)
    }

    /// Returns the pre-fetched colour bitmap for `unitig_id`, or `None`.
    #[inline]
    pub fn get_colors(&self, unitig_id: u32) -> Option<&roaring::RoaringBitmap> {
        self.unitig_colors.get(&unitig_id)
    }

    /// Returns the IDF weight for `unitig_id`, or 0.0 if absent.
    #[inline]
    pub fn get_idf(&self, unitig_id: u32) -> f32 {
        self.unitig_idf.get(&unitig_id).copied().unwrap_or(0.0)
    }

    /// Returns centrality score for `unitig_id`, falling back to IDF if centrality unavailable.
    #[inline]
    pub fn get_centrality(&self, unitig_id: u32) -> f32 {
        self.unitig_centrality
            .get(&unitig_id)
            .copied()
            .unwrap_or_else(|| self.get_idf(unitig_id))
    }

    pub fn len(&self) -> usize {
        self.kmers.len()
    }
}

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
    /// Enable batch query mode: pre-compute a per-shard KmerCache across all
    /// queries before searching, then process queries in parallel with rayon.
    /// Default: true.  Disable for single-query interactive use or when
    /// `dump_seeds_path` is set (which requires sequential processing).
    pub batch_queries: bool,
    /// Enable parallel shard loading: search multiple shards concurrently.
    /// Each shard loads its own index into memory; safe with mmap (OS demand-
    /// pages only what's accessed).  Default: true.
    pub parallel_shards: bool,
    /// Cross-species mode: use k=7-8 anchors for 15-30% divergence.
    /// Enables detection of homologs at 70-85% ANI (T2 cross-species tier).
    /// Default: false (within-species, short-read mode).
    pub cross_species: bool,
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
            batch_queries: true,
            parallel_shards: true,
            cross_species: false,
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

    /// Read VmRSS in MB from /proc/self/status. Linux-only; returns None elsewhere.
    fn rss_mb() -> Option<u64> {
        let s = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in s.lines() {
            if let Some(rest) = line.strip_prefix("VmRSS:") {
                let kb: u64 = rest.trim().split_whitespace().next()?.parse().ok()?;
                return Some(kb / 1024);
            }
        }
        None
    }
    macro_rules! log_rss {
        ($stage:expr) => {
            if let Some(mb) = rss_mb() {
                log::info!("[mem] {} RSS={} MB", $stage, mb);
            }
        };
    }
    log_rss!("after-load");

    // ── Batch mode: build KmerCache once for all queries, then run in parallel ─
    //
    // When batch_queries=true (default) and there is no seed-dump writer, we
    // pre-compute all FM-index lookups for unique k-mers across every query,
    // then process queries in parallel with rayon.  Each query thread reads the
    // shared KmerCache (immutable HashMap) and the shared index structures
    // (all Sync via Mmap / Vec).
    //
    // When batch_queries=false or seed_dump is Some, we fall back to the
    // original sequential loop so the mutable seed-dump writer is safe.

    let use_batch = config.batch_queries && seed_dump.is_none() && queries.len() > 1;

    let kmer_cache: Option<KmerCache> = if use_batch {
        log::info!(
            "Batch mode: building KmerCache for {} queries ...",
            queries.len()
        );
        Some(KmerCache::build(
            &fm_index,
            &color_index,
            &path_index,    // NEW
            &queries,
            metadata.kmer_size,
            config.max_seed_freq,
        ))
    } else {
        None
    };

    // ── Inner function: process one query against the already-loaded shard ────
    // Extracted so both the parallel and sequential paths share logic.
    let process_one = |query: &crate::io::fasta::Sequence,
                       cache: Option<&KmerCache>|
     -> QueryResult {
        // Stage 3: containment ranking with adaptive k-mer size.
        //
        // Dragon's FM-index is built at k=31, but fixed k=31 gives
        // P(exact seed) = (1-d)^31 ≈ 0.20 at d=5% divergence — too low for
        // reliable candidate identification.  We select the largest k such
        // that E[exact seeds] ≥ 50 at assumed 5% divergence, then run a
        // multi-k fallback if initial containment finds no candidates.
        let primary_k = adaptive_kmer::adaptive_k(query.seq.len(), metadata.kmer_size);

        let containment_hits = {
            // First attempt at primary (adaptive) k.
            let hits = containment::containment_rank(
                &query.seq,
                &fm_index,
                &color_index,
                primary_k,
                config.max_seed_freq,
                cache,
                config.cross_species,
            );

            if !hits.is_empty() {
                hits
            } else {
                // Multi-k fallback: try progressively smaller k values.
                // Each step increases expected seeds by ~sensitivity_gain(k, k-4).
                let fallback_ks = adaptive_kmer::fallback_k_sequence(primary_k);
                let mut best = Vec::new();
                for &k in fallback_ks.iter().skip(1) { // skip primary_k (already tried)
                    let gain = adaptive_kmer::sensitivity_gain(primary_k, k);
                    log::debug!(
                        "containment fallback k={} (primary={}, gain={:.1}×)",
                        k, primary_k, gain
                    );
                    let h = containment::containment_rank(
                        &query.seq,
                        &fm_index,
                        &color_index,
                        k,
                        config.max_seed_freq,
                        None, // cache built for primary_k; smaller k needs live lookups
                        config.cross_species,
                    );
                    if !h.is_empty() {
                        best = h;
                        break;
                    }
                }
                best
            }
        };

        // Stage 4: direct alignment against top candidates
        let mut alignments = direct_align::direct_align_candidates(
            &query.seq,
            &query.name,
            &containment_hits,
            &path_index,
            &unitigs,
            config.max_target_seqs,
        );

        // Stage 5: post-alignment filters
        alignments.retain(|r| r.identity() >= config.min_identity);
        alignments.retain(|r| {
            if query.seq.is_empty() {
                return false;
            }
            (r.query_end - r.query_start) as f64 / query.seq.len() as f64
                >= config.min_query_coverage
        });
        if !alignments.is_empty() && config.min_score_ratio > 0.0 {
            let best = alignments[0]
                .tags
                .iter()
                .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            let thr = best * config.min_score_ratio;
            alignments.retain(|r| {
                r.tags
                    .iter()
                    .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
                    .unwrap_or(0.0)
                    >= thr
            });
        }
        // Dedup per target genome + cap
        let mut seen: HashSet<String> = HashSet::new();
        alignments.retain(|r| seen.insert(r.target_name.clone()));
        alignments.truncate(config.max_target_seqs);

        QueryResult {
            query_name: query.name.clone(),
            query_len: query.seq.len(),
            alignments,
        }
    };

    // ── Execute: parallel (batch) or sequential ───────────────────────────────
    let results: Vec<QueryResult> = if use_batch {
        use rayon::prelude::*;
        let cache_ref = kmer_cache.as_ref();
        queries
            .par_iter()
            .map(|q| process_one(q, cache_ref))
            .collect()
    } else {
        // Sequential path — supports seed_dump writer
        let mut out = Vec::with_capacity(queries.len());
        for query in &queries {
            // Stage 1 (for seed dump only)
            let seeds_for_dump = if seed_dump.is_some() {
                Some(seed::find_seeds(
                    &query.seq,
                    &fm_index,
                    config.min_seed_len,
                    config.max_seed_freq,
                ))
            } else {
                None
            };
            // Stage 2: optional ML seed dump
            if let (Some(ref mut writer), Some(ref seeds)) = (seed_dump.as_mut(), &seeds_for_dump) {
                use std::io::Write;
                let base_min_votes = if query.seq.len() < 300 { 2u32 } else { 10u32 };
                let min_votes = if seeds.is_empty() {
                    base_min_votes
                } else {
                    base_min_votes.min((seeds.len() as u32 / 5).max(2))
                };
                let candidates = candidate::find_candidates(seeds, &color_index, min_votes);
                let all_positions: Vec<usize> = candidates
                    .iter()
                    .flat_map(|c| c.seeds.iter().map(|s| s.query_pos))
                    .collect();
                for cand in &candidates {
                    let genome_name = metadata
                        .genome_names
                        .get(cand.genome_id as usize)
                        .cloned()
                        .unwrap_or_else(|| format!("genome_{}", cand.genome_id));
                    let label = match &config.ground_truth_genome {
                        Some(gt) => {
                            if genome_name.contains(gt.as_str()) { 1i8 } else { 0i8 }
                        }
                        None => -1i8,
                    };
                    for s in &cand.seeds {
                        let features = ml_score::SeedFeatures::from_seed_with_context(
                            s,
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
            // Stages 3-5 (same logic as batch path, no cache)
            out.push(process_one(query, None));
        }
        out
    };

    log_rss!("after-direct_align");
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

    log::info!(
        "Searching across {} shard indices (parallel_shards={})",
        index_dirs.len(),
        base_config.parallel_shards,
    );

    // ── Per-shard search: parallel or sequential ──────────────────────────────
    //
    // Parallel shard mode: rayon work-steals over shards; each shard loads its
    // own mmap'd index.  Memory cost per shard is small in practice because the
    // OS demand-pages only accessed pages (FM-index backward search touches the
    // suffix array sequentially; colors are accessed sparsely).
    //
    // Sequential mode: load one shard at a time — lower peak RSS.

    let all_shard_results: Vec<Vec<QueryResult>> = if base_config.parallel_shards
        && index_dirs.len() > 1
    {
        use rayon::prelude::*;
        index_dirs
            .par_iter()
            .enumerate()
            .filter_map(|(i, index_dir)| {
                if !index_dir.exists() {
                    log::warn!("  Shard {:?} not found, skipping", index_dir);
                    return None;
                }
                log::info!("  [parallel] Shard {}/{}: {:?}", i + 1, index_dirs.len(), index_dir);
                let mut shard_cfg = base_config.clone();
                shard_cfg.index_dir = index_dir.clone().into();
                match search_with_overlays(query_file, &shard_cfg) {
                    Ok(r) => Some(r),
                    Err(e) => {
                        log::warn!("  Shard {:?} failed: {}, skipping", index_dir, e);
                        None
                    }
                }
            })
            .collect()
    } else {
        let mut out = Vec::with_capacity(index_dirs.len());
        for (i, index_dir) in index_dirs.iter().enumerate() {
            log::info!("  Shard {}/{}: {:?}", i + 1, index_dirs.len(), index_dir);
            if !index_dir.exists() {
                log::warn!("  Shard {:?} not found, skipping", index_dir);
                continue;
            }
            let mut shard_cfg = base_config.clone();
            shard_cfg.index_dir = index_dir.clone().into();
            match search_with_overlays(query_file, &shard_cfg) {
                Ok(r) => out.push(r),
                Err(e) => log::warn!("  Shard {:?} failed: {}, skipping", index_dir, e),
            }
        }
        out
    };

    // Merge all shard results by query name
    let mut merged: std::collections::HashMap<String, QueryResult> =
        std::collections::HashMap::new();

    for shard_results in all_shard_results {
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
