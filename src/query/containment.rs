/// Containment-based genome ranking.
///
/// Computes k-mer containment between a query and each candidate genome:
///   containment(Q, G) = |kmers(Q) ∩ kmers(G)| / |kmers(Q)|
///
/// This is mathematically related to ANI and provides near-perfect sensitivity
/// for queries with ANI > 90% to the target genome. It bypasses the chaining
/// pipeline entirely, using the color index to check genome membership directly.
///
/// For each query k-mer that matches a unitig in the FM-index, we check which
/// genomes contain that unitig via the color index. A genome's containment score
/// is the fraction of query k-mers found in it.

use std::collections::HashMap;

use roaring::RoaringBitmap;

use crate::index::color::ColorIndex;
use crate::index::fm::{DragonFmIndex, SeedHit};
use crate::query::KmerCache;
use crate::query::bayes::{bayesian_probs, bayesian_ani};
use crate::query::spaced_seed::{AnchorConfig, pigeonhole_search};
use crate::query::sprt::{SprtDecision, SprtState};

/// A genome ranked by k-mer containment.
#[derive(Clone, Debug)]
pub struct ContainmentHit {
    pub genome_id: u32,
    /// Fraction of query k-mers found in this genome (0.0 to 1.0).
    pub containment: f64,
    /// Number of query k-mers found in this genome.
    pub shared_kmers: usize,
    /// Total query k-mers searched.
    pub total_query_kmers: usize,
    /// Information-weighted score (sum of IC for matched unitigs).
    pub info_score: f64,
    /// Bayesian posterior P(true containment ≥ 0.5 | observed hits).
    /// Calibrated probability: accounts for sample size, unlike raw containment.
    pub bayes_prob: Option<f64>,
    /// High-confidence Bayesian P(true containment ≥ 0.9 | observed hits).
    /// For clinical/epidemiological applications requiring strong evidence.
    pub bayes_prob_hc: Option<f64>,
    /// Bayesian ANI estimate (posterior mean, Laplace-smoothed).
    pub bayes_ani: Option<f64>,
    /// Seeds that map to this genome (for optional chaining/alignment).
    pub seeds: Vec<SeedHit>,
}

/// Compute containment-based ranking for a query against all genomes.
///
/// Steps:
/// 1. Extract all k-mers from the query at stride 1
/// 2. Search each in the FM-index to find matching unitigs
/// 3. For each matching unitig, look up genome membership via color index
/// 4. Accumulate per-genome k-mer counts
/// 5. Rank by containment = shared_kmers / total_query_kmers
/// Main entry point — uses a pre-built `KmerCache` when available, otherwise
/// falls back to live FM-index lookups.  The cache is built once per shard
/// across all queries and shared read-only, so this function is `Sync`-safe
/// for parallel query processing.
pub fn containment_rank(
    query: &[u8],
    fm_index: &DragonFmIndex,
    color_index: &ColorIndex,
    kmer_size: usize,
    max_seed_freq: usize,
    kmer_cache: Option<&KmerCache>,
    cross_species: bool,
) -> Vec<ContainmentHit> {
    if query.len() < kmer_size {
        return Vec::new();
    }

    let total_genomes = color_index.num_genomes();
    let total_kmers = query.len() - kmer_size + 1;

    // Containment is a RANKING statistic, so we SAMPLE k-mers rather than
    // scan every position. A few hundred sampled k-mers estimate containment
    // as well as all ~N of them, at a fraction of the FM-index lookups and
    // genome-crediting — the dominant per-query cost on large shards. Exact
    // identity is re-derived downstream by direct_align.
    //
    // 384 samples keeps the perfect-hit rate on par with a full scan; a
    // 200-sample setting was faster still but occasionally missed a
    // localized match (perfect-rate 0.975 vs 0.980 on the saureus benchmark).
    const TARGET_SAMPLES: usize = 384;
    let stride = (total_kmers / TARGET_SAMPLES).max(1);

    /// One query position's FM-index hit: the distinct unitigs it maps to,
    /// with a representative seed per unitig.
    struct KmerHit {
        qpos: usize,
        is_reverse: bool,
        unitigs: Vec<u32>,
        seeds: Vec<SeedHit>,
    }

    // ---- Phase 1: collect sampled k-mer hits; deserialize each unitig's
    // color set exactly ONCE for the whole query (not once per k-mer). ----
    let mut kmer_hits: Vec<KmerHit> = Vec::new();
    let mut unitig_colors: HashMap<u32, RoaringBitmap> = HashMap::new();
    let mut n_sampled = 0usize;
    // Query positions where the solid k-mer missed (candidates for the
    // pigeonhole FALLBACK — see the gated pass after the sampling loop).
    let mut pigeon_positions: Vec<usize> = Vec::new();

    // SPRT: one state per containment_rank call — tests whether ANY k-mer hit
    // is observed (query-level: "does this query appear in ANY genome?").
    // Cross-species mode uses the fast (higher error) variant to avoid
    // early rejection when hit rates are low (15-30% divergence).
    let mut sprt = if cross_species {
        SprtState::fast_for_kmer31()   // higher α/β = less conservative
    } else {
        SprtState::default_for_kmer31()
    };
    // Set once the SPRT reaches a decision (H₀ or H₁); we then continue
    // processing normally but skip further SPRT updates.
    let mut sprt_decided = false;

    let mut p = 0usize;
    while p + kmer_size <= query.len() {
        n_sampled += 1;
        let fwd = &query[p..p + kmer_size];
        let has_ambig = fwd
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));
        if !has_ambig {
            // SPRT: treat any k-mer hit (on either strand) as a hit observation.
            // This is a query-level test — does the query appear in ANY genome?
            // We use a single hit/miss per query position (forward strand only,
            // to avoid double-counting the same position).

            // Search the k-mer on both strands (forward first, so the
            // forward hit marks the position before the reverse one).
            let rc = reverse_complement(fwd);
            // SPRT counts a hit if EITHER strand has a match at this position.
            // Using forward-only would miss queries that are entirely on the
            // reverse complement strand (e.g. reverse_strand_perfect_match test).
            let mut position_hit = false;
            for (kmer, is_reverse) in [(fwd.to_vec(), false), (rc, true)] {
                // Use pre-computed cache when available; fall back to live search.
                let positions: std::borrow::Cow<Vec<usize>> = match kmer_cache {
                    Some(cache) => match cache.get(&kmer) {
                        Some(cached) => std::borrow::Cow::Borrowed(cached),
                        None => continue, // cache miss means absent or too frequent
                    },
                    None => {
                        let p = fm_index.search(&kmer);
                        if p.is_empty() || p.len() > max_seed_freq { continue; }
                        std::borrow::Cow::Owned(p)
                    }
                };
                // Record a hit for SPRT if EITHER strand matched.
                if !positions.is_empty() {
                    position_hit = true;
                }
                let mut seen = RoaringBitmap::new();
                let mut unitigs: Vec<u32> = Vec::new();
                let mut seeds: Vec<SeedHit> = Vec::new();
                for &pos in positions.as_ref() {
                    if let Some((unitig_id, offset)) = fm_index.position_to_unitig(pos) {
                        if seen.insert(unitig_id) {
                            unitigs.push(unitig_id);
                            seeds.push(SeedHit {
                                unitig_id,
                                offset,
                                query_pos: p,
                                match_len: kmer_size,
                                is_reverse,
                                sa_count: positions.len(),
                            });
                            unitig_colors.entry(unitig_id).or_insert_with(|| {
                                // Use pre-fetched cache when available (avoids
                                // mmap deserialisation per unitig per query).
                                // Fall back to live color_index lookup for any
                                // unitig that wasn't in the pre-fetch set.
                                if let Some(cache) = kmer_cache {
                                    if let Some(bm) = cache.get_colors(unitig_id) {
                                        return bm.clone();
                                    }
                                }
                                color_index.get_colors(unitig_id).unwrap_or_default()
                            });
                        }
                    }
                }
                if !unitigs.is_empty() {
                    kmer_hits.push(KmerHit { qpos: p, is_reverse, unitigs, seeds });
                }
            }
            // Update SPRT once per query position using the forward-strand result.
            //
            // CRITICAL: AcceptH0 must NOT discard the query. The SPRT's H₁ hit
            // rate (p₁ = 0.98³¹ ≈ 0.53) is calibrated for ~2% divergence; a
            // divergent-but-present query (10–15%) has a solid-k-mer hit rate of
            // only ~3–9%, which the SPRT cannot distinguish from absence, so it
            // would fire AcceptH0 after ~4 early misses and `return` — throwing
            // away a real match (and any pigeonhole anchors already collected).
            // This was the dominant recall ceiling vs LexicMap (77% @ 15% div).
            //
            // Instead, treat EITHER decision as "stop updating the SPRT" and keep
            // sampling so the pigeonhole fallback (shorter k=7–10 anchors) can
            // seed the divergent regions. Truly-absent queries still terminate
            // correctly via the `kmer_hits.is_empty()` check below.
            if !sprt_decided {
                match sprt.update(position_hit) {
                    SprtDecision::AcceptH0 | SprtDecision::AcceptH1 => sprt_decided = true,
                    SprtDecision::Continue => {}
                }
            }

            // Record solid-miss positions for the pigeonhole FALLBACK pass.
            // We do NOT run pigeonhole inline per-position: doing so floods the
            // candidate pool with ~1M low-specificity core anchors on a 1 kb
            // query that has scattered edits (e.g. d3, ~13 edits → ~400 misses ×
            // ~2500 anchors). That noise (a) crowds the true region out of
            // candidate selection and (b) corrupts the per-genome anchor cluster
            // so even genomes that CONTAIN the query align at ~50% (observed:
            // q25 aligned its own 98.7%-identical source at 52%). The decision
            // to spend anchors is made ONCE, after we know the global solid-hit
            // strength — see the gated pass below.
            if !position_hit && kmer_size >= 7 {
                pigeon_positions.push(p);
            }
        }
        p += stride;
    }

    // ---- Pigeonhole FALLBACK (query-level, not per-position) ----
    // Solid 31-mer anchors are precise; pigeonhole (short k=7–10 anchors) is a
    // recall fallback for queries where solid seeding is too sparse to identify
    // candidates (short reads, high divergence). Running it only when solid
    // signal is WEAK keeps it from flooding/ corrupting queries that already
    // have ample solid anchors. `kmer_hits` here holds solid hits only.
    const MIN_SOLID_ANCHORS: usize = 30;
    let solid_weak = kmer_hits.len() < MIN_SOLID_ANCHORS;
    if solid_weak && !pigeon_positions.is_empty() {
        let anchor_cfg = if cross_species {
            crate::query::spaced_seed::AnchorConfig::for_cross_species()
        } else {
            AnchorConfig::default()
        };
        let min_window = anchor_cfg.min_window();
        for &p in &pigeon_positions {
            let window = &query[p..p.saturating_add(kmer_size).min(query.len())];
            if window.len() < min_window {
                continue;
            }
            let (_, anchor_hits) =
                pigeonhole_search(query, p, kmer_size, fm_index, max_seed_freq, &anchor_cfg);
            if anchor_hits.is_empty() {
                continue;
            }
            let mut seen = RoaringBitmap::new();
            let mut unitigs: Vec<u32> = Vec::new();
            let mut seeds: Vec<SeedHit> = Vec::new();
            for hit in &anchor_hits {
                if seen.insert(hit.unitig_id) {
                    unitigs.push(hit.unitig_id);
                    seeds.push(hit.clone());
                    unitig_colors.entry(hit.unitig_id).or_insert_with(|| {
                        if let Some(cache) = kmer_cache {
                            if let Some(bm) = cache.get_colors(hit.unitig_id) {
                                return bm.clone();
                            }
                        }
                        color_index.get_colors(hit.unitig_id).unwrap_or_default()
                    });
                }
            }
            if !unitigs.is_empty() {
                kmer_hits.push(KmerHit { qpos: p, is_reverse: false, unitigs, seeds });
            }
        }
    }

    if kmer_hits.is_empty() {
        return Vec::new();
    }

    // ---- Phase 2: select candidate genomes from the most SPECIFIC unitigs.
    //
    // A unitig shared by a large fraction of genomes (core genome) carries
    // almost no ranking signal, yet crediting it costs O(cardinality) — on
    // the 16K-genome saureus shards that meant crediting ~all genomes for
    // every k-mer, an O(genomes) pass per query and a ~16K-entry result.
    //
    // Instead, build the candidate set greedily from the lowest-cardinality
    // (most specific) unitigs until it is large enough. The true source
    // genome — and any genome that genuinely contains the query — shares the
    // query's most specific k-mers, so it is always captured. ----
    // Soft target for the greedy candidate set. 2000 capped recall on
    // population-scale surveillance: a gene carried by ~3600 genomes left a
    // sub-clade of carriers (sharing slightly less-specific unitigs) ranked
    // just below the cutoff and never aligned, even at 100% identity / full
    // coverage. 8000 captures all carriers of a per-species panel gene with
    // margin while staying bounded (the AND with each unitig's colours keeps
    // per-k-mer work proportional to this, not to the whole collection).
    const CANDIDATE_TARGET: u64 = 8000;
    // When a KmerCache is available, sort by IDF descending (most discriminative first).
    // Otherwise sort by cardinality ascending (rarest first) — same logic as before.
    let mut by_weight: Vec<(u32, f64)> = unitig_colors
        .iter()
        .map(|(uid, c)| {
            let weight = if let Some(cache) = kmer_cache {
                let idf = cache.get_centrality(*uid) as f64;
                // Negate so sort_by ascending puts highest-centrality first.
                if idf > 0.0 { -idf } else { c.len() as f64 }
            } else {
                c.len() as f64
            };
            (*uid, weight)
        })
        .collect();
    by_weight.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut candidate_set = RoaringBitmap::new();
    for (uid, _) in &by_weight {
        if candidate_set.len() >= CANDIDATE_TARGET {
            break;
        }
        if let Some(colors) = unitig_colors.get(uid) {
            candidate_set |= colors;
        }
    }
    if candidate_set.is_empty() {
        return Vec::new();
    }

    // ---- Phase 3: score ONLY candidate genomes. Each hit unitig's color set
    // is intersected with the candidate set (fast Roaring AND), so per-k-mer
    // work is bounded by the candidate set, not the whole genome collection. ----
    // Cap seeds kept per genome. A flat 64 is dense for ~1 kb AMR genes but far
    // too sparse for large genes: a 9.6 kb gene (e.g. clbB) filled the bucket
    // from its first ~1.5 kb and dropped every later seed, so the alignment
    // window missed the rest of the gene (and, for genes split across a contig
    // boundary, anchored on the wrong fragment). Scale ~1 seed per 16 bp of
    // query so seeds span the full gene, with a 64 floor for short queries.
    let max_seeds_per_genome: usize = (query.len() / 16).max(64);
    let mut genome_shared: HashMap<u32, usize> = HashMap::new();
    let mut genome_info: HashMap<u32, f64> = HashMap::new();
    let mut genome_seeds: HashMap<u32, Vec<SeedHit>> = HashMap::new();
    let mut counted_positions = vec![false; query.len()];

    for kh in &kmer_hits {
        for (idx, uid) in kh.unitigs.iter().enumerate() {
            let colors = match unitig_colors.get(uid) {
                Some(c) => c,
                None => continue,
            };
            let credited = colors & &candidate_set;
            if credited.is_empty() {
                continue;
            }
            let cardinality = colors.len();
            let ic = if cardinality > 0 && total_genomes > 0 {
                (total_genomes as f64 / cardinality as f64).log2()
            } else {
                0.0
            };
            let count_this = !counted_positions[kh.qpos];
            let seed = &kh.seeds[idx];
            for genome_id in credited.iter() {
                if count_this {
                    *genome_shared.entry(genome_id).or_insert(0) += 1;
                    // Gate IC with the same per-position guard as genome_shared;
                    // otherwise a position matched on both strands double-adds IC
                    // and skews the (info-score) ranking.
                    *genome_info.entry(genome_id).or_insert(0.0) += ic;
                }
                let bucket = genome_seeds.entry(genome_id).or_default();
                if bucket.len() < max_seeds_per_genome {
                    bucket.push(seed.clone());
                }
            }
        }
        // Mark this query position counted (forward strand runs first).
        if !kh.is_reverse {
            counted_positions[kh.qpos] = true;
        }
    }

    // Build ranked results
    let mut hits: Vec<ContainmentHit> = genome_shared
        .into_iter()
        .map(|(genome_id, shared)| {
            let containment = shared as f64 / n_sampled.max(1) as f64;
            let (bp, bp_hc) = bayesian_probs(shared, n_sampled);
            let bani = bayesian_ani(shared, n_sampled, kmer_size);
            ContainmentHit {
                genome_id,
                containment,
                shared_kmers: shared,
                total_query_kmers: n_sampled,
                info_score: genome_info.get(&genome_id).copied().unwrap_or(0.0),
                bayes_prob: bp,
                bayes_prob_hc: bp_hc,
                bayes_ani: bani,
                seeds: genome_seeds.remove(&genome_id).unwrap_or_default(),
            }
        })
        .collect();

    // Sort by IDF-weighted information content (info_score) descending, breaking
    // ties by Bayesian posterior P(match).
    //
    // WHY info_score FIRST (not raw containment / bayes_prob):
    // At high divergence the TRUE source shares only a few k-mers with the query
    // — but they include the query's RARE, discriminative k-mers (high IDF,
    // ic = log2(N/cardinality)). A WRONG genome that merely shares the query's
    // conserved/core k-mers accumulates a comparable raw shared-count (so a
    // similar bayes_prob) but its IDF contribution is ≈0 (core k-mers have
    // cardinality ≈ N → ic ≈ 0). Ranking by bayes_prob/containment therefore
    // let core-sharing impostors outrank the true source, pushing it below the
    // alignment cutoff (top max_target_seqs) so it was never aligned — the
    // dominant high-divergence recall gap vs LexicMap (d15: 46% vs 100%).
    //
    // info_score down-weights core k-mers to ~0 automatically, so the source's
    // rare-k-mer evidence dominates and it ranks into the aligned set.
    hits.sort_by(|a, b| {
        b.info_score.partial_cmp(&a.info_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                let a_score = a.bayes_prob.unwrap_or(a.containment);
                let b_score = b.bayes_prob.unwrap_or(b.containment);
                b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    hits
}

pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

/// Compute ANI 95% CI half-width using the delta method (Blanca et al. 2022).
pub fn ani_confidence_interval(
    shared_kmers: usize,
    total_kmers: usize,
    kmer_size: usize,
) -> Option<f64> {
    if total_kmers == 0 || kmer_size == 0 { return None; }
    let c = shared_kmers as f64 / total_kmers as f64;
    if c <= 0.0 { return None; }
    let ani = c.powf(1.0 / kmer_size as f64);
    if ani <= 0.0 { return None; }
    let var_c = c * (1.0 - c) / total_kmers as f64;
    let deriv = kmer_size as f64 * ani.powi(kmer_size as i32 - 1);
    if deriv <= f64::EPSILON { return None; }
    let var_ani = var_c / (deriv * deriv);
    Some((1.96 * var_ani.sqrt()).min(0.5)) // cap at 0.5 to avoid nonsensical values
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACGTACGT"), b"ACGTACGT");
    }

    #[test]
    fn ani_ci_zero_total_returns_none() {
        assert!(ani_confidence_interval(0, 0, 31).is_none());
    }

    #[test]
    fn ani_ci_zero_kmer_size_returns_none() {
        assert!(ani_confidence_interval(10, 100, 0).is_none());
    }

    #[test]
    fn ani_ci_zero_shared_returns_none() {
        // c = 0.0 → None
        assert!(ani_confidence_interval(0, 100, 31).is_none());
    }

    #[test]
    fn ani_ci_perfect_containment_small_ci() {
        // 1000/1000 shared k-mers of size 31 → ANI ≈ 1.0, CI should be very small
        let ci = ani_confidence_interval(1000, 1000, 31);
        // c=1.0 → var_c = 0 → CI = 0
        assert!(ci.is_some());
        assert!(ci.unwrap() < 1e-6, "Expected near-zero CI for perfect containment");
    }

    #[test]
    fn ani_ci_realistic_value() {
        // ~53% containment (ANI≈0.98, k=31), 384 sampled k-mers
        let shared = (0.533 * 384.0) as usize;
        let ci = ani_confidence_interval(shared, 384, 31);
        assert!(ci.is_some());
        let ci_val = ci.unwrap();
        // CI should be finite, positive, and well below 0.5
        assert!(ci_val > 0.0 && ci_val < 0.5,
            "Unexpected CI value: {}", ci_val);
    }

    #[test]
    fn ani_ci_capped_at_half() {
        // Very few k-mers → large variance → cap kicks in
        let ci = ani_confidence_interval(1, 2, 31);
        assert!(ci.is_some());
        assert!(ci.unwrap() <= 0.5, "CI must be capped at 0.5");
    }
}
