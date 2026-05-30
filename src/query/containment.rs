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

    let mut p = 0usize;
    while p + kmer_size <= query.len() {
        n_sampled += 1;
        let fwd = &query[p..p + kmer_size];
        let has_ambig = fwd
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));
        if !has_ambig {
            // Search the k-mer on both strands (forward first, so the
            // forward hit marks the position before the reverse one).
            let rc = reverse_complement(fwd);
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
        }
        p += stride;
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
    const CANDIDATE_TARGET: u64 = 2000;
    let mut by_card: Vec<(u32, u64)> = unitig_colors
        .iter()
        .map(|(uid, c)| (*uid, c.len()))
        .collect();
    by_card.sort_by_key(|&(_, card)| card);

    let mut candidate_set = RoaringBitmap::new();
    for (uid, _card) in &by_card {
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
    const MAX_SEEDS_PER_GENOME: usize = 64;
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
                }
                *genome_info.entry(genome_id).or_insert(0.0) += ic;
                let bucket = genome_seeds.entry(genome_id).or_default();
                if bucket.len() < MAX_SEEDS_PER_GENOME {
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
        .map(|(genome_id, shared)| ContainmentHit {
            genome_id,
            // containment estimated over the sampled k-mers
            containment: shared as f64 / n_sampled.max(1) as f64,
            shared_kmers: shared,
            total_query_kmers: n_sampled,
            info_score: genome_info.get(&genome_id).copied().unwrap_or(0.0),
            seeds: genome_seeds.remove(&genome_id).unwrap_or_default(),
        })
        .collect();

    // Sort by containment descending, break ties by info_score
    hits.sort_by(|a, b| {
        b.containment.partial_cmp(&a.containment)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.info_score.partial_cmp(&a.info_score)
                .unwrap_or(std::cmp::Ordering::Equal))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACGTACGT"), b"ACGTACGT");
    }
}
