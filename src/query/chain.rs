/// Colinear chaining algorithm with gap-sensitive scoring.
///
/// Uses a Fenwick tree for O(h log h) dynamic programming where
/// h is the number of seed anchors per candidate genome.

use crate::ds::fenwick::FenwickMax;
use crate::index::color::ColorIndex;
use crate::index::fm::SeedHit;
use crate::index::paths::PathIndex;
use crate::query::candidate::Candidate;
use crate::query::ml_score::{SeedFeatures, SeedScorer};

/// An anchor for chaining: a seed hit mapped to genome coordinates.
#[derive(Clone, Debug)]
pub struct Anchor {
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: u64,
    pub ref_end: u64,
    pub match_len: usize,
    pub is_reverse: bool,
    /// ML quality score (weighted by match_len). Used as DP weight.
    pub score: f64,
}

/// A chain of colinear anchors.
#[derive(Clone, Debug)]
pub struct Chain {
    pub genome_id: u32,
    pub score: f64,
    pub anchors: Vec<Anchor>,
    pub query_coverage: f64,
    pub is_reverse: bool,
}

/// Chain seeds for all candidate genomes.
pub fn chain_candidates(
    seeds: &[SeedHit],
    candidates: &[Candidate],
    path_index: &PathIndex,
    min_score: f64,
    scorer: Option<&SeedScorer>,
    color_index: &ColorIndex,
    query_bytes: &[u8],
) -> Vec<Chain> {
    chain_candidates_with_query_len(seeds, candidates, path_index, min_score, 0, scorer, color_index, query_bytes)
}

/// Chain seeds for all candidate genomes, computing coverage relative to full query length.
pub fn chain_candidates_with_query_len(
    seeds: &[SeedHit],
    candidates: &[Candidate],
    path_index: &PathIndex,
    min_score: f64,
    query_len: usize,
    scorer: Option<&SeedScorer>,
    color_index: &ColorIndex,
    query_bytes: &[u8],
) -> Vec<Chain> {
    let mut chains = Vec::new();
    let _ = seeds; // Seeds are accessed via candidates

    for candidate in candidates {
        if let Some(path) = path_index.get_path(candidate.genome_id) {
            // Map seed hits to genome coordinates using the path.
            // Only include seeds whose unitig is found in the path.
            let mut anchors: Vec<Anchor> = candidate
                .seeds
                .iter()
                .filter_map(|seed| {
                    // Look up the unitig in this genome's path
                    let step = path.steps.iter().find(|s| s.unitig_id == seed.unitig_id)?;
                    let genome_offset = step.genome_offset + seed.offset as u64;

                    // Compute ML score if scorer is provided
                    let anchor_score = if let Some(s) = scorer {
                        let features = SeedFeatures::from_seed(
                            seed,
                            query_len.max(query_bytes.len()),
                            query_bytes,
                            color_index,
                        );
                        s.weighted_score(&features)
                    } else {
                        seed.match_len as f64
                    };

                    Some(Anchor {
                        query_start: seed.query_pos,
                        query_end: seed.query_pos + seed.match_len,
                        ref_start: genome_offset,
                        ref_end: genome_offset + seed.match_len as u64,
                        match_len: seed.match_len,
                        is_reverse: seed.is_reverse,
                        score: anchor_score,
                    })
                })
                .collect();

            // Separate forward and reverse anchors
            let fwd_anchors: Vec<_> = anchors.iter().filter(|a| !a.is_reverse).cloned().collect();
            let rev_anchors: Vec<_> = anchors.iter().filter(|a| a.is_reverse).cloned().collect();

            // Chain forward anchors
            if let Some(chain) = compute_chain(candidate.genome_id, &fwd_anchors, false, query_len) {
                if chain.score >= min_score {
                    chains.push(chain);
                }
            }

            // Chain reverse anchors
            if let Some(chain) = compute_chain(candidate.genome_id, &rev_anchors, true, query_len) {
                if chain.score >= min_score {
                    chains.push(chain);
                }
            }

            anchors.clear();
        }
    }

    // Sort chains by total matching bases (containment) descending.
    // This ranks genomes by how many query bases are covered by exact-match
    // seeds, which is more discriminative than DP score for closely related genomes.
    chains.sort_by(|a, b| {
        let a_match: usize = a.anchors.iter().map(|anc| anc.match_len).sum();
        let b_match: usize = b.anchors.iter().map(|anc| anc.match_len).sum();
        b_match.cmp(&a_match)
    });

    chains
}

/// Compute the best colinear chain from a set of anchors using DP.
///
/// For small anchor sets (< 5000), uses O(h²) pairwise DP with proper gap penalty.
/// For larger sets, uses Fenwick tree for O(h log h) without gap penalty.
///
/// Gap cost between consecutive anchors:
///   gap_cost = 0.01 * (query_gap + ref_gap) + 0.5 * |query_gap - ref_gap|
///
/// `query_len` is the full length of the query sequence — if provided (>0),
/// query_coverage is computed as matched bases / query_len rather than
/// matched bases / chain span.
fn compute_chain(genome_id: u32, anchors: &[Anchor], is_reverse: bool, query_len: usize) -> Option<Chain> {
    if anchors.is_empty() {
        return None;
    }

    let mut sorted = anchors.to_vec();
    sorted.sort_by_key(|a| (a.ref_start, a.query_start));

    let n = sorted.len();

    let mut dp = vec![0i64; n];
    let mut parent = vec![usize::MAX; n];

    if n < 5000 {
        // O(h²) pairwise DP with gap penalty — accurate for small anchor sets
        for i in 0..n {
            dp[i] = (sorted[i].score * 1000.0) as i64;
            parent[i] = usize::MAX;

            for j in 0..i {
                // Colinearity check: both query and ref must be increasing
                if sorted[j].query_end > sorted[i].query_start {
                    continue;
                }
                if sorted[j].ref_end > sorted[i].ref_start {
                    continue;
                }

                let query_gap = (sorted[i].query_start - sorted[j].query_end) as i64;
                let ref_gap = (sorted[i].ref_start - sorted[j].ref_end) as i64;

                // Gap cost: linear in total gap + penalty for gap size mismatch (indels)
                let gap_cost = (0.01 * (query_gap + ref_gap) as f64
                    + 0.5 * (query_gap - ref_gap).unsigned_abs() as f64) as i64;

                let chain_score = dp[j] + (sorted[i].score * 1000.0) as i64 - gap_cost;
                if chain_score > dp[i] {
                    dp[i] = chain_score;
                    parent[i] = j;
                }
            }
        }
    } else {
        // Fenwick tree DP for large anchor sets (no gap penalty, but fast)
        let mut query_positions: Vec<usize> = sorted.iter().map(|a| a.query_end).collect();
        query_positions.sort_unstable();
        query_positions.dedup();

        let compress = |pos: usize| -> usize {
            query_positions.binary_search(&pos).unwrap_or_else(|i| i)
        };

        let mut fenwick = FenwickMax::new(query_positions.len() + 1);

        for i in 0..n {
            let match_score = (sorted[i].score * 1000.0) as i64;
            dp[i] = match_score;
            parent[i] = usize::MAX;

            let max_query_idx = compress(sorted[i].query_start);
            if max_query_idx > 0 {
                let best_prev = fenwick.prefix_max(max_query_idx - 1);
                if best_prev > i64::MIN {
                    let chain_score = best_prev + match_score;
                    if chain_score > dp[i] {
                        dp[i] = chain_score;
                        for j in (0..i).rev() {
                            if dp[j] == best_prev
                                && compress(sorted[j].query_end) < max_query_idx
                            {
                                parent[i] = j;
                                break;
                            }
                        }
                    }
                }
            }

            let comp_end = compress(sorted[i].query_end);
            if dp[i] > fenwick.prefix_max(comp_end) {
                fenwick.update(comp_end, dp[i]);
            }
        }
    }

    // Find best chain endpoint
    let (best_idx, &best_score) = dp.iter().enumerate().max_by_key(|(_, s)| *s)?;

    // Traceback using parent pointers
    let mut chain_anchors = Vec::new();
    let mut current = best_idx;
    while current != usize::MAX {
        chain_anchors.push(sorted[current].clone());
        current = parent[current];
    }

    chain_anchors.reverse();

    // Compute query coverage relative to full query length (if known) or chain span
    let total_match: usize = chain_anchors.iter().map(|a| a.match_len).sum();
    let coverage_denominator = if query_len > 0 {
        query_len
    } else if let (Some(first), Some(last)) = (chain_anchors.first(), chain_anchors.last()) {
        (last.query_end - first.query_start).max(1)
    } else {
        1
    };
    let coverage = total_match as f64 / coverage_denominator as f64;

    Some(Chain {
        genome_id,
        score: best_score as f64,
        anchors: chain_anchors,
        query_coverage: coverage,
        is_reverse,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_colinear_anchors() {
        let anchors = vec![
            Anchor {
                query_start: 0,
                query_end: 31,
                ref_start: 100,
                ref_end: 131,
                match_len: 31,
                is_reverse: false,
                score: 31.0,
            },
            Anchor {
                query_start: 50,
                query_end: 81,
                ref_start: 150,
                ref_end: 181,
                match_len: 31,
                is_reverse: false,
                score: 31.0,
            },
            Anchor {
                query_start: 100,
                query_end: 131,
                ref_start: 200,
                ref_end: 231,
                match_len: 31,
                is_reverse: false,
                score: 31.0,
            },
        ];

        let chain = compute_chain(0, &anchors, false, 200).unwrap();
        assert!(chain.score > 0.0);
        assert!(!chain.anchors.is_empty());
    }

    #[test]
    fn test_empty_anchors() {
        let chain = compute_chain(0, &[], false, 0);
        assert!(chain.is_none());
    }
}
