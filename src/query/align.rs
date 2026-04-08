/// Banded wavefront alignment (WFA) for base-level alignment.
///
/// Uses a simplified Needleman-Wunsch implementation for base-level
/// identity calculation. Production version should use libwfa2 for
/// optimal performance on very long sequences.

use crate::index::paths::PathIndex;
use crate::index::unitig::UnitigSet;
use crate::io::paf::PafRecord;
use crate::query::chain::Chain;

/// Align chains against their reference regions and produce PAF records.
///
/// For each chain:
/// 1. Extract reference sequence from path index + unitigs
/// 2. Run NW alignment for regions < 1M cells, otherwise estimate from seeds
/// 3. Estimate mapq from chain score, coverage, and score gap to second-best
pub fn align_chains(
    query: &[u8],
    chains: &[Chain],
    path_index: &PathIndex,
    unitigs: &UnitigSet,
) -> Vec<PafRecord> {
    let mut records = Vec::new();

    // Track best and second-best scores per genome for mapq calculation
    let best_score = chains.first().map(|c| c.score).unwrap_or(0.0);
    let second_best_score = chains.get(1).map(|c| c.score).unwrap_or(0.0);

    for chain in chains {
        if chain.anchors.is_empty() {
            continue;
        }

        let first = chain.anchors.first().unwrap();
        let last = chain.anchors.last().unwrap();

        // Compute alignment region
        let query_start = first.query_start;
        let query_end = last.query_end;
        let ref_start = first.ref_start;
        let ref_end = last.ref_end;

        let query_span = query_end - query_start;
        let ref_span = (ref_end - ref_start) as usize;

        // Determine if ref region is sane (within 3x of query span).
        // If ref_span >> query_span, the seeds are scattered across the genome
        // and NW alignment would be meaningless — use seed-density estimation.
        let ref_span_ratio = if query_span > 0 { ref_span as f64 / query_span as f64 } else { 100.0 };

        let (num_matches, alignment_len, real_ref_start, real_ref_end) = if query_span > 0 && ref_span_ratio <= 3.0 {
            // Extract the reference region from the path index
            let ref_seq = path_index.extract_sequence(
                chain.genome_id,
                ref_start,
                ref_end,
                unitigs,
            );

            if !ref_seq.is_empty() && query_span * ref_seq.len() < 1_000_000 {
                // Run real NW alignment
                let query_region = &query[query_start..query_end];
                let (matches, alen, _cigar) = banded_nw_align(query_region, &ref_seq, 100);
                (matches, alen, ref_start, ref_end)
            } else {
                let (matches, alen) = estimate_identity_from_seeds(chain, query_span);
                (matches, alen, ref_start, ref_end)
            }
        } else if query_span > 0 {
            // Ref region too spread out — estimate identity from seed anchors
            let (matches, alen) = estimate_identity_from_seeds(chain, query_span);
            (matches, alen, ref_start, ref_end)
        } else {
            (0, 1, ref_start, ref_end)
        };

        // Mapq estimation incorporating score gap
        let mapq = estimate_mapq(chain, best_score, second_best_score);

        let record = PafRecord {
            query_name: String::new(), // Filled in by caller
            query_len: query.len(),
            query_start,
            query_end,
            strand: if chain.is_reverse { '-' } else { '+' },
            target_name: format!("genome_{}", chain.genome_id),
            target_len: 0, // Filled in from metadata
            target_start: real_ref_start as usize,
            target_end: real_ref_end as usize,
            num_matches,
            alignment_len,
            mapq,
            tags: vec![
                format!("AS:i:{}", chain.anchors.iter().map(|a| a.match_len).sum::<usize>()),
                format!("cs:f:{:.4}", chain.query_coverage),
                format!("de:f:{:.4}", 1.0 - num_matches as f64 / alignment_len.max(1) as f64),
            ],
        };

        records.push(record);
    }

    records
}

/// Estimate identity from seed anchors when real alignment is not feasible.
///
/// Uses anchor coverage (fraction of query spanned by exact-match seeds) to estimate
/// overall alignment identity. With compacted unitigs, anchors can be much longer
/// than k, providing a strong lower bound on true identity.
fn estimate_identity_from_seeds(chain: &Chain, query_span: usize) -> (usize, usize) {
    let total_anchor_bases: usize = chain.anchors.iter().map(|a| a.match_len).sum();

    // Seed density = fraction of the query span covered by exact-match anchors
    let density = total_anchor_bases as f64 / query_span.max(1) as f64;

    // With density d of exact matches, the remaining (1-d) fraction has unknown identity.
    // For closely related genomes (typical use case), inter-anchor regions are likely
    // ~90% identical if density > 0.3, tapering to ~50% at low density.
    let inter_anchor_identity = if density >= 0.5 {
        0.95 // very high coverage → nearly identical
    } else if density >= 0.2 {
        0.7 + 0.5 * density // moderate coverage → likely similar
    } else if density > 0.0 {
        0.5 * density / 0.2 // low coverage → uncertain
    } else {
        0.0
    };

    let inter_anchor_bases = query_span.saturating_sub(total_anchor_bases);
    let inter_anchor_matches = (inter_anchor_bases as f64 * inter_anchor_identity) as usize;

    let total_matches = total_anchor_bases + inter_anchor_matches;
    (total_matches, query_span)
}

/// Estimate mapping quality from chain properties and score gap.
///
/// Uses the ratio between best and second-best chain score to estimate
/// confidence. A large gap means high confidence (high mapq); a small
/// gap means the hit is ambiguous.
fn estimate_mapq(chain: &Chain, best_score: f64, second_best_score: f64) -> u8 {
    // Component 1: score magnitude (log-scaled)
    let score_component = if chain.score > 0.0 {
        (chain.score.log2() * 5.0).min(30.0).max(0.0)
    } else {
        0.0
    };

    // Component 2: query coverage (higher = more confident)
    let coverage_component = chain.query_coverage * 20.0;

    // Component 3: score gap to second-best (larger gap = more unique)
    let gap_component = if best_score > 0.0 && chain.score == best_score {
        let ratio = second_best_score / best_score;
        // If second-best is much lower (ratio ~0), mapq bonus is high
        // If second-best is close (ratio ~1), mapq bonus is low
        ((1.0 - ratio) * 20.0).max(0.0)
    } else {
        0.0
    };

    let mapq = (score_component + coverage_component + gap_component).min(60.0).max(0.0);
    mapq as u8
}

/// Simplified banded Needleman-Wunsch alignment.
/// Returns (num_matches, alignment_length, cigar_string).
pub fn banded_nw_align(query: &[u8], reference: &[u8], bandwidth: usize) -> (usize, usize, String) {
    let n = query.len();
    let m = reference.len();

    if n == 0 || m == 0 {
        return (0, n.max(m), String::new());
    }

    // Simple non-banded NW for short sequences
    if n * m < 1_000_000 {
        return simple_nw(query, reference);
    }

    // For longer sequences, use banded approach
    let band = bandwidth.min(n).min(m);
    banded_nw_internal(query, reference, band)
}

fn simple_nw(query: &[u8], reference: &[u8]) -> (usize, usize, String) {
    let n = query.len();
    let m = reference.len();

    let match_score: i32 = 2;
    let mismatch_penalty: i32 = -3;
    let gap_penalty: i32 = -1;

    // DP matrix
    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    for i in 0..=n {
        dp[i][0] = i as i32 * gap_penalty;
    }
    for j in 0..=m {
        dp[0][j] = j as i32 * gap_penalty;
    }

    for i in 1..=n {
        for j in 1..=m {
            let score = if query[i - 1].to_ascii_uppercase() == reference[j - 1].to_ascii_uppercase()
            {
                match_score
            } else {
                mismatch_penalty
            };

            dp[i][j] = (dp[i - 1][j - 1] + score)
                .max(dp[i - 1][j] + gap_penalty)
                .max(dp[i][j - 1] + gap_penalty);
        }
    }

    // Traceback
    let mut matches = 0usize;
    let mut align_len = 0usize;
    let mut i = n;
    let mut j = m;

    while i > 0 && j > 0 {
        let score = if query[i - 1].to_ascii_uppercase() == reference[j - 1].to_ascii_uppercase()
        {
            match_score
        } else {
            mismatch_penalty
        };

        if dp[i][j] == dp[i - 1][j - 1] + score {
            if score == match_score {
                matches += 1;
            }
            i -= 1;
            j -= 1;
        } else if dp[i][j] == dp[i - 1][j] + gap_penalty {
            i -= 1;
        } else {
            j -= 1;
        }
        align_len += 1;
    }
    align_len += i + j;

    (matches, align_len, String::new())
}

fn banded_nw_internal(query: &[u8], reference: &[u8], band: usize) -> (usize, usize, String) {
    let n = query.len();
    let m = reference.len();

    let match_score: i32 = 2;
    let mismatch_penalty: i32 = -3;
    let gap_penalty: i32 = -1;

    // Banded DP: only compute cells within `band` of the diagonal
    let mut dp = vec![vec![i32::MIN; 2 * band + 1]; n + 1];

    // Initialize
    for offset in 0..=band.min(m) {
        dp[0][offset + band] = offset as i32 * gap_penalty;
    }
    for i in 0..=band.min(n) {
        dp[i][band] = i as i32 * gap_penalty;
    }

    for i in 1..=n {
        let j_center = (i as f64 * m as f64 / n as f64) as usize;
        let j_start = j_center.saturating_sub(band).max(1);
        let j_end = (j_center + band).min(m);

        for j in j_start..=j_end {
            let band_offset = j + band - j_center;
            if band_offset >= 2 * band + 1 {
                continue;
            }

            let score = if query[i - 1].to_ascii_uppercase() == reference[j - 1].to_ascii_uppercase()
            {
                match_score
            } else {
                mismatch_penalty
            };

            let mut best = i32::MIN;

            // Diagonal
            if j > j_start || (j == j_start && band_offset > 0) {
                let prev_offset = band_offset - 1;
                if prev_offset < 2 * band + 1 && dp[i - 1][prev_offset] != i32::MIN {
                    best = best.max(dp[i - 1][prev_offset] + score);
                }
            }

            // From above
            if dp[i - 1][band_offset] != i32::MIN {
                best = best.max(dp[i - 1][band_offset] + gap_penalty);
            }

            // From left
            if band_offset + 1 < 2 * band + 1 && dp[i][band_offset + 1] != i32::MIN {
                best = best.max(dp[i][band_offset + 1] + gap_penalty);
            }

            dp[i][band_offset] = best;
        }
    }

    // Traceback from bottom-right of the band
    let j_center_last = (n as f64 * m as f64 / n as f64) as usize;
    let band_offset_end = m + band - j_center_last;
    let final_score = if band_offset_end < 2 * band + 1 {
        dp[n][band_offset_end]
    } else {
        i32::MIN
    };

    // If we couldn't reach the end, fall back to approximate
    if final_score == i32::MIN {
        let matches = (n.min(m) as f64 * 0.9) as usize;
        let align_len = n.max(m);
        return (matches, align_len, String::new());
    }

    // Count matches by re-walking the band (simplified: use score-based estimate)
    // For a proper traceback we'd need O(n * band) backtrack array, which
    // is the same memory as the DP. Use ratio of actual score to perfect score.
    let perfect_score = n.min(m) as i32 * match_score;
    let score_ratio = if perfect_score > 0 {
        (final_score as f64 / perfect_score as f64).max(0.0).min(1.0)
    } else {
        0.0
    };
    let matches = (n.min(m) as f64 * score_ratio) as usize;
    let align_len = n.max(m);
    (matches, align_len, String::new())
}
