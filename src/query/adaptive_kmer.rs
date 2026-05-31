//! Adaptive k-mer size selection for containment ranking.
//!
//! # The core problem
//!
//! Dragon's FM-index is built at k=31.  For containment ranking, we sample
//! the query at stride and search 31-mers in the index.  The fraction of
//! query k-mers that survive as exact matches under divergence d is:
//!
//! ```text
//!   P(exact match) = (1 − d)^k
//!
//!   d=5%,  k=31 → 0.20   (only 20% of k-mers survive → low recall)
//!   d=5%,  k=21 → 0.34
//!   d=5%,  k=15 → 0.46
//!   d=10%, k=31 → 0.04   (catastrophic)
//!   d=10%, k=15 → 0.21
//! ```
//!
//! A fixed k=31 is optimal for specificity (minimising false-positive genome
//! candidates) but catastrophic for sensitivity when divergence exceeds ~3%
//! or when queries are too short for enough seeds to accumulate.
//!
//! # Solution: adaptive k
//!
//! Choose the *largest* k such that the expected number of exact-match seeds
//! from the query exceeds a minimum threshold `MIN_EXPECTED_SEEDS`.  This
//! ensures we always have enough seed evidence to reliably identify the source
//! genome, while still using the highest specificity k possible.
//!
//! Expected seeds = (query_len − k + 1) × (1 − d_assumed)^k
//!
//! We do not know the true divergence; we assume `D_ASSUMED = 0.05` (5%) as a
//! conservative upper bound for within-species surveillance.  If the query is
//! actually at 0% divergence, the larger k still works (we just use k=31 as
//! normal).  If the query is more divergent, the adaptive k gracefully
//! degrades further.
//!
//! # Multi-k fallback
//!
//! For maximum recall, `multi_k_containment` runs containment at the adaptive
//! k first.  If it finds fewer than `MIN_CANDIDATES` candidate genomes, it
//! falls back to progressively smaller k values until either enough candidates
//! are found or the minimum k (15) is reached.
//!
//! # The FM-index is k-agnostic
//!
//! The FM-index built at k=31 can search *any* pattern length ≤ 31.  Searching
//! with k=15 in a k=31 index is fully valid: we find all exact 15-mer matches
//! in the indexed text.  The colour bitmaps are at unitig level (independent
//! of k), so the containment score remains meaningful.  The tradeoff: shorter
//! k → more random matches → lower specificity → more candidate genomes →
//! higher recall at the cost of more alignment work in Stage 4.

/// Assumed divergence for k-mer size selection (conservative 5%).
const D_ASSUMED: f64 = 0.05;

/// Minimum expected exact-match seeds before we accept the chosen k.
/// Below this, containment estimates are unreliable.
const MIN_EXPECTED_SEEDS: f64 = 50.0;

/// Absolute minimum k regardless of query length.
const K_MIN: usize = 15;

/// Minimum candidate genomes before multi-k fallback is triggered.
pub const MIN_CANDIDATES_FOR_EARLY_EXIT: usize = 1;

/// Select the largest k ≤ `max_k` such that:
///   E[exact seeds] = (query_len − k + 1) × (1 − D_ASSUMED)^k ≥ MIN_EXPECTED_SEEDS
///
/// Returns a value in [K_MIN, max_k].
///
/// # Examples
/// ```
/// use dragon::query::adaptive_kmer::adaptive_k;
/// // Short read, 150 bp: k=31 gives only ~24 expected seeds; adaptive k drops to ~18
/// let k = adaptive_k(150, 31);
/// assert!(k < 31, "should reduce k for short reads at 5% divergence");
/// assert!(k >= 15, "never below minimum");
///
/// // Long contig, 5000 bp: k=31 gives ~940 expected seeds; keep k=31
/// let k = adaptive_k(5000, 31);
/// assert_eq!(k, 31);
/// ```
pub fn adaptive_k(query_len: usize, max_k: usize) -> usize {
    // Walk down from max_k to K_MIN, stopping at the first k that satisfies
    // the expected-seed threshold.
    for k in (K_MIN..=max_k).rev() {
        let n_positions = query_len.saturating_sub(k) + 1;
        if n_positions == 0 {
            continue;
        }
        let p_exact = (1.0 - D_ASSUMED).powi(k as i32);
        let expected = n_positions as f64 * p_exact;
        if expected >= MIN_EXPECTED_SEEDS {
            return k;
        }
    }
    K_MIN
}

/// Return the sequence of k values to try in multi-k fallback order.
///
/// Starts at `primary_k` (adaptive), then falls back to larger steps if
/// needed: [primary_k, primary_k-4, primary_k-8, ..., K_MIN].
/// Each step reduces k by 4, giving a geometric improvement in sensitivity.
pub fn fallback_k_sequence(primary_k: usize) -> Vec<usize> {
    let mut ks = vec![primary_k];
    let mut k = primary_k;
    while k > K_MIN + 4 {
        k = k.saturating_sub(4).max(K_MIN);
        if !ks.contains(&k) {
            ks.push(k);
        }
    }
    if !ks.contains(&K_MIN) {
        ks.push(K_MIN);
    }
    ks
}

/// Theoretical sensitivity gain from reducing k at assumed divergence D_ASSUMED.
///
/// Returns the ratio of expected exact-match seeds at `new_k` vs `old_k`.
/// Useful for logging and for deciding whether a fallback is worth the cost.
pub fn sensitivity_gain(old_k: usize, new_k: usize) -> f64 {
    if old_k == new_k || new_k == 0 {
        return 1.0;
    }
    let p_old = (1.0 - D_ASSUMED).powi(old_k as i32);
    let p_new = (1.0 - D_ASSUMED).powi(new_k as i32);
    if p_old <= 0.0 {
        return f64::INFINITY;
    }
    p_new / p_old
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_contig_keeps_max_k() {
        // 5000 bp: expected = 4970 × 0.20 = 994 seeds at k=31 → keep k=31
        assert_eq!(adaptive_k(5000, 31), 31);
    }

    #[test]
    fn short_read_reduces_k() {
        // 150 bp: at k=31, expected = 120 × 0.20 = 24 < 50 → reduce
        let k = adaptive_k(150, 31);
        assert!(k < 31, "short read should reduce k from 31, got {}", k);
        assert!(k >= K_MIN, "never below minimum {}", K_MIN);
    }

    #[test]
    fn medium_read_boundary() {
        // 500 bp: at k=31, expected = 470 × 0.20 = 94 > 50 → keep k=31
        assert_eq!(adaptive_k(500, 31), 31);
    }

    #[test]
    fn very_short_read_hits_minimum() {
        // 50 bp: even k=15 gives only 36 expected seeds — still below 50
        // but we cap at K_MIN
        let k = adaptive_k(50, 31);
        assert_eq!(k, K_MIN);
    }

    #[test]
    fn fallback_sequence_decreasing() {
        let seq = fallback_k_sequence(31);
        assert_eq!(seq[0], 31, "first element is primary_k");
        assert_eq!(*seq.last().unwrap(), K_MIN, "last element is K_MIN");
        for w in seq.windows(2) {
            assert!(w[0] > w[1], "fallback sequence must be strictly decreasing");
        }
    }

    #[test]
    fn fallback_sequence_no_duplicates() {
        let seq = fallback_k_sequence(21);
        let unique: std::collections::HashSet<_> = seq.iter().collect();
        assert_eq!(unique.len(), seq.len(), "no duplicates");
    }

    #[test]
    fn sensitivity_gain_positive() {
        // Going from k=31 to k=15 should increase expected seeds
        let gain = sensitivity_gain(31, 15);
        assert!(gain > 1.0, "smaller k should improve sensitivity, gain={}", gain);
    }

    #[test]
    fn adaptive_k_satisfies_threshold_when_possible() {
        // For 300 bp at 5% divergence, the chosen k should give ≥ 50 expected seeds
        let k = adaptive_k(300, 31);
        let n = 300usize.saturating_sub(k) + 1;
        let expected = n as f64 * (1.0 - D_ASSUMED).powi(k as i32);
        // Either expected >= threshold, OR we're at K_MIN (can't go lower)
        assert!(
            expected >= MIN_EXPECTED_SEEDS || k == K_MIN,
            "adaptive_k={} gives {:.1} expected seeds (threshold {})",
            k, expected, MIN_EXPECTED_SEEDS
        );
    }
}
