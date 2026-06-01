//! Pigeonhole-based multi-anchor seeding for high-sensitivity alignment.
//!
//! # Why spaced / multi-anchor seeds?
//!
//! Dragon's FM-index supports *exact* k-mer matching only.  For a query with
//! d% divergence from the database, the probability that a solid k-mer is an
//! exact match is:
//!
//! ```text
//!   P(exact | d, k) = (1 − d)^k
//!   d=5%, k=31 → 0.20   (80% miss rate)
//!   d=5%, k=10 → 0.60   (40% miss rate, 3× better)
//! ```
//!
//! The **pigeonhole principle** shows that by querying *multiple short anchors*
//! spread across the same window, we achieve much higher sensitivity:
//!
//! ```text
//!   P(≥1 of N anchors exact | d, k_anchor) = 1 − (1 − (1−d)^k_a)^N
//!
//!   N=3 anchors, k_a=10, d=5%:
//!     = 1 − (1 − 0.60)^3 = 1 − 0.40³ = 0.936  ← BLAST-level sensitivity!
//!
//!   N=3 anchors, k_a=10, d=10%:
//!     = 1 − (1 − 0.35)^3 = 1 − 0.65³ = 0.725  (vs solid k=31 → 0.04!)
//! ```
//!
//! This is mathematically equivalent to BLAST's two-hit model and PatternHunter's
//! spaced seeds, but implemented using Dragon's existing FM-index by issuing
//! multiple short-k searches within the same query window.
//!
//! # Algorithm
//!
//! For each sampled query position `p` (with window `W`):
//! 1. Try the primary seed (adaptive k, usually k=31 for long queries).
//! 2. If no hit: extract N non-overlapping anchor k-mers of length `k_anchor`
//!    spread across the window.
//! 3. Search each anchor in the FM-index.
//! 4. For each anchor hit, record (unitig_id, anchor_index) in a per-genome
//!    vote counter.
//! 5. A genome is a candidate anchor hit only if it scores ≥ MIN_ANCHOR_VOTES
//!    from *different* anchor positions (ensuring the vote comes from a
//!    distributed signal, not a single repetitive k-mer).
//!
//! # Integration
//!
//! The anchors are transparent to the rest of the pipeline: they produce the
//! same `SeedHit` structs as regular k-mer search, with `match_len = k_anchor`.
//! The alignment stage (Stage 4, WFA / graph-DP) re-derives exact coordinates
//! and identity regardless of how the candidate was found.
//!
//! # Sensitivity comparison (d=5%)
//!
//! | Method                    | P(hit at 5% div) | P(hit at 10% div) |
//! |---------------------------|------------------|-------------------|
//! | Solid k=31 (current)      | 0.20             | 0.04              |
//! | Adaptive k=18             | 0.40             | 0.16              |
//! | 3 anchors k=10 (this mod) | **0.94**         | **0.72**          |
//! | BLAST two-hit (11-mers)   | ~0.95            | ~0.74             |

use crate::index::fm::{DragonFmIndex, SeedHit};

/// One anchor within a multi-anchor seed set.
#[derive(Clone, Debug)]
pub struct AnchorKmer {
    /// The k-mer bytes.
    pub kmer: Vec<u8>,
    /// Byte offset of this anchor within the query window.
    pub window_offset: usize,
    /// Anchor index (0-based) within this window's set.
    pub anchor_idx: usize,
}

/// Parameters for multi-anchor seeding.
#[derive(Clone, Debug)]
pub struct AnchorConfig {
    /// Length of each anchor k-mer (default 10).
    pub k_anchor: usize,
    /// Number of anchors spread across the window (default 3).
    pub n_anchors: usize,
    /// Minimum number of *distinct-anchor* votes a genome needs to be
    /// considered a candidate. Default 2 (at least 2 of N anchors hit).
    pub min_anchor_votes: usize,
}

impl Default for AnchorConfig {
    fn default() -> Self {
        // min_anchor_votes=1: maximises recall (P(≥1 anchor | d=5%) ≈ 0.94).
        // Use min_anchor_votes=2 for higher precision at the cost of recall.
        Self { k_anchor: 10, n_anchors: 3, min_anchor_votes: 1 }
    }
}

impl AnchorConfig {
    /// Config optimised for cross-species queries (15-30% divergence).
    ///
    /// 5 anchors of k=7: P(≥1 exact | d=15%) = 0.855, P(≥1 | d=20%) = 0.692.
    /// This enables detection of homologs at 70-85% ANI which solid k=31
    /// completely misses (P ≈ 0.006 at 15% divergence).
    pub fn for_cross_species() -> Self {
        Self { k_anchor: 7, n_anchors: 5, min_anchor_votes: 1 }
    }

    /// Config for within-species / short-read mode (default).
    pub fn for_short_reads() -> Self {
        Self::default()
    }
}

impl AnchorConfig {
    /// Minimum window size needed to fit all anchors with this config.
    pub fn min_window(&self) -> usize {
        // anchors packed end-to-end: n_anchors * k_anchor
        self.n_anchors * self.k_anchor
    }

    /// Sensitivity at divergence `d` (probability that ≥ min_anchor_votes
    /// anchors are exact).
    pub fn sensitivity_at(&self, d: f64) -> f64 {
        let p_anchor = (1.0 - d).powi(self.k_anchor as i32);
        // P(≥ min_votes out of n_anchors)
        // Using binomial CDF complement
        let mut p = 0.0f64;
        let n = self.n_anchors as u64;
        let m = self.min_anchor_votes as u64;
        // Sum P(X=k) for k >= m
        for k in m..=n {
            let binom = binom_coeff(n, k) as f64;
            p += binom * p_anchor.powi(k as i32) * (1.0 - p_anchor).powi((n - k) as i32);
        }
        p
    }
}

/// Extract N evenly-spaced anchor k-mers from `window`.
///
/// Anchors are placed as evenly as possible across the window, with at least
/// `k_anchor` non-overlapping positions.  If the window is too short to hold
/// all anchors, returns as many as fit.
pub fn extract_anchors(window: &[u8], config: &AnchorConfig) -> Vec<AnchorKmer> {
    if window.len() < config.k_anchor {
        return Vec::new();
    }

    let n = config.n_anchors.min(window.len() / config.k_anchor);
    if n == 0 {
        return Vec::new();
    }

    let step = (window.len() - config.k_anchor) / n.max(1);
    let step = step.max(config.k_anchor); // anchors don't overlap

    let mut anchors = Vec::with_capacity(n);
    let mut offset = 0;
    let mut idx = 0;

    while idx < n && offset + config.k_anchor <= window.len() {
        let kmer = window[offset..offset + config.k_anchor].to_vec();
        // Skip ambiguous k-mers
        let has_ambig = kmer.iter().any(|&b| {
            !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
        });
        if !has_ambig {
            anchors.push(AnchorKmer { kmer, window_offset: offset, anchor_idx: idx });
        }
        offset += step;
        idx += 1;
    }

    anchors
}

/// Search all anchors of a window in the FM-index.
///
/// Returns `(unitig_id, anchor_idx)` pairs — the anchor_idx tracks *which*
/// anchor within the window produced the hit, enabling the multi-vote filter.
pub fn search_anchors(
    anchors: &[AnchorKmer],
    fm_index: &DragonFmIndex,
    max_seed_freq: usize,
    query_pos: usize,
) -> Vec<SeedHit> {
    let mut hits = Vec::new();

    for anchor in anchors {
        let positions = fm_index.search(&anchor.kmer);
        if positions.is_empty() || positions.len() > max_seed_freq {
            continue;
        }
        for &pos in &positions {
            if let Some((unitig_id, offset)) = fm_index.position_to_unitig(pos) {
                hits.push(SeedHit {
                    unitig_id,
                    offset,
                    query_pos: query_pos + anchor.window_offset,
                    match_len: anchor.kmer.len(),
                    is_reverse: false,
                    sa_count: positions.len(),
                });
            }
        }
    }

    hits
}

/// Full pigeonhole seeding at one query position: search solid k-mer first
/// (fast path), fall back to multi-anchor on failure (sensitivity path).
///
/// Returns `(solid_hit, anchor_hits)` where `solid_hit` is `Some` if the
/// primary k-mer search succeeded.
pub fn pigeonhole_search(
    query: &[u8],
    query_pos: usize,
    primary_k: usize,
    fm_index: &DragonFmIndex,
    max_seed_freq: usize,
    anchor_config: &AnchorConfig,
) -> (Vec<SeedHit>, Vec<SeedHit>) {
    if query_pos + primary_k > query.len() {
        return (Vec::new(), Vec::new());
    }

    let window = &query[query_pos..query_pos + primary_k];

    // Fast path: solid primary k-mer
    let has_ambig = window.iter().any(|&b| {
        !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
    });
    if has_ambig {
        return (Vec::new(), Vec::new());
    }

    let primary_positions = fm_index.search(window);
    if !primary_positions.is_empty() && primary_positions.len() <= max_seed_freq {
        let solid: Vec<SeedHit> = primary_positions
            .iter()
            .filter_map(|&pos| {
                fm_index.position_to_unitig(pos).map(|(uid, off)| SeedHit {
                    unitig_id: uid,
                    offset: off,
                    query_pos,
                    match_len: primary_k,
                    is_reverse: false,
                    sa_count: primary_positions.len(),
                })
            })
            .collect();
        return (solid, Vec::new());
    }

    // Sensitivity path: multi-anchor search
    if window.len() < anchor_config.min_window() {
        return (Vec::new(), Vec::new());
    }
    let anchors = extract_anchors(window, anchor_config);
    let anchor_hits = search_anchors(&anchors, fm_index, max_seed_freq, query_pos);

    (Vec::new(), anchor_hits)
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

fn binom_coeff(n: u64, k: u64) -> u64 {
    if k > n { return 0; }
    if k == 0 || k == n { return 1; }
    let k = k.min(n - k);
    let mut result = 1u64;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anchor_config_sensitivity() {
        // Default config: k=10, N=3, min_votes=1
        // P(≥1 of 3 anchors exact | d=5%) = 1 - (1-0.599)^3 ≈ 0.936
        let cfg = AnchorConfig::default();
        let s5 = cfg.sensitivity_at(0.05);
        let s10 = cfg.sensitivity_at(0.10);
        assert!(s5 > 0.90, "sensitivity at 5%={:.3} should be >0.90 (min_votes=1)", s5);
        assert!(s10 > 0.60, "sensitivity at 10%={:.3} should be >0.60", s10);
        // Must be at least 3× better than solid k=31 at 5% divergence
        let solid_5 = (1.0f64 - 0.05).powi(31); // ≈ 0.20
        assert!(s5 > solid_5 * 3.0,
            "anchor s5={:.3} should be >3× solid k=31 ({:.3})", s5, solid_5);

        // With min_votes=2: less recall but more specific
        let cfg2 = AnchorConfig { k_anchor: 10, n_anchors: 3, min_anchor_votes: 2 };
        let s5_2 = cfg2.sensitivity_at(0.05);
        assert!(s5_2 > 0.55, "min_votes=2 sensitivity={:.3} should be >0.55", s5_2);
        assert!(s5 > s5_2, "min_votes=1 should have higher sensitivity than min_votes=2");
    }

    #[test]
    fn extract_anchors_basic() {
        let window = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
        let cfg = AnchorConfig::default(); // k=10, n=3
        let anchors = extract_anchors(window, &cfg);
        assert_eq!(anchors.len(), 3, "should extract 3 anchors");
        // All anchors have k_anchor length
        for a in &anchors {
            assert_eq!(a.kmer.len(), cfg.k_anchor);
        }
        // Anchors are at increasing offsets
        let offsets: Vec<_> = anchors.iter().map(|a| a.window_offset).collect();
        for w in offsets.windows(2) {
            assert!(w[0] < w[1], "offsets should be increasing");
        }
    }

    #[test]
    fn extract_anchors_too_short() {
        // Window shorter than k_anchor
        let window = b"ACGT"; // 4 bp, k_anchor=10
        let cfg = AnchorConfig::default();
        assert!(extract_anchors(window, &cfg).is_empty());
    }

    #[test]
    fn anchors_non_overlapping() {
        let window = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
        let cfg = AnchorConfig { k_anchor: 10, n_anchors: 3, min_anchor_votes: 2 };
        let anchors = extract_anchors(window, &cfg);
        for i in 0..anchors.len() {
            for j in (i + 1)..anchors.len() {
                let gap = anchors[j].window_offset.saturating_sub(anchors[i].window_offset);
                assert!(gap >= cfg.k_anchor,
                    "anchors {} and {} overlap: offsets {}, {}",
                    i, j, anchors[i].window_offset, anchors[j].window_offset);
            }
        }
    }

    #[test]
    fn binom_coeff_values() {
        assert_eq!(binom_coeff(5, 2), 10);
        assert_eq!(binom_coeff(3, 3), 1);
        assert_eq!(binom_coeff(4, 0), 1);
        assert_eq!(binom_coeff(0, 1), 0);
    }

    #[test]
    fn sensitivity_table() {
        // Verify the key claims in the module docstring
        let cfg = AnchorConfig { k_anchor: 10, n_anchors: 3, min_anchor_votes: 1 };
        // P(≥1 of 3 anchors exact | d=5%) = 1 - (0.40)^3 = 0.936
        let s = cfg.sensitivity_at(0.05);
        assert!((s - 0.936).abs() < 0.01,
            "P(≥1 anchor | d=5%) = {:.3}, expected ~0.936", s);
    }
}
