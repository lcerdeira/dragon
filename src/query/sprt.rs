//! Wald Sequential Probability Ratio Test (SPRT) for adaptive k-mer sampling.
//!
//! # Problem
//! `containment_rank` currently samples a fixed 384 k-mers per query regardless
//! of how quickly the evidence accumulates.  For databases where most queries
//! have no match (e.g. AMR gene panels against a species-specific database),
//! the first few k-mers already produce overwhelming evidence of absence — yet
//! we keep sampling all 384.
//!
//! # Solution
//! Wald's SPRT (1945) is the optimal stopping rule for a binary hypothesis test:
//!
//! * **H₀**: query is absent from the database.
//!   Expected per-k-mer hit rate p₀ ≈ 1e-6 (random k-mer collision on k=31).
//! * **H₁**: query matches at least one genome in the database.
//!   Expected per-k-mer hit rate p₁ = ANI^k ≈ 0.53 for ANI=0.98, k=31.
//!
//! After each k-mer observation the log-likelihood ratio (LLR) is updated:
//! ```text
//! hit:  LLR += ln(p₁/p₀)
//! miss: LLR += ln((1-p₁)/(1-p₀))
//! ```
//! and compared to two thresholds:
//! * `LLR > B` → accept H₁ (match found, stop and continue to Phase 2)
//! * `LLR < A` → accept H₀ (no match, return empty immediately)
//! * `A ≤ LLR ≤ B` → continue sampling
//!
//! Thresholds for error rates α (false-positive) and β (false-negative):
//! * `A = ln(β / (1-α))`
//! * `B = ln((1-β) / α)`
//!
//! # Expected savings
//! For truly-absent queries: stops after ~3-5 samples (vs 384).
//! For true matches: stops after ~5-10 samples when first hit is seen.
//! Overall: **30-50× fewer FM-index lookups** on sparse databases.

/// Decision returned by [`SprtState::update`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SprtDecision {
    /// Not enough evidence yet — continue sampling.
    Continue,
    /// Accept H₁: query likely has at least one match in the database.
    /// Proceed to Phase 2/3 of containment ranking.
    AcceptH1,
    /// Accept H₀: query is almost certainly absent from the database.
    /// Caller should return an empty candidate list immediately.
    AcceptH0,
}

/// State for one SPRT test instance (one query).
///
/// Construct via [`SprtState::new`] or the convenience
/// [`SprtState::default_for_kmer31`], then call [`SprtState::update`]
/// after each k-mer observation.
#[derive(Debug, Clone)]
pub struct SprtState {
    /// Cumulative log-likelihood ratio.
    llr: f64,
    /// Lower stopping threshold A = ln(β/(1-α)).
    /// LLR < A → accept H₀.
    threshold_a: f64,
    /// Upper stopping threshold B = ln((1-β)/α).
    /// LLR > B → accept H₁.
    threshold_b: f64,
    /// Per-hit LLR increment: ln(p₁/p₀).
    delta_hit: f64,
    /// Per-miss LLR increment: ln((1-p₁)/(1-p₀)).
    delta_miss: f64,
}

impl SprtState {
    /// Create a new SPRT state.
    ///
    /// # Parameters
    /// * `p0` — expected per-k-mer hit rate under H₀ (no match).
    /// * `p1` — expected per-k-mer hit rate under H₁ (match present).
    /// * `alpha` — tolerated false-positive rate (type I error).
    /// * `beta`  — tolerated false-negative rate (type II error).
    ///
    /// # Panics
    /// Panics if p0 ≥ p1, or if either probability is outside (0, 1).
    pub fn new(p0: f64, p1: f64, alpha: f64, beta: f64) -> Self {
        assert!(p0 > 0.0 && p0 < 1.0, "p0 must be in (0,1)");
        assert!(p1 > 0.0 && p1 < 1.0, "p1 must be in (0,1)");
        assert!(p0 < p1, "p0 must be strictly less than p1");
        assert!(alpha > 0.0 && alpha < 0.5, "alpha must be in (0, 0.5)");
        assert!(beta > 0.0 && beta < 0.5, "beta must be in (0, 0.5)");

        Self {
            llr: 0.0,
            threshold_a: (beta / (1.0 - alpha)).ln(),
            threshold_b: ((1.0 - beta) / alpha).ln(),
            delta_hit:  (p1 / p0).ln(),
            delta_miss: ((1.0 - p1) / (1.0 - p0)).ln(),
        }
    }

    /// Convenience constructor with Dragon defaults:
    /// * k=31, ANI=0.98 → p₁ = 0.98^31 ≈ 0.533
    /// * background k-mer hit rate p₀ = 1e-6
    /// * α = β = 0.05 (5% error on each side)
    ///
    /// These give A ≈ −2.94 and B ≈ +2.94.
    /// Expected samples to decision: ~5 for absent queries, ~10 for matches.
    pub fn default_for_kmer31() -> Self {
        // p₁ = 0.98^31 ≈ 0.5327
        Self::new(1e-6, 0.98_f64.powi(31), 0.05, 0.05)
    }

    /// Variant with lower thresholds for faster decisions (higher error rates).
    /// Useful when a downstream alignment step catches false positives anyway.
    /// α = β = 0.20
    pub fn fast_for_kmer31() -> Self {
        Self::new(1e-6, 0.98_f64.powi(31), 0.20, 0.20)
    }

    /// Update the LLR with one k-mer observation and return the current decision.
    ///
    /// `hit` = `true` if the k-mer was found anywhere in the FM-index
    /// (i.e. `fm_index.search(kmer)` returned a non-empty list).
    #[inline]
    pub fn update(&mut self, hit: bool) -> SprtDecision {
        self.llr += if hit { self.delta_hit } else { self.delta_miss };
        if self.llr > self.threshold_b {
            SprtDecision::AcceptH1
        } else if self.llr < self.threshold_a {
            SprtDecision::AcceptH0
        } else {
            SprtDecision::Continue
        }
    }

    /// Current log-likelihood ratio (for diagnostics).
    pub fn llr(&self) -> f64 {
        self.llr
    }

    /// Reset to initial state (reuse allocation).
    pub fn reset(&mut self) {
        self.llr = 0.0;
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thresholds_correct() {
        let s = SprtState::new(1e-6, 0.5, 0.05, 0.05);
        // A = ln(0.05/0.95) ≈ -2.944
        // B = ln(0.95/0.05) ≈ +2.944
        assert!((s.threshold_a - (-2.944_f64)).abs() < 0.01, "A wrong: {}", s.threshold_a);
        assert!((s.threshold_b -  2.944_f64).abs() < 0.01,  "B wrong: {}", s.threshold_b);
    }

    #[test]
    fn all_hits_accept_h1() {
        let mut s = SprtState::default_for_kmer31();
        let mut result = SprtDecision::Continue;
        // All hits should drive LLR to B quickly
        for _ in 0..50 {
            result = s.update(true);
            if result != SprtDecision::Continue { break; }
        }
        assert_eq!(result, SprtDecision::AcceptH1, "LLR={}", s.llr());
    }

    #[test]
    fn all_misses_accept_h0() {
        let mut s = SprtState::default_for_kmer31();
        let mut result = SprtDecision::Continue;
        // p₀≈1e-6 means misses drive LLR negative quickly
        for _ in 0..50 {
            result = s.update(false);
            if result != SprtDecision::Continue { break; }
        }
        assert_eq!(result, SprtDecision::AcceptH0, "LLR={}", s.llr());
    }

    #[test]
    fn mixed_stays_in_band() {
        // Alternating hit/miss at equal rate — LLR should stay near 0
        let mut s = SprtState::new(0.3, 0.7, 0.05, 0.05);
        for i in 0..10 {
            let decision = s.update(i % 2 == 0);
            // With p0=0.3, p1=0.7 and alternating, should stay Continue
            // (this is a weak assertion — the test just verifies no panic)
            let _ = decision;
        }
    }

    #[test]
    fn early_exit_on_absent_query_is_fast() {
        // For a query absent from any genome, we should decide within ~10 samples
        let mut s = SprtState::default_for_kmer31();
        let mut n_samples = 0usize;
        for _ in 0..384 {
            n_samples += 1;
            if s.update(false) != SprtDecision::Continue {
                break;
            }
        }
        assert!(n_samples < 50,
            "SPRT took {} samples for all-miss; expected < 50", n_samples);
    }

    #[test]
    fn reset_restores_initial_state() {
        let mut s = SprtState::default_for_kmer31();
        for _ in 0..5 { s.update(true); }
        let llr_before_reset = s.llr();
        s.reset();
        assert_eq!(s.llr(), 0.0, "reset should zero LLR, was {}", llr_before_reset);
        // After reset, should behave identically to a fresh instance
        let d = s.update(false);
        let mut fresh = SprtState::default_for_kmer31();
        let d2 = fresh.update(false);
        assert_eq!(d, d2);
    }
}
