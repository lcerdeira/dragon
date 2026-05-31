//! Bayesian posterior probability of genomic match.
//!
//! # Why Bayesian scoring?
//!
//! The standard containment score c = n_hits / N_sampled is a point estimate
//! that discards sample-size information.  Two queries with c = 0.95 are treated
//! identically whether they have 19/20 or 190/200 sampled k-mers — but the
//! confidence in those scores is very different.
//!
//! The Bayesian approach maintains a full posterior distribution over the true
//! containment p:
//!
//! ```text
//! Prior:     p ~ Beta(1, 1)          (uniform — no prior knowledge)
//! Likelihood: n | p ~ Binomial(N, p)
//! Posterior:  p | n, N ~ Beta(n+1, N-n+1)
//! ```
//!
//! From this posterior, we compute:
//!
//! **P(p ≥ θ | n, N)** — the probability that the TRUE containment exceeds
//! a threshold θ.  This is the *calibrated* probability that the genome is a
//! genuine match: if Dragon reports `bp:f:0.95`, approximately 95% of such
//! calls are correct (in expectation under the Beta-Binomial model).
//!
//! # Threshold choice
//!
//! The default threshold `θ = 0.5` asks: "Is this genome sharing more than
//! half its k-mers with the query?" — a natural, interpretable boundary
//! between a match and a random coincidence.
//!
//! For clinical applications, a higher threshold (θ = 0.9) can be used to
//! obtain high-confidence calls at the cost of reduced recall.
//!
//! # Numerical method
//!
//! For moderate N (≥30), the Beta posterior is well approximated by a Normal:
//!
//! ```text
//! Beta(n+1, N-n+1) ≈ Normal(μ, σ²)
//!   μ = (n+1) / (N+2)
//!   σ² = μ(1−μ) / (N+2)
//!
//! P(p ≥ θ) ≈ Φ((μ − θ) / σ)
//! ```
//!
//! where Φ is the standard-Normal CDF.  The Horner's method erf approximation
//! (Abramowitz & Stegun 7.1.26) is accurate to 1.5 × 10⁻⁷ for all x.
//!
//! # Output PAF tags
//!
//! ```text
//! bp:f:0.9991   Bayesian posterior P(containment ≥ 0.5 | data)
//! bq:f:0.9500   Bayesian P(containment ≥ 0.9 | data) — high-confidence
//! ```

/// Default threshold: probability that true containment ≥ 50%.
pub const DEFAULT_MATCH_THRESHOLD: f64 = 0.5;

/// High-confidence threshold: probability that true containment ≥ 90%.
/// Useful for clinical/epidemiological applications.
pub const HIGH_CONF_THRESHOLD: f64 = 0.90;

/// Compute the Bayesian posterior probability that the true k-mer containment
/// between a query and a reference genome exceeds `threshold`.
///
/// Uses the Beta(n+1, N-n+1) posterior with a uniform Beta(1,1) prior and
/// the Normal approximation (valid for N ≥ 30).
///
/// # Parameters
/// - `n_hits`: number of sampled k-mers that matched the genome
/// - `n_total`: total number of k-mers sampled from the query
/// - `threshold`: minimum true containment to be considered a match (default: 0.5)
///
/// # Returns
/// `None` if `n_total == 0` (no data); otherwise `Some(prob)` in [0, 1].
///
/// # Example
/// ```
/// use dragon::query::bayes::bayesian_match_prob;
/// // 192/200 k-mers match at θ=0.5 → near-certain match
/// let p = bayesian_match_prob(192, 200, 0.5).unwrap();
/// assert!(p > 0.999, "high containment should give p≈1.0, got {}", p);
///
/// // 1/20 k-mers match → near-certain non-match
/// let p = bayesian_match_prob(1, 20, 0.5).unwrap();
/// assert!(p < 0.01, "low containment should give p≈0.0, got {}", p);
/// ```
pub fn bayesian_match_prob(n_hits: usize, n_total: usize, threshold: f64) -> Option<f64> {
    if n_total == 0 {
        return None;
    }

    // Posterior: Beta(n+1, N-n+1) — add pseudocounts for Laplace smoothing
    let n = n_hits as f64 + 1.0;
    let total = n_total as f64 + 2.0;   // N + alpha + beta = N + 1 + 1

    let mu = n / total;                  // posterior mean
    let variance = mu * (1.0 - mu) / total;

    // Edge cases: posterior is degenerate
    if variance <= 0.0 {
        return Some(if mu >= threshold { 1.0 } else { 0.0 });
    }

    // P(p ≥ threshold) = Φ((μ − θ) / σ)  ← CDF, NOT survival function.
    // When μ > θ: z > 0, Φ(z) > 0.5 ← genome is likely a match ✓
    // When μ < θ: z < 0, Φ(z) < 0.5 ← genome is likely not a match ✓
    let z = (mu - threshold) / variance.sqrt();
    Some(standard_normal_cdf(z))
}

/// Compute both standard (θ=0.5) and high-confidence (θ=0.9) Bayesian probs.
///
/// Returns `(p_standard, p_highconf)`.
pub fn bayesian_probs(n_hits: usize, n_total: usize) -> (Option<f64>, Option<f64>) {
    (
        bayesian_match_prob(n_hits, n_total, DEFAULT_MATCH_THRESHOLD),
        bayesian_match_prob(n_hits, n_total, HIGH_CONF_THRESHOLD),
    )
}

/// Bayesian ANI point estimate: the posterior mean of (containment)^(1/k).
///
/// Unlike the frequentist ANI = (n/N)^(1/k), this uses the Bayes shrinkage
/// estimator which pulls the estimate toward the centre when evidence is weak:
///
/// ```text
/// ANI_Bayes = ((n+1)/(N+2))^(1/k)
/// ```
///
/// For large N (≥384), the shrinkage is negligible (<0.001 ANI units).
/// For small N (short reads, few k-mers), shrinkage prevents extreme estimates.
pub fn bayesian_ani(n_hits: usize, n_total: usize, kmer_size: usize) -> Option<f64> {
    if n_total == 0 || kmer_size == 0 {
        return None;
    }
    let posterior_mean = (n_hits as f64 + 1.0) / (n_total as f64 + 2.0);
    if posterior_mean <= 0.0 {
        return Some(0.0);
    }
    Some(posterior_mean.powf(1.0 / kmer_size as f64))
}

// ─── Internal numerical routines ──────────────────────────────────────────────

/// Standard Normal CDF Φ(z) = P(Z ≤ z).
///
/// Uses the error function: Φ(z) = (1 + erf(z/√2)) / 2.
/// Accurate to ~1.5e-7 for all z via the A&S 7.1.26 erf approximation.
#[inline]
fn standard_normal_cdf(z: f64) -> f64 {
    let result = 0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2));
    result.clamp(0.0, 1.0)
}

/// Error function erf(x).
///
/// Horner's method approximation (Abramowitz & Stegun 7.1.26).
/// Maximum error: |ε| < 1.5 × 10⁻⁷ for all real x.
fn erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t * (0.254829592
        + t * (-0.284496736
            + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let correction = poly * (-(x * x)).exp();
    // erf(x) = sign(x) × (1 − correction)
    let result = 1.0 - correction;
    if x >= 0.0 { result } else { -result }
}

/// Complementary error function erfc(x) = 1 − erf(x).
#[allow(dead_code)]
fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_containment_gives_high_prob() {
        let p = bayesian_match_prob(384, 384, 0.5).unwrap();
        assert!(p > 0.9999, "100% containment → p≈1; got {:.6}", p);
    }

    #[test]
    fn zero_hits_gives_low_prob() {
        let p = bayesian_match_prob(0, 384, 0.5).unwrap();
        assert!(p < 1e-10, "0 hits → p≈0; got {:e}", p);
    }

    #[test]
    fn sample_size_matters() {
        // Use borderline containment (60%) close to threshold (0.5) to show
        // sample-size effect. With large N the posterior is tight; small N gives wider.
        let p_large = bayesian_match_prob(120, 200, 0.5).unwrap(); // 60%, N=200
        let p_small = bayesian_match_prob(12, 20, 0.5).unwrap();   // 60%, N=20
        assert!(p_large > p_small,
            "larger N should give higher confidence: {:.4} vs {:.4}", p_large, p_small);
        // Large N: tight posterior around 0.60, θ=0.50 → high confidence
        assert!(p_large > 0.99, "120/200 at θ=0.5 should be high-confidence: {:.4}", p_large);
        // Small N: wide posterior → more uncertainty
        assert!(p_small < 0.98, "12/20 at θ=0.5 should be less certain: {:.4}", p_small);
    }

    #[test]
    fn high_conf_threshold_is_stricter() {
        let (p_std, p_hc) = bayesian_probs(192, 384);  // ~50% containment
        let p_std = p_std.unwrap();
        let p_hc = p_hc.unwrap();
        assert!(p_std > p_hc,
            "standard threshold (0.5) should be easier to satisfy than high-conf (0.9)");
        // At 50% containment, θ=0.9 should be essentially impossible
        assert!(p_hc < 0.01,
            "50% containment can't meet θ=0.9 threshold: p_hc={:.4}", p_hc);
    }

    #[test]
    fn calibration_midpoint() {
        // At exactly the threshold, prob should be ~0.5
        // 50% containment, threshold=0.5 → prob ≈ 0.5
        let p = bayesian_match_prob(192, 384, 0.5).unwrap();
        // With N=384, posterior is tight: should be very close to 0.5
        assert!((p - 0.5).abs() < 0.05,
            "at threshold midpoint, prob should be ~0.5: got {:.4}", p);
    }

    #[test]
    fn no_data_returns_none() {
        assert!(bayesian_match_prob(0, 0, 0.5).is_none());
    }

    #[test]
    fn bayesian_ani_shrunken() {
        // All hits: frequentist ANI = 1.0, Bayesian ANI slightly < 1.0 (Laplace shrinkage)
        let ani = bayesian_ani(384, 384, 31).unwrap();
        assert!(ani < 1.0, "Bayesian ANI should be shrunk from 1.0 by pseudocount: {}", ani);
        assert!(ani > 0.999, "Shrinkage with N=384 is tiny (≤0.001 ANI): {}", ani);

        // 0 hits: Laplace pseudocount gives posterior mean 1/(N+2) — small but non-zero.
        // This is mathematically correct (don't assign zero probability to presence).
        // Use bp:f: tag for match probability, not ba: for threshold testing.
        let ani_zero = bayesian_ani(0, 384, 31).unwrap();
        assert!(ani_zero > 0.0, "Laplace pseudocount: ANI > 0 even with 0 hits: {}", ani_zero);
        // ba tag is meaningful alongside bp:f:~0 which correctly signals non-match.
        // The Bayesian ANI here ≈ (1/386)^(1/31) ≈ 0.825 — a shrinkage artefact.
        assert!(ani_zero < 0.9, "0 hits should give ani << 1.0: {}", ani_zero);
    }

    #[test]
    fn short_read_uncertainty() {
        // 15-mer sampled 50 times: fewer samples → more uncertainty
        let (p50, p10) = (
            bayesian_match_prob(45, 50, 0.5).unwrap(),   // 90% containment, N=50
            bayesian_match_prob(45, 50, 0.85).unwrap(),  // same, stricter threshold
        );
        assert!(p50 > 0.99, "90% containment at θ=0.5 should be very confident: {}", p50);
        assert!(p10 > 0.5,  "90% containment at θ=0.85 with N=50 should be moderately confident: {}", p10);
    }

    #[test]
    fn erfc_boundary_values() {
        // erfc(0) = 1.0, erfc(∞) → 0, erfc(-∞) → 2
        let e0 = erfc(0.0);
        assert!((e0 - 1.0).abs() < 1e-6, "erfc(0)=1: {}", e0);
        let e_large = erfc(6.0);
        assert!(e_large < 1e-7, "erfc(6)≈0: {}", e_large);
        let e_neg = erfc(-6.0);
        assert!((e_neg - 2.0).abs() < 1e-7, "erfc(-6)≈2: {}", e_neg);
    }
}
