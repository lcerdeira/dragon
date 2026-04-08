/// Learned seed scoring via logistic regression.
///
/// Replaces the raw `match_len` anchor weight with a quality score learned
/// from labeled training data. Features include seed frequency, match length,
/// positional bias, and genome specificity. Model is a logistic regression —
/// inference is a single dot product + sigmoid with zero ML dependencies.

use crate::index::color::ColorIndex;
use crate::index::fm::SeedHit;

/// Per-seed feature vector for ML scoring.
#[derive(Clone, Debug)]
pub struct SeedFeatures {
    pub match_len: f64,
    pub log_sa_count: f64,
    pub query_pos_frac: f64,
    pub match_frac: f64,
    pub color_cardinality: f64,
    pub gc_content: f64,
}

impl SeedFeatures {
    /// Compute features for a seed hit.
    pub fn from_seed(
        seed: &SeedHit,
        query_len: usize,
        query_bytes: &[u8],
        color_index: &ColorIndex,
    ) -> Self {
        let match_len = seed.match_len as f64;
        let log_sa_count = if seed.sa_count > 0 {
            (seed.sa_count as f64).log2()
        } else {
            0.0
        };
        let query_pos_frac = if query_len > 0 {
            seed.query_pos as f64 / query_len as f64
        } else {
            0.0
        };
        let match_frac = if query_len > 0 {
            seed.match_len as f64 / query_len as f64
        } else {
            0.0
        };

        // Color cardinality: how many genomes contain this unitig
        let color_cardinality = color_index
            .get_colors(seed.unitig_id)
            .map(|c| c.len() as f64)
            .unwrap_or(0.0);

        // GC content of the seed region in the query
        let gc_content = if seed.match_len > 0 && seed.query_pos + seed.match_len <= query_bytes.len() {
            let region = &query_bytes[seed.query_pos..seed.query_pos + seed.match_len];
            let gc = region.iter().filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c')).count();
            gc as f64 / seed.match_len as f64
        } else {
            0.5
        };

        Self {
            match_len,
            log_sa_count,
            query_pos_frac,
            match_frac,
            color_cardinality,
            gc_content,
        }
    }

    /// Convert to a fixed-size array for dot product.
    #[inline]
    fn as_array(&self) -> [f64; 6] {
        [
            self.match_len,
            self.log_sa_count,
            self.query_pos_frac,
            self.match_frac,
            self.color_cardinality,
            self.gc_content,
        ]
    }
}

/// Logistic regression seed scorer.
///
/// Weights: [bias, match_len, log_sa_count, query_pos_frac, match_frac, color_cardinality, gc_content]
pub struct SeedScorer {
    pub weights: [f64; 7],
}

/// Hand-tuned default weights. These encode domain knowledge:
/// - Longer seeds are more informative (positive weight)
/// - High-frequency seeds are likely repetitive (negative weight)
/// - Seeds in many genomes are less specific (negative weight)
/// - Extreme GC content correlates with repetitive regions (slight negative)
const DEFAULT_WEIGHTS: [f64; 7] = [
    0.5,    // bias
    0.15,   // match_len: longer = better
    -0.3,   // log_sa_count: frequent = likely repetitive
    0.0,    // query_pos_frac: neutral (no positional bias)
    0.1,    // match_frac: covers more query = better
    -0.2,   // color_cardinality: in many genomes = less specific
    -0.05,  // gc_content: extreme GC = slightly worse
];

impl Default for SeedScorer {
    fn default() -> Self {
        Self {
            weights: DEFAULT_WEIGHTS,
        }
    }
}

impl SeedScorer {
    /// Load weights from a JSON file, or use defaults if not found.
    pub fn load_or_default(index_dir: &std::path::Path) -> Self {
        let weights_path = index_dir.join("seed_scorer.json");
        if weights_path.exists() {
            match std::fs::read_to_string(&weights_path) {
                Ok(contents) => {
                    match serde_json::from_str::<Vec<f64>>(&contents) {
                        Ok(w) if w.len() == 7 => {
                            log::info!("Loaded ML seed scorer weights from {:?}", weights_path);
                            let mut weights = [0.0f64; 7];
                            weights.copy_from_slice(&w);
                            return Self { weights };
                        }
                        Ok(w) => {
                            log::warn!(
                                "seed_scorer.json has {} weights (expected 7), using defaults",
                                w.len()
                            );
                        }
                        Err(e) => {
                            log::warn!("Failed to parse seed_scorer.json: {}, using defaults", e);
                        }
                    }
                }
                Err(e) => {
                    log::warn!("Failed to read seed_scorer.json: {}, using defaults", e);
                }
            }
        }
        Self::default()
    }

    /// Score a seed using logistic regression.
    ///
    /// Returns a value in (0, 1) where higher means the seed is more likely
    /// to be a true positive. The score incorporates match length, so it
    /// can directly replace `match_len` in the chaining DP.
    #[inline]
    pub fn score(&self, features: &SeedFeatures) -> f64 {
        let feat = features.as_array();
        let mut z = self.weights[0]; // bias
        for i in 0..6 {
            z += self.weights[i + 1] * feat[i];
        }
        sigmoid(z)
    }

    /// Score a seed, scaled by match_len for use as chain anchor weight.
    ///
    /// This preserves the property that longer seeds contribute more to
    /// chain score, while down-weighting low-quality seeds.
    #[inline]
    pub fn weighted_score(&self, features: &SeedFeatures) -> f64 {
        let quality = self.score(features);
        features.match_len * quality
    }
}

#[inline]
fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sigmoid() {
        assert!((sigmoid(0.0) - 0.5).abs() < 1e-10);
        assert!(sigmoid(10.0) > 0.99);
        assert!(sigmoid(-10.0) < 0.01);
    }

    #[test]
    fn test_default_scorer() {
        let scorer = SeedScorer::default();

        // A long, rare seed should score high
        let good_seed = SeedFeatures {
            match_len: 50.0,
            log_sa_count: 1.0,  // only 2 occurrences
            query_pos_frac: 0.5,
            match_frac: 0.05,
            color_cardinality: 1.0,
            gc_content: 0.5,
        };

        // A short, frequent seed should score lower
        let bad_seed = SeedFeatures {
            match_len: 15.0,
            log_sa_count: 13.0,  // ~8000 occurrences
            query_pos_frac: 0.5,
            match_frac: 0.015,
            color_cardinality: 500.0,
            gc_content: 0.5,
        };

        let good_score = scorer.score(&good_seed);
        let bad_score = scorer.score(&bad_seed);

        assert!(good_score > bad_score,
            "good seed score ({:.4}) should be > bad seed score ({:.4})",
            good_score, bad_score);
    }

    #[test]
    fn test_weighted_score_scales_with_length() {
        let scorer = SeedScorer::default();

        let short = SeedFeatures {
            match_len: 15.0,
            log_sa_count: 3.0,
            query_pos_frac: 0.5,
            match_frac: 0.015,
            color_cardinality: 5.0,
            gc_content: 0.5,
        };

        let long = SeedFeatures {
            match_len: 50.0,
            log_sa_count: 3.0,
            query_pos_frac: 0.5,
            match_frac: 0.05,
            color_cardinality: 5.0,
            gc_content: 0.5,
        };

        assert!(scorer.weighted_score(&long) > scorer.weighted_score(&short));
    }
}
