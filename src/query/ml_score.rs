/// Learned seed scoring via logistic regression.
///
/// Replaces the raw `match_len` anchor weight with a quality score learned
/// from labeled training data. Features include seed frequency, match length,
/// positional bias, genome specificity, and information content.
///
/// Model is a logistic regression — inference is a single dot product + sigmoid
/// with zero ML dependencies.

use crate::index::color::ColorIndex;
use crate::index::fm::SeedHit;

/// Number of features in the model.
pub const NUM_FEATURES: usize = 10;
/// Total weights = 1 bias + NUM_FEATURES.
pub const NUM_WEIGHTS: usize = NUM_FEATURES + 1;

/// Feature names for TSV header (training data export).
pub const FEATURE_NAMES: [&str; NUM_FEATURES] = [
    "match_len",
    "log_sa_count",
    "query_pos_frac",
    "match_frac",
    "color_cardinality",
    "gc_content",
    "information_content",
    "inverse_sa_count",
    "local_seed_density",
    "seed_uniqueness",
];

/// Per-seed feature vector for ML scoring.
#[derive(Clone, Debug)]
pub struct SeedFeatures {
    pub match_len: f64,
    pub log_sa_count: f64,
    pub query_pos_frac: f64,
    pub match_frac: f64,
    pub color_cardinality: f64,
    pub gc_content: f64,
    /// Information content: log2(total_genomes / color_cardinality).
    /// High for genome-specific unitigs, near-zero for core genome.
    pub information_content: f64,
    /// 1.0 / sa_count — directly captures seed rarity (complements log_sa_count).
    pub inverse_sa_count: f64,
    /// Number of other seeds within ±500bp on the query, normalized by query length.
    /// Dense seed clusters indicate high-confidence genomic regions.
    pub local_seed_density: f64,
    /// Seed uniqueness: match_len / (sa_count * color_cardinality).
    /// Captures the combined specificity of length, frequency, and genome scope.
    pub seed_uniqueness: f64,
}

impl SeedFeatures {
    /// Compute features for a seed hit.
    ///
    /// `nearby_seeds` is the list of all seed query positions for computing
    /// local density. Pass an empty slice to skip density computation.
    pub fn from_seed(
        seed: &SeedHit,
        query_len: usize,
        query_bytes: &[u8],
        color_index: &ColorIndex,
    ) -> Self {
        Self::from_seed_with_context(seed, query_len, query_bytes, color_index, &[])
    }

    /// Compute features with full context (including nearby seeds for density).
    pub fn from_seed_with_context(
        seed: &SeedHit,
        query_len: usize,
        query_bytes: &[u8],
        color_index: &ColorIndex,
        all_seed_positions: &[usize],
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
        let color_card_raw = color_index
            .get_colors(seed.unitig_id)
            .map(|c| c.len())
            .unwrap_or(0);
        let color_cardinality = color_card_raw as f64;

        // GC content of the seed region in the query
        let gc_content = if seed.match_len > 0 && seed.query_pos + seed.match_len <= query_bytes.len() {
            let region = &query_bytes[seed.query_pos..seed.query_pos + seed.match_len];
            let gc = region.iter().filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c')).count();
            gc as f64 / seed.match_len as f64
        } else {
            0.5
        };

        // Information content: log2(total_genomes / color_cardinality)
        let total_genomes = color_index.num_genomes();
        let information_content = if color_card_raw > 0 && total_genomes > 0 {
            (total_genomes as f64 / color_card_raw as f64).log2()
        } else {
            0.0
        };

        // Inverse SA count: 1/sa_count (capped at 1.0 for sa_count=1)
        let inverse_sa_count = if seed.sa_count > 0 {
            1.0 / seed.sa_count as f64
        } else {
            0.0
        };

        // Local seed density: count seeds within ±500bp window
        let local_seed_density = if !all_seed_positions.is_empty() && query_len > 0 {
            let window = 500;
            let pos = seed.query_pos;
            let lo = pos.saturating_sub(window);
            let hi = (pos + window).min(query_len);
            let count = all_seed_positions.iter().filter(|&&p| p >= lo && p <= hi).count();
            // Normalize: seeds per kilobase
            let span = (hi - lo) as f64;
            if span > 0.0 { count as f64 / (span / 1000.0) } else { 0.0 }
        } else {
            0.0
        };

        // Seed uniqueness: match_len / (sa_count * color_cardinality)
        let denom = (seed.sa_count as f64) * color_cardinality.max(1.0);
        let seed_uniqueness = if denom > 0.0 { match_len / denom } else { 0.0 };

        Self {
            match_len,
            log_sa_count,
            query_pos_frac,
            match_frac,
            color_cardinality,
            gc_content,
            information_content,
            inverse_sa_count,
            local_seed_density,
            seed_uniqueness,
        }
    }

    /// Convert to a fixed-size array for dot product.
    #[inline]
    pub fn as_array(&self) -> [f64; NUM_FEATURES] {
        [
            self.match_len,
            self.log_sa_count,
            self.query_pos_frac,
            self.match_frac,
            self.color_cardinality,
            self.gc_content,
            self.information_content,
            self.inverse_sa_count,
            self.local_seed_density,
            self.seed_uniqueness,
        ]
    }

    /// Feature names (for TSV header).
    pub fn header() -> String {
        FEATURE_NAMES.join("\t")
    }

    /// Feature values as TSV string.
    pub fn as_tsv(&self) -> String {
        let arr = self.as_array();
        arr.iter()
            .map(|v| format!("{:.6}", v))
            .collect::<Vec<_>>()
            .join("\t")
    }
}

/// Logistic regression seed scorer.
///
/// Weights: [bias, match_len, log_sa_count, query_pos_frac, match_frac,
///           color_cardinality, gc_content, information_content,
///           inverse_sa_count, local_seed_density, seed_uniqueness]
pub struct SeedScorer {
    pub weights: [f64; NUM_WEIGHTS],
}

/// Default weights encode domain knowledge:
/// - Longer seeds are more informative (positive weight)
/// - High-frequency seeds are likely repetitive (negative weight on log_sa_count)
/// - Seeds in many genomes are less specific (negative weight on color_cardinality)
/// - High information content means genome-specific (strong positive)
/// - Dense seed clusters indicate reliable regions (positive)
/// - High seed uniqueness combines length + rarity (positive)
const DEFAULT_WEIGHTS: [f64; NUM_WEIGHTS] = [
    0.5,    // bias
    0.15,   // match_len: longer = better
    -0.3,   // log_sa_count: frequent = likely repetitive
    0.0,    // query_pos_frac: neutral (no positional bias)
    0.1,    // match_frac: covers more query = better
    -0.2,   // color_cardinality: in many genomes = less specific
    -0.05,  // gc_content: extreme GC = slightly worse
    0.3,    // information_content: genome-specific = much better
    0.1,    // inverse_sa_count: rare seeds = better
    0.05,   // local_seed_density: dense clusters = more reliable
    0.2,    // seed_uniqueness: long + rare + specific = best
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
    ///
    /// Supports both legacy 7-weight files (padded with zeros for new features)
    /// and current 11-weight files.
    pub fn load_or_default(index_dir: &std::path::Path) -> Self {
        let weights_path = index_dir.join("seed_scorer.json");
        if weights_path.exists() {
            match std::fs::read_to_string(&weights_path) {
                Ok(contents) => {
                    match serde_json::from_str::<Vec<f64>>(&contents) {
                        Ok(w) if w.len() == NUM_WEIGHTS => {
                            log::info!("Loaded ML seed scorer weights ({} weights) from {:?}", NUM_WEIGHTS, weights_path);
                            let mut weights = [0.0f64; NUM_WEIGHTS];
                            weights.copy_from_slice(&w);
                            return Self { weights };
                        }
                        Ok(w) if w.len() == 7 => {
                            // Legacy 7-weight format: pad new features with default weights
                            log::info!("Loaded legacy 7-weight seed scorer from {:?}, using defaults for new features", weights_path);
                            let mut weights = DEFAULT_WEIGHTS;
                            weights[..7].copy_from_slice(&w);
                            return Self { weights };
                        }
                        Ok(w) => {
                            log::warn!(
                                "seed_scorer.json has {} weights (expected {}), using defaults",
                                w.len(), NUM_WEIGHTS
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
        for i in 0..NUM_FEATURES {
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

        // A long, rare, genome-specific seed should score high
        let good_seed = SeedFeatures {
            match_len: 50.0,
            log_sa_count: 1.0,  // only 2 occurrences
            query_pos_frac: 0.5,
            match_frac: 0.05,
            color_cardinality: 1.0,
            gc_content: 0.5,
            information_content: 13.0,  // very specific
            inverse_sa_count: 0.5,
            local_seed_density: 10.0,
            seed_uniqueness: 25.0,
        };

        // A short, frequent, core-genome seed should score lower
        let bad_seed = SeedFeatures {
            match_len: 15.0,
            log_sa_count: 13.0,  // ~8000 occurrences
            query_pos_frac: 0.5,
            match_frac: 0.015,
            color_cardinality: 500.0,
            gc_content: 0.5,
            information_content: 1.0,  // very common
            inverse_sa_count: 0.000125,
            local_seed_density: 2.0,
            seed_uniqueness: 0.00003,
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
            information_content: 8.0,
            inverse_sa_count: 0.125,
            local_seed_density: 5.0,
            seed_uniqueness: 0.375,
        };

        let long = SeedFeatures {
            match_len: 50.0,
            log_sa_count: 3.0,
            query_pos_frac: 0.5,
            match_frac: 0.05,
            color_cardinality: 5.0,
            gc_content: 0.5,
            information_content: 8.0,
            inverse_sa_count: 0.125,
            local_seed_density: 5.0,
            seed_uniqueness: 1.25,
        };

        assert!(scorer.weighted_score(&long) > scorer.weighted_score(&short));
    }

    #[test]
    fn test_information_content_boosts_specific_seeds() {
        let scorer = SeedScorer::default();

        // Same seed length/frequency, but one is genome-specific
        let specific = SeedFeatures {
            match_len: 31.0,
            log_sa_count: 3.0,
            query_pos_frac: 0.5,
            match_frac: 0.03,
            color_cardinality: 1.0,
            gc_content: 0.5,
            information_content: 14.0,  // unique to 1 genome out of ~16K
            inverse_sa_count: 0.125,
            local_seed_density: 5.0,
            seed_uniqueness: 3.875,
        };

        let shared = SeedFeatures {
            match_len: 31.0,
            log_sa_count: 3.0,
            query_pos_frac: 0.5,
            match_frac: 0.03,
            color_cardinality: 500.0,
            gc_content: 0.5,
            information_content: 5.0,  // shared across many genomes
            inverse_sa_count: 0.125,
            local_seed_density: 5.0,
            seed_uniqueness: 0.00775,
        };

        assert!(scorer.score(&specific) > scorer.score(&shared),
            "genome-specific seed ({:.4}) should score higher than shared ({:.4})",
            scorer.score(&specific), scorer.score(&shared));
    }

    #[test]
    fn test_feature_array_length() {
        assert_eq!(DEFAULT_WEIGHTS.len(), NUM_WEIGHTS);
        assert_eq!(FEATURE_NAMES.len(), NUM_FEATURES);
    }
}
