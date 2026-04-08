/// Signal discretization: converts continuous pA current values into a finite
/// discrete alphabet suitable for FM-index construction and backward search.
///
/// The discretization uses equal-width bins over the normalized signal range,
/// mapping each current measurement to one of `num_levels` symbols (default 16).

use serde::{Deserialize, Serialize};

/// Defines the discrete alphabet for signal encoding.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SignalAlphabet {
    /// Number of discrete levels (symbols in the alphabet).
    pub num_levels: u8,
    /// Minimum normalized signal value (bin floor).
    pub min_val: f32,
    /// Maximum normalized signal value (bin ceiling).
    pub max_val: f32,
    /// Optional learned boundaries for non-uniform binning.
    /// When `Some`, uses binary search over sorted thresholds instead of equal-width bins.
    /// The vector should contain `num_levels - 1` boundary values in ascending order.
    #[serde(default)]
    pub boundaries: Option<Vec<f32>>,
}

impl Default for SignalAlphabet {
    fn default() -> Self {
        Self {
            num_levels: 16,
            // Typical normalized signal range after median-MAD normalization
            min_val: -4.0,
            max_val: 4.0,
            boundaries: None,
        }
    }
}

impl SignalAlphabet {
    /// Create a new alphabet with the given number of levels.
    pub fn new(num_levels: u8) -> Self {
        Self {
            num_levels,
            ..Default::default()
        }
    }

    /// Create an alphabet with custom range.
    pub fn with_range(num_levels: u8, min_val: f32, max_val: f32) -> Self {
        assert!(max_val > min_val, "max_val must be greater than min_val");
        assert!(num_levels >= 2, "need at least 2 levels");
        Self {
            num_levels,
            min_val,
            max_val,
            boundaries: None,
        }
    }

    /// Create an alphabet from learned boundary thresholds.
    ///
    /// `boundaries` should be a sorted vector of `num_levels - 1` threshold values.
    /// A value `v` is assigned to bin `i` if `boundaries[i-1] <= v < boundaries[i]`
    /// (with implicit -inf and +inf at the ends).
    pub fn from_boundaries(boundaries: Vec<f32>) -> Self {
        assert!(
            boundaries.len() >= 1,
            "need at least 1 boundary for 2+ levels"
        );
        // Verify sorted
        for i in 1..boundaries.len() {
            assert!(
                boundaries[i] >= boundaries[i - 1],
                "boundaries must be sorted"
            );
        }
        let num_levels = (boundaries.len() + 1) as u8;
        let min_val = boundaries.first().copied().unwrap_or(-4.0) - 1.0;
        let max_val = boundaries.last().copied().unwrap_or(4.0) + 1.0;
        Self {
            num_levels,
            min_val,
            max_val,
            boundaries: Some(boundaries),
        }
    }

    /// Width of each bin.
    #[inline]
    pub fn bin_width(&self) -> f32 {
        (self.max_val - self.min_val) / self.num_levels as f32
    }

    /// Discretize a single normalized value to a symbol in [0, num_levels).
    #[inline]
    pub fn discretize_value(&self, val: f32) -> u8 {
        if let Some(ref bounds) = self.boundaries {
            // Learned boundaries: binary search
            match bounds.binary_search_by(|b| b.partial_cmp(&val).unwrap_or(std::cmp::Ordering::Equal)) {
                Ok(i) => (i as u8 + 1).min(self.num_levels - 1),
                Err(i) => (i as u8).min(self.num_levels - 1),
            }
        } else {
            // Equal-width bins
            let clamped = val.clamp(self.min_val, self.max_val - f32::EPSILON);
            let bin = ((clamped - self.min_val) / self.bin_width()) as u8;
            bin.min(self.num_levels - 1)
        }
    }

    /// Map a symbol back to the center of its bin (for display/debug).
    #[inline]
    pub fn symbol_to_pa(&self, symbol: u8) -> f32 {
        let width = self.bin_width();
        self.min_val + (symbol as f32 + 0.5) * width
    }
}

/// Normalize a raw pA signal using median-MAD (Median Absolute Deviation).
///
/// This is the standard normalization for nanopore signals:
///   normalized = (raw - median) / (1.4826 * MAD)
///
/// where MAD = median(|x - median(x)|) and 1.4826 is the consistency constant
/// for estimating standard deviation from MAD under normality.
pub fn normalize_signal(raw: &[f32]) -> Vec<f32> {
    if raw.is_empty() {
        return Vec::new();
    }

    // Compute median
    let mut sorted = raw.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    // Compute MAD (Median Absolute Deviation)
    let mut abs_devs: Vec<f32> = raw.iter().map(|&x| (x - median).abs()).collect();
    abs_devs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = if abs_devs.len() % 2 == 0 {
        (abs_devs[abs_devs.len() / 2 - 1] + abs_devs[abs_devs.len() / 2]) / 2.0
    } else {
        abs_devs[abs_devs.len() / 2]
    };

    // Scale factor: 1.4826 * MAD (consistency constant for Gaussian)
    let scale = 1.4826 * mad;

    // Avoid division by zero for constant signals
    if scale < 1e-10 {
        return vec![0.0; raw.len()];
    }

    raw.iter().map(|&x| (x - median) / scale).collect()
}

/// Discretize a raw (unnormalized) signal into a discrete alphabet.
/// Applies median-MAD normalization first, then bins into `alphabet.num_levels` symbols.
pub fn discretize_signal(raw: &[f32], alphabet: &SignalAlphabet) -> Vec<u8> {
    let normalized = normalize_signal(raw);
    normalized
        .iter()
        .map(|&v| alphabet.discretize_value(v))
        .collect()
}

/// Discretize an already-normalized signal.
pub fn discretize_normalized(normalized: &[f32], alphabet: &SignalAlphabet) -> Vec<u8> {
    normalized
        .iter()
        .map(|&v| alphabet.discretize_value(v))
        .collect()
}

/// Extract overlapping k-mers from a discretized signal.
///
/// Each k-mer is a slice of `k` consecutive discrete symbols.
/// Returns the k-mers as byte vectors.
pub fn extract_signal_kmers(discretized: &[u8], k: usize) -> Vec<Vec<u8>> {
    if discretized.len() < k {
        return Vec::new();
    }
    (0..=discretized.len() - k)
        .map(|i| discretized[i..i + k].to_vec())
        .collect()
}

/// Extract signal k-mers as a flat concatenation with separators.
/// Each k-mer is separated by a sentinel symbol (num_levels, which is outside the alphabet).
/// This is useful for building a concatenated text for FM-index construction.
pub fn extract_signal_kmers_concat(
    discretized: &[u8],
    k: usize,
    sentinel: u8,
) -> Vec<u8> {
    if discretized.len() < k {
        return Vec::new();
    }
    let mut result = Vec::with_capacity(discretized.len());
    for i in 0..=discretized.len() - k {
        if i > 0 {
            result.push(sentinel);
        }
        result.extend_from_slice(&discretized[i..i + k]);
    }
    result.push(sentinel);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_alphabet() {
        let alpha = SignalAlphabet::default();
        assert_eq!(alpha.num_levels, 16);
        assert!((alpha.bin_width() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_discretize_value() {
        let alpha = SignalAlphabet::default();
        // Center of range -> middle bin
        assert_eq!(alpha.discretize_value(0.0), 8);
        // Very negative -> bin 0
        assert_eq!(alpha.discretize_value(-10.0), 0);
        // Very positive -> last bin
        assert_eq!(alpha.discretize_value(10.0), 15);
        // Just above min
        assert_eq!(alpha.discretize_value(-3.9), 0);
    }

    #[test]
    fn test_symbol_to_pa_roundtrip() {
        let alpha = SignalAlphabet::default();
        for sym in 0..alpha.num_levels {
            let pa = alpha.symbol_to_pa(sym);
            let back = alpha.discretize_value(pa);
            assert_eq!(back, sym, "symbol {} -> pA {} -> symbol {}", sym, pa, back);
        }
    }

    #[test]
    fn test_normalize_signal() {
        // Constant signal -> all zeros
        let constant = vec![50.0; 100];
        let norm = normalize_signal(&constant);
        assert!(norm.iter().all(|&v| v.abs() < 1e-6));

        // Linear signal
        let linear: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let norm = normalize_signal(&linear);
        // After normalization, median should map to ~0
        let mid_val = norm[49]; // close to median
        assert!(mid_val.abs() < 0.1);
    }

    #[test]
    fn test_discretize_signal() {
        let alpha = SignalAlphabet::new(4);
        let raw: Vec<f32> = (0..100).map(|i| 50.0 + i as f32).collect();
        let disc = discretize_signal(&raw, &alpha);
        assert_eq!(disc.len(), 100);
        // All values should be in [0, 4)
        assert!(disc.iter().all(|&v| v < 4));
    }

    #[test]
    fn test_extract_signal_kmers() {
        let disc = vec![0u8, 1, 2, 3, 4, 5];
        let kmers = extract_signal_kmers(&disc, 3);
        assert_eq!(kmers.len(), 4);
        assert_eq!(kmers[0], vec![0, 1, 2]);
        assert_eq!(kmers[1], vec![1, 2, 3]);
        assert_eq!(kmers[3], vec![3, 4, 5]);
    }

    #[test]
    fn test_extract_signal_kmers_too_short() {
        let disc = vec![0u8, 1];
        let kmers = extract_signal_kmers(&disc, 5);
        assert!(kmers.is_empty());
    }

    #[test]
    fn test_from_boundaries() {
        // 4 levels with 3 boundaries: [-1.0, 0.0, 1.0]
        let alpha = SignalAlphabet::from_boundaries(vec![-1.0, 0.0, 1.0]);
        assert_eq!(alpha.num_levels, 4);

        // Values below first boundary -> bin 0
        assert_eq!(alpha.discretize_value(-5.0), 0);
        assert_eq!(alpha.discretize_value(-1.5), 0);

        // Between -1.0 and 0.0 -> bin 1
        assert_eq!(alpha.discretize_value(-0.5), 1);

        // Between 0.0 and 1.0 -> bin 2
        assert_eq!(alpha.discretize_value(0.5), 2);

        // Above last boundary -> bin 3
        assert_eq!(alpha.discretize_value(2.0), 3);
    }

    #[test]
    fn test_from_boundaries_roundtrip_serde() {
        let alpha = SignalAlphabet::from_boundaries(vec![-2.0, -0.5, 0.5, 2.0]);
        let json = serde_json::to_string(&alpha).unwrap();
        let alpha2: SignalAlphabet = serde_json::from_str(&json).unwrap();
        assert_eq!(alpha2.num_levels, 5);
        assert!(alpha2.boundaries.is_some());
        assert_eq!(alpha2.boundaries.unwrap().len(), 4);
    }
}
