/// Pore model: maps DNA k-mers to expected pA (picoampere) current levels.
///
/// This provides an approximate R10.4.1 pore model mapping each 5-mer to
/// an expected current level. The model is used to convert reference genome
/// sequences into expected signal space for signal-level indexing.
///
/// The 1024 5-mer entries are generated from a simplified R10.4.1 model
/// with base-level contributions and nearest-neighbor interaction terms.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A pore model mapping DNA k-mers to expected pA current levels.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PoreModel {
    /// K-mer size (typically 5 or 6).
    pub kmer_size: usize,
    /// Lookup table: index = k-mer encoded as base-4 number (A=0,C=1,G=2,T=3),
    /// value = expected pA level.
    pub levels: Vec<f32>,
    /// Model name/description.
    pub name: String,
}

impl PoreModel {
    /// Look up expected pA level for a DNA k-mer (ASCII bytes: A,C,G,T).
    /// Returns None if the k-mer contains non-ACGT bases or wrong length.
    pub fn expected_signal(&self, kmer: &[u8]) -> Option<f32> {
        if kmer.len() != self.kmer_size {
            return None;
        }
        let idx = kmer_to_index(kmer)?;
        Some(self.levels[idx])
    }

    /// Create a pore model from a pre-computed lookup table.
    ///
    /// `levels` should contain `4^kmer_size` entries mapping each k-mer index
    /// to an expected pA level.
    pub fn from_table(kmer_size: usize, levels: Vec<f32>, name: String) -> Self {
        let expected = 4usize.pow(kmer_size as u32);
        assert_eq!(
            levels.len(),
            expected,
            "levels table must have {} entries for k={}, got {}",
            expected,
            kmer_size,
            levels.len()
        );
        Self {
            kmer_size,
            levels,
            name,
        }
    }

    /// Load a pore model from a JSON file.
    ///
    /// Expected format: `{"kmer_size": 5, "levels": [f32; 1024], "name": "..."}`
    /// or just a flat array of 1024 floats (assumes k=5).
    pub fn load_from_file(path: &std::path::Path) -> Option<Self> {
        let data = std::fs::read_to_string(path).ok()?;

        // Try full JSON object first
        if let Ok(model) = serde_json::from_str::<PoreModel>(&data) {
            log::info!("Loaded pore model '{}' from {:?}", model.name, path);
            return Some(model);
        }

        // Try flat array (assume k=5)
        if let Ok(levels) = serde_json::from_str::<Vec<f32>>(&data) {
            if levels.len() == 1024 {
                log::info!("Loaded flat pore model table ({} entries) from {:?}", levels.len(), path);
                return Some(Self::from_table(
                    5,
                    levels,
                    format!("custom-from-{}", path.display()),
                ));
            }
        }

        log::warn!("Failed to parse pore model from {:?}", path);
        None
    }

    /// Convert an entire DNA sequence to expected signal.
    /// Returns one pA value per k-mer position (len - kmer_size + 1 values).
    pub fn sequence_to_expected_signal(&self, seq: &[u8]) -> Vec<f32> {
        if seq.len() < self.kmer_size {
            return Vec::new();
        }
        let mut signal = Vec::with_capacity(seq.len() - self.kmer_size + 1);
        for i in 0..=seq.len() - self.kmer_size {
            let kmer = &seq[i..i + self.kmer_size];
            match self.expected_signal(kmer) {
                Some(level) => signal.push(level),
                None => signal.push(90.0), // default for N-containing k-mers
            }
        }
        signal
    }

    /// Number of entries in the model (4^k).
    pub fn num_entries(&self) -> usize {
        self.levels.len()
    }

    /// Get a hashmap from k-mer string to pA level (for debugging/export).
    pub fn to_hashmap(&self) -> HashMap<String, f32> {
        let mut map = HashMap::new();
        let num_kmers = 4usize.pow(self.kmer_size as u32);
        for idx in 0..num_kmers {
            let kmer = index_to_kmer(idx, self.kmer_size);
            let kmer_str = String::from_utf8(kmer).unwrap();
            map.insert(kmer_str, self.levels[idx]);
        }
        map
    }
}

/// Encode a DNA k-mer (ASCII) to a numeric index.
/// A=0, C=1, G=2, T=3; the first base is the most significant digit.
fn kmer_to_index(kmer: &[u8]) -> Option<usize> {
    let mut idx = 0usize;
    for &base in kmer {
        let val = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        idx = idx * 4 + val;
    }
    Some(idx)
}

/// Decode a numeric index back to a DNA k-mer.
fn index_to_kmer(mut idx: usize, k: usize) -> Vec<u8> {
    let mut kmer = vec![b'A'; k];
    for i in (0..k).rev() {
        kmer[i] = match idx % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        idx /= 4;
    }
    kmer
}

/// Load the default approximate R10.4.1 pore model (5-mer).
///
/// This model uses a simplified parametric approach:
///   level(kmer) = sum of base contributions + nearest-neighbor terms + noise
///
/// Base contributions (approximate pA offset from ~90 pA baseline):
///   A: -5.0, C: +3.0, G: +8.0, T: -2.0
///
/// Nearest-neighbor dinucleotide interaction terms are also included.
/// The resulting levels span approximately 60-130 pA, which matches the
/// typical range for R10.4.1 chemistry.
pub fn load_default_model() -> PoreModel {
    let k = 5;
    let num_kmers = 4usize.pow(k as u32); // 1024

    // Base-level contributions to current (pA)
    let base_contrib: [f32; 4] = [-5.0, 3.0, 8.0, -2.0]; // A, C, G, T

    // Nearest-neighbor dinucleotide interaction terms (4x4 matrix)
    // Rows = first base, Cols = second base (A=0,C=1,G=2,T=3)
    let nn_interact: [[f32; 4]; 4] = [
        // A->A  A->C  A->G  A->T
        [-1.2, 0.8, 2.1, -0.5],
        // C->A  C->C  C->G  C->T
        [0.5, -0.9, 1.5, 0.3],
        // G->A  G->C  G->G  G->T
        [1.8, 0.2, -1.1, 1.0],
        // T->A  T->C  T->G  T->T
        [-0.3, 0.7, 1.3, -1.5],
    ];

    // Position weights (central positions contribute more)
    let pos_weights: [f32; 5] = [0.6, 0.8, 1.0, 0.8, 0.6];

    let baseline = 90.0f32;

    let mut levels = Vec::with_capacity(num_kmers);

    for idx in 0..num_kmers {
        let kmer = index_to_kmer(idx, k);
        let mut level = baseline;

        // Add weighted base contributions
        for (pos, &base) in kmer.iter().enumerate() {
            let base_idx = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0,
            };
            level += base_contrib[base_idx] * pos_weights[pos];
        }

        // Add nearest-neighbor interactions
        for i in 0..k - 1 {
            let b1 = match kmer[i] {
                b'A' => 0usize,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0,
            };
            let b2 = match kmer[i + 1] {
                b'A' => 0usize,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0,
            };
            level += nn_interact[b1][b2];
        }

        levels.push(level);
    }

    PoreModel {
        kmer_size: k,
        levels,
        name: "R10.4.1-approx-5mer".to_string(),
    }
}

/// Compute the complement of a DNA k-mer.
pub fn complement_kmer(kmer: &[u8]) -> Vec<u8> {
    kmer.iter()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Compute the reverse complement of a DNA k-mer.
pub fn reverse_complement_kmer(kmer: &[u8]) -> Vec<u8> {
    kmer.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_index_roundtrip() {
        for k in 3..=6 {
            let num = 4usize.pow(k as u32);
            for idx in 0..num {
                let kmer = index_to_kmer(idx, k);
                let back = kmer_to_index(&kmer).unwrap();
                assert_eq!(idx, back, "k={}, idx={}, kmer={:?}", k, idx, kmer);
            }
        }
    }

    #[test]
    fn test_default_model_size() {
        let model = load_default_model();
        assert_eq!(model.kmer_size, 5);
        assert_eq!(model.num_entries(), 1024);
    }

    #[test]
    fn test_expected_signal_range() {
        let model = load_default_model();
        // All levels should be in a reasonable pA range (roughly 50-140 pA)
        for &level in &model.levels {
            assert!(
                level > 50.0 && level < 140.0,
                "level {} out of expected range",
                level
            );
        }
    }

    #[test]
    fn test_expected_signal_lookup() {
        let model = load_default_model();
        let level = model.expected_signal(b"ACGTG").unwrap();
        assert!(level > 50.0 && level < 140.0);

        // Wrong length
        assert!(model.expected_signal(b"ACG").is_none());
        // Invalid base
        assert!(model.expected_signal(b"ACNGT").is_none());
    }

    #[test]
    fn test_sequence_to_expected_signal() {
        let model = load_default_model();
        let seq = b"ACGTACGTACGT";
        let signal = model.sequence_to_expected_signal(seq);
        // For a 12-base sequence with k=5, expect 8 signal values
        assert_eq!(signal.len(), 8);
        // All values in range
        assert!(signal.iter().all(|&v| v > 50.0 && v < 140.0));
    }

    #[test]
    fn test_different_kmers_different_levels() {
        let model = load_default_model();
        let l1 = model.expected_signal(b"AAAAA").unwrap();
        let l2 = model.expected_signal(b"GGGGG").unwrap();
        // G-rich should be higher than A-rich
        assert!(l2 > l1, "GGGGG ({}) should be > AAAAA ({})", l2, l1);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement_kmer(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement_kmer(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement_kmer(b"AACG"), b"CGTT");
    }

    #[test]
    fn test_to_hashmap() {
        let model = load_default_model();
        let map = model.to_hashmap();
        assert_eq!(map.len(), 1024);
        assert!(map.contains_key("AAAAA"));
        assert!(map.contains_key("TTTTT"));
    }
}
