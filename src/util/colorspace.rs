/// SOLiD-style color-space encoding for DNA sequences.
///
/// Encodes dinucleotide transitions as one of 4 colors (0-3). A sequence of
/// length N becomes: 1 anchor base + (N-1) colors. This representation has
/// two key mathematical properties useful for sequence alignment:
///
/// 1. **Error-correcting geometry**: A single substitution in the original
///    DNA changes exactly 2 consecutive colors. A sequencing error changes
///    only 1 color. This lets us distinguish biological variation from noise.
///
/// 2. **Reduced seed ambiguity**: Divergent sequences (e.g., 10% mutated)
///    produce color-space seeds that often still match despite point
///    mutations, because a SNP locally perturbs only 2 adjacent colors.
///
/// The SOLiD color matrix (row = first base, column = second base):
///
/// ```text
///      A  C  G  T
///   A  0  1  2  3
///   C  1  0  3  2
///   G  2  3  0  1
///   T  3  2  1  0
/// ```
///
/// Properties:
/// - Each color has 4 dinucleotide pre-images (ambiguous without anchor)
/// - Color 0 = identity (AA, CC, GG, TT)
/// - Color 3 = transversion (AT, TA, CG, GC)
/// - Reverse complement preserves the color sequence (reversed)

/// Encode a nucleotide pair into a SOLiD color (0-3).
#[inline]
pub fn encode_pair(a: u8, b: u8) -> Option<u8> {
    let a_idx = base_to_idx(a)?;
    let b_idx = base_to_idx(b)?;
    Some(COLOR_MATRIX[a_idx][b_idx])
}

/// Decode a color + anchor base into the next nucleotide.
#[inline]
pub fn decode_pair(anchor: u8, color: u8) -> Option<u8> {
    let a_idx = base_to_idx(anchor)?;
    if color >= 4 {
        return None;
    }
    // Find b such that COLOR_MATRIX[a_idx][b] == color
    for b_idx in 0..4 {
        if COLOR_MATRIX[a_idx][b_idx] == color {
            return Some(IDX_TO_BASE[b_idx]);
        }
    }
    None
}

/// Encode a DNA sequence to color-space: returns (anchor_base, color_string).
///
/// The color_string has length (seq.len() - 1). N bases produce a gap (color 4).
pub fn encode_sequence(seq: &[u8]) -> Option<(u8, Vec<u8>)> {
    if seq.len() < 2 {
        return None;
    }

    let anchor = seq[0].to_ascii_uppercase();
    if !is_valid_base(anchor) {
        return None;
    }

    let mut colors = Vec::with_capacity(seq.len() - 1);
    for window in seq.windows(2) {
        let a = window[0].to_ascii_uppercase();
        let b = window[1].to_ascii_uppercase();
        match encode_pair(a, b) {
            Some(c) => colors.push(c),
            None => colors.push(4), // gap/N marker
        }
    }

    Some((anchor, colors))
}

/// Decode a color-space sequence back to DNA.
pub fn decode_sequence(anchor: u8, colors: &[u8]) -> Option<Vec<u8>> {
    let mut result = Vec::with_capacity(colors.len() + 1);
    result.push(anchor);

    let mut current = anchor;
    for &color in colors {
        if color == 4 {
            return None; // gap
        }
        match decode_pair(current, color) {
            Some(next) => {
                result.push(next);
                current = next;
            }
            None => return None,
        }
    }

    Some(result)
}

/// Count the number of color differences between two color-space sequences.
///
/// A single SNP in DNA produces exactly 2 consecutive color differences.
/// A sequencing error produces exactly 1. This lets us detect biological
/// variation vs. noise using the topology of differences.
pub fn color_differences(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Classify a set of color differences as SNPs or errors based on adjacency.
///
/// Returns (num_snps, num_errors): an isolated color difference is an error,
/// a pair of consecutive differences is a SNP.
pub fn classify_differences(a: &[u8], b: &[u8]) -> (usize, usize) {
    let diffs: Vec<bool> = a.iter().zip(b.iter()).map(|(x, y)| x != y).collect();
    let mut snps = 0;
    let mut errors = 0;
    let mut i = 0;
    while i < diffs.len() {
        if diffs[i] {
            if i + 1 < diffs.len() && diffs[i + 1] {
                snps += 1;
                i += 2;
            } else {
                errors += 1;
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    (snps, errors)
}

// ---- Internals ----

/// The SOLiD color matrix: rows = first base, columns = second base.
const COLOR_MATRIX: [[u8; 4]; 4] = [
    [0, 1, 2, 3], // A -> A, C, G, T
    [1, 0, 3, 2], // C -> A, C, G, T
    [2, 3, 0, 1], // G -> A, C, G, T
    [3, 2, 1, 0], // T -> A, C, G, T
];

const IDX_TO_BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_color_matrix_identity() {
        // Same base → color 0
        assert_eq!(encode_pair(b'A', b'A'), Some(0));
        assert_eq!(encode_pair(b'C', b'C'), Some(0));
        assert_eq!(encode_pair(b'G', b'G'), Some(0));
        assert_eq!(encode_pair(b'T', b'T'), Some(0));
    }

    #[test]
    fn test_color_matrix_transversions() {
        // Transversions → color 3
        assert_eq!(encode_pair(b'A', b'T'), Some(3));
        assert_eq!(encode_pair(b'T', b'A'), Some(3));
        assert_eq!(encode_pair(b'C', b'G'), Some(3));
        assert_eq!(encode_pair(b'G', b'C'), Some(3));
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        let dna = b"ACGTACGTACGT";
        let (anchor, colors) = encode_sequence(dna).unwrap();
        let decoded = decode_sequence(anchor, &colors).unwrap();
        assert_eq!(decoded, dna);
    }

    #[test]
    fn test_snp_produces_two_color_changes() {
        let dna1 = b"ACGTACGT";
        let dna2 = b"ACGTTCGT"; // single SNP at position 4 (A→T)
        let (_, colors1) = encode_sequence(dna1).unwrap();
        let (_, colors2) = encode_sequence(dna2).unwrap();

        // Count differences
        let diffs = color_differences(&colors1, &colors2);
        assert_eq!(diffs, 2, "Single SNP should produce exactly 2 color differences");
    }

    #[test]
    fn test_sequencing_error_produces_one_color_change() {
        // Simulate a single sequencing error (one color flipped)
        let (anchor, mut colors) = encode_sequence(b"ACGTACGT").unwrap();
        let original = colors.clone();
        colors[3] = (colors[3] + 1) % 4; // flip one color

        let (snps, errors) = classify_differences(&original, &colors);
        assert_eq!(errors, 1);
        assert_eq!(snps, 0);
    }

    #[test]
    fn test_classify_snp_vs_error() {
        // Two consecutive color changes = 1 SNP
        let colors1 = vec![0, 0, 0, 0, 0];
        let colors2 = vec![0, 1, 2, 0, 0]; // positions 1,2 changed = SNP

        let (snps, errors) = classify_differences(&colors1, &colors2);
        assert_eq!(snps, 1);
        assert_eq!(errors, 0);
    }

    #[test]
    fn test_classify_two_separate_errors() {
        let colors1 = vec![0, 0, 0, 0, 0];
        let colors2 = vec![1, 0, 2, 0, 0]; // two isolated changes

        let (snps, errors) = classify_differences(&colors1, &colors2);
        assert_eq!(snps, 0);
        assert_eq!(errors, 2);
    }

    #[test]
    fn test_decode_pair_coverage() {
        // Every (anchor, color) pair should decode uniquely
        for &anchor in &[b'A', b'C', b'G', b'T'] {
            for color in 0..4 {
                let b = decode_pair(anchor, color).unwrap();
                let c = encode_pair(anchor, b).unwrap();
                assert_eq!(c, color);
            }
        }
    }
}
