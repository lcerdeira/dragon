/// 2-bit DNA encoding: A=0, C=1, G=2, T=3
/// Packs 32 bases per u64.

/// Encode a single base to 2-bit representation.
#[inline]
pub fn encode_base(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0, // N and ambiguous bases default to A
    }
}

/// Decode a 2-bit value back to ASCII base.
#[inline]
pub fn decode_base(val: u8) -> u8 {
    match val & 0x03 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!(),
    }
}

/// Complement of a 2-bit encoded base.
#[inline]
pub fn complement_2bit(val: u8) -> u8 {
    3 - (val & 0x03)
}

/// A 2-bit packed DNA sequence.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct PackedSequence {
    pub data: Vec<u64>,
    pub len: usize,
}

impl PackedSequence {
    /// Encode an ASCII DNA sequence into 2-bit packed format.
    pub fn from_bytes(seq: &[u8]) -> Self {
        let num_words = (seq.len() + 31) / 32;
        let mut data = vec![0u64; num_words];

        for (i, &base) in seq.iter().enumerate() {
            let word_idx = i / 32;
            let bit_offset = (i % 32) * 2;
            data[word_idx] |= (encode_base(base) as u64) << bit_offset;
        }

        Self {
            data,
            len: seq.len(),
        }
    }

    /// Decode the packed sequence back to ASCII bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::with_capacity(self.len);
        for i in 0..self.len {
            let word_idx = i / 32;
            let bit_offset = (i % 32) * 2;
            let val = ((self.data[word_idx] >> bit_offset) & 0x03) as u8;
            result.push(decode_base(val));
        }
        result
    }

    /// Get the 2-bit encoded value at position i.
    #[inline]
    pub fn get(&self, i: usize) -> u8 {
        let word_idx = i / 32;
        let bit_offset = (i % 32) * 2;
        ((self.data[word_idx] >> bit_offset) & 0x03) as u8
    }

    /// Compute the reverse complement of this packed sequence.
    pub fn reverse_complement(&self) -> Self {
        let num_words = (self.len + 31) / 32;
        let mut data = vec![0u64; num_words];

        for i in 0..self.len {
            let src_val = self.get(i);
            let dest_pos = self.len - 1 - i;
            let word_idx = dest_pos / 32;
            let bit_offset = (dest_pos % 32) * 2;
            data[word_idx] |= (complement_2bit(src_val) as u64) << bit_offset;
        }

        Self {
            data,
            len: self.len,
        }
    }

    /// Extract a subsequence [start, end) as ASCII bytes.
    pub fn subsequence(&self, start: usize, end: usize) -> Vec<u8> {
        let mut result = Vec::with_capacity(end - start);
        for i in start..end.min(self.len) {
            result.push(decode_base(self.get(i)));
        }
        result
    }

    /// Extract a k-mer starting at position pos as a u64 (k <= 32).
    pub fn kmer_u64(&self, pos: usize, k: usize) -> u64 {
        assert!(k <= 32);
        let mut kmer: u64 = 0;
        for i in 0..k {
            kmer |= (self.get(pos + i) as u64) << (i * 2);
        }
        kmer
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// Compute canonical k-mer (lexicographically smaller of forward and revcomp).
pub fn canonical_kmer(kmer: u64, k: usize) -> u64 {
    let rc = revcomp_kmer(kmer, k);
    kmer.min(rc)
}

/// Reverse complement a k-mer stored as u64.
pub fn revcomp_kmer(kmer: u64, k: usize) -> u64 {
    let mut rc: u64 = 0;
    for i in 0..k {
        let base = (kmer >> (i * 2)) & 0x03;
        let comp = 3 - base;
        rc |= comp << ((k - 1 - i) * 2);
    }
    rc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let seq = b"ACGTACGTNNACGT";
        let packed = PackedSequence::from_bytes(seq);
        let decoded = packed.to_bytes();
        // N gets encoded as A
        assert_eq!(&decoded[..4], b"ACGT");
        assert_eq!(&decoded[8..10], b"AA"); // NN -> AA
        assert_eq!(packed.len, 14);
    }

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        let packed = PackedSequence::from_bytes(seq);
        let rc = packed.reverse_complement();
        assert_eq!(rc.to_bytes(), b"ACGT"); // ACGT is its own revcomp

        let seq2 = b"AACG";
        let packed2 = PackedSequence::from_bytes(seq2);
        let rc2 = packed2.reverse_complement();
        assert_eq!(rc2.to_bytes(), b"CGTT");
    }

    #[test]
    fn test_subsequence() {
        let seq = b"ACGTACGTACGT";
        let packed = PackedSequence::from_bytes(seq);
        assert_eq!(packed.subsequence(2, 6), b"GTAC");
    }

    #[test]
    fn test_canonical_kmer() {
        let kmer = 0b00_01_10_11u64; // ACGT
        let rc = revcomp_kmer(kmer, 4);
        assert_eq!(kmer, rc); // ACGT is palindromic
    }

    #[test]
    fn test_long_sequence() {
        let seq: Vec<u8> = (0..100).map(|i| b"ACGT"[i % 4]).collect();
        let packed = PackedSequence::from_bytes(&seq);
        assert_eq!(packed.to_bytes(), seq);
    }
}
