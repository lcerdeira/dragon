/// Elias-Fano encoded monotone sequence for mapping FM-index positions
/// back to unitig IDs. Supports O(1) predecessor queries.
///
/// For now, we use a simple sorted array with binary search.
/// This can be upgraded to a true Elias-Fano structure via the `sux` crate
/// for production use.

use serde::{Deserialize, Serialize};

/// Maps positions in the concatenated unitig text to (unitig_id, offset).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CumulativeLengthIndex {
    /// Cumulative start positions: cum_lengths[i] = start position of unitig i.
    /// cum_lengths[num_unitigs] = total text length.
    cum_lengths: Vec<u64>,
}

impl CumulativeLengthIndex {
    /// Build from a slice of unitig lengths.
    pub fn from_lengths(lengths: &[u64]) -> Self {
        let mut cum_lengths = Vec::with_capacity(lengths.len() + 1);
        let mut total = 0u64;
        cum_lengths.push(0);
        for &len in lengths {
            total += len;
            cum_lengths.push(total);
        }
        Self { cum_lengths }
    }

    /// Given a position in the concatenated FM-index text, return
    /// `(unitig_id, offset_within_unitig)`.
    ///
    /// **Separator awareness.** The FM-index text concatenates unitig
    /// sequences with a 1-byte `$` separator *after each unitig* (see
    /// [`crate::index::unitig::UnitigSet`]). `cum_lengths`, however, stores
    /// separator-*free* cumulative sums (it is built from raw unitig
    /// lengths). Therefore unitig `i` begins at text position
    /// `cum_lengths[i] + i` — there are `i` separators in front of it — and
    /// its `$` sits at `cum_lengths[i] + i + len_i`.
    ///
    /// Mapping a text position with the raw `cum_lengths` boundaries (as a
    /// previous version did) yields an offset inflated by the unitig index
    /// — e.g. `+1867` for unitig 1867 — and, for large indices, the wrong
    /// unitig entirely. That mis-anchored every FM seed and was the root
    /// cause of the saureus `mean_id 0.624` ceiling.
    pub fn unitig_at_position(&self, pos: u64) -> Option<(u32, u32)> {
        let n = self.num_unitigs();
        if n == 0 {
            return None;
        }
        // text_start(i) = cum_lengths[i] + i is strictly increasing, so we
        // can binary-search it for the rightmost start <= pos.
        let text_start = |i: usize| self.cum_lengths[i] + i as u64;
        let (mut lo, mut hi) = (0usize, n);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if text_start(mid) <= pos {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo == 0 {
            return None; // pos precedes the first unitig
        }
        let idx = lo - 1;
        let offset = pos - text_start(idx);
        let len = self.cum_lengths[idx + 1] - self.cum_lengths[idx];
        if offset >= len {
            // pos falls on the '$' separator after unitig `idx`.
            return None;
        }
        Some((idx as u32, offset as u32))
    }

    /// Number of unitigs in the index.
    pub fn num_unitigs(&self) -> usize {
        if self.cum_lengths.is_empty() {
            0
        } else {
            self.cum_lengths.len() - 1
        }
    }

    /// Total text length.
    pub fn total_length(&self) -> u64 {
        *self.cum_lengths.last().unwrap_or(&0)
    }

    /// Get the start position of a unitig.
    pub fn unitig_start(&self, unitig_id: u32) -> u64 {
        self.cum_lengths[unitig_id as usize]
    }

    /// Get the length of a unitig.
    pub fn unitig_length(&self, unitig_id: u32) -> u64 {
        let id = unitig_id as usize;
        self.cum_lengths[id + 1] - self.cum_lengths[id]
    }

    /// Reconstruct the individual lengths from cumulative sums.
    pub fn lengths(&self) -> Vec<u64> {
        (0..self.num_unitigs())
            .map(|i| self.cum_lengths[i + 1] - self.cum_lengths[i])
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The FM-index text concatenates unitigs with a 1-byte '$' separator
    // after each. For lengths [100, 200, 150] the text layout is:
    //   unitig 0 : text[0  .. 100)   '$' @ 100
    //   unitig 1 : text[101 .. 301)  '$' @ 301
    //   unitig 2 : text[302 .. 452)  '$' @ 452
    #[test]
    fn test_basic_lookup() {
        let idx = CumulativeLengthIndex::from_lengths(&[100, 200, 150]);
        assert_eq!(idx.num_unitigs(), 3);

        assert_eq!(idx.unitig_at_position(0), Some((0, 0)));
        assert_eq!(idx.unitig_at_position(99), Some((0, 99)));
        assert_eq!(idx.unitig_at_position(100), None); // '$' after unitig 0
        assert_eq!(idx.unitig_at_position(101), Some((1, 0)));
        assert_eq!(idx.unitig_at_position(300), Some((1, 199)));
        assert_eq!(idx.unitig_at_position(301), None); // '$' after unitig 1
        assert_eq!(idx.unitig_at_position(302), Some((2, 0)));
        assert_eq!(idx.unitig_at_position(451), Some((2, 149)));
        assert_eq!(idx.unitig_at_position(452), None); // '$' after unitig 2
        assert_eq!(idx.unitig_at_position(453), None); // past end
    }

    #[test]
    fn test_single_unitig() {
        let idx = CumulativeLengthIndex::from_lengths(&[500]);
        assert_eq!(idx.unitig_at_position(0), Some((0, 0)));
        assert_eq!(idx.unitig_at_position(499), Some((0, 499)));
        assert_eq!(idx.unitig_at_position(500), None); // '$' separator
    }

    /// Regression test for the saureus `mean_id 0.624` bug: with many
    /// unitigs, a position deep in the text must map to an offset that is
    /// NOT inflated by the unitig index. Here every unitig is 10 bp, so
    /// unitig i occupies text[11*i .. 11*i + 10) and unitig i's first base
    /// is at text position 11*i.
    #[test]
    fn test_many_unitigs_offset_not_inflated() {
        let idx = CumulativeLengthIndex::from_lengths(&vec![10u64; 2000]);
        for i in [0usize, 1, 500, 1867, 1999] {
            let base = (11 * i) as u64; // start of unitig i in separator-ful text
            assert_eq!(idx.unitig_at_position(base), Some((i as u32, 0)));
            assert_eq!(idx.unitig_at_position(base + 7), Some((i as u32, 7)));
            assert_eq!(idx.unitig_at_position(base + 9), Some((i as u32, 9)));
            assert_eq!(idx.unitig_at_position(base + 10), None); // '$'
        }
    }
}
