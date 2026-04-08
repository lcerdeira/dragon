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

    /// Given a position in the concatenated text, return (unitig_id, offset_within_unitig).
    pub fn unitig_at_position(&self, pos: u64) -> Option<(u32, u32)> {
        if self.cum_lengths.len() < 2 {
            return None;
        }
        // Binary search for the rightmost cumulative length <= pos
        let idx = match self.cum_lengths.binary_search(&pos) {
            Ok(i) => {
                // Exact match: pos is at the start of unitig i (or end sentinel)
                if i >= self.cum_lengths.len() - 1 {
                    return None; // past the end
                }
                i
            }
            Err(i) => {
                if i == 0 {
                    return None;
                }
                i - 1
            }
        };

        let offset = pos - self.cum_lengths[idx];
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

    #[test]
    fn test_basic_lookup() {
        // 3 unitigs of length 100, 200, 150
        let idx = CumulativeLengthIndex::from_lengths(&[100, 200, 150]);

        assert_eq!(idx.num_unitigs(), 3);
        assert_eq!(idx.total_length(), 450);

        // Position 0 -> unitig 0, offset 0
        assert_eq!(idx.unitig_at_position(0), Some((0, 0)));
        // Position 99 -> unitig 0, offset 99
        assert_eq!(idx.unitig_at_position(99), Some((0, 99)));
        // Position 100 -> unitig 1, offset 0
        assert_eq!(idx.unitig_at_position(100), Some((1, 0)));
        // Position 299 -> unitig 1, offset 199
        assert_eq!(idx.unitig_at_position(299), Some((1, 199)));
        // Position 300 -> unitig 2, offset 0
        assert_eq!(idx.unitig_at_position(300), Some((2, 0)));
        // Position 449 -> unitig 2, offset 149
        assert_eq!(idx.unitig_at_position(449), Some((2, 149)));
        // Position 450 -> past end
        assert_eq!(idx.unitig_at_position(450), None);
    }

    #[test]
    fn test_single_unitig() {
        let idx = CumulativeLengthIndex::from_lengths(&[500]);
        assert_eq!(idx.unitig_at_position(0), Some((0, 0)));
        assert_eq!(idx.unitig_at_position(499), Some((0, 499)));
        assert_eq!(idx.unitig_at_position(500), None);
    }
}
