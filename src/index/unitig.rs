/// Unitig representation and 2-bit encoding.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::io::fasta;
use crate::util::dna::PackedSequence;

/// A single unitig from the de Bruijn graph.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Unitig {
    pub id: u32,
    pub sequence: PackedSequence,
}

/// Collection of all unitigs in the index.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UnitigSet {
    pub unitigs: Vec<Unitig>,
    /// Concatenated text for FM-index construction (unitig sequences joined by separators).
    pub concatenated: Vec<u8>,
    /// Lengths of each unitig (for CumulativeLengthIndex).
    pub lengths: Vec<u64>,
}

impl UnitigSet {
    pub fn num_unitigs(&self) -> usize {
        self.unitigs.len()
    }

    pub fn total_bases(&self) -> u64 {
        self.lengths.iter().sum()
    }

    /// Get the ASCII sequence of a unitig.
    pub fn get_sequence(&self, unitig_id: u32) -> Vec<u8> {
        self.unitigs[unitig_id as usize].sequence.to_bytes()
    }

    /// Get a subsequence of a unitig.
    pub fn get_subsequence(&self, unitig_id: u32, start: usize, end: usize) -> Vec<u8> {
        self.unitigs[unitig_id as usize]
            .sequence
            .subsequence(start, end)
    }

    /// Reconstruct a UnitigSet from FM-index components.
    ///
    /// The FM-index stores concatenated unitig text with '$' separators.
    /// Combined with the CumulativeLengthIndex boundaries, we can reconstruct
    /// each unitig without serializing a separate file.
    pub fn from_fm_text(text: &[u8], lengths: &[u64]) -> Self {
        let mut unitigs = Vec::with_capacity(lengths.len());
        let mut concatenated = Vec::new();
        let mut offset = 0usize;

        for (id, &len) in lengths.iter().enumerate() {
            let end = offset + len as usize;
            let seq_bytes = &text[offset..end.min(text.len())];
            let packed = PackedSequence::from_bytes(seq_bytes);

            concatenated.extend_from_slice(seq_bytes);
            concatenated.push(b'$');

            unitigs.push(Unitig {
                id: id as u32,
                sequence: packed,
            });

            // Skip past the '$' separator
            offset = end + 1;
        }

        Self {
            unitigs,
            concatenated,
            lengths: lengths.to_vec(),
        }
    }
}

/// Parse unitigs from a FASTA file (GGCAT output) and encode in 2-bit format.
pub fn parse_and_encode_unitigs(unitig_file: &Path) -> Result<UnitigSet> {
    let sequences = fasta::read_sequences(unitig_file)?;

    let mut unitigs = Vec::with_capacity(sequences.len());
    let mut lengths = Vec::with_capacity(sequences.len());
    // Build concatenated text with '$' separators for FM-index
    let mut concatenated = Vec::new();

    for (id, seq) in sequences.into_iter().enumerate() {
        let packed = PackedSequence::from_bytes(&seq.seq);
        lengths.push(seq.seq.len() as u64);

        // Append to concatenated text
        concatenated.extend_from_slice(&seq.seq);
        concatenated.push(b'$'); // separator

        unitigs.push(Unitig {
            id: id as u32,
            sequence: packed,
        });
    }

    log::info!(
        "Parsed {} unitigs, {} total bases",
        unitigs.len(),
        lengths.iter().sum::<u64>()
    );

    Ok(UnitigSet {
        unitigs,
        concatenated,
        lengths,
    })
}
