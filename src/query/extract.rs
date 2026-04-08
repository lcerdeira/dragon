/// Reference sequence extraction by walking genome paths through the de Bruijn graph.

use crate::index::paths::PathIndex;
use crate::index::unitig::UnitigSet;
use crate::query::chain::Chain;

/// Region of a reference genome to align against.
pub struct ExtractedRegion {
    pub genome_id: u32,
    pub start: u64,
    pub end: u64,
    pub sequence: Vec<u8>,
}

/// Extract reference regions for a list of chains.
pub fn extract_regions(
    chains: &[Chain],
    path_index: &PathIndex,
    unitigs: &UnitigSet,
    flank: u64,
) -> Vec<ExtractedRegion> {
    let mut regions = Vec::new();

    for chain in chains {
        if chain.anchors.is_empty() {
            continue;
        }

        let first = chain.anchors.first().unwrap();
        let last = chain.anchors.last().unwrap();

        let start = first.ref_start.saturating_sub(flank);
        let end = last.ref_end + flank;

        // Clamp to genome length
        let genome_len = path_index
            .get_path(chain.genome_id)
            .map(|p| p.genome_length)
            .unwrap_or(0);
        let end = end.min(genome_len);

        let sequence = path_index.extract_sequence(chain.genome_id, start, end, unitigs);

        regions.push(ExtractedRegion {
            genome_id: chain.genome_id,
            start,
            end,
            sequence,
        });
    }

    regions
}
