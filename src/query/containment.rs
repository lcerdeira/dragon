/// Containment-based genome ranking.
///
/// Computes k-mer containment between a query and each candidate genome:
///   containment(Q, G) = |kmers(Q) ∩ kmers(G)| / |kmers(Q)|
///
/// This is mathematically related to ANI and provides near-perfect sensitivity
/// for queries with ANI > 90% to the target genome. It bypasses the chaining
/// pipeline entirely, using the color index to check genome membership directly.
///
/// For each query k-mer that matches a unitig in the FM-index, we check which
/// genomes contain that unitig via the color index. A genome's containment score
/// is the fraction of query k-mers found in it.

use std::collections::HashMap;

use roaring::RoaringBitmap;

use crate::index::color::ColorIndex;
use crate::index::fm::{DragonFmIndex, SeedHit};

/// A genome ranked by k-mer containment.
#[derive(Clone, Debug)]
pub struct ContainmentHit {
    pub genome_id: u32,
    /// Fraction of query k-mers found in this genome (0.0 to 1.0).
    pub containment: f64,
    /// Number of query k-mers found in this genome.
    pub shared_kmers: usize,
    /// Total query k-mers searched.
    pub total_query_kmers: usize,
    /// Information-weighted score (sum of IC for matched unitigs).
    pub info_score: f64,
    /// Seeds that map to this genome (for optional chaining/alignment).
    pub seeds: Vec<SeedHit>,
}

/// Compute containment-based ranking for a query against all genomes.
///
/// Steps:
/// 1. Extract all k-mers from the query at stride 1
/// 2. Search each in the FM-index to find matching unitigs
/// 3. For each matching unitig, look up genome membership via color index
/// 4. Accumulate per-genome k-mer counts
/// 5. Rank by containment = shared_kmers / total_query_kmers
pub fn containment_rank(
    query: &[u8],
    fm_index: &DragonFmIndex,
    color_index: &ColorIndex,
    kmer_size: usize,
    max_seed_freq: usize,
) -> Vec<ContainmentHit> {
    if query.len() < kmer_size {
        return Vec::new();
    }

    let total_genomes = color_index.num_genomes();
    let total_query_kmers = query.len() - kmer_size + 1;

    // Per-genome accumulators
    let mut genome_shared: HashMap<u32, usize> = HashMap::new();
    let mut genome_info: HashMap<u32, f64> = HashMap::new();
    let mut genome_seeds: HashMap<u32, Vec<SeedHit>> = HashMap::new();

    // Track which query positions already contributed (avoid double-counting overlapping seeds)
    let mut counted_positions = vec![false; query.len()];

    // Search both strands
    for (seq, is_reverse) in [(query.to_vec(), false), (reverse_complement(query), true)] {
        if seq.len() < kmer_size {
            continue;
        }

        for qpos in 0..=seq.len() - kmer_size {
            let kmer = &seq[qpos..qpos + kmer_size];

            // Skip k-mers with N
            if kmer.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
                continue;
            }

            // Search in FM-index
            let positions = fm_index.search(kmer);
            if positions.is_empty() || positions.len() > max_seed_freq {
                continue;
            }

            // Map each position to a unitig, collect unique unitigs
            let mut hit_unitigs = RoaringBitmap::new();
            let mut unitig_seeds: Vec<SeedHit> = Vec::new();

            for &pos in &positions {
                if let Some((unitig_id, offset)) = fm_index.position_to_unitig(pos) {
                    if hit_unitigs.insert(unitig_id) {
                        let hit_qpos = if is_reverse {
                            seq.len() - qpos - kmer_size
                        } else {
                            qpos
                        };
                        unitig_seeds.push(SeedHit {
                            unitig_id,
                            offset,
                            query_pos: hit_qpos,
                            match_len: kmer_size,
                            is_reverse,
                            sa_count: positions.len(),
                        });
                    }
                }
            }

            // For each unique unitig hit, credit the genomes that contain it
            for unitig_id in hit_unitigs.iter() {
                if let Ok(colors) = color_index.get_colors(unitig_id) {
                    let cardinality = colors.len();
                    let ic = if cardinality > 0 && total_genomes > 0 {
                        (total_genomes as f64 / cardinality as f64).log2()
                    } else {
                        0.0
                    };

                    // Only count this query position once per genome
                    let actual_qpos = if is_reverse {
                        seq.len() - qpos - kmer_size
                    } else {
                        qpos
                    };

                    for genome_id in colors.iter() {
                        if !counted_positions[actual_qpos] {
                            *genome_shared.entry(genome_id).or_insert(0) += 1;
                        }
                        *genome_info.entry(genome_id).or_insert(0.0) += ic;

                        // Store seeds for later alignment
                        for seed in &unitig_seeds {
                            if color_index.get_colors(seed.unitig_id)
                                .map(|c| c.contains(genome_id))
                                .unwrap_or(false)
                            {
                                genome_seeds.entry(genome_id).or_default().push(seed.clone());
                            }
                        }
                    }
                }
            }

            // Mark this position as counted (for first strand that matches)
            if !hit_unitigs.is_empty() && !is_reverse {
                counted_positions[qpos] = true;
            }
        }
    }

    // Build ranked results
    let mut hits: Vec<ContainmentHit> = genome_shared
        .into_iter()
        .map(|(genome_id, shared)| ContainmentHit {
            genome_id,
            containment: shared as f64 / total_query_kmers.max(1) as f64,
            shared_kmers: shared,
            total_query_kmers,
            info_score: genome_info.get(&genome_id).copied().unwrap_or(0.0),
            seeds: genome_seeds.remove(&genome_id).unwrap_or_default(),
        })
        .collect();

    // Sort by containment descending, break ties by info_score
    hits.sort_by(|a, b| {
        b.containment.partial_cmp(&a.containment)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.info_score.partial_cmp(&a.info_score)
                .unwrap_or(std::cmp::Ordering::Equal))
    });

    hits
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACGTACGT"), b"ACGTACGT");
    }
}
