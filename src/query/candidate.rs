/// Candidate genome filtering via information-weighted voting.
///
/// For each seed hit, look up which genomes contain the hit unitig,
/// and accumulate votes per genome weighted by information content.
/// Genome-specific unitigs (low color cardinality) contribute far more
/// than core-genome unitigs shared across many genomes.

use std::collections::HashMap;

use crate::index::color::ColorIndex;
use crate::index::fm::SeedHit;

/// A candidate genome with its accumulated vote score.
#[derive(Clone, Debug)]
pub struct Candidate {
    pub genome_id: u32,
    /// Information-weighted vote score. Higher means more specific evidence.
    pub vote_score: f64,
    /// Number of distinct unitig hits (unweighted count, for threshold filtering).
    pub vote_count: u32,
    /// Seed hits belonging to this genome.
    pub seeds: Vec<SeedHit>,
}

/// Compute information content of a unitig based on its color set size.
///
/// IC = -log2(|color_set| / total_genomes) = log2(total_genomes / |color_set|)
///
/// A genome-specific unitig (|color_set|=1) in a database of 10,000 genomes
/// scores ~13.3 bits. A core-genome unitig present in all genomes scores ~0 bits.
#[inline]
fn information_content(color_cardinality: u64, total_genomes: u64) -> f64 {
    if color_cardinality == 0 || total_genomes == 0 {
        return 0.0;
    }
    (total_genomes as f64 / color_cardinality as f64).log2()
}

/// Find candidate genomes that share unitigs with the query seeds.
///
/// Votes are weighted by information content: unitigs present in fewer genomes
/// contribute more to a genome's score than unitigs shared across many genomes.
pub fn find_candidates(
    seeds: &[SeedHit],
    color_index: &ColorIndex,
    min_votes: u32,
) -> Vec<Candidate> {
    let total_genomes = color_index.num_genomes();

    // Accumulate both weighted scores and raw counts per genome
    let mut vote_scores: HashMap<u32, f64> = HashMap::new();
    let mut vote_counts: HashMap<u32, u32> = HashMap::new();
    let mut genome_seeds: HashMap<u32, Vec<SeedHit>> = HashMap::new();

    // Track unique unitigs to avoid double-counting votes
    let mut seen_unitigs: std::collections::HashSet<u32> = std::collections::HashSet::new();

    for seed in seeds {
        if !seen_unitigs.insert(seed.unitig_id) {
            // Already counted this unitig, but still record the seed
            if let Ok(colors) = color_index.get_colors(seed.unitig_id) {
                for genome_id in colors.iter() {
                    genome_seeds
                        .entry(genome_id)
                        .or_default()
                        .push(seed.clone());
                }
            }
            continue;
        }

        // Look up which genomes contain this unitig
        if let Ok(colors) = color_index.get_colors(seed.unitig_id) {
            let cardinality = colors.len();
            let ic = information_content(cardinality, total_genomes);

            for genome_id in colors.iter() {
                *vote_scores.entry(genome_id).or_insert(0.0) += ic;
                *vote_counts.entry(genome_id).or_insert(0) += 1;
                genome_seeds
                    .entry(genome_id)
                    .or_default()
                    .push(seed.clone());
            }
        }
    }

    // Filter by minimum raw vote count (ensures sufficient evidence)
    let mut candidates: Vec<Candidate> = vote_counts
        .into_iter()
        .filter(|(_, count)| *count >= min_votes)
        .map(|(genome_id, vote_count)| Candidate {
            genome_id,
            vote_score: vote_scores.get(&genome_id).copied().unwrap_or(0.0),
            vote_count,
            seeds: genome_seeds.remove(&genome_id).unwrap_or_default(),
        })
        .collect();

    // Sort by information-weighted score descending (not raw count)
    candidates.sort_by(|a, b| b.vote_score.partial_cmp(&a.vote_score).unwrap_or(std::cmp::Ordering::Equal));

    candidates
}
