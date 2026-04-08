/// Candidate genome filtering via color-based voting.
///
/// For each seed hit, look up which genomes contain the hit unitig,
/// and accumulate votes per genome. Genomes exceeding a threshold are candidates.

use std::collections::HashMap;

use crate::index::color::ColorIndex;
use crate::index::fm::SeedHit;

/// A candidate genome with its vote count.
#[derive(Clone, Debug)]
pub struct Candidate {
    pub genome_id: u32,
    pub vote_count: u32,
    /// Seed hits belonging to this genome.
    pub seeds: Vec<SeedHit>,
}

/// Find candidate genomes that share unitigs with the query seeds.
pub fn find_candidates(
    seeds: &[SeedHit],
    color_index: &ColorIndex,
    min_votes: u32,
) -> Vec<Candidate> {
    // Accumulate votes per genome
    let mut votes: HashMap<u32, u32> = HashMap::new();
    let mut genome_seeds: HashMap<u32, Vec<SeedHit>> = HashMap::new();

    // Track unique unitigs to avoid double-counting
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
            for genome_id in colors.iter() {
                *votes.entry(genome_id).or_insert(0) += 1;
                genome_seeds
                    .entry(genome_id)
                    .or_default()
                    .push(seed.clone());
            }
        }
    }

    // Filter by minimum vote threshold
    let mut candidates: Vec<Candidate> = votes
        .into_iter()
        .filter(|(_, count)| *count >= min_votes)
        .map(|(genome_id, vote_count)| Candidate {
            genome_id,
            vote_count,
            seeds: genome_seeds.remove(&genome_id).unwrap_or_default(),
        })
        .collect();

    // Sort by vote count descending
    candidates.sort_by(|a, b| b.vote_count.cmp(&a.vote_count));

    candidates
}
