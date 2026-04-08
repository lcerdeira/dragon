/// Seed finding via FM-index backward search.

use crate::index::fm::{DragonFmIndex, SeedHit};

/// Find seeds in a query sequence using the FM-index.
///
/// Extracts overlapping k-mers and searches each in the FM-index.
/// Supports variable-length extension for longer, more specific matches.
/// For short queries (<300 bp), uses a denser seeding strategy with
/// reduced minimum seed length for improved sensitivity.
pub fn find_seeds(
    query: &[u8],
    fm_index: &DragonFmIndex,
    min_seed_len: usize,
    max_freq: usize,
) -> Vec<SeedHit> {
    let mut seeds = Vec::new();

    // Adaptive parameters for short reads
    let (k, stride, max_extend_cap) = if query.len() < 300 {
        // Short-read mode: lower min seed length, stride 1, cap extension
        // to avoid wasting time on short queries
        let short_k = min_seed_len.min(11).max(9);
        (short_k, 1, query.len())
    } else {
        (min_seed_len, 1, 500) // Normal mode: cap extension at 500 to avoid O(n^2)
    };

    if query.len() < k {
        return seeds;
    }

    find_seeds_on_strand(query, fm_index, k, stride, max_extend_cap, max_freq, false, &mut seeds);

    // Also search reverse complement
    let rc_query = reverse_complement(query);
    find_seeds_on_strand(&rc_query, fm_index, k, stride, max_extend_cap, max_freq, true, &mut seeds);

    seeds
}

/// Find seeds on a single strand, appending to the output vec.
fn find_seeds_on_strand(
    seq: &[u8],
    fm_index: &DragonFmIndex,
    k: usize,
    stride: usize,
    max_extend_cap: usize,
    max_freq: usize,
    is_reverse: bool,
    seeds: &mut Vec<SeedHit>,
) {
    let query_len = seq.len();

    for qpos in (0..=seq.len() - k).step_by(stride) {
        let kmer = &seq[qpos..qpos + k];

        // Skip k-mers containing N
        if kmer.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
            continue;
        }

        // Variable-length search: extend as far as possible (capped)
        let max_extend = (seq.len() - qpos).min(max_extend_cap);
        let pattern = &seq[qpos..qpos + max_extend];
        let (match_len, _sa_count) = fm_index.variable_length_search(pattern);

        if match_len < k {
            continue;
        }

        // Now search for the actual match length
        let search_pattern = &seq[qpos..qpos + match_len];
        let positions = fm_index.search(search_pattern);

        // Skip if too frequent (repetitive seed)
        if positions.len() > max_freq {
            continue;
        }

        // Convert positions to seed hits
        for &pos in &positions {
            if let Some((unitig_id, offset)) = fm_index.position_to_unitig(pos) {
                let hit_query_pos = if is_reverse {
                    query_len - qpos - match_len
                } else {
                    qpos
                };

                seeds.push(SeedHit {
                    unitig_id,
                    offset,
                    query_pos: hit_query_pos,
                    match_len,
                    is_reverse,
                    sa_count: positions.len(),
                });
            }
        }
    }
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
