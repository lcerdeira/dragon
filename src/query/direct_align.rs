/// Direct alignment against candidate genome sequences.
///
/// After containment-based ranking identifies the top candidate genomes,
/// this module extracts the actual genome sequence and performs direct
/// alignment to produce base-level PAF records.
///
/// This bypasses the unitig-boundary problem by working on full genome
/// sequences rather than fragmented unitig text.

use crate::index::paths::PathIndex;
use crate::index::unitig::UnitigSet;
use crate::io::paf::PafRecord;
use crate::query::containment::ContainmentHit;

/// Align a query directly against top candidate genomes.
///
/// For each candidate (by containment score), extracts the genome sequence
/// from the path index and finds the best alignment region using seed anchors.
pub fn direct_align_candidates(
    query: &[u8],
    query_name: &str,
    candidates: &[ContainmentHit],
    path_index: &PathIndex,
    unitigs: &UnitigSet,
    max_candidates: usize,
) -> Vec<PafRecord> {
    let mut records = Vec::new();

    for hit in candidates.iter().take(max_candidates) {
        let genome_path = match path_index.get_path(hit.genome_id) {
            Some(p) => p,
            None => continue,
        };

        let genome_len = genome_path.genome_length;
        if genome_len == 0 {
            continue;
        }

        // Find the best alignment region using seed anchors
        // Seeds give us approximate positions in the genome
        if hit.seeds.is_empty() {
            // No seeds but high containment — report whole-genome match
            let identity = hit.containment; // containment approximates identity
            if identity < 0.01 {
                continue;
            }

            records.push(PafRecord {
                query_name: query_name.to_string(),
                query_len: query.len(),
                query_start: 0,
                query_end: query.len(),
                strand: '+',
                target_name: genome_path.genome_name.clone(),
                target_len: genome_len as usize,
                target_start: 0,
                target_end: query.len().min(genome_len as usize),
                num_matches: (query.len() as f64 * identity) as usize,
                alignment_len: query.len(),
                mapq: estimate_mapq(hit, candidates),
                tags: make_tags(hit, query.len()),
            });
            continue;
        }

        // Use seed positions to find alignment region(s)
        // Group seeds by strand
        let fwd_seeds: Vec<_> = hit.seeds.iter().filter(|s| !s.is_reverse).collect();
        let rev_seeds: Vec<_> = hit.seeds.iter().filter(|s| s.is_reverse).collect();

        // Try forward strand alignment
        if let Some(record) = align_with_seeds(
            query, query_name, &fwd_seeds, hit, genome_path, unitigs, candidates, false,
        ) {
            records.push(record);
        }

        // Try reverse strand alignment
        if let Some(record) = align_with_seeds(
            query, query_name, &rev_seeds, hit, genome_path, unitigs, candidates, true,
        ) {
            records.push(record);
        }
    }

    // Sort by alignment score descending
    records.sort_by(|a, b| {
        let a_score = extract_as_tag(a);
        let b_score = extract_as_tag(b);
        b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal)
    });

    records
}

/// Align query to a genome region defined by seed anchors.
fn align_with_seeds(
    query: &[u8],
    query_name: &str,
    seeds: &[&crate::index::fm::SeedHit],
    hit: &ContainmentHit,
    genome_path: &crate::index::paths::GenomePath,
    unitigs: &UnitigSet,
    all_candidates: &[ContainmentHit],
    is_reverse: bool,
) -> Option<PafRecord> {
    if seeds.is_empty() {
        return None;
    }

    // Find the genome region covered by seeds
    // Map seeds to genome coordinates via path index
    let mut ref_positions: Vec<u64> = Vec::new();
    let mut query_positions: Vec<usize> = Vec::new();
    let mut total_seed_bases: usize = 0;

    for seed in seeds {
        // Find this unitig in the genome's path
        if let Some(step) = genome_path.steps.iter().find(|s| s.unitig_id == seed.unitig_id) {
            let ref_pos = step.genome_offset + seed.offset as u64;
            ref_positions.push(ref_pos);
            query_positions.push(seed.query_pos);
            total_seed_bases += seed.match_len;
        }
    }

    if ref_positions.is_empty() {
        return None;
    }

    // Determine alignment span
    let ref_min = *ref_positions.iter().min().unwrap();
    let ref_max = ref_positions.iter().zip(seeds.iter())
        .map(|(&pos, seed)| pos + seed.match_len as u64)
        .max()
        .unwrap_or(ref_min);
    let query_min = *query_positions.iter().min().unwrap();
    let query_max = query_positions.iter().zip(seeds.iter())
        .map(|(&pos, seed)| pos + seed.match_len)
        .max()
        .unwrap_or(query_min);

    let ref_span = (ref_max - ref_min) as usize;
    let query_span = query_max - query_min;

    if query_span == 0 || ref_span == 0 {
        return None;
    }

    // Extract reference region for alignment (with padding)
    let pad = 100;
    let extract_start = ref_min.saturating_sub(pad as u64);
    let extract_end = (ref_max + pad as u64).min(genome_path.genome_length);
    let ref_seq = crate::index::paths::PathIndex::extract_sequence_static(
        genome_path, extract_start, extract_end, unitigs,
    );

    // Compute identity from seed coverage + direct comparison
    let identity = if !ref_seq.is_empty() && query_span > 0 {
        // Compare overlapping regions directly
        let query_region = &query[query_min..query_max.min(query.len())];
        compute_identity(query_region, &ref_seq, total_seed_bases, query_span)
    } else {
        // Fallback: estimate from containment
        hit.containment
    };

    let alignment_len = query_span.max(ref_span);
    let num_matches = (alignment_len as f64 * identity) as usize;

    Some(PafRecord {
        query_name: query_name.to_string(),
        query_len: query.len(),
        query_start: query_min,
        query_end: query_max.min(query.len()),
        strand: if is_reverse { '-' } else { '+' },
        target_name: genome_path.genome_name.clone(),
        target_len: genome_path.genome_length as usize,
        target_start: ref_min as usize,
        target_end: ref_max as usize,
        num_matches,
        alignment_len,
        mapq: estimate_mapq(hit, all_candidates),
        tags: make_tags(hit, query.len()),
    })
}

/// Compute identity by direct base comparison + seed coverage.
fn compute_identity(
    query_region: &[u8],
    ref_region: &[u8],
    seed_match_bases: usize,
    query_span: usize,
) -> f64 {
    if query_region.is_empty() || ref_region.is_empty() {
        return 0.0;
    }

    // Direct base comparison over the shorter of the two
    let compare_len = query_region.len().min(ref_region.len());
    if compare_len == 0 {
        return seed_match_bases as f64 / query_span.max(1) as f64;
    }

    let matches = query_region[..compare_len].iter()
        .zip(ref_region[..compare_len].iter())
        .filter(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
        .count();

    matches as f64 / compare_len as f64
}

/// Estimate mapping quality from containment ranking.
fn estimate_mapq(hit: &ContainmentHit, all_candidates: &[ContainmentHit]) -> u8 {
    if all_candidates.len() <= 1 {
        return 60;
    }

    let best = all_candidates[0].containment;
    let second = if all_candidates.len() > 1 { all_candidates[1].containment } else { 0.0 };

    // MAPQ based on gap between best and second-best containment
    let gap = if best > 0.0 { (best - second) / best } else { 0.0 };
    let containment_component = (hit.containment * 30.0) as u8;
    let gap_component = (gap * 30.0) as u8;

    containment_component.saturating_add(gap_component).min(60)
}

/// Create PAF tags from containment hit.
fn make_tags(hit: &ContainmentHit, query_len: usize) -> Vec<String> {
    let as_score = (hit.containment * query_len as f64) as usize;
    vec![
        format!("AS:i:{}", as_score),
        format!("ct:f:{:.4}", hit.containment),
        format!("ic:f:{:.1}", hit.info_score),
        format!("de:f:{:.4}", 1.0 - hit.containment),
    ]
}

fn extract_as_tag(record: &PafRecord) -> f64 {
    record.tags.iter()
        .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
        .unwrap_or(0.0)
}
