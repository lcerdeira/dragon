//! Direct alignment against candidate genome sequences.
//!
//! After containment-based ranking identifies the top candidate genomes, this
//! module extracts the actual genome sequence around seed anchors and runs
//! Wavefront Alignment (WFA, Marcus & Marco-Rivas 2021) to produce real
//! base-level alignments with CIGAR strings.
//!
//! Output PAF rows include:
//!   * `cg:Z:` — CIGAR with `=` (match) and `X` (mismatch) operations.
//!   * `NM:i:` — edit distance (mismatches + indel bases).
//!   * `AS:i:` — total matched bases (length-normalised score).
//!   * `de:f:` — gap-compressed divergence (1 - identity).

// The `rust_wfa` crate exposes its library under the `lib::` namespace
// because its upstream Cargo.toml sets `[lib] name = "lib"`.
use lib::alignment_lib::{Alignment, Penalties};
use lib::wavefront_alignment::wavefront_align;

use crate::index::paths::PathIndex;
use crate::index::unitig::UnitigSet;
use crate::io::paf::PafRecord;
use crate::query::containment::ContainmentHit;

/// Padding (bp) added to either side of the seed-defined reference region
/// before WFA alignment, so the aligner can find the true endpoints.
const ALIGN_PAD: usize = 100;

/// WFA gap-affine penalties (mismatch=4, gap-open=6, gap-extend=2). These
/// match minimap2's `asm5` profile and are reasonable for divergent
/// prokaryotic homologs (≤10% divergence).
fn default_penalties() -> Penalties {
    Penalties {
        mismatch_pen: 4,
        open_pen: 6,
        extd_pen: 2,
    }
}

fn rss_mb() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in s.lines() {
        if let Some(rest) = line.strip_prefix("VmRSS:") {
            let kb: u64 = rest.trim().split_whitespace().next()?.parse().ok()?;
            return Some(kb / 1024);
        }
    }
    None
}

/// Align a query directly against top candidate genomes.
pub fn direct_align_candidates(
    query: &[u8],
    query_name: &str,
    candidates: &[ContainmentHit],
    path_index: &PathIndex,
    unitigs: &UnitigSet,
    max_candidates: usize,
) -> Vec<PafRecord> {
    let mut records = Vec::new();
    let to_process = candidates.len().min(max_candidates);
    log::info!("[direct_align] entering with {} candidates ({} considered)",
               candidates.len(), to_process);
    if let Some(mb) = rss_mb() { log::info!("[direct_align] RSS={} MB", mb); }

    for (i, hit) in candidates.iter().take(max_candidates).enumerate() {
        if i % 10 == 0 {
            if let Some(mb) = rss_mb() {
                log::info!("[direct_align] candidate {}/{} RSS={} MB", i, to_process, mb);
            }
        }
        let genome_path = match path_index.get_path(hit.genome_id) {
            Some(p) => p,
            None => continue,
        };
        let genome_len = genome_path.genome_length;
        if genome_len == 0 {
            continue;
        }

        if hit.seeds.is_empty() {
            // No seed anchors — fall back to a containment-only summary row.
            // No CIGAR is emitted (we have no aligned region).
            let identity = hit.containment;
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
                tags: containment_tags(hit, query.len()),
            });
            continue;
        }

        let fwd_seeds: Vec<_> = hit.seeds.iter().filter(|s| !s.is_reverse).collect();
        let rev_seeds: Vec<_> = hit.seeds.iter().filter(|s| s.is_reverse).collect();

        if let Some(record) = align_with_seeds(
            query, query_name, &fwd_seeds, hit, &genome_path, unitigs, candidates, false,
        ) {
            records.push(record);
        }
        if let Some(record) = align_with_seeds(
            query, query_name, &rev_seeds, hit, &genome_path, unitigs, candidates, true,
        ) {
            records.push(record);
        }
    }

    records.sort_by(|a, b| {
        let a_score = extract_as_tag(a);
        let b_score = extract_as_tag(b);
        b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal)
    });
    records
}

/// Build a PAF record by running WFA on the seed-anchored region.
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

    // Map seeds to genome coordinates via the path index. Collect
    // (ref_pos, query_pos, match_len) triples so we can cluster.
    let mut anchors: Vec<(u64, usize, usize)> = Vec::new();
    for seed in seeds {
        if let Some(step) = genome_path
            .steps
            .iter()
            .find(|s| s.unitig_id == seed.unitig_id)
        {
            anchors.push((
                step.genome_offset + seed.offset as u64,
                seed.query_pos,
                seed.match_len,
            ));
        }
    }
    if anchors.is_empty() {
        return None;
    }

    // Cluster anchors around the median ref_position to avoid one rogue seed
    // (that mapped to a duplicated region elsewhere in the genome) blowing
    // the alignment window up to megabases. WFA on 59 bp × 2.8 Mb allocates
    // O(s · n) cells = many GB and OOMs the process.
    //
    // Keep only anchors within MAX_CLUSTER_SPAN bp of the median; if that
    // cluster is too tight or empty, fall back to a window centred on the
    // median. Bounds the extracted reference to (3 · query.len() + 2 · pad).
    let max_cluster_span: u64 = (query.len() as u64 * 3).max(500);
    anchors.sort_by_key(|a| a.0);
    let median_pos = anchors[anchors.len() / 2].0;
    let half_span = max_cluster_span / 2;
    let cluster_lo = median_pos.saturating_sub(half_span);
    let cluster_hi = median_pos.saturating_add(half_span);
    anchors.retain(|(p, _, _)| *p >= cluster_lo && *p <= cluster_hi);
    if anchors.is_empty() {
        return None;
    }

    let ref_min = anchors.iter().map(|a| a.0).min().unwrap();
    let ref_max = anchors
        .iter()
        .map(|(p, _, ml)| *p + *ml as u64)
        .max()
        .unwrap_or(ref_min);
    let query_min = anchors.iter().map(|a| a.1).min().unwrap();
    let query_max = anchors
        .iter()
        .map(|(_, q, ml)| *q + *ml)
        .max()
        .unwrap_or(query_min);

    if query_max <= query_min || ref_max <= ref_min {
        return None;
    }

    // Pad the reference window. After clustering, this is bounded by
    // max_cluster_span + 2 · ALIGN_PAD ≈ a few hundred bp for short queries,
    // a few kb for longer ones — never megabases.
    let extract_start = ref_min.saturating_sub(ALIGN_PAD as u64);
    let extract_end = (ref_max + ALIGN_PAD as u64).min(genome_path.genome_length);
    let ref_seq = PathIndex::extract_sequence_static(
        genome_path,
        extract_start,
        extract_end,
        unitigs,
    );
    if ref_seq.is_empty() {
        return None;
    }

    // Slice + orient the query to match the seed strand.
    let query_slice_end = query_max.min(query.len());
    if query_slice_end <= query_min {
        return None;
    }
    let mut query_oriented: Vec<u8> = query[query_min..query_slice_end].to_vec();
    if is_reverse {
        query_oriented = reverse_complement(&query_oriented);
    }

    // WFA constraint: query.len() <= text.len(). If padding wasn't enough,
    // grow the reference window symmetrically.
    let q_len = query_oriented.len();
    let mut t = ref_seq;
    if q_len > t.len() {
        let need = q_len - t.len() + 16;
        let extra_left = (need / 2) as u64;
        let extra_right = (need - extra_left as usize) as u64;
        let new_start = extract_start.saturating_sub(extra_left);
        let new_end = (extract_end + extra_right).min(genome_path.genome_length);
        t = PathIndex::extract_sequence_static(genome_path, new_start, new_end, unitigs);
        if q_len > t.len() {
            // Cannot align — query truly longer than the genome region.
            return None;
        }
    }

    // Run WFA.
    let q_str = match std::str::from_utf8(&query_oriented) {
        Ok(s) => s,
        Err(_) => return None,
    };
    let t_str = match std::str::from_utf8(&t) {
        Ok(s) => s,
        Err(_) => return None,
    };
    let alignment = match wavefront_align(q_str, t_str, &default_penalties()) {
        Ok(a) => a,
        Err(_) => return None,
    };

    // Translate the WFA aligned strings into CIGAR + counts.
    let stats = AlignmentStats::from_wfa(&alignment);

    // Compute genome-coordinate target span from the leading/trailing
    // text-side gaps so target_start/target_end describe only the aligned
    // portion (not the padded extraction window).
    let text_lead_in_ref = leading_gap_consumed(&alignment.query_aligned, &alignment.text_aligned);

    // text-side leading skip = bases of t consumed before the alignment proper begins
    let aligned_target_len = stats.target_consumed;
    let target_start_in_window = text_lead_in_ref;
    let target_end_in_window = target_start_in_window + aligned_target_len;

    let extract_window_start = if q_len > t.len() {
        // Re-derived window already used above; recompute conservatively.
        let pad_lhs = ALIGN_PAD as u64;
        ref_min.saturating_sub(pad_lhs)
    } else {
        extract_start
    };
    let target_start = extract_window_start as usize + target_start_in_window;
    let target_end = extract_window_start as usize + target_end_in_window;

    // Query span on the original (un-RC'd) coordinates.
    let q_lead_gap = leading_gap_consumed(&alignment.text_aligned, &alignment.query_aligned);
    let q_aligned_len = stats.query_consumed;
    let (final_q_start, final_q_end) = if is_reverse {
        // We aligned RC(query[query_min..query_slice_end]) against t.
        // Map back: a RC slice's [a, b) corresponds to original [L-b, L-a)
        // where L = q_str.len().
        let l = q_str.len();
        let a = q_lead_gap;
        let b = a + q_aligned_len;
        let rc_start = query_min + (l - b);
        let rc_end = query_min + (l - a);
        (rc_start, rc_end)
    } else {
        (query_min + q_lead_gap, query_min + q_lead_gap + q_aligned_len)
    };

    let alignment_len = stats.alignment_len;
    let num_matches = stats.matches;

    let identity = if alignment_len == 0 {
        0.0
    } else {
        num_matches as f64 / alignment_len as f64
    };

    let mut tags = containment_tags(hit, query.len());
    tags.push(format!("NM:i:{}", stats.edits));
    tags.push(format!("de:f:{:.4}", 1.0 - identity));
    tags.push(format!("cg:Z:{}", stats.cigar));

    Some(PafRecord {
        query_name: query_name.to_string(),
        query_len: query.len(),
        query_start: final_q_start,
        query_end: final_q_end,
        strand: if is_reverse { '-' } else { '+' },
        target_name: genome_path.genome_name.clone(),
        target_len: genome_path.genome_length as usize,
        target_start,
        target_end,
        num_matches,
        alignment_len,
        mapq: estimate_mapq(hit, all_candidates),
        tags,
    })
}

/// Per-alignment statistics derived from the aligned strings emitted by WFA.
#[derive(Debug, Default, Clone)]
struct AlignmentStats {
    cigar: String,
    matches: usize,
    mismatches: usize,
    insertions: usize, // gaps in target → bases added from query
    deletions: usize,  // gaps in query → bases consumed from target only
    /// Bases consumed from the query (excluding leading/trailing gaps in query_aligned).
    query_consumed: usize,
    /// Bases consumed from the target (excluding leading/trailing gaps in text_aligned).
    target_consumed: usize,
    /// Total CIGAR length.
    alignment_len: usize,
    /// NM-equivalent edit distance.
    edits: usize,
}

impl AlignmentStats {
    fn from_wfa(a: &Alignment) -> Self {
        let q: Vec<u8> = a.query_aligned.bytes().collect();
        let t: Vec<u8> = a.text_aligned.bytes().collect();
        debug_assert_eq!(q.len(), t.len(), "WFA aligned strings must be equal length");

        // Trim leading and trailing positions where BOTH strings are gaps —
        // shouldn't happen with WFA but be defensive.
        // We separately handle leading/trailing single-sided gaps as soft
        // clips at the alignment endpoints.

        // Find first/last index where at least one side is a non-gap base —
        // this is the body of the alignment. Any single-sided gap region at
        // the ends is treated as part of the alignment (insertion/deletion).
        let n = q.len();
        let mut s = Self::default();
        if n == 0 {
            return s;
        }

        let mut cigar = String::with_capacity(n / 4);
        let mut run_op: Option<u8> = None;
        let mut run_len: usize = 0;

        let push_run = |cigar: &mut String, op: u8, len: usize| {
            if len > 0 {
                cigar.push_str(&len.to_string());
                cigar.push(op as char);
            }
        };

        for i in 0..n {
            let qb = q[i];
            let tb = t[i];
            let op: u8;
            if qb == b'-' && tb != b'-' {
                op = b'D';
                s.deletions += 1;
                s.target_consumed += 1;
            } else if tb == b'-' && qb != b'-' {
                op = b'I';
                s.insertions += 1;
                s.query_consumed += 1;
            } else if qb.eq_ignore_ascii_case(&tb) {
                op = b'=';
                s.matches += 1;
                s.query_consumed += 1;
                s.target_consumed += 1;
            } else {
                op = b'X';
                s.mismatches += 1;
                s.query_consumed += 1;
                s.target_consumed += 1;
            }

            if Some(op) == run_op {
                run_len += 1;
            } else {
                if let Some(prev) = run_op {
                    push_run(&mut cigar, prev, run_len);
                }
                run_op = Some(op);
                run_len = 1;
            }
        }
        if let Some(op) = run_op {
            push_run(&mut cigar, op, run_len);
        }

        s.cigar = cigar;
        s.alignment_len = n;
        s.edits = s.mismatches + s.insertions + s.deletions;
        s
    }
}

/// Number of bases consumed from `consumed` while `gapped` is leading-gap.
/// E.g. consumed = "ACGT...", gapped = "----..." → returns 4.
fn leading_gap_consumed(gapped: &str, consumed: &str) -> usize {
    gapped
        .bytes()
        .zip(consumed.bytes())
        .take_while(|&(g, _)| g == b'-')
        .filter(|&(_, c)| c != b'-')
        .count()
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len());
    for &b in s.iter().rev() {
        out.push(match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        });
    }
    out
}

fn estimate_mapq(hit: &ContainmentHit, all_candidates: &[ContainmentHit]) -> u8 {
    if all_candidates.len() <= 1 {
        return 60;
    }
    let best = all_candidates[0].containment;
    let second = all_candidates.get(1).map(|c| c.containment).unwrap_or(0.0);
    let gap = if best > 0.0 { (best - second) / best } else { 0.0 };
    let containment_component = (hit.containment * 30.0) as u8;
    let gap_component = (gap * 30.0) as u8;
    containment_component.saturating_add(gap_component).min(60)
}

fn containment_tags(hit: &ContainmentHit, query_len: usize) -> Vec<String> {
    let as_score = (hit.containment * query_len as f64) as usize;
    vec![
        format!("AS:i:{}", as_score),
        format!("ct:f:{:.4}", hit.containment),
        format!("ic:f:{:.1}", hit.info_score),
    ]
}

fn extract_as_tag(record: &PafRecord) -> f64 {
    record
        .tags
        .iter()
        .find_map(|t| t.strip_prefix("AS:i:").and_then(|v| v.parse::<f64>().ok()))
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_stats_perfect_match() {
        // Both strings identical → 4 matches, no edits.
        let a = Alignment {
            score: 0,
            query_aligned: "ACGT".to_string(),
            text_aligned: "ACGT".to_string(),
        };
        let s = AlignmentStats::from_wfa(&a);
        assert_eq!(s.matches, 4);
        assert_eq!(s.mismatches, 0);
        assert_eq!(s.insertions, 0);
        assert_eq!(s.deletions, 0);
        assert_eq!(s.cigar, "4=");
        assert_eq!(s.edits, 0);
        assert_eq!(s.alignment_len, 4);
    }

    #[test]
    fn alignment_stats_mismatch_in_middle() {
        // Single substitution.
        let a = Alignment {
            score: 4,
            query_aligned: "ACGT".to_string(),
            text_aligned: "ACAT".to_string(),
        };
        let s = AlignmentStats::from_wfa(&a);
        assert_eq!(s.matches, 3);
        assert_eq!(s.mismatches, 1);
        assert_eq!(s.cigar, "2=1X1=");
        assert_eq!(s.edits, 1);
    }

    #[test]
    fn alignment_stats_with_indels() {
        // Insertion in query (gap in text) + deletion (gap in query).
        let a = Alignment {
            score: 0,
            query_aligned: "AC-GT".to_string(),
            text_aligned: "ACAGT".to_string(),
        };
        let s = AlignmentStats::from_wfa(&a);
        assert_eq!(s.matches, 4);
        assert_eq!(s.deletions, 1);
        assert_eq!(s.insertions, 0);
        assert_eq!(s.cigar, "2=1D2=");
        assert_eq!(s.edits, 1);
    }

    #[test]
    fn wfa_end_to_end_smoke() {
        // Drives the rust_wfa crate end-to-end, then inspects the CIGAR.
        let q = "ACGTACGTACGT";
        let t = "ACGTACGTACGTAA"; // exact match in prefix, two extra bases
        let aln = wavefront_align(q, t, &default_penalties()).unwrap();
        let s = AlignmentStats::from_wfa(&aln);
        // The first 12 bases match exactly. WFA aligns to score=0 with all
        // matches. (Trailing bases of t are not consumed because alignment
        // can stop where the score is best.)
        assert!(s.matches >= 12);
        assert_eq!(s.mismatches, 0);
        assert!(s.cigar.contains('='));
    }

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"AGCN"), b"NGCT");
    }
}
