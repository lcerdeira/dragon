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
/// Kept small because rust_wfa does GLOBAL alignment — every byte of
/// reference padding past the actual match becomes a gap operation.
const ALIGN_PAD: usize = 20;

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

        // Partition seeds by EFFECTIVE strand of the query→genome match,
        // not by seed.is_reverse alone. The effective strand is the XOR
        // of seed.is_reverse (was the QUERY RC'd when searching the
        // FM-index?) and step.is_reverse (is the GENOME reading the
        // unitig in RC?). Only seeds whose XOR is false anchor a
        // forward-strand alignment; the others anchor a reverse-strand
        // alignment and need the query RC'd before WFA.
        let mut fwd_seeds: Vec<&crate::index::fm::SeedHit> = Vec::new();
        let mut rev_seeds: Vec<&crate::index::fm::SeedHit> = Vec::new();
        for s in &hit.seeds {
            let step_rev = genome_path
                .steps
                .iter()
                .find(|st| st.unitig_id == s.unitig_id)
                .map(|st| st.is_reverse)
                .unwrap_or(false);
            if s.is_reverse == step_rev {
                fwd_seeds.push(s);
            } else {
                rev_seeds.push(s);
            }
        }

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

    // Map seeds to genome coordinates via the path index. The genome
    // position of the seed match depends ONLY on step.is_reverse:
    //   * step.is_reverse=false: genome reads the unitig forward, so
    //     a seed at unitig offset O lands at genome step.genome_offset + O.
    //   * step.is_reverse=true: genome reads RC(unitig), so a seed at
    //     unitig offset O lands at genome step.genome_offset +
    //     unitig_len - (O + match_len).
    // The query orientation (seed.is_reverse) governs the alignment
    // strand — but that's already taken care of by the caller, who has
    // partitioned seeds into fwd_seeds / rev_seeds by EFFECTIVE strand
    // (seed.is_reverse XOR step.is_reverse). At this point all anchors
    // we accept share the same effective strand; we just need the
    // correct genome anchor coordinates.
    let mut anchors: Vec<(u64, usize, usize)> = Vec::new();
    for seed in seeds {
        if let Some(step) = genome_path
            .steps
            .iter()
            .find(|s| s.unitig_id == seed.unitig_id)
        {
            let unitig_len = unitigs
                .unitigs
                .get(seed.unitig_id as usize)
                .map(|u| u.sequence.len as u64)
                .unwrap_or(0);
            let pos_in_genome = if step.is_reverse {
                let off = seed.offset as u64 + seed.match_len as u64;
                step.genome_offset + unitig_len.saturating_sub(off)
            } else {
                step.genome_offset + seed.offset as u64
            };
            anchors.push((pos_in_genome, seed.query_pos, seed.match_len));
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

    // Project the reference window to span the FULL query.
    //
    // For FORWARD alignment, anchors map original query[X..X+M] → genome
    // [pos_in_genome..pos_in_genome+M]. So query[0] → ref_min - query_min,
    // and the window is [ref_min - query_min .. ref_min - query_min + L].
    //
    // For REVERSE alignment we align RC(query) (length L) against genome
    // forward. An anchor with original query_pos = X and match_len = M
    // corresponds to RC(query) positions [L-X-M .. L-X]. So RC(query)[0]
    // lands at pos_in_genome - (L - X - M), and the smallest such value
    // across seeds is ref_min - (L - max(X + M)) = ref_min - (L - query_max).
    // Using `query_min` here (the original-query-coord min) shifts the
    // window left by an arbitrary amount and was the source of the
    // saureus 0.74 identity ceiling on reverse-strand hits.
    let l = query.len() as u64;
    let est_ref_start_for_query = if is_reverse {
        ref_min.saturating_sub(l.saturating_sub(query_max as u64))
    } else {
        ref_min.saturating_sub(query_min as u64)
    };
    let est_ref_end_for_query = est_ref_start_for_query + l;
    let extract_start = est_ref_start_for_query.saturating_sub(ALIGN_PAD as u64);
    let extract_end = (est_ref_end_for_query + ALIGN_PAD as u64)
        .min(genome_path.genome_length);

    log::debug!(
        "[align_with_seeds] genome={} q_len={} anchors={} ref_min={} ref_max={} extract=[{}, {}) win={}",
        hit.genome_id,
        query.len(),
        anchors.len(),
        ref_min,
        ref_max,
        extract_start,
        extract_end,
        extract_end.saturating_sub(extract_start),
    );
    let ref_seq = PathIndex::extract_sequence_static(
        genome_path,
        extract_start,
        extract_end,
        unitigs,
    );
    if ref_seq.is_empty() {
        return None;
    }
    log::debug!(
        "[align_with_seeds] genome={} ref_seq.len()={}",
        hit.genome_id, ref_seq.len()
    );

    // Hard cap on the alignment text size: never feed WFA more than
    // 4 · query.len() of reference. If extract_sequence_static produced
    // more than that (e.g. because path steps overlap or expanded
    // unexpectedly), trim the reference to a window centred on the
    // anchor cluster's median in target-window coordinates.
    let max_ref_len = (query.len() * 4).max(400);
    let ref_seq = if ref_seq.len() > max_ref_len {
        let median_in_window = (ref_min - extract_start) as usize;
        let half = max_ref_len / 2;
        let lo = median_in_window.saturating_sub(half);
        let hi = (lo + max_ref_len).min(ref_seq.len());
        log::debug!(
            "[align_with_seeds] trimmed ref_seq from {} -> {} bp (window {}..{})",
            ref_seq.len(),
            hi - lo,
            lo,
            hi
        );
        ref_seq[lo..hi].to_vec()
    } else {
        ref_seq
    };

    // Align the FULL query (not just the seed-covered portion). Slicing to
    // query[query_min..query_max] artificially capped query coverage at
    // (query_max - query_min) / query.len() — typically ~50 % for short
    // queries — even when the rest of the query was present in the
    // genome. WFA does global alignment, so giving it the whole query
    // against a tightly-sized reference window produces meaningful
    // identity numbers.
    let mut query_oriented: Vec<u8> = query.to_vec();
    let query_min_in_align: usize = 0;
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

    // Trim alignment to the high-scoring core. WFA runs in global mode
    // and must consume the full reference window, so the surrounding pad
    // bases get scattered as D/X ops through the alignment. Without this
    // trim, identity is artificially capped at ~1 - 2·PAD/(q + 2·PAD).
    // Drop everything outside the first/last run of CORE_RUN consecutive
    // exact-match positions. If no such run exists, fall back to the
    // untrimmed alignment so we don't regress on highly-divergent hits.
    const CORE_RUN: usize = 3;
    let q_aln: Vec<u8> = alignment.query_aligned.bytes().collect();
    let t_aln: Vec<u8> = alignment.text_aligned.bytes().collect();
    debug_assert_eq!(q_aln.len(), t_aln.len());
    let n = q_aln.len();
    let is_match_at = |i: usize| -> bool {
        q_aln[i] != b'-' && t_aln[i] != b'-' && q_aln[i].eq_ignore_ascii_case(&t_aln[i])
    };
    let core_first = if n >= CORE_RUN {
        (0..=n - CORE_RUN).find(|&i| (i..i + CORE_RUN).all(|j| is_match_at(j)))
    } else {
        None
    };
    let (trim_start, trim_end) = match core_first {
        Some(s) => {
            let e = (CORE_RUN..=n)
                .rev()
                .find(|&i| (i - CORE_RUN..i).all(|j| is_match_at(j)))
                .unwrap_or(n);
            (s, e)
        }
        None => (0, n), // fall back: no clean core, use whole alignment
    };
    let pre_trim_target: usize = t_aln[..trim_start].iter().filter(|&&b| b != b'-').count();
    let pre_trim_query: usize = q_aln[..trim_start].iter().filter(|&&b| b != b'-').count();
    let stats = AlignmentStats::from_aligned_strings(
        &q_aln[trim_start..trim_end],
        &t_aln[trim_start..trim_end],
    );

    // text-side leading skip = bases of t consumed before the alignment proper begins
    let aligned_target_len = stats.target_consumed;
    let target_start_in_window = pre_trim_target;
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

    // Query span on the original (un-RC'd) coordinates. We now align the
    // full query (query_min_in_align == 0 after the slicing change above),
    // so coordinate mapping just needs to undo the reverse-complement when
    // applicable.
    let q_aligned_len = stats.query_consumed;
    let _ = query_min_in_align; // placeholder; full query is always 0..query.len()
    let (final_q_start, final_q_end) = if is_reverse {
        // We aligned RC(query) against t.
        // Map back: a RC slice's [a, b) corresponds to original [L-b, L-a)
        // where L = q_str.len().
        let l = q_str.len();
        let a = pre_trim_query;
        let b = a + q_aligned_len;
        let rc_start = l - b;
        let rc_end = l - a;
        (rc_start, rc_end)
    } else {
        (pre_trim_query, pre_trim_query + q_aligned_len)
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
        Self::from_aligned_strings(&q, &t)
    }

    /// Build stats + CIGAR from already-extracted gapped alignment strings.
    /// Caller passes equal-length byte slices where '-' marks gaps.
    fn from_aligned_strings(q: &[u8], t: &[u8]) -> Self {
        debug_assert_eq!(q.len(), t.len(), "aligned strings must be equal length");

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
