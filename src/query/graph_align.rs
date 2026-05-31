//! Graph-aware alignment for paths.bin v4.
//!
//! Instead of linearising a genome window and running WFA, this module aligns
//! a query against the local subgraph extracted from the v4 path index.  The
//! graph structure — alternative branches in the dBG — is captured as directed
//! edges between [`SubgraphNode`]s so future versions can exploit branching
//! paths.  For v1 the alignment itself is standard Needleman–Wunsch edit-
//! distance DP over the linearised `ref_seq`, with banded optimisation
//! (`O(n × max_edit)` instead of `O(n²)`).
//!
//! ## Why this helps over raw WFA
//!
//! * The extraction window is driven by the **graph** (anchor unitigs →
//!   genome positions) rather than raw coordinate arithmetic, so the window
//!   edges are always unitig-aligned.
//! * The subgraph structure is preserved for future branching-aware DP.
//! * Edit-distance DP handles affine-gap queries more uniformly than
//!   WFA's gap-affine scoring for very short (≤100 bp) queries.

use std::collections::HashMap;

use crate::index::paths::PathIndex;
use crate::index::unitig::UnitigSet;

// ─── Data structures ──────────────────────────────────────────────────────────

/// One node in the local subgraph — one traversal of a unitig in the genome
/// path.
#[derive(Debug, Clone)]
pub struct SubgraphNode {
    /// Unitig identifier.
    pub unitig_id: u32,
    /// Genome-coordinate offset where this traversal starts.
    pub genome_offset: u64,
    /// Length (bp) contributed by this traversal (after deduplication of the
    /// `k-1` boundary overlap with the preceding step).
    pub length: usize,
    /// Whether the genome reads this unitig on the reverse-complement strand.
    pub is_reverse: bool,
    /// Index of this node in [`LocalSubgraph::nodes`].
    pub node_idx: usize,
}

/// The extracted local subgraph around the anchor unitigs.
#[derive(Debug)]
pub struct LocalSubgraph {
    /// Ordered list of nodes (one per path step in the extraction window).
    pub nodes: Vec<SubgraphNode>,
    /// Directed edges: `(from_idx, to_idx)` pairs in `nodes`.
    pub edges: Vec<(usize, usize)>,
    /// Concatenated reference sequence assembled from the window steps.
    /// For v1 this is the linearised window sequence used as alignment target.
    pub ref_seq: Vec<u8>,
    /// Maps `unitig_id → node_idx` for the first occurrence in the window.
    pub unitig_to_node: HashMap<u32, usize>,
    /// Genome coordinate of the first base in `ref_seq`.
    pub window_start: u64,
}

/// Result of graph DP alignment.
#[derive(Debug, Clone)]
pub struct GraphAlignment {
    /// Number of matching bases.
    pub num_matches: usize,
    /// Total CIGAR length (matches + mismatches + insertions + deletions).
    pub alignment_len: usize,
    /// Edit distance (mismatches + indel bases).
    pub edit_distance: usize,
    /// CIGAR string using `=`, `X`, `I`, `D` operations.
    pub cigar: String,
    /// Indices into [`LocalSubgraph::nodes`] traversed by the alignment path.
    pub path: Vec<usize>,
    /// Number of bases consumed from the reference (target side).
    pub ref_consumed: usize,
    /// Number of bases consumed from the query.
    pub query_consumed: usize,
}

// ─── Integration check ────────────────────────────────────────────────────────

/// Returns true if graph alignment should be preferred over WFA for this
/// `PathIndex` variant.  Currently true only for `PathIndex::MmapV4`.
#[inline]
pub fn should_use_graph_align(path_index: &PathIndex) -> bool {
    matches!(path_index, PathIndex::MmapV4(_))
}

// ─── Subgraph extraction ──────────────────────────────────────────────────────

/// Extract a local subgraph around `anchor_unitig_ids` and assemble the
/// reference sequence.
///
/// # Algorithm (v1 — linear)
///
/// 1. Use [`PathIndex::find_unitig_steps`] to locate the anchor steps in the
///    genome path; compute the median anchor genome position.
/// 2. Compute an extraction window:
///    `[anchor_pos - query_len/2, anchor_pos + query_len * 1.5]`,
///    clamped to `[0, genome_len)`.  The window is capped at
///    `3 × query_len` to prevent quadratic DP blowup.
/// 3. Walk [`PathIndex::iter_window`] over the window to build
///    `SubgraphNode` entries and `edges`.
/// 4. Assemble `ref_seq` from the window steps via
///    [`PathIndex::extract_sequence_static`], which handles `k-1` boundary
///    deduplication.
///
/// Returns `None` if the subgraph is empty or no anchor steps were found.
pub fn extract_query_subgraph(
    path_index: &PathIndex,
    unitigs: &UnitigSet,
    genome_id: u32,
    genome_len: u64,
    anchor_unitig_ids: &[u32],
    query_len: usize,
) -> Option<LocalSubgraph> {
    if anchor_unitig_ids.is_empty() || genome_len == 0 {
        return None;
    }

    // Step 1 — locate anchor steps (early-exit walk for v4).
    let wanted: std::collections::HashSet<u32> = anchor_unitig_ids.iter().copied().collect();
    let anchor_steps = path_index.find_unitig_steps(genome_id, &wanted);
    if anchor_steps.is_empty() {
        return None;
    }

    // Step 2 — compute median anchor genome position, then the extraction
    // window.  We cap the window at 3 × query_len so the DP stays O(n · d).
    let mut positions: Vec<u64> = anchor_steps
        .values()
        .map(|s| s.genome_offset)
        .collect();
    positions.sort_unstable();
    let median_pos = positions[positions.len() / 2];

    let half_back = (query_len / 2) as u64;
    let window_start = median_pos.saturating_sub(half_back);

    let max_window = (query_len * 3).max(400) as u64;
    let window_end_unclamped = window_start + max_window;
    let window_end = window_end_unclamped.min(genome_len);

    if window_start >= window_end {
        return None;
    }

    // Step 3 — walk the window to build nodes + edges.
    let steps = path_index.iter_window(genome_id, window_start, window_end);
    if steps.is_empty() {
        return None;
    }

    let mut nodes: Vec<SubgraphNode> = Vec::with_capacity(steps.len());
    let mut edges: Vec<(usize, usize)> = Vec::with_capacity(steps.len().saturating_sub(1));
    let mut unitig_to_node: HashMap<u32, usize> = HashMap::with_capacity(steps.len());

    for (i, step) in steps.iter().enumerate() {
        let unitig_len = unitigs
            .unitigs
            .get(step.unitig_id as usize)
            .map(|u| u.sequence.len as usize)
            .unwrap_or(0);

        let node = SubgraphNode {
            unitig_id: step.unitig_id,
            genome_offset: step.genome_offset,
            length: unitig_len,
            is_reverse: step.is_reverse,
            node_idx: i,
        };
        nodes.push(node);

        // First occurrence in window → register in lookup table.
        unitig_to_node.entry(step.unitig_id).or_insert(i);

        // Sequential edge from previous node to this one.
        if i > 0 {
            edges.push((i - 1, i));
        }
    }

    // Step 4 — assemble the reference sequence.
    let ref_seq = PathIndex::extract_sequence_static(&steps, window_start, window_end, unitigs);

    if ref_seq.is_empty() {
        return None;
    }

    Some(LocalSubgraph {
        nodes,
        edges,
        ref_seq,
        unitig_to_node,
        window_start,
    })
}

// ─── DP alignment ─────────────────────────────────────────────────────────────

/// Align `query` against `subgraph.ref_seq` using banded Levenshtein DP.
///
/// Uses unit costs: match = 0, substitution = 1, insertion = 1, deletion = 1.
/// The band width is `max_edit + 1`, keeping the DP at `O(n × max_edit)`
/// instead of `O(n²)`.
///
/// Returns `None` if:
/// * The subgraph has no reference sequence.
/// * `max_edit == 0` and the sequences are not identical.
/// * The optimal edit distance exceeds `max_edit`.
pub fn graph_dp_align(
    query: &[u8],
    subgraph: &LocalSubgraph,
    max_edit: usize,
) -> Option<GraphAlignment> {
    let ref_seq = &subgraph.ref_seq;
    if ref_seq.is_empty() || query.is_empty() {
        return None;
    }

    let q_len = query.len();
    let r_len = ref_seq.len();

    // Fill the banded DP table.
    // dp[i][j] = min edit distance to align query[0..i] vs ref[0..j].
    // Band: only compute cells where |i - j| <= max_edit.
    // We use a flat Vec<u32> of size (q_len+1) * (r_len+1) for simplicity;
    // an optimised implementation would use only 2 rows.
    //
    // For the sizes we expect (query ≤ ~300 bp, ref ≤ 3 × 300 = 900 bp,
    // max_edit ≤ ~50), this is 901 × 301 × 4 = ~1 MB — acceptable.
    // Banding eliminates most cells from consideration.

    let max_e = max_edit as u32;
    let inf: u32 = max_e + 1;

    // Use two-row rolling DP to save memory.
    let mut prev: Vec<u32> = (0..=r_len as u32).collect();
    // Clamp row 0 to the band: cells further than max_edit from column 0
    // are unreachable.
    for j in (max_edit + 1)..=r_len {
        prev[j] = inf;
    }

    // For traceback we must store the full table.
    // Allocate (q_len+1) rows each of (r_len+1) cells.
    let mut table: Vec<Vec<u32>> = Vec::with_capacity(q_len + 1);
    table.push(prev.clone());

    for i in 1..=q_len {
        let qi = query[i - 1];
        let mut curr = vec![inf; r_len + 1];

        // First column: aligning query[0..i] against empty ref.
        curr[0] = if (i as u32) <= max_e { i as u32 } else { inf };

        // Column range for the band: j in [max(1, i - max_edit),
        //                                   (i + max_edit).min(r_len)].
        // j=0 is the first column and is already initialised above, so we
        // start at 1 to avoid a j-1 underflow.
        let j_lo = i.saturating_sub(max_edit).max(1);
        let j_hi = (i + max_edit).min(r_len);

        for j in j_lo..=j_hi {
            let rj = ref_seq[j - 1];
            let cost_sub = if qi.eq_ignore_ascii_case(&rj) { 0 } else { 1 };

            let from_diag = prev[j - 1].saturating_add(cost_sub);
            let from_del = prev[j].saturating_add(1); // gap in query (deletion)
            let from_ins = curr[j - 1].saturating_add(1); // gap in ref (insertion)

            curr[j] = from_diag.min(from_del).min(from_ins).min(inf);
        }

        table.push(curr.clone());
        prev = curr;
    }

    // Find the best ending cell: minimum in the last row (semi-global:
    // we allow the alignment to end at any column of the last query row,
    // effectively allowing free right-flank deletions in the reference).
    // This matches the WFA behaviour of aligning the full query against
    // a reference that may be longer on the right.
    let last_row = &table[q_len];
    let (best_j, best_dist) = last_row
        .iter()
        .enumerate()
        .min_by_key(|(_, &v)| v)
        .map(|(j, &v)| (j, v))?;

    if best_dist > max_e {
        return None;
    }

    // Traceback from (q_len, best_j) to (0, 0).
    let (cigar, num_matches, q_consumed, r_consumed, _mismatches) =
        traceback(&table, query, ref_seq, q_len, best_j);

    // Map reference positions to node indices for the path field.
    let path = nodes_covering_ref(&subgraph.nodes, r_consumed, subgraph.window_start);

    let edits = best_dist as usize;
    // CIGAR-column count (sum of all run lengths).
    let alignment_len = cigar_len(&cigar);

    Some(GraphAlignment {
        num_matches,
        alignment_len,
        edit_distance: edits,
        cigar,
        path,
        ref_consumed: r_consumed,
        query_consumed: q_consumed,
    })
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Traceback through the DP table from `(i, j)` to `(0, 0)`.
/// Returns `(cigar, num_matches, q_consumed, r_consumed, mismatches)`.
fn traceback(
    table: &[Vec<u32>],
    query: &[u8],
    ref_seq: &[u8],
    mut i: usize,
    mut j: usize,
) -> (String, usize, usize, usize, usize) {
    let mut ops: Vec<u8> = Vec::with_capacity(i + j);
    let mut matches = 0usize;
    let mut mismatches = 0usize;

    while i > 0 || j > 0 {
        if i == 0 {
            // Only deletions left (consume ref).
            ops.push(b'D');
            j -= 1;
            continue;
        }
        if j == 0 {
            // Only insertions left (consume query).
            ops.push(b'I');
            i -= 1;
            continue;
        }

        let qi = query[i - 1];
        let rj = ref_seq[j - 1];
        let cost_sub = if qi.eq_ignore_ascii_case(&rj) { 0u32 } else { 1 };

        let diag_val = table[i - 1][j - 1].saturating_add(cost_sub);
        let del_val = table[i - 1][j].saturating_add(1); // insertion in ref coords (query advances)
        let _ins_val = table[i][j - 1].saturating_add(1); // deletion in ref coords (ref advances)

        let cur = table[i][j];

        if cur == diag_val {
            if cost_sub == 0 {
                ops.push(b'=');
                matches += 1;
            } else {
                ops.push(b'X');
                mismatches += 1;
            }
            i -= 1;
            j -= 1;
        } else if cur == del_val {
            // query base consumes without a ref base (insertion in query = deletion in alignment)
            ops.push(b'I');
            i -= 1;
        } else {
            // ref base consumes without a query base (deletion in query = insertion in alignment)
            ops.push(b'D');
            j -= 1;
        }
    }

    ops.reverse();

    let q_consumed = ops.iter().filter(|&&o| o != b'D').count();
    let r_consumed = ops.iter().filter(|&&o| o != b'I').count();

    let cigar = encode_cigar(&ops);
    (cigar, matches, q_consumed, r_consumed, mismatches)
}

/// RLE-encode a flat ops vector into a compact CIGAR string.
fn encode_cigar(ops: &[u8]) -> String {
    if ops.is_empty() {
        return String::new();
    }
    let mut out = String::with_capacity(ops.len() / 2);
    let mut run_op = ops[0];
    let mut run_len: usize = 1;
    for &op in &ops[1..] {
        if op == run_op {
            run_len += 1;
        } else {
            out.push_str(&run_len.to_string());
            out.push(run_op as char);
            run_op = op;
            run_len = 1;
        }
    }
    out.push_str(&run_len.to_string());
    out.push(run_op as char);
    out
}

/// Count total CIGAR columns (sum of all run lengths).
fn cigar_len(cigar: &str) -> usize {
    let mut total = 0usize;
    let mut num_start: Option<usize> = None;
    for (idx, ch) in cigar.char_indices() {
        if ch.is_ascii_digit() {
            if num_start.is_none() {
                num_start = Some(idx);
            }
        } else {
            if let Some(ns) = num_start.take() {
                if let Ok(n) = cigar[ns..idx].parse::<usize>() {
                    total += n;
                }
            }
        }
    }
    total
}

/// Return the sequence of node indices in `nodes` that cover the first
/// `ref_consumed` bases of the reference starting from `window_start`.
fn nodes_covering_ref(
    nodes: &[SubgraphNode],
    ref_consumed: usize,
    window_start: u64,
) -> Vec<usize> {
    let window_end = window_start + ref_consumed as u64;
    nodes
        .iter()
        .filter(|n| {
            let node_end = n.genome_offset + n.length as u64;
            n.genome_offset < window_end && node_end > window_start
        })
        .map(|n| n.node_idx)
        .collect()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── encode_cigar ──────────────────────────────────────────────────────────

    #[test]
    fn encode_cigar_all_matches() {
        let ops = vec![b'=', b'=', b'=', b'='];
        assert_eq!(encode_cigar(&ops), "4=");
    }

    #[test]
    fn encode_cigar_mixed_ops() {
        let ops = vec![b'=', b'=', b'X', b'=', b'I', b'D', b'D'];
        assert_eq!(encode_cigar(&ops), "2=1X1=1I2D");
    }

    #[test]
    fn encode_cigar_empty() {
        assert_eq!(encode_cigar(&[]), "");
    }

    // ── cigar_len ─────────────────────────────────────────────────────────────

    #[test]
    fn cigar_len_basic() {
        assert_eq!(cigar_len("4="), 4);
        assert_eq!(cigar_len("2=1X1="), 4);
        assert_eq!(cigar_len("2=1X1=1I2D"), 7);
    }

    // ── graph_dp_align: perfect match ─────────────────────────────────────────

    #[test]
    fn dp_align_perfect_match() {
        // query == ref_seq → edit distance 0, all matches.
        let seq = b"ACGTACGT".to_vec();
        let subgraph = LocalSubgraph {
            nodes: vec![],
            edges: vec![],
            ref_seq: seq.clone(),
            unitig_to_node: HashMap::new(),
            window_start: 0,
        };
        let ga = graph_dp_align(&seq, &subgraph, 0).expect("perfect match should succeed");
        assert_eq!(ga.edit_distance, 0);
        assert_eq!(ga.num_matches, 8);
        assert_eq!(ga.cigar, "8=");
    }

    // ── graph_dp_align: single substitution ───────────────────────────────────

    #[test]
    fn dp_align_single_mismatch() {
        let query = b"ACTTACGT".to_vec(); // position 2: T instead of G
        let reference = b"ACGTACGT".to_vec();
        let subgraph = LocalSubgraph {
            nodes: vec![],
            edges: vec![],
            ref_seq: reference,
            unitig_to_node: HashMap::new(),
            window_start: 0,
        };
        let ga = graph_dp_align(&query, &subgraph, 2).expect("single mismatch within max_edit");
        assert_eq!(ga.edit_distance, 1);
        assert_eq!(ga.num_matches, 7);
        // CIGAR should contain one X.
        assert!(ga.cigar.contains('X'), "cigar='{}' should contain X", ga.cigar);
    }

    // ── graph_dp_align: exceeds max_edit → None ───────────────────────────────

    #[test]
    fn dp_align_exceeds_max_edit() {
        let query = b"AAAAAAAA".to_vec();
        let reference = b"TTTTTTTT".to_vec(); // 8 mismatches
        let subgraph = LocalSubgraph {
            nodes: vec![],
            edges: vec![],
            ref_seq: reference,
            unitig_to_node: HashMap::new(),
            window_start: 0,
        };
        // max_edit=3 → should return None (8 mismatches > 3).
        assert!(graph_dp_align(&query, &subgraph, 3).is_none());
    }

    // ── graph_dp_align: query shorter than ref (semi-global) ─────────────────

    #[test]
    fn dp_align_query_shorter_than_ref() {
        // Query is a prefix of the reference.  Semi-global alignment should
        // produce 0 edit distance (free right-flank on the reference).
        let query = b"ACGT".to_vec();
        let reference = b"ACGTTTTTTTTT".to_vec();
        let subgraph = LocalSubgraph {
            nodes: vec![],
            edges: vec![],
            ref_seq: reference,
            unitig_to_node: HashMap::new(),
            window_start: 0,
        };
        let ga = graph_dp_align(&query, &subgraph, 0).expect("prefix match should succeed");
        assert_eq!(ga.edit_distance, 0);
        assert_eq!(ga.num_matches, 4);
    }

    // ── graph_dp_align: insertion in query ───────────────────────────────────

    #[test]
    fn dp_align_insertion_in_query() {
        // query has an extra 'A' in the middle.
        let query = b"ACAAGT".to_vec();
        let reference = b"ACAGT".to_vec();
        let subgraph = LocalSubgraph {
            nodes: vec![],
            edges: vec![],
            ref_seq: reference,
            unitig_to_node: HashMap::new(),
            window_start: 0,
        };
        let ga = graph_dp_align(&query, &subgraph, 2).expect("1 insertion within max_edit");
        assert_eq!(ga.edit_distance, 1);
        assert!(ga.cigar.contains('I'), "cigar='{}' should contain I", ga.cigar);
    }

    // ── should_use_graph_align ────────────────────────────────────────────────

    #[test]
    fn should_use_graph_align_for_eager_is_false() {
        let eager = PathIndex::Eager(std::sync::Arc::new(vec![]));
        assert!(!should_use_graph_align(&eager));
    }

    // ── nodes_covering_ref ────────────────────────────────────────────────────

    #[test]
    fn nodes_covering_ref_basic() {
        // Two nodes: [0,10) and [10,20).  ref_consumed = 15, window_start = 0.
        let nodes = vec![
            SubgraphNode {
                unitig_id: 0,
                genome_offset: 0,
                length: 10,
                is_reverse: false,
                node_idx: 0,
            },
            SubgraphNode {
                unitig_id: 1,
                genome_offset: 10,
                length: 10,
                is_reverse: false,
                node_idx: 1,
            },
        ];
        let path = nodes_covering_ref(&nodes, 15, 0);
        assert_eq!(path, vec![0, 1]);
    }

    #[test]
    fn nodes_covering_ref_only_first() {
        let nodes = vec![
            SubgraphNode {
                unitig_id: 0,
                genome_offset: 0,
                length: 10,
                is_reverse: false,
                node_idx: 0,
            },
            SubgraphNode {
                unitig_id: 1,
                genome_offset: 10,
                length: 10,
                is_reverse: false,
                node_idx: 1,
            },
        ];
        // ref_consumed = 5 → only the first node is touched.
        let path = nodes_covering_ref(&nodes, 5, 0);
        assert_eq!(path, vec![0]);
    }
}
