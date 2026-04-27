/// Graph-context output: induced subgraph around alignment hits.
///
/// For each query hit, extracts the local neighbourhood of unitigs in the ccdBG
/// and outputs them in GFA (Graphical Fragment Assembly) format, enabling
/// visualisation of genomic context (e.g., plasmid vs chromosome, nearby
/// mobile elements).

use std::collections::{BTreeSet, HashSet};
use std::io::Write;

use crate::index::color::ColorIndex;
use crate::index::paths::PathIndex;
use crate::index::unitig::UnitigSet;
use crate::io::paf::PafRecord;

/// A node in the induced subgraph (a unitig).
#[derive(Clone, Debug)]
pub struct SubgraphNode {
    pub unitig_id: u32,
    pub genome_ids: Vec<u32>,
    pub is_hit: bool,
    pub sequence: Vec<u8>,
    /// Information content: log2(total_genomes / color_cardinality).
    /// High for genome-specific unitigs, low for core genome.
    pub information_content: f64,
    /// Whether this unitig is private (present in ≤3 genomes).
    pub is_private: bool,
}

/// An edge between two unitigs in the subgraph (adjacency in a genome path).
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct SubgraphEdge {
    pub from_unitig: u32,
    pub to_unitig: u32,
    pub from_orient: char,
    pub to_orient: char,
}

/// An induced subgraph around one or more alignment hits.
#[derive(Clone, Debug)]
pub struct HitSubgraph {
    pub query_name: String,
    pub target_name: String,
    pub nodes: Vec<SubgraphNode>,
    pub edges: Vec<SubgraphEdge>,
}

/// Extract induced subgraphs for alignment hits.
///
/// For each hit, walks the genome path to find the unitig containing the hit,
/// then expands `radius` steps in both directions along the path to capture
/// local genomic context.
pub fn extract_hit_subgraphs(
    records: &[PafRecord],
    path_index: &PathIndex,
    color_index: &ColorIndex,
    unitigs: &UnitigSet,
    radius: usize,
) -> Vec<HitSubgraph> {
    let mut subgraphs = Vec::new();

    for record in records {
        // Find genome in path index.
        // target_name is "genome_{id}" (from align.rs); match by genome_id or name.
        let genome_path = if let Some(id_str) = record.target_name.strip_prefix("genome_") {
            if let Ok(gid) = id_str.parse::<u32>() {
                path_index.get_path(gid)
            } else {
                path_index.iter().find(|p| p.genome_name == record.target_name)
            }
        } else {
            path_index.iter().find(|p| p.genome_name == record.target_name)
        };
        let genome_path = match genome_path {
            Some(p) => p,
            None => continue,
        };

        // Find which step(s) overlap the hit region
        let hit_start = record.target_start as u64;
        let hit_end = record.target_end as u64;

        let mut hit_step_indices: Vec<usize> = Vec::new();
        for (idx, step) in genome_path.steps.iter().enumerate() {
            let unitig_len = unitigs.unitigs.get(step.unitig_id as usize)
                .map(|u| u.sequence.len as u64)
                .unwrap_or(1000);
            let step_end = step.genome_offset + unitig_len;
            if step.genome_offset < hit_end && step_end > hit_start {
                hit_step_indices.push(idx);
            }
        }

        if hit_step_indices.is_empty() {
            continue;
        }

        // Expand by radius steps in both directions
        let min_idx = hit_step_indices
            .iter()
            .copied()
            .min()
            .unwrap_or(0)
            .saturating_sub(radius);
        let max_idx = hit_step_indices
            .iter()
            .copied()
            .max()
            .unwrap_or(0)
            .saturating_add(radius)
            .min(genome_path.steps.len().saturating_sub(1));

        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut seen_unitigs: HashSet<u32> = HashSet::new();
        let hit_unitigs: HashSet<u32> = hit_step_indices
            .iter()
            .map(|&i| genome_path.steps[i].unitig_id)
            .collect();

        for idx in min_idx..=max_idx {
            let step = &genome_path.steps[idx];

            if seen_unitigs.insert(step.unitig_id) {
                // Get genome colors for this unitig
                let colors = color_index
                    .get_colors(step.unitig_id)
                    .ok()
                    .map(|bm| bm.iter().collect::<Vec<u32>>())
                    .unwrap_or_default();

                let color_card = colors.len() as u64;
                let total_genomes = color_index.num_genomes();
                let ic = if color_card > 0 && total_genomes > 0 {
                    (total_genomes as f64 / color_card as f64).log2()
                } else {
                    0.0
                };

                // Get the actual unitig sequence
                let sequence = unitigs.get_sequence(step.unitig_id);

                nodes.push(SubgraphNode {
                    unitig_id: step.unitig_id,
                    genome_ids: colors,
                    is_hit: hit_unitigs.contains(&step.unitig_id),
                    sequence,
                    information_content: ic,
                    is_private: color_card <= 3,
                });
            }

            // Add edge to next step
            if idx < max_idx {
                let next_step = &genome_path.steps[idx + 1];
                let edge = SubgraphEdge {
                    from_unitig: step.unitig_id,
                    to_unitig: next_step.unitig_id,
                    from_orient: if step.is_reverse { '-' } else { '+' },
                    to_orient: if next_step.is_reverse { '-' } else { '+' },
                };
                edges.push(edge);
            }
        }

        // Deduplicate edges
        let edge_set: BTreeSet<(u32, u32, char, char)> = edges
            .iter()
            .map(|e| (e.from_unitig, e.to_unitig, e.from_orient, e.to_orient))
            .collect();
        let edges: Vec<SubgraphEdge> = edge_set
            .into_iter()
            .map(|(f, t, fo, to)| SubgraphEdge {
                from_unitig: f,
                to_unitig: t,
                from_orient: fo,
                to_orient: to,
            })
            .collect();

        subgraphs.push(HitSubgraph {
            query_name: record.query_name.clone(),
            target_name: record.target_name.clone(),
            nodes,
            edges,
        });
    }

    subgraphs
}

/// Write subgraphs in GFA format (one GFA per hit).
///
/// Output includes:
/// - H line: GFA header
/// - S lines: segments (unitigs) with sequence + annotation tags:
///   - LN:i: sequence length
///   - RC:i: genome count (color cardinality)
///   - CL:Z: genome IDs (comma-separated)
///   - TP:Z: HIT (direct alignment) or CTX (context neighbourhood)
///   - IC:f: information content (bits; high = genome-specific)
///   - PV:Z: YES if private unitig (≤3 genomes), NO otherwise
///   - CO:Z: Bandage-compatible hex color (#FF0000=hit, #00AA00=private, #AAAAAA=context)
/// - L lines: links (edges) between unitigs with overlap
/// - Comment lines with subgraph metadata
pub fn write_gfa<W: Write>(
    writer: &mut W,
    subgraphs: &[HitSubgraph],
) -> std::io::Result<()> {
    for (idx, sg) in subgraphs.iter().enumerate() {
        writeln!(writer, "H\tVN:Z:1.0")?;
        writeln!(
            writer,
            "# Subgraph {} for query={} target={}",
            idx, sg.query_name, sg.target_name
        )?;
        writeln!(
            writer,
            "# Nodes: {} ({} HIT, {} CTX, {} private)",
            sg.nodes.len(),
            sg.nodes.iter().filter(|n| n.is_hit).count(),
            sg.nodes.iter().filter(|n| !n.is_hit).count(),
            sg.nodes.iter().filter(|n| n.is_private).count(),
        )?;

        // Segment lines
        for node in &sg.nodes {
            let color_tag = node
                .genome_ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(",");
            let hit_tag = if node.is_hit { "HIT" } else { "CTX" };

            // Bandage-compatible color: red=hit, green=private context, grey=shared context
            let bandage_color = if node.is_hit {
                "#FF0000"
            } else if node.is_private {
                "#00AA00"
            } else {
                "#AAAAAA"
            };

            let seq_str = if node.sequence.is_empty() {
                "*".to_string()
            } else {
                String::from_utf8_lossy(&node.sequence).to_string()
            };

            writeln!(
                writer,
                "S\tunitig_{}\t{}\tLN:i:{}\tRC:i:{}\tCL:Z:{}\tTP:Z:{}\tIC:f:{:.2}\tPV:Z:{}\tCO:Z:{}",
                node.unitig_id,
                seq_str,
                node.sequence.len(),
                node.genome_ids.len(),
                color_tag,
                hit_tag,
                node.information_content,
                if node.is_private { "YES" } else { "NO" },
                bandage_color,
            )?;
        }

        // Link lines
        for edge in &sg.edges {
            writeln!(
                writer,
                "L\tunitig_{}\t{}\tunitig_{}\t{}\t0M",
                edge.from_unitig, edge.from_orient, edge.to_unitig, edge.to_orient,
            )?;
        }

        writeln!(writer)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_gfa_empty() {
        let mut buf = Vec::new();
        write_gfa(&mut buf, &[]).unwrap();
        assert!(buf.is_empty());
    }

    #[test]
    fn test_write_gfa_basic() {
        let sg = HitSubgraph {
            query_name: "query1".to_string(),
            target_name: "genome_A".to_string(),
            nodes: vec![
                SubgraphNode {
                    unitig_id: 10,
                    genome_ids: vec![0, 1, 2],
                    is_hit: true,
                    sequence: b"ACGTACGT".to_vec(),
                    information_content: 5.0,
                    is_private: true,
                },
                SubgraphNode {
                    unitig_id: 11,
                    genome_ids: vec![0, 1],
                    is_hit: false,
                    sequence: b"TTGCAA".to_vec(),
                    information_content: 2.0,
                    is_private: false,
                },
            ],
            edges: vec![SubgraphEdge {
                from_unitig: 10,
                to_unitig: 11,
                from_orient: '+',
                to_orient: '+',
            }],
        };

        let mut buf = Vec::new();
        write_gfa(&mut buf, &[sg]).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("S\tunitig_10\tACGTACGT\tLN:i:8"));
        assert!(output.contains("S\tunitig_11\tTTGCAA\tLN:i:6"));
        assert!(output.contains("L\tunitig_10\t+\tunitig_11\t+\t0M"));
        assert!(output.contains("TP:Z:HIT"));
        assert!(output.contains("TP:Z:CTX"));
        assert!(output.contains("IC:f:5.00"));
        assert!(output.contains("PV:Z:YES"));
        assert!(output.contains("CO:Z:#FF0000")); // hit = red
        assert!(output.contains("CO:Z:#AAAAAA")); // non-private context = grey
    }
}
