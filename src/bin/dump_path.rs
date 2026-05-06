//! Diagnostic tool: dump a single genome's path entry from `paths.bin`
//! and check whether its unitigs are coloured for that genome ID in
//! `colors.drgn`.
//!
//! Usage:
//!   dump_path <index_dir> <genome_name>            # one genome
//!   dump_path <index_dir> <name_a> <name_b>        # diff two
//!
//! Prints:
//!   - genome_id, num_steps, declared genome_length
//!   - first/last 5 steps (genome_offset, unitig_id, is_reverse)
//!   - implied coverage span and any gaps between consecutive steps
//!   - for the first 20 unique unitig ids in the path, whether the
//!     genome's id is set in colors.drgn (so seeding-via-FM can reach
//!     this genome through those unitigs)
//!
//! The two-name form is the saureus debug case: pass a SAMN that Dragon
//! finds and a SAMEA that it doesn't, and compare. If both look
//! structurally similar but only the SAMEA's unitigs are missing from
//! `colors.drgn`, the bug is in the build pipeline's color/path
//! ordering. If the SAMEA's path itself is short/empty/has gaps,
//! it's path-builder.

use std::collections::BTreeSet;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};

use dragon::index::color::load_color_index;
use dragon::index::paths::PathIndex;
use dragon::index::paths::{load_path_index, GenomePath};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 || args.len() > 4 {
        eprintln!("usage: dump_path <index_dir> <genome_name> [<other_genome_name>]");
        std::process::exit(2);
    }
    let index_dir = PathBuf::from(&args[1]);
    let name_a = &args[2];
    let name_b = args.get(3).cloned();

    let pi = load_path_index(&index_dir).context("load paths.bin")?;
    let ci = load_color_index(&index_dir).context("load colors.drgn")?;

    eprintln!(
        "loaded index: {} genomes in paths.bin, {} unitigs in colors.drgn",
        pi.num_genomes(),
        ci.num_unitigs()
    );

    let path_a = find_genome(&pi, name_a)?;
    report(&pi, &ci, &path_a, "A");

    if let Some(b) = name_b {
        let path_b = find_genome(&pi, &b)?;
        report(&pi, &ci, &path_b, "B");
    }
    Ok(())
}

fn find_genome(pi: &PathIndex, name: &str) -> Result<GenomePath> {
    // Linear scan — fine for diagnostic; eager and mmap iters both work.
    for g in pi.iter() {
        if g.genome_name == name {
            return Ok(g);
        }
    }
    bail!("genome {name:?} not found in paths.bin");
}

fn report(
    _pi: &PathIndex,
    ci: &dragon::index::color::ColorIndex,
    g: &GenomePath,
    label: &str,
) {
    println!("==== [{label}] {} ====", g.genome_name);
    println!("genome_id      = {}", g.genome_id);
    println!("genome_length  = {} bp (declared)", g.genome_length);
    println!("num_steps      = {}", g.steps.len());
    if g.steps.is_empty() {
        println!("(no steps — genome has zero unitig coverage)");
        println!();
        return;
    }
    let first_offset = g.steps.first().unwrap().genome_offset;
    let last_offset = g.steps.last().unwrap().genome_offset;
    println!("first step at  = {first_offset}");
    println!("last step at   = {last_offset}");

    println!("first 5 steps:");
    for s in g.steps.iter().take(5) {
        println!(
            "  off={:>10}  unitig={:>8}  rev={}",
            s.genome_offset, s.unitig_id, s.is_reverse
        );
    }
    if g.steps.len() > 5 {
        println!("  ...");
        println!("last 5 steps:");
        for s in g.steps.iter().rev().take(5).collect::<Vec<_>>().iter().rev() {
            println!(
                "  off={:>10}  unitig={:>8}  rev={}",
                s.genome_offset, s.unitig_id, s.is_reverse
            );
        }
    }

    // Coverage gap scan: any consecutive steps whose
    // (offset_{i+1} - offset_i) is unusually large?
    let mut deltas: Vec<u64> = Vec::with_capacity(g.steps.len());
    for w in g.steps.windows(2) {
        deltas.push(w[1].genome_offset.saturating_sub(w[0].genome_offset));
    }
    if !deltas.is_empty() {
        let max_d = *deltas.iter().max().unwrap();
        let big: Vec<(usize, u64)> = deltas
            .iter()
            .enumerate()
            .filter(|(_, &d)| d > 1000)
            .map(|(i, &d)| (i, d))
            .take(10)
            .collect();
        println!("max step delta = {max_d} bp; deltas > 1000 bp: {} (showing up to 10)", big.len());
        for (i, d) in big {
            let s0 = &g.steps[i];
            let s1 = &g.steps[i + 1];
            println!(
                "  gap@i={}: {} -> {} ({} bp), unitigs {} -> {}",
                i, s0.genome_offset, s1.genome_offset, d, s0.unitig_id, s1.unitig_id
            );
        }
    }

    // Color-index reachability: for each unique unitig in the first 20
    // steps, is this genome's id in the color bitmap? If a unitig from
    // this genome's path lacks the genome's id in colors.drgn, FM seeds
    // landing in that unitig will NOT reach this genome via containment
    // ranking — that's exactly the symptom on saureus.
    let unique_unitigs: BTreeSet<u32> = g
        .steps
        .iter()
        .take(20)
        .map(|s| s.unitig_id)
        .collect();
    println!(
        "color-index reachability for first {} unique unitigs in the path:",
        unique_unitigs.len()
    );
    let mut reachable = 0usize;
    let mut unreachable = 0usize;
    for u in &unique_unitigs {
        match ci.get_colors(*u) {
            Ok(bm) => {
                let has = bm.contains(g.genome_id);
                if has {
                    reachable += 1;
                } else {
                    unreachable += 1;
                    println!(
                        "  MISSING: unitig {} does NOT have genome_id {} in its color bitmap (cardinality={})",
                        u, g.genome_id, bm.len()
                    );
                }
            }
            Err(e) => {
                println!("  WARN: get_colors({u}) failed: {e:?}");
            }
        }
    }
    println!(
        "summary: {}/{} unitigs reachable via colors.drgn (unreachable={})",
        reachable,
        unique_unitigs.len(),
        unreachable
    );
    println!();
}
