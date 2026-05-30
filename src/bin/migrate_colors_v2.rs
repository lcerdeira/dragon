//! Migrate one or more index directories' `colors.drgn` from v1 to v2.
//!
//! v2 stores each distinct color set (RoaringBitmap) exactly once in a shared
//! pool; unitigs reference the pool by a 4-byte set_id.
//!
//! For large, diverse databases (e.g. 26 K-genome S. aureus) most unitigs
//! have near-unique color sets, so deduplication alone saves little.  The
//! real win comes from **high-cardinality truncation** (`--max-cardinality`):
//! unitigs in more than N genomes are replaced with an empty bitmap.
//! These high-cardinality unitigs are already ignored by `containment_rank`
//! (it stops building the candidate set at 2 000 genomes), so truncating them
//! is lossless for search quality while reducing disk/RAM by 3-6×.
//!
//! Recommended invocation for HPC 26 K-genome shards:
//!   migrate_colors_v2 --max-cardinality 2000 <index_dir> [<index_dir> ...]
//!
//! For each directory:
//!   * skips if `colors.drgn` is already v2,
//!   * streams one bitmap at a time (low RAM — no full in-memory collection),
//!   * writes `colors.drgn.v2.tmp`, atomically promotes to `colors.drgn`,
//!   * keeps the original as `colors.drgn.v1.bak`.
//!
//! Usage:
//!   migrate_colors_v2 [--max-cardinality N] <index_dir> [<index_dir> ...]

use std::path::PathBuf;
use anyhow::{Context, Result};

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    if raw_args.is_empty() {
        eprintln!(
            "usage: migrate_colors_v2 [--max-cardinality N] <index_dir> [<index_dir> ...]"
        );
        std::process::exit(2);
    }

    // Parse flags.
    let mut max_cardinality: Option<u32> = None;
    let mut dirs: Vec<PathBuf> = Vec::new();
    let mut i = 0;
    while i < raw_args.len() {
        match raw_args[i].as_str() {
            "--max-cardinality" | "-m" => {
                i += 1;
                let val = raw_args.get(i).ok_or_else(|| {
                    anyhow::anyhow!("--max-cardinality requires a value")
                })?;
                max_cardinality = Some(val.parse::<u32>().with_context(|| {
                    format!("--max-cardinality value '{}' is not a u32", val)
                })?);
            }
            arg if arg.starts_with("--max-cardinality=") => {
                let val = &arg["--max-cardinality=".len()..];
                max_cardinality = Some(val.parse::<u32>().with_context(|| {
                    format!("--max-cardinality value '{}' is not a u32", val)
                })?);
            }
            _ => dirs.push(PathBuf::from(&raw_args[i])),
        }
        i += 1;
    }

    if dirs.is_empty() {
        eprintln!("error: no index directories specified");
        std::process::exit(2);
    }

    if let Some(k) = max_cardinality {
        eprintln!(
            "Truncating bitmaps with cardinality > {} to empty \
             (safety check: will refuse if threshold < 90% of shard genome count).",
            k
        );
    } else {
        eprintln!(
            "Note: no --max-cardinality set. Truncation is only safe when \
             --max-cardinality >= 90% of the shard's genome count. \
             For sharded indexes with few genomes per shard (e.g. 7-shard 26K-genome \
             build, ~3700 genomes/shard) no safe truncation threshold exists — \
             colors.drgn deduplication without truncation has ~1× ratio for \
             diverse databases and is not recommended."
        );
    }

    let mut any_error = false;

    for dir in &dirs {
        let colors_path = dir.join("colors.drgn");
        if !colors_path.exists() {
            eprintln!("[skip] {:?}: no colors.drgn", dir);
            continue;
        }

        let t0 = std::time::Instant::now();
        match dragon::index::color::migrate_colors_v1_to_v2(dir, max_cardinality)
            .with_context(|| format!("migrating {:?}", dir))
        {
            Ok(stats) if stats.skipped => {
                eprintln!("[skip] {:?}: already v2", dir);
            }
            Ok(stats) => {
                let ratio = stats.v1_bytes as f64 / stats.v2_bytes.max(1) as f64;
                eprintln!(
                    "[done] {:?}: {} unitigs, {} unique sets, \
                     {:.2} GB → {:.2} GB ({:.2}×)  wall={:.0}s",
                    dir,
                    stats.num_unitigs,
                    stats.num_sets,
                    stats.v1_bytes as f64 / 1e9,
                    stats.v2_bytes as f64 / 1e9,
                    ratio,
                    t0.elapsed().as_secs_f64(),
                );
            }
            Err(e) => {
                eprintln!("[error] {:?}: {:#}", dir, e);
                any_error = true;
            }
        }
    }

    if any_error {
        std::process::exit(1);
    }
    eprintln!("done.");
    Ok(())
}
