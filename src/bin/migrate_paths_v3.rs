//! Migrate one or more index directories' `paths.bin` from v2 to v3
//! (graph-edge encoding). See `docs/design/paths-bin-v3.md`.
//!
//! Usage:
//!   migrate_paths_v3 <index_dir> [<index_dir> ...]
//!
//! For each directory:
//!   * skips if `paths.bin` is already v3,
//!   * writes `paths.bin.v3.tmp`, then atomically installs it as `paths.bin`,
//!   * keeps the original as `paths.bin.v2.bak` (delete manually once the
//!     v3 search has been validated).

use std::path::PathBuf;

use anyhow::{bail, Context, Result};

use dragon::index::paths_v3::{is_v3, migrate_v2_to_v3};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: migrate_paths_v3 <index_dir> [<index_dir> ...]");
        std::process::exit(2);
    }

    for dir in &args[1..] {
        let dir = PathBuf::from(dir);
        let paths = dir.join("paths.bin");
        if !paths.exists() {
            eprintln!("[skip] {dir:?}: no paths.bin");
            continue;
        }
        if is_v3(&paths)? {
            eprintln!("[skip] {dir:?}: paths.bin already v3");
            continue;
        }

        let tmp = dir.join("paths.bin.v3.tmp");
        let bak = dir.join("paths.bin.v2.bak");
        eprintln!("[migrate] {dir:?} ...");

        let stats = migrate_v2_to_v3(&paths, &tmp)
            .with_context(|| format!("migrate {paths:?}"))?;

        let ratio = stats.old_size as f64 / stats.new_size.max(1) as f64;
        eprintln!(
            "  {} genomes, {} unitigs, {} edges",
            stats.num_genomes, stats.num_unitigs, stats.num_edges
        );
        eprintln!(
            "  {:.2} GB -> {:.2} GB  ({:.2}x smaller)",
            stats.old_size as f64 / 1e9,
            stats.new_size as f64 / 1e9,
            ratio
        );

        if stats.new_size == 0 {
            bail!("migration produced an empty file for {dir:?}");
        }

        std::fs::rename(&paths, &bak)
            .with_context(|| format!("back up {paths:?} -> {bak:?}"))?;
        std::fs::rename(&tmp, &paths)
            .with_context(|| format!("install {tmp:?} -> {paths:?}"))?;
        eprintln!("  installed v3; v2 kept at {bak:?}");
    }
    eprintln!("done.");
    Ok(())
}
