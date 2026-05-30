//! Migrate one or more index directories' `paths.bin` from v3 to v4
//! (graph-edge encoding with per-genome checkpoints). See
//! `docs/design/paths-bin-v4-lazy-access.md`.
//!
//! Usage:
//!   migrate_paths_v4 <index_dir> [<index_dir> ...]
//!
//! For each directory:
//!   * skips if `paths.bin` is already v4,
//!   * writes `paths.bin.v4.tmp`, then atomically installs it as `paths.bin`,
//!   * keeps the existing file as `paths.bin.v3.bak` (delete manually once
//!     a v4 search has been validated end-to-end).

use std::path::PathBuf;

use anyhow::{bail, Context, Result};

use dragon::index::paths_v4::{is_v4, migrate_v3_to_v4};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: migrate_paths_v4 <index_dir> [<index_dir> ...]");
        std::process::exit(2);
    }

    for dir in &args[1..] {
        let dir = PathBuf::from(dir);
        let paths = dir.join("paths.bin");
        if !paths.exists() {
            eprintln!("[skip] {dir:?}: no paths.bin");
            continue;
        }
        if is_v4(&paths)? {
            eprintln!("[skip] {dir:?}: paths.bin already v4");
            continue;
        }

        let tmp = dir.join("paths.bin.v4.tmp");
        let bak = dir.join("paths.bin.v3.bak");
        eprintln!("[migrate] {dir:?} ...");

        let stats = migrate_v3_to_v4(&paths, &tmp)
            .with_context(|| format!("migrate {paths:?}"))?;

        let ratio = stats.old_size as f64 / stats.new_size.max(1) as f64;
        eprintln!("  {} genomes, {} unitigs", stats.num_genomes, stats.num_unitigs);
        eprintln!(
            "  {:.2} GB -> {:.2} GB  ({:.2}x)",
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
        eprintln!("  installed v4; v3 kept at {bak:?}");
    }
    eprintln!("done.");
    Ok(())
}
