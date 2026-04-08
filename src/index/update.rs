/// Incremental index updates: add new genomes without full rebuild.
///
/// Dragon supports appending new genomes to an existing index via an overlay
/// approach. New genomes are indexed into a separate "overlay" directory, and
/// queries are run against both the base index and all overlays, with results
/// merged by score.
///
/// When the overlay grows too large (e.g., >10% of base), a full rebuild
/// merges everything into a new base index.
///
/// Layout:
///   index/
///     metadata.json     — base index metadata
///     colors.drgn       — base color index
///     fm_index.bin      — base FM-index
///     paths.bin         — base path index
///     specificity.drgn  — base specificity index
///     overlays/
///       overlay_001/    — first overlay (same structure as base)
///       overlay_002/    — second overlay
///       ...
///     overlay_manifest.json  — tracks overlay count and genome offsets

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// Manifest tracking overlay state.
#[derive(Serialize, Deserialize, Debug)]
pub struct OverlayManifest {
    /// Number of overlays.
    pub overlay_count: usize,
    /// Total genomes across all overlays.
    pub overlay_genomes: usize,
    /// Genome ID offset for the next overlay (base num_genomes + cumulative overlay genomes).
    pub next_genome_offset: usize,
    /// Base index genome count.
    pub base_genomes: usize,
    /// Overlay directory names.
    pub overlays: Vec<OverlayEntry>,
}

/// A single overlay entry.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct OverlayEntry {
    pub name: String,
    pub num_genomes: usize,
    pub genome_offset: usize,
    pub created: String,
}

impl OverlayManifest {
    /// Load manifest from an index directory, or create a new one.
    pub fn load_or_create(index_dir: &Path) -> Result<Self> {
        let manifest_path = index_dir.join("overlay_manifest.json");
        if manifest_path.exists() {
            let contents = std::fs::read_to_string(&manifest_path)?;
            Ok(serde_json::from_str(&contents)?)
        } else {
            // Read base metadata to get genome count
            let metadata = crate::index::metadata::load_metadata(index_dir)?;
            Ok(Self {
                overlay_count: 0,
                overlay_genomes: 0,
                next_genome_offset: metadata.num_genomes,
                base_genomes: metadata.num_genomes,
                overlays: Vec::new(),
            })
        }
    }

    /// Save manifest to disk.
    pub fn save(&self, index_dir: &Path) -> Result<()> {
        let manifest_path = index_dir.join("overlay_manifest.json");
        let contents = serde_json::to_string_pretty(self)?;
        std::fs::write(&manifest_path, contents)?;
        Ok(())
    }

    /// Total genomes across base + all overlays.
    pub fn total_genomes(&self) -> usize {
        self.base_genomes + self.overlay_genomes
    }

    /// Whether the overlay is large enough to warrant a full rebuild.
    pub fn needs_compaction(&self, threshold_ratio: f64) -> bool {
        if self.base_genomes == 0 {
            return false;
        }
        (self.overlay_genomes as f64 / self.base_genomes as f64) > threshold_ratio
    }

    /// List all overlay directories.
    pub fn overlay_dirs(&self, index_dir: &Path) -> Vec<PathBuf> {
        let overlays_dir = index_dir.join("overlays");
        self.overlays
            .iter()
            .map(|e| overlays_dir.join(&e.name))
            .collect()
    }
}

/// Add new genomes to an existing index as an overlay.
///
/// Builds a mini-index from the new genomes and registers it in the manifest.
pub fn add_genomes(
    index_dir: &Path,
    new_genomes_dir: &Path,
    kmer_size: usize,
    threads: usize,
    max_ram_bytes: Option<usize>,
) -> Result<OverlayEntry> {
    log::info!("Adding genomes from {:?} to index {:?}", new_genomes_dir, index_dir);

    // Load or create manifest
    let mut manifest = OverlayManifest::load_or_create(index_dir)?;

    // Create overlay directory
    let overlay_name = format!("overlay_{:03}", manifest.overlay_count + 1);
    let overlays_dir = index_dir.join("overlays");
    let overlay_dir = overlays_dir.join(&overlay_name);
    std::fs::create_dir_all(&overlay_dir)?;

    log::info!(
        "Building overlay {} (genome offset: {})",
        overlay_name,
        manifest.next_genome_offset
    );

    // Build a mini-index for the new genomes
    crate::index::build_index_with_options(
        new_genomes_dir,
        &overlay_dir,
        kmer_size,
        threads,
        max_ram_bytes,
    )?;

    // Read the overlay's metadata to get genome count
    let overlay_metadata = crate::index::metadata::load_metadata(&overlay_dir)?;
    let new_genomes = overlay_metadata.num_genomes;

    // Create overlay entry
    let entry = OverlayEntry {
        name: overlay_name.clone(),
        num_genomes: new_genomes,
        genome_offset: manifest.next_genome_offset,
        created: chrono_now(),
    };

    // Update manifest
    manifest.overlays.push(entry.clone());
    manifest.overlay_count += 1;
    manifest.overlay_genomes += new_genomes;
    manifest.next_genome_offset += new_genomes;
    manifest.save(index_dir)?;

    log::info!(
        "Overlay {} added: {} new genomes (total: {} base + {} overlay = {})",
        overlay_name,
        new_genomes,
        manifest.base_genomes,
        manifest.overlay_genomes,
        manifest.total_genomes()
    );

    if manifest.needs_compaction(0.1) {
        log::warn!(
            "Overlay size ({} genomes) exceeds 10% of base ({} genomes). \
             Consider running `dragon index` to rebuild for optimal performance.",
            manifest.overlay_genomes,
            manifest.base_genomes
        );
    }

    Ok(entry)
}

/// Compact all overlays into the base index by rebuilding.
///
/// Collects all genome directories (base + overlays), rebuilds the full index,
/// and cleans up overlay directories.
pub fn compact(
    index_dir: &Path,
    all_genomes_dir: &Path,
    kmer_size: usize,
    threads: usize,
    max_ram_bytes: Option<usize>,
) -> Result<()> {
    log::info!("Compacting index: rebuilding from all genomes in {:?}", all_genomes_dir);

    let manifest = OverlayManifest::load_or_create(index_dir)?;
    log::info!(
        "Current state: {} base + {} overlay ({} overlays) = {} total genomes",
        manifest.base_genomes,
        manifest.overlay_genomes,
        manifest.overlay_count,
        manifest.total_genomes()
    );

    // Rebuild the full index
    crate::index::build_index_with_options(
        all_genomes_dir,
        index_dir,
        kmer_size,
        threads,
        max_ram_bytes,
    )?;

    // Clean up overlays
    let overlays_dir = index_dir.join("overlays");
    if overlays_dir.exists() {
        std::fs::remove_dir_all(&overlays_dir)?;
        log::info!("Removed overlay directories");
    }

    // Remove old manifest (the new index is a fresh base)
    let manifest_path = index_dir.join("overlay_manifest.json");
    if manifest_path.exists() {
        std::fs::remove_file(&manifest_path)?;
    }

    log::info!("Compaction complete: index rebuilt with all genomes");
    Ok(())
}

fn chrono_now() -> String {
    // Simple ISO 8601 timestamp without chrono dependency
    let duration = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = duration.as_secs();
    // Approximate date: good enough for a timestamp label
    format!("{}s_since_epoch", secs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_manifest_total_genomes() {
        let m = OverlayManifest {
            overlay_count: 2,
            overlay_genomes: 100,
            next_genome_offset: 1100,
            base_genomes: 1000,
            overlays: vec![],
        };
        assert_eq!(m.total_genomes(), 1100);
    }

    #[test]
    fn test_needs_compaction() {
        let m = OverlayManifest {
            overlay_count: 5,
            overlay_genomes: 150,
            next_genome_offset: 1150,
            base_genomes: 1000,
            overlays: vec![],
        };
        assert!(m.needs_compaction(0.1));
        assert!(!m.needs_compaction(0.2));
    }

    #[test]
    fn test_overlay_entry_serialize() {
        let entry = OverlayEntry {
            name: "overlay_001".into(),
            num_genomes: 50,
            genome_offset: 1000,
            created: "2026-04-08".into(),
        };
        let json = serde_json::to_string(&entry).unwrap();
        assert!(json.contains("overlay_001"));
        let parsed: OverlayEntry = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.num_genomes, 50);
    }
}
