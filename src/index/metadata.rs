/// Index metadata: header, genome names, statistics.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::index::dbg::DbgResult;
use crate::index::unitig::UnitigSet;

#[derive(Serialize, Deserialize, Debug)]
pub struct IndexMetadata {
    pub version: String,
    pub kmer_size: usize,
    pub num_genomes: usize,
    pub num_unitigs: usize,
    pub total_unitig_bases: u64,
    pub genome_names: Vec<String>,
}

/// Write index metadata to disk.
pub fn write_metadata(
    output_dir: &Path,
    dbg_result: &DbgResult,
    unitigs: &UnitigSet,
) -> Result<()> {
    let metadata = IndexMetadata {
        version: env!("CARGO_PKG_VERSION").to_string(),
        kmer_size: dbg_result.kmer_size,
        num_genomes: dbg_result.num_genomes,
        num_unitigs: dbg_result.num_unitigs,
        total_unitig_bases: unitigs.total_bases(),
        genome_names: Vec::new(), // Populated during path building
    };

    let meta_path = output_dir.join("metadata.json");
    let json = serde_json::to_string_pretty(&metadata)?;
    std::fs::write(&meta_path, json)?;

    log::info!("Metadata written to {:?}", meta_path);
    Ok(())
}

/// Load index metadata.
pub fn load_metadata(index_dir: &Path) -> Result<IndexMetadata> {
    let meta_path = index_dir.join("metadata.json");
    let json = std::fs::read_to_string(&meta_path)?;
    let metadata: IndexMetadata = serde_json::from_str(&json)?;
    Ok(metadata)
}
