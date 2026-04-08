/// Memory-mapped file helpers for low-RAM index access.

use anyhow::Result;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::path::Path;

/// Open a file as a read-only memory map.
pub fn mmap_open(path: &Path) -> Result<Mmap> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    Ok(mmap)
}

/// Read a bincode-serialized structure from a file.
pub fn read_bincode<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T> {
    let data = std::fs::read(path)?;
    let value: T = bincode::deserialize(&data)?;
    Ok(value)
}

/// Write a bincode-serialized structure to a file.
pub fn write_bincode<T: serde::Serialize>(path: &Path, value: &T) -> Result<()> {
    let data = bincode::serialize(value)?;
    std::fs::write(path, data)?;
    Ok(())
}
