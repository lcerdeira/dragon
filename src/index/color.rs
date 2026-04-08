/// Color index: maps each unitig to a set of genome IDs using Roaring Bitmaps.
/// Stored as a memory-mapped file for low-RAM on-demand access.

use anyhow::Result;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

const COLOR_INDEX_MAGIC: u32 = 0x4452_474E; // "DRGN"
const COLOR_INDEX_VERSION: u32 = 1;

/// Header for the on-disk color index.
#[derive(Serialize, Deserialize)]
struct ColorIndexHeader {
    magic: u32,
    version: u32,
    num_unitigs: u64,
    num_genomes: u64,
}

/// Runtime color index with memory-mapped access.
pub struct ColorIndex {
    data: memmap2::Mmap,
    num_unitigs: u64,
    num_genomes: u64,
    offsets_start: usize,
    data_start: usize,
}

impl ColorIndex {
    /// Look up which genomes contain a given unitig.
    pub fn get_colors(&self, unitig_id: u32) -> Result<RoaringBitmap> {
        let id = unitig_id as usize;
        if id >= self.num_unitigs as usize {
            return Ok(RoaringBitmap::new());
        }

        // Read offset for this unitig and next unitig
        let offset_pos = self.offsets_start + id * 8;
        let start = u64::from_le_bytes(
            self.data[offset_pos..offset_pos + 8].try_into().unwrap(),
        ) as usize;
        let end = u64::from_le_bytes(
            self.data[offset_pos + 8..offset_pos + 16].try_into().unwrap(),
        ) as usize;

        let bitmap_bytes = &self.data[self.data_start + start..self.data_start + end];
        let bitmap = RoaringBitmap::deserialize_from(bitmap_bytes)?;
        Ok(bitmap)
    }

    pub fn num_unitigs(&self) -> u64 {
        self.num_unitigs
    }

    pub fn num_genomes(&self) -> u64 {
        self.num_genomes
    }
}

/// Build the color index from a color file (TSV: unitig_id\tgenome_id1,genome_id2,...).
pub fn build_color_index(color_file: &Path, output_dir: &Path, num_genomes: usize) -> Result<()> {
    let file = std::fs::File::open(color_file)?;
    let reader = BufReader::new(file);

    // Parse color file into bitmaps
    let mut bitmaps: Vec<(u32, RoaringBitmap)> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        let unitig_id: u32 = parts[0].parse()?;
        let mut bitmap = RoaringBitmap::new();
        for genome_str in parts[1].split(',') {
            if let Ok(genome_id) = genome_str.trim().parse::<u32>() {
                bitmap.insert(genome_id);
            }
        }
        bitmaps.push((unitig_id, bitmap));
    }

    // Sort by unitig ID
    bitmaps.sort_by_key(|(id, _)| *id);

    // Fill gaps (unitigs without colors get empty bitmaps)
    let max_id = bitmaps.last().map(|(id, _)| *id).unwrap_or(0);
    let num_unitigs = (max_id + 1) as usize;
    let mut all_bitmaps: Vec<RoaringBitmap> = vec![RoaringBitmap::new(); num_unitigs];
    for (id, bitmap) in bitmaps {
        all_bitmaps[id as usize] = bitmap;
    }

    // Serialize to on-disk format
    let index_path = output_dir.join("colors.drgn");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&index_path)?);

    // Write header
    let header = ColorIndexHeader {
        magic: COLOR_INDEX_MAGIC,
        version: COLOR_INDEX_VERSION,
        num_unitigs: num_unitigs as u64,
        num_genomes: num_genomes as u64,
    };
    let header_bytes = bincode::serialize(&header)?;
    writer.write_all(&header_bytes)?;

    // Serialize all bitmaps to a buffer first to compute offsets
    let mut bitmap_data: Vec<Vec<u8>> = Vec::with_capacity(num_unitigs);
    for bitmap in &all_bitmaps {
        let mut buf = Vec::new();
        bitmap.serialize_into(&mut buf)?;
        bitmap_data.push(buf);
    }

    // Write offsets (num_unitigs + 1 entries for start/end pairs)
    let mut cumulative_offset = 0u64;
    for data in &bitmap_data {
        writer.write_all(&cumulative_offset.to_le_bytes())?;
        cumulative_offset += data.len() as u64;
    }
    writer.write_all(&cumulative_offset.to_le_bytes())?; // sentinel

    // Write bitmap data
    for data in &bitmap_data {
        writer.write_all(data)?;
    }

    log::info!(
        "Color index: {} unitigs, {} genomes, {:.1} MB on disk",
        num_unitigs,
        num_genomes,
        std::fs::metadata(&index_path)?.len() as f64 / 1_048_576.0
    );

    Ok(())
}

/// Load the color index via memory mapping.
pub fn load_color_index(index_dir: &Path) -> Result<ColorIndex> {
    let index_path = index_dir.join("colors.drgn");
    let mmap = crate::util::mmap::mmap_open(&index_path)?;

    // Parse header
    let header: ColorIndexHeader = bincode::deserialize(&mmap[..24])?;
    if header.magic != COLOR_INDEX_MAGIC {
        anyhow::bail!("Invalid color index magic number");
    }

    let header_size = 24; // bincode-serialized header size
    let offsets_start = header_size;
    let offsets_size = (header.num_unitigs as usize + 1) * 8;
    let data_start = offsets_start + offsets_size;

    Ok(ColorIndex {
        data: mmap,
        num_unitigs: header.num_unitigs,
        num_genomes: header.num_genomes,
        offsets_start,
        data_start,
    })
}
