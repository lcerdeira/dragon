/// Genome specificity index: identifies unitigs that are private or near-private
/// to individual genomes, enabling strain-level discrimination.
///
/// At build time, scans the color index to find unitigs with small color sets
/// (present in <= `max_sharing` genomes). These "private markers" carry
/// maximum discriminative power for within-species resolution.
///
/// At query time, counts how many of a genome's private unitigs are hit by seeds,
/// producing a specificity score that complements chain score for ranking.

use anyhow::Result;
use roaring::RoaringBitmap;
use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::index::color::ColorIndex;

const SPECIFICITY_MAGIC: u32 = 0x4452_5350; // "DRSP"
const SPECIFICITY_VERSION: u32 = 1;

/// Maximum color set size to consider a unitig as "private" to a genome.
/// Unitigs in <= this many genomes are treated as strain markers.
pub const DEFAULT_MAX_SHARING: u64 = 3;

/// Runtime specificity index: genome_id -> set of private unitig IDs.
pub struct SpecificityIndex {
    /// For each genome, the set of unitig IDs that are private/near-private to it.
    private_unitigs: HashMap<u32, RoaringBitmap>,
    /// Maximum sharing threshold used when building.
    pub max_sharing: u64,
    /// Total number of genomes in the database.
    pub num_genomes: u64,
}

impl SpecificityIndex {
    /// Build the specificity index by scanning the color index.
    ///
    /// A unitig is "private" to a genome if its color set size <= `max_sharing`.
    /// For each such unitig, it is added to every genome in its color set.
    pub fn build(color_index: &ColorIndex, max_sharing: u64) -> Result<Self> {
        let num_unitigs = color_index.num_unitigs();
        let num_genomes = color_index.num_genomes();
        let mut private_unitigs: HashMap<u32, RoaringBitmap> = HashMap::new();

        let mut total_private = 0u64;
        let mut total_exclusive = 0u64; // |color_set| == 1

        for unitig_id in 0..num_unitigs as u32 {
            if let Ok(colors) = color_index.get_colors(unitig_id) {
                let card = colors.len();
                if card >= 1 && card <= max_sharing {
                    total_private += 1;
                    if card == 1 {
                        total_exclusive += 1;
                    }
                    for genome_id in colors.iter() {
                        private_unitigs
                            .entry(genome_id)
                            .or_default()
                            .insert(unitig_id);
                    }
                }
            }
        }

        log::info!(
            "Specificity index: {} private unitigs (≤{} genomes), {} exclusive (=1 genome), across {} genomes",
            total_private, max_sharing, total_exclusive, private_unitigs.len()
        );

        Ok(Self {
            private_unitigs,
            max_sharing,
            num_genomes,
        })
    }

    /// Get the set of private unitigs for a genome. Returns None if genome has no private unitigs.
    pub fn get_private_unitigs(&self, genome_id: u32) -> Option<&RoaringBitmap> {
        self.private_unitigs.get(&genome_id)
    }

    /// Count how many of a genome's private unitigs appear in the given set of hit unitig IDs.
    ///
    /// Returns (private_hits, total_private) where:
    /// - private_hits: number of the genome's private unitigs that were hit
    /// - total_private: total number of private unitigs for this genome
    pub fn count_private_hits(&self, genome_id: u32, hit_unitigs: &RoaringBitmap) -> (u64, u64) {
        match self.private_unitigs.get(&genome_id) {
            Some(private) => {
                let hits = private.intersection_len(hit_unitigs);
                (hits, private.len())
            }
            None => (0, 0),
        }
    }

    /// Compute a specificity score for a genome given the query's hit unitigs.
    ///
    /// Score = private_hits * log2(num_genomes / max_sharing)
    ///
    /// This weights each private unitig hit by its information content.
    /// A genome with many private unitig hits gets a high specificity score.
    pub fn specificity_score(&self, genome_id: u32, hit_unitigs: &RoaringBitmap) -> f64 {
        let (hits, _total) = self.count_private_hits(genome_id, hit_unitigs);
        if hits == 0 {
            return 0.0;
        }
        // Each private hit carries information proportional to its exclusivity
        let ic_per_hit = if self.num_genomes > 0 && self.max_sharing > 0 {
            (self.num_genomes as f64 / self.max_sharing as f64).log2()
        } else {
            1.0
        };
        hits as f64 * ic_per_hit
    }

    /// Save the specificity index to disk.
    pub fn save(&self, index_dir: &Path) -> Result<()> {
        let path = index_dir.join("specificity.drgn");
        let mut writer = BufWriter::new(std::fs::File::create(&path)?);

        // Header
        writer.write_all(&SPECIFICITY_MAGIC.to_le_bytes())?;
        writer.write_all(&SPECIFICITY_VERSION.to_le_bytes())?;
        writer.write_all(&self.max_sharing.to_le_bytes())?;
        writer.write_all(&self.num_genomes.to_le_bytes())?;
        writer.write_all(&(self.private_unitigs.len() as u64).to_le_bytes())?;

        // Per-genome entries: genome_id (u32) + bitmap_len (u64) + bitmap_bytes
        for (&genome_id, bitmap) in &self.private_unitigs {
            writer.write_all(&genome_id.to_le_bytes())?;
            let mut buf = Vec::new();
            bitmap.serialize_into(&mut buf)?;
            writer.write_all(&(buf.len() as u64).to_le_bytes())?;
            writer.write_all(&buf)?;
        }

        let size = std::fs::metadata(&path)?.len();
        log::info!(
            "Specificity index saved: {:.1} MB",
            size as f64 / 1_048_576.0
        );

        Ok(())
    }

    /// Load the specificity index from disk.
    pub fn load(index_dir: &Path) -> Result<Self> {
        let path = index_dir.join("specificity.drgn");
        let file = std::fs::File::open(&path)?;
        let mut reader = BufReader::new(file);

        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf4)?;
        let magic = u32::from_le_bytes(buf4);
        if magic != SPECIFICITY_MAGIC {
            anyhow::bail!("Invalid specificity index magic number");
        }

        reader.read_exact(&mut buf4)?;
        let _version = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf8)?;
        let max_sharing = u64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let num_genomes = u64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let num_entries = u64::from_le_bytes(buf8);

        let mut private_unitigs = HashMap::with_capacity(num_entries as usize);

        for _ in 0..num_entries {
            reader.read_exact(&mut buf4)?;
            let genome_id = u32::from_le_bytes(buf4);

            reader.read_exact(&mut buf8)?;
            let bitmap_len = u64::from_le_bytes(buf8) as usize;

            let mut bitmap_bytes = vec![0u8; bitmap_len];
            reader.read_exact(&mut bitmap_bytes)?;
            let bitmap = RoaringBitmap::deserialize_from(&bitmap_bytes[..])?;

            private_unitigs.insert(genome_id, bitmap);
        }

        log::info!(
            "Loaded specificity index: {} genomes with private unitigs (max_sharing={})",
            private_unitigs.len(), max_sharing
        );

        Ok(Self {
            private_unitigs,
            max_sharing,
            num_genomes,
        })
    }

    /// Load from disk, or build from color index if not present.
    pub fn load_or_build(index_dir: &Path, color_index: &ColorIndex) -> Result<Self> {
        let path = index_dir.join("specificity.drgn");
        if path.exists() {
            Self::load(index_dir)
        } else {
            log::info!("Specificity index not found, building from color index...");
            let index = Self::build(color_index, DEFAULT_MAX_SHARING)?;
            index.save(index_dir)?;
            Ok(index)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_private_hits() {
        let mut private_unitigs = HashMap::new();
        let mut genome0 = RoaringBitmap::new();
        genome0.insert(10);
        genome0.insert(20);
        genome0.insert(30);
        private_unitigs.insert(0, genome0);

        let idx = SpecificityIndex {
            private_unitigs,
            max_sharing: 3,
            num_genomes: 1000,
        };

        // Query hits unitigs 10 and 30
        let mut hits = RoaringBitmap::new();
        hits.insert(10);
        hits.insert(30);
        hits.insert(99); // not private to genome 0

        let (private_hits, total) = idx.count_private_hits(0, &hits);
        assert_eq!(private_hits, 2);
        assert_eq!(total, 3);

        // Unknown genome
        let (private_hits, total) = idx.count_private_hits(999, &hits);
        assert_eq!(private_hits, 0);
        assert_eq!(total, 0);
    }

    #[test]
    fn test_specificity_score() {
        let mut private_unitigs = HashMap::new();
        let mut genome0 = RoaringBitmap::new();
        for i in 0..10 {
            genome0.insert(i);
        }
        private_unitigs.insert(0, genome0);

        let idx = SpecificityIndex {
            private_unitigs,
            max_sharing: 3,
            num_genomes: 10000,
        };

        let mut hits = RoaringBitmap::new();
        for i in 0..5 {
            hits.insert(i);
        }

        let score = idx.specificity_score(0, &hits);
        // 5 hits * log2(10000/3) ≈ 5 * 11.7 ≈ 58.5
        assert!(score > 50.0, "score = {}", score);

        // No hits
        let empty = RoaringBitmap::new();
        assert_eq!(idx.specificity_score(0, &empty), 0.0);
    }
}
