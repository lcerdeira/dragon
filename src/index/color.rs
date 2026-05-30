/// Color index: maps each unitig to a set of genome IDs using Roaring Bitmaps.
///
/// Two on-disk formats:
///
/// **v1** (`COLOR_INDEX_MAGIC_V1`): one serialised RoaringBitmap per unitig.
///   Layout: [header 24B] [offsets (num_unitigs+1)×8B] [bitmap data]
///   Problem: no deduplication — unitigs with identical color sets are stored
///   separately, inflating the file by num_unitigs / num_unique_sets (often 10-50×).
///
/// **v2** (`COLOR_INDEX_MAGIC_V2`): shared color-set pool + per-unitig 4-byte index.
///   Layout: [header 32B] [index num_unitigs×4B] [pool_offsets (num_sets+1)×4B] [pool data]
///   Pool stores each distinct RoaringBitmap exactly once; unitigs reference the pool
///   by a u32 set_id.  For a 26K-genome S. aureus shard this typically reduces
///   colors.drgn from ~7 GB to ~200-500 MB (10-50× reduction).
///
/// `load_color_index` detects the format automatically by magic number.

use std::collections::HashMap;

use anyhow::{bail, Context, Result};
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};
use std::io::{BufRead, Write};
use std::path::Path;

const COLOR_INDEX_MAGIC_V1: u32 = 0x4452_474E; // "DRGN"
const COLOR_INDEX_MAGIC_V2: u32 = 0x4452_474E; // same 4 bytes; version field differentiates
const COLOR_INDEX_VERSION_V1: u32 = 1;
const COLOR_INDEX_VERSION_V2: u32 = 2;

// ─── On-disk headers ────────────────────────────────────────────────────────

/// v1 header (bincode, 24 bytes).
#[derive(Serialize, Deserialize)]
struct ColorIndexHeaderV1 {
    magic: u32,
    version: u32,
    num_unitigs: u64,
    num_genomes: u64,
}

/// v2 header (bincode, 32 bytes).
#[derive(Serialize, Deserialize)]
struct ColorIndexHeaderV2 {
    magic: u32,
    version: u32,
    num_unitigs: u64,
    num_genomes: u64,
    num_sets: u32,
    _pad: u32,
}

// ─── Runtime index ──────────────────────────────────────────────────────────

enum ColorIndexInner {
    V1 {
        num_unitigs: u64,
        num_genomes: u64,
        offsets_start: usize,
        data_start: usize,
    },
    V2 {
        num_unitigs: u64,
        num_genomes: u64,
        num_sets: u32,
        /// Byte offset inside `data` where the per-unitig set_id array starts.
        index_start: usize,
        /// Byte offset where the (num_sets+1) u32 pool-offset entries start.
        pool_offsets_start: usize,
        /// Byte offset where the concatenated bitmap bytes start.
        pool_data_start: usize,
    },
}

/// Runtime color index backed by a memory-mapped file.
pub struct ColorIndex {
    data: memmap2::Mmap,
    inner: ColorIndexInner,
}

impl ColorIndex {
    /// Return the genome IDs that contain `unitig_id`.
    pub fn get_colors(&self, unitig_id: u32) -> Result<RoaringBitmap> {
        match &self.inner {
            ColorIndexInner::V1 {
                num_unitigs,
                offsets_start,
                data_start,
                ..
            } => {
                let id = unitig_id as usize;
                if id >= *num_unitigs as usize {
                    return Ok(RoaringBitmap::new());
                }
                let offset_pos = offsets_start + id * 8;
                let start = u64::from_le_bytes(
                    self.data[offset_pos..offset_pos + 8].try_into().unwrap(),
                ) as usize;
                let end = u64::from_le_bytes(
                    self.data[offset_pos + 8..offset_pos + 16].try_into().unwrap(),
                ) as usize;
                let bm = RoaringBitmap::deserialize_from(
                    &self.data[*data_start + start..*data_start + end],
                )?;
                Ok(bm)
            }
            ColorIndexInner::V2 {
                num_unitigs,
                num_sets,
                index_start,
                pool_offsets_start,
                pool_data_start,
                ..
            } => {
                let id = unitig_id as usize;
                if id >= *num_unitigs as usize {
                    return Ok(RoaringBitmap::new());
                }
                // Read set_id
                let ip = index_start + id * 4;
                let set_id =
                    u32::from_le_bytes(self.data[ip..ip + 4].try_into().unwrap()) as usize;
                if set_id >= *num_sets as usize {
                    return Ok(RoaringBitmap::new());
                }
                // Read pool offsets
                let op = pool_offsets_start + set_id * 4;
                let start =
                    u32::from_le_bytes(self.data[op..op + 4].try_into().unwrap()) as usize;
                let end =
                    u32::from_le_bytes(self.data[op + 4..op + 8].try_into().unwrap()) as usize;
                let bm = RoaringBitmap::deserialize_from(
                    &self.data[*pool_data_start + start..*pool_data_start + end],
                )?;
                Ok(bm)
            }
        }
    }

    pub fn num_unitigs(&self) -> u64 {
        match &self.inner {
            ColorIndexInner::V1 { num_unitigs, .. } => *num_unitigs,
            ColorIndexInner::V2 { num_unitigs, .. } => *num_unitigs,
        }
    }

    pub fn num_genomes(&self) -> u64 {
        match &self.inner {
            ColorIndexInner::V1 { num_genomes, .. } => *num_genomes,
            ColorIndexInner::V2 { num_genomes, .. } => *num_genomes,
        }
    }

    /// True if this index uses the v2 deduplicated format.
    pub fn is_v2(&self) -> bool {
        matches!(&self.inner, ColorIndexInner::V2 { .. })
    }
}

// ─── Building v1 ────────────────────────────────────────────────────────────

/// Build a v1 color index from a color TSV file.
///
/// TSV format: `unitig_id\tgenome_id1,genome_id2,...`
pub fn build_color_index(color_file: &Path, output_dir: &Path, num_genomes: usize) -> Result<()> {
    let file = std::fs::File::open(color_file)?;
    let reader = std::io::BufReader::new(file);

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

    bitmaps.sort_by_key(|(id, _)| *id);
    let max_id = bitmaps.last().map(|(id, _)| *id).unwrap_or(0);
    let num_unitigs = (max_id + 1) as usize;
    let mut all_bitmaps: Vec<RoaringBitmap> = vec![RoaringBitmap::new(); num_unitigs];
    for (id, bitmap) in bitmaps {
        all_bitmaps[id as usize] = bitmap;
    }

    write_color_index_v1(&all_bitmaps, output_dir, num_genomes)
}

fn write_color_index_v1(
    bitmaps: &[RoaringBitmap],
    output_dir: &Path,
    num_genomes: usize,
) -> Result<()> {
    let num_unitigs = bitmaps.len();
    let index_path = output_dir.join("colors.drgn");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&index_path)?);

    let header = ColorIndexHeaderV1 {
        magic: COLOR_INDEX_MAGIC_V1,
        version: COLOR_INDEX_VERSION_V1,
        num_unitigs: num_unitigs as u64,
        num_genomes: num_genomes as u64,
    };
    let header_bytes = bincode::serialize(&header)?;
    writer.write_all(&header_bytes)?;

    let mut bitmap_data: Vec<Vec<u8>> = Vec::with_capacity(num_unitigs);
    for bm in bitmaps {
        let mut buf = Vec::new();
        bm.serialize_into(&mut buf)?;
        bitmap_data.push(buf);
    }

    let mut cumulative_offset = 0u64;
    for data in &bitmap_data {
        writer.write_all(&cumulative_offset.to_le_bytes())?;
        cumulative_offset += data.len() as u64;
    }
    writer.write_all(&cumulative_offset.to_le_bytes())?;

    for data in &bitmap_data {
        writer.write_all(data)?;
    }

    log::info!(
        "Color index v1: {} unitigs, {} genomes, {:.1} MB on disk",
        num_unitigs,
        num_genomes,
        std::fs::metadata(&index_path)?.len() as f64 / 1_048_576.0
    );
    Ok(())
}

// ─── Building v2 ────────────────────────────────────────────────────────────

/// Write a v2 (deduplicated) color index directly from a slice of bitmaps.
///
/// Deduplicates bitmaps by serialisation hash, stores each unique bitmap once,
/// writes a per-unitig u32 set_id index, and a u32 pool-offset table.
pub fn write_color_index_v2(
    bitmaps: &[RoaringBitmap],
    output_dir: &Path,
    num_genomes: usize,
) -> Result<ColorV2Stats> {
    let num_unitigs = bitmaps.len();

    // --- Deduplication pass ---
    // Key: serialised bitmap bytes (exact hash).  Value: u32 pool index.
    let mut pool_map: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut pool: Vec<Vec<u8>> = Vec::new();
    let mut unitig_set_ids: Vec<u32> = Vec::with_capacity(num_unitigs);

    for bm in bitmaps {
        let mut buf = Vec::new();
        bm.serialize_into(&mut buf)?;
        let next_id = pool.len() as u32;
        let set_id = *pool_map.entry(buf.clone()).or_insert_with(|| {
            pool.push(buf);
            next_id
        });
        unitig_set_ids.push(set_id);
    }

    let num_sets = pool.len() as u32;

    // Total pool data size must fit in u32 (< 4 GB).
    let pool_total: u64 = pool.iter().map(|b| b.len() as u64).sum();
    if pool_total > u32::MAX as u64 {
        anyhow::bail!(
            "colors.drgn v2: pool data {:.1} GB exceeds 4 GB u32 offset limit; \
             use v1 for this index",
            pool_total as f64 / 1e9
        );
    }

    // --- Serialize ---
    let index_path = output_dir.join("colors.drgn");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&index_path)?);

    // Header (32 bytes)
    let header = ColorIndexHeaderV2 {
        magic: COLOR_INDEX_MAGIC_V2,
        version: COLOR_INDEX_VERSION_V2,
        num_unitigs: num_unitigs as u64,
        num_genomes: num_genomes as u64,
        num_sets,
        _pad: 0,
    };
    let header_bytes = bincode::serialize(&header)?;
    debug_assert_eq!(header_bytes.len(), 32, "v2 header must be 32 bytes");
    writer.write_all(&header_bytes)?;

    // Per-unitig index: num_unitigs × u32
    for &sid in &unitig_set_ids {
        writer.write_all(&sid.to_le_bytes())?;
    }

    // Pool offset table: (num_sets + 1) × u32
    let mut cum = 0u32;
    for bytes in &pool {
        writer.write_all(&cum.to_le_bytes())?;
        cum += bytes.len() as u32;
    }
    writer.write_all(&cum.to_le_bytes())?; // sentinel

    // Pool data
    for bytes in &pool {
        writer.write_all(bytes)?;
    }

    let on_disk = std::fs::metadata(&index_path)?.len();
    let stats = ColorV2Stats {
        num_unitigs: num_unitigs as u64,
        num_genomes: num_genomes as u64,
        num_sets,
        on_disk_bytes: on_disk,
    };
    log::info!(
        "Color index v2: {} unitigs, {} unique sets ({:.1}× dedup), \
         {} genomes, {:.1} MB on disk",
        num_unitigs,
        num_sets,
        num_unitigs as f64 / num_sets.max(1) as f64,
        num_genomes,
        on_disk as f64 / 1_048_576.0
    );
    Ok(stats)
}

/// Statistics returned by v2 writes.
pub struct ColorV2Stats {
    pub num_unitigs: u64,
    pub num_genomes: u64,
    pub num_sets: u32,
    pub on_disk_bytes: u64,
}

// ─── Loading ─────────────────────────────────────────────────────────────────

/// Load the color index from `index_dir/colors.drgn`.
///
/// Detects v1 vs v2 automatically by the `version` field in the header.
pub fn load_color_index(index_dir: &Path) -> Result<ColorIndex> {
    let index_path = index_dir.join("colors.drgn");
    let mmap = crate::util::mmap::mmap_open(&index_path)?;

    // Read the first 8 bytes to detect magic + version.
    if mmap.len() < 8 {
        bail!("colors.drgn too small to contain a valid header");
    }
    let magic = u32::from_le_bytes(mmap[0..4].try_into().unwrap());
    let version = u32::from_le_bytes(mmap[4..8].try_into().unwrap());

    // Both v1 and v2 share the same magic value; the version field
    // differentiates them (bincode encodes u32 fields at the same offsets).
    match (magic, version) {
        (m, 1) if m == COLOR_INDEX_MAGIC_V1 => load_v1(mmap),
        (m, 2) if m == COLOR_INDEX_MAGIC_V2 => load_v2(mmap),
        _ => bail!(
            "colors.drgn: unrecognised magic {:#010x} / version {}",
            magic,
            version
        ),
    }
}

fn load_v1(mmap: memmap2::Mmap) -> Result<ColorIndex> {
    let header: ColorIndexHeaderV1 = bincode::deserialize(&mmap[..24])?;
    if header.magic != COLOR_INDEX_MAGIC_V1 {
        bail!("colors.drgn v1: invalid magic {:#010x}", header.magic);
    }
    let header_size = 24usize;
    let offsets_start = header_size;
    let offsets_size = (header.num_unitigs as usize + 1) * 8;
    let data_start = offsets_start + offsets_size;
    Ok(ColorIndex {
        data: mmap,
        inner: ColorIndexInner::V1 {
            num_unitigs: header.num_unitigs,
            num_genomes: header.num_genomes,
            offsets_start,
            data_start,
        },
    })
}

fn load_v2(mmap: memmap2::Mmap) -> Result<ColorIndex> {
    if mmap.len() < 32 {
        bail!("colors.drgn v2: file too small for header");
    }
    let header: ColorIndexHeaderV2 = bincode::deserialize(&mmap[..32])?;
    let header_size = 32usize;
    let num_unitigs = header.num_unitigs as usize;
    let num_sets = header.num_sets as usize;

    let index_start = header_size;
    let index_size = num_unitigs * 4;

    let pool_offsets_start = index_start + index_size;
    let pool_offsets_size = (num_sets + 1) * 4;

    let pool_data_start = pool_offsets_start + pool_offsets_size;

    if mmap.len() < pool_data_start {
        bail!(
            "colors.drgn v2: file too small (need {} bytes for headers+index+pool_offsets, \
             got {})",
            pool_data_start,
            mmap.len()
        );
    }

    Ok(ColorIndex {
        data: mmap,
        inner: ColorIndexInner::V2 {
            num_unitigs: header.num_unitigs,
            num_genomes: header.num_genomes,
            num_sets: header.num_sets,
            index_start,
            pool_offsets_start,
            pool_data_start,
        },
    })
}

// ─── GGCAT-direct v2 build (no intermediate TSV) ─────────────────────────────

/// Build a v2 color index directly from GGCAT's colormap + unitigs.fa.
///
/// Reads the GGCAT subset pool, then assigns each unitig its set_id directly
/// from the `C:hex_id` tag in `unitigs.fa` — no per-unitig bitmap expansion.
/// This is the most efficient path for new index builds.
pub fn build_color_drgn_v2_from_ggcat(
    colormap_path: &Path,
    unitigs_path: &Path,
    output_dir: &Path,
) -> Result<()> {
    use crate::index::ggcat_colors;

    log::info!("Building colors.drgn v2 from GGCAT colormap {:?}", colormap_path);

    // Step 1: Build the pool from GGCAT subsets (each subset → RoaringBitmap).
    let (header, subsets_map) = ggcat_colors::parse_all(colormap_path)?;
    let num_genomes = header.colors_count as usize;
    let subsets_count = header.subsets_count as usize;
    log::info!(
        "  Colormap: {} genomes, {} color subsets",
        num_genomes,
        subsets_count
    );

    // Build pool: hex_id → pool_index, pool[pool_index] = serialised RoaringBitmap.
    // Pool is already deduplicated by GGCAT's own subset deduplication.
    let mut hex_to_pool: HashMap<u32, u32> = HashMap::with_capacity(subsets_count);
    let mut pool: Vec<Vec<u8>> = Vec::with_capacity(subsets_count);

    // Insert empty set at index 0 (unitigs with no C: tag map here).
    {
        let mut empty_buf = Vec::new();
        RoaringBitmap::new().serialize_into(&mut empty_buf)?;
        pool.push(empty_buf);
    }
    let empty_set_id: u32 = 0;

    for (&hex_id, genome_ids) in &subsets_map {
        let mut bm = RoaringBitmap::new();
        for &g in genome_ids {
            bm.insert(g);
        }
        let mut buf = Vec::new();
        bm.serialize_into(&mut buf)?;
        let idx = pool.len() as u32;
        hex_to_pool.insert(hex_id, idx);
        pool.push(buf);
    }

    // Step 2: Stream unitigs.fa to build the per-unitig set_id array.
    let unitig_reader = std::io::BufReader::new(std::fs::File::open(unitigs_path)?);
    let mut unitig_set_ids: Vec<u32> = Vec::new();
    let mut unitig_count = 0u64;

    for line in unitig_reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            continue;
        }
        let mut set_id = empty_set_id;
        for part in line[1..].split_whitespace().skip(1) {
            if let Some(rest) = part.strip_prefix("C:") {
                if let Some(hex_id_str) = rest.split(':').next() {
                    if let Ok(hex_id) = u32::from_str_radix(hex_id_str, 16) {
                        if let Some(&pid) = hex_to_pool.get(&hex_id) {
                            set_id = pid;
                        }
                    }
                }
                break;
            }
        }
        unitig_set_ids.push(set_id);
        unitig_count += 1;
        if unitig_count % 1_000_000 == 0 {
            log::info!("  Processed {} unitigs...", unitig_count);
        }
    }

    log::info!(
        "  {} unitigs → {} pool entries ({:.1}× GGCAT dedup factor)",
        unitig_count,
        pool.len(),
        unitig_count as f64 / pool.len().max(1) as f64
    );

    write_color_index_v2_raw(&unitig_set_ids, &pool, output_dir, num_genomes)
}

/// Build colors.drgn v2 from a bitmaps slice (used during new index builds
/// that already have per-unitig bitmaps in memory).
pub fn build_color_drgn_direct(
    colormap_path: &Path,
    unitigs_path: &Path,
    output_dir: &Path,
) -> Result<crate::index::ggcat_colors::GgcatColorHeader> {
    // Delegate to v2 builder for new builds (gets GGCAT-level dedup for free).
    build_color_drgn_v2_from_ggcat(colormap_path, unitigs_path, output_dir)?;
    crate::index::ggcat_colors::read_header(colormap_path)
}

// ─── Migration: v1 → v2 (with optional high-cardinality truncation) ──────────

/// Migrate an existing `colors.drgn` v1 file to v2 in-place.
///
/// # High-cardinality truncation (`max_cardinality`)
///
/// In a large, diverse genome database (e.g. 26 K S. aureus strains) most
/// unitigs have nearly unique color sets, so deduplication saves almost nothing.
/// The dominant disk cost is "core genome" unitigs that appear in most or all
/// genomes: each stores an 8 KB RoaringBitmap bitset container.
///
/// These high-cardinality unitigs are **useless for containment ranking**:
/// `containment_rank` sorts unitigs by cardinality (ascending) and stops
/// adding to the candidate set once `CANDIDATE_TARGET` (2 000) is reached.
/// A unitig in >2 000 genomes sits at the back of the list and is never used.
///
/// Passing `max_cardinality = Some(N)` replaces any bitmap with cardinality
/// > N with an empty set before writing v2.  All truncated unitigs share a
/// single empty pool entry (effective deduplication), and the surviving
/// low-cardinality bitmaps are small enough that the pool fits in u32 offsets.
///
/// Recommended value: `Some(2_000)` (matches `CANDIDATE_TARGET` in
/// `containment.rs`).  This trades perfect recall on queries whose every
/// k-mer is in a high-cardinality unitig (very rare for S. aureus) for a
/// 3–6× disk and RAM reduction.
///
/// # Memory
///
/// Bitmaps are streamed one at a time from the mmap — only the pool and the
/// per-unitig set-id array are kept in RAM.  With `max_cardinality = 2_000`
/// the pool is typically <1 GB even for 26 K-genome shards.
pub fn migrate_colors_v1_to_v2(
    index_dir: &Path,
    max_cardinality: Option<u32>,
) -> Result<MigrateColorsStats> {
    let path = index_dir.join("colors.drgn");
    let tmp  = index_dir.join("colors.drgn.v2.tmp");
    let bak  = index_dir.join("colors.drgn.v1.bak");

    let v1_size;
    let num_sets;
    {
        let v1 = load_color_index(index_dir)?;
        if v1.is_v2() {
            log::info!("[skip] {:?}: colors.drgn already v2", index_dir);
            return Ok(MigrateColorsStats { skipped: true, ..Default::default() });
        }

        let num_unitigs = v1.num_unitigs() as usize;
        let num_genomes = v1.num_genomes();
        v1_size = std::fs::metadata(&path)?.len();

        // Safety guard: truncating at a threshold below half the shard's genome
        // count removes core-genome unitigs that containment_rank relies on to
        // build its initial candidate set.  For a shard with N genomes and
        // CANDIDATE_TARGET K (default 2000), any unitig with cardinality > K is
        // still needed when the query is from a core-genome region — it pushes
        // the candidate set to K in one step. Only truncate when max_cardinality
        // is well above N/2 so the core genome is preserved.
        //
        // Rule: refuse if max_cardinality < num_genomes * 0.90.
        if let Some(k) = max_cardinality {
            let safe_min = (num_genomes as f64 * 0.90) as u64;
            if (k as u64) < safe_min {
                bail!(
                    "colors.drgn truncation unsafe for {:?}: \
                     --max-cardinality {} < 90% of shard genome count {} ({} min). \
                     Core-genome queries would lose all candidates. \
                     Use --max-cardinality >= {} or omit the flag.",
                    index_dir, k, num_genomes, safe_min, safe_min
                );
            }
        }

        log::info!(
            "[migrate] {:?}: {} unitigs, {} genomes{}, streaming ...",
            index_dir,
            num_unitigs,
            num_genomes,
            match max_cardinality {
                Some(k) => format!(", truncating cardinality > {}", k),
                None    => String::new(),
            }
        );

        // --- Streaming dedup pass ---
        // Process one bitmap at a time: deserialise from mmap, optionally
        // truncate, serialise, insert into pool.  Only pool + set_ids in RAM.
        let mut pool_map: HashMap<Vec<u8>, u32> = HashMap::new();
        let mut pool: Vec<Vec<u8>> = Vec::new();
        let mut unitig_set_ids: Vec<u32> = Vec::with_capacity(num_unitigs);

        // Pre-serialise an empty bitmap.  Only inserted into the pool if
        // truncation is active (avoids a spurious pool entry when max_cardinality
        // is None and no bitmap is naturally empty).
        let empty_buf = {
            let mut b = Vec::new();
            RoaringBitmap::new().serialize_into(&mut b)?;
            b
        };
        // If truncation is active, seed pool[0] = empty so all truncated
        // unitigs share a single lookup without repeated HashMap insertions.
        if max_cardinality.is_some() {
            pool_map.insert(empty_buf.clone(), 0u32);
            pool.push(empty_buf.clone());
        }

        let mut truncated = 0u64;

        for uid in 0..num_unitigs as u32 {
            let bm = v1.get_colors(uid).unwrap_or_default();

            // Truncate high-cardinality bitmaps to empty.
            let buf = if max_cardinality.is_some_and(|k| bm.len() > k as u64) {
                truncated += 1;
                empty_buf.clone()
            } else {
                let mut b = Vec::new();
                bm.serialize_into(&mut b)?;
                b
            };

            let next_id = pool.len() as u32;
            let sid = *pool_map.entry(buf.clone()).or_insert_with(|| {
                pool.push(buf);
                next_id
            });
            unitig_set_ids.push(sid);

            if (uid + 1) % 500_000 == 0 {
                log::info!(
                    "  {} / {} unitigs processed ({} pool entries so far) ...",
                    uid + 1, num_unitigs, pool.len()
                );
            }
        }
        num_sets = pool.len() as u32;

        if truncated > 0 {
            log::info!(
                "  truncated {}/{} unitigs (cardinality > {})",
                truncated, num_unitigs,
                max_cardinality.unwrap_or(0)
            );
        }

        // Check pool fits in u32 offsets; bail clearly if not.
        let pool_total: u64 = pool.iter().map(|b| b.len() as u64).sum();
        if pool_total > u32::MAX as u64 {
            bail!(
                "colors.drgn v2: pool {:.1} GB still exceeds 4 GB u32 limit after \
                 processing {:?}. \
                 Try a lower --max-cardinality value (currently {:?}).",
                pool_total as f64 / 1e9,
                index_dir,
                max_cardinality,
            );
        }

        write_color_index_v2_to_path(&unitig_set_ids, &pool, &tmp, num_genomes as usize)?;
    } // v1 mmap dropped

    let v2_size = std::fs::metadata(&tmp)?.len();
    log::info!(
        "  {:.2} GB → {:.2} GB ({:.2}×), {} unique sets",
        v1_size as f64 / 1e9,
        v2_size as f64 / 1e9,
        v1_size as f64 / v2_size.max(1) as f64,
        num_sets,
    );

    // Atomic install.
    std::fs::rename(&path, &bak)
        .with_context(|| format!("rename {:?} → {:?}", path, bak))?;
    std::fs::rename(&tmp, &path)
        .with_context(|| format!("rename {:?} → {:?}", tmp, path))?;
    log::info!("  installed v2; original kept as {:?}", bak);

    let installed = load_color_index(index_dir)?;
    Ok(MigrateColorsStats {
        skipped: false,
        num_unitigs: installed.num_unitigs(),
        num_genomes: installed.num_genomes(),
        num_sets,
        v1_bytes: v1_size,
        v2_bytes: v2_size,
    })
}

/// Statistics returned by `migrate_colors_v1_to_v2`.
#[derive(Default, Debug)]
pub struct MigrateColorsStats {
    pub skipped: bool,
    pub num_unitigs: u64,
    pub num_genomes: u64,
    pub num_sets: u32,
    pub v1_bytes: u64,
    pub v2_bytes: u64,
}

// ─── Internal helpers ────────────────────────────────────────────────────────

/// Write a v2 file to `output_dir/colors.drgn`.
fn write_color_index_v2_raw(
    unitig_set_ids: &[u32],
    pool: &[Vec<u8>],
    output_dir: &Path,
    num_genomes: usize,
) -> Result<()> {
    let out_path = output_dir.join("colors.drgn");
    write_color_index_v2_to_path(unitig_set_ids, pool, &out_path, num_genomes)
}

/// Write a v2 file to an explicit output path (used by migration to avoid
/// overwriting the source file).
fn write_color_index_v2_to_path(
    unitig_set_ids: &[u32],
    pool: &[Vec<u8>],
    out_path: &Path,
    num_genomes: usize,
) -> Result<()> {
    let num_unitigs = unitig_set_ids.len();
    let num_sets = pool.len() as u32;

    let pool_total: u64 = pool.iter().map(|b| b.len() as u64).sum();
    if pool_total > u32::MAX as u64 {
        anyhow::bail!(
            "colors.drgn v2: pool {:.1} GB exceeds u32 offset limit",
            pool_total as f64 / 1e9
        );
    }

    let mut writer = std::io::BufWriter::new(std::fs::File::create(out_path)?);

    let header = ColorIndexHeaderV2 {
        magic: COLOR_INDEX_MAGIC_V2,
        version: COLOR_INDEX_VERSION_V2,
        num_unitigs: num_unitigs as u64,
        num_genomes: num_genomes as u64,
        num_sets,
        _pad: 0,
    };
    let hdr_bytes = bincode::serialize(&header)?;
    debug_assert_eq!(hdr_bytes.len(), 32);
    writer.write_all(&hdr_bytes)?;

    for &sid in unitig_set_ids {
        writer.write_all(&sid.to_le_bytes())?;
    }

    let mut cum = 0u32;
    for bytes in pool {
        writer.write_all(&cum.to_le_bytes())?;
        cum += bytes.len() as u32;
    }
    writer.write_all(&cum.to_le_bytes())?;

    for bytes in pool {
        writer.write_all(bytes)?;
    }

    let on_disk = std::fs::metadata(out_path)?.len();
    log::info!(
        "colors.drgn v2: {} unitigs, {} unique sets ({:.1}× dedup), {:.1} MB",
        num_unitigs,
        num_sets,
        num_unitigs as f64 / num_sets.max(1) as f64,
        on_disk as f64 / 1_048_576.0
    );
    Ok(())
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn make_bitmap(ids: &[u32]) -> RoaringBitmap {
        let mut bm = RoaringBitmap::new();
        for &id in ids {
            bm.insert(id);
        }
        bm
    }

    #[test]
    fn v2_roundtrip_basic() {
        let dir = TempDir::new().unwrap();
        let bitmaps = vec![
            make_bitmap(&[0, 1, 2]),
            make_bitmap(&[0, 1, 2]), // duplicate
            make_bitmap(&[3, 4]),
            make_bitmap(&[0, 1, 2]), // duplicate again
            make_bitmap(&[]),        // empty
        ];
        let stats = write_color_index_v2(&bitmaps, dir.path(), 5).unwrap();
        // 3 unique sets: {0,1,2}, {3,4}, {}
        assert_eq!(stats.num_sets, 3);

        let idx = load_color_index(dir.path()).unwrap();
        assert!(idx.is_v2());
        assert_eq!(idx.num_unitigs(), 5);
        assert_eq!(idx.num_genomes(), 5);

        assert_eq!(idx.get_colors(0).unwrap(), make_bitmap(&[0, 1, 2]));
        assert_eq!(idx.get_colors(1).unwrap(), make_bitmap(&[0, 1, 2]));
        assert_eq!(idx.get_colors(2).unwrap(), make_bitmap(&[3, 4]));
        assert_eq!(idx.get_colors(3).unwrap(), make_bitmap(&[0, 1, 2]));
        assert_eq!(idx.get_colors(4).unwrap(), make_bitmap(&[]));
    }

    #[test]
    fn v2_out_of_range_returns_empty() {
        let dir = TempDir::new().unwrap();
        let bitmaps = vec![make_bitmap(&[1, 2])];
        write_color_index_v2(&bitmaps, dir.path(), 3).unwrap();
        let idx = load_color_index(dir.path()).unwrap();
        assert_eq!(idx.get_colors(999).unwrap(), RoaringBitmap::new());
    }

    #[test]
    fn v1_still_loads() {
        let dir = TempDir::new().unwrap();
        let bitmaps = vec![make_bitmap(&[0, 5]), make_bitmap(&[1])];
        write_color_index_v1(&bitmaps, dir.path(), 6).unwrap();
        let idx = load_color_index(dir.path()).unwrap();
        assert!(!idx.is_v2());
        assert_eq!(idx.get_colors(0).unwrap(), make_bitmap(&[0, 5]));
        assert_eq!(idx.get_colors(1).unwrap(), make_bitmap(&[1]));
    }

    #[test]
    fn migrate_v1_to_v2() {
        let dir = TempDir::new().unwrap();
        let bitmaps = vec![
            make_bitmap(&[0, 1]),
            make_bitmap(&[0, 1]), // dup
            make_bitmap(&[2]),
        ];
        write_color_index_v1(&bitmaps, dir.path(), 3).unwrap();

        let stats = migrate_colors_v1_to_v2(dir.path(), None).unwrap();
        assert!(!stats.skipped);
        assert_eq!(stats.num_sets, 2); // {0,1} and {2}
        assert!(stats.v2_bytes < stats.v1_bytes);

        let idx = load_color_index(dir.path()).unwrap();
        assert!(idx.is_v2());
        assert_eq!(idx.get_colors(0).unwrap(), make_bitmap(&[0, 1]));
        assert_eq!(idx.get_colors(1).unwrap(), make_bitmap(&[0, 1]));
        assert_eq!(idx.get_colors(2).unwrap(), make_bitmap(&[2]));

        // bak file must exist
        assert!(dir.path().join("colors.drgn.v1.bak").exists());
    }

    #[test]
    fn migrate_v1_to_v2_truncation_safety_guard() {
        // A shard with 10 genomes: threshold must be >= 9 (90%).
        // Passing threshold=5 (50%) should be rejected.
        let dir = TempDir::new().unwrap();
        let bitmaps: Vec<RoaringBitmap> = (0..5u32)
            .map(|i| make_bitmap(&(0..10u32).collect::<Vec<_>>()[..i as usize]))
            .collect();
        write_color_index_v1(&bitmaps, dir.path(), 10).unwrap();
        let result = migrate_colors_v1_to_v2(dir.path(), Some(5));
        assert!(result.is_err(), "should reject threshold < 90% of num_genomes");
        let msg = format!("{:#}", result.unwrap_err());
        assert!(msg.contains("unsafe"), "error should mention 'unsafe': {}", msg);
    }

    #[test]
    fn migrate_v2_skips_already_v2() {
        let dir = TempDir::new().unwrap();
        let bitmaps = vec![make_bitmap(&[0])];
        write_color_index_v2(&bitmaps, dir.path(), 1).unwrap();
        let stats = migrate_colors_v1_to_v2(dir.path(), None).unwrap();
        assert!(stats.skipped);
    }

    #[test]
    fn migrate_v1_to_v2_with_truncation() {
        // 3 unitigs: [0,1] (card=2), [0,1,2,3,4] (card=5, high), [3] (card=1)
        // num_genomes=5, safe_min = floor(5 * 0.90) = 4.
        // threshold=4 passes the safety guard (4 >= 4) and truncates card > 4.
        let dir = TempDir::new().unwrap();
        let bitmaps = vec![
            make_bitmap(&[0, 1]),
            make_bitmap(&[0, 1, 2, 3, 4]), // cardinality 5 > 4 → truncated
            make_bitmap(&[3]),
        ];
        write_color_index_v1(&bitmaps, dir.path(), 5).unwrap();

        let stats = migrate_colors_v1_to_v2(dir.path(), Some(4)).unwrap();
        assert!(!stats.skipped);

        let idx = load_color_index(dir.path()).unwrap();
        assert!(idx.is_v2());
        assert_eq!(idx.get_colors(0).unwrap(), make_bitmap(&[0, 1]));
        assert_eq!(idx.get_colors(1).unwrap(), RoaringBitmap::new()); // truncated
        assert_eq!(idx.get_colors(2).unwrap(), make_bitmap(&[3]));

        // 3 unique sets: empty (truncated), {0,1}, {3}
        assert_eq!(stats.num_sets, 3);
    }
}
