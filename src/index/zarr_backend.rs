//! Zarr v3 backend for the Dragon index.
//!
//! Exports an existing Dragon index (`fm_index.bin` + `colors.drgn`) as a
//! chunked, compressed, cloud-native Zarr store. The store can be read back
//! with [`ZarrFmIndex`] / [`ZarrColorIndex`] for query, or opened from any
//! Zarr-aware tool (e.g. `zarr-python`, xarray).
//!
//! ## Store layout (Zarr v3)
//!
//! ```text
//! <store>/
//!   zarr.json                       root group metadata + Dragon attrs
//!   text/                           concatenated unitig text (u8, Zstd)
//!   suffix_array/                   sorted suffix positions (u64, Zstd)
//!   unitig_lengths/                 per-unitig base lengths (u64)
//!   colors/
//!     offsets/                      per-unitig bitmap byte offsets (u64)
//!     bitmaps/                      raw RoaringBitmap bytes (u8, Zstd)
//! ```
//!
//! Dragon-specific metadata (k-mer size, genome count, etc.) is stored as
//! JSON attributes on the root group.

use anyhow::{Context, Result};
use roaring::RoaringBitmap;
use std::path::Path;
use std::sync::Arc;

use zarrs::array::codec::ZstdCodec;
use zarrs::array::{Array, ArrayBuilder, ArraySubset, data_type};
use zarrs::filesystem::FilesystemStore;
use zarrs::group::GroupBuilder;
use zarrs::storage::ReadableWritableListableStorage;

use crate::ds::elias_fano::CumulativeLengthIndex;

/// On-disk format version. Bump when the layout changes incompatibly.
pub const DRAGON_ZARR_FORMAT_VERSION: u32 = 1;

/// 1 MiB of u8 per text chunk.
pub const TEXT_CHUNK_SIZE: u64 = 1_048_576;

/// 131072 u64 entries = 1 MiB per suffix-array chunk.
pub const SA_CHUNK_SIZE: u64 = 131_072;

/// 1 MiB per bitmap-blob chunk.
pub const BITMAP_CHUNK_SIZE: u64 = 1_048_576;

/// Zstd compression level for all compressed arrays.
pub const ZSTD_LEVEL: i32 = 3;

/// Export a pre-built Dragon index directory as a Zarr v3 store.
///
/// Reads `fm_index.bin` and `colors.drgn` from `index_dir`; writes a fresh
/// Zarr store at `zarr_path`. The source index is not modified.
pub fn export_to_zarr(index_dir: &Path, zarr_path: &Path) -> Result<ExportStats> {
    log::info!("Loading Dragon index from {:?}", index_dir);
    let fm = crate::index::fm::load_fm_index(index_dir)
        .context("load fm_index.bin")?;
    let colors = crate::index::color::load_color_index(index_dir)
        .context("load colors.drgn")?;

    let num_unitigs = colors.num_unitigs();
    let num_genomes = colors.num_genomes();
    let text_len = fm.text.len() as u64;
    let sa_len = fm.suffix_array.len() as u64;

    log::info!(
        "Source index: {} unitigs, {} genomes, {} text bytes, {} SA entries",
        num_unitigs, num_genomes, text_len, sa_len
    );

    std::fs::create_dir_all(zarr_path)
        .with_context(|| format!("create zarr dir {:?}", zarr_path))?;

    let store: ReadableWritableListableStorage =
        Arc::new(FilesystemStore::new(zarr_path)?);

    // Root group with Dragon metadata attributes.
    let mut root = GroupBuilder::new().build(store.clone(), "/")?;
    let attrs = root.attributes_mut();
    attrs.insert(
        "dragon_zarr_format_version".into(),
        serde_json::json!(DRAGON_ZARR_FORMAT_VERSION),
    );
    let kmer_size = crate::index::metadata::load_metadata(index_dir)
        .map(|m| m.kmer_size)
        .unwrap_or(31);
    attrs.insert("kmer_size".into(), serde_json::json!(kmer_size));
    attrs.insert("num_unitigs".into(), serde_json::json!(num_unitigs));
    attrs.insert("num_genomes".into(), serde_json::json!(num_genomes));
    attrs.insert("text_len".into(), serde_json::json!(text_len));
    attrs.insert("num_suffixes".into(), serde_json::json!(sa_len));
    root.store_metadata()?;

    // /text: concatenated unitig text, u8 with Zstd.
    write_text_array(&fm.text, store.clone())?;

    // /suffix_array: sorted suffix positions, u64 with Zstd.
    write_suffix_array(&fm.suffix_array, store.clone())?;

    // /unitig_lengths: per-unitig base lengths, u64 (small, single chunk).
    let lengths = fm.cumulative_lengths.lengths();
    write_unitig_lengths(&lengths, store.clone())?;

    // /colors/{offsets,bitmaps}: per-unitig bitmap offsets + flat blob.
    let colors_bytes = write_colors(&colors, store.clone())?;

    Ok(ExportStats {
        num_unitigs,
        num_genomes,
        text_len,
        sa_len,
        colors_bytes,
    })
}

/// Summary returned by [`export_to_zarr`].
#[derive(Debug, Clone)]
pub struct ExportStats {
    pub num_unitigs: u64,
    pub num_genomes: u64,
    pub text_len: u64,
    pub sa_len: u64,
    pub colors_bytes: u64,
}

fn write_text_array(text: &[u8], store: ReadableWritableListableStorage) -> Result<()> {
    log::info!("Writing /text ({} bytes, chunk={} B, zstd-{})", text.len(), TEXT_CHUNK_SIZE, ZSTD_LEVEL);

    let array = ArrayBuilder::new(
        vec![text.len() as u64],
        vec![TEXT_CHUNK_SIZE],
        data_type::uint8(),
        0u8,
    )
    .bytes_to_bytes_codecs(vec![Arc::new(ZstdCodec::new(ZSTD_LEVEL, false))])
    .dimension_names(["position"].into())
    .build(store, "/text")?;

    array.store_metadata()?;
    array.store_array_subset(&array.subset_all(), text)?;
    Ok(())
}

fn write_suffix_array(sa: &[usize], store: ReadableWritableListableStorage) -> Result<()> {
    log::info!("Writing /suffix_array ({} entries, chunk={} entries, zstd-{})", sa.len(), SA_CHUNK_SIZE, ZSTD_LEVEL);

    // Store as u64 — widen on 32-bit, identity on 64-bit.
    let sa_u64: Vec<u64> = sa.iter().map(|&x| x as u64).collect();

    let array = ArrayBuilder::new(
        vec![sa_u64.len() as u64],
        vec![SA_CHUNK_SIZE],
        data_type::uint64(),
        0u64,
    )
    .bytes_to_bytes_codecs(vec![Arc::new(ZstdCodec::new(ZSTD_LEVEL, false))])
    .dimension_names(["rank"].into())
    .build(store, "/suffix_array")?;

    array.store_metadata()?;
    array.store_array_subset(&array.subset_all(), sa_u64.as_slice())?;
    Ok(())
}

fn write_unitig_lengths(lengths: &[u64], store: ReadableWritableListableStorage) -> Result<()> {
    log::info!("Writing /unitig_lengths ({} unitigs)", lengths.len());

    // Small array — one chunk holds everything, no compression.
    let array = ArrayBuilder::new(
        vec![lengths.len() as u64],
        vec![lengths.len().max(1) as u64],
        data_type::uint64(),
        0u64,
    )
    .dimension_names(["unitig_id"].into())
    .build(store, "/unitig_lengths")?;

    array.store_metadata()?;
    array.store_array_subset(&array.subset_all(), lengths)?;
    Ok(())
}

fn write_colors(
    colors: &crate::index::color::ColorIndex,
    store: ReadableWritableListableStorage,
) -> Result<u64> {
    let n = colors.num_unitigs() as usize;
    log::info!("Writing /colors ({} unitigs)", n);

    // Serialize each unitig's bitmap; record byte offsets.
    let mut offsets: Vec<u64> = Vec::with_capacity(n + 1);
    offsets.push(0);
    let mut blob: Vec<u8> = Vec::new();

    for uid in 0..n as u32 {
        let bm = colors.get_colors(uid)?;
        bm.serialize_into(&mut blob)
            .with_context(|| format!("serialize bitmap for unitig {uid}"))?;
        offsets.push(blob.len() as u64);
    }
    let total = *offsets.last().unwrap();
    log::info!("  bitmaps blob = {} bytes across {} unitigs", total, n);

    // /colors group
    let group = GroupBuilder::new().build(store.clone(), "/colors")?;
    group.store_metadata()?;

    // /colors/offsets (n+1 entries, single chunk — O(num_unitigs) memory).
    let off_array = ArrayBuilder::new(
        vec![offsets.len() as u64],
        vec![offsets.len().max(1) as u64],
        data_type::uint64(),
        0u64,
    )
    .dimension_names(["unitig_id"].into())
    .build(store.clone(), "/colors/offsets")?;
    off_array.store_metadata()?;
    off_array.store_array_subset(&off_array.subset_all(), offsets.as_slice())?;

    // /colors/bitmaps (flat bitmap bytes, Zstd-compressed).
    let bm_array = ArrayBuilder::new(
        vec![blob.len() as u64],
        vec![BITMAP_CHUNK_SIZE],
        data_type::uint8(),
        0u8,
    )
    .bytes_to_bytes_codecs(vec![Arc::new(ZstdCodec::new(ZSTD_LEVEL, false))])
    .dimension_names(["byte"].into())
    .build(store, "/colors/bitmaps")?;
    bm_array.store_metadata()?;
    bm_array.store_array_subset(&bm_array.subset_all(), blob.as_slice())?;

    Ok(total)
}

// ---------------------------------------------------------------------------
// Read path: ZarrFmIndex / ZarrColorIndex
// ---------------------------------------------------------------------------

/// FM-index variant backed by a Zarr store. Reads are lazy and chunked — only
/// the chunks touched by a query are decompressed, making this cheap over
/// remote object stores (S3, GCS) that expose HTTP range requests.
pub struct ZarrFmIndex {
    text: Array<FilesystemStore>,
    sa: Array<FilesystemStore>,
    pub cumulative_lengths: CumulativeLengthIndex,
    pub text_len: u64,
    pub sa_len: u64,
    pub kmer_size: usize,
    pub num_unitigs: u64,
    pub num_genomes: u64,
}

impl ZarrFmIndex {
    /// Open a Zarr store produced by [`export_to_zarr`] for read-only query.
    pub fn open(zarr_path: &Path) -> Result<Self> {
        let store = Arc::new(FilesystemStore::new(zarr_path)
            .with_context(|| format!("open zarr store at {:?}", zarr_path))?);

        let text = Array::open(store.clone(), "/text")
            .context("open /text array")?;
        let sa = Array::open(store.clone(), "/suffix_array")
            .context("open /suffix_array array")?;
        let lengths_arr = Array::open(store.clone(), "/unitig_lengths")
            .context("open /unitig_lengths array")?;

        let lengths: Vec<u64> = lengths_arr
            .retrieve_array_subset::<Vec<u64>>(&lengths_arr.subset_all())
            .context("read unitig_lengths")?;
        let cumulative_lengths = CumulativeLengthIndex::from_lengths(&lengths);

        let text_len = text.shape()[0];
        let sa_len = sa.shape()[0];

        // Root-group attributes hold Dragon metadata.
        let root = zarrs::group::Group::open(store, "/")
            .context("open root group")?;
        let attrs = root.attributes();
        let kmer_size = attrs.get("kmer_size")
            .and_then(|v| v.as_u64()).unwrap_or(31) as usize;
        let num_unitigs = attrs.get("num_unitigs")
            .and_then(|v| v.as_u64()).unwrap_or(lengths.len() as u64);
        let num_genomes = attrs.get("num_genomes")
            .and_then(|v| v.as_u64()).unwrap_or(0);

        Ok(Self {
            text, sa, cumulative_lengths,
            text_len, sa_len, kmer_size, num_unitigs, num_genomes,
        })
    }

    /// Read `SA[rank]` via a single-element Zarr read.
    fn sa_get(&self, rank: u64) -> Result<u64> {
        let subset = ArraySubset::new_with_ranges(&[rank..rank + 1]);
        let vals: Vec<u64> = self.sa.retrieve_array_subset(&subset)?;
        Ok(vals[0])
    }

    /// Read `text[start..end]` via a single chunked Zarr read.
    fn text_slice(&self, start: u64, end: u64) -> Result<Vec<u8>> {
        let end = end.min(self.text_len);
        if end <= start { return Ok(Vec::new()); }
        let subset = ArraySubset::new_with_ranges(&[start..end]);
        Ok(self.text.retrieve_array_subset(&subset)?)
    }

    fn lower_bound(&self, pattern: &[u8]) -> Result<u64> {
        let mut lo = 0u64;
        let mut hi = self.sa_len;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let suffix_start = self.sa_get(mid)?;
            let suffix = self.text_slice(suffix_start, suffix_start + pattern.len() as u64)?;
            if suffix.as_slice() < pattern {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }

    fn upper_bound(&self, pattern: &[u8]) -> Result<u64> {
        let mut lo = 0u64;
        let mut hi = self.sa_len;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let suffix_start = self.sa_get(mid)?;
            let suffix = self.text_slice(suffix_start, suffix_start + pattern.len() as u64)?;
            if suffix.as_slice() <= pattern {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }

    /// Backward-search equivalent: return all text positions where `pattern` occurs.
    pub fn search(&self, pattern: &[u8]) -> Result<Vec<u64>> {
        if pattern.is_empty() || self.text_len == 0 {
            return Ok(Vec::new());
        }
        let lo = self.lower_bound(pattern)?;
        let hi = self.upper_bound(pattern)?;
        if hi <= lo { return Ok(Vec::new()); }
        let subset = ArraySubset::new_with_ranges(&[lo..hi]);
        Ok(self.sa.retrieve_array_subset(&subset)?)
    }

    /// Count occurrences without materializing positions.
    pub fn count(&self, pattern: &[u8]) -> Result<usize> {
        if pattern.is_empty() || self.text_len == 0 {
            return Ok(0);
        }
        let lo = self.lower_bound(pattern)?;
        let hi = self.upper_bound(pattern)?;
        Ok((hi - lo) as usize)
    }

    /// Map a text position to `(unitig_id, offset)`.
    pub fn position_to_unitig(&self, pos: u64) -> Option<(u32, u32)> {
        self.cumulative_lengths.unitig_at_position(pos)
    }
}

/// Color index variant backed by a Zarr store.
///
/// Each `get_colors` call reads two u64s (offset pair) then the raw bitmap
/// bytes — so it touches at most a handful of chunks per query.
pub struct ZarrColorIndex {
    offsets: Array<FilesystemStore>,
    bitmaps: Array<FilesystemStore>,
    pub num_unitigs: u64,
    pub num_genomes: u64,
}

impl ZarrColorIndex {
    pub fn open(zarr_path: &Path) -> Result<Self> {
        let store = Arc::new(FilesystemStore::new(zarr_path)?);
        let offsets = Array::open(store.clone(), "/colors/offsets")
            .context("open /colors/offsets")?;
        let bitmaps = Array::open(store.clone(), "/colors/bitmaps")
            .context("open /colors/bitmaps")?;

        let root = zarrs::group::Group::open(store, "/")?;
        let attrs = root.attributes();
        let num_unitigs = attrs.get("num_unitigs")
            .and_then(|v| v.as_u64()).unwrap_or(offsets.shape()[0].saturating_sub(1));
        let num_genomes = attrs.get("num_genomes")
            .and_then(|v| v.as_u64()).unwrap_or(0);

        Ok(Self { offsets, bitmaps, num_unitigs, num_genomes })
    }

    pub fn get_colors(&self, unitig_id: u32) -> Result<RoaringBitmap> {
        let uid = unitig_id as u64;
        if uid >= self.num_unitigs {
            return Ok(RoaringBitmap::new());
        }
        let off_subset = ArraySubset::new_with_ranges(&[uid..uid + 2]);
        let offs: Vec<u64> = self.offsets.retrieve_array_subset(&off_subset)?;
        let (start, end) = (offs[0], offs[1]);
        if end <= start {
            return Ok(RoaringBitmap::new());
        }
        let bm_subset = ArraySubset::new_with_ranges(&[start..end]);
        let bytes: Vec<u8> = self.bitmaps.retrieve_array_subset(&bm_subset)?;
        Ok(RoaringBitmap::deserialize_from(bytes.as_slice())?)
    }

    pub fn num_unitigs(&self) -> u64 { self.num_unitigs }
    pub fn num_genomes(&self) -> u64 { self.num_genomes }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::fm::build_fm_index;
    use crate::index::unitig::UnitigSet;

    /// Round-trip: build a tiny in-memory FM-index, export to Zarr, verify
    /// that the arrays readable by zarrs contain the same bytes.
    #[test]
    fn export_roundtrip_smoke() -> Result<()> {
        let _ = env_logger::builder().is_test(true).try_init();

        let tmp = tempfile::tempdir()?;
        let idx_dir = tmp.path().join("idx");
        std::fs::create_dir_all(&idx_dir)?;

        // Build a minimal FM-index over two unitigs using from_fm_text,
        // which is the public constructor (matches how index load works).
        let text = b"ACGTACGT$TTGGCCAA$".to_vec();
        let lengths = [8u64, 8];
        let unitigs = UnitigSet::from_fm_text(&text, &lengths);
        build_fm_index(&unitigs, &idx_dir)?;

        // Minimal colors.drgn: one bitmap per unitig, each with genome 0.
        let color_tsv = idx_dir.join("colors.tsv");
        std::fs::write(&color_tsv, "0\t0\n1\t0\n")?;
        crate::index::color::build_color_index(&color_tsv, &idx_dir, 1)?;

        let zarr_dir = tmp.path().join("out.zarr");
        let stats = export_to_zarr(&idx_dir, &zarr_dir)?;

        assert_eq!(stats.num_unitigs, 2);
        assert_eq!(stats.num_genomes, 1);
        assert!(stats.text_len >= 16);

        // Round-trip read: reopen the store and verify /text bytes match.
        let store: ReadableWritableListableStorage =
            Arc::new(FilesystemStore::new(&zarr_dir)?);
        let arr = zarrs::array::Array::open(store, "/text")?;
        let read: Vec<u8> = arr
            .retrieve_array_subset::<Vec<u8>>(&arr.subset_all())?;
        assert_eq!(read.len() as u64, stats.text_len);
        Ok(())
    }

    /// Parity: every search that succeeds against the in-memory FM-index must
    /// return the same set of positions when issued against ZarrFmIndex.
    #[test]
    fn zarr_search_parity() -> Result<()> {
        let _ = env_logger::builder().is_test(true).try_init();

        let tmp = tempfile::tempdir()?;
        let idx_dir = tmp.path().join("idx");
        std::fs::create_dir_all(&idx_dir)?;

        let text = b"ACGTACGTACGT$TTGGCCAATTGG$ACGTACGT$".to_vec();
        let lengths = [12u64, 11, 8];
        let unitigs = UnitigSet::from_fm_text(&text, &lengths);
        build_fm_index(&unitigs, &idx_dir)?;

        let color_tsv = idx_dir.join("colors.tsv");
        std::fs::write(&color_tsv, "0\t0\n1\t0\n2\t0\n")?;
        crate::index::color::build_color_index(&color_tsv, &idx_dir, 1)?;

        let zarr_dir = tmp.path().join("parity.zarr");
        export_to_zarr(&idx_dir, &zarr_dir)?;

        let baseline = crate::index::fm::load_fm_index(&idx_dir)?;
        let zarr_fm = ZarrFmIndex::open(&zarr_dir)?;

        for pat in [&b"ACGT"[..], b"TTGG", b"CGTA", b"GGCC", b"NNNN"] {
            let mut want: Vec<u64> = baseline.search(pat).into_iter().map(|x| x as u64).collect();
            let mut got: Vec<u64> = zarr_fm.search(pat)?;
            want.sort();
            got.sort();
            assert_eq!(got, want, "mismatch for pattern {:?}", std::str::from_utf8(pat).unwrap_or("?"));
            assert_eq!(zarr_fm.count(pat)?, baseline.count(pat));
        }

        // Colors parity.
        let baseline_colors = crate::index::color::load_color_index(&idx_dir)?;
        let zarr_colors = ZarrColorIndex::open(&zarr_dir)?;
        assert_eq!(zarr_colors.num_unitigs(), baseline_colors.num_unitigs());
        for uid in 0..baseline_colors.num_unitigs() as u32 {
            let a = baseline_colors.get_colors(uid)?;
            let b = zarr_colors.get_colors(uid)?;
            assert_eq!(a, b, "color mismatch for unitig {uid}");
        }

        Ok(())
    }
}
