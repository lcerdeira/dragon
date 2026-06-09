//! HTTP-backed lazy Zarr reader for Dragon.
//!
//! Enables querying a Dragon Zarr index hosted on any HTTP(S) server (S3
//! pre-signed URLs, nginx, Apache, GitHub Releases, etc.) without downloading
//! the full store.  Only the array chunks actually touched by a query are
//! fetched — for a single 500 bp query, this amounts to O(log n) SA chunks
//! plus a handful of colour bitmap chunks: typically < 5 MB from a store
//! that can be 5 TB on disk.
//!
//! # Protocol
//!
//! Zarr v3 arrays store each chunk as a separate file at a predictable path:
//!
//! ```text
//! {store_url}/suffix_array/c/{chunk_idx}
//! {store_url}/text/c/{text_chunk_idx}
//! {store_url}/colors/offsets/c/0   (small, single chunk)
//! {store_url}/colors/bitmaps/c/{bitmap_chunk_idx}
//! ```
//!
//! Each chunk response is a raw zstd-compressed blob of the array elements
//! (little-endian u64 for SA/offsets, u8 for text/bitmaps), produced by
//! Dragon's `export-zarr` command.
//!
//! # Usage
//!
//! ```no_run
//! use dragon::index::zarr_http::HttpZarrIndex;
//!
//! let idx = HttpZarrIndex::open("https://example.com/saureus.zarr").unwrap();
//! let hits = idx.search(b"ACGTACGTACGT").unwrap();
//! println!("{} hits", hits.len());
//! ```
//!
//! For the `dragon search-zarr` CLI, pass the URL directly:
//! ```text
//! dragon search-zarr --zarr https://example.com/saureus.zarr -q queries.fa
//! ```

use anyhow::{bail, Context, Result};
use roaring::RoaringBitmap;
use serde::Deserialize;

use crate::ds::elias_fano::CumulativeLengthIndex;

/// Dragon Zarr index accessible via HTTP(S).
///
/// Metadata (unitig lengths, genome count, k-mer size) is loaded at
/// construction; array chunks are fetched lazily per query.
pub struct HttpZarrIndex {
    /// Base URL of the Zarr store (no trailing slash).
    base_url: String,
    /// Blocking HTTP client, shared across chunk fetches.
    client: reqwest::blocking::Client,
    /// Metadata from the root group.
    pub num_unitigs: u64,
    pub num_genomes: u64,
    pub kmer_size: usize,
    pub text_len: u64,
    pub sa_len: u64,
    /// Optional genome names (root attribute `genome_names`), indexed by genome id.
    /// Empty if the store predates name export.
    pub genome_names: Vec<String>,
    /// Cumulative unitig lengths — maps text positions to (unitig_id, offset).
    pub cumulative_lengths: CumulativeLengthIndex,
    /// Total entries in /colors/offsets (= num_unitigs + 1).
    color_offsets_len: u64,
    /// SA chunk size (entries per chunk).
    sa_chunk_size: u64,
    /// Text chunk size (bytes per chunk).
    text_chunk_size: u64,
    /// Bitmap chunk size (bytes per chunk).
    bitmap_chunk_size: u64,
    /// Decompressed-chunk cache, keyed by "{array}/{chunk_idx}".
    ///
    /// The binary search on the suffix array re-visits the same top-of-tree
    /// chunks for every k-mer (the root region of the SA is touched by every
    /// search).  Caching decompressed chunks turns the O(queries × log n)
    /// HTTP fetches into O(distinct chunks) — the first query warms the cache,
    /// subsequent queries hit it.  RefCell because search() takes &self.
    chunk_cache: std::cell::RefCell<std::collections::HashMap<String, std::rc::Rc<Vec<u8>>>>,
    /// Running counter of actual HTTP GETs (for diagnostics / cache-hit rate).
    http_gets: std::cell::Cell<u64>,
}

// ─── Zarr metadata deserialization ───────────────────────────────────────────

/// Root group metadata (`zarr.json`).
#[derive(Deserialize, Debug)]
struct ZarrRootMeta {
    attributes: serde_json::Value,
}

/// Zarr v3 array metadata (`{array}/zarr.json`).
#[derive(Deserialize, Debug)]
struct ZarrArrayMeta {
    shape: Vec<u64>,
    chunk_grid: serde_json::Value,
}

impl HttpZarrIndex {
    /// Open a Dragon Zarr store at `url`.
    ///
    /// Fetches the root group metadata and the unitig-lengths array to build
    /// the `CumulativeLengthIndex`.  No large data is downloaded.
    pub fn open(url: &str) -> Result<Self> {
        let base_url = url.trim_end_matches('/').to_owned();
        let client = reqwest::blocking::Client::builder()
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .context("build HTTP client")?;

        // ── Root metadata ─────────────────────────────────────────────────
        let root_meta: ZarrRootMeta = client
            .get(format!("{}/zarr.json", base_url))
            .send()
            .context("GET zarr.json")?
            .json()
            .context("parse zarr.json")?;

        let attrs = &root_meta.attributes;
        let kmer_size = attrs["kmer_size"].as_u64().unwrap_or(31) as usize;
        let num_unitigs = attrs["num_unitigs"].as_u64().unwrap_or(0);
        let num_genomes = attrs["num_genomes"].as_u64().unwrap_or(0);
        let text_len = attrs["text_len"].as_u64().unwrap_or(0);
        let sa_len = attrs["num_suffixes"].as_u64().unwrap_or(0);
        let genome_names: Vec<String> = attrs["genome_names"]
            .as_array()
            .map(|a| a.iter().filter_map(|v| v.as_str().map(str::to_owned)).collect())
            .unwrap_or_default();

        // ── Array metadata (chunk sizes) ──────────────────────────────────
        let sa_meta: ZarrArrayMeta = client
            .get(format!("{}/suffix_array/zarr.json", base_url))
            .send()
            .context("GET suffix_array/zarr.json")?
            .json()
            .context("parse suffix_array zarr.json")?;
        let sa_chunk_size = chunk_size_from_meta(&sa_meta)?;

        let text_meta: ZarrArrayMeta = client
            .get(format!("{}/text/zarr.json", base_url))
            .send()
            .context("GET text/zarr.json")?
            .json()
            .context("parse text zarr.json")?;
        let text_chunk_size = chunk_size_from_meta(&text_meta)?;

        let bm_meta: ZarrArrayMeta = client
            .get(format!("{}/colors/bitmaps/zarr.json", base_url))
            .send()
            .context("GET colors/bitmaps/zarr.json")?
            .json()
            .context("parse colors/bitmaps zarr.json")?;
        let bitmap_chunk_size = chunk_size_from_meta(&bm_meta)?;

        // ── Unitig lengths (small, single chunk) ─────────────────────────
        let lengths_bytes = fetch_and_decompress(
            &client,
            &format!("{}/unitig_lengths/c/0", base_url),
        )
        .context("fetch unitig_lengths chunk")?;
        let lengths: Vec<u64> = bytes_to_u64_le(&lengths_bytes);
        let cumulative_lengths = CumulativeLengthIndex::from_lengths(&lengths);

        log::info!(
            "HttpZarrIndex: {} unitigs, {} genomes, k={}, text={}, SA={} \
             [sa_chunk={}, text_chunk={}, bm_chunk={}]",
            num_unitigs, num_genomes, kmer_size, text_len, sa_len,
            sa_chunk_size, text_chunk_size, bitmap_chunk_size,
        );

        Ok(Self {
            base_url,
            client,
            num_unitigs,
            num_genomes,
            kmer_size,
            text_len,
            sa_len,
            genome_names,
            cumulative_lengths,
            color_offsets_len: num_unitigs + 1,
            sa_chunk_size,
            text_chunk_size,
            bitmap_chunk_size,
            chunk_cache: std::cell::RefCell::new(std::collections::HashMap::new()),
            http_gets: std::cell::Cell::new(0),
        })
    }

    /// Number of actual HTTP GETs performed (cache misses).
    pub fn http_get_count(&self) -> u64 {
        self.http_gets.get()
    }

    // ── SA access ────────────────────────────────────────────────────────────

    fn sa_get(&self, rank: u64) -> Result<u64> {
        if rank >= self.sa_len {
            bail!("SA rank {} out of bounds (len={})", rank, self.sa_len);
        }
        let chunk_idx = rank / self.sa_chunk_size;
        let offset = (rank % self.sa_chunk_size) as usize;
        let bytes = self.fetch_chunk("suffix_array", chunk_idx)?;
        let entries = bytes_to_u64_le(&bytes);
        entries
            .get(offset)
            .copied()
            .with_context(|| format!("SA offset {} out of chunk ({} entries)", offset, entries.len()))
    }

    fn sa_range(&self, lo: u64, hi: u64) -> Result<Vec<u64>> {
        if hi <= lo { return Ok(Vec::new()); }
        let mut result = Vec::with_capacity((hi - lo) as usize);
        // Group by chunk to minimise HTTP requests.
        let first_chunk = lo / self.sa_chunk_size;
        let last_chunk = (hi - 1) / self.sa_chunk_size;
        for chunk_idx in first_chunk..=last_chunk {
            let bytes = self.fetch_chunk("suffix_array", chunk_idx)?;
            let entries = bytes_to_u64_le(&bytes);
            let chunk_start = chunk_idx * self.sa_chunk_size;
            let from = (lo.max(chunk_start) - chunk_start) as usize;
            let to = ((hi.min(chunk_start + self.sa_chunk_size)) - chunk_start) as usize;
            result.extend_from_slice(&entries[from..to.min(entries.len())]);
        }
        Ok(result)
    }

    // ── Text access ──────────────────────────────────────────────────────────

    fn text_slice(&self, start: u64, end: u64) -> Result<Vec<u8>> {
        let end = end.min(self.text_len);
        if end <= start { return Ok(Vec::new()); }
        let first_chunk = start / self.text_chunk_size;
        let last_chunk = (end - 1) / self.text_chunk_size;
        let mut result = Vec::with_capacity((end - start) as usize);
        for chunk_idx in first_chunk..=last_chunk {
            let bytes = self.fetch_chunk("text", chunk_idx)?;
            let chunk_start = chunk_idx * self.text_chunk_size;
            let from = (start.max(chunk_start) - chunk_start) as usize;
            let to = ((end.min(chunk_start + self.text_chunk_size)) - chunk_start) as usize;
            result.extend_from_slice(&bytes[from..to.min(bytes.len())]);
        }
        Ok(result)
    }

    // ── Binary search ────────────────────────────────────────────────────────

    fn lower_bound(&self, pattern: &[u8]) -> Result<u64> {
        let mut lo = 0u64;
        let mut hi = self.sa_len;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let suffix_start = self.sa_get(mid)?;
            let suffix = self.text_slice(suffix_start, suffix_start + pattern.len() as u64)?;
            if suffix.as_slice() < pattern { lo = mid + 1; } else { hi = mid; }
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
            if suffix.as_slice() <= pattern { lo = mid + 1; } else { hi = mid; }
        }
        Ok(lo)
    }

    // ── Public search interface ───────────────────────────────────────────────

    /// Find all text positions where `pattern` occurs.
    ///
    /// Issues O(log n) HTTP requests for binary search bounds, then one
    /// request for the SA range.  Total bytes fetched ≈ chunks × chunk_size.
    pub fn search(&self, pattern: &[u8]) -> Result<Vec<u64>> {
        if pattern.is_empty() || self.text_len == 0 {
            return Ok(Vec::new());
        }
        let lo = self.lower_bound(pattern)?;
        let hi = self.upper_bound(pattern)?;
        self.sa_range(lo, hi)
    }

    /// Map a text position to `(unitig_id, offset_within_unitig)`.
    pub fn position_to_unitig(&self, pos: u64) -> Option<(u32, u32)> {
        self.cumulative_lengths.unitig_at_position(pos)
    }

    // ── Colour index ─────────────────────────────────────────────────────────

    /// Fetch the genome membership bitmap for `unitig_id`.
    ///
    /// Issues exactly 2 HTTP requests: one for the offset pair, one for the
    /// bitmap bytes.
    pub fn get_colors(&self, unitig_id: u32) -> Result<RoaringBitmap> {
        let uid = unitig_id as u64;
        if uid >= self.num_unitigs {
            return Ok(RoaringBitmap::new());
        }

        // /colors/offsets is stored as a single chunk (small).
        let off_bytes = self.fetch_chunk("colors/offsets", 0)?;
        let offsets = bytes_to_u64_le(&off_bytes);
        let (start, end) = (
            offsets.get(uid as usize).copied().unwrap_or(0),
            offsets.get(uid as usize + 1).copied().unwrap_or(0),
        );
        if end <= start {
            return Ok(RoaringBitmap::new());
        }

        // /colors/bitmaps chunked.
        let first_chunk = start / self.bitmap_chunk_size;
        let last_chunk = (end - 1) / self.bitmap_chunk_size;
        let mut blob = Vec::new();
        for chunk_idx in first_chunk..=last_chunk {
            let bytes = self.fetch_chunk("colors/bitmaps", chunk_idx)?;
            let chunk_start = chunk_idx * self.bitmap_chunk_size;
            let from = (start.max(chunk_start) - chunk_start) as usize;
            let to = ((end.min(chunk_start + self.bitmap_chunk_size)) - chunk_start) as usize;
            blob.extend_from_slice(&bytes[from..to.min(bytes.len())]);
        }

        Ok(RoaringBitmap::deserialize_from(blob.as_slice())
            .context("deserialize RoaringBitmap from HTTP store")?)
    }

    // ── Internal HTTP chunk fetcher ───────────────────────────────────────────

    fn fetch_chunk(&self, array_path: &str, chunk_idx: u64) -> Result<std::rc::Rc<Vec<u8>>> {
        let key = format!("{}/{}", array_path, chunk_idx);

        // Cache hit: return the shared Rc without any HTTP traffic.
        if let Some(cached) = self.chunk_cache.borrow().get(&key) {
            return Ok(cached.clone());
        }

        // Cache miss: fetch + decompress, then store.
        let url = format!("{}/{}/c/{}", self.base_url, array_path, chunk_idx);
        let resp = self.client
            .get(&url)
            .send()
            .with_context(|| format!("GET {}", url))?;
        if !resp.status().is_success() {
            bail!("HTTP {} for {}", resp.status(), url);
        }
        self.http_gets.set(self.http_gets.get() + 1);

        let raw = resp.bytes().context("read HTTP response body")?.to_vec();
        let decoded = std::rc::Rc::new(
            maybe_zstd_decode(raw).with_context(|| format!("decode chunk {}", url))?,
        );
        self.chunk_cache.borrow_mut().insert(key, decoded.clone());
        Ok(decoded)
    }
}

/// Decompress a chunk only if it is a zstd frame.
///
/// Dragon's Zarr export compresses `/text`, `/suffix_array`, `/colors/bitmaps`
/// and `/paths/blobs` with zstd, but stores `/unitig_lengths`, `/colors/offsets`
/// and `/paths/offsets` UNCOMPRESSED.  Rather than track which array is which,
/// we detect the zstd magic number (0x28 0xB5 0x2F 0xFD) at the start of the
/// chunk and only decode when present; otherwise the bytes are returned as-is.
fn maybe_zstd_decode(bytes: Vec<u8>) -> Result<Vec<u8>> {
    const ZSTD_MAGIC: [u8; 4] = [0x28, 0xB5, 0x2F, 0xFD];
    if bytes.len() >= 4 && bytes[0..4] == ZSTD_MAGIC {
        Ok(zstd::decode_all(bytes.as_slice())?)
    } else {
        Ok(bytes)
    }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn bytes_to_u64_le(bytes: &[u8]) -> Vec<u64> {
    bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn chunk_size_from_meta(meta: &ZarrArrayMeta) -> Result<u64> {
    // Zarr v3: chunk_grid = { "name": "regular", "configuration": { "chunk_shape": [N] } }
    let size = meta.chunk_grid["configuration"]["chunk_shape"][0]
        .as_u64()
        .with_context(|| format!("parse chunk_shape from {:?}", meta.chunk_grid))?;
    Ok(size)
}

fn fetch_and_decompress(
    client: &reqwest::blocking::Client,
    url: &str,
) -> Result<Vec<u8>> {
    let resp = client.get(url).send().with_context(|| format!("GET {}", url))?;
    if !resp.status().is_success() {
        bail!("HTTP {} for {}", resp.status(), url);
    }
    let raw = resp.bytes().context("read body")?.to_vec();
    // /unitig_lengths is stored uncompressed; only decode if it's a zstd frame.
    maybe_zstd_decode(raw)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maybe_zstd_decode_passes_through_raw() {
        // Raw (uncompressed) bytes — no zstd magic — returned as-is
        let raw = vec![1u8, 0, 0, 0, 0, 0, 0, 0]; // a u64=1 LE, uncompressed
        let out = maybe_zstd_decode(raw.clone()).unwrap();
        assert_eq!(out, raw, "non-zstd bytes must pass through unchanged");
    }

    #[test]
    fn maybe_zstd_decode_handles_zstd() {
        // A real zstd frame round-trips
        let original = b"ACGTACGTACGT".to_vec();
        let compressed = zstd::encode_all(original.as_slice(), 3).unwrap();
        assert_eq!(&compressed[0..4], &[0x28, 0xB5, 0x2F, 0xFD], "zstd magic present");
        let out = maybe_zstd_decode(compressed).unwrap();
        assert_eq!(out, original, "zstd frame must decompress");
    }

    #[test]
    fn bytes_to_u64_le_basic() {
        let bytes: Vec<u8> = (0u64..4).flat_map(|i| i.to_le_bytes()).collect();
        let vals = bytes_to_u64_le(&bytes);
        assert_eq!(vals, vec![0, 1, 2, 3]);
    }

    #[test]
    fn bytes_to_u64_le_partial_ignored() {
        // 11 bytes → 1 complete u64 + 3 ignored
        let bytes: Vec<u8> = vec![1, 0, 0, 0, 0, 0, 0, 0, 0xFF, 0xFF, 0xFF];
        let vals = bytes_to_u64_le(&bytes);
        assert_eq!(vals, vec![1u64]);
    }
}
