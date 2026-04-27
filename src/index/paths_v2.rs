//! `paths.bin` v2 — fixed-layout, mmap-friendly format for the genome path index.
//!
//! ## Why
//!
//! The legacy `bincode::serialize(&PathIndex)` format slurps the entire
//! `Vec<GenomePath>` into RAM at load time. For 32K-genome shards that's a
//! 480 GB read served from NFS — currently ~13 minutes per shard during
//! multi-shard search.
//!
//! v2 stores a fixed-position offset table after the header so a reader can
//! `mmap()` the file and jump directly to genome `i`'s blob in O(1) without
//! reading any of the surrounding data. With OS page cache + lazy fault-in,
//! cold-start search latency drops from minutes to milliseconds.
//!
//! ## Wire format
//!
//! ```text
//! +0   magic[8]            = b"DRGNPTH2"
//! +8   version: u32 LE     = 1
//! +12  num_genomes: u64 LE
//! +20  reserved: u64       = 0   (must be 0 in v1)
//! +28  offsets[num_genomes + 1]: u64 LE
//!         offsets[i]   = absolute byte offset of genome i's blob
//!         offsets[N]   = end-of-file (one past last byte)
//! +28+8(N+1)
//!      blob_0, blob_1, ...
//! ```
//!
//! Per-genome blob (LEB128 / varint encoded throughout):
//!
//! ```text
//!   varint name_length         (bytes)
//!   bytes  name                (UTF-8, name_length bytes)
//!   varint genome_length       (bp)
//!   varint num_steps
//!   repeated num_steps times:
//!     varint unitig_with_rev   (low bit = is_reverse, rest = unitig_id)
//!     varint delta_offset      (genome_offset - previous step's genome_offset;
//!                               first step's delta is relative to 0)
//! ```
//!
//! Delta encoding shrinks per-step storage to typically ~3-5 bytes
//! (1 byte for the ~k bp stride, 2-4 bytes for the unitig id) versus
//! 8 bytes/step for the bincode varint encoding of `Vec<PathStep>`.

use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use crate::index::paths::{GenomePath, PathStep};

pub const MAGIC: &[u8; 8] = b"DRGNPTH2";
pub const VERSION: u32 = 1;

/// Header size: magic(8) + version(4) + num_genomes(8) + reserved(8) = 28
pub const HEADER_SIZE: u64 = 28;

// ---------------------------------------------------------------------------
// LEB128 varint
// ---------------------------------------------------------------------------

/// Append a LEB128-encoded `value` to `out`. 1–10 bytes for u64.
fn write_varint(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push((value as u8) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

/// Decode a LEB128 varint at `bytes[pos..]`. Returns (value, new_pos).
/// Errors if the varint runs off the end of the buffer.
fn read_varint(bytes: &[u8], pos: usize) -> Result<(u64, usize)> {
    let mut value: u64 = 0;
    let mut shift: u32 = 0;
    let mut p = pos;
    loop {
        if p >= bytes.len() {
            bail!("varint truncated at offset {p}");
        }
        let byte = bytes[p];
        p += 1;
        value |= u64::from(byte & 0x7F) << shift;
        if byte & 0x80 == 0 {
            return Ok((value, p));
        }
        shift += 7;
        if shift >= 64 {
            bail!("varint exceeds u64 range at offset {pos}");
        }
    }
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Streaming writer. Per-genome RAM stays bounded (one blob in flight at a
/// time); the only growth is the offset table itself (8 × num_genomes bytes).
pub struct PathV2Writer {
    file: BufWriter<File>,
    num_genomes: u64,
    offsets: Vec<u64>,
    written: u64,
}

impl PathV2Writer {
    pub fn create(path: &Path, num_genomes: u64) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("create {path:?}"))?;
        let mut w = BufWriter::with_capacity(16 * 1024 * 1024, file);

        // Header
        w.write_all(MAGIC)?;
        w.write_all(&VERSION.to_le_bytes())?;
        w.write_all(&num_genomes.to_le_bytes())?;
        w.write_all(&0u64.to_le_bytes())?; // reserved

        // Reserve offset table — we'll seek back to fill it in `finish`.
        let table_bytes = (num_genomes as usize + 1) * 8;
        let zeros = vec![0u8; table_bytes];
        w.write_all(&zeros)?;

        let body_start = HEADER_SIZE + table_bytes as u64;
        Ok(Self {
            file: w,
            num_genomes,
            offsets: {
                let mut v = Vec::with_capacity(num_genomes as usize + 1);
                v.push(body_start);
                v
            },
            written: body_start,
        })
    }

    pub fn write_genome(&mut self, g: &GenomePath) -> Result<()> {
        if self.offsets.len() as u64 > self.num_genomes {
            bail!(
                "write_genome called {} times; declared {}",
                self.offsets.len(),
                self.num_genomes
            );
        }

        let mut blob: Vec<u8> = Vec::with_capacity(64 + g.steps.len() * 4);
        write_varint(&mut blob, g.genome_name.len() as u64);
        blob.extend_from_slice(g.genome_name.as_bytes());
        write_varint(&mut blob, g.genome_length);
        write_varint(&mut blob, g.steps.len() as u64);

        let mut prev: u64 = 0;
        for s in &g.steps {
            let mixed = (u64::from(s.unitig_id) << 1) | u64::from(s.is_reverse);
            write_varint(&mut blob, mixed);
            // genome_offset is monotonic non-decreasing within a genome
            // (build_path_index emits steps in increasing offset order).
            let delta = s.genome_offset.checked_sub(prev).ok_or_else(|| {
                anyhow::anyhow!(
                    "non-monotonic genome_offset in genome {}: {} < prev {}",
                    g.genome_id,
                    s.genome_offset,
                    prev
                )
            })?;
            write_varint(&mut blob, delta);
            prev = s.genome_offset;
        }

        self.file.write_all(&blob)?;
        self.written += blob.len() as u64;
        self.offsets.push(self.written);
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        if (self.offsets.len() as u64) != self.num_genomes + 1 {
            bail!(
                "expected {} genomes, wrote {}",
                self.num_genomes,
                self.offsets.len() - 1
            );
        }
        self.file.flush()?;

        // Seek back and patch the offset table.
        let mut file = self.file.into_inner()?;
        file.seek(SeekFrom::Start(HEADER_SIZE))?;
        let mut bytes: Vec<u8> = Vec::with_capacity(self.offsets.len() * 8);
        for &o in &self.offsets {
            bytes.extend_from_slice(&o.to_le_bytes());
        }
        file.write_all(&bytes)?;
        file.flush()?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Mmap reader
// ---------------------------------------------------------------------------

/// Lazy, mmap-backed view over a v2 paths file. Open is O(1); per-genome
/// access decodes that genome's blob on demand.
pub struct MmapPathIndex {
    mmap: Mmap,
    num_genomes: u64,
    offsets_start: usize,
}

impl MmapPathIndex {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {path:?}"))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("mmap {path:?}"))?;

        if mmap.len() < HEADER_SIZE as usize {
            bail!("file too small for header: {} bytes", mmap.len());
        }
        if &mmap[..8] != MAGIC {
            bail!("not a v2 paths file (magic mismatch)");
        }
        let version = u32::from_le_bytes(mmap[8..12].try_into().unwrap());
        if version != VERSION {
            bail!("unsupported v2 version: {version}");
        }
        let num_genomes = u64::from_le_bytes(mmap[12..20].try_into().unwrap());

        let offsets_start = HEADER_SIZE as usize;
        let offsets_end = offsets_start + (num_genomes as usize + 1) * 8;
        if mmap.len() < offsets_end {
            bail!("file truncated within offset table");
        }

        Ok(Self {
            mmap,
            num_genomes,
            offsets_start,
        })
    }

    pub fn num_genomes(&self) -> u64 {
        self.num_genomes
    }

    /// Read the i-th offset table entry.
    fn offset_at(&self, i: u64) -> u64 {
        let pos = self.offsets_start + (i as usize) * 8;
        u64::from_le_bytes(self.mmap[pos..pos + 8].try_into().unwrap())
    }

    /// Decode genome `genome_id`'s full path. Returns `None` if out of range.
    pub fn get_path(&self, genome_id: u32) -> Result<Option<GenomePath>> {
        let gid = u64::from(genome_id);
        if gid >= self.num_genomes {
            return Ok(None);
        }
        let start = self.offset_at(gid) as usize;
        let end = self.offset_at(gid + 1) as usize;
        if end > self.mmap.len() || start > end {
            bail!("corrupt offset table for genome {gid}: [{start}, {end})");
        }
        let blob = &self.mmap[start..end];

        let mut p = 0usize;
        let (name_len, np) = read_varint(blob, p)?;
        p = np;
        if p + name_len as usize > blob.len() {
            bail!("genome {gid}: name overruns blob");
        }
        let name = std::str::from_utf8(&blob[p..p + name_len as usize])
            .with_context(|| format!("genome {gid}: name not UTF-8"))?
            .to_owned();
        p += name_len as usize;

        let (genome_length, np) = read_varint(blob, p)?;
        p = np;
        let (num_steps, np) = read_varint(blob, p)?;
        p = np;

        let mut steps: Vec<PathStep> = Vec::with_capacity(num_steps as usize);
        let mut prev: u64 = 0;
        for _ in 0..num_steps {
            let (mixed, np) = read_varint(blob, p)?;
            p = np;
            let (delta, np) = read_varint(blob, p)?;
            p = np;
            prev = prev
                .checked_add(delta)
                .ok_or_else(|| anyhow::anyhow!("genome {gid}: offset overflow"))?;
            steps.push(PathStep {
                unitig_id: (mixed >> 1) as u32,
                is_reverse: (mixed & 1) == 1,
                genome_offset: prev,
            });
        }

        Ok(Some(GenomePath {
            genome_id,
            genome_name: name,
            genome_length,
            steps,
        }))
    }
}

/// Magic-peeking helper: returns true if `path` starts with the v2 magic.
pub fn is_v2(path: &Path) -> Result<bool> {
    use std::io::Read;
    let mut f = File::open(path)?;
    let mut head = [0u8; 8];
    if f.read(&mut head)? < 8 {
        return Ok(false);
    }
    Ok(&head == MAGIC)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_genome(id: u32, name: &str, n_steps: u32) -> GenomePath {
        let mut steps = Vec::new();
        let mut off: u64 = 0;
        for i in 0..n_steps {
            steps.push(PathStep {
                unitig_id: i % 7,
                is_reverse: i % 2 == 0,
                genome_offset: off,
            });
            off += (i as u64 % 31) + 1;
        }
        GenomePath {
            genome_id: id,
            genome_name: name.to_owned(),
            genome_length: off,
            steps,
        }
    }

    #[test]
    fn roundtrip_three_genomes() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let path = tmp.path().join("paths_v2.bin");

        let genomes = vec![
            make_genome(0, "g_alpha", 10),
            make_genome(1, "g_bravo", 200),
            make_genome(2, "g_charlie_with_a_long_name_to_exercise_varint", 1234),
        ];

        let mut w = PathV2Writer::create(&path, genomes.len() as u64)?;
        for g in &genomes {
            w.write_genome(g)?;
        }
        w.finish()?;

        assert!(is_v2(&path)?);

        let r = MmapPathIndex::open(&path)?;
        assert_eq!(r.num_genomes(), 3);

        for orig in &genomes {
            let got = r.get_path(orig.genome_id)?.unwrap();
            assert_eq!(got.genome_id, orig.genome_id);
            assert_eq!(got.genome_name, orig.genome_name);
            assert_eq!(got.genome_length, orig.genome_length);
            assert_eq!(got.steps.len(), orig.steps.len());
            for (a, b) in got.steps.iter().zip(orig.steps.iter()) {
                assert_eq!(a.unitig_id, b.unitig_id);
                assert_eq!(a.is_reverse, b.is_reverse);
                assert_eq!(a.genome_offset, b.genome_offset);
            }
        }

        // Out-of-range
        assert!(r.get_path(99)?.is_none());
        Ok(())
    }

    #[test]
    fn empty_path_roundtrips() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let path = tmp.path().join("empty.bin");

        let g = GenomePath {
            genome_id: 0,
            genome_name: "lonely".to_owned(),
            genome_length: 0,
            steps: vec![],
        };
        let mut w = PathV2Writer::create(&path, 1)?;
        w.write_genome(&g)?;
        w.finish()?;

        let r = MmapPathIndex::open(&path)?;
        let got = r.get_path(0)?.unwrap();
        assert_eq!(got.steps.len(), 0);
        assert_eq!(got.genome_name, "lonely");
        Ok(())
    }

    #[test]
    fn varint_roundtrip() {
        for &v in &[0u64, 1, 127, 128, 16384, u32::MAX as u64, u64::MAX] {
            let mut buf = Vec::new();
            write_varint(&mut buf, v);
            let (out, n) = read_varint(&buf, 0).unwrap();
            assert_eq!(out, v);
            assert_eq!(n, buf.len());
        }
    }

    #[test]
    fn rejects_wrong_magic() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let path = tmp.path().join("bad.bin");
        std::fs::write(&path, b"NOTDRGN!\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00")?;
        assert!(MmapPathIndex::open(&path).is_err());
        assert!(!is_v2(&path)?);
        Ok(())
    }
}
