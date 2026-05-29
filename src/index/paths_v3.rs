//! `paths.bin` v3 — graph-edge encoding for the genome path index.
//!
//! See `docs/design/paths-bin-v3.md` for the full rationale. In short: v2
//! stores each genome's de Bruijn walk as `num_steps` full unitig IDs
//! (~3.2 bytes each — 76% of a 41 GB shard). v3 stores the unitig
//! *adjacency* once in a shared CSR table and encodes each genome path as
//! a start node plus per-step **edge indices** into that table.
//!
//! `get_path` reconstructs the identical [`GenomePath`] as v2, so every
//! downstream consumer is unchanged.
//!
//! ## Wire format (`DRGNPTH3`, all integers little-endian)
//!
//! ```text
//! Header (56 bytes):
//!   +0   magic[8]            = b"DRGNPTH3"
//!   +8   version: u32        = 1
//!   +12  _pad: u32           = 0
//!   +16  num_genomes: u64
//!   +24  num_unitigs: u64
//!   +32  succ_offsets_pos: u64
//!   +40  succ_values_pos: u64
//!   +48  genome_table_pos: u64
//!
//! Successor table (CSR) — node = (unitig_id << 1) | is_reverse:
//!   succ_offsets: [u64; 2*num_unitigs + 1]   at succ_offsets_pos
//!   succ_values:  [u32; total_edges]         at succ_values_pos
//!     successors of `node` = succ_values[succ_offsets[node] .. succ_offsets[node+1]]
//!     each value = (next_unitig_id << 1) | next_is_reverse, ascending
//!
//! Genome offset table:
//!   genome_offsets: [u64; num_genomes + 1]   at genome_table_pos
//!
//! Per-genome blob (LEB128 varints):
//!   varint name_length, bytes name
//!   varint genome_length
//!   varint num_steps
//!   if num_steps >= 1:
//!     varint first_node           ((unitig_id<<1)|is_reverse of step 0)
//!     varint first_genome_offset
//!   repeat (num_steps - 1):
//!     varint edge_index           (index into current node's successor slice)
//!     varint delta_offset
//! ```

use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use crate::index::paths::{GenomePath, PathStep};

pub const MAGIC: &[u8; 8] = b"DRGNPTH3";
pub const VERSION: u32 = 1;
pub const HEADER_SIZE: u64 = 56;

#[inline]
fn node_of(unitig_id: u32, is_reverse: bool) -> u32 {
    (unitig_id << 1) | (is_reverse as u32)
}

// ---------------------------------------------------------------------------
// LEB128 varint
// ---------------------------------------------------------------------------

fn write_varint(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push((value as u8) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

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
// CSR successor table
// ---------------------------------------------------------------------------

/// Compressed-sparse-row adjacency over de Bruijn nodes.
struct SuccTable {
    /// `offsets[node] .. offsets[node+1]` indexes into `values`.
    offsets: Vec<u64>,
    /// successor node IDs, each slice ascending.
    values: Vec<u32>,
}

impl SuccTable {
    /// Build from observed transitions across a set of genome paths.
    fn build<'a>(paths: impl Iterator<Item = &'a GenomePath>, num_nodes: u64) -> Self {
        let mut succ: HashMap<u32, BTreeSet<u32>> = HashMap::new();
        for path in paths {
            for w in path.steps.windows(2) {
                let a = node_of(w[0].unitig_id, w[0].is_reverse);
                let b = node_of(w[1].unitig_id, w[1].is_reverse);
                succ.entry(a).or_default().insert(b);
            }
        }
        let mut offsets: Vec<u64> = Vec::with_capacity(num_nodes as usize + 1);
        let mut values: Vec<u32> = Vec::new();
        offsets.push(0);
        for node in 0..num_nodes as u32 {
            if let Some(set) = succ.get(&node) {
                values.extend(set.iter().copied()); // BTreeSet iterates ascending
            }
            offsets.push(values.len() as u64);
        }
        SuccTable { offsets, values }
    }

    /// Successor slice for `node`.
    fn successors(&self, node: u32) -> &[u32] {
        let lo = self.offsets[node as usize] as usize;
        let hi = self.offsets[node as usize + 1] as usize;
        &self.values[lo..hi]
    }
}

// ---------------------------------------------------------------------------
// Encoding
// ---------------------------------------------------------------------------

/// Encode one genome path into a v3 blob using the successor table.
fn encode_blob(path: &GenomePath, succ: &SuccTable) -> Result<Vec<u8>> {
    let mut blob: Vec<u8> = Vec::with_capacity(32 + path.steps.len() * 2);
    write_varint(&mut blob, path.genome_name.len() as u64);
    blob.extend_from_slice(path.genome_name.as_bytes());
    write_varint(&mut blob, path.genome_length);
    write_varint(&mut blob, path.steps.len() as u64);

    if let Some(first) = path.steps.first() {
        let mut prev_node = node_of(first.unitig_id, first.is_reverse);
        let mut prev_off = first.genome_offset;
        write_varint(&mut blob, u64::from(prev_node));
        write_varint(&mut blob, prev_off);

        for s in &path.steps[1..] {
            let target = node_of(s.unitig_id, s.is_reverse);
            let slice = succ.successors(prev_node);
            let edge_index = slice.binary_search(&target).map_err(|_| {
                anyhow::anyhow!(
                    "transition {prev_node}->{target} absent from successor table"
                )
            })?;
            write_varint(&mut blob, edge_index as u64);
            let delta = s.genome_offset.checked_sub(prev_off).ok_or_else(|| {
                anyhow::anyhow!("non-monotonic genome_offset in genome {}", path.genome_id)
            })?;
            write_varint(&mut blob, delta);
            prev_node = target;
            prev_off = s.genome_offset;
        }
    }
    Ok(blob)
}

/// Write a complete v3 file from a slice of genome paths (in-memory; used by
/// tests and small indexes). For large shards use [`migrate_v2_to_v3`].
pub fn write_v3_from_paths(path: &Path, genomes: &[GenomePath]) -> Result<()> {
    let num_genomes = genomes.len() as u64;
    let max_unitig = genomes
        .iter()
        .flat_map(|g| g.steps.iter())
        .map(|s| s.unitig_id)
        .max()
        .unwrap_or(0);
    let num_unitigs = u64::from(max_unitig) + 1;
    let num_nodes = num_unitigs * 2;

    let succ = SuccTable::build(genomes.iter(), num_nodes);
    write_v3_file(
        path,
        num_genomes,
        num_unitigs,
        &succ,
        genomes.iter().cloned().map(Ok),
    )
}

/// Core writer: header + CSR table + offset table + blobs. Consumes genomes
/// lazily (one `GenomePath` resident at a time) so it scales to huge shards.
fn write_v3_file(
    path: &Path,
    num_genomes: u64,
    num_unitigs: u64,
    succ: &SuccTable,
    genomes: impl Iterator<Item = Result<GenomePath>>,
) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {path:?}"))?;
    let mut w = BufWriter::with_capacity(16 * 1024 * 1024, file);

    let succ_offsets_pos = HEADER_SIZE;
    let succ_values_pos = succ_offsets_pos + (succ.offsets.len() as u64) * 8;
    let genome_table_pos = succ_values_pos + (succ.values.len() as u64) * 4;
    let body_start = genome_table_pos + (num_genomes + 1) * 8;

    // Header
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;
    w.write_all(&0u32.to_le_bytes())?; // pad
    w.write_all(&num_genomes.to_le_bytes())?;
    w.write_all(&num_unitigs.to_le_bytes())?;
    w.write_all(&succ_offsets_pos.to_le_bytes())?;
    w.write_all(&succ_values_pos.to_le_bytes())?;
    w.write_all(&genome_table_pos.to_le_bytes())?;

    // Successor table
    for &o in &succ.offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for &v in &succ.values {
        w.write_all(&v.to_le_bytes())?;
    }

    // Reserve the genome offset table; patched in after the blobs.
    w.write_all(&vec![0u8; (num_genomes as usize + 1) * 8])?;

    let mut genome_offsets: Vec<u64> = Vec::with_capacity(num_genomes as usize + 1);
    genome_offsets.push(body_start);
    let mut written = body_start;
    let mut count = 0u64;
    for g in genomes {
        let g = g?;
        let blob = encode_blob(&g, succ)?;
        w.write_all(&blob)?;
        written += blob.len() as u64;
        genome_offsets.push(written);
        count += 1;
    }
    if count != num_genomes {
        bail!("expected {num_genomes} genomes, wrote {count}");
    }

    w.flush()?;
    let mut file = w.into_inner()?;
    file.seek(SeekFrom::Start(genome_table_pos))?;
    let mut bytes = Vec::with_capacity(genome_offsets.len() * 8);
    for &o in &genome_offsets {
        bytes.extend_from_slice(&o.to_le_bytes());
    }
    file.write_all(&bytes)?;
    file.flush()?;
    Ok(())
}

/// Summary returned by [`migrate_v2_to_v3`].
#[derive(Debug, Clone)]
pub struct MigrationV3Stats {
    pub num_genomes: u64,
    pub num_unitigs: u64,
    pub num_edges: u64,
    pub old_size: u64,
    pub new_size: u64,
}

/// Migrate an existing v2 `paths.bin` to v3. Two streaming passes over the
/// v2 file (pass 1 builds the successor table, pass 2 encodes blobs); never
/// holds more than one `GenomePath` in RAM at a time.
pub fn migrate_v2_to_v3(v2_path: &Path, v3_path: &Path) -> Result<MigrationV3Stats> {
    use crate::index::paths_v2::MmapPathIndex;

    let src = MmapPathIndex::open(v2_path)
        .with_context(|| format!("open v2 source {v2_path:?}"))?;
    let num_genomes = src.num_genomes();

    // Pass 1: observed transitions + max unitig id.
    let mut succ_sets: HashMap<u32, BTreeSet<u32>> = HashMap::new();
    let mut max_unitig: u32 = 0;
    for gid in 0..num_genomes {
        let path = src
            .get_path(gid as u32)?
            .ok_or_else(|| anyhow::anyhow!("v2 genome {gid} missing"))?;
        for w in path.steps.windows(2) {
            let a = node_of(w[0].unitig_id, w[0].is_reverse);
            let b = node_of(w[1].unitig_id, w[1].is_reverse);
            succ_sets.entry(a).or_default().insert(b);
        }
        for s in &path.steps {
            max_unitig = max_unitig.max(s.unitig_id);
        }
    }
    let num_unitigs = u64::from(max_unitig) + 1;
    let num_nodes = num_unitigs * 2;

    let mut offsets: Vec<u64> = Vec::with_capacity(num_nodes as usize + 1);
    let mut values: Vec<u32> = Vec::new();
    offsets.push(0);
    for node in 0..num_nodes as u32 {
        if let Some(set) = succ_sets.get(&node) {
            values.extend(set.iter().copied());
        }
        offsets.push(values.len() as u64);
    }
    let num_edges = values.len() as u64;
    let succ = SuccTable { offsets, values };
    drop(succ_sets);

    // Pass 2: re-decode each genome from the v2 mmap on demand and stream
    // it straight into the v3 writer — one GenomePath resident at a time.
    let genome_iter = (0..num_genomes).map(|gid| {
        src.get_path(gid as u32)?
            .ok_or_else(|| anyhow::anyhow!("v2 genome {gid} missing in pass 2"))
    });
    write_v3_file(v3_path, num_genomes as u64, num_unitigs, &succ, genome_iter)?;

    let old_size = std::fs::metadata(v2_path)?.len();
    let new_size = std::fs::metadata(v3_path)?.len();
    Ok(MigrationV3Stats {
        num_genomes: num_genomes as u64,
        num_unitigs,
        num_edges,
        old_size,
        new_size,
    })
}

// ---------------------------------------------------------------------------
// Mmap reader
// ---------------------------------------------------------------------------

/// Mmap-backed view over a v3 paths file.
///
/// The per-genome blob section stays mmap'd (each `get_path` reads one
/// blob sequentially). The successor CSR table and the genome offset
/// table are pulled into RAM at open time (~49 MB/shard total): they are
/// random-accessed ~1.3M times per path decode, and serving that through
/// the mmap faults scattered 4 KB pages from (often networked) storage —
/// a measured ~2x search slowdown. 49 MB resident is negligible against
/// the laptop RAM budget.
pub struct MmapPathIndexV3 {
    mmap: Mmap,
    num_genomes: u64,
    num_nodes: u64,
    succ_offsets: Vec<u64>,
    succ_values: Vec<u32>,
    genome_offsets: Vec<u64>,
}

impl MmapPathIndexV3 {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {path:?}"))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("mmap {path:?}"))?;
        if mmap.len() < HEADER_SIZE as usize {
            bail!("file too small for v3 header");
        }
        if &mmap[..8] != MAGIC {
            bail!("not a v3 paths file (magic mismatch)");
        }
        let version = u32::from_le_bytes(mmap[8..12].try_into().unwrap());
        if version != VERSION {
            bail!("unsupported v3 version: {version}");
        }
        let num_genomes = u64::from_le_bytes(mmap[16..24].try_into().unwrap());
        let num_unitigs = u64::from_le_bytes(mmap[24..32].try_into().unwrap());
        let succ_offsets_pos = u64::from_le_bytes(mmap[32..40].try_into().unwrap()) as usize;
        let succ_values_pos = u64::from_le_bytes(mmap[40..48].try_into().unwrap()) as usize;
        let genome_table_pos = u64::from_le_bytes(mmap[48..56].try_into().unwrap()) as usize;
        let num_nodes = num_unitigs * 2;

        let table_end = genome_table_pos + (num_genomes as usize + 1) * 8;
        if mmap.len() < table_end || succ_values_pos < succ_offsets_pos
            || genome_table_pos < succ_values_pos
        {
            bail!("v3 file truncated or malformed within tables");
        }

        // Pull the randomly-accessed tables into RAM (one sequential read each).
        let read_u64_table = |start: usize, end: usize| -> Result<Vec<u64>> {
            if (end - start) % 8 != 0 {
                bail!("v3 u64 table not 8-byte aligned");
            }
            Ok((start..end)
                .step_by(8)
                .map(|p| u64::from_le_bytes(mmap[p..p + 8].try_into().unwrap()))
                .collect())
        };
        let succ_offsets = read_u64_table(succ_offsets_pos, succ_values_pos)?;
        if succ_offsets.len() != num_nodes as usize + 1 {
            bail!(
                "v3 succ_offsets length {} != num_nodes+1 {}",
                succ_offsets.len(),
                num_nodes + 1
            );
        }
        let succ_values: Vec<u32> = (succ_values_pos..genome_table_pos)
            .step_by(4)
            .map(|p| u32::from_le_bytes(mmap[p..p + 4].try_into().unwrap()))
            .collect();
        let genome_offsets = read_u64_table(genome_table_pos, table_end)?;

        Ok(Self {
            mmap,
            num_genomes,
            num_nodes,
            succ_offsets,
            succ_values,
            genome_offsets,
        })
    }

    pub fn num_genomes(&self) -> u64 {
        self.num_genomes
    }

    #[inline]
    fn succ_offset(&self, node: u64) -> u64 {
        self.succ_offsets[node as usize]
    }

    #[inline]
    fn succ_value(&self, idx: u64) -> u32 {
        self.succ_values[idx as usize]
    }

    #[inline]
    fn genome_offset(&self, i: u64) -> u64 {
        self.genome_offsets[i as usize]
    }

    /// Decode genome `genome_id`'s path. Reconstructs the same `GenomePath`
    /// a v2 reader would have produced.
    pub fn get_path(&self, genome_id: u32) -> Result<Option<GenomePath>> {
        let gid = u64::from(genome_id);
        if gid >= self.num_genomes {
            return Ok(None);
        }
        let start = self.genome_offset(gid) as usize;
        let end = self.genome_offset(gid + 1) as usize;
        if end > self.mmap.len() || start > end {
            bail!("corrupt genome offset table for genome {gid}");
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
        if num_steps >= 1 {
            let (first_node, np) = read_varint(blob, p)?;
            p = np;
            let (first_off, np) = read_varint(blob, p)?;
            p = np;
            let mut node = first_node as u32;
            let mut offset = first_off;
            steps.push(PathStep {
                unitig_id: node >> 1,
                is_reverse: (node & 1) == 1,
                genome_offset: offset,
            });
            for _ in 1..num_steps {
                let (edge_index, np) = read_varint(blob, p)?;
                p = np;
                let (delta, np) = read_varint(blob, p)?;
                p = np;
                if u64::from(node) >= self.num_nodes {
                    bail!("genome {gid}: node {node} out of range");
                }
                let lo = self.succ_offset(u64::from(node));
                let hi = self.succ_offset(u64::from(node) + 1);
                if edge_index >= hi - lo {
                    bail!(
                        "genome {gid}: edge_index {edge_index} out of range \
                         (node {node} has {} successors)",
                        hi - lo
                    );
                }
                node = self.succ_value(lo + edge_index);
                offset = offset
                    .checked_add(delta)
                    .ok_or_else(|| anyhow::anyhow!("genome {gid}: offset overflow"))?;
                steps.push(PathStep {
                    unitig_id: node >> 1,
                    is_reverse: (node & 1) == 1,
                    genome_offset: offset,
                });
            }
        }

        Ok(Some(GenomePath {
            genome_id,
            genome_name: name,
            genome_length,
            steps,
        }))
    }
}

/// True if `path` starts with the v3 magic.
pub fn is_v3(path: &Path) -> Result<bool> {
    let mut f = match File::open(path) {
        Ok(f) => f,
        Err(_) => return Ok(false),
    };
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
        // Deterministic walk: unitig ids cycle, offsets strictly increasing.
        let mut steps = Vec::new();
        let mut off: u64 = 0;
        for i in 0..n_steps {
            steps.push(PathStep {
                unitig_id: (i * 7 + 3) % 97,
                is_reverse: i % 3 == 0,
                genome_offset: off,
            });
            off += u64::from(i % 13) + 1;
        }
        GenomePath {
            genome_id: id,
            genome_name: name.to_owned(),
            genome_length: off,
            steps,
        }
    }

    fn assert_same(a: &GenomePath, b: &GenomePath) {
        assert_eq!(a.genome_id, b.genome_id);
        assert_eq!(a.genome_name, b.genome_name);
        assert_eq!(a.genome_length, b.genome_length);
        assert_eq!(a.steps.len(), b.steps.len(), "step count");
        for (i, (x, y)) in a.steps.iter().zip(b.steps.iter()).enumerate() {
            assert_eq!(x.unitig_id, y.unitig_id, "step {i} unitig_id");
            assert_eq!(x.is_reverse, y.is_reverse, "step {i} is_reverse");
            assert_eq!(x.genome_offset, y.genome_offset, "step {i} genome_offset");
        }
    }

    #[test]
    fn roundtrip_multiple_genomes() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let path = tmp.path().join("v3.bin");
        let genomes = vec![
            make_genome(0, "alpha", 1),
            make_genome(1, "bravo", 250),
            make_genome(2, "charlie_long_name", 4096),
        ];
        write_v3_from_paths(&path, &genomes)?;
        assert!(is_v3(&path)?);

        let r = MmapPathIndexV3::open(&path)?;
        assert_eq!(r.num_genomes(), 3);
        for g in &genomes {
            let got = r.get_path(g.genome_id)?.unwrap();
            assert_same(g, &got);
        }
        assert!(r.get_path(99)?.is_none());
        Ok(())
    }

    #[test]
    fn empty_path_roundtrips() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let path = tmp.path().join("empty_v3.bin");
        let g = GenomePath {
            genome_id: 0,
            genome_name: "lonely".to_owned(),
            genome_length: 0,
            steps: vec![],
        };
        write_v3_from_paths(&path, std::slice::from_ref(&g))?;
        let r = MmapPathIndexV3::open(&path)?;
        let got = r.get_path(0)?.unwrap();
        assert_same(&g, &got);
        Ok(())
    }

    #[test]
    fn migrate_from_v2_matches() -> Result<()> {
        use crate::index::paths_v2::{MmapPathIndex, PathV2Writer};
        let tmp = tempfile::tempdir()?;
        let v2 = tmp.path().join("paths_v2.bin");
        let v3 = tmp.path().join("paths_v3.bin");

        let genomes = vec![
            make_genome(0, "g0", 5),
            make_genome(1, "g1", 1000),
            make_genome(2, "g2", 333),
        ];
        let mut w = PathV2Writer::create(&v2, genomes.len() as u64)?;
        for g in &genomes {
            w.write_genome(g)?;
        }
        w.finish()?;

        let stats = migrate_v2_to_v3(&v2, &v3)?;
        assert_eq!(stats.num_genomes, 3);

        // v3 must decode byte-identically to v2.
        let r2 = MmapPathIndex::open(&v2)?;
        let r3 = MmapPathIndexV3::open(&v3)?;
        for g in &genomes {
            let from_v2 = r2.get_path(g.genome_id)?.unwrap();
            let from_v3 = r3.get_path(g.genome_id)?.unwrap();
            assert_same(&from_v2, &from_v3);
        }
        Ok(())
    }
}
