//! `paths.bin` v4 — graph-edge encoding with per-genome checkpoints.
//!
//! v3 reconstructs every step by chasing edges through a shared CSR
//! successor table. The table is ~60 MB/shard — larger than L3 cache — so
//! decoding a 650 k-step path costs ~100 ms (1.3M random table reads × 80 ns
//! cache-miss to DRAM). v3 was ~30× slower per `get_path` than v2.
//!
//! v4 keeps the disk savings of v3 (same CSR successor table; same
//! edge-index + delta step stream) and adds a per-genome **checkpoint
//! table** so callers can decode only the steps inside a small genome-offset
//! window. `direct_align` only ever needs steps near the query anchor —
//! typically a few hundred steps out of 650 k — so `iter_window` decodes
//! ~1 ms instead of ~100 ms.
//!
//! ## Wire format (`DRGNPTH4`)
//!
//! ```text
//! Header (56 bytes):
//!   +0   magic[8]            = b"DRGNPTH4"
//!   +8   version: u32        = 1
//!   +12  _pad: u32           = 0
//!   +16  num_genomes: u64
//!   +24  num_unitigs: u64
//!   +32  succ_offsets_pos: u64
//!   +40  succ_values_pos: u64
//!   +48  genome_table_pos: u64
//!
//! Successor table (CSR) — node = (unitig_id << 1) | is_reverse:
//!   succ_offsets: [u64; 2*num_unitigs + 1]
//!   succ_values:  [u32; total_edges]
//!
//! Genome offset table: [u64; num_genomes + 1]
//!
//! Per-genome blob (LEB128 varints):
//!   varint name_length, bytes name
//!   varint genome_length
//!   varint num_steps
//!   varint num_checkpoints
//!   varint stream_size                (bytes; lets the reader find the end)
//!   repeated num_checkpoints:
//!     varint step_idx                 (0-based, in [0, num_steps))
//!     varint genome_offset_at_step    (the step's genome_offset)
//!     varint node_at_step             ((unitig_id<<1)|is_reverse at step_idx)
//!     varint stream_pos               (byte offset within step stream where
//!                                       step `step_idx` BEGINS — for step 0
//!                                       this is the first_node varint; for
//!                                       step >0 this is the edge_index varint)
//!   step stream (stream_size bytes):
//!     if num_steps >= 1:
//!       varint first_node
//!       varint first_genome_offset
//!     repeated (num_steps - 1):
//!       varint edge_index
//!       varint delta_offset
//! ```
//!
//! The checkpoint at step_idx records the STATE entering that step (its node
//! and genome_offset). `stream_pos` points to the byte where that step's
//! record starts. So to resume at a checkpoint: seek to stream + stream_pos,
//! set node = checkpoint.node, offset = checkpoint.genome_offset, and that
//! IS step `step_idx`. To decode step_idx+1: read (edge, delta) following
//! the current node and accumulate.
//!
//! By convention checkpoint[0] is always step 0 (the genome's first step,
//! whose record is the first_node + first_offset pair).

use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use crate::index::paths::{GenomePath, PathStep};

pub const MAGIC: &[u8; 8] = b"DRGNPTH4";
pub const VERSION: u32 = 1;
pub const HEADER_SIZE: u64 = 56;

/// Checkpoint every N steps. ~635 checkpoints per saureus genome (~650 k
/// steps), so a query window covers ≤1024 step decodes instead of 650 k —
/// ~3 orders of magnitude faster.
pub const CHECKPOINT_INTERVAL: u32 = 1024;

#[inline]
fn node_of(unitig_id: u32, is_reverse: bool) -> u32 {
    (unitig_id << 1) | (is_reverse as u32)
}

// ---------------------------------------------------------------------------
// LEB128 varint (same as v2/v3)
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
// Encoding helpers
// ---------------------------------------------------------------------------

/// One row of a genome's checkpoint table.
#[derive(Clone, Copy, Debug)]
struct Checkpoint {
    step_idx: u32,
    genome_offset: u64,
    node: u32,
    stream_pos: u64, // byte offset within the step stream
}

/// Encode one genome's blob using a borrowed CSR table.
fn encode_blob_with_table(
    path: &GenomePath,
    succ_offsets: &[u64],
    succ_values: &[u32],
) -> Result<Vec<u8>> {
    // ---- Build the step stream and record checkpoint positions ----
    let mut stream: Vec<u8> = Vec::with_capacity(path.steps.len() * 2);
    let mut checkpoints: Vec<Checkpoint> = Vec::new();

    if !path.steps.is_empty() {
        let first = &path.steps[0];
        let first_node = node_of(first.unitig_id, first.is_reverse);
        // Checkpoint 0 = step 0 at stream position 0.
        checkpoints.push(Checkpoint {
            step_idx: 0,
            genome_offset: first.genome_offset,
            node: first_node,
            stream_pos: 0,
        });
        write_varint(&mut stream, u64::from(first_node));
        write_varint(&mut stream, first.genome_offset);

        let mut prev_node = first_node;
        let mut prev_off = first.genome_offset;
        for (i, s) in path.steps[1..].iter().enumerate() {
            let step_idx = (i + 1) as u32;
            if step_idx % CHECKPOINT_INTERVAL == 0 {
                // Record state ENTERING step_idx; stream_pos is the byte where
                // step_idx's (edge, delta) record will begin.
                let stream_pos_for_this_step = stream.len() as u64;
                // But the checkpoint records the state AT step_idx, which we
                // only know after decoding (edge, delta). We need the
                // checkpoint's stored node/offset to BE step_idx's node/offset,
                // and stream_pos to point AFTER step_idx's record so the
                // decoder doesn't re-read it. Reorder: emit the (edge, delta)
                // FIRST, compute the new node/offset, THEN record the
                // checkpoint pointing to the next read position.
                //
                // Practical choice: define stream_pos = position of (edge, delta)
                // for step_idx, and checkpoint records (prev_node, prev_off)
                // because that's what the decoder needs to seek to and resume.
                // Adjust: store (step_idx, prev_off, prev_node, stream_pos_for_this_step).
                checkpoints.push(Checkpoint {
                    step_idx,
                    genome_offset: prev_off,
                    node: prev_node,
                    stream_pos: stream_pos_for_this_step,
                });
            }
            let target = node_of(s.unitig_id, s.is_reverse);
            let lo = succ_offsets[prev_node as usize] as usize;
            let hi = succ_offsets[prev_node as usize + 1] as usize;
            let slice = &succ_values[lo..hi];
            let edge_index = slice.binary_search(&target).map_err(|_| {
                anyhow::anyhow!("transition {prev_node}->{target} absent from successor table")
            })?;
            write_varint(&mut stream, edge_index as u64);
            let delta = s.genome_offset.checked_sub(prev_off).ok_or_else(|| {
                anyhow::anyhow!("non-monotonic genome_offset in genome {}", path.genome_id)
            })?;
            write_varint(&mut stream, delta);
            prev_node = target;
            prev_off = s.genome_offset;
        }
    }

    // ---- Assemble the blob: header + checkpoint table + stream ----
    let mut blob: Vec<u8> = Vec::with_capacity(32 + checkpoints.len() * 12 + stream.len());
    write_varint(&mut blob, path.genome_name.len() as u64);
    blob.extend_from_slice(path.genome_name.as_bytes());
    write_varint(&mut blob, path.genome_length);
    write_varint(&mut blob, path.steps.len() as u64);
    write_varint(&mut blob, checkpoints.len() as u64);
    write_varint(&mut blob, stream.len() as u64);
    for cp in &checkpoints {
        write_varint(&mut blob, u64::from(cp.step_idx));
        write_varint(&mut blob, cp.genome_offset);
        write_varint(&mut blob, u64::from(cp.node));
        write_varint(&mut blob, cp.stream_pos);
    }
    blob.extend_from_slice(&stream);
    Ok(blob)
}

// ---------------------------------------------------------------------------
// Writer (in-memory)
// ---------------------------------------------------------------------------

/// Write a complete v4 file from a slice of genome paths.
pub fn write_v4_from_paths(path: &Path, genomes: &[GenomePath]) -> Result<()> {
    let num_genomes = genomes.len() as u64;
    let max_unitig = genomes
        .iter()
        .flat_map(|g| g.steps.iter())
        .map(|s| s.unitig_id)
        .max()
        .unwrap_or(0);
    let num_unitigs = u64::from(max_unitig) + 1;
    let num_nodes = num_unitigs * 2;

    // Build CSR successor table from observed transitions.
    let mut succ_sets: HashMap<u32, BTreeSet<u32>> = HashMap::new();
    for g in genomes {
        for w in g.steps.windows(2) {
            let a = node_of(w[0].unitig_id, w[0].is_reverse);
            let b = node_of(w[1].unitig_id, w[1].is_reverse);
            succ_sets.entry(a).or_default().insert(b);
        }
    }
    let (succ_offsets, succ_values) = freeze_csr(&succ_sets, num_nodes);
    drop(succ_sets);

    write_v4_file(
        path,
        num_genomes,
        num_unitigs,
        &succ_offsets,
        &succ_values,
        genomes.iter().cloned().map(Ok),
    )
}

fn freeze_csr(succ: &HashMap<u32, BTreeSet<u32>>, num_nodes: u64) -> (Vec<u64>, Vec<u32>) {
    let mut offsets: Vec<u64> = Vec::with_capacity(num_nodes as usize + 1);
    let mut values: Vec<u32> = Vec::new();
    offsets.push(0);
    for node in 0..num_nodes as u32 {
        if let Some(set) = succ.get(&node) {
            values.extend(set.iter().copied());
        }
        offsets.push(values.len() as u64);
    }
    (offsets, values)
}

/// Core writer: header + CSR table + offset table + per-genome blobs.
fn write_v4_file(
    path: &Path,
    num_genomes: u64,
    num_unitigs: u64,
    succ_offsets: &[u64],
    succ_values: &[u32],
    genomes: impl Iterator<Item = Result<GenomePath>>,
) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {path:?}"))?;
    let mut w = BufWriter::with_capacity(16 * 1024 * 1024, file);

    let succ_offsets_pos = HEADER_SIZE;
    let succ_values_pos = succ_offsets_pos + (succ_offsets.len() as u64) * 8;
    let genome_table_pos = succ_values_pos + (succ_values.len() as u64) * 4;
    let body_start = genome_table_pos + (num_genomes + 1) * 8;

    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;
    w.write_all(&0u32.to_le_bytes())?; // pad
    w.write_all(&num_genomes.to_le_bytes())?;
    w.write_all(&num_unitigs.to_le_bytes())?;
    w.write_all(&succ_offsets_pos.to_le_bytes())?;
    w.write_all(&succ_values_pos.to_le_bytes())?;
    w.write_all(&genome_table_pos.to_le_bytes())?;
    for &o in succ_offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for &v in succ_values {
        w.write_all(&v.to_le_bytes())?;
    }
    w.write_all(&vec![0u8; (num_genomes as usize + 1) * 8])?;

    let mut genome_offsets: Vec<u64> = Vec::with_capacity(num_genomes as usize + 1);
    genome_offsets.push(body_start);
    let mut written = body_start;
    let mut count = 0u64;
    for g in genomes {
        let g = g?;
        let blob = encode_blob_with_table(&g, succ_offsets, succ_values)?;
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

// ---------------------------------------------------------------------------
// Migration v3 -> v4 (re-use the successor table, add checkpoints)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct MigrationV4Stats {
    pub num_genomes: u64,
    pub num_unitigs: u64,
    pub old_size: u64,
    pub new_size: u64,
}

/// Migrate a v3 paths.bin to v4. Two streaming passes:
/// pass 1 collects every transition observed across all genomes (the v3
/// successor table is *already* in the source file but we re-derive to be
/// fully self-contained — pass 1 is fast); pass 2 re-encodes each genome
/// with checkpoints. Only one `GenomePath` is resident at a time.
pub fn migrate_v3_to_v4(v3_path: &Path, v4_path: &Path) -> Result<MigrationV4Stats> {
    use crate::index::paths_v3::MmapPathIndexV3;

    let src = MmapPathIndexV3::open(v3_path)
        .with_context(|| format!("open v3 source {v3_path:?}"))?;
    let num_genomes = src.num_genomes();

    let mut succ_sets: HashMap<u32, BTreeSet<u32>> = HashMap::new();
    let mut max_unitig: u32 = 0;
    for gid in 0..num_genomes {
        let p = src
            .get_path(gid as u32)?
            .ok_or_else(|| anyhow::anyhow!("v3 genome {gid} missing"))?;
        for w in p.steps.windows(2) {
            let a = node_of(w[0].unitig_id, w[0].is_reverse);
            let b = node_of(w[1].unitig_id, w[1].is_reverse);
            succ_sets.entry(a).or_default().insert(b);
        }
        for s in &p.steps {
            max_unitig = max_unitig.max(s.unitig_id);
        }
    }
    let num_unitigs = u64::from(max_unitig) + 1;
    let num_nodes = num_unitigs * 2;
    let (succ_offsets, succ_values) = freeze_csr(&succ_sets, num_nodes);
    drop(succ_sets);

    let genome_iter = (0..num_genomes).map(|gid| {
        src.get_path(gid as u32)?
            .ok_or_else(|| anyhow::anyhow!("v3 genome {gid} missing in pass 2"))
    });
    write_v4_file(
        v4_path,
        num_genomes as u64,
        num_unitigs,
        &succ_offsets,
        &succ_values,
        genome_iter,
    )?;

    Ok(MigrationV4Stats {
        num_genomes: num_genomes as u64,
        num_unitigs,
        old_size: std::fs::metadata(v3_path)?.len(),
        new_size: std::fs::metadata(v4_path)?.len(),
    })
}

// ---------------------------------------------------------------------------
// Mmap reader
// ---------------------------------------------------------------------------

pub struct MmapPathIndexV4 {
    mmap: Mmap,
    num_genomes: u64,
    succ_offsets: Vec<u64>,
    succ_values: Vec<u32>,
    genome_offsets: Vec<u64>,
}

impl MmapPathIndexV4 {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {path:?}"))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("mmap {path:?}"))?;
        if mmap.len() < HEADER_SIZE as usize {
            bail!("file too small for v4 header");
        }
        if &mmap[..8] != MAGIC {
            bail!("not a v4 paths file (magic mismatch)");
        }
        let version = u32::from_le_bytes(mmap[8..12].try_into().unwrap());
        if version != VERSION {
            bail!("unsupported v4 version: {version}");
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
            bail!("v4 file truncated or malformed within tables");
        }

        let read_u64_table = |start: usize, end: usize| -> Vec<u64> {
            (start..end)
                .step_by(8)
                .map(|p| u64::from_le_bytes(mmap[p..p + 8].try_into().unwrap()))
                .collect()
        };
        let succ_offsets = read_u64_table(succ_offsets_pos, succ_values_pos);
        if succ_offsets.len() != num_nodes as usize + 1 {
            bail!(
                "v4 succ_offsets length {} != num_nodes+1 {}",
                succ_offsets.len(),
                num_nodes + 1
            );
        }
        let succ_values: Vec<u32> = (succ_values_pos..genome_table_pos)
            .step_by(4)
            .map(|p| u32::from_le_bytes(mmap[p..p + 4].try_into().unwrap()))
            .collect();
        let genome_offsets = read_u64_table(genome_table_pos, table_end);

        Ok(Self {
            mmap,
            num_genomes,
            succ_offsets,
            succ_values,
            genome_offsets,
        })
    }

    pub fn num_genomes(&self) -> u64 {
        self.num_genomes
    }

    /// Return the number of distinct successors of unitig `uid` in the shared
    /// CSR successor table.  O(1): reads two entries from succ_offsets.
    /// Returns 0 if `uid` is out of range or the v4 index has no successor table.
    pub fn unitig_successor_degree(&self, unitig_id: u32) -> u32 {
        let uid = unitig_id as usize;
        let n_unitigs = self.succ_offsets.len().saturating_sub(1);
        if uid >= n_unitigs {
            return 0;
        }
        let lo = self.succ_offsets[uid];
        let hi = self.succ_offsets[uid + 1];
        hi.saturating_sub(lo) as u32
    }

    /// Decode the entire genome path. Provided for compatibility; for hot
    /// paths use [`iter_window`].
    pub fn get_path(&self, genome_id: u32) -> Result<Option<GenomePath>> {
        let gid = u64::from(genome_id);
        if gid >= self.num_genomes {
            return Ok(None);
        }
        let (name, genome_length, steps) = self.decode_genome(gid, 0, u64::MAX)?;
        Ok(Some(GenomePath {
            genome_id,
            genome_name: name,
            genome_length,
            steps,
        }))
    }

    /// Decode only the steps with `genome_offset < end` and seek straight to
    /// the checkpoint closest before `start`. Returns owned `PathStep`s in
    /// order. Drops the genome header/metadata.
    pub fn iter_window(
        &self,
        genome_id: u32,
        start: u64,
        end: u64,
    ) -> Result<Vec<PathStep>> {
        let gid = u64::from(genome_id);
        if gid >= self.num_genomes {
            return Ok(Vec::new());
        }
        let (_n, _gl, steps) = self.decode_genome(gid, start, end)?;
        Ok(steps)
    }

    /// Lookup the first PathStep matching each `wanted` unitig_id by walking
    /// the genome's step stream from the start, with **early-exit** as soon
    /// as every wanted unitig has been found. For real alignments the seeds
    /// cluster near one genome region, so this terminates well before the
    /// full ~650k-step decode that `get_path` performs.
    ///
    /// Returns an empty map if `wanted` is empty or the genome is absent.
    pub fn find_unitig_steps(
        &self,
        genome_id: u32,
        wanted: &std::collections::HashSet<u32>,
    ) -> Result<HashMap<u32, PathStep>> {
        let mut found: HashMap<u32, PathStep> = HashMap::with_capacity(wanted.len());
        if wanted.is_empty() {
            return Ok(found);
        }
        let gid = u64::from(genome_id);
        if gid >= self.num_genomes {
            return Ok(found);
        }
        let need = wanted.len();
        self.walk_genome(gid, 0, u64::MAX, |step| {
            if wanted.contains(&step.unitig_id) {
                found.entry(step.unitig_id).or_insert(step);
                if found.len() == need {
                    return false; // early-exit: all wanted unitigs located
                }
            }
            true
        })?;
        Ok(found)
    }

    /// Read just the genome's name and total length without decoding any
    /// steps. Useful for direct_align which now anchors via
    /// [`find_unitig_steps`] + [`iter_window`] and only needs metadata.
    pub fn genome_meta(&self, genome_id: u32) -> Result<Option<(String, u64)>> {
        let gid = u64::from(genome_id);
        if gid >= self.num_genomes {
            return Ok(None);
        }
        let lo = self.genome_offsets[gid as usize] as usize;
        let hi = self.genome_offsets[gid as usize + 1] as usize;
        if hi > self.mmap.len() || lo > hi {
            bail!("corrupt genome offset table for genome {gid}");
        }
        let blob = &self.mmap[lo..hi];
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
        let (genome_length, _np) = read_varint(blob, p)?;
        Ok(Some((name, genome_length)))
    }

    /// Shared decoder used by `get_path` (start=0, end=MAX), `iter_window`,
    /// and (via [`walk_genome`]) `find_unitig_steps`. Returns
    /// `(genome_name, genome_length, steps)`.
    ///
    /// The checkpoint table embedded in the blob is searched for the entry
    /// with `genome_offset <= start`. Decoding starts from there. For
    /// `start = 0` this is checkpoint 0 (the first step), so the full path
    /// is decoded — identical output to v3 `get_path`.
    fn decode_genome(
        &self,
        gid: u64,
        start: u64,
        end: u64,
    ) -> Result<(String, u64, Vec<PathStep>)> {
        let mut steps: Vec<PathStep> = Vec::new();
        let (name, genome_length) = self.walk_genome(gid, start, end, |s| {
            steps.push(s);
            true
        })?;
        Ok((name, genome_length, steps))
    }

    /// Core walker: from the latest checkpoint at or before `start`, emit
    /// every step whose `genome_offset < end` via `emit(step) -> bool`.
    /// Walking stops when `emit` returns false, `step.genome_offset >= end`,
    /// or the path ends. Returns `(name, genome_length)`.
    fn walk_genome<F: FnMut(PathStep) -> bool>(
        &self,
        gid: u64,
        start: u64,
        end: u64,
        mut emit: F,
    ) -> Result<(String, u64)> {
        let lo = self.genome_offsets[gid as usize] as usize;
        let hi = self.genome_offsets[gid as usize + 1] as usize;
        if hi > self.mmap.len() || lo > hi {
            bail!("corrupt genome offset table for genome {gid}");
        }
        let blob = &self.mmap[lo..hi];

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
        let (num_checkpoints, np) = read_varint(blob, p)?;
        p = np;
        let (stream_size, np) = read_varint(blob, p)?;
        p = np;

        if num_steps == 0 {
            return Ok((name, genome_length));
        }
        if num_checkpoints == 0 {
            bail!("genome {gid}: non-empty path with zero checkpoints");
        }

        // Decode checkpoint table.
        struct CP {
            step_idx: u32,
            genome_offset: u64,
            node: u32,
            stream_pos: u64,
        }
        let mut ckpts: Vec<CP> = Vec::with_capacity(num_checkpoints as usize);
        for _ in 0..num_checkpoints {
            let (si, np) = read_varint(blob, p)?;
            p = np;
            let (go, np) = read_varint(blob, p)?;
            p = np;
            let (nd, np) = read_varint(blob, p)?;
            p = np;
            let (sp, np) = read_varint(blob, p)?;
            p = np;
            ckpts.push(CP {
                step_idx: si as u32,
                genome_offset: go,
                node: nd as u32,
                stream_pos: sp,
            });
        }
        let stream_start = p;
        let stream_end = stream_start + stream_size as usize;
        if stream_end > blob.len() {
            bail!("genome {gid}: stream overruns blob");
        }

        // Find the latest checkpoint with genome_offset <= start.
        // partition_point returns the first idx where pred is FALSE; we want
        // the last where it's TRUE.
        let cp_idx = ckpts
            .partition_point(|c| c.genome_offset <= start)
            .saturating_sub(1);
        let cp = &ckpts[cp_idx];

        // The checkpoint's `stream_pos` is where step `step_idx`'s record
        // starts. For step 0 (the first step) this is the first_node/first_off
        // pair; for step >0 this is the (edge, delta) pair to step_idx.
        let mut pp = stream_start + cp.stream_pos as usize;
        let mut node;
        let mut offset;
        let mut step_idx;

        if cp.step_idx == 0 {
            // Read first step from full encoding.
            let (fn_, np) = read_varint(blob, pp)?;
            pp = np;
            let (fo, np) = read_varint(blob, pp)?;
            pp = np;
            node = fn_ as u32;
            offset = fo;
            step_idx = 0u32;
        } else {
            // Resume from the recorded (prev_node, prev_off) and consume the
            // (edge, delta) at stream_pos to advance to step_idx.
            node = cp.node;
            offset = cp.genome_offset;
            // To reach step_idx we'd consume one (edge, delta) — but the
            // checkpoint already stores (node, offset) AT step_idx, with
            // stream_pos pointing to step_idx's record. So actually the state
            // recorded IS step_idx; stream_pos is where step_idx's data lives.
            // The "checkpoint records prev_node/prev_off" path was the
            // build-time view. At read time we treat the checkpoint as: "step
            // step_idx has node=cp.node, offset=cp.genome_offset, and its
            // record is at stream_pos" — same as step 0's special encoding.
            //
            // The writer's chosen layout: for step_idx >= 1 the record at
            // stream_pos is (edge, delta) that takes (prev_node, prev_off) TO
            // (node, offset). The checkpoint stores prev_node/prev_off so we
            // must apply the record to obtain step_idx's data.
            let (edge, np) = read_varint(blob, pp)?;
            pp = np;
            let (delta, np) = read_varint(blob, pp)?;
            pp = np;
            let slo = self.succ_offsets[node as usize] as usize;
            let shi = self.succ_offsets[node as usize + 1] as usize;
            if (edge as usize) >= shi - slo {
                bail!(
                    "genome {gid} ckpt {} edge {} out of range (node {} has {} successors)",
                    cp.step_idx,
                    edge,
                    node,
                    shi - slo
                );
            }
            node = self.succ_values[slo + edge as usize];
            offset = offset
                .checked_add(delta)
                .ok_or_else(|| anyhow::anyhow!("offset overflow genome {gid}"))?;
            step_idx = cp.step_idx;
        }

        // Emit the resumed step IF in window. Stop if the caller says so.
        if offset < end
            && !emit(PathStep {
                unitig_id: node >> 1,
                is_reverse: (node & 1) == 1,
                genome_offset: offset,
            })
        {
            return Ok((name, genome_length));
        }
        // Walk forward.
        step_idx += 1;
        while step_idx < num_steps as u32 && pp < stream_end {
            let (edge, np) = read_varint(blob, pp)?;
            pp = np;
            let (delta, np) = read_varint(blob, pp)?;
            pp = np;
            let slo = self.succ_offsets[node as usize] as usize;
            let shi = self.succ_offsets[node as usize + 1] as usize;
            if (edge as usize) >= shi - slo {
                bail!(
                    "genome {gid} step {} edge {} out of range",
                    step_idx,
                    edge
                );
            }
            node = self.succ_values[slo + edge as usize];
            offset = offset
                .checked_add(delta)
                .ok_or_else(|| anyhow::anyhow!("offset overflow genome {gid}"))?;
            if offset >= end {
                break;
            }
            if !emit(PathStep {
                unitig_id: node >> 1,
                is_reverse: (node & 1) == 1,
                genome_offset: offset,
            }) {
                return Ok((name, genome_length));
            }
            step_idx += 1;
        }

        Ok((name, genome_length))
    }
}

/// True if `path` starts with the v4 magic.
pub fn is_v4(path: &Path) -> Result<bool> {
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
            assert_eq!(x.genome_offset, y.genome_offset, "step {i} offset");
        }
    }

    #[test]
    fn get_path_roundtrips() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let p = tmp.path().join("v4.bin");
        let genomes = vec![
            make_genome(0, "g0", 1),
            make_genome(1, "g1_short", 500),
            // Span multiple checkpoint intervals.
            make_genome(2, "g2_long", CHECKPOINT_INTERVAL * 3 + 17),
        ];
        write_v4_from_paths(&p, &genomes)?;
        assert!(is_v4(&p)?);
        let r = MmapPathIndexV4::open(&p)?;
        assert_eq!(r.num_genomes(), genomes.len() as u64);
        for g in &genomes {
            let got = r.get_path(g.genome_id)?.unwrap();
            assert_same(g, &got);
        }
        Ok(())
    }

    #[test]
    fn iter_window_matches_full_decode() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let p = tmp.path().join("v4w.bin");
        // Long path spanning several checkpoint intervals so iter_window
        // actually exercises the seek.
        let g = make_genome(0, "long", CHECKPOINT_INTERVAL * 5 + 123);
        write_v4_from_paths(&p, std::slice::from_ref(&g))?;
        let r = MmapPathIndexV4::open(&p)?;

        let full = r.get_path(0)?.unwrap();
        // For a variety of windows, iter_window must yield the SAME steps as
        // filtering `full` by genome_offset < end.
        let total_len = g.genome_length;
        let cases = [
            (0u64, total_len),
            (0, total_len / 10),
            (total_len / 2, total_len),
            (total_len / 3, total_len * 2 / 3),
            (total_len - 50, total_len),
            (0, 50),
        ];
        for (start, end) in cases {
            let win = r.iter_window(0, start, end)?;
            if end == 0 {
                assert!(win.is_empty(), "empty window for end=0");
                continue;
            }
            assert!(!win.is_empty(), "window [{start},{end}) unexpectedly empty");
            // iter_window seeks to the latest checkpoint <= start, so its
            // first step's offset may be ≤ start. Assert it matches the
            // contiguous tail of `full` starting at that same step.
            let first_off = win[0].genome_offset;
            let first_idx = full
                .steps
                .iter()
                .position(|s| s.genome_offset == first_off)
                .expect("first step must be in full path");
            let expected: Vec<PathStep> = full.steps[first_idx..]
                .iter()
                .take_while(|s| s.genome_offset < end)
                .copied()
                .collect();
            assert_eq!(win.len(), expected.len(), "window [{start},{end}) length");
            for (i, (a, b)) in win.iter().zip(expected.iter()).enumerate() {
                assert_eq!(a.unitig_id, b.unitig_id, "win[{i}] unitig");
                assert_eq!(a.is_reverse, b.is_reverse, "win[{i}] rev");
                assert_eq!(a.genome_offset, b.genome_offset, "win[{i}] offset");
            }
        }
        Ok(())
    }

    #[test]
    fn find_unitig_steps_returns_first_match_and_exits_early() -> Result<()> {
        use std::collections::HashSet;
        let tmp = tempfile::tempdir()?;
        let p = tmp.path().join("v4f.bin");
        // Long path so the early-exit path is exercised.
        let g = make_genome(0, "long", CHECKPOINT_INTERVAL * 3 + 17);
        write_v4_from_paths(&p, std::slice::from_ref(&g))?;
        let r = MmapPathIndexV4::open(&p)?;

        // Pick three unitig_ids that DO appear in the path (the make_genome
        // formula cycles unitig_id = (i*7 + 3) % 97 so we know all of these
        // will appear) and one that does NOT (unitig 200 > 97).
        let mut wanted: HashSet<u32> = HashSet::new();
        wanted.insert(3); // step 0
        wanted.insert(10); // step 1
        wanted.insert(17); // step 2
        let found = r.find_unitig_steps(0, &wanted)?;
        assert_eq!(found.len(), wanted.len(), "all wanted unitigs found");
        for &uid in &wanted {
            let step = found.get(&uid).expect("unitig found");
            assert_eq!(step.unitig_id, uid);
            // The match should be the FIRST occurrence in the path.
            let first = g.steps.iter().find(|s| s.unitig_id == uid).unwrap();
            assert_eq!(step.genome_offset, first.genome_offset);
            assert_eq!(step.is_reverse, first.is_reverse);
        }

        // Absent unitig — find_unitig_steps walks the whole path looking for
        // it, but never finds it.
        let mut not_found: HashSet<u32> = HashSet::new();
        not_found.insert(200);
        let res = r.find_unitig_steps(0, &not_found)?;
        assert!(res.is_empty());

        // Empty wanted set returns empty.
        let res = r.find_unitig_steps(0, &HashSet::new())?;
        assert!(res.is_empty());
        Ok(())
    }

    #[test]
    fn genome_meta_reads_without_decoding_steps() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let p = tmp.path().join("v4m.bin");
        let genomes = vec![
            make_genome(0, "alpha", 0),
            make_genome(1, "bravo_long_name", 4096),
        ];
        write_v4_from_paths(&p, &genomes)?;
        let r = MmapPathIndexV4::open(&p)?;
        let (n0, l0) = r.genome_meta(0)?.unwrap();
        assert_eq!(n0, "alpha");
        assert_eq!(l0, genomes[0].genome_length);
        let (n1, l1) = r.genome_meta(1)?.unwrap();
        assert_eq!(n1, "bravo_long_name");
        assert_eq!(l1, genomes[1].genome_length);
        assert!(r.genome_meta(99)?.is_none());
        Ok(())
    }

    #[test]
    fn migrate_v3_to_v4_matches() -> Result<()> {
        use crate::index::paths_v3::{write_v3_from_paths, MmapPathIndexV3};
        let tmp = tempfile::tempdir()?;
        let v3 = tmp.path().join("v3.bin");
        let v4 = tmp.path().join("v4.bin");
        let genomes = vec![
            make_genome(0, "a", 1),
            make_genome(1, "b", CHECKPOINT_INTERVAL * 2 + 33),
            make_genome(2, "c", 700),
        ];
        write_v3_from_paths(&v3, &genomes)?;
        let _stats = migrate_v3_to_v4(&v3, &v4)?;

        let r3 = MmapPathIndexV3::open(&v3)?;
        let r4 = MmapPathIndexV4::open(&v4)?;
        for g in &genomes {
            let a = r3.get_path(g.genome_id)?.unwrap();
            let b = r4.get_path(g.genome_id)?.unwrap();
            assert_same(&a, &b);
        }
        Ok(())
    }
}
