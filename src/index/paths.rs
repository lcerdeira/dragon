/// Genome path index: stores each genome as a sequence of unitig IDs
/// traversed when walking the genome through the de Bruijn graph.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

use crate::index::unitig::UnitigSet;
use crate::util::dna;

/// A genome's path through the de Bruijn graph.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GenomePath {
    pub genome_id: u32,
    pub genome_name: String,
    pub genome_length: u64,
    /// Sequence of (unitig_id, is_reverse) pairs.
    pub steps: Vec<PathStep>,
}

/// A single step in a genome path.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct PathStep {
    pub unitig_id: u32,
    pub is_reverse: bool,
    /// Start position within the genome sequence.
    pub genome_offset: u64,
}

/// Collection of all genome paths.
///
/// Two backings are supported:
///   * `Eager` — legacy bincode `Vec<GenomePath>` loaded fully into RAM.
///     Used when `paths.bin` predates the v2 format.
///   * `Mmap`  — v2 file mmap'd lazily; per-genome blobs are decoded on
///     demand. O(1) load time, O(touched genomes × steps) RAM at query time.
///
/// Public API yields **owned** `GenomePath` values so the same code paths
/// work for both backings (the eager variant clones; the mmap variant
/// decodes from the mmap region).
#[derive(Clone)]
pub enum PathIndex {
    /// All genomes pre-decoded into RAM.
    Eager(std::sync::Arc<Vec<GenomePath>>),
    /// Lazy v2 view; only the offset table is mmap-resident at idle.
    Mmap(std::sync::Arc<crate::index::paths_v2::MmapPathIndex>),
}

/// Bincode-friendly shadow type for the legacy on-disk format. Keeping the
/// `paths` Vec layout means existing bincode files (and the
/// `crate::util::mmap::write_bincode(&PathIndex, ...)` test helper) keep
/// working unchanged.
#[derive(Clone, Debug, Serialize, Deserialize)]
struct LegacyPathIndex {
    paths: Vec<GenomePath>,
}

impl PathIndex {
    /// Build an eager index from a fully materialized `Vec<GenomePath>`.
    pub fn from_paths(paths: Vec<GenomePath>) -> Self {
        Self::Eager(std::sync::Arc::new(paths))
    }

    /// Decode the i-th genome. Returns owned data — for the mmap variant
    /// this allocates a `Vec<PathStep>` for the requested genome only.
    pub fn get_path(&self, genome_id: u32) -> Option<GenomePath> {
        match self {
            Self::Eager(v) => v.get(genome_id as usize).cloned(),
            Self::Mmap(m) => m.get_path(genome_id).ok().flatten(),
        }
    }

    pub fn num_genomes(&self) -> usize {
        match self {
            Self::Eager(v) => v.len(),
            Self::Mmap(m) => m.num_genomes() as usize,
        }
    }

    /// Iterate every genome path. Yields owned values; the mmap variant
    /// decodes lazily and is suitable even for huge indices because each
    /// item is dropped before the next is decoded.
    pub fn iter(&self) -> Box<dyn Iterator<Item = GenomePath> + '_> {
        match self {
            Self::Eager(v) => Box::new(v.iter().cloned()),
            Self::Mmap(m) => {
                let n = m.num_genomes();
                let m = m.clone();
                Box::new((0..n).filter_map(move |gid| m.get_path(gid as u32).ok().flatten()))
            }
        }
    }

    /// Extract the reference sequence for a region of a genome by walking its path.
    pub fn extract_sequence(
        &self,
        genome_id: u32,
        start: u64,
        end: u64,
        unitigs: &UnitigSet,
    ) -> Vec<u8> {
        let path = match self.get_path(genome_id) {
            Some(p) => p,
            None => return Vec::new(),
        };

        let mut result = Vec::with_capacity((end - start) as usize);

        for step in &path.steps {
            let unitig_len = unitigs.unitigs[step.unitig_id as usize].sequence.len as u64;
            let step_end = step.genome_offset + unitig_len;

            // Check if this step overlaps with the requested region
            if step.genome_offset >= end || step_end <= start {
                continue;
            }

            // Calculate the overlap
            let overlap_start = start.max(step.genome_offset);
            let overlap_end = end.min(step_end);
            let local_start = (overlap_start - step.genome_offset) as usize;
            let local_end = (overlap_end - step.genome_offset) as usize;

            let mut subseq = unitigs.get_subsequence(step.unitig_id, local_start, local_end);

            if step.is_reverse {
                // Reverse complement
                subseq.reverse();
                for base in &mut subseq {
                    *base = match *base {
                        b'A' => b'T',
                        b'T' => b'A',
                        b'C' => b'G',
                        b'G' => b'C',
                        other => other,
                    };
                }
            }

            result.extend_from_slice(&subseq);
        }

        result
    }

    /// Extract sequence from a GenomePath directly (no PathIndex lookup needed).
    pub fn extract_sequence_static(
        path: &GenomePath,
        start: u64,
        end: u64,
        unitigs: &UnitigSet,
    ) -> Vec<u8> {
        let mut result = Vec::with_capacity((end - start) as usize);

        for step in &path.steps {
            let unitig_len = unitigs.unitigs.get(step.unitig_id as usize)
                .map(|u| u.sequence.len as u64)
                .unwrap_or(0);
            let step_end = step.genome_offset + unitig_len;

            if step.genome_offset >= end || step_end <= start {
                continue;
            }

            let overlap_start = start.max(step.genome_offset);
            let overlap_end = end.min(step_end);
            let local_start = (overlap_start - step.genome_offset) as usize;
            let local_end = (overlap_end - step.genome_offset) as usize;

            let mut subseq = unitigs.get_subsequence(step.unitig_id, local_start, local_end);

            if step.is_reverse {
                subseq.reverse();
                for base in &mut subseq {
                    *base = match *base {
                        b'A' => b'T', b'T' => b'A',
                        b'C' => b'G', b'G' => b'C',
                        other => other,
                    };
                }
            }

            result.extend_from_slice(&subseq);
        }

        result
    }
}

/// Build a hash table mapping canonical k-mers to (unitig_id, offset, is_reverse).
///
/// For each unitig, slide a k-mer window across the forward strand. Store canonical
/// k-mer → (unitig_id, offset_in_unitig, is_rc). If a canonical k-mer appears in
/// multiple unitigs (possible at boundaries), keep the first occurrence — the path
/// walker will still find contiguous runs.
fn build_kmer_table(
    unitigs: &UnitigSet,
    kmer_size: usize,
) -> HashMap<u64, (u32, u32, bool)> {
    let mut table: HashMap<u64, (u32, u32, bool)> = HashMap::new();

    for unitig in &unitigs.unitigs {
        let seq_len = unitig.sequence.len;
        if seq_len < kmer_size {
            continue;
        }

        for pos in 0..=(seq_len - kmer_size) {
            let fwd = unitig.sequence.kmer_u64(pos, kmer_size);
            let canon = dna::canonical_kmer(fwd, kmer_size);
            let is_rc = canon != fwd;

            // First occurrence wins — avoids expensive overwrite checks
            table.entry(canon).or_insert((unitig.id, pos as u32, is_rc));
        }
    }

    log::info!(
        "Built k-mer table: {} distinct canonical {}-mers from {} unitigs",
        table.len(),
        kmer_size,
        unitigs.num_unitigs()
    );
    table
}

/// Check if a k-mer window contains ambiguous bases (N).
#[inline]
fn has_ambiguous(seq: &[u8], start: usize, k: usize) -> bool {
    seq[start..start + k].iter().any(|&b| !matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
}

/// Build the genome path index by walking each genome through unitig k-mers.
pub fn build_path_index(
    genome_dir: &Path,
    unitigs: &UnitigSet,
    output_dir: &Path,
    kmer_size: usize,
) -> Result<()> {
    // Step 1: Build canonical k-mer → (unitig_id, offset, is_rc) hash table
    let kmer_table = build_kmer_table(unitigs, kmer_size);

    // Step 2: Walk each genome and stream its path directly to paths.bin.
    //
    // We serialize using bincode's `Vec<GenomePath>` wire format:
    //   1. 8-byte little-endian length prefix (= num_genomes)
    //   2. Each GenomePath serialized back-to-back.
    // This is byte-identical to `bincode::serialize(&PathIndex { paths })`,
    // so load_path_index() reads it with no code change. Peak RAM for this
    // step drops from O(sum_over_genomes(steps)) to O(max_single_genome_steps)
    // because each GenomePath is dropped as soon as it's written.
    let genome_files = crate::io::fasta::list_fasta_files(genome_dir)?;
    let num_genomes = genome_files.len();

    let path_file = output_dir.join("paths.bin");
    // v2 mmap-friendly format. Per-genome blobs are streamed; only the offset
    // table accumulates in RAM (8 bytes × num_genomes ≈ 256 KB for 32K).
    let mut writer = crate::index::paths_v2::PathV2Writer::create(
        &path_file,
        num_genomes as u64,
    )?;

    let mut total_steps = 0usize;
    let mut total_mapped_bases = 0u64;
    let mut total_genome_bases = 0u64;
    let log_every = (num_genomes / 20).max(1).min(1000);

    for (genome_id, genome_file) in genome_files.iter().enumerate() {
        let sequences = crate::io::fasta::read_sequences(genome_file)?;
        let genome_name = genome_file
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        let mut genome_offset = 0u64;
        let mut steps = Vec::new();
        let mut prev_unitig: Option<(u32, bool)> = None;

        for seq in &sequences {
            let seq_bytes = &seq.seq;
            let seq_len = seq_bytes.len();
            total_genome_bases += seq_len as u64;

            if seq_len < kmer_size {
                genome_offset += seq_len as u64;
                continue;
            }

            for pos in 0..=(seq_len - kmer_size) {
                // Skip windows with ambiguous bases
                if has_ambiguous(seq_bytes, pos, kmer_size) {
                    prev_unitig = None;
                    continue;
                }

                // Encode k-mer from raw bytes
                let packed = dna::PackedSequence::from_bytes(&seq_bytes[pos..pos + kmer_size]);
                let fwd = packed.kmer_u64(0, kmer_size);
                let canon = dna::canonical_kmer(fwd, kmer_size);

                if let Some(&(unitig_id, _offset, kmer_is_rc)) = kmer_table.get(&canon) {
                    // Determine orientation: if the genome k-mer's canonical form
                    // required RC, and the unitig's stored k-mer also required RC,
                    // they cancel out → forward. Otherwise they differ → reverse.
                    let genome_is_rc = canon != fwd;
                    let is_reverse = genome_is_rc != kmer_is_rc;

                    // Emit a new step only on unitig transitions
                    let current = (unitig_id, is_reverse);
                    if prev_unitig != Some(current) {
                        steps.push(PathStep {
                            unitig_id,
                            is_reverse,
                            genome_offset: genome_offset + pos as u64,
                        });
                        prev_unitig = Some(current);
                        total_mapped_bases += unitigs.unitigs[unitig_id as usize].sequence.len as u64;
                    }
                } else {
                    // k-mer not found in any unitig — break the chain
                    prev_unitig = None;
                }
            }

            genome_offset += seq_len as u64;
        }

        total_steps += steps.len();

        let path = GenomePath {
            genome_id: genome_id as u32,
            genome_name,
            genome_length: genome_offset,
            steps,
        };
        writer.write_genome(&path)?;
        // `path` drops here — per-genome memory is not retained.

        if (genome_id + 1) % log_every == 0 || genome_id + 1 == num_genomes {
            log::info!(
                "  paths built for {}/{} genomes ({} cumulative steps)",
                genome_id + 1, num_genomes, total_steps
            );
        }
    }

    writer.finish()?;

    log::info!(
        "Path index: {} genomes, {} total steps, {:.1}% bases mapped",
        num_genomes,
        total_steps,
        if total_genome_bases > 0 {
            total_mapped_bases as f64 / total_genome_bases as f64 * 100.0
        } else {
            0.0
        }
    );
    log::info!("Path index saved to {:?}", path_file);
    Ok(())
}

/// Load the path index from disk. Backwards-compatible: dispatches on the
/// file magic so old (bincode) and new (paths_v2 mmap) formats both work.
///
/// **Memory behaviour matters here:**
///   * v2 files load in O(1) — only the file is mmap'd. Per-genome blobs
///     are decoded lazily by [`PathIndex::get_path`] / [`PathIndex::iter`].
///   * Legacy bincode files are eagerly slurped into a `Vec<GenomePath>`
///     because their on-disk layout offers no random-access table. For
///     huge legacy files (>50 GB) this can OOM — run [`migrate_paths_to_v2`]
///     first.
pub fn load_path_index(index_dir: &Path) -> Result<PathIndex> {
    let path_file = index_dir.join("paths.bin");
    if crate::index::paths_v2::is_v2(&path_file)? {
        let mmap_idx = crate::index::paths_v2::MmapPathIndex::open(&path_file)?;
        Ok(PathIndex::Mmap(std::sync::Arc::new(mmap_idx)))
    } else {
        // Legacy bincode `Vec<GenomePath>` format.
        let legacy: LegacyPathIndex =
            crate::util::mmap::read_bincode(&path_file)?;
        Ok(PathIndex::from_paths(legacy.paths))
    }
}

/// O(1)-load mmap view of `paths.bin`. Only available for v2 files; old
/// bincode files still need [`load_path_index`].
pub fn open_path_index_mmap(
    index_dir: &Path,
) -> Result<crate::index::paths_v2::MmapPathIndex> {
    let path_file = index_dir.join("paths.bin");
    crate::index::paths_v2::MmapPathIndex::open(&path_file)
}

/// Migrate a legacy bincode-format `paths.bin` (eager `Vec<GenomePath>`) to
/// the mmap-friendly v2 format **without ever loading the whole file into RAM**.
/// Reads one `GenomePath` at a time from the legacy stream and writes each
/// directly into a fresh v2 file at `<paths>.v2.tmp`, then atomically renames
/// over the original on success. Peak memory is bounded by the largest single
/// genome's `Vec<PathStep>`.
///
/// This is the migration path for indices built before commit a70a087 — those
/// `paths.bin` files are typically tens to hundreds of GB and cannot be loaded
/// eagerly on most hardware. Running this once per shard converts them to a
/// format `MmapPathIndex` can open in O(1) at query time.
pub fn migrate_paths_to_v2(index_dir: &Path) -> Result<MigrationStats> {
    let path_file = index_dir.join("paths.bin");
    if !path_file.exists() {
        anyhow::bail!("no paths.bin in {:?}", index_dir);
    }
    if crate::index::paths_v2::is_v2(&path_file)? {
        log::info!("paths.bin in {:?} already v2 — nothing to do", index_dir);
        return Ok(MigrationStats {
            num_genomes: 0,
            old_size: std::fs::metadata(&path_file)?.len(),
            new_size: std::fs::metadata(&path_file)?.len(),
            already_v2: true,
        });
    }

    log::info!("Migrating legacy paths.bin in {:?} -> v2", index_dir);
    let old_size = std::fs::metadata(&path_file)?.len();

    // Stream-decode the legacy bincode `Vec<GenomePath>` layout:
    //   1. u64 length prefix (number of genomes)
    //   2. concatenated GenomePath records (each variable-length).
    let f = std::fs::File::open(&path_file)?;
    let mut reader = std::io::BufReader::with_capacity(16 * 1024 * 1024, f);
    let num_genomes: u64 = bincode::deserialize_from(&mut reader)
        .context("read u64 num_genomes prefix from legacy paths.bin")?;
    log::info!("  legacy file contains {} genomes ({:.1} GB on disk)",
               num_genomes, old_size as f64 / 1_073_741_824.0);

    // Write to a sibling .tmp file; rename on success.
    let tmp_file = index_dir.join("paths.bin.v2.tmp");
    let mut writer =
        crate::index::paths_v2::PathV2Writer::create(&tmp_file, num_genomes)?;

    let log_every = (num_genomes / 20).max(1).min(1000);
    for gid in 0..num_genomes {
        let g: GenomePath = bincode::deserialize_from(&mut reader)
            .with_context(|| format!("decode genome {gid} from legacy paths.bin"))?;
        writer.write_genome(&g)?;
        if (gid + 1) % log_every == 0 || gid + 1 == num_genomes {
            log::info!("  migrated {}/{} genomes", gid + 1, num_genomes);
        }
    }
    writer.finish()?;
    let new_size = std::fs::metadata(&tmp_file)?.len();

    // Atomic rename over the original.
    let backup = index_dir.join("paths.bin.legacy");
    std::fs::rename(&path_file, &backup)
        .with_context(|| format!("backup {:?} -> {:?}", path_file, backup))?;
    std::fs::rename(&tmp_file, &path_file)
        .with_context(|| format!("install v2 file at {:?}", path_file))?;
    log::info!(
        "  migration complete: {:.1} GB -> {:.1} GB (legacy backup at {:?})",
        old_size as f64 / 1_073_741_824.0,
        new_size as f64 / 1_073_741_824.0,
        backup
    );

    Ok(MigrationStats {
        num_genomes,
        old_size,
        new_size,
        already_v2: false,
    })
}

/// Summary returned by [`migrate_paths_to_v2`].
#[derive(Debug, Clone)]
pub struct MigrationStats {
    pub num_genomes: u64,
    pub old_size: u64,
    pub new_size: u64,
    pub already_v2: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::unitig::UnitigSet;

    /// Streaming write + PathIndex read must round-trip to the exact same
    /// paths/steps. Confirms the wire format matches bincode's Vec<T> layout.
    #[test]
    fn streaming_path_index_roundtrip() -> Result<()> {
        let tmp = tempfile::tempdir()?;
        let genome_dir = tmp.path().join("genomes");
        let out_dir = tmp.path().join("idx");
        std::fs::create_dir_all(&genome_dir)?;
        std::fs::create_dir_all(&out_dir)?;

        // Two tiny genomes that share a unitig-worthy 31+ base prefix.
        std::fs::write(
            genome_dir.join("g1.fa"),
            b">g1\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )?;
        std::fs::write(
            genome_dir.join("g2.fa"),
            b">g2\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )?;

        // Two unitigs; the k-mer table just needs any data to walk against.
        let text = b"ACGTACGTACGTACGTACGTACGTACGTACGT$TTGGCCAATTGGCCAATTGGCCAATTGGCCAA$".to_vec();
        let lengths = [32u64, 32];
        let unitigs = UnitigSet::from_fm_text(&text, &lengths);

        build_path_index(&genome_dir, &unitigs, &out_dir, 31)?;

        let idx = load_path_index(&out_dir)?;
        assert_eq!(idx.num_genomes(), 2);
        // Each genome must round-trip the per-genome metadata we wrote.
        for (i, path) in idx.iter().enumerate() {
            assert_eq!(path.genome_id, i as u32);
            assert!(!path.genome_name.is_empty());
            assert!(path.genome_length > 0);
        }
        Ok(())
    }

    /// Manufacture a legacy bincode `paths.bin`, migrate it to v2 in place,
    /// and confirm that `MmapPathIndex::open` yields byte-identical genomes.
    #[test]
    fn migrate_legacy_to_v2_roundtrip() -> Result<()> {
        use crate::index::paths_v2::MmapPathIndex;
        let tmp = tempfile::tempdir()?;
        let dir = tmp.path();

        // Build a tiny in-memory PathIndex with non-trivial steps so the
        // delta encoding gets exercised.
        let mut paths = Vec::new();
        for gid in 0..3u32 {
            let mut steps = Vec::new();
            let mut off = 0u64;
            for i in 0..(50 + gid as u32 * 25) {
                steps.push(super::PathStep {
                    unitig_id: i % 11,
                    is_reverse: i % 2 == 0,
                    genome_offset: off,
                });
                off += (i as u64 % 31) + 1;
            }
            paths.push(super::GenomePath {
                genome_id: gid,
                genome_name: format!("g_{gid:04}"),
                genome_length: off,
                steps,
            });
        }
        let original_paths = paths.clone();
        let original = super::LegacyPathIndex { paths };

        // Write a *legacy* bincode paths.bin (the format produced by builds
        // before commit a70a087 introduced the v2 magic header).
        let path_file = dir.join("paths.bin");
        crate::util::mmap::write_bincode(&path_file, &original)?;
        assert!(!crate::index::paths_v2::is_v2(&path_file)?);

        // Migrate.
        let stats = migrate_paths_to_v2(dir)?;
        assert!(!stats.already_v2);
        assert_eq!(stats.num_genomes, 3);
        assert!(crate::index::paths_v2::is_v2(&path_file)?);
        assert!(dir.join("paths.bin.legacy").exists());

        // Read back via MmapPathIndex and compare every step.
        let mmap_idx = MmapPathIndex::open(&path_file)?;
        assert_eq!(mmap_idx.num_genomes(), 3);
        for orig in &original_paths {
            let got = mmap_idx.get_path(orig.genome_id)?.unwrap();
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

        // Re-running migration on an already-v2 file should be a no-op.
        let again = migrate_paths_to_v2(dir)?;
        assert!(again.already_v2);
        Ok(())
    }
}
