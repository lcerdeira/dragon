/// Genome path index: stores each genome as a sequence of unitig IDs
/// traversed when walking the genome through the de Bruijn graph.

use anyhow::Result;
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
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PathIndex {
    pub paths: Vec<GenomePath>,
}

impl PathIndex {
    pub fn get_path(&self, genome_id: u32) -> Option<&GenomePath> {
        self.paths.get(genome_id as usize)
    }

    pub fn num_genomes(&self) -> usize {
        self.paths.len()
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
/// For new code that wants O(1) load time and per-genome lazy decoding, use
/// [`open_path_index_mmap`] instead — this function eagerly materializes
/// every genome into RAM for compatibility with legacy callers.
pub fn load_path_index(index_dir: &Path) -> Result<PathIndex> {
    let path_file = index_dir.join("paths.bin");
    if crate::index::paths_v2::is_v2(&path_file)? {
        // Eagerly decode all genomes — preserves the legacy in-memory shape.
        // Hot-path consumers should switch to open_path_index_mmap().
        let mmap_idx = crate::index::paths_v2::MmapPathIndex::open(&path_file)?;
        let n = mmap_idx.num_genomes();
        let mut paths = Vec::with_capacity(n as usize);
        for gid in 0..n {
            let p = mmap_idx
                .get_path(gid as u32)?
                .ok_or_else(|| anyhow::anyhow!("missing genome {gid}"))?;
            paths.push(p);
        }
        Ok(PathIndex { paths })
    } else {
        // Legacy bincode `Vec<GenomePath>` format.
        crate::util::mmap::read_bincode(&path_file)
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
        for (i, path) in idx.paths.iter().enumerate() {
            assert_eq!(path.genome_id, i as u32);
            assert!(!path.genome_name.is_empty());
            assert!(path.genome_length > 0);
        }
        Ok(())
    }
}
