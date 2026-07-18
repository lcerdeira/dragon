//! Zarr-backed full alignment → BLAST output (path B).
//!
//! Reuses the binary-index aligner (`direct_align_candidates`, with its
//! seed-cap / contig-split / candidate fixes) but sources the reference from a
//! Zarr store. The store's `text`, `suffix_array`, `colors` and `paths` arrays
//! are loaded into the in-memory index types once, so seeding and alignment run
//! at binary-index speed (no per-k-mer Zarr chunk reads). Requires `/paths` in
//! the store (export after `migrate-paths`). Local store only for now; the
//! HTTP-lazy variant is future work (path A).

use anyhow::Result;
use std::collections::HashMap;
use std::path::Path;

use roaring::RoaringBitmap;

use crate::ds::elias_fano::CumulativeLengthIndex;
use crate::index::fm::{DragonFmIndex, SeedHit};
use crate::index::paths::{GenomePath, PathIndex};
use crate::index::unitig::UnitigSet;
use crate::index::zarr_backend::{ZarrColorIndex, ZarrFmIndex, ZarrPathIndex};
use crate::io::paf::PafRecord;
use crate::query::containment::ContainmentHit;
use crate::query::direct_align::direct_align_candidates;

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            o => o,
        })
        .collect()
}

/// In-memory reference structures materialised once from a Zarr store.
pub struct ZarrAlignRef {
    pub fm: DragonFmIndex,
    /// Colour bitmap per unitig (indexed by unitig id).
    pub colors: Vec<RoaringBitmap>,
    pub path_index: PathIndex,
    pub unitigs: UnitigSet,
    pub k: usize,
    /// Effective database size (text bytes) for BLAST E-values.
    pub db_size: u64,
}

/// Load a filesystem Zarr store into in-memory alignment structures.
/// Requires `/paths` in the store.
pub fn load_zarr_ref(zarr: &Path) -> Result<ZarrAlignRef> {
    let zfm = ZarrFmIndex::open(zarr)?;
    let zcolors = ZarrColorIndex::open(zarr)?;
    let zpaths = ZarrPathIndex::open(zarr)?;

    let k = zfm.kmer_size;
    let db_size = zfm.text_len;
    let lengths: Vec<u64> = zfm.cumulative_lengths.lengths().to_vec();
    let text = zfm.text_slice(0, zfm.text_len)?;
    let unitigs = UnitigSet::from_fm_text(&text, &lengths);
    let suffix_array = zfm.read_all_sa()?;

    let fm = DragonFmIndex {
        text,
        suffix_array,
        eytzinger_sa: Vec::new(),
        cumulative_lengths: CumulativeLengthIndex::from_lengths(&lengths),
    };

    // Materialise colours per unitig (one Roaring bitmap each).
    let nu = zcolors.num_unitigs() as u32;
    let mut colors: Vec<RoaringBitmap> = Vec::with_capacity(nu as usize);
    for u in 0..nu {
        colors.push(zcolors.get_colors(u)?);
    }

    // Materialise genome paths (index i == genome i).
    let mut gpaths: Vec<GenomePath> = Vec::with_capacity(zpaths.num_genomes as usize);
    for g in 0..zpaths.num_genomes as u32 {
        gpaths.push(zpaths.get_path(g)?.unwrap_or(GenomePath {
            genome_id: g,
            genome_name: String::new(),
            genome_length: 0,
            steps: Vec::new(),
        }));
    }
    let path_index = PathIndex::from_paths(gpaths);

    Ok(ZarrAlignRef {
        fm,
        colors,
        path_index,
        unitigs,
        k,
        db_size,
    })
}

/// Align one query against a Zarr-backed reference, returning PAF records.
/// `max_target_seqs == 0` means no cap.
/// Align one query against the Zarr-sourced reference.
///
/// `align_once` enables the align-once aligner: candidates sharing a
/// byte-identical reference window are aligned with a single WFA call and the
/// result is projected onto every genome that shares it. Output is identical to
/// the per-genome path; only the WFA count changes. This matters most for the
/// published cloud catalog, where a single gene query can match tens of
/// thousands of genomes.
#[allow(clippy::too_many_arguments)]
pub fn align_query(
    r: &ZarrAlignRef,
    query: &[u8],
    query_name: &str,
    max_target_seqs: usize,
    max_seed_freq: usize,
    min_identity: f64,
    min_query_coverage: f64,
    align_once: bool,
) -> Result<Vec<PafRecord>> {
    let k = r.k;
    if query.len() < k {
        return Ok(Vec::new());
    }
    let total_kmers = query.len() - k + 1;
    // Spread ~500 sampled seed positions across the query so seeds span the
    // whole gene (the split-contig fix needs this) while bounding work.
    let stride = (total_kmers / 500).max(1);
    // ~1 seed per 16 bp per genome (mirrors the binary-index seed-cap fix).
    let max_seeds_per_genome = (query.len() / 16).max(64);
    const CANDIDATE_TARGET: usize = 8000;

    let mut genome_seeds: HashMap<u32, Vec<SeedHit>> = HashMap::new();
    let mut qpos = 0usize;
    while qpos < total_kmers {
        let kmer = &query[qpos..qpos + k];
        for (is_rev, pat) in [(false, kmer.to_vec()), (true, revcomp(kmer))] {
            let positions = r.fm.search(&pat);
            if positions.is_empty() || positions.len() > max_seed_freq {
                continue;
            }
            for pos in positions {
                if let Some((uid, off)) = r.fm.position_to_unitig(pos) {
                    if let Some(bm) = r.colors.get(uid as usize) {
                        for g in bm.iter() {
                            let bucket = genome_seeds.entry(g).or_default();
                            if bucket.len() < max_seeds_per_genome {
                                bucket.push(SeedHit {
                                    unitig_id: uid,
                                    offset: off,
                                    query_pos: qpos,
                                    match_len: k,
                                    is_reverse: is_rev,
                                    sa_count: 0,
                                });
                            }
                        }
                    }
                }
            }
        }
        qpos += stride;
    }
    if genome_seeds.is_empty() {
        return Ok(Vec::new());
    }

    let mut hits: Vec<ContainmentHit> = genome_seeds
        .into_iter()
        .map(|(g, seeds)| {
            // Distinct matched query positions, not seed instances: a k-mer can
            // seed several unitig positions and both strands, so seeds.len() would
            // overcount and push containment above 1.0.
            let shared = seeds
                .iter()
                .map(|s| s.query_pos)
                .collect::<std::collections::HashSet<_>>()
                .len();
            ContainmentHit {
                genome_id: g,
                containment: shared as f64 / total_kmers.max(1) as f64,
                shared_kmers: shared,
                total_query_kmers: total_kmers,
                info_score: shared as f64,
                bayes_prob: None,
                bayes_prob_hc: None,
                bayes_ani: None,
                seeds,
            }
        })
        .collect();
    // Keep the strongest candidates (bounds alignment work).
    hits.sort_by(|a, b| b.shared_kmers.cmp(&a.shared_kmers));
    hits.truncate(CANDIDATE_TARGET);

    let cap = if max_target_seqs == 0 {
        usize::MAX
    } else {
        max_target_seqs
    };
    let mut records =
        direct_align_candidates(query, query_name, &hits, &r.path_index, &r.unitigs, cap, align_once);
    // Same post-alignment filters as the binary `dragon search` pipeline.
    let ql = query.len() as f64;
    records.retain(|rec| {
        rec.identity() >= min_identity
            && ql > 0.0
            && (rec.query_end.saturating_sub(rec.query_start)) as f64 / ql >= min_query_coverage
    });
    Ok(records)
}
