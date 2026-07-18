//! Correctness of the Zarr-backed aligner's candidate scoring.
//!
//! `zarr_align::align_query` re-implements the seeding/ranking that
//! `containment.rs` does for binary indices, and the two drifted apart. These
//! tests pin the two properties that drift broke:
//!
//!   1. **Containment denominator.** Seeds are sampled with a stride once a
//!      query exceeds ~500 k-mers, but containment was divided by *every*
//!      k-mer in the query, so a perfect match could never report 1.0 — a
//!      2 kb self-match reported ≈0.33 and a 6 kb one ≈0.09.
//!   2. **IDF-weighted ranking.** Candidates were ranked by raw shared-k-mer
//!      count, so genomes sharing only conserved/core sequence could outrank
//!      the true source and crowd it out of the candidate cap.

use std::path::Path;

use dragon::index::build_index_with_options;
use dragon::index::zarr_backend::export_to_zarr;
use dragon::query::zarr_align::{align_query, load_zarr_ref};

/// Deterministic pseudo-random DNA (no `rand` dependency in tests).
fn random_dna(n: usize, seed: u64) -> Vec<u8> {
    let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let bases = b"ACGT";
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[(s >> 60) as usize & 3]
        })
        .collect()
}

/// Build a Zarr store from one FASTA file per genome; returns (tmp, zarr_path).
fn build_zarr(genomes: &[(&str, Vec<u8>)]) -> (tempfile::TempDir, std::path::PathBuf) {
    let _ = env_logger::builder().is_test(true).try_init();
    let tmp = tempfile::tempdir().expect("tempdir");
    let g_dir = tmp.path().join("genomes");
    std::fs::create_dir_all(&g_dir).unwrap();
    for (name, seq) in genomes {
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(b">");
        out.extend_from_slice(name.as_bytes());
        out.push(b'\n');
        for chunk in seq.chunks(80) {
            out.extend_from_slice(chunk);
            out.push(b'\n');
        }
        std::fs::write(g_dir.join(format!("{name}.fa")), &out).unwrap();
    }
    let idx = tmp.path().join("idx");
    build_index_with_options(&g_dir, &idx, 31, 2, None).expect("build_index");
    let zarr = tmp.path().join("store.zarr");
    export_to_zarr(&idx, &zarr).expect("export_to_zarr");
    (tmp, zarr)
}

fn ct_tag(rec: &dragon::io::paf::PafRecord) -> f64 {
    rec.tags
        .iter()
        .find_map(|t| t.strip_prefix("ct:f:").and_then(|v| v.parse::<f64>().ok()))
        .expect("ct:f: tag present")
}

fn run(zarr: &Path, name: &str, query: &[u8]) -> Vec<dragon::io::paf::PafRecord> {
    let aref = load_zarr_ref(zarr).expect("load_zarr_ref");
    align_query(&aref, query, name, 0, 10_000, 0.0, 0.0, false).expect("align_query")
}

/// A perfect self-match must report containment ≈ 1.0 regardless of query
/// length. Before the fix, any query past ~530 bp was sampled with stride > 1
/// yet divided by the full k-mer count, capping containment at ~1/stride.
#[test]
fn containment_is_one_for_perfect_match_of_any_length() {
    let genome = random_dna(20_000, 0xC0FFEE01);
    let (_tmp, zarr) = build_zarr(&[("genomeA", genome.clone())]);

    // 400 bp → stride 1 (was already fine); 3 kb and 8 kb → stride > 1.
    for len in [400usize, 3_000, 8_000] {
        let q = &genome[2_000..2_000 + len];
        let recs = run(&zarr, &format!("self{len}"), q);
        assert!(!recs.is_empty(), "no hit for a {len} bp self-query");
        let ct = ct_tag(&recs[0]);
        assert!(
            ct > 0.95,
            "{len} bp perfect self-match reported containment {ct:.3}; \
             expected ≈1.0 (stride-sampled numerator must use the sampled \
             k-mer count as denominator)"
        );
    }
}

/// The true source must rank first even when many other genomes share a large
/// conserved block with the query. Ranking by raw shared-k-mer count let those
/// core-sharing genomes tie or win; IDF weighting drives conserved sequence to
/// ≈0 information so the source's unique sequence dominates.
#[test]
fn true_source_outranks_core_sharing_genomes() {
    // A conserved block every genome carries, plus a segment unique to the source.
    let core = random_dna(4_000, 0xABCD_0001);
    let unique = random_dna(1_200, 0xABCD_0002);

    let mut genomes: Vec<(String, Vec<u8>)> = Vec::new();
    // The true source: core + its unique segment.
    let mut src = core.clone();
    src.extend_from_slice(&unique);
    genomes.push(("source".to_string(), src));
    // Decoys: the same core plus filler that is unrelated to the query.
    for i in 0..8 {
        let mut d = core.clone();
        d.extend_from_slice(&random_dna(1_200, 0xDEAD_0000 + i as u64));
        genomes.push((format!("decoy{i}"), d));
    }
    let refs: Vec<(&str, Vec<u8>)> =
        genomes.iter().map(|(n, s)| (n.as_str(), s.clone())).collect();
    let (_tmp, zarr) = build_zarr(&refs);

    // Query = conserved block + the source-unique segment. Every decoy matches
    // the conserved 4 kb; only the source matches the final 1.2 kb.
    let mut query = core.clone();
    query.extend_from_slice(&unique);

    let recs = run(&zarr, "q", &query);
    assert!(!recs.is_empty(), "no hits at all");
    assert_eq!(
        recs[0].target_name, "source",
        "true source must rank first; got {} (ranking must use IDF-weighted \
         information, not raw shared-k-mer count). Order: {:?}",
        recs[0].target_name,
        recs.iter().map(|r| &r.target_name).collect::<Vec<_>>()
    );
}
