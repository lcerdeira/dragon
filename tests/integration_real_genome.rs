//! Real-data regression test: reproduces the saureus `mean_id 0.624` bug
//! with a single real genome and a GGCAT-built production index.
//!
//! `SAMEA110247553` is a 49-contig S. aureus draft assembly (2.8 Mbp). The
//! query is a verbatim 500 bp slice of its own sequence at concatenated
//! genome offset 1999..2499 (awk `substr(seq,2000,500)`, 1-indexed).
//!
//! Ground truth is exact: searching this query against an index built from
//! this one genome MUST return SAMEA110247553, '+' strand, target_start ≈
//! 1999, num_matches == 500, identity == 1.0.
//!
//! On HPC this produced target_start=112, num_matches=277, identity≈0.55 —
//! the reference window was extracted from the wrong genome coordinates.
//!
//! The index in `tests/data/sg_index` is a real GGCAT-built index (the
//! production build path), copied from HPC, so this exercises the same
//! paths.bin format the saureus shards use.
//!
//! Run: cargo test --test integration_real_genome --release -- --nocapture

use std::path::Path;

use dragon::query::{search, SearchConfig};

const QUERY_OFFSET: usize = 1999; // 0-based; awk substr(seq,2000,500)
const QUERY_LEN: usize = 500;

/// Concatenate all non-header lines of a FASTA into one sequence (matches
/// the awk `!/^>/{seq=seq $0}` that generated the original query).
fn concat_fasta(fasta_bytes: &[u8]) -> Vec<u8> {
    let mut seq = Vec::new();
    for line in fasta_bytes.split(|&b| b == b'\n') {
        if line.first() == Some(&b'>') || line.is_empty() {
            continue;
        }
        seq.extend_from_slice(line);
    }
    seq
}

#[test]
fn real_genome_self_slice_aligns_perfectly() {
    let _ = env_logger::builder().is_test(true).try_init();

    let data_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let fasta_path = data_dir.join("SAMEA110247553.fa");
    let idx_dir = data_dir.join("sg_index");

    // The fixture (a real 2.8 Mbp genome + its GGCAT-built index, ~30 MB)
    // is intentionally not committed — it lives in `tests/data/`, which is
    // git-ignored. Skip gracefully when it is absent (e.g. in CI) rather
    // than failing; the separator-offset regression is also covered by the
    // self-contained unit test `elias_fano::tests::test_many_unitigs_*`.
    if !fasta_path.exists() || !idx_dir.join("paths.bin").exists() {
        eprintln!("SKIP real_genome_self_slice_aligns_perfectly: \
                   fixture tests/data/ not present");
        return;
    }
    let fasta = std::fs::read(&fasta_path).expect("read SAMEA110247553.fa");

    let seq = concat_fasta(&fasta);
    assert!(seq.len() >= QUERY_OFFSET + QUERY_LEN, "genome too short");
    let query: Vec<u8> = seq[QUERY_OFFSET..QUERY_OFFSET + QUERY_LEN].to_vec();

    let tmp = tempfile::tempdir().expect("tempdir");
    let q_path = tmp.path().join("query.fa");
    let mut q_bytes: Vec<u8> = b">unique3\n".to_vec();
    q_bytes.extend_from_slice(&query);
    q_bytes.push(b'\n');
    std::fs::write(&q_path, &q_bytes).unwrap();

    let config = SearchConfig {
        index_dir: idx_dir.into_boxed_path(),
        min_seed_len: 15,
        max_seed_freq: 10_000,
        min_chain_score: 0.0,
        max_target_seqs: 10,
        threads: 1,
        max_ram_gb: 4.0,
        min_identity: 0.0,
        min_query_coverage: 0.0,
        min_score_ratio: 0.0,
        no_ml: true,
        ml_weights_path: None,
        dump_seeds_path: None,
        ground_truth_genome: None,
        batch_queries: true,
        parallel_shards: false,
        cross_species: false,
    };
    let results = search(&q_path, &config).expect("search");
    assert_eq!(results.len(), 1);
    let r = &results[0];

    println!("{} alignments for '{}':", r.alignments.len(), r.query_name);
    for a in &r.alignments {
        println!(
            "  q={}..{} strand={} t={}..{} matches={} alen={} id={:.3}",
            a.query_start, a.query_end, a.strand,
            a.target_start, a.target_end, a.num_matches,
            a.alignment_len, a.identity()
        );
    }
    assert!(!r.alignments.is_empty(), "expected at least one alignment");

    let best = &r.alignments[0];
    assert_eq!(best.strand, '+', "expected forward-strand hit");
    assert!(
        (best.target_start as i64 - QUERY_OFFSET as i64).abs() <= 40,
        "target_start={} but query is at genome offset {} — reference window \
         extracted from wrong coordinates",
        best.target_start, QUERY_OFFSET
    );
    assert_eq!(
        best.num_matches, QUERY_LEN,
        "every base of a verbatim self-slice must match (got {}/{})",
        best.num_matches, QUERY_LEN
    );
    assert!(
        best.identity() >= 0.99,
        "expected identity ≈1.0 for a verbatim self-slice, got {:.3}",
        best.identity()
    );
}
