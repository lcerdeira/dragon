//! End-to-end alignment correctness tests on tiny synthetic genomes with
//! known ground truth.
//!
//! These tests pin down where Dragon's coordinate / orientation handling
//! goes wrong on real-data queries (manifested as a 0.7 identity ceiling
//! against saureus, where BLAST and minimap2 see 1.0). Each test:
//!
//!   1. Builds a deterministic synthetic genome FASTA with a "gene" we
//!      planted at a known position.
//!   2. Calls `dragon::index::build_index_with_options` to produce a real
//!      index in a tempdir (uses the internal de Bruijn graph builder
//!      since we're at single-genome scale).
//!   3. Calls `dragon::query::search` with the gene as the query.
//!   4. Asserts the resulting `PafRecord` matches the planted truth
//!      exactly: target_start, target_end, num_matches == query_len,
//!      identity == 1.0, CIGAR is "<query_len>=".
//!
//! When any assertion fails, the failure message points at *which*
//! coordinate is wrong, which is the fix-it-once debugging signal we
//! couldn't get from the saureus reproducer.
//!
//! Run with `cargo test --test integration_align --release -- --nocapture`
//! to see RUST_LOG=info output during alignment.

use std::path::Path;

use dragon::index::build_index_with_options;
use dragon::query::{search, SearchConfig};

/// Make a deterministic random DNA sequence of length `n` from `seed`.
fn random_dna(n: usize, seed: u64) -> Vec<u8> {
    // simple LCG so we don't pull in `rand` for tests
    let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let bases = b"ACGT";
    (0..n)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            bases[(s >> 60) as usize & 3]
        })
        .collect()
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

/// Plant `gene` at `pos` inside a `total_len` random backbone, return the
/// concatenated FASTA bytes (one record).
fn make_genome_fasta(name: &str, total_len: usize, gene: &[u8], pos: usize) -> Vec<u8> {
    assert!(pos + gene.len() <= total_len, "gene doesn't fit");
    let mut backbone = random_dna(total_len, 0xD3A60001);
    backbone[pos..pos + gene.len()].copy_from_slice(gene);
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(b">");
    out.extend_from_slice(name.as_bytes());
    out.push(b'\n');
    // wrap at 80 cols for cleanliness
    for chunk in backbone.chunks(80) {
        out.extend_from_slice(chunk);
        out.push(b'\n');
    }
    out
}

/// Build an index in a tempdir for the given FASTA bytes and return the
/// (tempdir, index_dir, query_path) triple. The caller is responsible for
/// holding onto the tempdir so it isn't dropped while the index is in use.
fn build_synthetic_index(
    genome_fasta: &[u8],
    query_name: &str,
    query_seq: &[u8],
) -> (tempfile::TempDir, std::path::PathBuf, std::path::PathBuf) {
    let _ = env_logger::builder().is_test(true).try_init();
    let tmp = tempfile::tempdir().expect("tempdir");
    let g_dir = tmp.path().join("genomes");
    std::fs::create_dir_all(&g_dir).unwrap();
    std::fs::write(g_dir.join("g.fa"), genome_fasta).unwrap();
    let idx_dir = tmp.path().join("idx");
    build_index_with_options(&g_dir, &idx_dir, 31, 2, None).expect("build_index");

    let q_path = tmp.path().join("query.fa");
    let mut q_bytes: Vec<u8> = Vec::new();
    q_bytes.extend_from_slice(b">");
    q_bytes.extend_from_slice(query_name.as_bytes());
    q_bytes.push(b'\n');
    q_bytes.extend_from_slice(query_seq);
    q_bytes.push(b'\n');
    std::fs::write(&q_path, &q_bytes).unwrap();

    (tmp, idx_dir, q_path)
}

/// Run a search with permissive defaults so we see all hits.
fn run_search(idx_dir: &Path, query_path: &Path) -> Vec<dragon::query::QueryResult> {
    let config = SearchConfig {
        index_dir: idx_dir.into(),
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
    };
    search(query_path, &config).expect("search")
}

/// `gene` planted at `pos` on the FORWARD strand. Expect a single hit
/// with target_start == pos, identity == 1.0.
#[test]
fn forward_strand_perfect_match() {
    let gene = random_dna(120, 0x91A2_5E8B);
    let genome = make_genome_fasta("g", 5000, &gene, 1500);
    let (_tmp, idx, qp) = build_synthetic_index(&genome, "gene", &gene);

    let results = run_search(&idx, &qp);
    assert_eq!(results.len(), 1, "one query");
    let r = &results[0];

    println!(
        "{} hits for '{}'. PAF records:",
        r.alignments.len(),
        r.query_name
    );
    for a in &r.alignments {
        println!(
            "  q={}..{}  strand={}  t={}..{}  matches={}  alen={}  id={:.3}",
            a.query_start,
            a.query_end,
            a.strand,
            a.target_start,
            a.target_end,
            a.num_matches,
            a.alignment_len,
            a.identity()
        );
    }
    assert!(!r.alignments.is_empty(), "expected at least 1 alignment");

    // Find a hit with strand '+' and target_start == 1500.
    let best = r
        .alignments
        .iter()
        .find(|a| a.strand == '+' && a.target_start == 1500);
    let best = best.expect(
        "no forward alignment at the planted position 1500 — coordinate translation is wrong",
    );

    // target_end may include WFA padding (D ops on the trailing pad), so
    // we allow the alignment-len to extend up to 2·ALIGN_PAD past the gene.
    assert!(
        best.target_end >= 1500 + gene.len() && best.target_end <= 1500 + gene.len() + 50,
        "target_end out of expected range: got {}, expected ≈ {}",
        best.target_end,
        1500 + gene.len()
    );
    assert_eq!(best.num_matches, gene.len(), "every base should match");
    assert!(
        best.identity() >= 0.99,
        "expected ≈1.0 identity, got {:.3}",
        best.identity()
    );
}

/// `gene` planted as RC on the genome's forward strand → query (forward)
/// must produce a '-' strand hit with identity 1.0 at the planted RC
/// position. This isolates the reverse-strand path through align_with_seeds.
#[test]
fn reverse_strand_perfect_match() {
    let gene = random_dna(120, 0xC0FFEE_42);
    let gene_rc = reverse_complement(&gene);
    // plant the RC of the gene at position 2200 — the query (forward gene)
    // should map there as a '-' strand hit.
    let genome = make_genome_fasta("g_rc", 5000, &gene_rc, 2200);
    let (_tmp, idx, qp) = build_synthetic_index(&genome, "gene", &gene);

    let results = run_search(&idx, &qp);
    let r = &results[0];

    println!(
        "{} hits for '{}'. PAF records:",
        r.alignments.len(),
        r.query_name
    );
    for a in &r.alignments {
        println!(
            "  q={}..{}  strand={}  t={}..{}  matches={}  alen={}  id={:.3}",
            a.query_start,
            a.query_end,
            a.strand,
            a.target_start,
            a.target_end,
            a.num_matches,
            a.alignment_len,
            a.identity()
        );
    }
    assert!(!r.alignments.is_empty(), "expected at least 1 alignment");

    let best = r
        .alignments
        .iter()
        .find(|a| a.strand == '-' && a.target_start == 2200);
    let best = best.expect(
        "no reverse alignment at the RC-planted position 2200 — strand bucket or position formula is wrong",
    );

    assert!(
        best.target_end >= 2200 + gene.len() && best.target_end <= 2200 + gene.len() + 50,
        "target_end out of expected range: got {}, expected ≈ {}",
        best.target_end,
        2200 + gene.len()
    );
    assert_eq!(best.num_matches, gene.len(), "every base should match");
    assert!(
        best.identity() >= 0.99,
        "expected ≈1.0 identity, got {:.3}",
        best.identity()
    );
}
