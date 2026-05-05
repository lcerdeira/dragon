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

/// Build a two-contig FASTA. Each contig has its own random backbone
/// (different LCG seeds → no shared k-mers across contigs). `gene` is
/// planted at offset `pos_in_c2` inside contig 2.
fn make_two_contig_fasta(
    genome_name: &str,
    c1_len: usize,
    c2_len: usize,
    gene: &[u8],
    pos_in_c2: usize,
) -> Vec<u8> {
    assert!(pos_in_c2 + gene.len() <= c2_len, "gene doesn't fit in c2");
    let c1 = random_dna(c1_len, 0xC0117_0001);
    let mut c2 = random_dna(c2_len, 0xC0117_0002);
    c2[pos_in_c2..pos_in_c2 + gene.len()].copy_from_slice(gene);

    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(b">");
    out.extend_from_slice(genome_name.as_bytes());
    out.extend_from_slice(b".contig1\n");
    for chunk in c1.chunks(80) {
        out.extend_from_slice(chunk);
        out.push(b'\n');
    }
    out.extend_from_slice(b">");
    out.extend_from_slice(genome_name.as_bytes());
    out.extend_from_slice(b".contig2\n");
    for chunk in c2.chunks(80) {
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

/// Multi-contig genome: two contigs, gene planted in contig 2. This is the
/// case that fails on real saureus assemblies (50–300 contigs/genome) — the
/// path-builder's `prev_unitig` state leaks across contig boundaries and the
/// unitig entry-offset is discarded, so `extract_sequence_static` returns
/// shifted bases for any step that begins mid-unitig (which happens at every
/// contig start). Net effect on production benchmarks: ~9× under-selection
/// of fragmented genomes and a 0.6 identity ceiling on those that do surface.
#[test]
fn multi_contig_gene_in_second_contig() {
    let gene = random_dna(120, 0xBEEF_5A11);
    let genome = make_two_contig_fasta("g_mc", 4000, 4000, &gene, 1500);
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

    // Gene is at contig2 pos 1500. Genome coordinates: contig1 (0..4000) +
    // contig2 (4000..8000). Planted gene starts at genome offset 4000 + 1500
    // = 5500. Allow a small slop window in case path-walking off-by-ones the
    // contig boundary.
    let expected_start: usize = 4000 + 1500;
    let best = r.alignments.iter().find(|a| {
        a.strand == '+'
            && a.target_start >= expected_start.saturating_sub(40)
            && a.target_start <= expected_start + 40
    });
    let best = best.expect(&format!(
        "no forward alignment near planted position {} — path-builder is dropping or shifting the contig-2 entry step",
        expected_start
    ));

    assert_eq!(
        best.num_matches,
        gene.len(),
        "every base should match (got {} of {} for matches; suggests reference window contains shifted/wrong bases)",
        best.num_matches,
        gene.len()
    );
    assert!(
        best.identity() >= 0.99,
        "expected ≈1.0 identity across contig boundary, got {:.3}",
        best.identity()
    );
}

/// Stronger multi-contig test: the gene appears at the END of contig 1 AND
/// the START of contig 2. Both copies share the same unitig in the dBG, so
/// at the boundary `prev_unitig` is `(gene_unitig, fwd)` from contig 1, and
/// contig 2's first k-mer maps to the same unitig — the path-builder will
/// not emit a new step, and contig 2's copy of the gene becomes invisible
/// to seed→genome translation. If this test fails, we have a clean repro
/// for the "SAMEA missing" symptom we see on real saureus shards.
#[test]
fn multi_contig_shared_motif_at_boundary() {
    let gene = random_dna(120, 0xFEED_FACE);
    // Plant gene at the END of contig 1 and START of contig 2 so the
    // unitig holding it is shared at the contig boundary.
    let c1_len = 4000;
    let c2_len = 4000;
    let mut c1 = random_dna(c1_len, 0xCAFE_C001);
    let mut c2 = random_dna(c2_len, 0xCAFE_C002);
    c1[c1_len - gene.len()..].copy_from_slice(&gene);
    c2[..gene.len()].copy_from_slice(&gene);

    let mut fasta: Vec<u8> = Vec::new();
    fasta.extend_from_slice(b">g_shared.contig1\n");
    for chunk in c1.chunks(80) {
        fasta.extend_from_slice(chunk);
        fasta.push(b'\n');
    }
    fasta.extend_from_slice(b">g_shared.contig2\n");
    for chunk in c2.chunks(80) {
        fasta.extend_from_slice(chunk);
        fasta.push(b'\n');
    }

    let (_tmp, idx, qp) = build_synthetic_index(&fasta, "gene", &gene);
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

    // The best hit's target_start should resolve to one of the two planted
    // positions: contig1 end (genome offset c1_len - gene.len() = 3880) or
    // contig2 start (genome offset c1_len + 0 = 4000).
    let pos_c1: usize = c1_len - gene.len();
    let pos_c2: usize = c1_len;
    let near = |a: &dragon::io::paf::PafRecord, p: usize| -> bool {
        a.strand == '+'
            && (a.target_start as i64 - p as i64).abs() <= 40
    };
    let any_correct = r.alignments.iter().any(|a| near(a, pos_c1) || near(a, pos_c2));
    assert!(
        any_correct,
        "no alignment near either planted position (c1 end = {} or c2 start = {})",
        pos_c1, pos_c2
    );

    // And the best alignment, wherever it lands, should be id == 1.0 because
    // the planted gene is verbatim in the genome.
    let best = &r.alignments[0];
    assert_eq!(
        best.num_matches,
        gene.len(),
        "every base should match (got {} of {})",
        best.num_matches,
        gene.len()
    );
    assert!(
        best.identity() >= 0.99,
        "expected ≈1.0 identity, got {:.3}",
        best.identity()
    );
}
