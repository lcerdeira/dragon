//! Round-trip + end-to-end validation of the v3 (graph-edge) path format
//! against the real GGCAT-built single-genome index in `tests/data/sg_index`.
//!
//! Gate: a wrong v3 decode silently corrupts every alignment, so before any
//! shard is migrated on HPC these must pass:
//!   1. v3.get_path == v2.get_path for every genome (exact).
//!   2. a search over a v3-migrated index reproduces the v2 result.
//!
//! Fixture is git-ignored; tests skip gracefully when absent.
//!
//! Run: cargo test --test integration_paths_v3 --release -- --nocapture

use std::path::Path;

use dragon::index::paths::load_path_index;
use dragon::index::paths_v2::MmapPathIndex;
use dragon::index::paths_v3::{migrate_v2_to_v3, MmapPathIndexV3};
use dragon::query::{search, SearchConfig};

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
fn v3_migration_roundtrips_real_index() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let v2_path = data.join("sg_index/paths.bin");
    if !v2_path.exists() {
        eprintln!("SKIP v3_migration_roundtrips_real_index: fixture absent");
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let v3_path = tmp.path().join("paths_v3.bin");

    let stats = migrate_v2_to_v3(&v2_path, &v3_path).expect("migrate v2->v3");
    println!(
        "migrated {} genomes, {} unitigs, {} edges: {} B -> {} B ({:.2}x)",
        stats.num_genomes, stats.num_unitigs, stats.num_edges,
        stats.old_size, stats.new_size,
        stats.old_size as f64 / stats.new_size.max(1) as f64,
    );

    let v2 = MmapPathIndex::open(&v2_path).expect("open v2");
    let v3 = MmapPathIndexV3::open(&v3_path).expect("open v3");
    assert_eq!(v2.num_genomes(), v3.num_genomes());

    for gid in 0..v2.num_genomes() {
        let a = v2.get_path(gid as u32).unwrap().expect("v2 genome");
        let b = v3.get_path(gid as u32).unwrap().expect("v3 genome");
        assert_eq!(a.genome_name, b.genome_name, "genome {gid} name");
        assert_eq!(a.genome_length, b.genome_length, "genome {gid} length");
        assert_eq!(a.steps.len(), b.steps.len(), "genome {gid} step count");
        for (i, (x, y)) in a.steps.iter().zip(b.steps.iter()).enumerate() {
            assert_eq!(x.unitig_id, y.unitig_id, "genome {gid} step {i} unitig_id");
            assert_eq!(x.is_reverse, y.is_reverse, "genome {gid} step {i} is_reverse");
            assert_eq!(
                x.genome_offset, y.genome_offset,
                "genome {gid} step {i} genome_offset"
            );
        }
    }
    println!("v3 decodes byte-identically to v2 for all {} genomes", v2.num_genomes());
}

#[test]
fn search_works_on_v3_migrated_index() {
    let _ = env_logger::builder().is_test(true).try_init();
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let fasta_path = data.join("SAMEA110247553.fa");
    let idx_src = data.join("sg_index");
    if !fasta_path.exists() || !idx_src.join("paths.bin").exists() {
        eprintln!("SKIP search_works_on_v3_migrated_index: fixture absent");
        return;
    }

    // Copy the index to a tempdir and replace paths.bin with a v3 migration.
    let tmp = tempfile::tempdir().unwrap();
    let idx = tmp.path().join("idx");
    std::fs::create_dir_all(&idx).unwrap();
    for entry in std::fs::read_dir(&idx_src).unwrap() {
        let e = entry.unwrap();
        std::fs::copy(e.path(), idx.join(e.file_name())).unwrap();
    }
    let v2_tmp = idx.join("paths.bin.v2");
    std::fs::rename(idx.join("paths.bin"), &v2_tmp).unwrap();
    migrate_v2_to_v3(&v2_tmp, &idx.join("paths.bin")).expect("migrate");

    // load_path_index must dispatch to the v3 reader.
    assert!(
        matches!(
            load_path_index(&idx).unwrap(),
            dragon::index::paths::PathIndex::MmapV3(_)
        ),
        "load_path_index should select the v3 backing"
    );

    // The uniq3 self-slice must still align perfectly through the v3 path.
    let fasta = std::fs::read(&fasta_path).unwrap();
    let seq = concat_fasta(&fasta);
    let query = &seq[1999..2499];
    let q_path = tmp.path().join("q.fa");
    let mut q = b">unique3\n".to_vec();
    q.extend_from_slice(query);
    q.push(b'\n');
    std::fs::write(&q_path, &q).unwrap();

    let config = SearchConfig {
        index_dir: idx.into_boxed_path(),
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
    let results = search(&q_path, &config).expect("search");
    let best = &results[0].alignments[0];
    println!(
        "v3 search: t={}..{} matches={} id={:.3}",
        best.target_start, best.target_end, best.num_matches, best.identity()
    );
    assert_eq!(best.num_matches, 500, "every base must match through v3 path");
    assert!(best.identity() >= 0.99, "v3 identity {:.3}", best.identity());
    assert!(
        (best.target_start as i64 - 1999).abs() <= 40,
        "v3 target_start {} off",
        best.target_start
    );
}
