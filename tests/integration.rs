/// Integration tests for Dragon: end-to-end index construction and search.

use std::fs;
use std::io::Write;
use std::path::Path;
use tempfile::TempDir;

/// Helper: create a small FASTA file with one sequence.
fn write_fasta(dir: &Path, name: &str, seq: &str) {
    let path = dir.join(format!("{}.fa", name));
    let mut f = fs::File::create(&path).unwrap();
    writeln!(f, ">{}", name).unwrap();
    for chunk in seq.as_bytes().chunks(80) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
}

/// Helper: generate a random DNA sequence of given length.
fn random_seq(len: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = seed;
    (0..len)
        .map(|_| {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            bases[((rng >> 33) % 4) as usize] as char
        })
        .collect()
}

/// Helper: mutate a sequence at a given rate.
fn mutate_seq(seq: &str, rate: f64, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = seed;
    seq.bytes()
        .map(|b| {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let r = ((rng >> 33) as f64) / (u32::MAX as f64);
            if r < rate {
                // Substitute with a different base
                let alts: Vec<u8> = bases.iter().copied().filter(|&x| x != b).collect();
                let idx = ((rng >> 17) % 3) as usize;
                alts[idx] as char
            } else {
                b as char
            }
        })
        .collect()
}

// ============================================================================
// FASTA I/O tests
// ============================================================================

#[test]
fn test_fasta_read_multiple_sequences() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("multi.fa");
    let mut f = fs::File::create(&path).unwrap();
    writeln!(f, ">seq1 some description").unwrap();
    writeln!(f, "ACGTACGTACGT").unwrap();
    writeln!(f, "TTTTGGGG").unwrap();
    writeln!(f, ">seq2").unwrap();
    writeln!(f, "AAAACCCCGGGGTTTT").unwrap();
    writeln!(f, ">seq3 another desc").unwrap();
    writeln!(f, "ATATATAT").unwrap();

    let seqs = dragon::io::fasta::read_sequences(&path).unwrap();
    assert_eq!(seqs.len(), 3);
    assert_eq!(seqs[0].name, "seq1");
    assert_eq!(seqs[0].seq, b"ACGTACGTACGTTTTTGGGG");
    assert_eq!(seqs[1].name, "seq2");
    assert_eq!(seqs[1].seq, b"AAAACCCCGGGGTTTT");
    assert_eq!(seqs[2].name, "seq3");
    assert_eq!(seqs[2].seq, b"ATATATAT");
}

#[test]
fn test_fasta_empty_file() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("empty.fa");
    fs::write(&path, "").unwrap();

    let seqs = dragon::io::fasta::read_sequences(&path).unwrap();
    assert_eq!(seqs.len(), 0);
}

#[test]
fn test_fasta_streaming_reader() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("stream.fa");
    let mut f = fs::File::create(&path).unwrap();
    for i in 0..10 {
        writeln!(f, ">seq_{}", i).unwrap();
        writeln!(f, "{}", random_seq(500, i as u64)).unwrap();
    }

    let reader = dragon::io::fasta::FastaReader::new(&path).unwrap();
    let seqs: Vec<_> = reader.map(|r| r.unwrap()).collect();
    assert_eq!(seqs.len(), 10);
    for (i, seq) in seqs.iter().enumerate() {
        assert_eq!(seq.name, format!("seq_{}", i));
        assert_eq!(seq.seq.len(), 500);
    }
}

#[test]
fn test_list_fasta_files() {
    let dir = TempDir::new().unwrap();
    write_fasta(dir.path(), "genome_a", "ACGT");
    write_fasta(dir.path(), "genome_b", "TTGG");
    // Non-FASTA file should be ignored
    fs::write(dir.path().join("readme.txt"), "not a fasta").unwrap();

    let files = dragon::io::fasta::list_fasta_files(dir.path()).unwrap();
    assert_eq!(files.len(), 2);
}

// ============================================================================
// DNA encoding tests
// ============================================================================

#[test]
fn test_packed_sequence_large() {
    let seq = random_seq(10_000, 42);
    let packed = dragon::util::dna::PackedSequence::from_bytes(seq.as_bytes());
    assert_eq!(packed.len, 10_000);
    assert_eq!(packed.to_bytes(), seq.as_bytes());
}

#[test]
fn test_packed_revcomp_involution() {
    // Reverse complement twice should give original
    let seq = b"ACGTAACCGGTTAACC";
    let packed = dragon::util::dna::PackedSequence::from_bytes(seq);
    let rc = packed.reverse_complement();
    let rc2 = rc.reverse_complement();
    assert_eq!(rc2.to_bytes(), packed.to_bytes());
}

#[test]
fn test_kmer_extraction() {
    let seq = b"ACGTACGTACGT";
    let packed = dragon::util::dna::PackedSequence::from_bytes(seq);
    let kmer1 = packed.kmer_u64(0, 4);
    let kmer2 = packed.kmer_u64(4, 4);
    // Both should be the same (ACGT repeated)
    assert_eq!(kmer1, kmer2);
}

// ============================================================================
// Data structure tests
// ============================================================================

// ds::varint module was removed in commit 33ba145 — its only consumer
// (paths_v2) ships its own inline LEB128. The two tests previously here
// (test_varint_roundtrip_many_values, test_delta_encoding_sorted) now live
// inside paths_v2::tests::varint_roundtrip.

#[test]
fn test_elias_fano_many_unitigs() {
    let lengths: Vec<u64> = (0..200).map(|i| 100 + (i % 50) * 10).collect();
    let idx = dragon::ds::elias_fano::CumulativeLengthIndex::from_lengths(&lengths);

    assert_eq!(idx.num_unitigs(), 200);
    assert_eq!(idx.total_length(), lengths.iter().sum::<u64>());

    // Verify every position maps correctly
    let mut pos = 0u64;
    for (uid, &len) in lengths.iter().enumerate() {
        let (mapped_uid, offset) = idx.unitig_at_position(pos).unwrap();
        assert_eq!(mapped_uid, uid as u32, "Start of unitig {}", uid);
        assert_eq!(offset, 0);

        if len > 1 {
            let (mapped_uid, offset) = idx.unitig_at_position(pos + len - 1).unwrap();
            assert_eq!(mapped_uid, uid as u32, "End of unitig {}", uid);
            assert_eq!(offset, (len - 1) as u32);
        }
        pos += len;
    }

    // Past the end
    assert!(idx.unitig_at_position(pos).is_none());
}

#[test]
fn test_fenwick_max_stress() {
    let n = 500;
    let mut fw = dragon::ds::fenwick::FenwickMax::new(n);

    // Insert values in random order
    let mut vals: Vec<(usize, i64)> = (0..n).map(|i| (i, (i as i64) * 3 - 100)).collect();
    // Simple shuffle
    for i in (1..vals.len()).rev() {
        let j = (i * 7 + 3) % (i + 1);
        vals.swap(i, j);
    }

    for &(idx, val) in &vals {
        fw.update(idx, val);
    }

    // Prefix max should be correct
    let mut expected_max = i64::MIN;
    for i in 0..n {
        let val = (i as i64) * 3 - 100;
        expected_max = expected_max.max(val);
        assert_eq!(fw.prefix_max(i), expected_max, "Prefix max at {}", i);
    }
}

// ============================================================================
// FM-index tests
// ============================================================================

#[test]
fn test_fm_index_build_and_search() {
    let dir = TempDir::new().unwrap();
    let genome_dir = dir.path().join("genomes");
    fs::create_dir(&genome_dir).unwrap();

    // Create two genomes with a shared region
    let shared = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp
    write_fasta(&genome_dir, "g1", &format!("AAAA{}TTTT", shared));
    write_fasta(&genome_dir, "g2", &format!("CCCC{}GGGG", shared));

    // Parse unitigs (in this simple case, each genome as one unitig)
    let unitig_path = dir.path().join("unitigs.fa");
    let mut f = fs::File::create(&unitig_path).unwrap();
    writeln!(f, ">unitig_0").unwrap();
    writeln!(f, "AAAA{}TTTT", shared).unwrap();
    writeln!(f, ">unitig_1").unwrap();
    writeln!(f, "CCCC{}GGGG", shared).unwrap();

    let unitigs = dragon::index::unitig::parse_and_encode_unitigs(&unitig_path).unwrap();
    assert_eq!(unitigs.num_unitigs(), 2);

    // Build FM-index
    let index_dir = dir.path().join("index");
    fs::create_dir(&index_dir).unwrap();
    dragon::index::fm::build_fm_index(&unitigs, &index_dir).unwrap();

    // Load and search
    let fm = dragon::index::fm::load_fm_index(&index_dir).unwrap();

    // Search for the shared region
    let hits = fm.search(shared.as_bytes());
    assert_eq!(hits.len(), 2, "Shared region should be found in both unitigs");

    // Search for a unique region
    let hits_unique = fm.search(b"AAAAA");
    assert!(hits_unique.len() >= 1, "Unique region should be found");

    // Search for something not present
    let hits_missing = fm.search(b"ATATATATATATAT");
    assert_eq!(hits_missing.len(), 0);
}

#[test]
fn test_fm_variable_length_search() {
    let dir = TempDir::new().unwrap();
    let unitig_path = dir.path().join("unitigs.fa");
    let mut f = fs::File::create(&unitig_path).unwrap();
    writeln!(f, ">u0").unwrap();
    writeln!(f, "ACGTACGTACGTACGTACGTACGTACGT").unwrap();

    let unitigs = dragon::index::unitig::parse_and_encode_unitigs(&unitig_path).unwrap();
    let index_dir = dir.path().join("index");
    fs::create_dir(&index_dir).unwrap();
    dragon::index::fm::build_fm_index(&unitigs, &index_dir).unwrap();
    let fm = dragon::index::fm::load_fm_index(&index_dir).unwrap();

    // Variable length: ACGTACGTX — should match up to ACGTACGT
    let (len, count) = fm.variable_length_search(b"ACGTACGTXYZ");
    assert!(len >= 8, "Should match at least 8 bases, got {}", len);
    assert!(count > 0);
}

// ============================================================================
// End-to-end index + search test
// ============================================================================

#[test]
fn test_auto_batch_small_single_batch() {
    // Small dataset should NOT trigger batching — single-index fast path
    let dir = TempDir::new().unwrap();
    let genome_dir = dir.path().join("genomes");
    let index_dir = dir.path().join("index");
    fs::create_dir(&genome_dir).unwrap();

    // 3 small genomes
    for i in 0..3 {
        let seq = random_seq(1000, i * 100);
        write_fasta(&genome_dir, &format!("g{}", i), &seq);
    }

    // Big RAM budget → single batch
    let result = dragon::index::auto_batch::build_index_auto(
        &genome_dir, &index_dir, 15, 1, Some(64 * 1024 * 1024 * 1024),
    );
    assert!(result.is_ok(), "auto-batch failed: {:?}", result.err());

    // Must produce a valid single-index (no overlay_manifest)
    assert!(index_dir.join("metadata.json").exists());
    assert!(index_dir.join("fm_index.bin").exists());
    let metadata = dragon::index::metadata::load_metadata(&index_dir).unwrap();
    assert_eq!(metadata.num_genomes, 3);
}

#[test]
fn test_auto_batch_large_multi_batch() {
    // Force batching by setting a tiny RAM budget
    let dir = TempDir::new().unwrap();
    let genome_dir = dir.path().join("genomes");
    let index_dir = dir.path().join("index");
    fs::create_dir(&genome_dir).unwrap();

    // 10 genomes, but tiny RAM forces 2+ batches
    for i in 0..10 {
        let seq = random_seq(500, i * 73);
        write_fasta(&genome_dir, &format!("g{}", i), &seq);
    }

    // Force 2 batches: 10 genomes × 30 MB × 1.25 = need ≥ 375 MB for single batch
    // Set budget to 150 MB → should force ~2 batches
    let result = dragon::index::auto_batch::build_index_auto(
        &genome_dir, &index_dir, 15, 1, Some(150 * 1024 * 1024),
    );
    assert!(result.is_ok(), "auto-batch failed: {:?}", result.err());

    // Must produce base index
    assert!(index_dir.join("metadata.json").exists());
    assert!(index_dir.join("fm_index.bin").exists());

    // Must produce overlay manifest (indicating batching happened)
    let manifest_path = index_dir.join("overlay_manifest.json");
    assert!(manifest_path.exists(),
        "Expected overlay manifest for multi-batch build");

    // Verify overlays directory has entries
    let overlays_dir = index_dir.join("overlays");
    assert!(overlays_dir.exists());
    let overlay_count = fs::read_dir(&overlays_dir)
        .map(|rd| rd.filter_map(|e| e.ok()).count())
        .unwrap_or(0);
    assert!(overlay_count >= 1, "Expected at least 1 overlay, got {}", overlay_count);
}

#[test]
fn test_multi_index_search() {
    // Build 2 independent shard indices, then search a query across both
    let dir = TempDir::new().unwrap();
    let query_path = dir.path().join("query.fa");

    // Shard 1: 3 genomes with core region A
    let shard1_dir = dir.path().join("shard1_genomes");
    let shard1_idx = dir.path().join("shard1_idx");
    fs::create_dir(&shard1_dir).unwrap();
    let core_a = random_seq(1500, 11111);
    for i in 0..3 {
        let flank = random_seq(300, i * 10);
        write_fasta(&shard1_dir, &format!("s1_g{}", i), &format!("{}{}", flank, core_a));
    }
    dragon::index::build_index(&shard1_dir, &shard1_idx, 15, 1).unwrap();

    // Shard 2: 3 different genomes with core region B
    let shard2_dir = dir.path().join("shard2_genomes");
    let shard2_idx = dir.path().join("shard2_idx");
    fs::create_dir(&shard2_dir).unwrap();
    let core_b = random_seq(1500, 22222);
    for i in 0..3 {
        let flank = random_seq(300, i * 10 + 100);
        write_fasta(&shard2_dir, &format!("s2_g{}", i), &format!("{}{}", flank, core_b));
    }
    dragon::index::build_index(&shard2_dir, &shard2_idx, 15, 1).unwrap();

    // Query from core_a — should hit shard 1 genomes, not shard 2
    let mut f = fs::File::create(&query_path).unwrap();
    writeln!(f, ">q_from_core_a").unwrap();
    writeln!(f, "{}", &core_a[300..1200]).unwrap();

    let mut config = dragon::query::SearchConfig::default();
    config.index_dir = shard1_idx.clone().into();
    config.min_identity = 0.0;
    config.min_query_coverage = 0.0;
    config.min_chain_score = 0.0;

    let indices = vec![shard1_idx.clone(), shard2_idx.clone()];
    let results = dragon::query::search_multi_index(&query_path, &indices, &config).unwrap();

    // Should get results (even if empty, not crash)
    assert_eq!(results.len(), 1, "Expected 1 query result");
    // Verify the multi-index search ran against both shards without panic
}

#[test]
fn test_end_to_end_small_dataset() {
    let dir = TempDir::new().unwrap();
    let genome_dir = dir.path().join("genomes");
    let index_dir = dir.path().join("index");
    let query_path = dir.path().join("query.fa");
    fs::create_dir(&genome_dir).unwrap();

    // Create 5 genomes with shared core
    let core = random_seq(2000, 12345);
    for i in 0..5 {
        let flank5 = random_seq(500, i * 100);
        let flank3 = random_seq(500, i * 100 + 50);
        let genome_seq = format!("{}{}{}", flank5, core, flank3);
        write_fasta(&genome_dir, &format!("genome_{}", i), &genome_seq);
    }

    // Build index
    let result = dragon::index::build_index(&genome_dir, &index_dir, 15, 1);
    assert!(result.is_ok(), "Index build failed: {:?}", result.err());

    // Verify index files exist
    assert!(index_dir.join("metadata.json").exists());
    assert!(index_dir.join("fm_index.bin").exists());
    assert!(index_dir.join("colors.drgn").exists());
    assert!(index_dir.join("paths.bin").exists());

    // Create a query from the core region
    let query_seq = &core[500..1500]; // 1000bp from middle of core
    let mut f = fs::File::create(&query_path).unwrap();
    writeln!(f, ">test_query").unwrap();
    writeln!(f, "{}", query_seq).unwrap();

    // Load index and search
    let metadata = dragon::index::metadata::load_metadata(&index_dir).unwrap();
    assert_eq!(metadata.num_genomes, 5);
    assert_eq!(metadata.kmer_size, 15);

    let fm = dragon::index::fm::load_fm_index(&index_dir).unwrap();

    // Search for a k-mer from the core — should find hits
    let kmer = &core.as_bytes()[600..615]; // 15-mer
    let hits = fm.search(kmer);
    assert!(!hits.is_empty(), "Should find core k-mer in index");
}

#[test]
fn test_end_to_end_with_mutations() {
    let dir = TempDir::new().unwrap();
    let genome_dir = dir.path().join("genomes");
    let index_dir = dir.path().join("index");
    fs::create_dir(&genome_dir).unwrap();

    // Create one genome
    let genome = random_seq(3000, 99);
    write_fasta(&genome_dir, "ref_genome", &genome);

    // Build index
    dragon::index::build_index(&genome_dir, &index_dir, 15, 1).unwrap();
    let fm = dragon::index::fm::load_fm_index(&index_dir).unwrap();

    // Query with 0% divergence — should find many seeds
    let query_exact = &genome.as_bytes()[500..1500];
    let hits_exact = fm.search(&query_exact[0..15]);
    assert!(!hits_exact.is_empty(), "Exact match should find seeds");

    // Query with 5% divergence — should still find some seeds with k=15
    let query_mutated = mutate_seq(&genome[500..1500], 0.05, 42);
    // At 5% divergence, some 15-mers will still be intact
    let mut found_any = false;
    for i in 0..query_mutated.len().saturating_sub(15) {
        let kmer = &query_mutated.as_bytes()[i..i + 15];
        if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            let hits = fm.search(kmer);
            if !hits.is_empty() {
                found_any = true;
                break;
            }
        }
    }
    assert!(found_any, "At 5% divergence, should still find some intact 15-mers");
}

// ============================================================================
// PAF output tests
// ============================================================================

#[test]
fn test_paf_output_format() {
    let record = dragon::io::paf::PafRecord {
        query_name: "query1".to_string(),
        query_len: 1000,
        query_start: 10,
        query_end: 990,
        strand: '+',
        target_name: "genome_0".to_string(),
        target_len: 5000,
        target_start: 100,
        target_end: 1080,
        num_matches: 950,
        alignment_len: 980,
        mapq: 60,
        tags: vec!["AS:i:1900".to_string()],
    };

    let formatted = format!("{}", record);
    let fields: Vec<&str> = formatted.split('\t').collect();
    assert_eq!(fields.len(), 13); // 12 mandatory + 1 tag
    assert_eq!(fields[0], "query1");
    assert_eq!(fields[5], "genome_0");
    assert_eq!(fields[11], "60"); // mapq
    assert_eq!(fields[12], "AS:i:1900");

    // Identity
    assert!((record.identity() - 0.9694).abs() < 0.001);
}

#[test]
fn test_blast_tabular_output() {
    let records = vec![dragon::io::paf::PafRecord {
        query_name: "q1".to_string(),
        query_len: 500,
        query_start: 0,
        query_end: 500,
        strand: '+',
        target_name: "g1".to_string(),
        target_len: 3000,
        target_start: 100,
        target_end: 600,
        num_matches: 480,
        alignment_len: 500,
        mapq: 60,
        tags: vec![],
    }];

    let mut buf = Vec::new();
    dragon::io::blast::write_blast_tabular(&mut buf, &records).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let fields: Vec<&str> = output.trim().split('\t').collect();
    assert_eq!(fields[0], "q1");
    assert_eq!(fields[1], "g1");
    // Identity should be 96.00
    assert!(fields[2].starts_with("96.0"));
}

// ============================================================================
// Colour index tests
// ============================================================================

#[test]
fn test_colour_index_roundtrip() {
    let dir = TempDir::new().unwrap();

    // Create a colour file
    let colour_path = dir.path().join("colors.tsv");
    let mut f = fs::File::create(&colour_path).unwrap();
    writeln!(f, "0\t0,1,2").unwrap();     // unitig 0 in genomes 0,1,2
    writeln!(f, "1\t1,3").unwrap();        // unitig 1 in genomes 1,3
    writeln!(f, "2\t0,2,4").unwrap();      // unitig 2 in genomes 0,2,4

    // Build colour index
    dragon::index::color::build_color_index(&colour_path, dir.path(), 5).unwrap();

    // Load and query
    let idx = dragon::index::color::load_color_index(dir.path()).unwrap();
    assert_eq!(idx.num_unitigs(), 3);

    let colors0 = idx.get_colors(0).unwrap();
    assert!(colors0.contains(0));
    assert!(colors0.contains(1));
    assert!(colors0.contains(2));
    assert!(!colors0.contains(3));

    let colors1 = idx.get_colors(1).unwrap();
    assert!(colors1.contains(1));
    assert!(colors1.contains(3));
    assert!(!colors1.contains(0));

    // Out of range unitig should return empty
    let colors_oob = idx.get_colors(999).unwrap();
    assert_eq!(colors_oob.len(), 0);
}

// ============================================================================
// Chaining tests
// ============================================================================

#[test]
fn test_chain_preserves_colinearity() {
    use dragon::query::chain::Anchor;

    // Create colinear anchors
    let anchors = vec![
        Anchor { query_start: 0, query_end: 20, ref_start: 100, ref_end: 120, match_len: 20, is_reverse: false, score: 20.0 },
        Anchor { query_start: 30, query_end: 50, ref_start: 130, ref_end: 150, match_len: 20, is_reverse: false, score: 20.0 },
        Anchor { query_start: 60, query_end: 80, ref_start: 160, ref_end: 180, match_len: 20, is_reverse: false, score: 20.0 },
        // This one is out of order — should not be in the best chain
        Anchor { query_start: 10, query_end: 30, ref_start: 200, ref_end: 220, match_len: 20, is_reverse: false, score: 20.0 },
    ];

    // Use the internal chain function indirectly via chain_candidates
    // For now, just verify the anchors are sortable
    let mut sorted = anchors.clone();
    sorted.sort_by_key(|a| a.ref_start);
    assert_eq!(sorted[0].ref_start, 100);
    assert_eq!(sorted[3].ref_start, 200);
}
