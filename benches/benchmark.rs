use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_dna_encoding(c: &mut Criterion) {
    let seq: Vec<u8> = (0..10000).map(|i| b"ACGT"[i % 4]).collect();

    c.bench_function("pack_10k_bases", |b| {
        b.iter(|| {
            dragon::util::dna::PackedSequence::from_bytes(black_box(&seq))
        })
    });

    let packed = dragon::util::dna::PackedSequence::from_bytes(&seq);
    c.bench_function("unpack_10k_bases", |b| {
        b.iter(|| black_box(&packed).to_bytes())
    });

    c.bench_function("revcomp_10k_bases", |b| {
        b.iter(|| black_box(&packed).reverse_complement())
    });
}

fn bench_varint(c: &mut Criterion) {
    let values: Vec<u64> = (0..10000).map(|i| i * 137).collect();

    c.bench_function("varint_encode_10k", |b| {
        b.iter(|| dragon::ds::varint::encode_varint_slice(black_box(&values)))
    });

    let encoded = dragon::ds::varint::encode_varint_slice(&values);
    c.bench_function("varint_decode_10k", |b| {
        b.iter(|| dragon::ds::varint::decode_varint_slice(black_box(&encoded), 10000))
    });
}

fn bench_fenwick(c: &mut Criterion) {
    c.bench_function("fenwick_max_10k_ops", |b| {
        b.iter(|| {
            let mut fw = dragon::ds::fenwick::FenwickMax::new(10000);
            for i in 0..10000 {
                fw.update(i, i as i64);
            }
            for i in 0..10000 {
                black_box(fw.prefix_max(i));
            }
        })
    });
}

criterion_group!(benches, bench_dna_encoding, bench_varint, bench_fenwick);
criterion_main!(benches);
