//! Micro-benchmark: isolate `get_path` decode cost for v2 vs v3 path
//! formats. Both are opened in the same process so node/cache variance
//! is controlled for.
//!
//! Usage:
//!   bench_get_path <index_dir>
//!
//! Expects `<index_dir>/paths.bin` (v3) and `<index_dir>/paths.bin.v2.bak`
//! (v2). Calls get_path over the same pseudo-random genome IDs for each
//! and reports wall time + steps decoded.

use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;

use dragon::index::paths_v2::MmapPathIndex;
use dragon::index::paths_v3::MmapPathIndexV3;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: bench_get_path <index_dir>");
        std::process::exit(2);
    }
    let dir = PathBuf::from(&args[1]);
    let v2_path = dir.join("paths.bin.v2.bak");
    let v3_path = dir.join("paths.bin");

    let n_calls: usize = 5000;

    let v2 = MmapPathIndex::open(&v2_path)?;
    let v3 = MmapPathIndexV3::open(&v3_path)?;
    let num_genomes = v2.num_genomes();
    println!("genomes: v2={} v3={}", num_genomes, v3.num_genomes());

    // Pseudo-random genome IDs (fixed, identical for both formats).
    let ids: Vec<u32> = (0..n_calls)
        .map(|i| ((i as u64 * 7919) % num_genomes) as u32)
        .collect();

    for round in 0..3 {
        let t = Instant::now();
        let mut total_steps: u64 = 0;
        for &g in &ids {
            let p = v2.get_path(g)?.expect("v2 genome");
            total_steps += p.steps.len() as u64;
        }
        let el = t.elapsed();
        println!(
            "v2 round{round}: {:?}  {:.1} us/call  ({} steps total)",
            el,
            el.as_micros() as f64 / n_calls as f64,
            total_steps
        );
    }

    for round in 0..3 {
        let t = Instant::now();
        let mut total_steps: u64 = 0;
        for &g in &ids {
            let p = v3.get_path(g)?.expect("v3 genome");
            total_steps += p.steps.len() as u64;
        }
        let el = t.elapsed();
        println!(
            "v3 round{round}: {:?}  {:.1} us/call  ({} steps total)",
            el,
            el.as_micros() as f64 / n_calls as f64,
            total_steps
        );
    }
    Ok(())
}
