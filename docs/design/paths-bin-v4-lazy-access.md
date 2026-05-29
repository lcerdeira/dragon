# Design: `paths.bin` v4 — lazy path access (Phase 2/3)

**Status:** Phase 1 implemented and committed (`paths_v4.rs`, tests green).
Phase 2 (consumer rewrite) and Phase 3 (HPC migration + benchmark) are
specified below for a focused fresh-session implementation.

**Parent design:** `docs/design/paths-bin-v3.md` (the shipped v3 graph-edge
format). This doc is the speed follow-up.

---

## Why

v3 shrank `paths.bin` 41 GB → 21 GB per shard (index 343 → 204 GB) with
byte-identical accuracy. But search slowed 2.5× on the 200-query batch
(1041 s → 2586 s). Micro-benchmark on a real shard pinned the cause
exactly: warm `MmapPathIndexV3::get_path` takes **98 ms/call** (vs
**3.4 ms/call** for v2) because reconstructing each step random-accesses
the **60 MB** CSR successor table — larger than L3 cache, so ~1.3 M cache
misses per call × 80 ns ≈ 100 ms. Loading the table into RAM at open did
not help (the bottleneck is cache, not network).

For the laptop mission per-query cost ≈13 s is already acceptable, but
the batch number costs Dragon its speed lead over LexicMap (2586 s vs
2897 s — Dragon barely faster). Phase 2 closes that.

## The puzzle, and the fix

`direct_align` uses the full path for two things, both **before** the
alignment window is known:

1. Strand partition — needs `step.is_reverse` for every seed's unitig.
2. Anchor coordinates — needs `step.genome_offset` for every seed's
   unitig.

So `iter_window(gid, start, end)` alone is not enough — you don't know
the window until you've anchored the seeds. The fix is a different new
method that queries **by unitig set**:

```rust
impl MmapPathIndexV4 {
    /// Walk the genome's step stream from the start and collect, for each
    /// unitig_id in `wanted`, the FIRST PathStep that matches it. Returns
    /// once every wanted unitig has been found OR the path ends. For true
    /// matches the seeds cluster near one genome region, so the walk
    /// terminates well before the full 650k-step decode.
    pub fn find_unitig_steps(
        &self,
        genome_id: u32,
        wanted: &HashSet<u32>,
    ) -> Result<HashMap<u32, PathStep>>;
}
```

`PathIndex` (the enum) grows a corresponding method that dispatches to
the v4 variant; v2/v3/Eager variants implement it via their existing
`get_path` + a linear scan (no speed regression for them).

## Phase 2 — consumer rewrite

### `direct_align_candidates` (`src/query/direct_align.rs`)

Currently:
```rust
let genome_path = path_index.get_path(hit.genome_id)?;            // ~100 ms in v3
let mut unitig_step: HashMap<u32, &PathStep> = ...;               // built from ALL 650k steps
// strand partition, then:
align_with_seeds(query, ..., &genome_path, &unitig_step, ...);
```

After Phase 2:
```rust
let wanted: HashSet<u32> = hit.seeds.iter().map(|s| s.unitig_id).collect();
let unitig_step = path_index.find_unitig_steps(hit.genome_id, &wanted)?;
let (genome_name, genome_length) = path_index.genome_meta(hit.genome_id)?;
// strand partition (uses unitig_step.get(&s.unitig_id)) — unchanged
align_with_seeds(query, ..., genome_name, genome_length, &unitig_step, path_index, hit.genome_id, ...);
```

`get_path` is no longer on the hot path.

### `align_with_seeds`

Currently calls `extract_sequence_static(genome_path, start, end, unitigs)`
— iterates `path.steps` to find steps overlapping `[start, end)`.

After Phase 2:
```rust
extract_sequence(path_index, genome_id, start, end, unitigs) {
    let steps = path_index.iter_window(genome_id, start, end)?;
    // existing emitted_up_to logic over those steps
}
```

`extract_sequence_static` becomes a thin wrapper around `iter_window`. For
v4 this decodes ~window/4.3 bp/step ≈ a few hundred steps (well under
1 ms). For v2/v3 it falls back to filtering the full path (no regression).

### `PathIndex` enum & dispatch

Add the methods on the enum:

```rust
impl PathIndex {
    pub fn find_unitig_steps(&self, genome_id: u32, wanted: &HashSet<u32>) -> ...;
    pub fn iter_window(&self, genome_id: u32, start: u64, end: u64) -> ...;
    pub fn genome_meta(&self, genome_id: u32) -> Option<(String, u64)>;
}
```

For variants that lack native support, implement via existing `get_path`.

### Tests

- Unit: `find_unitig_steps` returns the same `PathStep`s that
  `get_path`-then-`HashMap::get` would, for synthetic and real paths.
- Integration: `integration_real_genome` continues to pass with v4 paths
  (search returns identity 1.000 / target_start 1999).
- A new bench `bench_direct_align`: time direct_align's per-candidate
  cost on a real shard for v2 vs v3 vs v4; target v4 ≈ v2.

## Phase 3 — HPC deployment

1. Add `src/bin/migrate_paths_v4` (mirrors `migrate_paths_v3`): in-place
   migrate `paths.bin` → v4, backup as `paths.bin.v3.bak`.
2. SLURM job: migrate the 7 saureus shards (one-pass, faster than v3
   migration since v3 → v4 is a simple re-encode with checkpoints).
3. Run the saureus benchmark on v4. Targets:
   - `mean_id` ≥ 0.996, perfect-rate 0.980 (must hold).
   - search wall ≤ 1300 s (back near v2's 1041 s).
   - index disk: same as v3 (~204 GB) plus per-blob checkpoint overhead
     (~100 MB/shard, negligible).
4. If all targets hit: delete the .v2.bak and .v3.bak files. Index is
   v4 going forward.

## Risks & mitigations

| risk | mitigation |
|---|---|
| `find_unitig_steps` early-exit misses a relevant seed | wanted set is built from `hit.seeds`; if a seed's unitig is NOT in the genome's path the lookup correctly returns nothing (no false anchor); pathological case = silent unitig absence, which v3 also has |
| Pathological queries with seeds spread across the genome | early-exit degrades to full walk — same as v3's `get_path` worst case, no new regression |
| `align_with_seeds` change breaks alignment correctness | `integration_real_genome` + `integration_align` + the saureus benchmark gate every commit; never push without all green |
| v4 migration on HPC takes hours | one-pass re-encode; estimate ~30 min/shard × 7 = ~3.5 h; runs as a background SLURM job |

## Out of scope (separate)

- `colors.drgn` shrink (the other 51 GB) — needed to fully reach <100 GB.
- A true BWT+rank FM-index — would shrink `fm_index.bin` (6 GB) and cut
  search RAM further.
