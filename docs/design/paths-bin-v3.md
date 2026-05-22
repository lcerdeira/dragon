# Design: `paths.bin` v3 — graph-edge encoding

**Status:** proposed, not yet implemented
**Author:** prepared for a focused implementation session
**Goal:** shrink the genome path index so the full 104k-genome *S. aureus*
database fits the "runs on a laptop" target — 8–16 GB RAM, <100 GB disk,
4–8 CPUs.

---

## 1. Motivation

Dragon's mission: be the **low-resource alternative** to LexicMap so
researchers without HPC access can search large genome databases on a
laptop.

Current standing on the laptop target (104k-genome saureus index, 7 shards):

| Resource | Laptop target | Dragon today | Status |
|---|---|---|---|
| Search RAM | 8–16 GB | ~2–3 GB | ✅ met (one shard resident at a time) |
| CPUs | 4–8 | scales down cleanly | ✅ met |
| **Index on disk** | <100 GB | **343 GB** | ❌ **the blocker** |

Per-shard disk breakdown:

| file | per shard | ×7 |
|---|---|---|
| `paths.bin` | 41 GB | **287 GB** |
| `colors.drgn` | 7.3 GB | 51 GB |
| `fm_index.bin` | 0.86 GB | 6 GB |
| `specificity.drgn` | 3.5 MB | negligible |

`paths.bin` is 84% of the index. This document specifies how to shrink it.
A separate follow-up (`colors.drgn`) is needed to fully clear <100 GB.

## 2. Why `paths.bin` is 41 GB (measured)

`paths.bin` v2 stores, per genome, the explicit walk through the colored
de Bruijn graph. Measured on `saureus_index_batch1` (16,000 genomes):

| quantity | value |
|---|---|
| steps per genome | **~650,000** |
| genome length | ~2.8 Mbp |
| **bp advanced per step** | **~4.3** |
| bytes per step | 4.24 |
| of which: unitig-id bytes | **76%** (~3.24 B/step) |
| of which: delta-offset bytes | ~24% (~1.0 B/step) |

A *step* is one unitig in the genome's graph walk. With 16k *S. aureus*
strains in one graph, inter-strain SNPs fragment the graph so heavily that
the average walked unitig is only ~34 bp — hence ~650k steps per 2.8 Mbp
genome. This is **inherent to the colored-dBG representation**, not a bug.

The **unitig IDs dominate** (76%): each step stores a full
`varint(unitig_id << 1 | is_reverse)` ≈ 3.24 bytes, because unitig IDs run
to ~2.07M (21 bits).

v2 wire format (see `src/index/paths_v2.rs`) per genome blob:

```
varint name_length, bytes name
varint genome_length
varint num_steps
repeat num_steps:
  varint unitig_with_rev   (unitig_id << 1 | is_reverse)
  varint delta_offset      (this genome_offset - previous)
```

## 3. Key idea: store graph edges, not unitig IDs

Consecutive unitigs in a genome walk are connected by an **edge** in the
de Bruijn graph. A unitig end has out-degree ≤ 4 (one per next nucleotide);
in practice the *observed* successor set of a `(unitig, orientation)` across
all 16k genomes is small (typically 1–4, rarely up to ~8).

So instead of storing the next unitig's full 21-bit ID (~3.24 B), store the
**edge index** — which successor of the current node to follow (~1 byte
varint, value usually 0–3; or 4 bits packed).

A **shared successor table** (built once per shard) lets the reader
reconstruct full unitig IDs by following edges. The table is tiny relative
to the savings.

### Expected size

| encoding | bytes/step | paths/shard | index total (7) |
|---|---|---|---|
| v2 (today) | 4.24 | 41 GB | 343 GB |
| **v3 simple** (1 B edge varint + 1 B delta varint) | ~2.0 | ~21 GB | ~200 GB |
| **v3 packed** (4-bit edge + 4-bit delta, escape for outliers) | ~1.0 | ~10 GB | ~127 GB |

Recommendation: **implement v3-simple first** (correctness-critical, keep it
boring), land it, validate, then optionally add the nibble-packing as a
follow-up. v3-simple already halves the index.

`colors.drgn` work (separate doc) is required to cross <100 GB on the full
104k-genome database. With v3-packed + a colors pass the projection is
~90–100 GB.

## 4. v3 wire format (`DRGNPTH3`)

All multi-byte integers little-endian. Varints are LEB128 (as in v2).

```
Header (48 bytes):
  +0   magic[8]              = b"DRGNPTH3"
  +8   version: u32          = 1
  +12  num_genomes: u64
  +20  num_nodes: u64        = 2 * num_unitigs   (one node per orientation)
  +28  succ_offsets_pos: u64 (absolute byte offset of successor offset table)
  +36  succ_values_pos: u64  (absolute byte offset of successor value array)
  +44  genome_table_pos: u32 ... NOTE: keep header 8-byte aligned; use u64.
       (use +44 reserved u32 + put genome_table_pos earlier, OR widen header
        to 56 bytes — see implementation note below)

Successor table (CSR):
  succ_offsets: [u64; num_nodes + 1]   at succ_offsets_pos
  succ_values:  [u32; total_edges]     at succ_values_pos
     node id      = (unitig_id << 1) | is_reverse
     succ of node = succ_values[succ_offsets[node] .. succ_offsets[node+1]]
     each value   = (next_unitig_id << 1) | next_is_reverse

Genome offset table:
  genome_offsets: [u64; num_genomes + 1]   (absolute offsets of blobs)

Per-genome blob:
  varint name_length, bytes name
  varint genome_length
  varint num_steps
  if num_steps >= 1:
    varint first_node          ((unitig_id<<1)|is_reverse of step 0)
    varint first_genome_offset (genome_offset of step 0)
  repeat (num_steps - 1):
    varint edge_index          (index into current node's successor slice)
    varint delta_offset        (genome_offset - previous genome_offset)
```

**Implementation note on the header:** keep it simple — use a fixed
**56-byte** header so every field is 8-byte aligned:
`magic[8], version u32, _pad u32, num_genomes u64, num_unitigs u64,
succ_offsets_pos u64, succ_values_pos u64, genome_table_pos u64`.

## 5. Successor table construction (migration, no GGCAT rerun)

The table is built by scanning the existing v2 paths — every transition any
genome needs is observed by construction, so no graph recomputation from
unitig sequences is required.

```
# Pass 1: collect observed transitions
succ: HashMap<u32 node, BTreeSet<u32 next_node>>
for each genome g (decode v2 GenomePath):
    for i in 0 .. g.steps.len()-1:
        a = (g.steps[i].unitig_id   << 1) | g.steps[i].is_reverse
        b = (g.steps[i+1].unitig_id << 1) | g.steps[i+1].is_reverse
        succ[a].insert(b)

# Freeze to CSR: succ_offsets[node], succ_values[]
# (sorted successor lists -> edge_index is the position in the sorted slice)
```

`edge_index` for a transition `a -> b` = position of `b` in the **sorted**
successor slice of `a`. Sorted order makes encode/decode agree
deterministically.

Sanity check during migration: every observed `b` MUST be findable in
`succ[a]` (it will be — same scan). Assert it.

## 6. Reader / decode (`get_path`)

```
decode genome blob:
  read name, genome_length, num_steps
  if num_steps == 0: return GenomePath{ steps: [] }
  node   = first_node
  offset = first_genome_offset
  steps  = [ PathStep{ unitig_id: node>>1, is_reverse: node&1, genome_offset: offset } ]
  for _ in 1 .. num_steps:
    edge_index = read varint
    delta      = read varint
    slice      = succ_values[ succ_offsets[node] .. succ_offsets[node+1] ]
    node       = slice[edge_index]          # bounds-check!
    offset    += delta
    steps.push( PathStep{ unitig_id: node>>1, is_reverse: node&1, genome_offset: offset } )
  return GenomePath{ ..., steps }
```

`get_path` returns the **identical `GenomePath`** as v2, so every downstream
consumer (`extract_sequence_static`, `align_with_seeds`, `dump_path`, …) is
**unchanged**. The abstraction boundary holds at `PathIndex::get_path`.

## 7. Integration points

- `src/index/paths_v3.rs` — new module: `PathV3Writer`, `MmapPathIndexV3`,
  `migrate_v2_to_v3`, `is_v3`.
- `src/index/paths.rs`:
  - `enum PathIndex` — add `MmapV3(Arc<MmapPathIndexV3>)`.
  - `get_path`, `num_genomes`, `iter` — add the `MmapV3` arm.
  - `load_path_index` — dispatch on magic: `DRGNPTH3` → v3, `DRGNPTH2` → v2,
    else legacy bincode.
- `build_path_index` — optionally emit v3 directly (later; migration covers
  existing indexes first).

## 8. Testing plan (the safety net — do this BEFORE any HPC reindex)

A wrong path decode silently corrupts **every** alignment. Tests are
non-negotiable:

1. **Unit round-trip** (`paths_v3::tests`): synthetic genomes → v3 write →
   v3 read → assert `GenomePath` equality (unitig_id, is_reverse,
   genome_offset for every step).
2. **Migration round-trip on the real single-genome index**
   (`tests/data/sg_index`, already in the repo, git-ignored):
   - `migrate_v2_to_v3(sg_index/paths.bin)` → `paths.bin.v3`.
   - For every genome, assert `v3.get_path(g) == v2.get_path(g)` exactly.
3. **End-to-end**: `integration_real_genome` must still pass with the v3
   file in place (`SAMEA110247553` self-slice → identity 1.000,
   target_start ≈ 1999).
4. **On HPC, before trusting it**: migrate ONE shard, run the saureus
   benchmark, confirm `mean_id` and perfect-rate are unchanged vs the v2
   baseline (0.9966 / 0.980).

## 9. Reindex plan (HPC)

- `migrate_v2_to_v3` reads a v2 `paths.bin` (~41 GB), writes v3 (~10–21 GB).
  No GGCAT, no genome rewalk.
- SLURM job, one task per shard or sequential; ~287 GB read + ~100 GB write
  total. Estimate 1–3 h wall.
- Keep the v2 file as `paths.bin.v2.bak` until the v3 benchmark passes;
  then reclaim.

## 10. Risks & mitigations

| risk | mitigation |
|---|---|
| Wrong decode → silent alignment corruption | Round-trip tests (§8) gate every step; HPC single-shard benchmark before full rollout |
| Successor slice index out of range | Bounds-check in decode; assert during migration that every transition is in the table |
| `edge_index` varint larger than expected (high-degree nodes) | varint handles it; measure max out-degree during migration and log it |
| Header alignment / format drift | Fixed 56-byte aligned header; magic + version checked on open |
| Disk during migration (need v2 + v3 simultaneously) | ~140 GB extra transient; verify free space, migrate shard-by-shard |

## 11. Out of scope (separate follow-up)

- `colors.drgn` (51 GB total) — needed to fully reach <100 GB on the 104k
  database. Likely approach: evaluate Roaring run-containers / a denser
  encoding for high-cardinality unitig color sets.
- `fm_index.bin` suffix array (744 MB/shard, slurped into RAM) — a true
  BWT+rank FM-index would cut both this file and search RAM; large change.
- Emitting v3 directly from `build_path_index` (migration covers existing
  indexes first).
