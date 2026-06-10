# Dragon Indexes — Manifest

> **DRAFT.** Regenerate the data columns with [`scripts/gen_manifest.sh`](scripts/gen_manifest.sh).
> Checksums and final genome counts are filled at publication time.

- **Bucket:** `s3://dragon-zarr/` · **Region:** `eu-west-2` · **Format:** Zarr (chunked)
- **Snapshot date:** «YYYY-MM-DD»
- **Dragon version:** «dragon --version» · **k-mer:** 31 (canonical) · **build preset:** «default/amr»

---

## Genome-collection indexes

### AllTheBacteria (ATB) — ✅ Published
- **Description:** Index over the full AllTheBacteria v«0.2» corpus of quality-controlled bacterial/archaeal assemblies.
- **Genomes:** 2,292,989
- **S3:** `s3://dragon-zarr/atb/` — 1,147 shards (`shard_0000.zarr` … `shard_1146.zarr`)
- **Objects:** 9,999,899 · **Size:** 8,099,470,677,679 bytes (≈ 8.1 TB)
- **Source data:** AllTheBacteria (Hunt et al.) — <https://allthebacteria.org/>
- **Build:** `dragon index -k 31 --threads 32 --low-memory` then `dragon export-zarr`
- **Checksum manifest:** `s3://dragon-zarr/atb/SHA256SUMS` — [TBD generate]
- **License:** open (INSDC-derived) «confirmar»

### GTDB — 🟡 Partial
- **Genomes:** «[TBD] (release r«220»)»
- **S3:** `s3://dragon-zarr/gtdb_batch1.zarr` (batch 1 only; batches 2–«N» pending upload)
- **Objects:** 32,749 · **Size:** 19,136,407,977 bytes (≈ 19.1 GB) *(batch 1)*
- **Note:** local raw indexes `gtdb_index_batch2..6` exist on HPC but are **not yet exported/uploaded**.
- **License:** «CC-BY-4.0 confirmar»

### Klebsiella pneumoniae — 🟡 Partial
- **Genomes:** «[TBD]»
- **S3:** `s3://dragon-zarr/kpneumo/b1` (batch 2 pending)
- **Objects:** 27,082 · **Size:** 6,898,920,254 bytes (≈ 6.9 GB)

### Staphylococcus aureus — 🟡 Partial
- **Genomes:** «[TBD]»
- **S3:** `s3://dragon-zarr/saureus/b1` + `s3://dragon-zarr/saureus_batch7.zarr`
- **Objects:** 8,199 + 5,470 = 13,669 · **Size:** 1,218,172,654 + 801,767,419 ≈ 2.0 GB

### Escherichia coli — 🔄 Building
- **Genomes:** ~364,000 (13 batches × ~28,000)
- **S3:** `s3://dragon-zarr/ecoli/` (prefix currently empty)
- **Build job:** SLURM array `6138675` `[1-13]%3` on `loginhpc` (started 2026-06-09)
- **Status:** 0/13 batches uploaded as of «timestamp»

### Shigella · Salmonella enterica · Salmonella Typhi — ⏳ Planned
- **S3 (proposed):** `shigella/`, `salmonella_enterica/`, `salmonella_typhi/`
- **Plan:** per-species ATB subset → `dragon index` → `export-zarr` → upload.

---

## Reference gene-database indexes

| DB | Records | S3 (proposed) | Source | License note |
|---|---|---|---|---|
| **AMRFinderPlus** | 25.6M hit rows / 2,440,377 genomes | `amrfinderplus.parquet` | ATB (OSF), NCBI AMRFinderPlus | NCBI public domain |
| **CARD** | «~5,000 genes» | `card/` | <https://card.mcmaster.ca> | ⚠️ **restricted redistribution** — may need to ship derived-only |
| **VFDB** | «~28,000 genes» | `vfdb/` | <http://www.mgc.ac.cn/VFs/> | academic use «confirmar» |
| **PlasmidFinder** | «~2,000 genes» | `plasmidfinder/` | DTU-CGE | open «confirmar» |

---

## Open items before publication
1. Fill all `«[TBD]»` genome counts and the GTDB/«…» release versions.
2. Generate and upload `SHA256SUMS` per index.
3. Finish + upload pending uploads: GTDB batches 2–N, kpneumo batch 2, E. coli (#9), and build Shigella/Salmonella×2.
4. Resolve **CARD** redistribution licensing.
5. Decide whether to also publish per-index `.tar.zst` for one-click download.
