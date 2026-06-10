# NOTICE — sources, attribution & redistribution status

The Dragon Indexes are **derived works** built from third-party databases. Each source keeps
its own license. This file records attribution and whether each derived index is
**redistributed** here (under CC-BY-4.0) or must be **built locally** by the user.

| Source | Used for | Source license | Redistributed here? | How to get the index |
|---|---|---|---|---|
| **AllTheBacteria** (Hunt et al.) | ATB + per-species genome indexes | open (INSDC-derived) | ✅ Yes — CC-BY-4.0 | `s3://dragon-zarr/atb/`, `…/ecoli/`, etc. |
| **INSDC** (ENA/NCBI/DDBJ) | underlying assemblies | public / unrestricted | ✅ Yes — CC-BY-4.0 | as above |
| **GTDB** | GTDB index | **CC-BY-4.0** | ✅ Yes — CC-BY-4.0 (attribution to GTDB) | `s3://dragon-zarr/gtdb_*` |
| **NCBI AMRFinderPlus** | AMR index / annotations | U.S. Gov. **public domain** | ✅ Yes — CC-BY-4.0 | `s3://dragon-zarr/amrfinderplus.parquet` |
| **PlasmidFinder** (DTU-CGE) | plasmid index | free / ~Apache-2.0 (verify) | ✅ Yes — CC-BY-4.0 + attribution | `s3://dragon-zarr/plasmidfinder/` |
| **VFDB** (Chen et al.) | virulence index | "free for academic use", no clear open redistribution | ⚠️ **Not yet** — pending author confirmation | `scripts/build_vfdb_index.sh` (build locally) |
| **CARD** (Alcock et al., McMaster) | AMR index | **Data Usage Agreement — redistribution prohibited; non-commercial** | ❌ **No** | `scripts/build_card_index.sh` (build locally after accepting CARD terms) |

## Attribution / citations for sources
- **AllTheBacteria** — Hunt M, et al. *AllTheBacteria — all bacterial genomes assembled,
  available and searchable.* bioRxiv 2024. doi:10.1101/2024.03.08.584059
- **GTDB** — Parks DH, et al. *GTDB: an ongoing census of bacterial and archaeal diversity…*
  Nucleic Acids Research. <https://gtdb.ecogenomic.org/> (CC-BY-4.0)
- **AMRFinderPlus** — Feldgarden M, et al. NCBI. <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>
- **CARD** — Alcock BP, et al. *CARD 2023…* Nucleic Acids Research. <https://card.mcmaster.ca/>
- **VFDB** — Chen L, et al. *VFDB…* <http://www.mgc.ac.cn/VFs/>
- **PlasmidFinder** — Carattoli A, et al. DTU-CGE. <https://bitbucket.org/genomicepidemiology/plasmidfinder_db>

## Why CARD/VFDB are not redistributed
CARD's Data Usage Agreement prohibits redistribution of the database (and restricts
commercial use), so we cannot ship a CARD-derived index under CC-BY-4.0. VFDB is "free for
academic use" but lacks an explicit open-redistribution license. To stay clearly compliant,
those indexes are **not distributed**: users run the provided scripts to download the source
DB (accepting its terms) and build the Dragon index locally. This keeps the published
collection cleanly CC-BY-4.0.

> Not legal advice — LSHTM Research/Legal should confirm CARD and VFDB handling before release.
