#!/usr/bin/env python3
"""
download_genomes.py -- Download real bacterial genomes from NCBI for Dragon benchmarking.

Tier 1: 50 complete E. coli genomes (within-species alignment benchmark)
Tier 2: 30 genomes from 6 species, 5 each (cross-species alignment benchmark)

Usage:
    python3 benchmark/download_genomes.py [--tier1-only | --tier2-only] [--output-dir DIR]

Downloads are re-runnable: already-downloaded files (verified by checksum) are skipped.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import os
import shutil
import subprocess
import sys
import textwrap
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NCBI_DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# Taxonomy IDs for target species
SPECIES = {
    "Escherichia_coli":          562,
    "Klebsiella_pneumoniae":     573,
    "Salmonella_enterica":       28901,
    "Pseudomonas_aeruginosa":    287,
    "Staphylococcus_aureus":     1280,
    "Streptococcus_pneumoniae":  1313,
}

# How many genomes per species in each tier
TIER1_SPECIES = "Escherichia_coli"
TIER1_COUNT = 50
TIER2_PER_SPECIES = 5

REQUEST_HEADERS = {
    "Accept": "application/json",
    "User-Agent": "Dragon-benchmark/1.0 (https://github.com/Dragon; mailto:benchmark@example.com)",
}

# Retry / rate-limit settings
MAX_RETRIES = 4
RETRY_BACKOFF = 2.0        # seconds, doubles each retry
RATE_LIMIT_SLEEP = 0.35    # minimum pause between NCBI requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("download_genomes")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class GenomeRecord:
    accession: str
    species: str
    ftp_path: Optional[str] = None
    filename: Optional[str] = None
    total_length: int = 0
    num_contigs: int = 0
    md5: str = ""


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _urlopen_with_retry(url: str, retries: int = MAX_RETRIES) -> bytes:
    """Fetch a URL with retries and exponential backoff."""
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=REQUEST_HEADERS)
            with urllib.request.urlopen(req, timeout=120) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            last_err = exc
            wait = RETRY_BACKOFF * (2 ** attempt)
            log.warning("  Attempt %d/%d failed for %s: %s -- retrying in %.1fs",
                        attempt + 1, retries, url, exc, wait)
            time.sleep(wait)
    raise RuntimeError(f"Failed to fetch {url} after {retries} attempts: {last_err}")


def md5_file(path: Path) -> str:
    """Compute MD5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fasta_stats(path: Path) -> tuple[int, int]:
    """Return (total_bases, num_contigs) for a FASTA file."""
    total = 0
    contigs = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs += 1
            else:
                total += len(line.strip())
    return total, contigs


def sanitize_species(name: str) -> str:
    """Convert species name to filesystem-safe form."""
    return name.replace(" ", "_").replace("/", "_")


# ---------------------------------------------------------------------------
# Strategy 1: NCBI Datasets CLI (fastest if installed)
# ---------------------------------------------------------------------------

def _datasets_cli_available() -> bool:
    """Check whether the NCBI `datasets` CLI is on PATH."""
    return shutil.which("datasets") is not None


def download_via_cli(
    taxid: int,
    species: str,
    limit: int,
    dest_dir: Path,
) -> list[GenomeRecord]:
    """
    Use `datasets download genome taxon ...` to fetch genomes.
    Returns list of GenomeRecord with filenames populated.
    """
    log.info("Strategy: NCBI Datasets CLI for %s (taxid=%d, limit=%d)", species, taxid, limit)
    tmp_zip = dest_dir / f"_ncbi_dataset_{taxid}.zip"

    cmd = [
        "datasets", "download", "genome", "taxon", str(taxid),
        "--assembly-level", "complete",
        "--assembly-source", "refseq",
        "--include", "genome",
        "--limit", str(limit),
        "--filename", str(tmp_zip),
    ]
    log.info("  Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        log.warning("  datasets CLI failed: %s", result.stderr.strip())
        raise RuntimeError(f"datasets CLI returned {result.returncode}")

    records = _extract_zip(tmp_zip, species, dest_dir)
    tmp_zip.unlink(missing_ok=True)
    return records


def _extract_zip(zip_path: Path, species: str, dest_dir: Path) -> list[GenomeRecord]:
    """Extract FASTA files from an NCBI datasets zip archive."""
    records = []
    safe_species = sanitize_species(species)

    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            # Typical path: ncbi_dataset/data/GCF_xxx/GCF_xxx_genomic.fna
            if not info.filename.endswith(".fna") and not info.filename.endswith(".fna.gz"):
                continue
            parts = info.filename.split("/")
            accession = None
            for p in parts:
                if p.startswith("GCF_") or p.startswith("GCA_"):
                    accession = p
                    break
            if accession is None:
                continue

            out_name = f"{safe_species}_{accession}.fa"
            out_path = dest_dir / out_name

            if out_path.exists():
                log.info("  Already exists, skipping: %s", out_name)
            else:
                data = zf.read(info.filename)
                # Decompress if gzipped inside the zip
                if info.filename.endswith(".gz"):
                    data = gzip.decompress(data)
                out_path.write_bytes(data)
                log.info("  Extracted: %s", out_name)

            total_len, n_contigs = fasta_stats(out_path)
            records.append(GenomeRecord(
                accession=accession,
                species=species,
                filename=out_name,
                total_length=total_len,
                num_contigs=n_contigs,
                md5=md5_file(out_path),
            ))
    return records


# ---------------------------------------------------------------------------
# Strategy 2: NCBI Datasets REST API v2 -> individual FTP downloads
# ---------------------------------------------------------------------------

def _fetch_accessions_rest(taxid: int, limit: int) -> list[dict]:
    """
    Query the NCBI Datasets v2 REST API for complete RefSeq genomes,
    returning a list of dicts with 'accession' and 'ftp_path' keys.
    """
    accessions: list[dict] = []
    page_token = None

    while len(accessions) < limit:
        page_size = min(limit - len(accessions), 100)
        url = (
            f"{NCBI_DATASETS_API}/genome/taxon/{taxid}/dataset_report"
            f"?filters.assembly_level=complete"
            f"&filters.assembly_source=refseq"
            f"&page_size={page_size}"
        )
        if page_token:
            url += f"&page_token={page_token}"

        raw = _urlopen_with_retry(url)
        data = json.loads(raw)
        time.sleep(RATE_LIMIT_SLEEP)

        reports = data.get("reports", [])
        if not reports:
            break

        for rpt in reports:
            acc = rpt.get("accession", "")
            # Build FTP path from the report if available
            ftp = ""
            assembly_info = rpt.get("assembly_info", {})
            refseq_cat = rpt.get("source_database", "")
            # Try to get the FTP path from various locations in the response
            paired_assembly = rpt.get("paired_accession", "")

            # Construct FTP path from accession
            if acc.startswith("GCF_") or acc.startswith("GCA_"):
                ftp = _accession_to_ftp_url(acc)

            accessions.append({"accession": acc, "ftp_path": ftp})
            if len(accessions) >= limit:
                break

        page_token = data.get("next_page_token")
        if not page_token:
            break

    return accessions[:limit]


def _accession_to_ftp_url(accession: str) -> str:
    """
    Convert a GenBank/RefSeq accession to an NCBI FTP directory URL.
    E.g. GCF_000005845.2 -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_*
    We return the directory URL; the caller will list it to find the FASTA.
    """
    prefix = accession[:3]    # GCF or GCA
    digits = accession[4:].split(".")[0]  # numeric part before version
    # Pad to at least 9 digits
    digits = digits.zfill(9)
    p1, p2, p3 = digits[0:3], digits[3:6], digits[6:9]
    return f"{NCBI_FTP_BASE}/{prefix}/{p1}/{p2}/{p3}/"


def _find_fasta_url_from_ftp_dir(dir_url: str, accession: str) -> Optional[str]:
    """
    Given an FTP directory URL, find the *_genomic.fna.gz FASTA URL.
    We list the directory via HTTPS and parse the HTML.
    """
    try:
        html = _urlopen_with_retry(dir_url).decode("utf-8", errors="replace")
    except Exception:
        return None

    # Look for directories that start with the accession
    import re
    # Find subdirectory matching the accession (e.g. GCF_000005845.2_ASM584v2)
    pattern = re.compile(rf'href="({re.escape(accession)}[^"/]*)/?"', re.IGNORECASE)
    matches = pattern.findall(html)
    if not matches:
        # Try a looser match
        acc_base = accession.split(".")[0]
        pattern = re.compile(rf'href="({re.escape(acc_base)}[^"/]*)/?"', re.IGNORECASE)
        matches = pattern.findall(html)

    if not matches:
        return None

    asm_dir = matches[0]
    asm_url = f"{dir_url}{asm_dir}/"

    # Now look for *_genomic.fna.gz in that directory
    try:
        html2 = _urlopen_with_retry(asm_url).decode("utf-8", errors="replace")
    except Exception:
        return None

    fna_pattern = re.compile(rf'href="([^"]*_genomic\.fna\.gz)"')
    fna_matches = fna_pattern.findall(html2)

    # Filter out _from_genomic and _rna_from_genomic variants
    fna_matches = [
        m for m in fna_matches
        if "_from_genomic" not in m
        and "_rna_from_genomic" not in m
        and "_cds_from_genomic" not in m
    ]

    if not fna_matches:
        return None

    return f"{asm_url}{fna_matches[0]}"


def download_via_rest_api(
    taxid: int,
    species: str,
    limit: int,
    dest_dir: Path,
) -> list[GenomeRecord]:
    """
    Use the NCBI Datasets v2 REST API to get accession list, then download
    each genome's FASTA individually from the FTP mirror.
    """
    log.info("Strategy: NCBI Datasets REST API for %s (taxid=%d, limit=%d)", species, taxid, limit)

    accession_info = _fetch_accessions_rest(taxid, limit)
    if not accession_info:
        raise RuntimeError(f"REST API returned no accessions for taxid {taxid}")

    log.info("  Found %d accession(s) for %s", len(accession_info), species)
    records = []
    safe_species = sanitize_species(species)

    for i, info in enumerate(accession_info, 1):
        acc = info["accession"]
        ftp_dir = info["ftp_path"]
        out_name = f"{safe_species}_{acc}.fa"
        out_path = dest_dir / out_name

        if out_path.exists():
            log.info("  [%d/%d] Already exists: %s", i, len(accession_info), out_name)
            total_len, n_contigs = fasta_stats(out_path)
            records.append(GenomeRecord(
                accession=acc, species=species, filename=out_name,
                total_length=total_len, num_contigs=n_contigs, md5=md5_file(out_path),
            ))
            continue

        log.info("  [%d/%d] Downloading %s ...", i, len(accession_info), acc)

        # Find the actual FASTA URL
        fasta_url = _find_fasta_url_from_ftp_dir(ftp_dir, acc) if ftp_dir else None

        if fasta_url is None:
            # Try the NCBI Datasets download endpoint as individual fallback
            fasta_url = _try_datasets_download_url(acc)

        if fasta_url is None:
            log.warning("  Could not resolve FASTA URL for %s, skipping", acc)
            continue

        try:
            raw = _urlopen_with_retry(fasta_url)
            # Decompress if gzipped
            if fasta_url.endswith(".gz"):
                raw = gzip.decompress(raw)
            out_path.write_bytes(raw)
            log.info("  Saved: %s (%.1f MB)", out_name, len(raw) / 1e6)
        except Exception as exc:
            log.warning("  Failed to download %s: %s", acc, exc)
            continue

        time.sleep(RATE_LIMIT_SLEEP)

        total_len, n_contigs = fasta_stats(out_path)
        records.append(GenomeRecord(
            accession=acc, species=species, filename=out_name,
            total_length=total_len, num_contigs=n_contigs, md5=md5_file(out_path),
        ))

    return records


def _try_datasets_download_url(accession: str) -> Optional[str]:
    """
    Attempt to build a direct download URL via the NCBI Datasets v2 download
    endpoint for a single genome by accession.
    Returns a URL string on success or None.
    """
    url = (
        f"{NCBI_DATASETS_API}/genome/accession/{accession}/download"
        f"?include_annotation_type=GENOME_FASTA"
    )
    try:
        # We don't actually download here -- just return the URL for the
        # zip-based approach below
        return url
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Strategy 3: NCBI Datasets API zip download per accession (last resort)
# ---------------------------------------------------------------------------

def download_via_api_zip(
    taxid: int,
    species: str,
    limit: int,
    dest_dir: Path,
) -> list[GenomeRecord]:
    """
    Fetch accession list from the REST API, then download each genome as an
    NCBI Datasets zip via the download endpoint and extract the FASTA.
    """
    log.info("Strategy: NCBI Datasets API (zip per accession) for %s (taxid=%d, limit=%d)",
             species, taxid, limit)

    accession_info = _fetch_accessions_rest(taxid, limit)
    if not accession_info:
        raise RuntimeError(f"API returned no accessions for taxid {taxid}")

    records = []
    safe_species = sanitize_species(species)

    for i, info in enumerate(accession_info, 1):
        acc = info["accession"]
        out_name = f"{safe_species}_{acc}.fa"
        out_path = dest_dir / out_name

        if out_path.exists():
            log.info("  [%d/%d] Already exists: %s", i, len(accession_info), out_name)
            total_len, n_contigs = fasta_stats(out_path)
            records.append(GenomeRecord(
                accession=acc, species=species, filename=out_name,
                total_length=total_len, num_contigs=n_contigs, md5=md5_file(out_path),
            ))
            continue

        log.info("  [%d/%d] Downloading zip for %s ...", i, len(accession_info), acc)
        zip_url = (
            f"{NCBI_DATASETS_API}/genome/accession/{acc}/download"
            f"?include_annotation_type=GENOME_FASTA"
        )

        try:
            raw_zip = _urlopen_with_retry(zip_url)
        except Exception as exc:
            log.warning("  Failed to download zip for %s: %s", acc, exc)
            continue

        time.sleep(RATE_LIMIT_SLEEP)

        # Extract FASTA from the zip
        try:
            with zipfile.ZipFile(io.BytesIO(raw_zip)) as zf:
                fasta_data = None
                for zinfo in zf.infolist():
                    name = zinfo.filename
                    if name.endswith(".fna") or name.endswith(".fna.gz") or name.endswith(".fa"):
                        data = zf.read(name)
                        if name.endswith(".gz"):
                            data = gzip.decompress(data)
                        fasta_data = data
                        break

                if fasta_data is None:
                    log.warning("  No FASTA found in zip for %s", acc)
                    continue

                out_path.write_bytes(fasta_data)
                log.info("  Saved: %s (%.1f MB)", out_name, len(fasta_data) / 1e6)
        except (zipfile.BadZipFile, Exception) as exc:
            log.warning("  Failed to process zip for %s: %s", acc, exc)
            continue

        total_len, n_contigs = fasta_stats(out_path)
        records.append(GenomeRecord(
            accession=acc, species=species, filename=out_name,
            total_length=total_len, num_contigs=n_contigs, md5=md5_file(out_path),
        ))

    return records


# ---------------------------------------------------------------------------
# Orchestration: try strategies in order
# ---------------------------------------------------------------------------

def download_genomes(
    taxid: int,
    species: str,
    limit: int,
    dest_dir: Path,
) -> list[GenomeRecord]:
    """
    Download `limit` complete RefSeq genomes for a given species/taxid.
    Tries three strategies in order:
      1. NCBI Datasets CLI
      2. NCBI Datasets REST API + FTP download
      3. NCBI Datasets REST API + per-accession zip download
    Returns a list of GenomeRecord for successfully downloaded genomes.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)

    strategies = []

    if _datasets_cli_available():
        strategies.append(("CLI", download_via_cli))
    strategies.append(("REST+FTP", download_via_rest_api))
    strategies.append(("REST+ZIP", download_via_api_zip))

    for name, fn in strategies:
        try:
            records = fn(taxid, species, limit, dest_dir)
            if records:
                log.info("  %s strategy yielded %d genome(s) for %s", name, len(records), species)
                return records
            else:
                log.warning("  %s strategy returned 0 genomes, trying next ...", name)
        except Exception as exc:
            log.warning("  %s strategy failed for %s: %s", name, species, exc)

    log.error("All download strategies failed for %s (taxid=%d)", species, taxid)
    return []


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def write_manifest(records: list[GenomeRecord], manifest_path: Path) -> None:
    """Write a TSV manifest of all downloaded genomes."""
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with open(manifest_path, "w") as fh:
        fh.write("filename\taccession\tspecies\ttotal_length\tnum_contigs\tmd5\n")
        for r in sorted(records, key=lambda x: (x.species, x.accession)):
            fh.write(f"{r.filename}\t{r.accession}\t{r.species}\t"
                     f"{r.total_length}\t{r.num_contigs}\t{r.md5}\n")

    log.info("Manifest written: %s (%d entries)", manifest_path, len(records))


# ---------------------------------------------------------------------------
# Checksum-based skip logic
# ---------------------------------------------------------------------------

def load_existing_manifest(manifest_path: Path) -> dict[str, str]:
    """Load an existing manifest, returning {filename: md5}."""
    if not manifest_path.exists():
        return {}
    result = {}
    with open(manifest_path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                result[parts[0]] = parts[5]
    return result


def verify_existing_files(dest_dir: Path, manifest: dict[str, str]) -> set[str]:
    """Return set of filenames that already exist and match their checksum."""
    verified = set()
    for fname, expected_md5 in manifest.items():
        fpath = dest_dir / fname
        if fpath.exists():
            actual_md5 = md5_file(fpath)
            if actual_md5 == expected_md5:
                verified.add(fname)
            else:
                log.warning("  Checksum mismatch for %s, will re-download", fname)
                fpath.unlink()
    return verified


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download real bacterial genomes from NCBI for Dragon benchmarking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python3 benchmark/download_genomes.py
              python3 benchmark/download_genomes.py --tier1-only
              python3 benchmark/download_genomes.py --tier2-only
              python3 benchmark/download_genomes.py --output-dir /tmp/genomes
        """),
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).resolve().parent / "data" / "real_genomes",
        help="Root output directory (default: benchmark/data/real_genomes/)",
    )
    parser.add_argument("--tier1-only", action="store_true", help="Download only Tier 1 genomes")
    parser.add_argument("--tier2-only", action="store_true", help="Download only Tier 2 genomes")
    parser.add_argument("--tier1-count", type=int, default=TIER1_COUNT,
                        help=f"Number of E. coli genomes for Tier 1 (default: {TIER1_COUNT})")
    parser.add_argument("--tier2-per-species", type=int, default=TIER2_PER_SPECIES,
                        help=f"Genomes per species for Tier 2 (default: {TIER2_PER_SPECIES})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Only query accessions, don't download")
    args = parser.parse_args()

    root = args.output_dir.resolve()
    do_tier1 = not args.tier2_only
    do_tier2 = not args.tier1_only

    all_records: list[GenomeRecord] = []

    # ------------------------------------------------------------------
    # Tier 1: 50 complete E. coli genomes
    # ------------------------------------------------------------------
    if do_tier1:
        tier1_dir = root / "tier1"
        tier1_manifest = tier1_dir / "manifest.tsv"
        log.info("=" * 60)
        log.info("TIER 1: %d complete %s genomes", args.tier1_count, TIER1_SPECIES)
        log.info("  Output: %s", tier1_dir)
        log.info("=" * 60)

        existing = load_existing_manifest(tier1_manifest)
        verified = verify_existing_files(tier1_dir, existing)
        if len(verified) >= args.tier1_count:
            log.info("  All %d Tier 1 genomes already present and verified.", len(verified))
            # Reload records from manifest
            with open(tier1_manifest) as fh:
                fh.readline()  # skip header
                for line in fh:
                    p = line.strip().split("\t")
                    if len(p) >= 6:
                        all_records.append(GenomeRecord(
                            accession=p[1], species=p[2], filename=p[0],
                            total_length=int(p[3]), num_contigs=int(p[4]), md5=p[5],
                        ))
        else:
            if not args.dry_run:
                records = download_genomes(
                    SPECIES[TIER1_SPECIES], TIER1_SPECIES, args.tier1_count, tier1_dir,
                )
                all_records.extend(records)
                write_manifest(records, tier1_manifest)
                log.info("Tier 1 complete: %d / %d genomes downloaded",
                         len(records), args.tier1_count)
            else:
                log.info("  [DRY RUN] Would download %d %s genomes",
                         args.tier1_count, TIER1_SPECIES)

    # ------------------------------------------------------------------
    # Tier 2: 5 genomes each from 6 species
    # ------------------------------------------------------------------
    if do_tier2:
        tier2_dir = root / "tier2"
        tier2_manifest = tier2_dir / "manifest.tsv"
        log.info("=" * 60)
        log.info("TIER 2: %d genomes each from %d species", args.tier2_per_species, len(SPECIES))
        log.info("  Output: %s", tier2_dir)
        log.info("=" * 60)

        existing = load_existing_manifest(tier2_manifest)
        verified = verify_existing_files(tier2_dir, existing)
        expected_total = args.tier2_per_species * len(SPECIES)

        if len(verified) >= expected_total:
            log.info("  All %d Tier 2 genomes already present and verified.", len(verified))
            with open(tier2_manifest) as fh:
                fh.readline()
                for line in fh:
                    p = line.strip().split("\t")
                    if len(p) >= 6:
                        all_records.append(GenomeRecord(
                            accession=p[1], species=p[2], filename=p[0],
                            total_length=int(p[3]), num_contigs=int(p[4]), md5=p[5],
                        ))
        else:
            tier2_records: list[GenomeRecord] = []
            for species_name, taxid in SPECIES.items():
                log.info("-" * 40)
                log.info("Species: %s (taxid=%d)", species_name, taxid)

                # Count how many of this species are already verified
                already_have = sum(
                    1 for f in verified
                    if f.startswith(sanitize_species(species_name) + "_")
                )
                need = args.tier2_per_species - already_have
                if need <= 0:
                    log.info("  Already have %d genome(s), skipping", already_have)
                    # Add existing records from manifest
                    for fname, md5_val in existing.items():
                        if fname.startswith(sanitize_species(species_name) + "_"):
                            fpath = tier2_dir / fname
                            if fpath.exists():
                                tl, nc = fasta_stats(fpath)
                                acc = fname.replace(sanitize_species(species_name) + "_", "").replace(".fa", "")
                                tier2_records.append(GenomeRecord(
                                    accession=acc, species=species_name, filename=fname,
                                    total_length=tl, num_contigs=nc, md5=md5_val,
                                ))
                    continue

                if args.dry_run:
                    log.info("  [DRY RUN] Would download %d genomes", need)
                    continue

                records = download_genomes(
                    taxid, species_name, args.tier2_per_species, tier2_dir,
                )
                tier2_records.extend(records)
                log.info("  Downloaded %d / %d for %s",
                         len(records), args.tier2_per_species, species_name)

            all_records.extend(tier2_records)
            if tier2_records and not args.dry_run:
                write_manifest(tier2_records, tier2_manifest)
            log.info("Tier 2 complete: %d / %d genomes downloaded",
                     len(tier2_records), expected_total)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log.info("=" * 60)
    log.info("SUMMARY")
    log.info("  Total genomes: %d", len(all_records))
    species_counts = {}
    for r in all_records:
        species_counts[r.species] = species_counts.get(r.species, 0) + 1
    for sp, cnt in sorted(species_counts.items()):
        log.info("    %s: %d", sp, cnt)
    total_bp = sum(r.total_length for r in all_records)
    log.info("  Total base pairs: %s", f"{total_bp:,}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
