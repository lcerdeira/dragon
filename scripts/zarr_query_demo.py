#!/usr/bin/env python3
"""Time a single-pattern backward search against a Dragon Zarr store.

Reads a query FASTA, opens a Zarr-backed index (local path or s3:// URI),
and times naive binary-search lookup over the suffix array. Designed for
small reproducible queries (≤ a few thousand bp) that fit in one chunk.

Usage:
    python scripts/zarr_query_demo.py <store> <query.fa>

    # public S3 (anonymous read, no AWS credentials needed):
    python scripts/zarr_query_demo.py s3://dragon-zarr/saureus/b1 ermC.fa

    # local-disk Zarr exported by `dragon export-zarr`:
    python scripts/zarr_query_demo.py ~/dragon_work/saureus_index_batch1.zarr ermC.fa

Reports for each query record: open latency, search latency, hit count,
and (if any) the unitig IDs covering the first few hits.
"""

import sys
import time
from pathlib import Path

try:
    import zarr
except ImportError:
    print("ERROR: pip install 'zarr>=3.0' s3fs numcodecs", file=sys.stderr)
    sys.exit(1)


def open_store(uri: str):
    """Open a Zarr group from a local path or remote URI."""
    if uri.startswith(("s3://", "gs://", "az://", "abfs://", "http://", "https://")):
        import fsspec
        proto = uri.split("://", 1)[0]
        path = uri[len(proto) + 3:]
        kwargs = {"anon": True} if proto in ("s3", "gs") else {}
        if proto == "s3":
            kwargs["client_kwargs"] = {"region_name": "eu-west-2"}
        fs = fsspec.filesystem(proto, **kwargs)
        return zarr.open_group(zarr.storage.FsspecStore(fs, path=path), mode="r")
    return zarr.open(uri, mode="r")


def parse_fasta(path: str):
    """Yield (name, sequence_bytes) for each record in a FASTA file."""
    name, buf = None, []
    with open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, b"".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.encode("ascii", errors="ignore"))
    if name is not None:
        yield name, b"".join(buf).upper()


def search_pattern(text, sa, pattern: bytes):
    """Return all text positions where `pattern` occurs (sorted ascending).

    Pure binary search over the suffix array. Each step touches one SA
    entry and one len(pattern)-byte slice of /text — both Zarr-chunked,
    so cold cache means a couple of HTTP range requests per step over S3.
    """
    n = sa.shape[0]
    plen = len(pattern)

    def lower(p):
        lo, hi = 0, n
        while lo < hi:
            mid = (lo + hi) // 2
            ss = int(sa[mid])
            suf = bytes(text[ss:ss + plen])
            if suf < p:
                lo = mid + 1
            else:
                hi = mid
        return lo

    def upper(p):
        lo, hi = 0, n
        while lo < hi:
            mid = (lo + hi) // 2
            ss = int(sa[mid])
            suf = bytes(text[ss:ss + plen])
            if suf <= p:
                lo = mid + 1
            else:
                hi = mid
        return lo

    a = lower(pattern)
    b = upper(pattern)
    if b <= a:
        return []
    return sorted(int(p) for p in sa[a:b])


def main(store_uri: str, query_path: str) -> None:
    print(f"=== store: {store_uri} ===")
    print(f"=== query: {query_path} ===\n")

    t0 = time.perf_counter()
    g = open_store(store_uri)
    text = g["text"]
    sa = g["suffix_array"]
    lengths = g["unitig_lengths"][:]
    attrs = dict(g.attrs)
    open_ms = (time.perf_counter() - t0) * 1000

    print(f"index opened in {open_ms:6.0f} ms")
    print(
        f"  num_genomes={attrs.get('num_genomes')}  "
        f"num_unitigs={attrs.get('num_unitigs')}  "
        f"text_len={attrs.get('text_len')}  "
        f"k={attrs.get('kmer_size')}\n"
    )

    # Cumulative unitig boundaries for position -> unitig mapping.
    import numpy as np
    cum = np.concatenate([[0], np.cumsum(lengths, dtype=np.int64)])

    for name, seq in parse_fasta(query_path):
        print(f"--- query: {name}  ({len(seq)} bp) ---")
        if len(seq) == 0:
            print("  empty sequence, skipping")
            continue

        t0 = time.perf_counter()
        positions = search_pattern(text, sa, seq)
        search_ms = (time.perf_counter() - t0) * 1000

        print(f"  {len(positions)} occurrence(s) found in {search_ms:.0f} ms")
        for pos in positions[:5]:
            uid = int(np.searchsorted(cum, pos, side="right") - 1)
            offset = pos - int(cum[uid])
            print(f"    unitig {uid:>10} offset {offset:>8}  text_pos {pos:>10}")
        if len(positions) > 5:
            print(f"    ... and {len(positions) - 5} more")
        print()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
