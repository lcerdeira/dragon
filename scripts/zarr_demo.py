#!/usr/bin/env python3
"""Interop demo: open a Dragon Zarr store from Python and recover unitig text.

This script demonstrates the paper claim that a Dragon index exported with
`dragon export-zarr` is readable from any Zarr-aware tool — no Dragon code or
Rust runtime needed. Works against local paths and remote object stores.

Requires:
    pip install 'zarr>=3.0' numcodecs                       # local
    pip install 'zarr>=3.0' numcodecs s3fs                  # also for s3://
    pip install 'zarr>=3.0' numcodecs gcsfs                 # also for gs://

Usage:
    python scripts/zarr_demo.py /path/to/store.zarr
    python scripts/zarr_demo.py s3://dragon-zarr/smoke/idx.zarr
    python scripts/zarr_demo.py gs://bucket/path/idx.zarr
"""

import sys
import time

try:
    import zarr
    import numpy as np
except ImportError:
    print("ERROR: pip install 'zarr>=3.0' numpy", file=sys.stderr)
    sys.exit(1)


def open_store(uri: str):
    """Return a zarr group for either a local path or remote URI."""
    if uri.startswith(("s3://", "gs://", "az://", "abfs://", "http://", "https://")):
        try:
            import fsspec
        except ImportError:
            print(
                "ERROR: remote URIs require fsspec + a backend\n"
                "       s3fs (s3://), gcsfs (gs://), adlfs (az://), aiohttp (http://)",
                file=sys.stderr,
            )
            sys.exit(1)
        proto = uri.split("://", 1)[0]
        path = uri[len(proto) + 3:]
        # Public buckets: anonymous read for reproducibility.
        kwargs = {"anon": True} if proto in ("s3", "gs") else {}
        if proto == "s3":
            kwargs["client_kwargs"] = {"region_name": "eu-west-2"}
        fs = fsspec.filesystem(proto, **kwargs)
        store = zarr.storage.FsspecStore(fs, path=path)
        return zarr.open_group(store, mode="r")
    return zarr.open(uri, mode="r")


def main(store_path: str) -> None:
    t0 = time.perf_counter()
    g = open_store(store_path)
    t_open = time.perf_counter() - t0

    attrs = dict(g.attrs)
    print(f"=== Dragon Zarr store: {store_path} ===")
    print(f"  open latency    : {t_open*1000:.1f} ms")
    print(f"  format_version  : {attrs.get('dragon_zarr_format_version')}")
    print(f"  kmer_size       : {attrs.get('kmer_size')}")
    print(f"  num_unitigs     : {attrs.get('num_unitigs')}")
    print(f"  num_genomes     : {attrs.get('num_genomes')}")
    print(f"  text_len        : {attrs.get('text_len')}")
    print(f"  num_suffixes    : {attrs.get('num_suffixes')}")
    print()

    text = g["text"]
    sa = g["suffix_array"]
    lengths = g["unitig_lengths"][:]
    print("=== Array descriptors ===")
    for name, arr in (("/text", text), ("/suffix_array", sa), ("/unitig_lengths", lengths)):
        if hasattr(arr, "shape"):
            print(f"  {name:<24} shape={arr.shape} dtype={arr.dtype}")
        else:
            print(f"  {name:<24} numpy array shape={arr.shape} dtype={arr.dtype}")
    print()

    # Chunked random read latency — paper-relevant number.
    n = min(128, int(attrs.get("text_len", 0)))
    if n > 0:
        t0 = time.perf_counter()
        chunk = bytes(text[:n])
        t_read = time.perf_counter() - t0
        try:
            preview = chunk.decode("ascii", errors="replace")
        except Exception:
            preview = repr(chunk)
        print(f"=== First {n} bases of /text (chunked decompressing read) ===")
        print(preview)
        print(f"  read latency    : {t_read*1000:.1f} ms")
        print()

    # Walk the first few unitigs using the cumulative-length array.
    cum = np.concatenate([[0], np.cumsum(lengths, dtype=np.int64)])
    preview_n = min(3, len(lengths))
    print(f"=== First {preview_n} unitigs ===")
    for uid in range(preview_n):
        s, e = int(cum[uid]), int(cum[uid] + lengths[uid])
        seq = bytes(text[s:e]).decode("ascii", errors="replace")
        print(f"  unitig {uid}: {seq[:80]}{'…' if len(seq) > 80 else ''}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: zarr_demo.py <store.zarr | s3://... | gs://...>", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1])
