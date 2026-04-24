#!/usr/bin/env python3
"""Interop demo: open a Dragon Zarr store from Python and recover unitig text.

This script demonstrates the paper claim that a Dragon index exported with
`dragon export-zarr` is readable from any Zarr-aware tool — no Dragon code or
Rust runtime needed.

Requires:
    pip install zarr>=3.0 numcodecs

Usage:
    dragon export-zarr -i /path/to/index -o /path/to/store.zarr
    python scripts/zarr_demo.py /path/to/store.zarr
"""

import sys
from pathlib import Path

try:
    import zarr
    import numpy as np
except ImportError:
    print("ERROR: install zarr and numpy first: pip install 'zarr>=3.0' numpy", file=sys.stderr)
    sys.exit(1)


def main(store_path: str) -> None:
    store = zarr.open(store_path, mode="r")

    attrs = dict(store.attrs)
    print(f"=== Dragon Zarr store: {store_path} ===")
    print(f"  format_version  : {attrs.get('dragon_zarr_format_version')}")
    print(f"  kmer_size       : {attrs.get('kmer_size')}")
    print(f"  num_unitigs     : {attrs.get('num_unitigs')}")
    print(f"  num_genomes     : {attrs.get('num_genomes')}")
    print(f"  text_len        : {attrs.get('text_len')}")
    print(f"  num_suffixes    : {attrs.get('num_suffixes')}")
    print()

    text = store["text"]
    sa = store["suffix_array"]
    lengths = store["unitig_lengths"][:]
    print("=== Array descriptors ===")
    for name, arr in (("/text", text), ("/suffix_array", sa), ("/unitig_lengths", lengths)):
        if hasattr(arr, "shape"):
            print(f"  {name:<24} shape={arr.shape} dtype={arr.dtype} chunks={getattr(arr, 'chunks', '-')}")
        else:
            print(f"  {name:<24} numpy array shape={arr.shape} dtype={arr.dtype}")
    print()

    # Demonstrate chunked random access: pull the first 128 bases of the
    # concatenated text without decompressing the whole array.
    n = min(128, int(attrs.get("text_len", 0)))
    if n > 0:
        chunk = bytes(text[:n])
        try:
            preview = chunk.decode("ascii", errors="replace")
        except Exception:
            preview = repr(chunk)
        print(f"=== First {n} bases of /text (chunked read) ===")
        print(preview)
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
        print("Usage: zarr_demo.py <store.zarr>", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1])
