# Data structures

## 2-bit DNA encoding

DNA bases are packed 32 per `u64` word:

| Base | Encoding |
|------|----------|
| A | `00` |
| C | `01` |
| G | `10` |
| T | `11` |

Operations:
- **Encode**: 10,000 bases in ~2 microseconds
- **Decode**: 10,000 bases in ~3 microseconds
- **Reverse complement**: bitwise NOT + reverse pairs
- **K-mer extraction**: shift and mask operations

Ambiguous bases (N, R, Y, etc.) are encoded as A (`00`).

## Roaring Bitmaps

Used for the colour index (mapping unitigs to genome sets). Roaring Bitmaps partition the integer universe into blocks of 2^16, using the optimal representation per block:

- **Array container**: for sparse blocks (<4,096 integers)
- **Bitmap container**: for dense blocks (>4,096 integers)
- **Run container**: for blocks with long consecutive runs

This is ideal for genome ID sets, which tend to be clustered by species.

## Elias-Fano cumulative length index

Maps FM-index text positions to unitig IDs via binary search on a sorted array of cumulative start positions. Current implementation uses a simple `Vec<u64>` with binary search; production version can upgrade to a true Elias-Fano structure for O(1) predecessor queries.

## Fenwick tree (Binary Indexed Tree)

Used in colinear chaining for prefix maximum queries:

- **Update**: set position i to max(current, value) in O(log n)
- **Query**: find maximum in prefix [0, i] in O(log n)
- **Space**: O(n)

Two variants:
- `FenwickMax`: prefix maximum (used in chaining)
- `FenwickSum`: prefix sum (used in vote counting)

## Variable-length integers (varint)

LEB128 encoding for genome path compression:

| Value range | Bytes |
|------------|-------|
| 0-127 | 1 |
| 128-16,383 | 2 |
| 16,384-2,097,151 | 3 |
| 2,097,152-268,435,455 | 4 |

Also supports:
- **Zigzag encoding**: efficient encoding of signed integers
- **Delta encoding**: for sorted sequences (store differences)
- **Slice encoding/decoding**: batch operations
