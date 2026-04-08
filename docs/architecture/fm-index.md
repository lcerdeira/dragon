# Run-length FM-index

## Background

The FM-index is a compressed full-text index based on the Burrows-Wheeler Transform (BWT). Given a text T of length n, it supports:
- **Count**: how many times pattern P occurs in T, in O(|P|) time
- **Locate**: the positions of all occurrences, in O(|P| + occ) time
- **Space**: O(n) for standard FM-index, O(r) for run-length variant

## Why run-length?

The **run-length FM-index** (r-index) stores the BWT in run-length encoded form. The number of BWT runs **r** is the key parameter:

- For random text: r ~ n (no compression)
- For repetitive text: r << n (massive compression)

Prokaryotic genome collections are **highly repetitive**: related species share >95% of their sequence. When their unitig sequences are concatenated, the BWT contains long runs of identical characters.

| Database | Text length n | BWT runs r | r/n | Compression |
|----------|--------------|-----------|-----|-------------|
| 500 E. coli | ~2.5 Gbp | ~50M | 0.02 | 50x |
| 85K GTDB | ~250 Gbp | ~2.5B | 0.01 | 100x |
| 2.34M GenBank | ~5 Gbp unitigs | ~250M | 0.05 | 20x |

## Construction

1. **Concatenate** all unitig sequences with separator characters (`$`)
2. **Build** the suffix array using SA-IS (linear time)
3. **Derive** the BWT from the suffix array
4. **Run-length encode** the BWT
5. **Sample** the suffix array at run boundaries

## Backward search

To search for pattern P[1..m]:

1. Initialise suffix array interval [lo, hi) to the full range
2. Process P from right to left: for each character P[i], update [lo, hi) using the LF-mapping
3. After m steps, `hi - lo` = number of occurrences

## Variable-length seed matching

Unlike fixed k-mer lookup, Dragon extends each backward search beyond k characters as long as the SA interval remains non-empty:

```
Query:  ACGTACGTACGTNNNN...
Search: ACGTACGTACGT       (12 chars, SA interval = [50, 55))
        ACGTACGTACGTN      (13 chars, SA interval = empty -> stop)
Result: 12-mer match with 5 occurrences
```

This produces **super-maximal exact matches (SMEMs)** of variable length, which are more informative than fixed-length k-mers.

## Position-to-unitig mapping

FM-index positions are in the concatenated text. Dragon maps them back to unitig IDs using a **cumulative length index** (Elias-Fano encoded monotone sequence):

```
Position 12345 in text
  -> binary search on cumulative lengths
  -> unitig_id = 42, offset = 789
```
