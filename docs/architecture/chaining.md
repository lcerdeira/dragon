# Graph-aware colinear chaining

## The chaining problem

After seed finding, Dragon has a set of seed hits: `{(query_pos, unitig_id, offset, match_len)}`. The goal is to find the **best colinear subset** — a chain of seeds that are consistent in both query and reference order, representing a contiguous alignment.

## Why graph-aware?

Traditional tools (Minimap2, LexicMap) chain seeds in **linear coordinate space** — assuming the reference is a flat string. Dragon chains seeds along the **genome's path through the de Bruijn graph**, which:

1. Naturally handles structural variation (alternative paths in the graph)
2. Respects the actual genome structure
3. Avoids spurious chains across unrelated genomic regions

## Algorithm

### Step 1: Map seeds to genome coordinates

For each candidate genome, Dragon loads its path (ordered sequence of unitig IDs) and maps each seed hit to a genome-level coordinate:

```
seed at (unitig_42, offset_100)
  + genome path: [..., unitig_42 at genome_pos_50000, ...]
  = anchor at (query_pos, genome_pos_50100, match_len)
```

### Step 2: Sort anchors by reference position

Anchors are sorted by their genome coordinate, separating forward and reverse-complement hits.

### Step 3: Fenwick tree DP

The colinear chaining DP finds the highest-scoring subset of anchors where each anchor starts after the previous one ends in both query and reference:

```
For each anchor a[i] (left to right by ref position):
    score[i] = match_len[i]  +  max over valid predecessors j of score[j]
    where: a[j].query_end <= a[i].query_start
           a[j].ref_end   <= a[i].ref_start
```

Using a Fenwick tree (binary indexed tree) indexed by query position, this runs in **O(h log h)** time where h is the number of anchors.

### Step 4: Gap-sensitive scoring

The chain score includes a gap penalty proportional to the difference between reference and query gap lengths:

```
gap_cost(a[i], a[j]) = alpha * |ref_gap - query_gap|
```

This penalises indels and structural variations proportional to their size.

## Complexity

| Step | Time | Space |
|------|------|-------|
| Map to genome coordinates | O(h) | O(h) |
| Sort | O(h log h) | O(h) |
| Fenwick tree DP | O(h log h) | O(h) |
| Traceback | O(chain_len) | O(chain_len) |
| **Total** | **O(h log h)** | **O(h)** |

For a typical gene query with h ~ 100 anchors, chaining takes <1 ms per genome.
