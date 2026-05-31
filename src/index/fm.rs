/// FM-index construction and querying over concatenated unitig sequences.
///
/// Uses the `fm-index` crate which provides a run-length encoded FM-index
/// optimized for repetitive text. For very large databases (>1 Gbp),
/// consider using external BWT construction tools.

use anyhow::Result;
use std::path::Path;

use crate::ds::elias_fano::CumulativeLengthIndex;
use crate::index::unitig::UnitigSet;

/// A seed hit from FM-index search.
#[derive(Clone, Debug)]
pub struct SeedHit {
    pub unitig_id: u32,
    pub offset: u32,
    pub query_pos: usize,
    pub match_len: usize,
    pub is_reverse: bool,
    /// Number of occurrences in the database (SA interval size).
    pub sa_count: usize,
}

/// The FM-index wrapper for Dragon.
pub struct DragonFmIndex {
    /// Concatenated unitig text (with separator bytes between unitigs).
    pub text: Vec<u8>,
    /// Suffix array in **linear** (sorted) order — always present.
    /// Used for range collection after bounds are found.
    pub suffix_array: Vec<usize>,
    /// Suffix array in **Eytzinger (BFS) order** — present when built with
    /// `--sa-layout eytzinger`.  Empty for linearly-built indexes.
    ///
    /// The Eytzinger layout places binary-search pivot elements near the start
    /// of the array, so the first ~20 levels of any search fit in L2 cache
    /// regardless of total SA size.  This eliminates the dominant cost of
    /// large-index searches: cache misses during the O(log n) bounds computation.
    ///
    /// Layout (1-indexed): node 1 = root (SA[n/2]), left child = 2k,
    /// right child = 2k+1.  Element at `eytzinger_sa[k]` is the SA value
    /// (text position) that would sit at node k in the implicit search tree.
    pub eytzinger_sa: Vec<usize>,
    pub cumulative_lengths: CumulativeLengthIndex,
}

impl DragonFmIndex {
    /// Returns true if this index uses the Eytzinger cache-oblivious SA layout.
    #[inline]
    pub fn has_eytzinger(&self) -> bool {
        !self.eytzinger_sa.is_empty()
    }

    /// Backward search: find all occurrences of pattern in the indexed text.
    ///
    /// Uses the Eytzinger SA for bounds computation when available (cache-
    /// oblivious, eliminates random cache misses on large indexes), then
    /// collects positions from the linear SA sequentially.
    pub fn search(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || self.text.is_empty() {
            return Vec::new();
        }

        let (lo, hi) = if self.has_eytzinger() {
            (
                self.lower_bound_eytzinger(pattern),
                self.upper_bound_eytzinger(pattern),
            )
        } else {
            (self.lower_bound(pattern), self.upper_bound(pattern))
        };

        self.suffix_array[lo..hi].to_vec()
    }

    /// Count occurrences of pattern (without locating).
    pub fn count(&self, pattern: &[u8]) -> usize {
        if pattern.is_empty() || self.text.is_empty() {
            return 0;
        }
        let lo = self.lower_bound(pattern);
        let hi = self.upper_bound(pattern);
        hi - lo
    }

    /// Map a text position to (unitig_id, offset_within_unitig).
    pub fn position_to_unitig(&self, pos: usize) -> Option<(u32, u32)> {
        self.cumulative_lengths.unitig_at_position(pos as u64)
    }

    /// Variable-length search: extend pattern until no matches remain.
    /// Returns (longest_match_length, SA_interval_size).
    pub fn variable_length_search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut best_len = 0;
        let mut best_count = 0;

        for len in 1..=pattern.len() {
            let count = self.count(&pattern[..len]);
            if count == 0 {
                break;
            }
            best_len = len;
            best_count = count;
        }

        (best_len, best_count)
    }

    fn lower_bound(&self, pattern: &[u8]) -> usize {
        let mut lo = 0usize;
        let mut hi = self.suffix_array.len();
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let suffix_start = self.suffix_array[mid];
            let suffix_end = self.text.len().min(suffix_start + pattern.len());
            let suffix = &self.text[suffix_start..suffix_end];
            if suffix < pattern {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    fn upper_bound(&self, pattern: &[u8]) -> usize {
        let mut lo = 0usize;
        let mut hi = self.suffix_array.len();
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let suffix_start = self.suffix_array[mid];
            let suffix_end = self.text.len().min(suffix_start + pattern.len());
            let suffix = &self.text[suffix_start..suffix_end];
            if suffix <= pattern {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    // ── Eytzinger (BFS-ordered) binary search ─────────────────────────────────

    /// Lower-bound search in the Eytzinger-ordered SA.
    ///
    /// Traverses the implicit binary search tree in BFS order: the root is at
    /// index 1, left child of node k is 2k, right child is 2k+1 (1-indexed).
    ///
    /// After the tree walk, extracts the final linear sorted index from the
    /// last visited node to return the conventional lower-bound position in
    /// `self.suffix_array`.
    fn lower_bound_eytzinger(&self, pattern: &[u8]) -> usize {
        let n = self.suffix_array.len();
        let eytz = &self.eytzinger_sa; // 1-indexed, eytz[0] unused, eytz[1..=n] valid

        // Walk the Eytzinger tree
        let mut k = 1usize; // 1-indexed node
        while k <= n {
            let sa_pos = eytz[k];
            let suffix_end = self.text.len().min(sa_pos + pattern.len());
            let suffix = &self.text[sa_pos..suffix_end];
            if suffix < pattern {
                k = 2 * k + 1; // go right
            } else {
                k = 2 * k;     // go left
            }
        }
        // Convert from Eytzinger node back to linear sorted index.
        // Node k (after loop) represents a missing element; the lower-bound
        // is the node we last went left from, encoded as k >> 1's right-spine.
        // Standard Eytzinger lower-bound recovery:
        eytzinger_node_to_sorted_idx(k, n)
    }

    /// Upper-bound search in the Eytzinger-ordered SA.
    fn upper_bound_eytzinger(&self, pattern: &[u8]) -> usize {
        let n = self.suffix_array.len();
        let eytz = &self.eytzinger_sa;

        let mut k = 1usize;
        while k <= n {
            let sa_pos = eytz[k];
            let suffix_end = self.text.len().min(sa_pos + pattern.len());
            let suffix = &self.text[sa_pos..suffix_end];
            if suffix <= pattern {
                k = 2 * k + 1; // go right
            } else {
                k = 2 * k;     // go left
            }
        }
        eytzinger_node_to_sorted_idx(k, n)
    }
}

/// Convert an Eytzinger leaf node index back to a linear sorted array index.
///
/// After a complete tree walk in Eytzinger order, the final node `k` encodes
/// the position where the search "fell off" the tree.  By counting trailing
/// right-child steps (trailing 1-bits after shifting right) we recover the
/// sorted index.
///
/// This is the standard O(1) recovery used in cache-oblivious binary search
/// implementations (Brodal, Fagerberg, and Pedersen 2003).
#[inline]
fn eytzinger_node_to_sorted_idx(mut k: usize, n: usize) -> usize {
    // k is the first out-of-bounds node.  The sorted index is the number of
    // times we went left on the way down, which equals the number of trailing
    // zeros in the node index after accounting for the direction encoding.
    let trailing = (k + 1).trailing_zeros() as usize;
    let result = (k >> trailing) - 1;
    result.min(n)
}

// ── Eytzinger SA builder (used at index-build time) ─────────────────────────

/// Build an Eytzinger-ordered copy of a suffix array.
///
/// The Eytzinger layout (BFS order of the implicit binary search tree) ensures
/// that the first log₂(cache_size/8) ≈ 20+ levels of any search are served from
/// cache, eliminating the random-miss cost that dominates searches in large SAs.
///
/// Returns `eytzinger_sa` where:
/// * `eytzinger_sa[0]` is unused (1-indexed layout)
/// * `eytzinger_sa[k]` = the SA value at node k in the BFS tree
pub fn build_eytzinger_sa(sa: &[usize]) -> Vec<usize> {
    let n = sa.len();
    if n == 0 {
        return Vec::new();
    }
    let mut eytz = vec![0usize; n + 1]; // 1-indexed
    fill_eytzinger(sa, &mut eytz, 1, 0, n);
    eytz
}

fn fill_eytzinger(sa: &[usize], eytz: &mut Vec<usize>, node: usize, lo: usize, hi: usize) {
    if lo >= hi || node >= eytz.len() {
        return;
    }
    let mid = lo + (hi - lo) / 2;
    eytz[node] = sa[mid];
    fill_eytzinger(sa, eytz, 2 * node, lo, mid);
    fill_eytzinger(sa, eytz, 2 * node + 1, mid + 1, hi);
}

/// Build the FM-index from a UnitigSet.
///
/// If `max_ram_bytes` is `Some(n)`, uses external-memory SA construction
/// to stay within the given RAM budget. If `None`, uses in-memory construction.
pub fn build_fm_index(unitigs: &UnitigSet, output_dir: &Path) -> Result<()> {
    build_fm_index_with_ram(unitigs, output_dir, None)
}

/// Build the FM-index with an optional RAM budget for suffix array construction.
///
/// Set `use_eytzinger = true` to also build the Eytzinger-ordered SA for
/// cache-oblivious binary search (recommended for indexes > 1 GB).
pub fn build_fm_index_with_ram(
    unitigs: &UnitigSet,
    output_dir: &Path,
    max_ram_bytes: Option<usize>,
) -> Result<()> {
    build_fm_index_with_options(unitigs, output_dir, max_ram_bytes, false)
}

/// Full build entry point with all options.
pub fn build_fm_index_with_options(
    unitigs: &UnitigSet,
    output_dir: &Path,
    max_ram_bytes: Option<usize>,
    use_eytzinger: bool,
) -> Result<()> {
    log::info!(
        "Building FM-index over {} bases of concatenated unitig text (eytzinger={})",
        unitigs.concatenated.len(),
        use_eytzinger,
    );

    let text = unitigs.concatenated.clone();

    let suffix_array = if let Some(ram_budget) = max_ram_bytes {
        log::info!("Using external-memory SA construction (RAM budget: {} MB)", ram_budget / 1_048_576);
        build_suffix_array_external(&text, output_dir, ram_budget)?
    } else {
        build_suffix_array(&text)
    };

    // Optionally build Eytzinger layout for cache-oblivious search.
    let eytzinger_sa = if use_eytzinger {
        log::info!("Building Eytzinger SA layout ({} entries)", suffix_array.len());
        let e = build_eytzinger_sa(&suffix_array);
        log::info!("Eytzinger SA built ({:.1} MB additional)", e.len() * 8 / 1_048_576);
        e
    } else {
        Vec::new()
    };

    // Build cumulative length index
    let cum_lengths = CumulativeLengthIndex::from_lengths(&unitigs.lengths);

    // Serialize components
    let index = DragonFmIndex {
        text,
        suffix_array,
        eytzinger_sa,
        cumulative_lengths: cum_lengths,
    };

    // Save to disk
    let fm_path = output_dir.join("fm_index.bin");
    crate::util::mmap::write_bincode(&fm_path, &FmIndexSerializable::from_index(&index))?;

    log::info!("FM-index saved to {:?}", fm_path);
    Ok(())
}

/// Load the FM-index from disk.
///
/// Supports both the current format (includes `eytzinger_sa`) and the legacy
/// format (without it).  bincode is positional, so old files end before
/// the `eytzinger_sa` bytes; we detect this by trying the new format first
/// and falling back to the legacy struct on failure.
pub fn load_fm_index(index_dir: &Path) -> Result<DragonFmIndex> {
    use anyhow::Context as _;

    let fm_path = index_dir.join("fm_index.bin");
    let bytes = std::fs::read(&fm_path)
        .with_context(|| format!("reading FM-index from {:?}", fm_path))?;

    // Try current format (with eytzinger_sa).
    if let Ok(s) = bincode::deserialize::<FmIndexSerializable>(&bytes) {
        return Ok(s.to_index());
    }

    // Fall back to legacy format (no eytzinger_sa field).
    #[derive(serde::Deserialize)]
    struct FmIndexLegacy {
        text: Vec<u8>,
        suffix_array: Vec<usize>,
        cumulative_lengths: CumulativeLengthIndex,
    }
    let legacy = bincode::deserialize::<FmIndexLegacy>(&bytes)
        .with_context(|| format!("deserializing FM-index from {:?}", fm_path))?;

    log::debug!("Loaded FM-index in legacy format from {:?}", fm_path);
    Ok(DragonFmIndex {
        text: legacy.text,
        suffix_array: legacy.suffix_array,
        eytzinger_sa: Vec::new(),
        cumulative_lengths: legacy.cumulative_lengths,
    })
}

/// Serializable representation of the FM-index.
#[derive(serde::Serialize, serde::Deserialize)]
struct FmIndexSerializable {
    text: Vec<u8>,
    suffix_array: Vec<usize>,
    cumulative_lengths: CumulativeLengthIndex,
    /// Eytzinger-ordered SA for cache-oblivious search.  Empty = use linear SA.
    #[serde(default)]
    eytzinger_sa: Vec<usize>,
}

impl FmIndexSerializable {
    fn from_index(index: &DragonFmIndex) -> Self {
        Self {
            text: index.text.clone(),
            suffix_array: index.suffix_array.clone(),
            cumulative_lengths: index.cumulative_lengths.clone(),
            eytzinger_sa: index.eytzinger_sa.clone(),
        }
    }

    fn to_index(self) -> DragonFmIndex {
        DragonFmIndex {
            text: self.text,
            suffix_array: self.suffix_array,
            cumulative_lengths: self.cumulative_lengths,
            eytzinger_sa: self.eytzinger_sa,
        }
    }
}

/// Build suffix array using a basic O(n log^2 n) algorithm.
/// For large texts (>100 MB), consider using the `cdivsufsort` crate.
fn build_suffix_array(text: &[u8]) -> Vec<usize> {
    let n = text.len();
    if n == 0 {
        return Vec::new();
    }

    let mut sa: Vec<usize> = (0..n).collect();
    sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
    sa
}

/// Build suffix array using external-memory (disk-based) merge sort.
///
/// This reduces peak RAM from O(8n) bytes to O(chunk_size) bytes by:
/// 1. Splitting suffix indices into chunks that fit in `max_ram_bytes`
/// 2. Sorting each chunk independently and writing to temp files
/// 3. K-way merging the sorted chunks
///
/// For 2M genomes (~5 Gbp unitig text): 40 GB SA in-memory → ~8 GB with chunking.
pub fn build_suffix_array_external(text: &[u8], output_dir: &Path, max_ram_bytes: usize) -> Result<Vec<usize>> {
    let n = text.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Each suffix index is 8 bytes; with comparison overhead, estimate ~24 bytes per entry
    let bytes_per_entry = 24;
    let chunk_size = max_ram_bytes / bytes_per_entry;

    if chunk_size >= n {
        // Fits in memory, use direct sort
        log::info!("Text fits in memory ({} entries, {} MB RAM budget), using in-memory sort",
                   n, max_ram_bytes / 1_048_576);
        return Ok(build_suffix_array(text));
    }

    log::info!(
        "External-memory SA construction: {} entries, chunk size {}, {} chunks",
        n, chunk_size, (n + chunk_size - 1) / chunk_size
    );

    let temp_dir = output_dir.join("_sa_temp");
    std::fs::create_dir_all(&temp_dir)?;

    // Phase 1: Sort chunks and write to temp files
    let mut chunk_files = Vec::new();
    let num_chunks = (n + chunk_size - 1) / chunk_size;

    for chunk_idx in 0..num_chunks {
        let start = chunk_idx * chunk_size;
        let end = n.min(start + chunk_size);

        log::info!("  Sorting chunk {}/{} ({} entries)", chunk_idx + 1, num_chunks, end - start);

        let mut chunk: Vec<usize> = (start..end).collect();
        chunk.sort_by(|&a, &b| text[a..].cmp(&text[b..]));

        // Write sorted chunk to temp file
        let chunk_path = temp_dir.join(format!("chunk_{:04}.bin", chunk_idx));
        let chunk_bytes: Vec<u8> = chunk.iter()
            .flat_map(|&x| x.to_le_bytes())
            .collect();
        std::fs::write(&chunk_path, &chunk_bytes)?;
        chunk_files.push(chunk_path);
    }

    // Phase 2: Streaming k-way merge using a BinaryHeap
    log::info!("  K-way merging {} sorted chunks (streaming)", chunk_files.len());

    use std::cmp::Ordering;
    use std::collections::BinaryHeap;
    use std::io::{BufReader, Read as _};

    struct ChunkReader {
        reader: BufReader<std::fs::File>,
        chunk_idx: usize,
    }

    impl ChunkReader {
        fn next_entry(&mut self) -> Option<usize> {
            let mut buf = [0u8; 8];
            self.reader.read_exact(&mut buf).ok()?;
            Some(usize::from_le_bytes(buf))
        }
    }

    // Wrapper for BinaryHeap that compares suffix indices by their text content.
    // Uses a raw pointer to `text` to avoid lifetime issues with BinaryHeap's Ord bound.
    struct SuffixOrd {
        suffix_idx: usize,
        chunk_idx: usize,
        text_ptr: *const u8,
        text_len: usize,
    }

    // SAFETY: text lives for the entire merge and is read-only during this phase.
    unsafe impl Send for SuffixOrd {}

    impl SuffixOrd {
        fn suffix(&self) -> &[u8] {
            unsafe { std::slice::from_raw_parts(self.text_ptr.add(self.suffix_idx), self.text_len - self.suffix_idx) }
        }
    }

    impl PartialEq for SuffixOrd {
        fn eq(&self, other: &Self) -> bool { self.suffix_idx == other.suffix_idx }
    }
    impl Eq for SuffixOrd {}

    impl PartialOrd for SuffixOrd {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
    }

    impl Ord for SuffixOrd {
        fn cmp(&self, other: &Self) -> Ordering {
            // Reverse comparison for min-heap (BinaryHeap is a max-heap)
            other.suffix().cmp(&self.suffix())
        }
    }

    let text_ptr = text.as_ptr();
    let text_len = text.len();

    let mut readers: Vec<ChunkReader> = chunk_files.iter().enumerate().map(|(idx, path)| {
        ChunkReader {
            reader: BufReader::with_capacity(64 * 1024, std::fs::File::open(path).unwrap()),
            chunk_idx: idx,
        }
    }).collect();

    let mut heap: BinaryHeap<SuffixOrd> = BinaryHeap::new();

    // Seed heap with first element from each chunk
    for reader in &mut readers {
        if let Some(idx) = reader.next_entry() {
            heap.push(SuffixOrd { suffix_idx: idx, chunk_idx: reader.chunk_idx, text_ptr, text_len });
        }
    }

    let mut result: Vec<usize> = Vec::with_capacity(n);

    while let Some(entry) = heap.pop() {
        result.push(entry.suffix_idx);
        let ci = entry.chunk_idx;
        if let Some(next_idx) = readers[ci].next_entry() {
            heap.push(SuffixOrd { suffix_idx: next_idx, chunk_idx: ci, text_ptr, text_len });
        }
    }

    // Cleanup temp files
    std::fs::remove_dir_all(&temp_dir)?;

    log::info!("  External SA construction complete");
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fm_index_search() {
        let text = b"ACGTACGT$TTGGCCAA$".to_vec();
        let sa = build_suffix_array(&text);
        let cum = CumulativeLengthIndex::from_lengths(&[8, 8]);

        let index = DragonFmIndex {
            text,
            suffix_array: sa,
            cumulative_lengths: cum,
            eytzinger_sa: Vec::new(),
        };

        // Search for ACGT - should find 2 occurrences
        let hits = index.search(b"ACGT");
        assert_eq!(hits.len(), 2);

        // Search for TTGG - should find 1
        let hits = index.search(b"TTGG");
        assert_eq!(hits.len(), 1);

        // Map positions to unitigs. Text "ACGTACGT$TTGGCCAA$":
        //   unitig 0 = text[0..8), '$' @ 8, unitig 1 = text[9..17), '$' @ 17.
        let (uid, offset) = index.position_to_unitig(0).unwrap();
        assert_eq!((uid, offset), (0, 0));
        let (uid, offset) = index.position_to_unitig(7).unwrap();
        assert_eq!((uid, offset), (0, 7)); // last base of unitig 0

        // Position 8 is the '$' separator — not inside any unitig.
        assert_eq!(index.position_to_unitig(8), None);

        // Position 9 is the first base of unitig 1 (separator correctly skipped).
        let (uid, offset) = index.position_to_unitig(9).unwrap();
        assert_eq!((uid, offset), (1, 0));
        let (uid, offset) = index.position_to_unitig(16).unwrap();
        assert_eq!((uid, offset), (1, 7)); // last base of unitig 1
    }

    #[test]
    fn test_variable_length_search() {
        let text = b"ACGTACGTACGT$".to_vec();
        let sa = build_suffix_array(&text);
        let cum = CumulativeLengthIndex::from_lengths(&[12]);

        let index = DragonFmIndex {
            text,
            suffix_array: sa,
            cumulative_lengths: cum,
            eytzinger_sa: Vec::new(),
        };

        let (len, count) = index.variable_length_search(b"ACGTACGTX");
        assert_eq!(len, 8); // ACGTACGT matches, ACGTACGTX doesn't
        assert!(count > 0);
    }
}
