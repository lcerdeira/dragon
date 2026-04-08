# Module reference

## Index modules (`src/index/`)

### `index::dbg`

De Bruijn graph construction via GGCAT or internal builder.

- `build_cdbg(genome_dir, output_dir, kmer_size, threads)` &mdash; Build a coloured compacted de Bruijn graph. Uses GGCAT if available, otherwise falls back to internal builder.
- `DbgResult` &mdash; Result struct containing paths to unitig and colour files.

### `index::unitig`

Unitig parsing and 2-bit encoding.

- `parse_and_encode_unitigs(path)` &mdash; Parse a unitig FASTA file and encode sequences.
- `UnitigSet` &mdash; Collection of all unitigs with concatenated text and length metadata.
- `Unitig` &mdash; Single unitig with ID and 2-bit packed sequence.

### `index::color`

Roaring Bitmap colour index.

- `build_color_index(color_file, output_dir, num_genomes)` &mdash; Build and serialise the colour index.
- `load_color_index(index_dir)` &mdash; Load colour index via memory mapping.
- `ColorIndex::get_colors(unitig_id)` &mdash; Look up which genomes contain a unitig.

### `index::fm`

FM-index construction and querying.

- `build_fm_index(unitigs, output_dir)` &mdash; Build FM-index from a UnitigSet.
- `load_fm_index(index_dir)` &mdash; Load FM-index from disk.
- `DragonFmIndex::search(pattern)` &mdash; Find all occurrences of a pattern.
- `DragonFmIndex::count(pattern)` &mdash; Count occurrences without locating.
- `DragonFmIndex::variable_length_search(pattern)` &mdash; Extend search to maximum match length.

### `index::paths`

Genome path index.

- `build_path_index(genome_dir, unitigs, output_dir)` &mdash; Build path index from genomes.
- `load_path_index(index_dir)` &mdash; Load path index from disk.
- `PathIndex::extract_sequence(genome_id, start, end, unitigs)` &mdash; Reconstruct a genome region.

### `index::metadata`

Index statistics and metadata.

- `write_metadata(output_dir, dbg_result, unitigs)` &mdash; Write metadata JSON.
- `load_metadata(index_dir)` &mdash; Load metadata.

---

## Query modules (`src/query/`)

### `query::seed`

FM-index seed finding.

- `find_seeds(query, fm_index, min_seed_len, max_freq)` &mdash; Find all seeds in a query using backward search with variable-length extension. Searches both forward and reverse complement.

### `query::candidate`

Candidate genome filtering.

- `find_candidates(seeds, color_index, min_votes)` &mdash; Identify genomes sharing unitigs with query seeds. Returns candidates sorted by vote count.

### `query::chain`

Colinear chaining.

- `chain_candidates(seeds, candidates, path_index, min_score)` &mdash; Compute optimal colinear chains for each candidate genome using Fenwick tree DP.
- `Chain` &mdash; A scored chain of colinear anchors with coverage information.

### `query::align`

Wavefront alignment.

- `align_chains(query, chains, path_index)` &mdash; Align chains and produce PAF records.
- `banded_nw_align(query, reference, bandwidth)` &mdash; Banded Needleman-Wunsch alignment.

---

## Data structures (`src/ds/`)

### `ds::fenwick`

- `FenwickMax` &mdash; Prefix maximum queries in O(log n).
- `FenwickSum` &mdash; Prefix sum queries in O(log n).

### `ds::elias_fano`

- `CumulativeLengthIndex` &mdash; Maps text positions to unitig IDs via binary search on cumulative lengths.

### `ds::varint`

- `encode_varint / decode_varint` &mdash; LEB128 variable-length integer encoding.
- `encode_zigzag / decode_zigzag` &mdash; Zigzag encoding for signed integers.
- `delta_encode / delta_decode` &mdash; Delta + varint encoding for sorted sequences.

---

## Utilities (`src/util/`)

### `util::dna`

- `PackedSequence` &mdash; 2-bit packed DNA sequence (32 bases per u64).
- `canonical_kmer(kmer, k)` &mdash; Lexicographically smaller of forward and reverse complement.

### `util::mmap`

- `mmap_open(path)` &mdash; Memory-map a file for read-only access.
- `read_bincode / write_bincode` &mdash; Serialise/deserialise via bincode.

---

## I/O modules (`src/io/`)

### `io::fasta`

- `read_sequences(path)` &mdash; Read all sequences from a FASTA file.
- `FastaReader` &mdash; Streaming iterator over FASTA records.
- `list_fasta_files(dir)` &mdash; List FASTA files in a directory.

### `io::paf`

- `PafRecord` &mdash; PAF alignment record with Display formatting.
- `write_paf(writer, records)` &mdash; Write PAF records.

### `io::blast`

- `write_blast_tabular(writer, records)` &mdash; Write BLAST outfmt 6 records.
