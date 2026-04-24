/// Direct parser for GGCAT's `unitigs.colors.dat` binary colormap format.
///
/// Bypasses `ggcat query --colors` (which produces 500+ GB JSONL for large
/// indices) by reading the compact binary colormap directly. For a 26K genome
/// S. aureus build, this is 11 GB of binary data vs 665 GB JSONL — a 60x
/// reduction in intermediate storage.
///
/// File format (authoritative source: ggcat/crates/colors/src/storage/):
///
/// ```text
/// [0..64)                    ColorsFileHeader (Desse, packed LE)
/// [64..color_names_end)      LZ4 frame: bincode(Vec<String>) color names
/// [color_names_end..idx)     Independent LZ4 frames, each a chunk of subsets
/// [index_offset..EOF)        bincode(ColorsIndexMap): chunk start offsets
/// ```
///
/// Each chunk's LZ4 frame decodes to a concatenation of RunLength-encoded
/// color subsets. A color subset is a sorted list of genome IDs encoded with
/// varint + 2nd-order RLE (see decode_subset below).
///
/// The `C:hex_id:count` tags in `unitigs.fa` are `u32` subset indices in
/// hexadecimal, range `[0, subsets_count)`. To look up a unitig's genomes:
///
/// 1. Parse `C:hex_id` → decimal u32
/// 2. Binary-search `index.pairs` for the chunk containing this hex_id
/// 3. LZ4-decompress that chunk; skip (hex_id - chunk.start_index) subsets
/// 4. Decode the target subset → Vec<genome_id>

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};

use anyhow::{anyhow, bail, Context, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use serde::Deserialize;

const EXPECTED_MAGIC: &[u8; 16] = b"GGCAT_CMAP_RNLEN";
const HEADER_SIZE: usize = 64; // 16-byte magic + 6 × u64

/// Parsed GGCAT colormap header.
#[derive(Debug, Clone)]
pub struct GgcatColorHeader {
    pub version: u64,
    pub index_offset: u64,
    pub colors_count: u64,    // Number of genomes
    pub subsets_count: u64,   // Number of distinct color subsets
    pub total_size: u64,
    pub total_uncompressed_size: u64,
}

#[derive(Deserialize, Debug)]
struct IndexEntry {
    start_index: u32,
    file_offset: u64,
}

#[derive(Deserialize, Debug)]
struct IndexMap {
    pairs: Vec<IndexEntry>,
    subsets_count: u64,
}

/// Read and verify the GGCAT colormap header.
pub fn read_header(path: &std::path::Path) -> Result<GgcatColorHeader> {
    let mut f = File::open(path).with_context(|| format!("opening {:?}", path))?;
    let mut buf = [0u8; HEADER_SIZE];
    f.read_exact(&mut buf)?;

    if &buf[0..16] != EXPECTED_MAGIC {
        bail!(
            "Bad GGCAT colormap magic: expected {:?}, got {:?}",
            std::str::from_utf8(EXPECTED_MAGIC).unwrap(),
            std::str::from_utf8(&buf[0..16]).unwrap_or("<non-utf8>")
        );
    }

    let mut c = std::io::Cursor::new(&buf[16..]);
    let version = c.read_u64::<LittleEndian>()?;
    let index_offset = c.read_u64::<LittleEndian>()?;
    let colors_count = c.read_u64::<LittleEndian>()?;
    let subsets_count = c.read_u64::<LittleEndian>()?;
    let total_size = c.read_u64::<LittleEndian>()?;
    let total_uncompressed_size = c.read_u64::<LittleEndian>()?;

    if version != 1 {
        bail!("Unsupported GGCAT colormap version {} (expected 1)", version);
    }

    Ok(GgcatColorHeader {
        version,
        index_offset,
        colors_count,
        subsets_count,
        total_size,
        total_uncompressed_size,
    })
}

/// Read a varint (LEB128 unsigned, little-endian base-128).
#[inline]
fn read_varint<R: Read>(r: &mut R) -> std::io::Result<u64> {
    let mut v = 0u64;
    let mut shift = 0u32;
    loop {
        let b = r.read_u8()?;
        v |= ((b & 0x7f) as u64) << shift;
        if b & 0x80 == 0 {
            return Ok(v);
        }
        shift += 7;
        if shift >= 64 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "varint too large",
            ));
        }
    }
}

/// Decode one RunLength-encoded color subset from the decompressed stream.
/// Appends genome IDs to `out`. Returns on success or error.
fn decode_subset<R: Read>(r: &mut R, out: &mut Vec<u32>) -> std::io::Result<()> {
    // First genome ID (biased by +2)
    let first = read_varint(r)?;
    if first < 2 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid first varint in subset",
        ));
    }
    let mut last = (first - 2) as u32;
    out.push(last);

    loop {
        let x = read_varint(r)? as u32;
        if x == 0 {
            return Ok(()); // terminator
        }
        if x == 1 {
            // Run-length encoded gap
            let gap = read_varint(r)? as u32;
            let count = read_varint(r)? as u32;
            for _ in 0..count {
                last += gap;
                out.push(last);
            }
        } else {
            last += x - 1;
            out.push(last);
        }
    }
}

/// Read the index map at the given offset.
fn read_index_map(f: &mut File, offset: u64) -> Result<IndexMap> {
    f.seek(SeekFrom::Start(offset))?;
    let reader = BufReader::new(f);
    bincode::deserialize_from(reader)
        .map_err(|e| anyhow!("Failed to deserialize IndexMap: {}", e))
}

/// Parse all color subsets, calling `cb(hex_id, genome_ids)` for each.
/// Memory-efficient streaming parser — no large HashMap materialized.
pub fn for_each_subset<F>(path: &std::path::Path, mut cb: F) -> Result<GgcatColorHeader>
where
    F: FnMut(u32, &[u32]),
{
    let header = read_header(path)?;
    let mut f = File::open(path)?;
    let index = read_index_map(&mut f, header.index_offset)?;

    let pairs = index.pairs;
    let total_subsets = index.subsets_count as u32;

    let mut genomes = Vec::new();

    for (ci, entry) in pairs.iter().enumerate() {
        let end = pairs
            .get(ci + 1)
            .map(|p| p.start_index)
            .unwrap_or(total_subsets);
        let n = (end - entry.start_index) as usize;

        f.seek(SeekFrom::Start(entry.file_offset))?;

        // Read enough bytes for this chunk (until next chunk or index_offset)
        let chunk_end = pairs
            .get(ci + 1)
            .map(|p| p.file_offset)
            .unwrap_or(header.index_offset);
        let chunk_size = (chunk_end - entry.file_offset) as usize;
        let mut chunk_bytes = vec![0u8; chunk_size];
        f.read_exact(&mut chunk_bytes)?;

        // Decompress with lz4_flex frame
        let mut decoder = lz4_flex::frame::FrameDecoder::new(&chunk_bytes[..]);

        for i in 0..n {
            let hex_id = entry.start_index + i as u32;
            genomes.clear();
            decode_subset(&mut decoder, &mut genomes)?;
            cb(hex_id, &genomes);
        }
    }

    Ok(header)
}

/// Parse every color subset into a HashMap.
///
/// **WARNING**: For large colormaps (millions of subsets), this materializes
/// everything in memory. Use `for_each_subset` for streaming instead.
pub fn parse_all(path: &std::path::Path) -> Result<(GgcatColorHeader, HashMap<u32, Vec<u32>>)> {
    let mut out = HashMap::new();
    let header = for_each_subset(path, |hex_id, genomes| {
        out.insert(hex_id, genomes.to_vec());
    })?;
    Ok((header, out))
}

/// Write a colors.tsv file directly from a GGCAT colormap + unitigs.fa.
///
/// Reads `unitigs.fa` to get the `C:hex_id:count` tags per unitig, decodes
/// each hex_id from the colormap, and writes TSV: `unitig_id\tgenome1,genome2,...`
///
/// This is the main entry point to replace `ggcat query --colors` for Dragon's
/// color index construction.
pub fn write_colors_tsv(
    colormap_path: &std::path::Path,
    unitigs_path: &std::path::Path,
    output_tsv: &std::path::Path,
) -> Result<GgcatColorHeader> {
    use std::io::{BufRead, Write};

    log::info!("Parsing GGCAT colormap {:?}", colormap_path);
    let (header, subsets) = parse_all(colormap_path)?;
    log::info!(
        "  Colormap: {} genomes, {} color subsets, {} MB disk",
        header.colors_count,
        header.subsets_count,
        header.total_size / 1_048_576
    );

    log::info!("Parsing unitig FASTA headers {:?}", unitigs_path);
    let unitig_reader = BufReader::new(File::open(unitigs_path)?);
    let mut tsv_writer = std::io::BufWriter::new(File::create(output_tsv)?);

    let mut unitig_count = 0u64;
    for line in unitig_reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            continue;
        }
        // Parse: >unitig_id LN:i:len C:hex_id:count ...
        let parts: Vec<&str> = line[1..].split_whitespace().collect();
        let unitig_id: u64 = parts.first()
            .and_then(|s| s.parse().ok())
            .unwrap_or(unitig_count);

        // Find first C:hex_id:count tag (most unitigs have exactly one
        // dominant color set)
        let mut written = false;
        for part in &parts[1..] {
            if let Some(rest) = part.strip_prefix("C:") {
                if let Some(hex_id_str) = rest.split(':').next() {
                    if let Ok(hex_id) = u32::from_str_radix(hex_id_str, 16) {
                        if let Some(genomes) = subsets.get(&hex_id) {
                            let ids: Vec<String> = genomes.iter().map(|g| g.to_string()).collect();
                            writeln!(tsv_writer, "{}\t{}", unitig_id, ids.join(","))?;
                            written = true;
                            break;
                        }
                    }
                }
            }
        }
        if !written {
            writeln!(tsv_writer, "{}\t", unitig_id)?;
        }
        unitig_count += 1;

        if unitig_count % 1_000_000 == 0 {
            log::info!("  Processed {} unitigs...", unitig_count);
        }
    }

    log::info!("Wrote {} unitig color entries to {:?}", unitig_count, output_tsv);
    Ok(header)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_varint_roundtrip() {
        // varint 0 → 0x00
        let mut r = &[0u8][..];
        assert_eq!(read_varint(&mut r).unwrap(), 0);

        // varint 127 → 0x7f
        let mut r = &[0x7f][..];
        assert_eq!(read_varint(&mut r).unwrap(), 127);

        // varint 128 → 0x80 0x01
        let mut r = &[0x80, 0x01][..];
        assert_eq!(read_varint(&mut r).unwrap(), 128);

        // varint 300 → 0xac 0x02
        let mut r = &[0xac, 0x02][..];
        assert_eq!(read_varint(&mut r).unwrap(), 300);
    }

    #[test]
    fn test_decode_subset_single_genome() {
        // Subset with just genome ID 5: varint(5+2=7), varint(0)
        let data = vec![7u8, 0u8];
        let mut out = Vec::new();
        decode_subset(&mut &data[..], &mut out).unwrap();
        assert_eq!(out, vec![5]);
    }

    #[test]
    fn test_decode_subset_multiple_genomes() {
        // Subset [3, 7, 12]: first=3+2=5, gap1=(7-3)+1=5, gap2=(12-7)+1=6, term=0
        let data = vec![5u8, 5, 6, 0];
        let mut out = Vec::new();
        decode_subset(&mut &data[..], &mut out).unwrap();
        assert_eq!(out, vec![3, 7, 12]);
    }

    #[test]
    fn test_magic_validation() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"WRONG_MAGIC_HERE").unwrap();
        let result = read_header(tmp.path());
        assert!(result.is_err());
    }
}
