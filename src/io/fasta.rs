/// Streaming FASTA/FASTQ parser with transparent gzip support.

use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// A named DNA sequence.
#[derive(Clone, Debug)]
pub struct Sequence {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Open a path as a buffered reader, transparently decompressing gzip.
///
/// Detects gzip by the `.gz` extension OR the magic bytes (0x1f 0x8b), so it
/// works for `.fa.gz`, `.fasta.gz`, `.fna.gz`, and also plain files renamed
/// without the extension.  AllTheBacteria / NCBI / EBI genomes are distributed
/// as gzip, so this is the common case for real-world inputs.
///
/// `flate2::MultiGzDecoder` handles multi-member gzip streams (concatenated
/// `.gz` blocks), which some bulk genome archives use.
pub fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;

    // Peek the first two bytes for the gzip magic number.
    let mut magic = [0u8; 2];
    let is_gzip = match file.read_exact(&mut magic) {
        Ok(()) => magic == [0x1f, 0x8b],
        Err(_) => false, // file shorter than 2 bytes — treat as plain
    };
    // Rewind to the start so the chosen reader sees the whole file.
    use std::io::Seek;
    file.seek(std::io::SeekFrom::Start(0))?;

    if is_gzip {
        let decoder = flate2::read::MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Read all sequences from a FASTA file (plain or gzipped).
pub fn read_sequences(path: &Path) -> Result<Vec<Sequence>> {
    let reader = open_maybe_gzip(path)?;
    let mut sequences = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }

        if let Some(header) = line.strip_prefix('>') {
            // Save previous sequence
            if let Some(name) = current_name.take() {
                if !current_seq.is_empty() {
                    sequences.push(Sequence {
                        name,
                        seq: current_seq.clone(),
                    });
                }
            }
            // Start new sequence
            let name = header.split_whitespace().next().unwrap_or(header);
            current_name = Some(name.to_string());
            current_seq.clear();
        } else if line.starts_with('@') {
            // FASTQ header - handle similarly
            if let Some(name) = current_name.take() {
                if !current_seq.is_empty() {
                    sequences.push(Sequence {
                        name,
                        seq: current_seq.clone(),
                    });
                }
            }
            let name = line[1..].split_whitespace().next().unwrap_or(&line[1..]);
            current_name = Some(name.to_string());
            current_seq.clear();
        } else if line.starts_with('+') {
            // FASTQ quality header - skip next line
            continue;
        } else {
            // Filter to only valid DNA bases (strip GGCAT's trailing '$' and other non-DNA chars)
            current_seq.extend(
                line.bytes().filter(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'N' | b'n'))
            );
        }
    }

    // Save last sequence
    if let Some(name) = current_name {
        if !current_seq.is_empty() {
            sequences.push(Sequence { name, seq: current_seq });
        }
    }

    Ok(sequences)
}

/// Iterator-based FASTA reader for streaming large files (plain or gzipped).
pub struct FastaReader {
    reader: Box<dyn BufRead>,
    buf: String,
    peeked_header: Option<String>,
}

impl FastaReader {
    pub fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            reader: open_maybe_gzip(path)?,
            buf: String::new(),
            peeked_header: None,
        })
    }
}

impl Iterator for FastaReader {
    type Item = Result<Sequence>;

    fn next(&mut self) -> Option<Self::Item> {
        let name = if let Some(header) = self.peeked_header.take() {
            header
        } else {
            // Find next header
            loop {
                self.buf.clear();
                match self.reader.read_line(&mut self.buf) {
                    Ok(0) => return None,
                    Ok(_) => {
                        let line = self.buf.trim();
                        if let Some(header) = line.strip_prefix('>') {
                            break header.split_whitespace().next().unwrap_or(header).to_string();
                        }
                    }
                    Err(e) => return Some(Err(e.into())),
                }
            }
        };

        let mut seq = Vec::new();
        loop {
            self.buf.clear();
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => break,
                Ok(_) => {
                    let line = self.buf.trim();
                    if let Some(header) = line.strip_prefix('>') {
                        self.peeked_header = Some(
                            header.split_whitespace().next().unwrap_or(header).to_string(),
                        );
                        break;
                    }
                    seq.extend_from_slice(line.as_bytes());
                }
                Err(e) => return Some(Err(e.into())),
            }
        }

        Some(Ok(Sequence { name, seq }))
    }
}

/// List all FASTA files in a directory.
pub fn list_fasta_files(dir: &Path) -> Result<Vec<std::path::PathBuf>> {
    let mut files = Vec::new();
    list_fasta_files_recursive(dir, &mut files)?;
    files.sort();
    Ok(files)
}

fn list_fasta_files_recursive(dir: &Path, files: &mut Vec<std::path::PathBuf>) -> Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            list_fasta_files_recursive(&path, files)?;
        } else {
            // Accept plain FASTA (.fa/.fasta/.fna/.fsa) AND gzipped variants
            // (.fa.gz etc.).  We match on the full filename so that the
            // double-extension `.fa.gz` is recognised (path.extension() would
            // only return "gz").
            let name = path.file_name()
                .map(|n| n.to_string_lossy().to_lowercase())
                .unwrap_or_default();
            let is_fasta = name.ends_with(".fa")
                || name.ends_with(".fasta")
                || name.ends_with(".fna")
                || name.ends_with(".fsa")
                || name.ends_with(".fa.gz")
                || name.ends_with(".fasta.gz")
                || name.ends_with(".fna.gz")
                || name.ends_with(".fsa.gz");
            if is_fasta {
                files.push(path);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_fasta() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fa");
        let mut f = File::create(&path).unwrap();
        writeln!(f, ">seq1 description").unwrap();
        writeln!(f, "ACGTACGT").unwrap();
        writeln!(f, "TTGGCCAA").unwrap();
        writeln!(f, ">seq2").unwrap();
        writeln!(f, "AAAACCCC").unwrap();

        let seqs = read_sequences(&path).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "seq1");
        assert_eq!(seqs[0].seq, b"ACGTACGTTTGGCCAA");
        assert_eq!(seqs[1].name, "seq2");
        assert_eq!(seqs[1].seq, b"AAAACCCC");
    }

    #[test]
    fn test_read_gzipped_fasta() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fa.gz");
        let f = File::create(&path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        writeln!(enc, ">seq1 description").unwrap();
        writeln!(enc, "ACGTACGT").unwrap();
        writeln!(enc, "TTGGCCAA").unwrap();
        writeln!(enc, ">seq2").unwrap();
        writeln!(enc, "AAAACCCC").unwrap();
        enc.finish().unwrap();

        // read_sequences must transparently decompress .fa.gz
        let seqs = read_sequences(&path).unwrap();
        assert_eq!(seqs.len(), 2, "gzip FASTA should yield 2 sequences");
        assert_eq!(seqs[0].name, "seq1");
        assert_eq!(seqs[0].seq, b"ACGTACGTTTGGCCAA");
        assert_eq!(seqs[1].seq, b"AAAACCCC");
    }

    #[test]
    fn test_list_includes_gzipped() {
        let dir = tempfile::tempdir().unwrap();
        // Mix of plain and gzipped FASTA + a non-FASTA file
        File::create(dir.path().join("a.fa")).unwrap();
        File::create(dir.path().join("b.fa.gz")).unwrap();
        File::create(dir.path().join("c.fasta.gz")).unwrap();
        File::create(dir.path().join("d.fna")).unwrap();
        File::create(dir.path().join("readme.txt")).unwrap();

        let files = list_fasta_files(dir.path()).unwrap();
        assert_eq!(files.len(), 4, "should list 4 FASTA (plain + gz), not the .txt");
        let names: Vec<String> = files
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        assert!(names.contains(&"b.fa.gz".to_string()), "must include .fa.gz");
        assert!(names.contains(&"c.fasta.gz".to_string()), "must include .fasta.gz");
        assert!(!names.iter().any(|n| n.ends_with(".txt")), "must exclude .txt");
    }

    #[test]
    fn test_gzip_magic_byte_detection() {
        // A gzipped file WITHOUT .gz extension must still be detected by magic bytes
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nogzext.fa"); // .fa but actually gzipped
        let f = File::create(&path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        writeln!(enc, ">s\nACGT").unwrap();
        enc.finish().unwrap();

        let seqs = read_sequences(&path).unwrap();
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].seq, b"ACGT");
    }
}
