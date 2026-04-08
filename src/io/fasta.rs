/// Streaming FASTA/FASTQ parser.

use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A named DNA sequence.
#[derive(Clone, Debug)]
pub struct Sequence {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Read all sequences from a FASTA file.
pub fn read_sequences(path: &Path) -> Result<Vec<Sequence>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
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
            current_seq.extend_from_slice(line.as_bytes());
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

/// Iterator-based FASTA reader for streaming large files.
pub struct FastaReader {
    reader: BufReader<File>,
    buf: String,
    peeked_header: Option<String>,
}

impl FastaReader {
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::new(file),
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
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if let Some(ext) = path.extension() {
            let ext = ext.to_string_lossy().to_lowercase();
            if matches!(ext.as_str(), "fa" | "fasta" | "fna" | "fsa") {
                files.push(path);
            }
        }
    }
    files.sort();
    Ok(files)
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
}
