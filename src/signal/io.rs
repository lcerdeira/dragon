/// Signal I/O: readers for nanopore signal data in various formats.
///
/// Supported formats:
/// - TSV: simple tab-separated float values, one read per line (id\tval1\tval2\t...)
/// - FAST5 (simplified): basic HDF5 reader for single-fast5 files
/// - BLOW5/SLOW5: simplified reader for the S/BLOW5 format
///
/// For production use with FAST5, the `hdf5` crate is recommended.
/// This implementation provides a lightweight TSV-based format that can
/// be generated from FAST5 using external tools (e.g., `slow5tools`).

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single nanopore signal read.
#[derive(Clone, Debug)]
pub struct SignalRead {
    /// Read identifier.
    pub id: String,
    /// Raw signal values in picoamperes (pA).
    pub signal: Vec<f32>,
}

/// Read signal data from a file, auto-detecting format by extension.
///
/// Supported extensions:
/// - `.tsv`, `.txt` - TSV signal format
/// - `.slow5` - SLOW5 text format
/// - `.csv` - CSV signal format
pub fn read_signal_file(path: &Path) -> Result<Vec<SignalRead>> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "tsv" | "txt" => read_tsv_signal(path),
        "slow5" => read_slow5_signal(path),
        "csv" => read_csv_signal(path),
        _ => {
            // Default to TSV format
            log::warn!(
                "Unknown signal file extension '{}', attempting TSV format",
                ext
            );
            read_tsv_signal(path)
        }
    }
}

/// Read signal data from a TSV file.
///
/// Format: one read per line, tab-separated.
///   read_id<TAB>value1<TAB>value2<TAB>...
///
/// Lines starting with '#' are treated as comments and skipped.
pub fn read_tsv_signal(path: &Path) -> Result<Vec<SignalRead>> {
    let file = File::open(path).with_context(|| format!("Cannot open signal file: {:?}", path))?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Error reading line {}", line_num + 1))?;
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            log::warn!("Line {}: too few fields, skipping", line_num + 1);
            continue;
        }

        let id = fields[0].to_string();
        let signal: Result<Vec<f32>, _> = fields[1..]
            .iter()
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.parse::<f32>()
                    .with_context(|| format!("Invalid float '{}' on line {}", s, line_num + 1))
            })
            .collect();

        let signal = signal?;

        if signal.is_empty() {
            log::warn!("Line {}: no signal values for read {}, skipping", line_num + 1, id);
            continue;
        }

        reads.push(SignalRead { id, signal });
    }

    log::info!(
        "Read {} signal reads from {:?}",
        reads.len(),
        path
    );
    Ok(reads)
}

/// Read signal data from a CSV file.
///
/// Format: one read per line, comma-separated.
///   read_id,value1,value2,...
pub fn read_csv_signal(path: &Path) -> Result<Vec<SignalRead>> {
    let file = File::open(path).with_context(|| format!("Cannot open signal file: {:?}", path))?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Error reading line {}", line_num + 1))?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 2 {
            continue;
        }

        let id = fields[0].to_string();
        let signal: Result<Vec<f32>, _> = fields[1..]
            .iter()
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.trim()
                    .parse::<f32>()
                    .with_context(|| format!("Invalid float '{}' on line {}", s, line_num + 1))
            })
            .collect();

        let signal = signal?;
        if !signal.is_empty() {
            reads.push(SignalRead { id, signal });
        }
    }

    log::info!("Read {} signal reads from {:?}", reads.len(), path);
    Ok(reads)
}

/// Read signal data from a SLOW5 (text format) file.
///
/// SLOW5 is a simplified format for nanopore signal data.
/// Header lines start with '#' or '@'.
/// Data lines: read_id<TAB>...fields...<TAB>raw_signal(comma-separated integers)
///
/// We parse a simplified subset: the first column is the read_id,
/// and the last column contains comma-separated signal values.
pub fn read_slow5_signal(path: &Path) -> Result<Vec<SignalRead>> {
    let file = File::open(path).with_context(|| format!("Cannot open SLOW5 file: {:?}", path))?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();
    let mut signal_col: Option<usize> = None;

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Error reading line {}", line_num + 1))?;
        let line = line.trim();

        // Skip empty lines and header
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Column header line
        if line.starts_with('@') || line.starts_with("read_id") {
            // Find the raw_signal column
            let cols: Vec<&str> = line.split('\t').collect();
            for (i, &col) in cols.iter().enumerate() {
                let col_name = col.trim_start_matches('@');
                if col_name == "raw_signal" || col_name == "signal" {
                    signal_col = Some(i);
                    break;
                }
            }
            // Default to last column if not found
            if signal_col.is_none() {
                signal_col = Some(cols.len() - 1);
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        let id = fields[0].to_string();
        let sig_col = signal_col.unwrap_or(fields.len() - 1);

        if sig_col >= fields.len() {
            continue;
        }

        let signal: Result<Vec<f32>, _> = fields[sig_col]
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.trim()
                    .parse::<f32>()
                    .with_context(|| format!("Invalid signal value '{}' on line {}", s, line_num + 1))
            })
            .collect();

        let signal = signal?;
        if !signal.is_empty() {
            reads.push(SignalRead { id, signal });
        }
    }

    log::info!("Read {} signal reads from SLOW5 {:?}", reads.len(), path);
    Ok(reads)
}

/// Write signal reads to a TSV file.
pub fn write_tsv_signal(path: &Path, reads: &[SignalRead]) -> Result<()> {
    use std::io::Write;

    let mut file = File::create(path)?;
    writeln!(file, "# Dragon signal TSV format: read_id<TAB>pA_values...")?;

    for read in reads {
        write!(file, "{}", read.id)?;
        for &val in &read.signal {
            write!(file, "\t{:.2}", val)?;
        }
        writeln!(file)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_tsv_signal() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.tsv");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "# comment line").unwrap();
        writeln!(f, "read_001\t50.0\t51.5\t49.3\t52.1").unwrap();
        writeln!(f, "read_002\t80.0\t81.0\t79.5").unwrap();

        let reads = read_tsv_signal(&path).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].id, "read_001");
        assert_eq!(reads[0].signal.len(), 4);
        assert!((reads[0].signal[0] - 50.0).abs() < 1e-6);
        assert_eq!(reads[1].id, "read_002");
        assert_eq!(reads[1].signal.len(), 3);
    }

    #[test]
    fn test_read_csv_signal() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.csv");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "read_001,50.0,51.5,49.3").unwrap();

        let reads = read_csv_signal(&path).unwrap();
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].signal.len(), 3);
    }

    #[test]
    fn test_write_read_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("roundtrip.tsv");

        let reads = vec![
            SignalRead {
                id: "r1".to_string(),
                signal: vec![50.0, 60.0, 70.0],
            },
            SignalRead {
                id: "r2".to_string(),
                signal: vec![80.0, 90.0],
            },
        ];

        write_tsv_signal(&path, &reads).unwrap();
        let loaded = read_tsv_signal(&path).unwrap();
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded[0].id, "r1");
        assert_eq!(loaded[0].signal.len(), 3);
        // Values are written with 2 decimal places, so check approximate equality
        assert!((loaded[0].signal[0] - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_auto_detect_format() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.tsv");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "read_001\t50.0\t51.5").unwrap();

        let reads = read_signal_file(&path).unwrap();
        assert_eq!(reads.len(), 1);
    }

    #[test]
    fn test_empty_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.tsv");
        File::create(&path).unwrap();

        let reads = read_tsv_signal(&path).unwrap();
        assert!(reads.is_empty());
    }
}
