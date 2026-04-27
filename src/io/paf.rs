/// PAF (Pairwise Alignment Format) output writer.

use std::fmt;
use std::io::Write;

/// A single PAF record.
#[derive(Clone, Debug)]
pub struct PafRecord {
    pub query_name: String,
    pub query_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub target_name: String,
    pub target_len: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub num_matches: usize,
    pub alignment_len: usize,
    pub mapq: u8,
    /// Optional tags (e.g., AS:i:score, NM:i:mismatches)
    pub tags: Vec<String>,
}

impl PafRecord {
    /// BLAST-style query identity: fraction of the query that matched
    /// identically. This excludes reference-padding deletions that would
    /// otherwise depress the score for short queries embedded in a larger
    /// reference window (Dragon's WFA does global alignment, so any
    /// padding bp past the actual hit becomes a D op).
    pub fn identity(&self) -> f64 {
        if self.query_len == 0 {
            0.0
        } else {
            self.num_matches as f64 / self.query_len as f64
        }
    }
}

impl fmt::Display for PafRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            self.strand,
            self.target_name,
            self.target_len,
            self.target_start,
            self.target_end,
            self.num_matches,
            self.alignment_len,
            self.mapq,
        )?;
        for tag in &self.tags {
            write!(f, "\t{}", tag)?;
        }
        Ok(())
    }
}

/// Write PAF records to a writer.
pub fn write_paf<W: Write>(writer: &mut W, records: &[PafRecord]) -> std::io::Result<()> {
    for record in records {
        writeln!(writer, "{}", record)?;
    }
    Ok(())
}
