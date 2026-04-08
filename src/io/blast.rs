/// BLAST-tabular (outfmt 6) output writer.

use std::io::Write;

use super::paf::PafRecord;

/// Write records in BLAST tabular format (outfmt 6).
/// Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
pub fn write_blast_tabular<W: Write>(writer: &mut W, records: &[PafRecord]) -> std::io::Result<()> {
    for r in records {
        let identity = r.identity() * 100.0;
        let mismatches = r.alignment_len.saturating_sub(r.num_matches);
        // Gap opens approximation (not tracked precisely in PAF)
        let gap_opens = 0;
        // E-value and bit-score are placeholders (would need database size for real calculation)
        let evalue = 0.0f64;
        let bitscore = r.num_matches as f64 * 2.0; // rough approximation

        writeln!(
            writer,
            "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
            r.query_name,
            r.target_name,
            identity,
            r.alignment_len,
            mismatches,
            gap_opens,
            r.query_start + 1, // 1-based
            r.query_end,
            r.target_start + 1, // 1-based
            r.target_end,
            evalue,
            bitscore,
        )?;
    }
    Ok(())
}
