/// BLAST-tabular output writer (outfmt 6 / outfmt 7).
///
/// Produces output a BLAST user can read directly: the same 12 columns as
/// `blastn -outfmt 6`, with a real percent identity, mismatch/gap counts
/// recovered from the CIGAR, and a Karlin–Altschul E-value / bit-score.
///
/// `blast6` is the bare, machine-parseable form (no header). `blast7` adds a
/// `# Fields:` comment header like `blastn -outfmt 7`, for humans.

use std::io::Write;

use super::paf::PafRecord;

/// Karlin–Altschul parameters for ungapped nucleotide alignment with the
/// common BLASTN scoring (match reward +1, mismatch penalty −2). These are the
/// standard published values (see NCBI BLAST statistics). We re-score the
/// reported alignment under this model so the E-value is comparable to what a
/// BLASTN user expects, rather than exposing Dragon's raw WFA score.
const KA_LAMBDA: f64 = 1.374;
const KA_K: f64 = 0.711;
const MATCH_REWARD: f64 = 1.0;
const MISMATCH_PENALTY: f64 = 2.0;

/// Column header in `blastn -outfmt 7` style.
pub const BLAST_FIELDS: &str = "# Fields: query id\tsubject id\t% identity\t\
alignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\t\
evalue\tbit score";

/// CIGAR-derived alignment stats, with terminal gaps trimmed.
struct AlnStats {
    aln_len: usize,
    mismatches: usize,
    gap_opens: usize,
}

/// Parse a minimap2-style CIGAR (`cg:Z:` value, ops `=`/`X`/`M`/`I`/`D`),
/// trimming leading and trailing indel runs (reference padding from Dragon's
/// global WFA, which would otherwise inflate length and gap counts for a short
/// query embedded in a larger reference window). Returns alignment length over
/// the trimmed span, substitution count, and internal gap-open count.
fn parse_cigar(cigar: &str) -> Option<AlnStats> {
    let mut ops: Vec<(usize, char)> = Vec::new();
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            let len: usize = num.parse().ok()?;
            num.clear();
            ops.push((len, ch));
        }
    }
    if ops.is_empty() {
        return None;
    }
    // Trim leading / trailing indel runs (terminal gaps = reference padding).
    let is_indel = |c: char| c == 'I' || c == 'D';
    let start = ops.iter().position(|&(_, c)| !is_indel(c))?;
    let end = ops.iter().rposition(|&(_, c)| !is_indel(c))?;
    let core = &ops[start..=end];

    let mut aln_len = 0;
    let mut mismatches = 0;
    let mut gap_opens = 0;
    for &(len, op) in core {
        match op {
            '=' | 'M' => aln_len += len,
            'X' => {
                aln_len += len;
                mismatches += len;
            }
            'I' | 'D' => {
                aln_len += len;
                gap_opens += 1;
            }
            _ => {}
        }
    }
    Some(AlnStats {
        aln_len,
        mismatches,
        gap_opens,
    })
}

fn tag_str<'a>(r: &'a PafRecord, prefix: &str) -> Option<&'a str> {
    r.tags.iter().find_map(|t| t.strip_prefix(prefix))
}

/// Bit-score and E-value under the BLASTN +1/−2 ungapped model.
///
/// raw score S = matches·reward − mismatches·penalty; bit-score
/// S' = (λ·S − ln K) / ln 2; E = m·n·2^(−S'), where m is the query length and
/// n is the effective database size (total searched bases). For strong, long
/// alignments S' is large and E underflows to 0.0 — exactly as BLAST reports.
fn evalue_bitscore(matches: usize, mismatches: usize, qlen: usize, db_size: u64) -> (f64, f64) {
    let raw = matches as f64 * MATCH_REWARD - mismatches as f64 * MISMATCH_PENALTY;
    if raw <= 0.0 {
        return (f64::INFINITY, 0.0);
    }
    let bitscore = (KA_LAMBDA * raw - KA_K.ln()) / std::f64::consts::LN_2;
    let search_space = qlen.max(1) as f64 * db_size.max(1) as f64;
    let evalue = search_space * 2f64.powf(-bitscore);
    (evalue, bitscore)
}

/// Write records in BLAST tabular format.
///
/// Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart
/// send evalue bitscore. `db_size` is the effective database size (sum of
/// `total_unitig_bases` across the searched index and shards), used for the
/// E-value. `header` emits the outfmt-7 `# Fields:` line.
pub fn write_blast_tabular<W: Write>(
    writer: &mut W,
    records: &[PafRecord],
    db_size: u64,
    header: bool,
) -> std::io::Result<()> {
    if header {
        writeln!(writer, "{}", BLAST_FIELDS)?;
    }
    for r in records {
        // Prefer CIGAR-derived stats (terminal-gap-trimmed); fall back to PAF fields.
        let (aln_len, mismatches, gap_opens) = match tag_str(r, "cg:Z:").and_then(parse_cigar) {
            Some(s) => (s.aln_len, s.mismatches, s.gap_opens),
            None => {
                // No CIGAR: NM:i: is the edit distance; span is the aligned query window.
                let nm = tag_str(r, "NM:i:")
                    .and_then(|v| v.parse::<usize>().ok())
                    .unwrap_or_else(|| r.alignment_len.saturating_sub(r.num_matches));
                let span = r.query_end.saturating_sub(r.query_start).max(r.num_matches);
                (span, nm, 0)
            }
        };

        // BLAST percent identity: matched bases over the alignment length, so
        // pident x length reconciles with the match count (BLAST column
        // semantics). Clamped to 100 in case terminal-gap trimming leaves
        // aln_len slightly below the matched span.
        let pident = if aln_len > 0 {
            (r.num_matches as f64 / aln_len as f64 * 100.0).min(100.0)
        } else {
            r.identity() * 100.0
        };

        let (evalue, bitscore) = evalue_bitscore(r.num_matches, mismatches, r.query_len, db_size);

        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
            r.query_name,
            r.target_name,
            pident,
            aln_len,
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
