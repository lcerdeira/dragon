/// Per-query prevalence summary and surveillance output.
///
/// Groups alignment hits by species/lineage and computes:
/// - Species-level prevalence (fraction of hit genomes per species)
/// - Number of unique sequence variants (distinct identity values)
/// - ANI distribution statistics (mean, median, stdev, min, max)
/// - Hit count and genome count per species
/// - Lineage-level rollup (genus, family if parseable)
///
/// Output formats: TSV (default) and JSON (for programmatic consumption).

use std::collections::HashMap;
use std::io::Write;

use super::paf::PafRecord;

/// A species-level summary for one query.
#[derive(Clone, Debug)]
pub struct SpeciesSummary {
    pub species: String,
    pub hit_count: usize,
    pub genome_count: usize,
    pub prevalence: f64,
    pub mean_identity: f64,
    pub median_identity: f64,
    pub stdev_identity: f64,
    pub min_identity: f64,
    pub max_identity: f64,
    pub unique_variants: usize,
    /// Individual genome hits with identities (for detailed output).
    pub genome_identities: Vec<(String, f64)>,
}

/// A multi-query surveillance report across all queries.
#[derive(Clone, Debug)]
pub struct SurveillanceReport {
    /// Per-query summaries.
    pub query_summaries: Vec<QuerySummary>,
    /// Aggregated species prevalence across all queries.
    pub aggregate_species: Vec<AggregateSpecies>,
    /// Total queries processed.
    pub total_queries: usize,
    /// Total genomes in the database.
    pub total_genomes: usize,
}

/// Summary for a single query.
#[derive(Clone, Debug)]
pub struct QuerySummary {
    pub query_name: String,
    pub query_len: usize,
    pub total_hits: usize,
    pub species: Vec<SpeciesSummary>,
}

/// Aggregated species presence across multiple queries.
#[derive(Clone, Debug)]
pub struct AggregateSpecies {
    pub species: String,
    pub queries_with_hits: usize,
    pub total_hits: usize,
    pub total_genomes: usize,
    pub mean_identity_across_queries: f64,
    pub prevalence: f64,
}

/// Compute species-level prevalence summaries for a set of query results.
///
/// Species is inferred from genome names using a configurable delimiter.
/// Default: split on '_' and take the first two tokens as "Genus_species".
pub fn summarise_hits(
    records: &[PafRecord],
    total_genomes_in_db: usize,
) -> Vec<SpeciesSummary> {
    if records.is_empty() {
        return Vec::new();
    }

    // Group by inferred species
    let mut species_hits: HashMap<String, Vec<&PafRecord>> = HashMap::new();

    for record in records {
        let species = infer_species(&record.target_name);
        species_hits.entry(species).or_default().push(record);
    }

    let mut summaries: Vec<SpeciesSummary> = species_hits
        .into_iter()
        .map(|(species, hits)| {
            let identities: Vec<f64> = hits.iter().map(|r| r.identity()).collect();

            // Count unique genomes
            let unique_genomes: std::collections::HashSet<&str> =
                hits.iter().map(|r| r.target_name.as_str()).collect();

            // Genome-level identities for detailed output
            let genome_identities: Vec<(String, f64)> = hits
                .iter()
                .map(|r| (r.target_name.clone(), r.identity()))
                .collect();

            // Count unique identity bins (rounded to 0.1%) as proxy for sequence variants
            let unique_variants: std::collections::HashSet<u32> = identities
                .iter()
                .map(|&id| (id * 1000.0) as u32)
                .collect();

            let mean_id = identities.iter().sum::<f64>() / identities.len() as f64;
            let median_id = median(&identities);
            let stdev_id = stdev(&identities, mean_id);
            let min_id = identities.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_id = identities.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

            let prevalence = if total_genomes_in_db > 0 {
                unique_genomes.len() as f64 / total_genomes_in_db as f64
            } else {
                0.0
            };

            SpeciesSummary {
                species,
                hit_count: hits.len(),
                genome_count: unique_genomes.len(),
                prevalence,
                mean_identity: mean_id,
                median_identity: median_id,
                stdev_identity: stdev_id,
                min_identity: min_id,
                max_identity: max_id,
                unique_variants: unique_variants.len(),
                genome_identities,
            }
        })
        .collect();

    // Sort by hit count descending
    summaries.sort_by(|a, b| b.hit_count.cmp(&a.hit_count));
    summaries
}

/// Build a surveillance report from multiple query results.
pub fn build_surveillance_report(
    query_results: &[(String, usize, Vec<PafRecord>)], // (name, len, records)
    total_genomes: usize,
) -> SurveillanceReport {
    let mut query_summaries = Vec::new();
    let mut species_agg: HashMap<String, (usize, usize, usize, Vec<f64>)> = HashMap::new();

    for (name, len, records) in query_results {
        let summaries = summarise_hits(records, total_genomes);

        for s in &summaries {
            let entry = species_agg
                .entry(s.species.clone())
                .or_insert((0, 0, 0, Vec::new()));
            entry.0 += 1; // queries_with_hits
            entry.1 += s.hit_count;
            entry.2 += s.genome_count;
            entry.3.push(s.mean_identity);
        }

        query_summaries.push(QuerySummary {
            query_name: name.clone(),
            query_len: *len,
            total_hits: records.len(),
            species: summaries,
        });
    }

    let aggregate_species: Vec<AggregateSpecies> = {
        let mut agg: Vec<_> = species_agg
            .into_iter()
            .map(|(species, (queries, hits, genomes, identities))| {
                let mean_id = if identities.is_empty() {
                    0.0
                } else {
                    identities.iter().sum::<f64>() / identities.len() as f64
                };
                AggregateSpecies {
                    species,
                    queries_with_hits: queries,
                    total_hits: hits,
                    total_genomes: genomes,
                    mean_identity_across_queries: mean_id,
                    prevalence: if total_genomes > 0 {
                        genomes as f64 / total_genomes as f64
                    } else {
                        0.0
                    },
                }
            })
            .collect();
        agg.sort_by(|a, b| b.total_hits.cmp(&a.total_hits));
        agg
    };

    SurveillanceReport {
        query_summaries,
        aggregate_species,
        total_queries: query_results.len(),
        total_genomes,
    }
}

/// Write prevalence summary as TSV.
pub fn write_summary<W: Write>(
    writer: &mut W,
    query_name: &str,
    summaries: &[SpeciesSummary],
) -> std::io::Result<()> {
    writeln!(
        writer,
        "#query\tspecies\thit_count\tgenome_count\tprevalence\tmean_identity\tmin_identity\tmax_identity\tunique_variants"
    )?;

    for s in summaries {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{}",
            query_name,
            s.species,
            s.hit_count,
            s.genome_count,
            s.prevalence,
            s.mean_identity,
            s.min_identity,
            s.max_identity,
            s.unique_variants,
        )?;
    }

    Ok(())
}

/// Write a full surveillance report as TSV.
pub fn write_surveillance_tsv<W: Write>(
    writer: &mut W,
    report: &SurveillanceReport,
) -> std::io::Result<()> {
    // Section 1: Per-query species breakdown
    writeln!(writer, "## Per-query species breakdown")?;
    writeln!(
        writer,
        "#query\tquery_len\tspecies\thit_count\tgenome_count\tprevalence\tmean_ANI\tmedian_ANI\tstdev_ANI\tmin_ANI\tmax_ANI\tunique_variants"
    )?;

    for qs in &report.query_summaries {
        if qs.species.is_empty() {
            writeln!(writer, "{}\t{}\tNO_HITS\t0\t0\t0\t0\t0\t0\t0\t0\t0", qs.query_name, qs.query_len)?;
            continue;
        }
        for s in &qs.species {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}",
                qs.query_name,
                qs.query_len,
                s.species,
                s.hit_count,
                s.genome_count,
                s.prevalence,
                s.mean_identity,
                s.median_identity,
                s.stdev_identity,
                s.min_identity,
                s.max_identity,
                s.unique_variants,
            )?;
        }
    }

    // Section 2: Aggregate species prevalence
    writeln!(writer)?;
    writeln!(writer, "## Aggregate species prevalence across {} queries ({} genomes in database)",
        report.total_queries, report.total_genomes)?;
    writeln!(
        writer,
        "#species\tqueries_with_hits\ttotal_hits\ttotal_genomes\tprevalence\tmean_ANI"
    )?;

    for a in &report.aggregate_species {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.4}",
            a.species,
            a.queries_with_hits,
            a.total_hits,
            a.total_genomes,
            a.prevalence,
            a.mean_identity_across_queries,
        )?;
    }

    Ok(())
}

/// Write a surveillance report as JSON.
pub fn write_surveillance_json<W: Write>(
    writer: &mut W,
    report: &SurveillanceReport,
) -> std::io::Result<()> {
    writeln!(writer, "{{")?;
    writeln!(writer, "  \"total_queries\": {},", report.total_queries)?;
    writeln!(writer, "  \"total_genomes_in_db\": {},", report.total_genomes)?;

    // Aggregate
    writeln!(writer, "  \"aggregate_species\": [")?;
    for (i, a) in report.aggregate_species.iter().enumerate() {
        let comma = if i + 1 < report.aggregate_species.len() { "," } else { "" };
        writeln!(writer, "    {{")?;
        writeln!(writer, "      \"species\": \"{}\",", a.species)?;
        writeln!(writer, "      \"queries_with_hits\": {},", a.queries_with_hits)?;
        writeln!(writer, "      \"total_hits\": {},", a.total_hits)?;
        writeln!(writer, "      \"total_genomes\": {},", a.total_genomes)?;
        writeln!(writer, "      \"prevalence\": {:.6},", a.prevalence)?;
        writeln!(writer, "      \"mean_ANI\": {:.4}", a.mean_identity_across_queries)?;
        writeln!(writer, "    }}{}", comma)?;
    }
    writeln!(writer, "  ],")?;

    // Per-query
    writeln!(writer, "  \"queries\": [")?;
    for (qi, qs) in report.query_summaries.iter().enumerate() {
        let qcomma = if qi + 1 < report.query_summaries.len() { "," } else { "" };
        writeln!(writer, "    {{")?;
        writeln!(writer, "      \"query_name\": \"{}\",", qs.query_name)?;
        writeln!(writer, "      \"query_len\": {},", qs.query_len)?;
        writeln!(writer, "      \"total_hits\": {},", qs.total_hits)?;
        writeln!(writer, "      \"species\": [")?;
        for (si, s) in qs.species.iter().enumerate() {
            let scomma = if si + 1 < qs.species.len() { "," } else { "" };
            writeln!(writer, "        {{")?;
            writeln!(writer, "          \"species\": \"{}\",", s.species)?;
            writeln!(writer, "          \"hit_count\": {},", s.hit_count)?;
            writeln!(writer, "          \"genome_count\": {},", s.genome_count)?;
            writeln!(writer, "          \"prevalence\": {:.6},", s.prevalence)?;
            writeln!(writer, "          \"mean_ANI\": {:.4},", s.mean_identity)?;
            writeln!(writer, "          \"median_ANI\": {:.4},", s.median_identity)?;
            writeln!(writer, "          \"stdev_ANI\": {:.4},", s.stdev_identity)?;
            writeln!(writer, "          \"min_ANI\": {:.4},", s.min_identity)?;
            writeln!(writer, "          \"max_ANI\": {:.4},", s.max_identity)?;
            writeln!(writer, "          \"unique_variants\": {},", s.unique_variants)?;
            // Genome-level details
            writeln!(writer, "          \"genomes\": [")?;
            for (gi, (gname, gid)) in s.genome_identities.iter().enumerate() {
                let gcomma = if gi + 1 < s.genome_identities.len() { "," } else { "" };
                writeln!(writer, "            {{\"name\": \"{}\", \"ANI\": {:.4}}}{}", gname, gid, gcomma)?;
            }
            writeln!(writer, "          ]")?;
            writeln!(writer, "        }}{}", scomma)?;
        }
        writeln!(writer, "      ]")?;
        writeln!(writer, "    }}{}", qcomma)?;
    }
    writeln!(writer, "  ]")?;
    writeln!(writer, "}}")?;

    Ok(())
}

/// Parse PAF records from a reader (for standalone summarize command).
pub fn parse_paf<R: std::io::BufRead>(reader: R) -> Vec<PafRecord> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let record = PafRecord {
            query_name: fields[0].to_string(),
            query_len: fields[1].parse().unwrap_or(0),
            query_start: fields[2].parse().unwrap_or(0),
            query_end: fields[3].parse().unwrap_or(0),
            strand: fields[4].chars().next().unwrap_or('+'),
            target_name: fields[5].to_string(),
            target_len: fields[6].parse().unwrap_or(0),
            target_start: fields[7].parse().unwrap_or(0),
            target_end: fields[8].parse().unwrap_or(0),
            num_matches: fields[9].parse().unwrap_or(0),
            alignment_len: fields[10].parse().unwrap_or(0),
            mapq: fields[11].parse().unwrap_or(0),
            tags: fields[12..].iter().map(|s| s.to_string()).collect(),
        };

        records.push(record);
    }

    records
}

/// Infer species from a genome name.
///
/// Heuristic: split on '_' or '.' and take the first two tokens as "Genus_species".
/// Falls back to the full name if fewer than 2 tokens.
pub fn infer_species(genome_name: &str) -> String {
    // Strip common extensions
    let name = genome_name
        .trim_end_matches(".fa")
        .trim_end_matches(".fasta")
        .trim_end_matches(".fna");

    // Try underscore split first (e.g., "Escherichia_coli_K12_MG1655")
    let parts: Vec<&str> = name.split('_').collect();
    if parts.len() >= 2 {
        return format!("{}_{}", parts[0], parts[1]);
    }

    // Try dot split
    let parts: Vec<&str> = name.split('.').collect();
    if parts.len() >= 2 {
        return format!("{}_{}", parts[0], parts[1]);
    }

    name.to_string()
}

// ---- Statistics helpers ----

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

fn stdev(values: &[f64], mean: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let variance = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(target: &str, matches: usize, align_len: usize) -> PafRecord {
        PafRecord {
            query_name: "query1".to_string(),
            query_len: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target_name: target.to_string(),
            target_len: 5_000_000,
            target_start: 0,
            target_end: 1000,
            num_matches: matches,
            alignment_len: align_len,
            mapq: 60,
            tags: vec![],
        }
    }

    #[test]
    fn test_infer_species() {
        assert_eq!(infer_species("Escherichia_coli_K12_MG1655.fa"), "Escherichia_coli");
        assert_eq!(infer_species("Klebsiella_pneumoniae_strain42"), "Klebsiella_pneumoniae");
        assert_eq!(infer_species("unknown"), "unknown");
    }

    #[test]
    fn test_summarise_hits() {
        let records = vec![
            make_record("Escherichia_coli_K12", 950, 1000),
            make_record("Escherichia_coli_O157", 920, 1000),
            make_record("Klebsiella_pneumoniae_01", 880, 1000),
        ];

        let summaries = summarise_hits(&records, 100);
        assert_eq!(summaries.len(), 2);
        assert_eq!(summaries[0].species, "Escherichia_coli");
        assert_eq!(summaries[0].hit_count, 2);
        assert_eq!(summaries[0].genome_count, 2);
        assert!((summaries[0].prevalence - 0.02).abs() < 1e-6);
    }

    #[test]
    fn test_write_summary() {
        let summaries = vec![SpeciesSummary {
            species: "Escherichia_coli".to_string(),
            hit_count: 5,
            genome_count: 3,
            prevalence: 0.001,
            mean_identity: 0.95,
            median_identity: 0.95,
            stdev_identity: 0.02,
            min_identity: 0.90,
            max_identity: 0.99,
            unique_variants: 4,
            genome_identities: vec![
                ("Ecoli_K12".into(), 0.99),
                ("Ecoli_O157".into(), 0.95),
                ("Ecoli_DH5a".into(), 0.90),
            ],
        }];

        let mut buf = Vec::new();
        write_summary(&mut buf, "query1", &summaries).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("Escherichia_coli"));
        assert!(output.contains("query1"));
    }

    #[test]
    fn test_surveillance_report() {
        let q1_records = vec![
            make_record("Escherichia_coli_K12", 950, 1000),
            make_record("Escherichia_coli_O157", 920, 1000),
        ];
        let q2_records = vec![
            make_record("Escherichia_coli_K12", 940, 1000),
            make_record("Klebsiella_pneumoniae_01", 880, 1000),
        ];

        let query_data = vec![
            ("gene_A".into(), 1000usize, q1_records),
            ("gene_B".into(), 800usize, q2_records),
        ];

        let report = build_surveillance_report(&query_data, 10000);
        assert_eq!(report.total_queries, 2);
        assert_eq!(report.aggregate_species.len(), 2);
        // E. coli should appear in both queries
        let ecoli = report.aggregate_species.iter().find(|a| a.species == "Escherichia_coli").unwrap();
        assert_eq!(ecoli.queries_with_hits, 2);
    }

    #[test]
    fn test_parse_paf() {
        let paf_line = "query1\t1000\t0\t1000\t+\tEcoli_K12\t5000000\t0\t1000\t950\t1000\t60\tAS:i:950\n";
        let records = parse_paf(std::io::BufReader::new(paf_line.as_bytes()));
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name, "query1");
        assert_eq!(records[0].target_name, "Ecoli_K12");
        assert_eq!(records[0].num_matches, 950);
    }

    #[test]
    fn test_median_stdev() {
        assert!((median(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-10);
        assert!((median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
        // Sample stdev: deviations^2 = 9+1+1+1+0+0+4+16 = 32, variance = 32/7 ≈ 4.571
        assert!((stdev(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0], 5.0) - 2.138).abs() < 0.01);
    }

    #[test]
    fn test_json_output() {
        let records = vec![
            make_record("Escherichia_coli_K12", 950, 1000),
        ];
        let query_data = vec![("query1".into(), 1000usize, records)];
        let report = build_surveillance_report(&query_data, 100);

        let mut buf = Vec::new();
        write_surveillance_json(&mut buf, &report).unwrap();
        let json = String::from_utf8(buf).unwrap();
        assert!(json.contains("\"species\": \"Escherichia_coli\""));
        assert!(json.contains("\"total_queries\": 1"));
    }
}
