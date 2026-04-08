/// Per-query prevalence summary output.
///
/// Groups alignment hits by species/lineage and computes:
/// - Species-level prevalence (fraction of hit genomes per species)
/// - Number of unique sequence variants (distinct identity values)
/// - ANI distribution statistics (mean, min, max identity along hits)
/// - Hit count per species

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
    pub min_identity: f64,
    pub max_identity: f64,
    pub unique_variants: usize,
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

            // Count unique identity bins (rounded to 0.1%) as proxy for sequence variants
            let unique_variants: std::collections::HashSet<u32> = identities
                .iter()
                .map(|&id| (id * 1000.0) as u32)
                .collect();

            let mean_id = identities.iter().sum::<f64>() / identities.len() as f64;
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
                min_identity: min_id,
                max_identity: max_id,
                unique_variants: unique_variants.len(),
            }
        })
        .collect();

    // Sort by hit count descending
    summaries.sort_by(|a, b| b.hit_count.cmp(&a.hit_count));
    summaries
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

/// Infer species from a genome name.
///
/// Heuristic: split on '_' or '.' and take the first two tokens as "Genus_species".
/// Falls back to the full name if fewer than 2 tokens.
fn infer_species(genome_name: &str) -> String {
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
            min_identity: 0.90,
            max_identity: 0.99,
            unique_variants: 4,
        }];

        let mut buf = Vec::new();
        write_summary(&mut buf, "query1", &summaries).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("Escherichia_coli"));
        assert!(output.contains("query1"));
    }
}
