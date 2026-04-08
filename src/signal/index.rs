/// Signal-level index construction.
///
/// Converts genome reference sequences into expected nanopore signal space
/// using a pore model, discretizes the signal, and builds an FM-index over
/// the discretized signal alphabet for backward search.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::ds::elias_fano::CumulativeLengthIndex;
use crate::index::fm::DragonFmIndex;
use crate::io::fasta;
use crate::signal::discretize::SignalAlphabet;
use crate::signal::model::{self, PoreModel};

/// Configuration for signal index construction.
#[derive(Clone, Debug)]
pub struct SignalIndexConfig {
    /// Signal discretization alphabet.
    pub alphabet: SignalAlphabet,
    /// Pore model to use.
    pub pore_model: PoreModel,
    /// Number of threads for parallel processing.
    pub threads: usize,
}

impl Default for SignalIndexConfig {
    fn default() -> Self {
        Self {
            alphabet: SignalAlphabet::default(),
            pore_model: model::load_default_model(),
            threads: 4,
        }
    }
}

/// Metadata for a signal index, saved alongside the index data.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SignalIndexMetadata {
    /// Dragon version that built this index.
    pub version: String,
    /// Number of genomes indexed.
    pub num_genomes: usize,
    /// Genome names in order.
    pub genome_names: Vec<String>,
    /// Total discretized signal length.
    pub total_signal_length: u64,
    /// Pore model name used.
    pub pore_model_name: String,
    /// Pore model k-mer size.
    pub pore_model_kmer_size: usize,
    /// Number of discrete levels in alphabet.
    pub num_levels: u8,
    /// Signal alphabet min value.
    pub alphabet_min: f32,
    /// Signal alphabet max value.
    pub alphabet_max: f32,
    /// Optional learned discretization boundaries.
    /// When present, the alphabet uses non-uniform bins via binary search.
    #[serde(default)]
    pub alphabet_boundaries: Option<Vec<f32>>,
}

/// Build a signal-level FM-index from a directory of genome FASTA files.
///
/// Pipeline:
/// 1. Read each genome FASTA
/// 2. Convert DNA sequences to expected pA signal via the pore model
/// 3. Normalize and discretize the expected signal
/// 4. Concatenate all discretized signals with sentinel separators
/// 5. Build suffix array / FM-index over the concatenated text
/// 6. Save index + metadata to output_dir
pub fn build_signal_index(
    genomes_dir: &Path,
    output_dir: &Path,
    config: &SignalIndexConfig,
) -> Result<()> {
    std::fs::create_dir_all(output_dir)?;

    // Step 1: Find all genome FASTA files
    let fasta_files = fasta::list_fasta_files(genomes_dir)?;
    if fasta_files.is_empty() {
        anyhow::bail!("No FASTA files found in {:?}", genomes_dir);
    }
    log::info!("Found {} genome FASTA files", fasta_files.len());

    // Step 2: Convert each genome to discretized signal and concatenate
    let sentinel = config.alphabet.num_levels; // sentinel is outside the alphabet
    let mut concatenated: Vec<u8> = Vec::new();
    let mut genome_signal_lengths: Vec<u64> = Vec::new();
    let mut genome_names: Vec<String> = Vec::new();

    for fasta_path in &fasta_files {
        let genome_name = fasta_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        log::info!("Processing genome: {}", genome_name);

        let sequences = fasta::read_sequences(fasta_path)
            .with_context(|| format!("Failed to read {:?}", fasta_path))?;

        let mut genome_signal: Vec<u8> = Vec::new();

        for seq_record in &sequences {
            // Convert DNA to expected pA signal
            let expected_signal = config.pore_model.sequence_to_expected_signal(&seq_record.seq);

            if expected_signal.is_empty() {
                continue;
            }

            // Normalize the expected signal (even though it's already "ideal",
            // normalization ensures consistency with query signal processing)
            let normalized = normalize_expected_signal(&expected_signal, &config.pore_model);

            // Discretize
            let discretized =
                crate::signal::discretize::discretize_normalized(&normalized, &config.alphabet);

            genome_signal.extend_from_slice(&discretized);

            // Also process reverse complement
            let rc_seq: Vec<u8> = seq_record
                .seq
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' | b'a' => b'T',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    b'T' | b't' => b'A',
                    _ => b'N',
                })
                .collect();

            let rc_signal = config.pore_model.sequence_to_expected_signal(&rc_seq);
            if !rc_signal.is_empty() {
                let rc_norm = normalize_expected_signal(&rc_signal, &config.pore_model);
                let rc_disc =
                    crate::signal::discretize::discretize_normalized(&rc_norm, &config.alphabet);
                // Add separator between forward and reverse complement
                genome_signal.push(sentinel);
                genome_signal.extend_from_slice(&rc_disc);
            }
        }

        let sig_len = genome_signal.len() as u64;
        genome_signal_lengths.push(sig_len);
        genome_names.push(genome_name);

        // Append to concatenated text with sentinel separator
        concatenated.extend_from_slice(&genome_signal);
        concatenated.push(sentinel);
    }

    log::info!(
        "Total discretized signal length: {} symbols across {} genomes",
        concatenated.len(),
        genome_names.len()
    );

    // Step 3: Build suffix array
    log::info!("Building suffix array over signal text...");
    let suffix_array = build_signal_suffix_array(&concatenated);

    // Step 4: Build cumulative length index for genome mapping
    let cum_lengths = CumulativeLengthIndex::from_lengths(&genome_signal_lengths);

    // Step 5: Create the FM-index structure
    let fm_index = DragonFmIndex {
        text: concatenated,
        suffix_array,
        cumulative_lengths: cum_lengths,
    };

    // Step 6: Save to disk
    let index_path = output_dir.join("signal_fm_index.bin");
    crate::util::mmap::write_bincode(&index_path, &SignalFmIndexSerializable::from_index(&fm_index))?;
    log::info!("Signal FM-index saved to {:?}", index_path);

    // Step 7: Save metadata
    let metadata = SignalIndexMetadata {
        version: env!("CARGO_PKG_VERSION").to_string(),
        num_genomes: genome_names.len(),
        genome_names,
        total_signal_length: fm_index.text.len() as u64,
        pore_model_name: config.pore_model.name.clone(),
        pore_model_kmer_size: config.pore_model.kmer_size,
        num_levels: config.alphabet.num_levels,
        alphabet_min: config.alphabet.min_val,
        alphabet_max: config.alphabet.max_val,
        alphabet_boundaries: config.alphabet.boundaries.clone(),
    };

    let meta_path = output_dir.join("signal_metadata.json");
    let json = serde_json::to_string_pretty(&metadata)?;
    std::fs::write(&meta_path, json)?;
    log::info!("Signal index metadata saved to {:?}", meta_path);

    // Step 8: Save the pore model and alphabet config
    let model_path = output_dir.join("signal_pore_model.bin");
    crate::util::mmap::write_bincode(&model_path, &config.pore_model)?;

    let alpha_path = output_dir.join("signal_alphabet.bin");
    crate::util::mmap::write_bincode(&alpha_path, &config.alphabet)?;

    log::info!("Signal index construction complete");
    Ok(())
}

/// Load a signal FM-index from disk.
pub fn load_signal_index(index_dir: &Path) -> Result<DragonFmIndex> {
    let index_path = index_dir.join("signal_fm_index.bin");
    let serialized: SignalFmIndexSerializable = crate::util::mmap::read_bincode(&index_path)?;
    Ok(serialized.to_index())
}

/// Load signal index metadata.
pub fn load_signal_metadata(index_dir: &Path) -> Result<SignalIndexMetadata> {
    let meta_path = index_dir.join("signal_metadata.json");
    let json = std::fs::read_to_string(&meta_path)?;
    let metadata: SignalIndexMetadata = serde_json::from_str(&json)?;
    Ok(metadata)
}

/// Load the pore model from a signal index directory.
pub fn load_signal_pore_model(index_dir: &Path) -> Result<PoreModel> {
    let model_path = index_dir.join("signal_pore_model.bin");
    crate::util::mmap::read_bincode(&model_path)
}

/// Load the signal alphabet from a signal index directory.
pub fn load_signal_alphabet(index_dir: &Path) -> Result<SignalAlphabet> {
    let alpha_path = index_dir.join("signal_alphabet.bin");
    crate::util::mmap::read_bincode(&alpha_path)
}

/// Normalize expected signal from pore model to zero-mean, unit-variance.
///
/// For expected signal (which has no noise), we use the model's own statistics
/// rather than per-read median-MAD normalization.
fn normalize_expected_signal(signal: &[f32], _model: &PoreModel) -> Vec<f32> {
    if signal.is_empty() {
        return Vec::new();
    }

    let n = signal.len() as f32;
    let mean: f32 = signal.iter().sum::<f32>() / n;
    let var: f32 = signal.iter().map(|&x| (x - mean) * (x - mean)).sum::<f32>() / n;
    let std_dev = var.sqrt();

    if std_dev < 1e-10 {
        return vec![0.0; signal.len()];
    }

    signal.iter().map(|&x| (x - mean) / std_dev).collect()
}

/// Build a suffix array for the signal text.
/// Uses the same simple approach as the DNA FM-index.
fn build_signal_suffix_array(text: &[u8]) -> Vec<usize> {
    let n = text.len();
    if n == 0 {
        return Vec::new();
    }
    let mut sa: Vec<usize> = (0..n).collect();
    sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
    sa
}

/// Serializable form of the signal FM-index.
#[derive(Serialize, Deserialize)]
struct SignalFmIndexSerializable {
    text: Vec<u8>,
    suffix_array: Vec<usize>,
    cumulative_lengths: CumulativeLengthIndex,
}

impl SignalFmIndexSerializable {
    fn from_index(index: &DragonFmIndex) -> Self {
        Self {
            text: index.text.clone(),
            suffix_array: index.suffix_array.clone(),
            cumulative_lengths: index.cumulative_lengths.clone(),
        }
    }

    fn to_index(self) -> DragonFmIndex {
        DragonFmIndex {
            text: self.text,
            suffix_array: self.suffix_array,
            cumulative_lengths: self.cumulative_lengths,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn create_test_genome_dir() -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        let genome_path = dir.path().join("genome1.fa");
        let mut f = std::fs::File::create(&genome_path).unwrap();
        writeln!(f, ">chr1").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        writeln!(f, ">chr2").unwrap();
        writeln!(f, "TTGGCCAATTGGCCAATTGGCCAATTGGCCAA").unwrap();
        dir
    }

    #[test]
    fn test_build_signal_index() {
        let genome_dir = create_test_genome_dir();
        let output_dir = tempfile::tempdir().unwrap();

        let config = SignalIndexConfig::default();
        build_signal_index(genome_dir.path(), output_dir.path(), &config).unwrap();

        // Verify files were created
        assert!(output_dir.path().join("signal_fm_index.bin").exists());
        assert!(output_dir.path().join("signal_metadata.json").exists());
        assert!(output_dir.path().join("signal_pore_model.bin").exists());
        assert!(output_dir.path().join("signal_alphabet.bin").exists());
    }

    #[test]
    fn test_load_signal_index() {
        let genome_dir = create_test_genome_dir();
        let output_dir = tempfile::tempdir().unwrap();

        let config = SignalIndexConfig::default();
        build_signal_index(genome_dir.path(), output_dir.path(), &config).unwrap();

        let index = load_signal_index(output_dir.path()).unwrap();
        assert!(!index.text.is_empty());
        assert_eq!(index.text.len(), index.suffix_array.len());

        let metadata = load_signal_metadata(output_dir.path()).unwrap();
        assert_eq!(metadata.num_genomes, 1);
        assert_eq!(metadata.genome_names, vec!["genome1"]);
    }

    #[test]
    fn test_normalize_expected_signal() {
        let model = model::load_default_model();
        let signal = vec![80.0, 90.0, 100.0, 85.0, 95.0];
        let norm = normalize_expected_signal(&signal, &model);
        assert_eq!(norm.len(), 5);

        // Mean of normalized should be ~0
        let mean: f32 = norm.iter().sum::<f32>() / norm.len() as f32;
        assert!(mean.abs() < 1e-5);
    }

    #[test]
    fn test_signal_index_searchable() {
        let genome_dir = create_test_genome_dir();
        let output_dir = tempfile::tempdir().unwrap();

        let config = SignalIndexConfig::default();
        build_signal_index(genome_dir.path(), output_dir.path(), &config).unwrap();

        let index = load_signal_index(output_dir.path()).unwrap();

        // The index should be searchable: look for a short pattern from the text itself
        if index.text.len() >= 3 {
            let pattern = &index.text[0..3];
            let hits = index.search(pattern);
            assert!(!hits.is_empty(), "Should find at least one hit for pattern from the text");
        }
    }
}
