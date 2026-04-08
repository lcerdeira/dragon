/// Signal-level nanopore search module for Dragon.
///
/// This module enables searching raw nanopore current signals directly against
/// a genome index, bypassing basecalling. The approach:
///
/// 1. **Pore Model**: Maps DNA k-mers to expected pA current levels (R10.4.1 approximate model)
/// 2. **Discretization**: Converts continuous pA values to a finite alphabet (default 16 levels)
/// 3. **Signal Index**: Converts genome sequences to expected signal, discretizes,
///    and builds an FM-index in signal space
/// 4. **Signal Search**: Normalizes query signal, discretizes, and performs FM-index
///    backward search for signal pattern matching
///
/// # Example workflow
///
/// ```no_run
/// use dragon::signal;
///
/// // Build signal index from genomes
/// let config = signal::index::SignalIndexConfig::default();
/// signal::index::build_signal_index(
///     &std::path::Path::new("genomes/"),
///     &std::path::Path::new("signal_index/"),
///     &config,
/// ).unwrap();
///
/// // Search signal reads against the index
/// let search_config = signal::search::SignalSearchConfig {
///     index_dir: std::path::Path::new("signal_index/").into(),
///     signal_kmer_size: 10,
///     ..Default::default()
/// };
/// let results = signal::search::search_signal_file(
///     &std::path::Path::new("reads.tsv"),
///     &search_config,
/// ).unwrap();
/// ```

pub mod discretize;
pub mod index;
pub mod io;
pub mod model;
pub mod search;

// Re-export key types for convenience.
pub use discretize::SignalAlphabet;
pub use index::{SignalIndexConfig, SignalIndexMetadata};
pub use io::SignalRead;
pub use model::PoreModel;
pub use search::{SignalHit, SignalSearchConfig, SignalSearchResult};
