use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "dragon",
    about = "Dragon: resource-efficient sequence alignment against millions of prokaryotic genomes",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a Dragon index from a directory of genome FASTA files
    Index {
        /// Directory containing genome FASTA files (.fa, .fasta, .fna)
        #[arg(short, long)]
        input: PathBuf,

        /// Output directory for the Dragon index
        #[arg(short, long)]
        output: PathBuf,

        /// K-mer size for the de Bruijn graph (default: 31)
        #[arg(short, long, default_value = "31")]
        kmer_size: usize,

        /// Number of threads for parallel processing
        #[arg(short = 'j', long, default_value = "4")]
        threads: usize,

        /// Use external-memory construction to reduce peak RAM (e.g., 8 GB instead of 64 GB)
        #[arg(long)]
        low_memory: bool,

        /// Maximum RAM in GB for low-memory mode (default: 8)
        #[arg(long, default_value = "8.0")]
        max_ram: f64,
    },

    /// Search query sequences against a Dragon index
    Search {
        /// Path to the Dragon index directory
        #[arg(short, long)]
        index: PathBuf,

        /// Query FASTA/FASTQ file
        #[arg(short, long)]
        query: PathBuf,

        /// Output file (PAF format by default, use --format for others)
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Output format: paf, blast6, summary, or gfa
        #[arg(short, long, default_value = "paf")]
        format: String,

        /// Number of threads
        #[arg(short = 'j', long, default_value = "4")]
        threads: usize,

        /// Maximum RAM usage in GB
        #[arg(long, default_value = "4.0")]
        max_ram: f64,

        /// Minimum seed length
        #[arg(long, default_value = "15")]
        min_seed_len: usize,

        /// Maximum seed frequency (skip seeds more frequent than this)
        #[arg(long, default_value = "10000")]
        max_seed_freq: usize,

        /// Minimum chain score to report
        #[arg(long, default_value = "50")]
        min_chain_score: f64,

        /// Maximum number of target genomes to report per query
        #[arg(long, default_value = "10")]
        max_target_seqs: usize,

        /// Minimum alignment identity to report (0.0-1.0)
        #[arg(long, default_value = "0.7")]
        min_identity: f64,

        /// Minimum query coverage to report (0.0-1.0)
        #[arg(long, default_value = "0.3")]
        min_query_coverage: f64,

        /// Minimum score ratio: only keep hits with score >= ratio * best_score (0.0-1.0)
        #[arg(long, default_value = "0.1")]
        min_score_ratio: f64,

        /// Hardware profile: laptop (≤8 GB, conservative) or workstation (full resources)
        #[arg(long, default_value = "workstation")]
        profile: String,

        /// GFA context radius: number of unitig steps around each hit (for --format gfa)
        #[arg(long, default_value = "5")]
        gfa_radius: usize,

        /// Disable ML seed scoring (use raw match_len as anchor weight)
        #[arg(long)]
        no_ml: bool,

        /// Path to custom ML seed scorer weights (JSON array of 7 floats)
        #[arg(long)]
        ml_weights: Option<String>,

        /// Dump all seeds with features to a TSV file (for ML training)
        #[arg(long)]
        dump_seeds: Option<String>,
    },

    /// Display information about a Dragon index
    Info {
        /// Path to the Dragon index directory
        #[arg(short, long)]
        index: PathBuf,
    },

    /// Build a signal-level index from genome FASTA files for nanopore signal search
    SignalIndex {
        /// Directory containing genome FASTA files (.fa, .fasta, .fna)
        #[arg(short, long)]
        input: PathBuf,

        /// Output directory for the signal index
        #[arg(short, long)]
        output: PathBuf,

        /// Number of discrete signal levels (alphabet size, default: 16)
        #[arg(long, default_value = "16")]
        num_levels: u8,

        /// Number of threads for parallel processing
        #[arg(short = 'j', long, default_value = "4")]
        threads: usize,

        /// Path to learned signal discretization boundaries (JSON array of floats)
        #[arg(long)]
        signal_boundaries: Option<String>,

        /// Path to a custom pore model (JSON with kmer_size, levels, name)
        #[arg(long)]
        pore_model: Option<String>,
    },

    /// Download genomes or a pre-built Dragon index
    ///
    /// Databases:
    ///   gtdb-r220, allthebacteria-v2, refseq-bacteria  — download pre-built index
    ///   allthebacteria            — download genomes from EBI AllTheBacteria
    ///   refseq                    — download genomes from NCBI RefSeq (bacteria)
    ///   refseq-representative     — download only representative RefSeq genomes
    ///   http://...                — custom URL to a pre-built index tarball
    Download {
        /// Database name or URL (see --help for full list)
        #[arg(short, long)]
        database: String,

        /// Output directory for downloaded genomes or index
        #[arg(short, long)]
        output: PathBuf,

        /// Use external-memory construction if building locally (reduces RAM to ~8 GB)
        #[arg(long)]
        low_memory: bool,

        /// K-mer size for index build (used with allthebacteria/refseq genome downloads)
        #[arg(short, long, default_value = "31")]
        kmer_size: usize,

        /// Number of threads (used with allthebacteria/refseq genome downloads)
        #[arg(short = 'j', long, default_value = "4")]
        threads: usize,

        /// Maximum RAM in GB for low-memory mode (default: 8)
        #[arg(long, default_value = "8.0")]
        max_ram: f64,

        /// Parallel download connections (used with refseq genome downloads)
        #[arg(long, default_value = "8")]
        parallel_downloads: usize,
    },

    /// Search nanopore signal reads against a signal-level index
    SignalSearch {
        /// Path to the signal index directory
        #[arg(short, long)]
        index: PathBuf,

        /// Signal file (TSV, CSV, or SLOW5 format)
        #[arg(short, long)]
        query: PathBuf,

        /// Output file (TSV format, use "-" for stdout)
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Signal k-mer size for search (default: 10)
        #[arg(long, default_value = "10")]
        signal_kmer_size: usize,

        /// Minimum number of k-mer hits to report a genome (default: 3)
        #[arg(long, default_value = "3")]
        min_hits: usize,

        /// Maximum seed frequency (skip signal k-mers more frequent than this)
        #[arg(long, default_value = "10000")]
        max_seed_freq: usize,

        /// Maximum number of results per query read
        #[arg(long, default_value = "50")]
        max_results: usize,

        /// Number of threads
        #[arg(short = 'j', long, default_value = "4")]
        threads: usize,
    },
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Index {
            input,
            output,
            kmer_size,
            threads,
            low_memory,
            max_ram,
        } => {
            log::info!("Dragon index construction");
            log::info!("Input: {:?}", input);
            log::info!("Output: {:?}", output);

            let max_ram_bytes = if low_memory {
                log::info!(
                    "Low-memory mode enabled: RAM budget {:.1} GB \
                     (uses disk-based sorting for suffix array construction)",
                    max_ram
                );
                Some((max_ram * 1_073_741_824.0) as usize)
            } else {
                None
            };

            dragon::index::build_index_with_options(&input, &output, kmer_size, threads, max_ram_bytes)?;
        }

        Commands::Search {
            index,
            query,
            output,
            format,
            threads,
            max_ram,
            min_seed_len,
            max_seed_freq,
            min_chain_score,
            max_target_seqs,
            min_identity,
            min_query_coverage,
            min_score_ratio,
            profile,
            gfa_radius,
            no_ml,
            ml_weights,
            dump_seeds,
        } => {
            log::info!("Dragon search");

            // Apply hardware profile overrides
            let (effective_ram, effective_threads) = match profile.as_str() {
                "laptop" => {
                    log::info!("Using laptop profile (≤8 GB RAM, conservative threading)");
                    (max_ram.min(8.0), threads.min(4))
                }
                _ => (max_ram, threads),
            };

            let config = dragon::query::SearchConfig {
                index_dir: index.clone().into(),
                min_seed_len,
                max_seed_freq,
                min_chain_score,
                max_target_seqs,
                threads: effective_threads,
                max_ram_gb: effective_ram,
                min_identity,
                min_query_coverage,
                min_score_ratio,
                no_ml,
                ml_weights_path: ml_weights,
                dump_seeds_path: dump_seeds,
            };

            let results = dragon::query::search(&query, &config)?;

            // Write output
            let mut writer: Box<dyn std::io::Write> = if output == "-" {
                Box::new(std::io::stdout())
            } else {
                Box::new(std::fs::File::create(&output)?)
            };

            // Load path and color indices for summary/gfa modes (lazy)
            let _need_indices = matches!(format.as_str(), "summary" | "gfa");

            for result in &results {
                let mut records = result.alignments.clone();
                for record in &mut records {
                    record.query_name = result.query_name.clone();
                }

                match format.as_str() {
                    "paf" => dragon::io::paf::write_paf(&mut writer, &records)?,
                    "blast6" => dragon::io::blast::write_blast_tabular(&mut writer, &records)?,
                    "summary" => {
                        // Load metadata to get total genome count
                        let metadata = dragon::index::metadata::load_metadata(&index)?;
                        let summaries = dragon::io::summary::summarise_hits(
                            &records,
                            metadata.num_genomes,
                        );
                        dragon::io::summary::write_summary(
                            &mut writer,
                            &result.query_name,
                            &summaries,
                        )?;
                    }
                    "gfa" => {
                        let path_index = dragon::index::paths::load_path_index(&index)?;
                        let color_index = dragon::index::color::load_color_index(&index)?;
                        let fm_index = dragon::index::fm::load_fm_index(&index)?;
                        let unitig_lengths: Vec<u64> = fm_index.cumulative_lengths.lengths().to_vec();
                        let unitigs = dragon::index::unitig::UnitigSet::from_fm_text(&fm_index.text, &unitig_lengths);
                        let subgraphs = dragon::io::graph_context::extract_hit_subgraphs(
                            &records,
                            &path_index,
                            &color_index,
                            &unitigs,
                            gfa_radius,
                        );
                        dragon::io::graph_context::write_gfa(&mut writer, &subgraphs)?;
                    }
                    _ => anyhow::bail!("Unknown output format: {}. Use paf, blast6, summary, or gfa", format),
                }
            }

            log::info!(
                "Search complete: {} queries, {} total alignments",
                results.len(),
                results.iter().map(|r| r.alignments.len()).sum::<usize>()
            );
        }

        Commands::Info { index } => {
            let metadata = dragon::index::metadata::load_metadata(&index)?;
            println!("Dragon Index Information");
            println!("========================");
            println!("Version:         {}", metadata.version);
            println!("K-mer size:      {}", metadata.kmer_size);
            println!("Genomes:         {}", metadata.num_genomes);
            println!("Unitigs:         {}", metadata.num_unitigs);
            println!("Total bases:     {}", metadata.total_unitig_bases);

            // Show disk usage
            let index_path = std::path::Path::new(&index);
            if index_path.exists() {
                let mut total_size = 0u64;
                for entry in std::fs::read_dir(index_path)? {
                    let entry = entry?;
                    total_size += entry.metadata()?.len();
                }
                println!(
                    "Index size:      {:.2} GB",
                    total_size as f64 / 1_073_741_824.0
                );
            }
        }

        Commands::Download {
            database,
            output,
            low_memory,
            kmer_size,
            threads,
            max_ram,
            parallel_downloads,
        } => {
            log::info!("Dragon download");

            std::fs::create_dir_all(&output)?;

            // Detect download tools
            let has_curl = std::process::Command::new("curl")
                .arg("--version")
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .status()
                .map(|s| s.success())
                .unwrap_or(false);

            let has_wget = std::process::Command::new("wget")
                .arg("--version")
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .status()
                .map(|s| s.success())
                .unwrap_or(false);

            if !has_curl && !has_wget {
                anyhow::bail!(
                    "Neither curl nor wget found in PATH. Please install one and retry."
                );
            }

            // Helper: download a single URL to a file
            let download_file = |url: &str, dest: &std::path::Path| -> Result<()> {
                log::info!("Downloading {} -> {:?}", url, dest);
                let status = if has_curl {
                    std::process::Command::new("curl")
                        .args(["-L", "--retry", "3", "--retry-delay", "5", "-o"])
                        .arg(dest)
                        .arg(url)
                        .arg("--progress-bar")
                        .status()?
                } else {
                    std::process::Command::new("wget")
                        .args(["--tries=3", "--waitretry=5", "-O"])
                        .arg(dest)
                        .arg(url)
                        .arg("--show-progress")
                        .status()?
                };
                if !status.success() {
                    anyhow::bail!("Download failed: {}", url);
                }
                Ok(())
            };

            match database.as_str() {
                // ===========================================================
                // Pre-built index downloads (tarball → extract)
                // ===========================================================
                "gtdb-r220" | "allthebacteria-v2" | "refseq-bacteria" => {
                    let url = match database.as_str() {
                        "gtdb-r220" => {
                            log::info!("Downloading GTDB r220 representative genomes index (~15 GB)");
                            "https://zenodo.org/records/dragon-indices/files/dragon-gtdb-r220.tar.gz"
                        }
                        "allthebacteria-v2" => {
                            log::info!("Downloading AllTheBacteria v2 index (~100 GB)");
                            "https://zenodo.org/records/dragon-indices/files/dragon-allthebacteria-v2.tar.gz"
                        }
                        "refseq-bacteria" => {
                            log::info!("Downloading RefSeq bacteria index (~50 GB)");
                            "https://zenodo.org/records/dragon-indices/files/dragon-refseq-bacteria.tar.gz"
                        }
                        _ => unreachable!(),
                    };

                    let tarball = output.join("dragon-index.tar.gz");
                    download_file(url, &tarball)?;

                    log::info!("Extracting {:?} -> {:?}", tarball, output);
                    let extract_status = std::process::Command::new("tar")
                        .arg("xzf")
                        .arg(&tarball)
                        .arg("-C")
                        .arg(&output)
                        .status()?;
                    if !extract_status.success() {
                        anyhow::bail!("Extraction failed");
                    }
                    std::fs::remove_file(&tarball)?;

                    log::info!(
                        "Download complete. Search with:\n  dragon search --index {:?} --query <reads.fasta>",
                        output
                    );
                }

                // ===========================================================
                // AllTheBacteria: download genomes → build index → package
                // ===========================================================
                "allthebacteria" => {
                    let atb_release = "0.2";
                    let atb_base = format!(
                        "https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/{}/assembly",
                        atb_release
                    );
                    log::info!(
                        "AllTheBacteria release {}: downloading genomes from EBI",
                        atb_release
                    );
                    log::info!("Source: {}", atb_base);

                    let genome_dir = output.join("genomes");
                    let index_dir = output.join("index");
                    let downloads_dir = output.join("downloads");
                    std::fs::create_dir_all(&genome_dir)?;
                    std::fs::create_dir_all(&index_dir)?;
                    std::fs::create_dir_all(&downloads_dir)?;

                    // Step 1: Fetch the file listing
                    let listing_file = downloads_dir.join("listing.html");
                    download_file(&atb_base, &listing_file)?;

                    let listing = std::fs::read_to_string(&listing_file)?;

                    // Parse tar archive URLs from directory listing
                    let archive_urls: Vec<String> = listing
                        .lines()
                        .filter_map(|line| {
                            // Match href="something.tar.{xz,gz,zst}"
                            let href_start = line.find("href=\"")?;
                            let rest = &line[href_start + 6..];
                            let href_end = rest.find('"')?;
                            let href = &rest[..href_end];
                            if href.ends_with(".tar.xz")
                                || href.ends_with(".tar.gz")
                                || href.ends_with(".tar.zst")
                            {
                                let url = if href.starts_with("http") {
                                    href.to_string()
                                } else {
                                    format!("{}/{}", atb_base, href)
                                };
                                Some(url)
                            } else {
                                None
                            }
                        })
                        .collect();

                    if archive_urls.is_empty() {
                        anyhow::bail!(
                            "No assembly archives found at {}. \
                             Check AllTheBacteria release and try:\n  \
                             bash aws/build_allthebacteria_index.sh",
                            atb_base
                        );
                    }

                    log::info!("Found {} assembly archive(s)", archive_urls.len());

                    // Step 2: Download + extract each archive
                    for (i, url) in archive_urls.iter().enumerate() {
                        let filename = url.rsplit('/').next().unwrap_or("archive.tar.xz");
                        let local_path = downloads_dir.join(filename);

                        log::info!(
                            "[{}/{}] Downloading {}",
                            i + 1,
                            archive_urls.len(),
                            filename
                        );
                        download_file(url, &local_path)?;

                        log::info!("Extracting {} -> {:?}", filename, genome_dir);
                        let tar_flag = if filename.ends_with(".tar.gz") {
                            "xzf"
                        } else if filename.ends_with(".tar.xz") {
                            "xJf"
                        } else {
                            // .tar.zst — extract via zstd pipe
                            let status = std::process::Command::new("sh")
                                .arg("-c")
                                .arg(format!(
                                    "zstd -d {:?} --stdout | tar -xf - -C {:?}",
                                    local_path, genome_dir
                                ))
                                .status()?;
                            if !status.success() {
                                log::warn!("Failed to extract {}", filename);
                            }
                            // Remove archive to save disk
                            let _ = std::fs::remove_file(&local_path);
                            continue;
                        };

                        let status = std::process::Command::new("tar")
                            .arg(tar_flag)
                            .arg(&local_path)
                            .arg("-C")
                            .arg(&genome_dir)
                            .status()?;
                        if !status.success() {
                            log::warn!("Failed to extract {}", filename);
                        }

                        // Remove archive to save disk
                        let _ = std::fs::remove_file(&local_path);
                    }

                    // Decompress any .gz FASTAs
                    log::info!("Decompressing gzipped FASTA files...");
                    let _ = std::process::Command::new("sh")
                        .arg("-c")
                        .arg(format!(
                            "find {:?} -name '*.fna.gz' -exec gunzip -f {{}} \\; 2>/dev/null; \
                             find {:?} -name '*.fa.gz' -exec gunzip -f {{}} \\; 2>/dev/null; \
                             find {:?} -name '*.fasta.gz' -exec gunzip -f {{}} \\; 2>/dev/null",
                            genome_dir, genome_dir, genome_dir
                        ))
                        .status();

                    // Flatten nested directories
                    let _ = std::process::Command::new("sh")
                        .arg("-c")
                        .arg(format!(
                            "find {:?} -mindepth 2 \\( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \\) \
                             -exec mv -n {{}} {:?}/ \\; 2>/dev/null; \
                             find {:?} -mindepth 1 -type d -empty -delete 2>/dev/null",
                            genome_dir, genome_dir, genome_dir
                        ))
                        .status();

                    let genome_count = dragon::io::fasta::list_fasta_files(&genome_dir)
                        .map(|v: Vec<std::path::PathBuf>| v.len())
                        .unwrap_or(0);
                    log::info!("AllTheBacteria: {} genome files ready", genome_count);

                    if genome_count == 0 {
                        anyhow::bail!("No FASTA files found after extraction");
                    }

                    // Step 3: Build index
                    log::info!("Building Dragon index (k={}, threads={})...", kmer_size, threads);
                    let max_ram_bytes = if low_memory {
                        log::info!("Low-memory mode: {:.1} GB RAM budget", max_ram);
                        Some((max_ram * 1_073_741_824.0) as usize)
                    } else {
                        None
                    };
                    dragon::index::build_index_with_options(
                        &genome_dir, &index_dir, kmer_size, threads, max_ram_bytes,
                    )?;

                    // Step 4: Package tarball
                    let tarball = output.join(format!(
                        "dragon_allthebacteria_v{}_k{}.tar.gz",
                        atb_release, kmer_size
                    ));
                    log::info!("Packaging index -> {:?}", tarball);
                    let status = std::process::Command::new("tar")
                        .arg("-czf")
                        .arg(&tarball)
                        .arg("-C")
                        .arg(&output)
                        .arg("index")
                        .status()?;
                    if !status.success() {
                        log::warn!("Tarball creation failed — index is still usable at {:?}", index_dir);
                    } else {
                        log::info!(
                            "Tarball created: {:?} ({:.2} GB)",
                            tarball,
                            std::fs::metadata(&tarball)?.len() as f64 / 1_073_741_824.0
                        );
                    }

                    log::info!(
                        "Done! {} genomes indexed.\n  \
                         Search: dragon search --index {:?} --query <reads.fasta>\n  \
                         Distribute: upload {:?} to Zenodo/S3",
                        genome_count, index_dir, tarball
                    );
                }

                // ===========================================================
                // RefSeq: download genomes → build index → package
                // ===========================================================
                db @ ("refseq" | "refseq-representative") => {
                    let representative_only = db == "refseq-representative";
                    let assembly_summary_url =
                        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt";

                    log::info!(
                        "RefSeq bacteria{}: downloading genomes from NCBI",
                        if representative_only { " (representative only)" } else { "" }
                    );

                    let genome_dir = output.join("genomes");
                    let index_dir = output.join("index");
                    let downloads_dir = output.join("downloads");
                    std::fs::create_dir_all(&genome_dir)?;
                    std::fs::create_dir_all(&index_dir)?;
                    std::fs::create_dir_all(&downloads_dir)?;

                    // Step 1: Download assembly summary
                    let summary_file = downloads_dir.join("assembly_summary.txt");
                    log::info!("Downloading assembly summary...");
                    download_file(assembly_summary_url, &summary_file)?;

                    // Step 2: Parse FTP URLs from assembly summary
                    let summary_text = std::fs::read_to_string(&summary_file)?;
                    let mut ftp_urls: Vec<String> = Vec::new();

                    for line in summary_text.lines() {
                        if line.starts_with('#') {
                            continue;
                        }
                        let cols: Vec<&str> = line.split('\t').collect();
                        if cols.len() < 20 {
                            continue;
                        }

                        let refseq_category = cols[4];
                        let assembly_level = cols[11];
                        let ftp_path = cols[19];

                        if ftp_path == "na" {
                            continue;
                        }

                        let include = if representative_only {
                            refseq_category == "representative genome"
                                || refseq_category == "reference genome"
                        } else {
                            assembly_level == "Complete Genome"
                                || assembly_level == "Chromosome"
                        };

                        if include {
                            // Build URL: ftp_path/basename_genomic.fna.gz
                            let basename = ftp_path.rsplit('/').next().unwrap_or("");
                            let url = format!(
                                "{}/{}_genomic.fna.gz",
                                ftp_path.replace("ftp://", "https://"),
                                basename
                            );
                            ftp_urls.push(url);
                        }
                    }

                    log::info!("Found {} genomes to download", ftp_urls.len());

                    if ftp_urls.is_empty() {
                        anyhow::bail!("No genomes found in assembly summary");
                    }

                    // Step 3: Download genomes in parallel using xargs
                    //
                    // Write URL list to file, then use xargs + curl for parallel downloads.
                    let url_list_file = downloads_dir.join("urls.txt");
                    std::fs::write(
                        &url_list_file,
                        ftp_urls
                            .iter()
                            .map(|u| u.as_str())
                            .collect::<Vec<_>>()
                            .join("\n"),
                    )?;

                    log::info!(
                        "Downloading {} genomes ({} parallel connections)...",
                        ftp_urls.len(),
                        parallel_downloads
                    );
                    log::info!("This may take several hours.");

                    // Download each genome: curl → gunzip → genome_dir/basename.fna
                    // We use a shell loop with background jobs for parallelism.
                    let download_script = format!(
                        r#"
download_one() {{
    URL="$1"
    BASENAME=$(basename "$URL" .gz)
    DEST="{genome_dir}/$BASENAME"
    [ -f "$DEST" ] && return 0
    curl -sL --retry 3 --retry-delay 2 "$URL" 2>/dev/null | gunzip -c > "$DEST" 2>/dev/null
    [ ! -s "$DEST" ] && rm -f "$DEST"
}}
export -f download_one 2>/dev/null

TOTAL=$(wc -l < "{url_file}")
DONE=0
while IFS= read -r URL; do
    while [ "$(jobs -r 2>/dev/null | wc -l)" -ge {parallel} ]; do
        sleep 0.2
    done
    download_one "$URL" &
    DONE=$((DONE + 1))
    if [ $((DONE % 5000)) -eq 0 ]; then
        echo "  Progress: $DONE / $TOTAL genomes..."
    fi
done < "{url_file}"
wait
echo "  Download complete."
"#,
                        genome_dir = genome_dir.display(),
                        url_file = url_list_file.display(),
                        parallel = parallel_downloads,
                    );

                    let status = std::process::Command::new("bash")
                        .arg("-c")
                        .arg(&download_script)
                        .status()?;

                    if !status.success() {
                        log::warn!("Some genome downloads may have failed — continuing with what we have");
                    }

                    let genome_count = dragon::io::fasta::list_fasta_files(&genome_dir)
                        .map(|v: Vec<std::path::PathBuf>| v.len())
                        .unwrap_or(0);
                    log::info!("RefSeq: {} genome files ready", genome_count);

                    if genome_count == 0 {
                        anyhow::bail!("No FASTA files found after download");
                    }

                    // Step 4: Build index
                    log::info!("Building Dragon index (k={}, threads={})...", kmer_size, threads);
                    let max_ram_bytes = if low_memory {
                        log::info!("Low-memory mode: {:.1} GB RAM budget", max_ram);
                        Some((max_ram * 1_073_741_824.0) as usize)
                    } else {
                        None
                    };
                    dragon::index::build_index_with_options(
                        &genome_dir, &index_dir, kmer_size, threads, max_ram_bytes,
                    )?;

                    // Step 5: Package tarball
                    let mode = if representative_only {
                        "representative"
                    } else {
                        "complete"
                    };
                    let tarball = output.join(format!(
                        "dragon_refseq_{}_k{}.tar.gz",
                        mode, kmer_size
                    ));
                    log::info!("Packaging index -> {:?}", tarball);
                    let status = std::process::Command::new("tar")
                        .arg("-czf")
                        .arg(&tarball)
                        .arg("-C")
                        .arg(&output)
                        .arg("index")
                        .status()?;
                    if !status.success() {
                        log::warn!("Tarball creation failed — index is still usable at {:?}", index_dir);
                    } else {
                        log::info!(
                            "Tarball created: {:?} ({:.2} GB)",
                            tarball,
                            std::fs::metadata(&tarball)?.len() as f64 / 1_073_741_824.0
                        );
                    }

                    log::info!(
                        "Done! {} genomes indexed.\n  \
                         Search: dragon search --index {:?} --query <reads.fasta>\n  \
                         Distribute: upload {:?} to Zenodo/S3",
                        genome_count, index_dir, tarball
                    );
                }

                // ===========================================================
                // Custom URL: download pre-built index tarball
                // ===========================================================
                custom if custom.starts_with("http") => {
                    log::info!("Downloading custom index from: {}", custom);

                    let tarball = output.join("dragon-index.tar.gz");
                    download_file(custom, &tarball)?;

                    log::info!("Extracting {:?} -> {:?}", tarball, output);
                    let extract_status = std::process::Command::new("tar")
                        .arg("xzf")
                        .arg(&tarball)
                        .arg("-C")
                        .arg(&output)
                        .status()?;
                    if !extract_status.success() {
                        anyhow::bail!("Extraction failed");
                    }
                    std::fs::remove_file(&tarball)?;

                    log::info!(
                        "Download complete. Search with:\n  dragon search --index {:?} --query <reads.fasta>",
                        output
                    );
                }

                _ => {
                    anyhow::bail!(
                        "Unknown database: {}.\n\n\
                         Pre-built indices (download only):\n  \
                           gtdb-r220            GTDB r220 representative genomes (~15 GB)\n  \
                           allthebacteria-v2    AllTheBacteria v2 index (~100 GB)\n  \
                           refseq-bacteria      RefSeq bacteria index (~50 GB)\n\n\
                         Download genomes + build index:\n  \
                           allthebacteria       Download from EBI + build index\n  \
                           refseq               Download complete genomes from NCBI + build\n  \
                           refseq-representative Download representative genomes only + build\n\n\
                         Custom:\n  \
                           http://...           Any URL to a .tar.gz index tarball",
                        database
                    );
                }
            }
        }

        Commands::SignalIndex {
            input,
            output,
            num_levels,
            threads,
            signal_boundaries,
            pore_model,
        } => {
            log::info!("Dragon signal index construction");
            log::info!("Input: {:?}", input);
            log::info!("Output: {:?}", output);
            log::info!("Discrete levels: {}", num_levels);

            // Load alphabet: use learned boundaries if provided, else equal-width
            let alphabet = if let Some(ref bounds_path) = signal_boundaries {
                let data = std::fs::read_to_string(bounds_path)?;
                let boundaries: Vec<f32> = serde_json::from_str(&data)?;
                log::info!("Loaded {} learned signal boundaries from {}", boundaries.len(), bounds_path);
                dragon::signal::SignalAlphabet::from_boundaries(boundaries)
            } else {
                dragon::signal::SignalAlphabet::new(num_levels)
            };

            // Load pore model: use custom if provided, else default
            let model = if let Some(ref model_path) = pore_model {
                let path = std::path::Path::new(model_path);
                dragon::signal::model::PoreModel::load_from_file(path)
                    .unwrap_or_else(|| {
                        log::warn!("Failed to load custom pore model, using default");
                        dragon::signal::model::load_default_model()
                    })
            } else {
                dragon::signal::model::load_default_model()
            };

            let config = dragon::signal::SignalIndexConfig {
                alphabet,
                pore_model: model,
                threads,
            };

            dragon::signal::index::build_signal_index(&input, &output, &config)?;

            log::info!("Signal index construction complete");
        }

        Commands::SignalSearch {
            index,
            query,
            output,
            signal_kmer_size,
            min_hits,
            max_seed_freq,
            max_results,
            threads,
        } => {
            log::info!("Dragon signal search");
            log::info!("Index: {:?}", index);
            log::info!("Query: {:?}", query);

            let config = dragon::signal::SignalSearchConfig {
                index_dir: index.into(),
                signal_kmer_size,
                min_hits,
                max_seed_freq,
                max_results,
                threads,
            };

            let results = dragon::signal::search::search_signal_file(&query, &config)?;

            // Write output
            let mut writer: Box<dyn std::io::Write> = if output == "-" {
                Box::new(std::io::stdout())
            } else {
                Box::new(std::fs::File::create(&output)?)
            };

            dragon::signal::search::write_signal_results(&results, &mut writer)?;

            log::info!(
                "Signal search complete: {} reads, {} total hits",
                results.len(),
                results.iter().map(|r| r.hits.len()).sum::<usize>()
            );
        }
    }

    Ok(())
}
