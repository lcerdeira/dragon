/// Hardware profiles for resource-constrained deployment.
///
/// Dragon supports two hardware profiles that control resource usage across
/// all operations (index build, search, summarize):
///
/// - **Laptop** (≤8 GB RAM, ≤4 threads, single SSD):
///   Enables external-memory suffix array construction, limits concurrent
///   operations, caps candidate genomes per query. Designed for field
///   deployments and resource-constrained environments.
///
/// - **Workstation** (64-128 GB RAM, full threading, NVMe):
///   Full resources, in-memory construction, no artificial limits.
///   Designed for lab servers and cloud instances.

use std::fmt;

/// A hardware profile controlling resource limits.
#[derive(Clone, Debug)]
pub struct HardwareProfile {
    pub name: String,
    /// Maximum RAM in bytes for index construction.
    pub max_ram_bytes: usize,
    /// Maximum threads for parallel operations.
    pub max_threads: usize,
    /// Whether to use external-memory (disk-based) suffix array construction.
    pub external_memory: bool,
    /// Maximum candidate genomes to chain per query (limits computation).
    pub max_candidates: usize,
    /// Maximum target sequences to report per query.
    pub max_target_seqs: usize,
    /// Maximum seed frequency (skip seeds more frequent than this).
    pub max_seed_freq: usize,
}

impl HardwareProfile {
    /// Laptop profile: ≤8 GB RAM, ≤4 threads, conservative limits.
    pub fn laptop() -> Self {
        Self {
            name: "laptop".to_string(),
            max_ram_bytes: 8 * 1_073_741_824, // 8 GB
            max_threads: 4,
            external_memory: true,
            max_candidates: 50,
            max_target_seqs: 10,
            max_seed_freq: 5_000,
        }
    }

    /// Workstation profile: full resources, no artificial limits.
    pub fn workstation() -> Self {
        let cpus = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8);
        Self {
            name: "workstation".to_string(),
            max_ram_bytes: 64 * 1_073_741_824, // 64 GB default
            max_threads: cpus,
            external_memory: false,
            max_candidates: 500,
            max_target_seqs: 50,
            max_seed_freq: 10_000,
        }
    }

    /// Parse a profile from a string name, with optional RAM/thread overrides.
    pub fn from_name(name: &str) -> Self {
        match name {
            "laptop" => Self::laptop(),
            "workstation" | _ => Self::workstation(),
        }
    }

    /// Apply user overrides (CLI flags take precedence over profile defaults).
    pub fn with_overrides(
        mut self,
        max_ram_gb: Option<f64>,
        threads: Option<usize>,
    ) -> Self {
        if let Some(ram) = max_ram_gb {
            self.max_ram_bytes = (ram * 1_073_741_824.0) as usize;
        }
        if let Some(t) = threads {
            self.max_threads = t;
        }
        self
    }

    /// Get RAM budget for external-memory construction, or None for in-memory.
    pub fn ram_budget(&self) -> Option<usize> {
        if self.external_memory {
            Some(self.max_ram_bytes)
        } else {
            None
        }
    }

    /// Log the active profile settings.
    pub fn log_settings(&self) {
        log::info!("Hardware profile: {}", self.name);
        log::info!(
            "  RAM: {:.1} GB, Threads: {}, External-memory: {}",
            self.max_ram_bytes as f64 / 1_073_741_824.0,
            self.max_threads,
            self.external_memory,
        );
        log::info!(
            "  Max candidates: {}, Max targets: {}, Max seed freq: {}",
            self.max_candidates,
            self.max_target_seqs,
            self.max_seed_freq,
        );
    }
}

impl fmt::Display for HardwareProfile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} ({:.0} GB RAM, {} threads{})",
            self.name,
            self.max_ram_bytes as f64 / 1_073_741_824.0,
            self.max_threads,
            if self.external_memory { ", disk-based SA" } else { "" },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_laptop_profile() {
        let p = HardwareProfile::laptop();
        assert_eq!(p.max_ram_bytes, 8 * 1_073_741_824);
        assert_eq!(p.max_threads, 4);
        assert!(p.external_memory);
    }

    #[test]
    fn test_workstation_profile() {
        let p = HardwareProfile::workstation();
        assert!(!p.external_memory);
        assert!(p.max_threads >= 1);
    }

    #[test]
    fn test_overrides() {
        let p = HardwareProfile::laptop().with_overrides(Some(4.0), Some(2));
        assert_eq!(p.max_ram_bytes, 4 * 1_073_741_824);
        assert_eq!(p.max_threads, 2);
    }

    #[test]
    fn test_from_name() {
        let laptop = HardwareProfile::from_name("laptop");
        assert_eq!(laptop.name, "laptop");
        let ws = HardwareProfile::from_name("workstation");
        assert_eq!(ws.name, "workstation");
    }
}
