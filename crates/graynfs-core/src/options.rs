use std::time::Duration;

pub const DEFAULT_MEMORY_BUDGET_BYTES: u64 = 1 << 30;
pub const MEMORY_BYTES_PER_COLUMN: u64 = 64;
pub const MEMORY_BYTES_PER_SIEVE_POSITION: u64 = 68;

pub fn estimate_collect_bytes(columns: usize, half_width: i64, workers: usize) -> u64 {
    let width = half_width.max(0) as u64 * 2 + 1;
    columns as u64 * MEMORY_BYTES_PER_COLUMN
        + width * MEMORY_BYTES_PER_SIEVE_POSITION * workers.max(1) as u64
}

pub fn estimate_rank_bytes(columns: usize, relations: usize) -> u64 {
    let mask_words = columns.div_ceil(64) as u64;
    let relation_words = relations.div_ceil(64) as u64;
    columns as u64 * (mask_words + relation_words) * 8
}

pub fn estimate_total_bytes(columns: usize, relations: usize) -> u64 {
    let mask_words = columns.div_ceil(64) as u64;
    let per_relation = mask_words.saturating_mul(8).saturating_add(256);
    estimate_rank_bytes(columns, relations)
        .saturating_add((relations as u64).saturating_mul(per_relation))
        .saturating_add((columns as u64).saturating_mul(8))
}

pub fn estimate_relation_matrix_megabytes(column_count: usize) -> usize {
    let columns = column_count.max(1);
    let rows = (columns.saturating_mul(3) / 2 + 16).max(16);
    let bytes = estimate_total_bytes(columns, rows);
    usize::try_from(bytes.div_ceil(1024 * 1024)).unwrap_or(usize::MAX)
}

#[derive(Clone, Debug)]
pub struct FactorOptions {
    pub degree: Option<usize>,
    pub rational_base_bound: usize,
    pub algebraic_base_bound: usize,
    pub sieve_a_bound: usize,
    pub sieve_b_bound: usize,
    pub large_prime_bound: usize,
    pub polynomial_candidates: usize,
    pub seed: u64,
    pub threads: Option<usize>,
    pub max_duration: Option<Duration>,
    pub max_memory_mb: Option<usize>,
    pub allow_fallback: bool,
    pub use_gnfs: bool,
    pub skip_preprocessing: bool,
    pub enable_double_large_primes: bool,
    pub degree_race_fraction: f64,
    pub use_fermat: bool,
    pub use_squfof: bool,
    pub use_pm1: bool,
    pub use_pp1: bool,
    pub use_ecm: bool,
    pub use_mpqs: bool,
    pub perf_mode: bool,
}

impl Default for FactorOptions {
    fn default() -> Self {
        Self {
            degree: None,
            rational_base_bound: 2_000,
            algebraic_base_bound: 2_000,
            sieve_a_bound: 700,
            sieve_b_bound: 120,
            large_prime_bound: 50_000,
            polynomial_candidates: 8,
            seed: 0x9e3779b97f4a7c15,
            threads: None,
            max_duration: None,
            max_memory_mb: None,
            allow_fallback: true,
            use_gnfs: true,
            skip_preprocessing: false,
            enable_double_large_primes: false,
            degree_race_fraction: 0.2,
            use_fermat: true,
            use_squfof: true,
            use_pm1: true,
            use_pp1: true,
            use_ecm: true,
            use_mpqs: true,
            perf_mode: false,
        }
    }
}

impl FactorOptions {
    pub fn memory_budget_bytes(&self) -> u64 {
        match self.max_memory_mb {
            Some(mb) => (mb as u64).saturating_mul(1024 * 1024),
            None => DEFAULT_MEMORY_BUDGET_BYTES,
        }
    }

    pub fn validate(&self) -> Result<(), crate::FactorError> {
        if let Some(degree) = self.degree
            && !(2..=5).contains(&degree)
        {
            return Err(crate::FactorError::InvalidInput(format!(
                "unsupported degree {degree}: supported degrees are 2, 3, 4 and 5"
            )));
        }
        if !self.use_squfof {
            return Err(crate::FactorError::InvalidInput(
                "use_squfof is reserved for a future SQUFOF stage and has no effect".to_string(),
            ));
        }
        if !(0.05..=0.5).contains(&self.degree_race_fraction) {
            return Err(crate::FactorError::InvalidInput(
                "degree race fraction must be between 0.05 and 0.5".to_string(),
            ));
        }
        if self.rational_base_bound < 2 || self.algebraic_base_bound < 2 {
            return Err(crate::FactorError::InvalidInput(
                "factor base bounds must be at least two".to_string(),
            ));
        }
        if self.sieve_a_bound == 0 || self.sieve_b_bound == 0 {
            return Err(crate::FactorError::InvalidInput(
                "sieve bounds must be greater than zero".to_string(),
            ));
        }
        Ok(())
    }

    pub fn with_base_bounds(mut self, rational: usize, algebraic: usize) -> Self {
        self.rational_base_bound = rational;
        self.algebraic_base_bound = algebraic;
        self
    }

    pub fn with_sieve_bounds(mut self, a: usize, b: usize) -> Self {
        self.sieve_a_bound = a;
        self.sieve_b_bound = b;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_budget_defaults_to_one_gibibyte() {
        let options = FactorOptions::default();
        assert_eq!(options.max_memory_mb, None);
        assert_eq!(options.memory_budget_bytes(), DEFAULT_MEMORY_BUDGET_BYTES);
        assert_eq!(DEFAULT_MEMORY_BUDGET_BYTES, 1 << 30);
    }

    #[test]
    fn memory_budget_scales_megabytes() {
        let options = FactorOptions {
            max_memory_mb: Some(512),
            ..FactorOptions::default()
        };
        assert_eq!(options.memory_budget_bytes(), 512 * 1024 * 1024);
        let zero = FactorOptions {
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        assert_eq!(zero.memory_budget_bytes(), 0);
    }

    #[test]
    fn estimators_grow_with_problem_size() {
        let small = estimate_collect_bytes(100, 700, 1);
        let large = estimate_collect_bytes(1_000, 7_000, 1);
        assert!(large > small);
        let small_rank = estimate_rank_bytes(100, 100);
        let large_rank = estimate_rank_bytes(1_000, 1_000);
        assert!(large_rank > small_rank);
        assert!(estimate_relation_matrix_megabytes(200_000) > 1_024);
    }

    #[test]
    fn relation_matrix_estimate_covers_rank_slabs() {
        for columns in [600_usize, 2000, 200_000] {
            let rows = columns * 3 / 2;
            let rank_bytes = estimate_rank_bytes(columns, rows);
            let matrix_bytes = estimate_relation_matrix_megabytes(columns) as u64 * 1024 * 1024;
            assert!(
                matrix_bytes >= rank_bytes,
                "columns={columns} rows={rows} matrix={matrix_bytes} rank={rank_bytes}"
            );
        }
    }

    #[test]
    fn total_bytes_grows_monotonically() {
        assert!(estimate_total_bytes(600, 900) < estimate_total_bytes(2000, 3000));
        assert!(estimate_total_bytes(2000, 100) < estimate_total_bytes(2000, 3000));
        assert!(estimate_total_bytes(2000, 3000) >= estimate_rank_bytes(2000, 3000));
    }

    #[test]
    fn large_columns_relation_storage_scales_with_mask_width() {
        let columns = 200_000_usize;
        let relations = 300_000_usize;
        let mask_words = columns.div_ceil(64) as u64;
        let floor = estimate_rank_bytes(columns, relations)
            .saturating_add((relations as u64).saturating_mul(mask_words.saturating_mul(8)));
        assert!(estimate_total_bytes(columns, relations) >= floor);
    }

    #[test]
    fn degree_accepts_supported_two_through_five() {
        for degree in [2_usize, 3, 4, 5] {
            let options = FactorOptions {
                degree: Some(degree),
                ..FactorOptions::default()
            };
            assert!(options.validate().is_ok());
        }
        assert!(FactorOptions::default().validate().is_ok());
    }

    #[test]
    fn degree_outside_supported_range_is_invalid() {
        for degree in [0_usize, 1, 6, 7, 100] {
            let options = FactorOptions {
                degree: Some(degree),
                ..FactorOptions::default()
            };
            assert!(matches!(
                options.validate(),
                Err(crate::FactorError::InvalidInput(_))
            ));
        }
    }

    #[test]
    fn reserved_squfof_toggle_is_rejected() {
        let options = FactorOptions {
            use_squfof: false,
            ..FactorOptions::default()
        };
        assert!(matches!(
            options.validate(),
            Err(crate::FactorError::InvalidInput(_))
        ));
        assert!(FactorOptions::default().validate().is_ok());
    }

    #[test]
    fn degree_race_fraction_boundaries() {
        for fraction in [0.0, 0.04, 0.51] {
            let options = FactorOptions {
                degree_race_fraction: fraction,
                ..FactorOptions::default()
            };
            assert!(
                matches!(options.validate(), Err(crate::FactorError::InvalidInput(_))),
                "fraction {fraction} should be rejected"
            );
        }
        for fraction in [0.05, 0.2, 0.5] {
            let options = FactorOptions {
                degree_race_fraction: fraction,
                ..FactorOptions::default()
            };
            assert!(
                options.validate().is_ok(),
                "fraction {fraction} should be accepted"
            );
        }
    }
}
