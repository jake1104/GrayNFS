use std::time::Duration;

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
    pub fn validate(&self) -> Result<(), crate::FactorError> {
        if self.degree == Some(0) {
            return Err(crate::FactorError::InvalidInput(
                "degree must be greater than zero".to_string(),
            ));
        }
        if !(0.0..=0.5).contains(&self.degree_race_fraction) {
            return Err(crate::FactorError::InvalidInput(
                "degree race fraction must be between zero and one half".to_string(),
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
