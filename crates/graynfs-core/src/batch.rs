use crate::options::{FactorOptions, estimate_collect_bytes};

pub fn normalize_batch_jobs(jobs: usize) -> usize {
    jobs.max(1)
}

pub fn batch_job_threads(total_workers: usize, jobs: usize) -> usize {
    (total_workers.max(1) / normalize_batch_jobs(jobs)).max(1)
}

pub fn batch_job_options(options: &FactorOptions, jobs: usize) -> FactorOptions {
    let total = options.effective_workers();
    let mut scoped = options.clone();
    scoped.threads = Some(batch_job_threads(total, jobs));
    scoped
}

pub fn batch_per_job_bytes(options: &FactorOptions) -> u64 {
    let columns = options
        .rational_base_bound
        .saturating_add(options.algebraic_base_bound)
        .max(1);
    estimate_collect_bytes(
        columns,
        options.sieve_a_bound as i64,
        options.effective_workers().max(1),
    )
}

pub fn batch_concurrency(budget_bytes: u64, per_job_bytes: u64, jobs: usize) -> usize {
    let jobs = normalize_batch_jobs(jobs);
    if per_job_bytes == 0 {
        return jobs;
    }
    if per_job_bytes > budget_bytes {
        return 1;
    }
    let slots = (budget_bytes / per_job_bytes.max(1)).min(usize::MAX as u64) as usize;
    slots.min(jobs).max(1)
}

pub fn batch_slots(options: &FactorOptions, jobs: usize) -> usize {
    let scoped = batch_job_options(options, jobs);
    batch_concurrency(
        options.memory_budget_bytes(),
        batch_per_job_bytes(&scoped),
        jobs,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_batch_jobs_clamps_zero_to_one() {
        assert_eq!(normalize_batch_jobs(0), 1);
        assert_eq!(normalize_batch_jobs(1), 1);
        assert_eq!(normalize_batch_jobs(4), 4);
    }

    #[test]
    fn job_threads_divide_total_with_floor_one() {
        assert_eq!(batch_job_threads(16, 1), 16);
        assert_eq!(batch_job_threads(16, 4), 4);
        assert_eq!(batch_job_threads(7, 2), 3);
        assert_eq!(batch_job_threads(3, 4), 1);
        assert_eq!(batch_job_threads(16, 0), 16);
    }

    #[test]
    fn job_options_preserve_everything_except_threads() {
        let options = FactorOptions {
            threads: Some(16),
            seed: 111,
            max_memory_mb: Some(512),
            ..FactorOptions::default()
        };
        let scoped = batch_job_options(&options, 4);
        assert_eq!(scoped.threads, Some(4));
        assert_eq!(scoped.seed, 111);
        assert_eq!(scoped.max_memory_mb, Some(512));
        assert_eq!(scoped.effective_workers(), 4);
    }

    #[test]
    fn job_options_cover_unset_thread_counts() {
        let options = FactorOptions::default();
        let scoped = batch_job_options(&options, 1);
        assert_eq!(scoped.effective_workers(), options.effective_workers());
    }

    #[test]
    fn concurrency_caps_at_requested_jobs() {
        assert_eq!(batch_concurrency(u64::MAX, 1, 4), 4);
        assert_eq!(batch_concurrency(100, 30, 4), 3);
        assert_eq!(batch_concurrency(100, 30, 2), 2);
    }

    #[test]
    fn concurrency_queues_over_budget_instead_of_aborting() {
        assert_eq!(batch_concurrency(0, 351_268, 8), 1);
        assert_eq!(batch_concurrency(100, 101, 8), 1);
        assert_eq!(batch_concurrency(100, 0, 8), 8);
    }

    #[test]
    fn slots_pin_default_options_and_over_budget_gate() {
        let options = FactorOptions {
            threads: Some(16),
            ..FactorOptions::default()
        };
        assert_eq!(batch_slots(&options, 4), 4);
        let starved = FactorOptions {
            threads: Some(4),
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        assert_eq!(batch_slots(&starved, 8), 1);
    }
}
