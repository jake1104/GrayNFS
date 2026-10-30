use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, PoisonError};
use std::time::{Duration, Instant};

use crate::integer::{Integer, IntegerExt};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::algebraic::QuadraticRing;
use crate::base::FactorBase;
use crate::error::FactorError;
use crate::general;
use crate::linear_algebra::rank_and_dependencies;
use crate::options::FactorOptions;
use crate::polynomial::{PolynomialPair, select_polynomial};
use crate::relations::collect_relations;
use crate::square_root::extract_factor;

#[derive(Clone, Debug)]
pub struct FactorReport {
    pub input: String,
    pub factors: Vec<String>,
    pub method: String,
    pub degree: Option<usize>,
    pub polynomial: Option<String>,
    pub root: Option<String>,
    pub relations: usize,
    pub rank: usize,
    pub elapsed_ms: u128,
    pub timings: Vec<StageTiming>,
}

#[derive(Clone, Debug)]
pub struct ProgressEvent {
    pub stage: String,
    pub current: u64,
    pub total: u64,
    pub message: String,
}

#[derive(Clone, Debug)]
pub struct StageTiming {
    pub stage: String,
    pub degree: Option<usize>,
    pub elapsed_ms: u128,
    pub detail: String,
}

#[derive(Default)]
struct StageOutcome {
    winners: Vec<String>,
}

impl StageOutcome {
    fn record(&mut self, stage: &str) {
        if !self.winners.iter().any(|existing| existing == stage) {
            self.winners.push(stage.to_string());
        }
    }

    fn merge(&mut self, other: StageOutcome) {
        for stage in other.winners {
            self.record(&stage);
        }
    }

    fn method(&self) -> String {
        if self.winners.is_empty() {
            return "preprocessing".to_string();
        }
        self.winners.join("+")
    }
}

#[derive(Default)]
pub struct TimingCollector {
    entries: Mutex<Vec<StageTiming>>,
    cpu_ms: Mutex<u128>,
}

impl TimingCollector {
    pub(crate) fn record(
        &self,
        stage: &str,
        degree: Option<usize>,
        started: Instant,
        detail: String,
    ) -> u128 {
        let elapsed_ms = started.elapsed().as_millis();
        if let Ok(mut entries) = self.entries.lock() {
            entries.push(StageTiming {
                stage: stage.to_string(),
                degree,
                elapsed_ms,
                detail,
            });
        }
        elapsed_ms
    }

    pub(crate) fn record_elapsed(
        &self,
        stage: &str,
        degree: Option<usize>,
        elapsed_ms: u128,
        detail: String,
    ) {
        if let Ok(mut entries) = self.entries.lock() {
            entries.push(StageTiming {
                stage: stage.to_string(),
                degree,
                elapsed_ms,
                detail,
            });
        }
    }

    pub(crate) fn snapshot(&self) -> Vec<StageTiming> {
        self.entries
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .clone()
    }

    pub(crate) fn add_cpu_ms(&self, elapsed_ms: u128, workers: usize) {
        let charge = elapsed_ms.saturating_mul(workers.max(1) as u128);
        if let Ok(mut used) = self.cpu_ms.lock() {
            *used = used.saturating_add(charge);
        }
    }

    pub(crate) fn cpu_used_ms(&self) -> u128 {
        self.cpu_ms.lock().map(|used| *used).unwrap_or(0)
    }

    pub(crate) fn check_cpu_budget(&self, options: &FactorOptions) -> Result<(), FactorError> {
        if let Some(cap_ms) = options.cpu_budget_ms()
            && self.cpu_used_ms() >= cap_ms
        {
            return Err(FactorError::ResourceLimit(format!(
                "cpu time budget exceeded: used {} ms across {} threads vs budget {} ms (cap {} s)",
                self.cpu_used_ms(),
                options.effective_workers(),
                cap_ms,
                options.max_cpu_seconds.unwrap_or(0),
            )));
        }
        Ok(())
    }

    pub(crate) fn merge(&self, entries: Vec<StageTiming>) {
        if entries.is_empty() {
            return;
        }
        if let Ok(mut current) = self.entries.lock() {
            current.extend(entries);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn record_memory_shadow(
    collector: Option<&TimingCollector>,
    site: &str,
    degree: Option<usize>,
    columns: usize,
    relations: usize,
    partial_singles: usize,
    partial_doubles: usize,
    max_a: i64,
    workers: usize,
    budget: u64,
    verdict_allow: bool,
) {
    let Some(collector) = collector else {
        return;
    };
    let collect = crate::options::estimate_collect_bytes(columns, max_a, workers);
    let rank = crate::options::estimate_rank_bytes(columns, relations);
    let total = crate::options::estimate_total_bytes(columns, relations);
    let live =
        crate::options::estimate_live_bytes(columns, relations, partial_singles, partial_doubles);
    let verdict = if verdict_allow { "allow" } else { "reject" };
    collector.record_elapsed(
        "memory_shadow",
        degree,
        0,
        format!(
            "site={site} columns={columns} relations={relations} singles={partial_singles} doubles={partial_doubles} max_a={max_a} workers={workers} collect={collect} rank={rank} total={total} live={live} budget={budget} verdict={verdict}"
        ),
    );
}

pub(crate) struct Timer<'a> {
    collector: Option<&'a TimingCollector>,
    stage: &'static str,
    degree: Option<usize>,
    started: Instant,
}

impl<'a> Timer<'a> {
    pub(crate) fn new(
        collector: Option<&'a TimingCollector>,
        stage: &'static str,
        degree: Option<usize>,
    ) -> Self {
        Self {
            collector,
            stage,
            degree,
            started: Instant::now(),
        }
    }

    pub(crate) fn finish(self, detail: String) -> u128 {
        if let Some(collector) = self.collector {
            collector.record(self.stage, self.degree, self.started, detail)
        } else {
            self.started.elapsed().as_millis()
        }
    }
}

struct GnfsStats {
    polynomial: PolynomialPair,
    relations: usize,
    rank: usize,
}

pub fn factor(n: &Integer, options: &FactorOptions) -> Result<FactorReport, FactorError> {
    factor_with_progress(n, options, |_| Ok(()))
}

pub fn factor_with_progress<F>(
    n: &Integer,
    options: &FactorOptions,
    progress: F,
) -> Result<FactorReport, FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let collector = TimingCollector::default();
    factor_collected(n, options, progress, &collector)
}

pub fn factor_with_timings<F>(
    n: &Integer,
    options: &FactorOptions,
    progress: F,
) -> (Result<FactorReport, FactorError>, Vec<StageTiming>)
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let collector = TimingCollector::default();
    let result = factor_collected(n, options, progress, &collector);
    let timings = collector.snapshot();
    (result, timings)
}

fn factor_collected<F>(
    n: &Integer,
    options: &FactorOptions,
    mut progress: F,
    collector: &TimingCollector,
) -> Result<FactorReport, FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let started = Instant::now();
    options.validate()?;
    check_cpu(Some(collector), options)?;
    if n < &Integer::from(1) {
        return Err(FactorError::InvalidInput(
            "N must be a positive integer".to_string(),
        ));
    }
    if n == &Integer::from(1) {
        return Ok(FactorReport {
            input: n.to_string(),
            factors: Vec::new(),
            method: "trivial".to_string(),
            degree: None,
            polynomial: None,
            root: None,
            relations: 0,
            rank: 0,
            elapsed_ms: started.elapsed().as_millis(),
            timings: Vec::new(),
        });
    }
    let cancelled = AtomicBool::new(false);
    let deadline = options.effective_deadline(Instant::now());
    emit(
        &mut progress,
        "preprocessing",
        0,
        7,
        "checking primality and small factors",
    )?;
    let mut factors = Vec::new();
    let mut outcome = StageOutcome::default();
    let mut polynomial = None;
    let mut root = None;
    let mut relation_count = 0;
    let mut rank = 0;
    let mut rng = ChaCha8Rng::seed_from_u64(options.seed);
    factor_recursive(
        n,
        options,
        &cancelled,
        deadline,
        &mut rng,
        &mut factors,
        &mut outcome,
        &mut polynomial,
        &mut root,
        &mut relation_count,
        &mut rank,
        &mut progress,
        Some(collector),
    )?;
    factors.sort();
    let mut product = Integer::from(1);
    for value in &factors {
        product *= value;
    }
    if product != *n {
        return Err(FactorError::Factorization(
            "reported factors do not multiply back to the input".to_string(),
        ));
    }
    let method = outcome.method();
    let degree = polynomial.as_ref().map(|value| value.degree);
    Ok(FactorReport {
        input: n.to_string(),
        factors: factors.iter().map(ToString::to_string).collect(),
        method,
        degree,
        polynomial: polynomial.map(|value| value.describe()),
        root,
        relations: relation_count,
        rank,
        elapsed_ms: started.elapsed().as_millis(),
        timings: collector.snapshot(),
    })
}

#[allow(clippy::too_many_arguments)]
fn factor_recursive<F>(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    rng: &mut ChaCha8Rng,
    factors: &mut Vec<Integer>,
    outcome: &mut StageOutcome,
    polynomial: &mut Option<PolynomialPair>,
    root: &mut Option<String>,
    relation_count: &mut usize,
    rank: &mut usize,
    progress: &mut F,
    timings: Option<&TimingCollector>,
) -> Result<(), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    if n == &Integer::from(1) {
        return Ok(());
    }
    check_cpu(timings, options)?;
    if cancelled.load(Ordering::Relaxed) {
        return Err(FactorError::Cancelled);
    }
    if crate::prime::is_probable_prime_with_options(n, options) {
        factors.push(n.clone());
        if outcome.winners.is_empty() {
            outcome.record("prime");
        }
        return Ok(());
    }
    if let Some(square_root) = perfect_square(n) {
        outcome.record("perfect_square");
        let mut square_factors = Vec::new();
        let mut square_outcome = StageOutcome::default();
        let mut square_polynomial = None;
        let mut square_root_value = None;
        let mut square_relations = 0;
        let mut square_rank = 0;
        factor_recursive(
            &square_root,
            options,
            cancelled,
            deadline,
            rng,
            &mut square_factors,
            &mut square_outcome,
            &mut square_polynomial,
            &mut square_root_value,
            &mut square_relations,
            &mut square_rank,
            progress,
            timings,
        )?;
        for value in square_factors {
            factors.push(value.clone());
            factors.push(value);
        }
        outcome.merge(square_outcome);
        if polynomial.is_none() {
            *polynomial = square_polynomial;
        }
        if root.is_none() {
            *root = square_root_value;
        }
        *relation_count += square_relations;
        *rank += square_rank;
        return Ok(());
    }
    if let Some((base, exponent)) = perfect_power_root(n) {
        outcome.record("perfect_power");
        let mut power_factors = Vec::new();
        let mut power_outcome = StageOutcome::default();
        let mut power_polynomial = None;
        let mut power_root_value = None;
        let mut power_relations = 0;
        let mut power_rank = 0;
        factor_recursive(
            &base,
            options,
            cancelled,
            deadline,
            rng,
            &mut power_factors,
            &mut power_outcome,
            &mut power_polynomial,
            &mut power_root_value,
            &mut power_relations,
            &mut power_rank,
            progress,
            timings,
        )?;
        for value in power_factors {
            for _ in 0..exponent {
                factors.push(value.clone());
            }
        }
        outcome.merge(power_outcome);
        if polynomial.is_none() {
            *polynomial = power_polynomial;
        }
        if root.is_none() {
            *root = power_root_value;
        }
        *relation_count += power_relations;
        *rank += power_rank;
        return Ok(());
    }
    let snfs_note = match detect_special_form(n) {
        Some(form) => form.describe(),
        None => "none".to_string(),
    };
    let trial_checked = if !options.skip_preprocessing {
        let (hit, checked) = trial_division(n);
        if let Some(divisor) = hit {
            outcome.record("trial_division");
            let cofactor = n / &divisor;
            return factor_recursive(
                &divisor,
                options,
                cancelled,
                deadline,
                rng,
                factors,
                outcome,
                polynomial,
                root,
                relation_count,
                rank,
                progress,
                timings,
            )
            .and_then(|_| {
                factor_recursive(
                    &cofactor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
            });
        }
        checked
    } else {
        0
    };
    let mut signals = crate::cascade::CascadeSignals::collect_with_options(n, options, deadline);
    if !options.skip_preprocessing {
        signals.observe_trial_division(trial_checked);
    }
    let (fermat_budget, pm1_budget, pp1_budget, ecm_budget) =
        crate::cascade::allocate_slices(&signals, signals.remaining_ms(deadline));
    let fermat_slice = budget_to_deadline(fermat_budget);
    let pm1_slice = budget_to_deadline(pm1_budget);
    let pp1_slice = budget_to_deadline(pp1_budget);
    let mut cheap_budget: Option<u128> = None;
    let workers = options.effective_workers();
    signals.observe_workers(workers);
    let trial_full_scan = trial_checked >= crate::prime::small_primes().len() as u64;
    let stage_one = crate::pm1::pm1_bound_for_slice(
        signals.target_digits,
        signals.remaining_ms(deadline),
        options.perf_mode,
        trial_full_scan,
    );
    if !options.skip_preprocessing && options.use_fermat {
        let timer = Timer::new(timings, "fermat", None);
        let slice = fermat_slice;
        if let Some(divisor) = fermat_probe(n, slice) {
            let cofactor = n / &divisor;
            outcome.record("fermat");
            let elapsed = timer.finish(format!("divisor={divisor}"));
            bill_cpu(timings, options, elapsed)?;
            return factor_recursive(
                &divisor,
                options,
                cancelled,
                deadline,
                rng,
                factors,
                outcome,
                polynomial,
                root,
                relation_count,
                rank,
                progress,
                timings,
            )
            .and_then(|_| {
                factor_recursive(
                    &cofactor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
            });
        }
        let elapsed = timer.finish("no divisor in slice".to_string());
        bill_cpu(timings, options, elapsed)?;
    }
    if !options.skip_preprocessing && options.use_pm1 {
        let timer = Timer::new(timings, "pm1", None);
        let slice = pm1_slice;
        if let Some(divisor) = crate::pm1::pollard_pm1(n, stage_one, slice) {
            let cofactor = n / &divisor;
            outcome.record("pm1");
            let elapsed = timer.finish(format!("divisor={divisor}"));
            bill_cpu(timings, options, elapsed)?;
            return factor_recursive(
                &divisor,
                options,
                cancelled,
                deadline,
                rng,
                factors,
                outcome,
                polynomial,
                root,
                relation_count,
                rank,
                progress,
                timings,
            )
            .and_then(|_| {
                factor_recursive(
                    &cofactor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
            });
        }
        let elapsed = timer.finish("no divisor".to_string());
        bill_cpu(timings, options, elapsed)?;
    }
    if !options.skip_preprocessing && options.use_pp1 {
        let timer = Timer::new(timings, "pp1", None);
        let slice = pp1_slice;
        if let Some(divisor) = crate::pm1::williams_pp1(n, crate::pm1::STAGE_ONE_BOUND, slice) {
            let cofactor = n / &divisor;
            outcome.record("pp1");
            let elapsed = timer.finish(format!("divisor={divisor}"));
            bill_cpu(timings, options, elapsed)?;
            return factor_recursive(
                &divisor,
                options,
                cancelled,
                deadline,
                rng,
                factors,
                outcome,
                polynomial,
                root,
                relation_count,
                rank,
                progress,
                timings,
            )
            .and_then(|_| {
                factor_recursive(
                    &cofactor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
            });
        }
        let elapsed = timer.finish("no divisor".to_string());
        bill_cpu(timings, options, elapsed)?;
    }
    if !options.skip_preprocessing && options.use_ecm {
        let timer = Timer::new(timings, "ecm_cheap", None);
        let cheap_remaining = signals.remaining_ms(deadline);
        let live_budget = crate::cascade::cheap_slice_ms(cheap_remaining);
        cheap_budget = live_budget;
        let cheap_slice = budget_to_deadline(live_budget);
        let (cheap_one, cheap_two, cheap_curves) = crate::cascade::cheap_ecm_plan(cheap_remaining);
        let cheap = crate::ecm::ecm_cheap_report(
            n,
            cheap_one,
            cheap_two,
            cheap_curves,
            cheap_slice,
            workers,
        );
        if let Some(divisor) = cheap.divisor {
            let cofactor = n / &divisor;
            outcome.record("ecm");
            let elapsed = timer.finish(format!(
                "divisor={divisor} curves={} ecm_ms={}",
                cheap.curves_run, cheap.elapsed_ms
            ));
            bill_cpu(timings, options, elapsed)?;
            return factor_recursive(
                &divisor,
                options,
                cancelled,
                deadline,
                rng,
                factors,
                outcome,
                polynomial,
                root,
                relation_count,
                rank,
                progress,
                timings,
            )
            .and_then(|_| {
                factor_recursive(
                    &cofactor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
            });
        }
        let elapsed = timer.finish(format!(
            "no divisor curves={} ecm_ms={}",
            cheap.curves_run, cheap.elapsed_ms
        ));
        bill_cpu(timings, options, elapsed)?;
    }
    let mut ecm_tranche_note = String::from("ecm_t1=none ecm_cap=open");
    if !options.skip_preprocessing && options.use_ecm {
        let timer = Timer::new(timings, "ecm", None);
        let joint = crate::cascade::run_joint_ecm_siqs(
            n,
            &mut signals,
            options,
            cancelled,
            deadline,
            timings,
            None,
            crate::siqs::SiqsBank::new(),
        );
        let t1_text = joint
            .ecm_t1
            .map_or_else(|| "none".to_string(), |value| value.to_string());
        ecm_tranche_note = match joint.ecm_cap {
            Some(cap) => format!("ecm_t1={t1_text} ecm_cap={cap}"),
            None => format!("ecm_t1={t1_text} ecm_cap=open"),
        };
        let ecm_box_ratio = match joint.ecm_box_ms {
            Some(box_ms) if box_ms > 0 => {
                format!(" ecm_box_ratio={:.2}", joint.ecm_ms as f64 / box_ms as f64)
            }
            _ => String::new(),
        };
        if let Some(collector) = timings {
            let relations_needed = crate::cascade::estimate_siqs_relations_needed(n, options);
            let share_text = match crate::cascade::governor_slices(
                &signals,
                signals.target_digits,
                relations_needed,
            ) {
                Some((ecm_share, siqs_share)) => {
                    format!("ecm_share={ecm_share:.4} siqs_share={siqs_share:.4}")
                }
                None => "shares=fallback".to_string(),
            };
            collector.record_elapsed(
                "joint",
                None,
                0,
                format!(
                    "{} ecm_workers={} siqs_workers={} fallback={} winner={} banked={} next_poly={} ecm_resume={}:{}",
                    share_text,
                    joint.ecm_workers,
                    joint.siqs_workers,
                    joint.used_fallback,
                    joint.winner.unwrap_or("none"),
                    joint.siqs_bank.len(),
                    joint.siqs_bank.next_poly,
                    joint.ecm_resume.map_or(0, |(digits, _)| digits),
                    joint.ecm_resume.map_or(0, |(_, curve)| curve),
                ),
            );
        }
        if let Some(divisor) = joint.divisor {
            let cofactor = n / &divisor;
            outcome.record(joint.winner.unwrap_or("ecm"));
            let elapsed = timer.finish(format!(
                "divisor={divisor} curves={} ecm_ms={} stage1_ms={} stage2_ms={}{ecm_box_ratio}",
                joint.ecm_curves, joint.ecm_ms, joint.ecm_stage_one_ms, joint.ecm_stage_two_ms
            ));
            bill_cpu(timings, options, elapsed)?;
            return factor_recursive(
                &divisor,
                options,
                cancelled,
                deadline,
                rng,
                factors,
                outcome,
                polynomial,
                root,
                relation_count,
                rank,
                progress,
                timings,
            )
            .and_then(|_| {
                factor_recursive(
                    &cofactor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
            });
        }
        let elapsed = timer.finish(format!(
            "no divisor curves={} ecm_ms={} stage1_ms={} stage2_ms={}{ecm_box_ratio}",
            joint.ecm_curves, joint.ecm_ms, joint.ecm_stage_one_ms, joint.ecm_stage_two_ms
        ));
        bill_cpu(timings, options, elapsed)?;
    }
    {
        let timer = Timer::new(timings, "cascade", None);
        let remaining = match signals.remaining_ms(deadline) {
            Some(ms) => ms.to_string(),
            None => "none".to_string(),
        };
        let ecm_ms_per_curve = match signals.measured_ms_per_curve {
            Some(ms) => ms.to_string(),
            None => "none".to_string(),
        };
        let ecm_projected_ms = match crate::cascade::project_ecm_cost_ms(&signals) {
            Some(ms) => ms.to_string(),
            None => "none".to_string(),
        };
        timer.finish(format!(
            "bits={} target_digits={} remaining_ms={} ecm_ms_per_curve={} ecm_failed_curves={} ecm_projected_ms={} slices={}/{}/{}/{}/{} {ecm_tranche_note} snfs={snfs_note}{}",
            signals.bits,
            signals.target_digits,
            remaining,
            ecm_ms_per_curve,
            signals.failed_curves,
            ecm_projected_ms,
            optional_u128_text(cheap_budget),
            optional_u128_text(fermat_budget),
            optional_u128_text(pm1_budget),
            optional_u128_text(pp1_budget),
            optional_u128_text(ecm_budget),
            if options.perf_mode { " perf_mode=true" } else { "" }
        ));
    }
    if options.use_gnfs {
        check_cpu(timings, options)?;
        emit(
            progress,
            "gnfs",
            1,
            7,
            "selecting a polynomial and collecting relations",
        )?;
        match gnfs_divisor(
            n,
            options,
            cancelled,
            deadline,
            progress,
            timings,
            &mut signals,
        ) {
            Ok((divisor, stats, strategy)) => {
                outcome.record(strategy);
                *relation_count += stats.relations;
                *rank += stats.rank;
                let polynomial_root = stats.polynomial.m.to_string();
                if polynomial.is_none() {
                    *polynomial = Some(stats.polynomial);
                }
                if root.is_none() {
                    *root = Some(polynomial_root);
                }
                let cofactor = n / &divisor;
                return factor_recursive(
                    &divisor,
                    options,
                    cancelled,
                    deadline,
                    rng,
                    factors,
                    outcome,
                    polynomial,
                    root,
                    relation_count,
                    rank,
                    progress,
                    timings,
                )
                .and_then(|_| {
                    factor_recursive(
                        &cofactor,
                        options,
                        cancelled,
                        deadline,
                        rng,
                        factors,
                        outcome,
                        polynomial,
                        root,
                        relation_count,
                        rank,
                        progress,
                        timings,
                    )
                });
            }
            Err(error) if options.allow_fallback => {
                if matches!(
                    error,
                    FactorError::ResourceLimit(_) | FactorError::Cancelled
                ) {
                    return Err(error);
                }
                outcome.record("fallback");
            }
            Err(error) => return Err(error),
        }
    }
    if !options.allow_fallback {
        return Err(FactorError::Factorization(
            "no factorization strategy succeeded".to_string(),
        ));
    }
    outcome.record("fallback");
    check_cpu(timings, options)?;
    emit(progress, "fallback", 3, 7, "running Pollard-Rho fallback")?;
    let divisor = pollard_rho(n, rng, cancelled, deadline)?;
    let cofactor = n / &divisor;
    factor_recursive(
        &divisor,
        options,
        cancelled,
        deadline,
        rng,
        factors,
        outcome,
        polynomial,
        root,
        relation_count,
        rank,
        progress,
        timings,
    )?;
    factor_recursive(
        &cofactor,
        options,
        cancelled,
        deadline,
        rng,
        factors,
        outcome,
        polynomial,
        root,
        relation_count,
        rank,
        progress,
        timings,
    )
}

fn optional_u128_text(value: Option<u128>) -> String {
    match value {
        Some(number) => number.to_string(),
        None => "none".to_string(),
    }
}

fn optional_rate_text(value: Option<f64>) -> String {
    match value {
        Some(rate) => format!("{rate:.4}"),
        None => "none".to_string(),
    }
}

fn candidate_degrees(n: &Integer) -> Vec<usize> {
    let bits = n.bits();
    let mut candidates = Vec::new();
    if bits <= 32 {
        return vec![2];
    }
    candidates.push(2);
    if bits >= 36 {
        candidates.push(3);
    }
    if bits >= 44 {
        candidates.push(4);
    }
    if bits >= 60 {
        candidates.push(5);
    }
    candidates
}

const SCALED_RATIONAL_CAP: usize = 50_000;
const SCALED_ALGEBRAIC_CAP: usize = 250_000;
const SCALED_SIEVE_A_CAP: usize = 500_000;
const SCALED_SIEVE_B_CAP: usize = 60_000;
const SCALED_LARGE_PRIME_CAP: usize = 1_000_000;

fn scaled_anchor(degree: usize) -> Option<(u64, usize, usize, usize, usize)> {
    match degree {
        3 => Some((27, 2_000, 500, 1_600, 700)),
        4 => Some((27, 2_000, 2_000, 6_000, 1_500)),
        5 => Some((43, 4_000, 20_000, 72_000, 8_000)),
        _ => None,
    }
}

fn scaled_growth(anchor_bits: u64, bits: u64, step: u32) -> usize {
    let delta = bits.saturating_sub(anchor_bits);
    let shift = u64::from(step).saturating_mul(delta) / 32;
    1_usize << shift.min(20)
}

fn resolve_scaled_bounds(n: &Integer, degree: usize, options: &mut FactorOptions) {
    let Some((anchor_bits, rational, algebraic, sieve_a, sieve_b)) = scaled_anchor(degree) else {
        return;
    };
    let bits = n.bits();
    if bits < anchor_bits {
        return;
    }
    let rational_growth = scaled_growth(anchor_bits, bits, 1);
    let algebraic_growth = scaled_growth(anchor_bits, bits, 5);
    let sieve_a_growth = scaled_growth(anchor_bits, bits, 3);
    let sieve_b_growth = scaled_growth(anchor_bits, bits, 2);
    let large_prime_growth = scaled_growth(anchor_bits, bits, 1);
    let defaults = FactorOptions::default();
    if options.rational_base_bound == defaults.rational_base_bound {
        options.rational_base_bound = rational
            .saturating_mul(rational_growth)
            .min(SCALED_RATIONAL_CAP);
    }
    if options.algebraic_base_bound == defaults.algebraic_base_bound {
        options.algebraic_base_bound = algebraic
            .saturating_mul(algebraic_growth)
            .min(SCALED_ALGEBRAIC_CAP);
    }
    if options.sieve_a_bound == defaults.sieve_a_bound {
        options.sieve_a_bound = sieve_a
            .saturating_mul(sieve_a_growth)
            .min(SCALED_SIEVE_A_CAP);
    }
    if options.sieve_b_bound == defaults.sieve_b_bound {
        options.sieve_b_bound = sieve_b
            .saturating_mul(sieve_b_growth)
            .min(SCALED_SIEVE_B_CAP);
    }
    if options.large_prime_bound == defaults.large_prime_bound {
        options.large_prime_bound = 50_000_usize
            .saturating_mul(large_prime_growth)
            .min(SCALED_LARGE_PRIME_CAP);
    }
}

fn preflight_collect_exceeds(n: &Integer, degree: usize, options: &FactorOptions) -> bool {
    if degree < 3 {
        return false;
    }
    let mut attempt = options.clone();
    resolve_scaled_bounds(n, degree, &mut attempt);
    let Ok(max_a) = i64::try_from(attempt.sieve_a_bound) else {
        return false;
    };
    let workers = attempt.effective_workers();
    crate::options::estimate_collect_bytes(1, max_a, workers) > attempt.memory_budget_bytes()
}

#[allow(clippy::too_many_arguments)]
fn run_degree<F>(
    n: &Integer,
    options: &FactorOptions,
    degree: usize,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
    timings: Option<&TimingCollector>,
    pilot: bool,
) -> Result<(Integer, GnfsStats, &'static str), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let mut attempt = options.clone();
    attempt.degree = Some(degree);
    resolve_scaled_bounds(n, degree, &mut attempt);
    if degree == 2 {
        quadratic_divisor(n, &attempt, cancelled, deadline, progress, timings, pilot)
    } else {
        match general::gnfs_divisor(n, &attempt, cancelled, deadline, timings) {
            Ok((divisor, polynomial, relations, rank, entries)) => {
                if let Some(collector) = timings {
                    collector.merge(entries);
                }
                match emit(progress, "gnfs", 7, 7, "verified a congruence of squares") {
                    Ok(()) => Ok((
                        divisor,
                        GnfsStats {
                            polynomial,
                            relations,
                            rank,
                        },
                        "gnfs",
                    )),
                    Err(error) => Err(error),
                }
            }
            Err(error) => Err(error),
        }
    }
}

fn pilot_error_kind(
    outcome: &Result<(Integer, GnfsStats, &'static str), FactorError>,
) -> &'static str {
    match outcome {
        Ok(_) => "none",
        Err(FactorError::ResourceLimit(_)) => "ResourceLimit",
        Err(FactorError::Cancelled) => "Cancelled",
        Err(FactorError::InsufficientRelations { .. }) => "InsufficientRelations",
        Err(FactorError::NoPolynomial) => "NoPolynomial",
        Err(FactorError::SearchExhausted) => "SearchExhausted",
        Err(FactorError::AlgebraicSquareRootFailed) => "AlgebraicSquareRootFailed",
        Err(FactorError::LinearAlgebra(_)) => "LinearAlgebra",
        Err(FactorError::Factorization(_)) => "Factorization",
        Err(FactorError::InvalidInput(_)) => "InvalidInput",
        Err(FactorError::UnsupportedDegree(_)) => "UnsupportedDegree",
    }
}

type RatedPilot = (usize, f64, Option<u128>);

fn measured_select_ms(timings: Option<&TimingCollector>, degree: usize) -> Option<u128> {
    let collector = timings?;
    collector
        .snapshot()
        .iter()
        .rev()
        .find(|entry| entry.stage == "select_polynomial" && entry.degree == Some(degree))
        .map(|entry| entry.elapsed_ms)
}

#[allow(clippy::too_many_arguments)]
fn note_pilot_for_escalation(
    best_rated: Option<RatedPilot>,
    first_limited: Option<usize>,
    degree: usize,
    pilot_yield: Option<usize>,
    pilot_ms: u128,
    timed_out: bool,
    relations_needed: Option<usize>,
    remaining_ms: u128,
    select_ms: Option<u128>,
) -> (Option<RatedPilot>, Option<usize>) {
    let mut best = best_rated;
    let mut limited = first_limited;
    let rate = pilot_yield.map(|found| found as f64 / pilot_ms.max(1) as f64);
    let projected = match (relations_needed, rate) {
        (Some(needed), Some(rate)) if degree >= 3 => {
            crate::cascade::project_gnfs_cost_ms(needed, rate, select_ms)
        }
        (Some(needed), Some(rate)) => crate::cascade::project_siqs_cost_ms(needed, rate),
        _ => None,
    };
    let affordable = projected.is_none_or(|cost| cost <= remaining_ms);
    if timed_out && limited.is_none() && pilot_yield != Some(0) && affordable {
        limited = Some(degree);
    }
    if let Some(rate) = rate {
        best = better_pilot_pick(best, degree, rate, projected, remaining_ms);
    }
    (best, limited)
}

fn better_pilot_pick(
    current: Option<RatedPilot>,
    degree: usize,
    rate: f64,
    projected_cost_ms: Option<u128>,
    remaining_ms: u128,
) -> Option<RatedPilot> {
    if rate <= 0.0 {
        return current;
    }
    if projected_cost_ms.is_some_and(|cost| cost > remaining_ms) {
        return current;
    }
    let cheaper = match (current, projected_cost_ms) {
        (Some((_, _, Some(best))), Some(cost)) => cost < best,
        (Some((_, _, None)), Some(_)) => true,
        _ => current.is_none_or(|(_, best, _)| rate > best),
    };
    if cheaper {
        Some((degree, rate, projected_cost_ms))
    } else {
        current
    }
}

fn escalation_pick(
    best_rated: Option<RatedPilot>,
    first_limited: Option<usize>,
) -> Option<(usize, &'static str, Option<f64>)> {
    best_rated
        .map(|(degree, rate, _)| (degree, "best_rated", Some(rate)))
        .or(first_limited.map(|degree| (degree, "first_limited", None)))
}

fn unaffordable_pilot_cost(
    needed: Option<usize>,
    rate: Option<f64>,
    select_ms: Option<u128>,
    remaining_ms: u128,
) -> Option<u128> {
    let (Some(needed), Some(rate)) = (needed, rate) else {
        return None;
    };
    let cost = crate::cascade::project_gnfs_cost_ms(needed, rate, select_ms)?;
    (cost > remaining_ms).then_some(cost)
}

fn escalation_has_no_rated_pick(
    best_rated: Option<RatedPilot>,
    first_limited: Option<usize>,
) -> bool {
    best_rated.is_none() && first_limited.is_some()
}

fn probe_budget(n: &Integer) -> Duration {
    let ln_n = (n.bits() as f64 * std::f64::consts::LN_2).max(1.0);
    let ln_ln_n = ln_n.ln().max(f64::EPSILON);
    let exponent = ln_n.powf(1.0 / 3.0) * ln_ln_n.powf(2.0 / 3.0) * (64.0f64 / 9.0).powf(1.0 / 3.0);
    let millis = exponent.exp().clamp(1_000.0, 3_600_000.0);
    Duration::from_millis(millis as u64)
}

fn gnfs_divisor<F>(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
    timings: Option<&TimingCollector>,
    signals: &mut crate::cascade::CascadeSignals,
) -> Result<(Integer, GnfsStats, &'static str), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    if let Some(degree) = options.degree {
        return run_degree(
            n, options, degree, cancelled, deadline, progress, timings, false,
        );
    }
    let candidates = candidate_degrees(n);
    if candidates.len() == 1 {
        return run_degree(
            n,
            options,
            candidates[0],
            cancelled,
            deadline,
            progress,
            timings,
            false,
        );
    }
    let fraction = options.degree_race_fraction.clamp(0.05, 0.5);
    let bound = deadline.unwrap_or_else(|| Instant::now() + probe_budget(n));
    let mut first_limited: Option<usize> = None;
    let mut best_rated: Option<RatedPilot> = None;
    let mut floor_required: Option<usize> = None;
    let mut floor_rate: Option<f64> = None;
    for &degree in candidates.iter() {
        if cancelled.load(Ordering::Relaxed) {
            return Err(FactorError::Cancelled);
        }
        let now = Instant::now();
        if now >= bound {
            break;
        }
        let pilot_start = Instant::now();
        if degree >= 3 {
            let remaining_ms = bound.saturating_duration_since(pilot_start).as_millis();
            if let Some(cost) =
                unaffordable_pilot_cost(floor_required, floor_rate, None, remaining_ms)
            {
                if let Some(collector) = timings {
                    collector.record(
                        "race_pilot",
                        Some(degree),
                        pilot_start,
                        format!(
                            "slice_ms=0 skipped relations=none err=Skipped projected_gnfs_ms={cost} remaining_ms={remaining_ms}"
                        ),
                    );
                }
                continue;
            }
        }
        let (outcome, slice_ms) = if degree >= 3 && preflight_collect_exceeds(n, degree, options) {
            (
                Err(general::collect_memory_error(options.memory_budget_bytes())),
                0,
            )
        } else {
            let slice = now + (bound - now).mul_f64(fraction);
            let slice_ms = slice.saturating_duration_since(now).as_millis();
            let outcome = run_degree(
                n,
                options,
                degree,
                cancelled,
                Some(slice),
                progress,
                timings,
                true,
            );
            (outcome, slice_ms)
        };
        let pilot_ms = pilot_start.elapsed().as_millis();
        let pilot_yield = match &outcome {
            Ok((_, stats, _)) => Some(stats.relations),
            Err(FactorError::InsufficientRelations { found, .. }) => Some(*found),
            _ => None,
        };
        let relations_needed = match &outcome {
            Err(FactorError::InsufficientRelations { required, .. }) => Some(*required),
            _ => None,
        };
        let select_ms = measured_select_ms(timings, degree);
        if degree >= 3 {
            if let Some(required) = relations_needed {
                floor_required = Some(floor_required.map_or(required, |floor| floor.min(required)));
            }
            if let Some(found) = pilot_yield {
                let rate = found as f64 / pilot_ms.max(1) as f64;
                if rate > 0.0 {
                    floor_rate = Some(floor_rate.map_or(rate, |floor| floor.max(rate)));
                }
            }
        }
        if let Some(found) = pilot_yield {
            let rate = found as f64 / pilot_ms.max(1) as f64;
            if let Some(collector) = timings {
                collector.record_elapsed(
                    "pilot_rate",
                    Some(degree),
                    pilot_ms,
                    format!("relations={found} rate={rate:.4}"),
                );
            }
            if rate > 0.0
                && signals
                    .measured_relations_per_sec()
                    .is_none_or(|best| rate * 1000.0 > best)
            {
                signals.report_pilot(found, pilot_ms);
            }
        }
        if let Some(collector) = timings {
            let pilot_rate = pilot_yield.map(|found| found as f64 / pilot_ms.max(1) as f64);
            let cost_ecm_ms = crate::cascade::project_ecm_cost_ms(signals);
            let cost_ecm_next_ms = crate::cascade::project_ecm_cost_ms_for_digits(
                signals,
                crate::cascade::escalate_digits(signals.target_digits),
            );
            let cost_siqs_ms = match (relations_needed, pilot_rate) {
                (Some(needed), Some(rate)) => crate::cascade::project_siqs_cost_ms(needed, rate),
                _ => None,
            };
            let cost_gnfs_ms = match (relations_needed, pilot_rate) {
                (Some(needed), Some(rate)) if degree >= 3 => {
                    crate::cascade::project_gnfs_cost_ms(needed, rate, select_ms)
                }
                _ => None,
            };
            let lower_cost = crate::cascade::projected_lower_cost(cost_ecm_next_ms, cost_siqs_ms);
            collector.record_elapsed(
                "cascade_projection",
                Some(degree),
                0,
                format!(
                    "cost_ecm_ms={} cost_ecm_next_ms={} cost_siqs_ms={} relations_needed={} pilot_rate={} lower_cost={} cost_gnfs_ms={}",
                    optional_u128_text(cost_ecm_ms),
                    optional_u128_text(cost_ecm_next_ms),
                    optional_u128_text(cost_siqs_ms),
                    optional_u128_text(relations_needed.map(|needed| needed as u128)),
                    optional_rate_text(pilot_rate),
                    lower_cost.unwrap_or("none"),
                    optional_u128_text(cost_gnfs_ms),
                ),
            );
            if let Some(needed) = relations_needed {
                let _governor = crate::cascade::governor_slices(
                    signals,
                    crate::cascade::escalate_digits(signals.target_digits),
                    needed,
                );
            }
        }
        if let Some(collector) = timings {
            let status = match &outcome {
                Ok(_) => "won".to_string(),
                Err(FactorError::ResourceLimit(_)) => "out of time".to_string(),
                Err(FactorError::Cancelled) => "cancelled".to_string(),
                Err(_) => "failed".to_string(),
            };
            let elapsed = collector.record(
                "race_pilot",
                Some(degree),
                now,
                format!(
                    "slice_ms={slice_ms} {status} relations={} err={}",
                    pilot_yield.map_or("none".to_string(), |found| found.to_string()),
                    pilot_error_kind(&outcome),
                ),
            );
            bill_cpu(timings, options, elapsed)?;
        }
        match outcome {
            Ok(result) => return Ok(result),
            Err(FactorError::Cancelled) => return Err(FactorError::Cancelled),
            Err(FactorError::ResourceLimit(_)) => {
                (best_rated, first_limited) = note_pilot_for_escalation(
                    best_rated,
                    first_limited,
                    degree,
                    pilot_yield,
                    pilot_ms,
                    true,
                    relations_needed,
                    bound.saturating_duration_since(Instant::now()).as_millis(),
                    select_ms,
                );
            }
            Err(_) => {
                (best_rated, first_limited) = note_pilot_for_escalation(
                    best_rated,
                    first_limited,
                    degree,
                    pilot_yield,
                    pilot_ms,
                    false,
                    relations_needed,
                    bound.saturating_duration_since(Instant::now()).as_millis(),
                    select_ms,
                );
            }
        }
    }
    let escalation = escalation_pick(best_rated, first_limited);
    match escalation {
        Some((degree, pick, rate)) => {
            if escalation_has_no_rated_pick(best_rated, first_limited) {
                if let Some(collector) = timings {
                    collector.record(
                        "race_escalation",
                        Some(degree),
                        Instant::now(),
                        format!(
                            "full budget skipped pick={pick} rate={} failed_degree={degree}",
                            optional_rate_text(rate)
                        ),
                    );
                }
                return Err(FactorError::ResourceLimit(format!(
                    "GNFS race escalation skipped: degree {degree} pilot failed with no rated affordable pick"
                )));
            }
            let started = Instant::now();
            let result = run_degree(
                n, options, degree, cancelled, deadline, progress, timings, false,
            );
            if let Some(collector) = timings {
                let status = if result.is_ok() { "won" } else { "failed" };
                let elapsed = collector.record(
                    "race_escalation",
                    Some(degree),
                    started,
                    format!(
                        "full budget {status} pick={pick} rate={}",
                        optional_rate_text(rate)
                    ),
                );
                bill_cpu(timings, options, elapsed)?;
            }
            result
        }
        None => Err(FactorError::SearchExhausted),
    }
}

#[allow(clippy::too_many_arguments)]
fn quadratic_divisor<F>(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
    timings: Option<&TimingCollector>,
    pilot: bool,
) -> Result<(Integer, GnfsStats, &'static str), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let total = Instant::now();
    if options.use_mpqs {
        let siqs_deadline = deadline.map(|limit| {
            let remaining = limit.saturating_duration_since(Instant::now());
            Instant::now() + remaining / 4
        });
        let siqs_timer = Timer::new(timings, "siqs", Some(2));
        match crate::siqs::siqs_divisor(n, options, cancelled, siqs_deadline, progress, timings) {
            Ok((divisor, outcome)) => {
                siqs_timer.finish(format!(
                    "relations={} rank={}",
                    outcome.relations, outcome.rank
                ));
                return Ok((
                    divisor,
                    GnfsStats {
                        polynomial: outcome.display,
                        relations: outcome.relations,
                        rank: outcome.rank,
                    },
                    "siqs",
                ));
            }
            Err(err) => {
                siqs_timer.finish(format!("miss: {err}"));
                if matches!(err, FactorError::Cancelled) {
                    return Err(err);
                }
            }
        }
    }
    let polynomial = {
        let timer = Timer::new(timings, "select_polynomial", Some(2));
        let value = select_polynomial(n, options)?;
        timer.finish(format!("m={}", value.m));
        value
    };
    emit(progress, "gnfs", 2, 7, "polynomial selected")?;
    let mut working = options.clone();
    let mut last_error = FactorError::SearchExhausted;
    let attempts: u32 = if pilot { 1 } else { QUADRATIC_ATTEMPTS };
    for attempt in 0..attempts {
        if cancelled.load(Ordering::Relaxed) {
            return Err(FactorError::Cancelled);
        }
        if deadline.is_some_and(|deadline| Instant::now() >= deadline) {
            return Err(FactorError::ResourceLimit(
                "GNFS deadline exceeded".to_string(),
            ));
        }
        check_cpu(timings, options)?;
        let base = {
            let timer = Timer::new(timings, "factor_base", Some(2));
            let value = FactorBase::build(
                &polynomial,
                n,
                working.rational_base_bound,
                working.algebraic_base_bound,
            )?;
            timer.finish(format!(
                "attempt={attempt} columns={}",
                value.column_count()
            ));
            value
        };
        emit(progress, "gnfs", 3, 7, "factor base built")?;
        let relations = {
            let timer = Timer::new(timings, "relation_collection", Some(2));
            match collect_relations(
                n,
                &polynomial,
                &base,
                &working,
                cancelled,
                deadline,
                progress,
                timings,
            ) {
                Ok(relations) => {
                    timer.finish(format!("attempt={attempt} relations={}", relations.len()));
                    relations
                }
                Err(error) => {
                    timer.finish(format!("attempt={attempt} failed {error}"));
                    last_error = error;
                    if !expand_search(&mut working) {
                        break;
                    }
                    continue;
                }
            }
        };
        emit(progress, "gnfs", 4, 7, "relations collected")?;
        let workers = options.effective_workers();
        let (rank, dependencies) = {
            let timer = Timer::new(timings, "matrix_rank", Some(2));
            let value = rank_and_dependencies(&relations, base.column_count(), 128, 1, workers);
            timer.finish(format!(
                "attempt={attempt} rank={} of {}",
                value.0,
                relations.len()
            ));
            value
        };
        emit(progress, "gnfs", 5, 7, "matrix rank computed")?;
        if relations.len() <= rank {
            last_error = FactorError::InsufficientRelations {
                found: relations.len(),
                required: rank + 1,
            };
            if !expand_search(&mut working) {
                break;
            }
            continue;
        }
        {
            let timer = Timer::new(timings, "dependency_search", Some(2));
            timer.finish(format!(
                "attempt={attempt} dependencies={}",
                dependencies.len()
            ));
        }
        emit(progress, "gnfs", 6, 7, "relation dependencies found")?;
        let ring = QuadraticRing::new(polynomial.constant());
        let attempt_order = crate::linear_algebra::order_dependencies_by_cost(&dependencies);
        for (di, index) in attempt_order.iter().enumerate() {
            let timer = Timer::new(timings, "square_root", Some(2));
            let attempt_result =
                extract_factor(n, &ring, &polynomial.m, &relations, &dependencies[*index]);
            let found = attempt_result.is_ok();
            timer.finish(format!("dependency={di} found={found}"));
            if let Ok(divisor) = attempt_result {
                emit(progress, "gnfs", 7, 7, "verified a congruence of squares")?;
                if let Some(collector) = timings {
                    collector.record(
                        "total",
                        Some(2),
                        total,
                        format!(
                            "relations={} rank={} dependencies={}",
                            relations.len(),
                            rank,
                            dependencies.len()
                        ),
                    );
                }
                return Ok((
                    divisor,
                    GnfsStats {
                        polynomial,
                        relations: relations.len(),
                        rank,
                    },
                    "quadratic-qs",
                ));
            }
        }
        last_error = FactorError::AlgebraicSquareRootFailed;
        if !expand_search(&mut working) {
            break;
        }
    }
    match last_error {
        FactorError::Cancelled => Err(last_error),
        _ if pilot => Err(FactorError::ResourceLimit(
            "quadratic pilot yielded to the degree race".to_string(),
        )),
        _ => Err(last_error),
    }
}

const QUADRATIC_ATTEMPTS: u32 = 6;

fn expand_search(options: &mut FactorOptions) -> bool {
    let mut grew = false;
    if options.rational_base_bound < SCALED_RATIONAL_CAP {
        options.rational_base_bound = (options.rational_base_bound * 2).min(SCALED_RATIONAL_CAP);
        grew = true;
    }
    if options.algebraic_base_bound < SCALED_ALGEBRAIC_CAP {
        options.algebraic_base_bound = (options.algebraic_base_bound * 2).min(SCALED_ALGEBRAIC_CAP);
        grew = true;
    }
    if options.sieve_a_bound < SCALED_SIEVE_A_CAP {
        options.sieve_a_bound = (options.sieve_a_bound * 2).min(SCALED_SIEVE_A_CAP);
        grew = true;
    }
    if options.sieve_b_bound < SCALED_SIEVE_B_CAP {
        options.sieve_b_bound = (options.sieve_b_bound * 2).min(SCALED_SIEVE_B_CAP);
        grew = true;
    }
    grew
}

fn trial_division(n: &Integer) -> (Option<Integer>, u64) {
    let digits: Vec<u64> = n.magnitude().iter_u64_digits().collect();
    let small_n = n.to_u64_saturating();
    let mut checked = 0_u64;
    for prime in crate::prime::small_primes() {
        checked += 1;
        let mut residue: u128 = 0;
        for limb in digits.iter().rev() {
            residue = ((residue << 64) | *limb as u128) % u128::from(*prime);
        }
        if residue == 0 {
            return (Some(Integer::from(*prime)), checked);
        }
        if prime
            .checked_mul(*prime)
            .is_some_and(|square| square > small_n)
        {
            break;
        }
    }
    (None, checked)
}

fn perfect_power_root(n: &Integer) -> Option<(Integer, u32)> {
    if n <= &Integer::from(1) {
        return None;
    }
    let mut exponent = 3_u64;
    while exponent <= n.bits() {
        let mut root = crate::general::integer_root(n, exponent as usize);
        while root.pow(exponent as u32) > *n {
            root -= Integer::from(1);
        }
        while (&root + Integer::from(1)).pow(exponent as u32) <= *n {
            root += Integer::from(1);
        }
        if root.pow(exponent as u32) == *n {
            return Some((root, exponent as u32));
        }
        exponent += 1;
    }
    None
}

pub(crate) struct SpecialForm {
    pub(crate) a: u64,
    pub(crate) b: u64,
    pub(crate) k: u32,
    pub(crate) c: u64,
    pub(crate) sign: i8,
    pub(crate) low_weight: bool,
}

impl SpecialForm {
    pub(crate) fn describe(&self) -> String {
        if self.k == 0 {
            return "heuristic,low_weight=true".to_string();
        }
        format!(
            "a={},b={},k={},c={},sign={},low_weight={}",
            self.a, self.b, self.k, self.c, self.sign, self.low_weight
        )
    }
}

const SPECIAL_BASE_LIMIT: u64 = 1_000_000;
const LOW_HAMMING_WEIGHT_LIMIT: u64 = 12;
const SMALL_COEFFICIENT_SCAN_MIN_BITS: u64 = 128;
const SMALL_COEFFICIENT_SCAN_MAX: u64 = 16;

pub(crate) fn popcount_weight(n: &Integer) -> u64 {
    let mut total = 0_u64;
    for limb in n.magnitude().iter_u64_digits() {
        total += u64::from(limb.count_ones());
    }
    total
}

pub(crate) fn detect_special_form(n: &Integer) -> Option<SpecialForm> {
    if n <= &Integer::from(3) {
        return None;
    }
    let low_weight = popcount_weight(n) <= LOW_HAMMING_WEIGHT_LIMIT;
    if let Some(mut form) = exact_unit_coefficient_form(n) {
        form.low_weight = low_weight;
        return Some(form);
    }
    if n.bits() >= SMALL_COEFFICIENT_SCAN_MIN_BITS
        && let Some(mut form) = exact_small_coefficient_form(n)
    {
        form.low_weight = low_weight;
        return Some(form);
    }
    if low_weight {
        return Some(SpecialForm {
            a: 0,
            b: 0,
            k: 0,
            c: 0,
            sign: 0,
            low_weight: true,
        });
    }
    None
}

fn exact_unit_coefficient_form(n: &Integer) -> Option<SpecialForm> {
    let one = Integer::from(1);
    if let Some((base, exponent)) = unit_power_params(&(n - &one)) {
        return Some(SpecialForm {
            a: 1,
            b: base,
            k: exponent,
            c: 1,
            sign: 1,
            low_weight: false,
        });
    }
    if let Some((base, exponent)) = unit_power_params(&(n + &one)) {
        return Some(SpecialForm {
            a: 1,
            b: base,
            k: exponent,
            c: 1,
            sign: -1,
            low_weight: false,
        });
    }
    None
}

fn exact_small_coefficient_form(n: &Integer) -> Option<SpecialForm> {
    let one = Integer::from(1);
    let mut coefficient = 2_u64;
    while coefficient <= SMALL_COEFFICIENT_SCAN_MAX {
        let divisor = Integer::from(coefficient);
        if (n - &one) % &divisor == Integer::from(0)
            && let Some((base, exponent)) = unit_power_params(&((n - &one) / &divisor))
        {
            return Some(SpecialForm {
                a: coefficient,
                b: base,
                k: exponent,
                c: 1,
                sign: 1,
                low_weight: false,
            });
        }
        if (n + &one) % &divisor == Integer::from(0)
            && let Some((base, exponent)) = unit_power_params(&((n + &one) / &divisor))
        {
            return Some(SpecialForm {
                a: coefficient,
                b: base,
                k: exponent,
                c: 1,
                sign: -1,
                low_weight: false,
            });
        }
        coefficient += 1;
    }
    None
}

fn unit_power_params(target: &Integer) -> Option<(u64, u32)> {
    if target < &Integer::from(4) {
        return None;
    }
    let mut degree = target.bits();
    let floor = target.bits().div_ceil(20).max(3);
    while degree >= floor {
        let root = snfs_floor_root(target, degree);
        if root.pow(degree as u32) == *target {
            let base = root.to_u64_saturating();
            if (2..=SPECIAL_BASE_LIMIT).contains(&base) {
                return Some((base, degree as u32));
            }
            return None;
        }
        if degree == floor {
            break;
        }
        degree -= 1;
    }
    let square = crate::integer::integer_sqrt(target);
    if &square * &square == *target {
        let base = square.to_u64_saturating();
        if (2..=SPECIAL_BASE_LIMIT).contains(&base) {
            return Some((base, 2));
        }
    }
    None
}

fn snfs_floor_root(value: &Integer, degree: u64) -> Integer {
    if value <= &Integer::from(1) {
        return value.clone();
    }
    let one = Integer::from(1);
    let mut guess = Integer::from(1) << value.bits().div_ceil(degree);
    loop {
        let previous = guess.pow(degree as u32 - 1);
        let next = (&guess * Integer::from(degree - 1) + value / &previous) / Integer::from(degree);
        if next >= guess {
            break;
        }
        guess = next;
    }
    while (&guess + &one).pow(degree as u32) <= *value {
        guess += &one;
    }
    while guess.pow(degree as u32) > *value {
        guess -= &one;
    }
    guess
}

fn fermat_probe(n: &Integer, deadline: Option<Instant>) -> Option<Integer> {
    if n <= &Integer::from(2) {
        return None;
    }
    let cap = crate::cascade::fermat_cap_for(n);
    let mut base = crate::integer::integer_sqrt(n);
    if &base * &base < *n {
        base += Integer::from(1);
    }
    let block_start = Instant::now();
    let mut iterations = 0_u64;
    let mut iter_cap = cap;
    let mut calibrated = deadline.is_none();
    loop {
        if iterations & 0x3ff == 0 {
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            if !calibrated && iterations > 0 {
                let elapsed_ns = block_start.elapsed().as_nanos();
                let ns_per_iter = elapsed_ns / u128::from(iterations);
                if let Some(limit) = deadline {
                    let remaining_ms = limit.saturating_duration_since(Instant::now()).as_millis();
                    iter_cap =
                        crate::cascade::fermat_iter_budget(cap, Some(remaining_ms), ns_per_iter);
                }
                calibrated = true;
            }
        }
        if iterations > iter_cap {
            return None;
        }
        let square = &base * &base - n;
        let root = crate::integer::integer_sqrt(&square);
        if &root * &root == square && root != Integer::from(0) {
            let divisor = &base - &root;
            if divisor > Integer::from(1) && divisor < *n {
                return Some(divisor);
            }
        }
        base += Integer::from(1);
        iterations += 1;
    }
}

fn budget_to_deadline(budget_ms: Option<u128>) -> Option<Instant> {
    let ms = budget_ms?;
    let duration = Duration::from_millis(u64::try_from(ms).ok()?);
    Instant::now().checked_add(duration)
}

fn bill_cpu(
    timings: Option<&TimingCollector>,
    options: &FactorOptions,
    elapsed_ms: u128,
) -> Result<(), FactorError> {
    if let Some(collector) = timings {
        collector.add_cpu_ms(elapsed_ms, options.effective_workers());
        collector.check_cpu_budget(options)?;
    }
    Ok(())
}

fn check_cpu(
    timings: Option<&TimingCollector>,
    options: &FactorOptions,
) -> Result<(), FactorError> {
    if let Some(collector) = timings {
        collector.check_cpu_budget(options)?;
    }
    Ok(())
}

fn perfect_square(n: &Integer) -> Option<Integer> {
    let root = crate::prime::integer_sqrt(n);
    if &root * &root == *n {
        Some(root)
    } else {
        None
    }
}

fn pollard_rho(
    n: &Integer,
    rng: &mut ChaCha8Rng,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
) -> Result<Integer, FactorError> {
    if n <= &Integer::from(3) {
        return Err(FactorError::Factorization(
            "Pollard-Rho received a non-composite value".to_string(),
        ));
    }
    if n % Integer::from(2) == Integer::from(0) {
        return Ok(Integer::from(2));
    }
    let range = n - 2;
    for _ in 0..16 {
        let constant = Integer::from(rng.random::<u64>()) % &range + 1;
        let mut outer = Integer::from(2);
        let mut radius = 1_u64;
        let mut total = 0_u64;
        let mut answer: Option<Integer> = None;
        while answer.is_none() {
            let mut inner = outer.clone();
            let mut advanced = 0_u64;
            while advanced < radius {
                inner = ((&inner * &inner) + &constant) % n;
                advanced += 1;
                total += 1;
                if total & 0x3ff == 0 {
                    if cancelled.load(Ordering::Relaxed) {
                        return Err(FactorError::Cancelled);
                    }
                    if deadline.is_some_and(|deadline| Instant::now() >= deadline) {
                        return Err(FactorError::ResourceLimit(
                            "Pollard-Rho deadline exceeded".to_string(),
                        ));
                    }
                }
                if total > 2_000_000 {
                    break;
                }
            }
            if total > 2_000_000 {
                break;
            }
            let mut product = Integer::from(1);
            let checkpoint = inner.clone();
            let mut used = 0_u64;
            let mut divisor = Integer::from(1);
            while used < radius {
                let chunk = (radius - used).min(128);
                let mut stepped = 0_u64;
                while stepped < chunk {
                    inner = ((&inner * &inner) + &constant) % n;
                    total += 1;
                    if total & 0x3ff == 0 {
                        if cancelled.load(Ordering::Relaxed) {
                            return Err(FactorError::Cancelled);
                        }
                        if deadline.is_some_and(|deadline| Instant::now() >= deadline) {
                            return Err(FactorError::ResourceLimit(
                                "Pollard-Rho deadline exceeded".to_string(),
                            ));
                        }
                    }
                    if total > 2_000_000 {
                        break;
                    }
                    let difference = if inner > outer {
                        &inner - &outer
                    } else {
                        &outer - &inner
                    };
                    product = (product * difference) % n;
                    stepped += 1;
                }
                if total > 2_000_000 {
                    break;
                }
                used += chunk;
                divisor = crate::integer::integer_gcd(&product, n);
                if divisor != Integer::from(1) {
                    break;
                }
            }
            if divisor == Integer::from(1) {
                outer = inner;
                radius *= 2;
                continue;
            }
            if divisor == *n {
                let mut walk = checkpoint;
                let mut guard = 0_u64;
                loop {
                    walk = ((&walk * &walk) + &constant) % n;
                    let difference = if walk > outer {
                        &walk - &outer
                    } else {
                        &outer - &walk
                    };
                    if difference == Integer::from(0) {
                        break;
                    }
                    divisor = crate::integer::integer_gcd(&difference, n);
                    guard += 1;
                    if divisor != Integer::from(1) || guard > radius * 2 {
                        break;
                    }
                }
            }
            if divisor > Integer::from(1) && divisor < *n {
                answer = Some(divisor);
            } else {
                break;
            }
        }
        if let Some(divisor) = answer {
            return Ok(divisor);
        }
    }
    Err(FactorError::Factorization(
        "Pollard-Rho did not find a divisor".to_string(),
    ))
}

fn emit<F>(
    progress: &mut F,
    stage: &str,
    current: u64,
    total: u64,
    message: &str,
) -> Result<(), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    progress(&ProgressEvent {
        stage: stage.to_string(),
        current,
        total,
        message: message.to_string(),
    })
}

pub fn parse_decimal(value: &str) -> Result<Integer, FactorError> {
    crate::integer::parse_decimal(value)
        .map_err(|error| FactorError::InvalidInput(format!("invalid decimal integer: {error}")))
}

pub fn probable_prime(value: &Integer) -> bool {
    crate::prime::is_probable_prime(value)
}

pub fn probable_prime_fast(value: &Integer, rounds: u64) -> bool {
    crate::prime::is_probable_prime_fast(value, rounds)
}

#[cfg(test)]
mod factor_tests {
    use super::better_pilot_pick;

    #[test]
    fn frobenius_gated_gate_preserves_default_verdicts() {
        let corpus = ["15", "91", "101", "10403", "100160063"];
        for digits in corpus {
            let value = super::parse_decimal(digits).expect("valid decimal");
            let default_report =
                super::factor(&value, &super::FactorOptions::default()).expect("factors");
            let gated_report = super::factor(
                &value,
                &super::FactorOptions {
                    use_frobenius: true,
                    ..super::FactorOptions::default()
                },
            )
            .expect("factors");
            assert_eq!(default_report.factors, gated_report.factors);
        }
        let prime = super::parse_decimal("101").expect("valid decimal");
        let report =
            super::factor(&prime, &super::FactorOptions::default()).expect("prime factors");
        assert_eq!(report.factors, vec!["101".to_string()]);
    }

    #[test]
    fn pilot_pick_prefers_measured_rate() {
        assert_eq!(better_pilot_pick(None, 2, 0.0, None, 60_000), None);
        assert_eq!(
            better_pilot_pick(None, 2, 1.5, None, 60_000),
            Some((2, 1.5, None))
        );
        assert_eq!(
            better_pilot_pick(Some((2, 1.5, None)), 3, 1.5, None, 60_000),
            Some((2, 1.5, None))
        );
        assert_eq!(
            better_pilot_pick(Some((2, 1.5, None)), 3, 2.5, None, 60_000),
            Some((3, 2.5, None))
        );
        assert_eq!(
            better_pilot_pick(Some((3, 2.5, None)), 4, 1.0, None, 60_000),
            Some((3, 2.5, None))
        );
    }

    #[test]
    fn prefilter_rejects_recorded_hw140_deg3_estimate() {
        let n = super::parse_decimal("825344960245717185263001325382831186370797")
            .expect("valid decimal");
        let recorded = super::FactorOptions {
            threads: Some(16),
            ..super::FactorOptions::default()
        };
        assert!(super::preflight_collect_exceeds(&n, 3, &recorded));
        let tiny = super::FactorOptions {
            max_memory_mb: Some(1),
            threads: Some(2),
            ..super::FactorOptions::default()
        };
        assert!(super::preflight_collect_exceeds(&n, 3, &tiny));
        assert!(super::preflight_collect_exceeds(&n, 4, &tiny));
        assert!(super::preflight_collect_exceeds(&n, 5, &tiny));
        let fitting = super::FactorOptions {
            threads: Some(8),
            ..super::FactorOptions::default()
        };
        let m = super::parse_decimal("770817968652704265841785071059").expect("valid decimal");
        assert!(!super::preflight_collect_exceeds(&m, 3, &fitting));
        assert!(!super::preflight_collect_exceeds(&n, 2, &tiny));
    }

    #[test]
    fn prefilter_feeds_escalation_identically() {
        let fed =
            super::note_pilot_for_escalation(None, None, 3, None, 0, true, None, 60_000, None);
        assert_eq!(fed, (None, Some(3)));
    }

    #[test]
    fn estimation_reject_skips_pilot_slice_without_collect() {
        let n = super::parse_decimal("825344960245717185263001325382831186370797")
            .expect("valid decimal");
        let options = super::FactorOptions {
            max_memory_mb: Some(1),
            max_duration: Some(std::time::Duration::from_millis(500)),
            threads: Some(2),
            skip_preprocessing: true,
            ..super::FactorOptions::default()
        };
        let collector = super::TimingCollector::default();
        let _ = super::factor_collected(&n, &options, |_| Ok(()), &collector);
        let rows = collector.snapshot();
        assert!(
            rows.iter().any(|entry| entry.stage == "race_pilot"
                && entry.degree == Some(3)
                && entry.detail.contains("err=ResourceLimit")),
            "degree 3 must still feed the escalation pick"
        );
        assert!(
            rows.iter().any(|entry| entry.stage == "race_pilot"
                && entry.degree == Some(3)
                && entry.detail.contains("slice_ms=0")),
            "degree 3 must be skipped before any pilot slice runs"
        );
        for entry in rows.iter().filter(|entry| entry.degree == Some(3)) {
            assert!(
                entry.stage == "race_pilot" || entry.stage == "cascade_projection",
                "degree 3 must do zero pre-slice work, got {}",
                entry.stage
            );
        }
    }

    #[test]
    fn budget_deadline_overflows_to_open() {
        assert!(super::budget_to_deadline(None).is_none());
        assert!(super::budget_to_deadline(Some(100)).is_some());
        assert!(super::budget_to_deadline(Some(u128::MAX)).is_none());
    }

    #[cfg(not(debug_assertions))]
    #[test]
    fn tight_budget_win_holds_in_release() {
        let n = super::parse_decimal(
            "166871515843890537614148105919789085907512212182703922830333764223351416812060867",
        )
        .expect("valid decimal");
        let options = super::FactorOptions {
            max_duration: Some(std::time::Duration::from_secs(1)),
            ..super::FactorOptions::default()
        };
        let (outcome, _) = super::factor_with_timings(&n, &options, |_| Ok(()));
        let report = outcome.expect("demo267 must factor inside a 1 s budget in release");
        assert_eq!(
            report.factors,
            vec![
                "121680281986369904489".to_string(),
                "1371393237423486234387340775961402466614689310566958802486603".to_string(),
            ]
        );
    }

    #[test]
    fn tight_budget_dispatch_runs_cheap_tier_with_recorded_slices() {
        let n = super::parse_decimal(
            "139607997073317332591812715527151593909624701340257947590834623479484508836409470911348221",
        )
        .expect("valid decimal");
        let options = super::FactorOptions {
            max_duration: Some(std::time::Duration::from_millis(2500)),
            ..super::FactorOptions::default()
        };
        let (outcome, timings) = super::factor_with_timings(&n, &options, |_| Ok(()));
        assert!(
            matches!(outcome, Err(crate::error::FactorError::ResourceLimit(_))),
            "90-digit input must stay out of reach on a 2.5 s slice, got {outcome:?}"
        );
        let cheap: Vec<_> = timings
            .iter()
            .filter(|entry| entry.stage == "ecm_cheap")
            .collect();
        assert!(
            !cheap.is_empty(),
            "cheap ECM tier must engage on the tight-budget path"
        );
        assert!(
            cheap.iter().any(|entry| entry.detail.contains("curves=")),
            "cheap tier must report curve progress"
        );
        let cascade = timings
            .iter()
            .find(|entry| entry.stage == "cascade")
            .expect("cascade record");
        assert!(
            cascade.detail.contains("slices="),
            "cascade record must carry deadline-derived slices: {}",
            cascade.detail
        );
    }

    #[test]
    fn cascade_projection_reports_after_ecm_miss() {
        let n = super::parse_decimal(
            "139607997073317332591812715527151593909624701340257947590834623479484508836409470911348221",
        )
        .expect("valid decimal");
        let options = super::FactorOptions {
            max_duration: Some(std::time::Duration::from_millis(2500)),
            ..super::FactorOptions::default()
        };
        let collector = super::TimingCollector::default();
        let outcome = super::factor_collected(&n, &options, |_| Ok(()), &collector);
        assert!(
            matches!(outcome, Err(crate::error::FactorError::ResourceLimit(_))),
            "tight-budget 90-digit run must end in ResourceLimit, got {outcome:?}"
        );
        let snapshot = collector.snapshot();
        let projections: Vec<&super::StageTiming> = snapshot
            .iter()
            .filter(|entry| entry.stage == "cascade_projection")
            .collect();
        assert!(
            !projections.is_empty(),
            "expected cascade_projection records after the race pilots"
        );
        for entry in &projections {
            assert!(
                entry.detail.contains("cost_ecm_ms=") && !entry.detail.contains("cost_ecm_ms=none"),
                "ECM cost must be measurable after a missed ECM stage: {}",
                entry.detail
            );
        }
    }

    #[test]
    fn failed_pilot_yield_feeds_escalation_pick() {
        let (best, limited) = super::note_pilot_for_escalation(
            None,
            None,
            3,
            Some(100),
            50,
            false,
            Some(400),
            60_000,
            None,
        );
        assert_eq!(best, Some((3, 2.0, Some(200))));
        assert_eq!(limited, None);
        let (best, limited) =
            super::note_pilot_for_escalation(None, None, 5, None, 50, true, None, 60_000, None);
        assert_eq!(best, None);
        assert_eq!(limited, Some(5));
        let (best, limited) = super::note_pilot_for_escalation(
            Some((2, 3.0, None)),
            Some(5),
            3,
            Some(100),
            50,
            false,
            None,
            60_000,
            None,
        );
        assert_eq!(best, Some((2, 3.0, None)));
        assert_eq!(limited, Some(5));
        assert_eq!(
            super::escalation_pick(Some((3, 2.0, Some(200))), None),
            Some((3, "best_rated", Some(2.0)))
        );
        assert_eq!(
            super::escalation_pick(None, Some(5)),
            Some((5, "first_limited", None))
        );
        assert_eq!(
            super::escalation_pick(Some((2, 3.0, None)), Some(5)),
            Some((2, "best_rated", Some(3.0)))
        );
        assert_eq!(super::escalation_pick(None, None), None);
        assert_eq!(super::optional_rate_text(None), "none");
        assert_eq!(super::optional_rate_text(Some(2.0)), "2.0000");
    }

    #[test]
    fn pilot_cost_link_skips_dominated_escalation() {
        let (best, limited) = super::note_pilot_for_escalation(
            None,
            None,
            3,
            Some(0),
            50,
            false,
            Some(500),
            60_000,
            None,
        );
        assert_eq!(best, None);
        assert_eq!(limited, None);
        let (best, limited) = super::note_pilot_for_escalation(
            None,
            None,
            3,
            Some(0),
            50,
            true,
            Some(500),
            60_000,
            None,
        );
        assert_eq!(best, None);
        assert_eq!(limited, None);
        assert_eq!(super::escalation_pick(best, limited), None);
        let (best, limited) = super::note_pilot_for_escalation(
            None,
            None,
            3,
            Some(100),
            50,
            false,
            Some(400),
            60_000,
            None,
        );
        assert_eq!(best, Some((3, 2.0, Some(200))));
        let (best, limited) = super::note_pilot_for_escalation(
            best,
            limited,
            2,
            Some(300),
            50,
            false,
            Some(300),
            60_000,
            None,
        );
        assert_eq!(best, Some((2, 6.0, Some(50))));
        assert_eq!(limited, None);
        assert_eq!(
            super::escalation_pick(best, limited),
            Some((2, "best_rated", Some(6.0)))
        );
        let (best, limited) = super::note_pilot_for_escalation(
            None,
            None,
            4,
            Some(100),
            50,
            false,
            Some(400),
            100,
            None,
        );
        assert_eq!(best, None);
        assert_eq!(limited, None);
        assert_eq!(super::escalation_pick(best, limited), None);
    }

    #[test]
    fn pilot_gnfs_cost_adds_measured_select() {
        let (best, _) = super::note_pilot_for_escalation(
            None,
            None,
            3,
            Some(100),
            50,
            false,
            Some(400),
            60_000,
            Some(30),
        );
        assert_eq!(best, Some((3, 2.0, Some(230))));
        let (best, _) = super::note_pilot_for_escalation(
            None,
            None,
            2,
            Some(300),
            50,
            false,
            Some(300),
            60_000,
            Some(30),
        );
        assert_eq!(best, Some((2, 6.0, Some(50))));
    }

    #[test]
    fn pilot_skip_needs_a_measured_unaffordable_projection() {
        assert_eq!(
            super::unaffordable_pilot_cost(None, Some(2.0), None, 60_000),
            None
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(400), None, None, 60_000),
            None
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(0), Some(2.0), None, 100),
            None
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(400), Some(0.0), None, 100),
            None
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(400), Some(2.0), None, 60_000),
            None
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(400), Some(2.0), None, 100),
            Some(200)
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(400), Some(2.0), Some(50), 200),
            Some(250)
        );
        assert_eq!(
            super::unaffordable_pilot_cost(Some(400), Some(2.0), Some(50), 250),
            None
        );
    }

    #[test]
    fn escalation_guard_fires_only_without_a_rated_pick() {
        assert!(super::escalation_has_no_rated_pick(None, Some(2)));
        assert!(!super::escalation_has_no_rated_pick(
            Some((2, 3.0, None)),
            Some(5)
        ));
        assert!(!super::escalation_has_no_rated_pick(
            Some((2, 3.0, Some(50))),
            None
        ));
        assert!(!super::escalation_has_no_rated_pick(None, None));
    }

    #[test]
    fn unknown_cost_never_displaces_known_affordable() {
        let known = Some((2, 6.0, Some(50)));
        let (best, _) = super::note_pilot_for_escalation(
            known,
            None,
            3,
            Some(300),
            50,
            false,
            None,
            60_000,
            None,
        );
        assert_eq!(best, known);
        let (incumbent, _) = super::note_pilot_for_escalation(
            None,
            None,
            3,
            Some(500),
            50,
            false,
            None,
            60_000,
            None,
        );
        assert_eq!(incumbent, Some((3, 10.0, None)));
        let (best, _) = super::note_pilot_for_escalation(
            incumbent,
            None,
            2,
            Some(100),
            50,
            false,
            Some(400),
            60_000,
            None,
        );
        assert_eq!(best, Some((2, 2.0, Some(200))));
        assert_eq!(
            super::escalation_pick(best, None),
            Some((2, "best_rated", Some(2.0)))
        );
        let (best, limited) = super::note_pilot_for_escalation(
            None,
            None,
            2,
            Some(100),
            50,
            true,
            Some(400),
            60_000,
            None,
        );
        let (best, limited) = super::note_pilot_for_escalation(
            best,
            limited,
            3,
            Some(100),
            50,
            true,
            None,
            60_000,
            None,
        );
        assert_eq!(best, Some((2, 2.0, Some(200))));
        assert_eq!(
            super::escalation_pick(best, limited),
            Some((2, "best_rated", Some(2.0)))
        );
        assert_eq!(
            super::escalation_pick(None, limited),
            Some((2, "first_limited", None))
        );
    }

    #[test]
    fn retained_timings_survive_failure() {
        let n = super::parse_decimal("2255755539706396516527333142242253921536071")
            .expect("valid decimal");
        let options = super::FactorOptions {
            max_duration: Some(std::time::Duration::from_millis(50)),
            ..super::FactorOptions::default()
        };
        let (outcome, timings) = super::factor_with_timings(&n, &options, |_| Ok(()));
        assert!(
            outcome.is_err(),
            "a 50 ms budget on a 141-bit input must fail, got {outcome:?}"
        );
        assert!(
            !timings.is_empty(),
            "stage timings collected before the error must be retained on the error path"
        );
    }

    #[test]
    fn ecm_miss_bound_spills_sieves_to_collect() {
        let n = super::parse_decimal(
            "139607997073317332591812715527151593909624701340257947590834623479484508836409470911348221",
        )
        .expect("valid decimal");
        let options = super::FactorOptions {
            max_duration: Some(std::time::Duration::from_millis(2500)),
            ..super::FactorOptions::default()
        };
        let (outcome, timings) = super::factor_with_timings(&n, &options, |_| Ok(()));
        assert!(
            matches!(outcome, Err(crate::error::FactorError::ResourceLimit(_))),
            "tight-budget 90-digit run must end in ResourceLimit, got {outcome:?}"
        );
        let cascade = timings
            .iter()
            .find(|entry| entry.stage == "cascade")
            .expect("cascade record");
        let cap_text = cascade
            .detail
            .split("ecm_cap=")
            .nth(1)
            .and_then(|tail| tail.split_whitespace().next())
            .expect("cascade record must carry the ecm burn bound");
        if cap_text != "open" {
            let cap_ms: u128 = cap_text.parse().expect("burn bound must be numeric");
            assert!(
                cap_ms > 0,
                "burn bound must leave a live ecm box, got {cap_ms}"
            );
        }
        let ecm_miss = timings
            .iter()
            .any(|entry| entry.stage == "ecm" && entry.detail.contains("no divisor"));
        assert!(ecm_miss, "miss input must reach main ECM and miss");
        let sieve_reached = timings
            .iter()
            .any(|entry| entry.stage == "cascade_projection");
        assert!(
            sieve_reached,
            "sieve pilots must still run under the burn bound"
        );
    }

    #[test]
    fn ecm_burn_bound_reports_numeric_cap() {
        let n = super::parse_decimal(
            "139607997073317332591812715527151593909624701340257947590834623479484508836409470911348221",
        )
        .expect("valid decimal");
        let mut sightings: Vec<String> = Vec::new();
        for _ in 0..3 {
            let options = super::FactorOptions {
                max_duration: Some(std::time::Duration::from_millis(1800)),
                ..super::FactorOptions::default()
            };
            let (outcome, timings) = super::factor_with_timings(&n, &options, |_| Ok(()));
            assert!(
                matches!(outcome, Err(crate::error::FactorError::ResourceLimit(_))),
                "tight-budget 90-digit run must end in ResourceLimit, got {outcome:?}"
            );
            let ecm_miss = timings
                .iter()
                .any(|entry| entry.stage == "ecm" && entry.detail.contains("no divisor"));
            assert!(ecm_miss, "miss input must reach main ECM and miss");
            let cap_text = timings
                .iter()
                .find(|entry| entry.stage == "cascade")
                .and_then(|entry| {
                    entry
                        .detail
                        .split("ecm_cap=")
                        .nth(1)
                        .and_then(|tail| tail.split_whitespace().next())
                })
                .expect("cascade record must carry the ecm burn bound");
            if cap_text != "open" && cap_text.parse::<u128>().is_ok_and(|cap_ms| cap_ms > 0) {
                return;
            }
            sightings.push(cap_text.to_string());
        }
        panic!("burn bound never produced a numeric cap sighting in 3 attempts, got {sightings:?}");
    }

    #[test]
    fn joint_ecm_win_records_ecm_with_joint_and_reclaim_rows() {
        let n = super::parse_decimal("3994424150216394531668178182146365185669")
            .expect("valid decimal");
        let options = super::FactorOptions::default();
        let (outcome, snapshot) = super::factor_with_timings(&n, &options, |_| Ok(()));
        let report = outcome.expect("joint ECM input factors");
        assert_eq!(report.method, "ecm");
        assert!(!report.method.contains("joint"));
        assert!(!report.method.contains("reclaim"));
        let joint = report
            .timings
            .iter()
            .find(|entry| entry.stage == "joint")
            .expect("joint decision row present");
        assert!(
            joint.detail.contains("winner=ecm"),
            "joint row must name the ECM winner: {}",
            joint.detail
        );
        assert!(
            report.timings.iter().any(|entry| entry.stage == "reclaim"),
            "reclaim decision row present"
        );
        assert!(
            report
                .timings
                .iter()
                .any(|entry| entry.stage == "ecm" && entry.detail.contains("divisor=")),
            "winning ECM tranche row present"
        );
        assert_eq!(report.timings.len(), snapshot.len());
        for (report_entry, snapshot_entry) in report.timings.iter().zip(snapshot.iter()) {
            assert_eq!(report_entry.stage, snapshot_entry.stage);
            assert_eq!(report_entry.degree, snapshot_entry.degree);
            assert_eq!(report_entry.elapsed_ms, snapshot_entry.elapsed_ms);
            assert_eq!(report_entry.detail, snapshot_entry.detail);
        }
        let mut product = super::Integer::from(1);
        for factor in &report.factors {
            product *= &super::parse_decimal(factor).expect("valid factor");
        }
        assert_eq!(product, n);
    }

    #[test]
    fn siqs_win_records_siqs_with_siqs_rows() {
        let n = super::parse_decimal("100460333").expect("valid decimal");
        let options = super::FactorOptions {
            degree: Some(2),
            skip_preprocessing: true,
            threads: Some(2),
            max_duration: Some(std::time::Duration::from_secs(300)),
            allow_fallback: false,
            ..super::FactorOptions::default()
        };
        let (outcome, snapshot) = super::factor_with_timings(&n, &options, |_| Ok(()));
        let report = outcome.expect("SIQS factorization succeeds");
        assert_eq!(report.method, "siqs");
        assert!(!report.method.contains("joint"));
        assert!(!report.method.contains("reclaim"));
        assert!(
            report
                .timings
                .iter()
                .any(|entry| entry.stage.starts_with("siqs")),
            "SIQS stage rows present"
        );
        assert_eq!(report.timings.len(), snapshot.len());
    }

    #[test]
    fn method_vocab_never_leaks_scheduler_names() {
        let cases: Vec<(String, super::FactorOptions, String)> = vec![
            (
                "3000099".to_string(),
                super::FactorOptions::default(),
                "trial_division".to_string(),
            ),
            (
                "97".to_string(),
                super::FactorOptions::default(),
                "prime".to_string(),
            ),
            (
                "1000036000099".to_string(),
                super::FactorOptions::default(),
                "fermat".to_string(),
            ),
            (
                "17497577401".to_string(),
                super::FactorOptions {
                    use_fermat: false,
                    use_pp1: false,
                    use_ecm: false,
                    use_gnfs: false,
                    allow_fallback: false,
                    ..super::FactorOptions::default()
                },
                "pm1".to_string(),
            ),
            (
                "100460333".to_string(),
                super::FactorOptions {
                    degree: Some(2),
                    skip_preprocessing: true,
                    threads: Some(2),
                    max_duration: Some(std::time::Duration::from_secs(300)),
                    allow_fallback: false,
                    ..super::FactorOptions::default()
                },
                "siqs".to_string(),
            ),
        ];
        let allowed = [
            "trial_division",
            "fermat",
            "pm1",
            "pp1",
            "ecm",
            "prime",
            "perfect_square",
            "perfect_power",
            "siqs",
            "quadratic-qs",
            "gnfs",
            "fallback",
        ];
        let scheduler = [
            "joint",
            "reclaim",
            "governor",
            "cascade_projection",
            "pilot_rate",
            "race_pilot",
            "race_escalation",
        ];
        for (decimal, options, expected) in &cases {
            let n = super::parse_decimal(decimal).expect("valid decimal");
            let (outcome, _) = super::factor_with_timings(&n, options, |_| Ok(()));
            let report = outcome.expect("factorization succeeds");
            assert_eq!(report.method, *expected, "input {decimal}");
            for part in report.method.split('+') {
                assert!(
                    allowed.contains(&part),
                    "method part {part} outside pre-existing vocabulary on input {decimal}"
                );
            }
            for name in &scheduler {
                assert!(
                    !report.method.contains(name),
                    "scheduler stage {name} leaked into method on input {decimal}"
                );
            }
        }
        let n = super::parse_decimal("3994424150216394531668178182146365185669")
            .expect("valid decimal");
        let (outcome, _) =
            super::factor_with_timings(&n, &super::FactorOptions::default(), |_| Ok(()));
        let report = outcome.expect("joint ECM input factors");
        assert_eq!(report.method, "ecm");
        for name in &scheduler {
            assert!(
                !report.method.contains(name),
                "scheduler stage {name} leaked into method on joint input"
            );
        }
    }

    fn trial_oracle(n: &super::Integer) -> (Option<super::Integer>, u64) {
        let mut checked = 0_u64;
        let small_n = super::IntegerExt::to_u64_saturating(n);
        let big_hit = |prime: u64| n % &super::Integer::from(prime) == super::Integer::from(0);
        for prime in crate::prime::small_primes() {
            checked += 1;
            if big_hit(*prime) {
                return (Some(super::Integer::from(*prime)), checked);
            }
            if prime
                .checked_mul(*prime)
                .is_some_and(|square| square > small_n)
            {
                break;
            }
        }
        (None, checked)
    }

    #[test]
    fn trial_division_matches_bigint_oracle_on_recorded_vectors() {
        let decimals = [
            "2",
            "3",
            "4",
            "9",
            "25",
            "49",
            "121",
            "240",
            "145",
            "926",
            "1177",
            "9999",
            "10000",
            "10007",
            "10403",
            "999983",
            "1000003",
            "18446744073709551557",
            "18446744073709551615",
            "18446744073709551614",
            "340282366920938463463374607431768211455",
            "770817968652704265841785071059",
            "3426790303412814276981061",
            "1155285037654380637",
        ];
        for decimal in decimals {
            let n = super::parse_decimal(decimal).expect("valid decimal");
            assert_eq!(
                super::trial_division(&n),
                trial_oracle(&n),
                "trial split at {decimal}"
            );
        }
        for small in [0_u64, 1, 2, 9_999, 10_000, 10_001] {
            let n = super::Integer::from(small);
            assert_eq!(
                super::trial_division(&n),
                trial_oracle(&n),
                "trial split at u64 {small}"
            );
        }
        let wide = super::Integer::from(u64::MAX);
        assert_eq!(super::trial_division(&wide), trial_oracle(&wide));
    }
}

#[cfg(test)]
mod snfs_tests {
    use super::{detect_special_form, exact_small_coefficient_form};

    fn expect_exact(decimal: &str, a: u64, b: u64, k: u32, c: u64, sign: i8, low_weight: bool) {
        let n = super::parse_decimal(decimal).expect("valid input");
        let form = detect_special_form(&n).expect("special form detected");
        assert_eq!(form.a, a, "a mismatch for {decimal}");
        assert_eq!(form.b, b, "b mismatch for {decimal}");
        assert_eq!(form.k, k, "k mismatch for {decimal}");
        assert_eq!(form.c, c, "c mismatch for {decimal}");
        assert_eq!(form.sign, sign, "sign mismatch for {decimal}");
        assert_eq!(
            form.low_weight, low_weight,
            "low_weight mismatch for {decimal}"
        );
    }

    #[test]
    fn detects_small_corpus_with_exact_parameters() {
        expect_exact("2047", 1, 2, 11, 1, -1, true);
        expect_exact("32767", 1, 2, 15, 1, -1, false);
        expect_exact("2097151", 1, 2, 21, 1, -1, false);
        expect_exact("730", 1, 3, 6, 1, 1, true);
    }

    #[test]
    fn detects_large_corpus_with_exact_parameters() {
        expect_exact(
            "51422017416287688817342786954917203280710495801049370729644031",
            1,
            2,
            205,
            1,
            -1,
            false,
        );
        expect_exact(
            "65185151242703554760590262029100101153646988597309960020356494379340201592426774597868716033",
            1,
            2,
            305,
            1,
            1,
            true,
        );
        expect_exact(
            "18511095575145533338447403836109581505353226740782032803932380476024584374357028238763044",
            1,
            3,
            185,
            1,
            1,
            false,
        );
    }

    #[test]
    fn prefers_primitive_exponent_over_composite_power() {
        let form =
            detect_special_form(&super::parse_decimal("32767").expect("valid")).expect("detected");
        assert_eq!((form.b, form.k), (2, 15));
        let form =
            detect_special_form(&super::parse_decimal("730").expect("valid")).expect("detected");
        assert_eq!((form.b, form.k), (3, 6));
        let form =
            detect_special_form(&super::parse_decimal("9").expect("valid")).expect("detected");
        assert_eq!((form.a, form.b, form.k, form.c, form.sign), (1, 2, 3, 1, 1));
    }

    #[test]
    fn rejects_canonical_inputs() {
        for decimal in [
            "770817968652704265841785071059",
            "3426790303412814276981061",
            "1155285037654380637",
            "166871515843890537614148105919789085907512212182703922830333764223351416812060867",
        ] {
            let n = super::parse_decimal(decimal).expect("valid input");
            assert!(
                detect_special_form(&n).is_none(),
                "canonical input must not be detected: {decimal}"
            );
        }
    }

    #[test]
    fn base_limit_boundary() {
        expect_exact("999966000290", 1, 999983, 2, 1, 1, false);
        let n = super::parse_decimal("1000006000008").expect("valid input");
        assert!(
            detect_special_form(&n).is_none(),
            "base above the limit must not be detected"
        );
    }

    #[test]
    fn low_weight_constant_is_twelve() {
        assert_eq!(super::LOW_HAMMING_WEIGHT_LIMIT, 12);
        let heleven = super::parse_decimal("1000000192513").expect("valid input");
        let form = detect_special_form(&heleven).expect("heuristic detection");
        assert_eq!(form.k, 0);
        assert!(form.low_weight);
        let thirteen = super::parse_decimal("1000000061441").expect("valid input");
        assert!(
            detect_special_form(&thirteen).is_none(),
            "popcount thirteen must not be detected"
        );
    }

    fn load_column_opt(relative: &str, column: usize) -> Option<Vec<String>> {
        let path = format!("{}/../../{}", env!("CARGO_MANIFEST_DIR"), relative);
        let text = std::fs::read_to_string(&path).ok()?;
        Some(
            text.lines()
                .filter(|line| !line.trim().is_empty())
                .map(|line| {
                    line.split_whitespace()
                        .nth(column)
                        .expect("column present")
                        .to_string()
                })
                .collect(),
        )
    }

    #[test]
    fn sweep_general_corpora_has_zero_exact_hits() {
        let files = [
            ("far16", "logs/sq-cases2.txt", 1),
            ("rs35", "logs/rs-cases.txt", 0),
            ("r9", "probe/r9-cases.txt", 0),
            ("cases-100120", "probe/cases-100120.txt", 0),
            ("cases-105", "probe/cases-105.txt", 0),
            ("cases-110", "probe/cases-110.txt", 0),
            ("cases-8486", "probe/cases-8486.txt", 0),
            ("cases-8788", "probe/cases-8788.txt", 0),
            ("cases-8889", "probe/cases-8889.txt", 0),
            ("cases-8990", "probe/cases-8990.txt", 0),
            ("m2", "probe/m2-cases.txt", 0),
        ];
        let mut checked = 0_usize;
        let mut skipped = 0_usize;
        for (name, relative, column) in files {
            let Some(entries) = load_column_opt(relative, column) else {
                skipped += 1;
                println!("SNFS_SWEEP_SKIP name={name} relative={relative}");
                continue;
            };
            checked += 1;
            let mut heuristic = 0_usize;
            for decimal in &entries {
                let n = super::parse_decimal(decimal).expect("valid corpus input");
                match detect_special_form(&n) {
                    None => {}
                    Some(form) if form.k == 0 => heuristic += 1,
                    Some(form) => panic!(
                        "exact false positive on {name} input {decimal}: {}",
                        form.describe()
                    ),
                }
            }
            println!(
                "SNFS_SWEEP name={name} count={} exact=0 heuristic={heuristic}",
                entries.len()
            );
        }
        println!("SNFS_SWEEP_SUMMARY checked={checked} skipped={skipped}");
    }

    #[test]
    fn detector_cost_stays_within_budget() {
        let mid = super::parse_decimal("3426790303412814276981061").expect("valid");
        let big = super::parse_decimal(
            "166871515843890537614148105919789085907512212182703922830333764223351416812060867",
        )
        .expect("valid");
        let start = std::time::Instant::now();
        let _ = detect_special_form(&mid);
        let mid_ms = start.elapsed().as_millis();
        let start = std::time::Instant::now();
        let _ = detect_special_form(&big);
        let big_ms = start.elapsed().as_millis();
        let start = std::time::Instant::now();
        let _ = exact_small_coefficient_form(&mid);
        let inner_ms = start.elapsed().as_millis();
        println!("SNFS_COST mid_ms={mid_ms} big_ms={big_ms} ungated_small_a_on_mid_ms={inner_ms}");
        assert!(mid_ms < 1000, "82-bit detector cost exploded: {mid_ms} ms");
        assert!(big_ms < 5000, "demo267 detector cost exploded: {big_ms} ms");
    }

    #[test]
    fn dumps_gnfs_polynomials_for_projection() {
        let inputs = [
            (
                "L1",
                "51422017416287688817342786954917203280710495801049370729644031",
            ),
            (
                "L2",
                "65185151242703554760590262029100101153646988597309960020356494379340201592426774597868716033",
            ),
            (
                "L3",
                "18511095575145533338447403836109581505353226740782032803932380476024584374357028238763044",
            ),
        ];
        for (tag, decimal) in inputs {
            let n = super::parse_decimal(decimal).expect("valid input");
            for degree in [2_usize, 3, 4, 5] {
                let pair = if degree == 2 {
                    let attempt = super::FactorOptions {
                        degree: Some(2),
                        ..super::FactorOptions::default()
                    };
                    crate::polynomial::select_polynomial(&n, &attempt)
                        .expect("quadratic polynomial")
                } else {
                    crate::general::inspect_polynomial(&n, degree).expect("general polynomial")
                };
                assert_eq!(pair.degree, degree);
                assert_eq!(pair.evaluate(&pair.m), n);
                let coeffs: Vec<String> =
                    pair.coefficients.iter().map(ToString::to_string).collect();
                println!(
                    "SNFS_POLY tag={tag} degree={degree} m={} coeffs={}",
                    pair.m,
                    coeffs.join(",")
                );
            }
        }
    }
}

#[cfg(test)]
mod cpu_budget_tests {
    use super::{TimingCollector, parse_decimal};
    use crate::options::FactorOptions;

    #[test]
    fn accumulator_charges_wall_times_threads() {
        let collector = TimingCollector::default();
        assert_eq!(collector.cpu_used_ms(), 0);
        collector.add_cpu_ms(500, 16);
        assert_eq!(collector.cpu_used_ms(), 8_000);
        collector.add_cpu_ms(250, 4);
        assert_eq!(collector.cpu_used_ms(), 9_000);
    }

    #[test]
    fn accumulator_saturates_instead_of_wrapping() {
        let collector = TimingCollector::default();
        collector.add_cpu_ms(u128::MAX, 16);
        assert_eq!(collector.cpu_used_ms(), u128::MAX);
        collector.add_cpu_ms(1, 1);
        assert_eq!(collector.cpu_used_ms(), u128::MAX);
    }

    #[test]
    fn unset_cap_never_trips() {
        let collector = TimingCollector::default();
        collector.add_cpu_ms(60_000, 16);
        assert!(
            collector
                .check_cpu_budget(&FactorOptions::default())
                .is_ok()
        );
    }

    #[test]
    fn trip_reports_cpu_against_budget() {
        let collector = TimingCollector::default();
        let options = FactorOptions {
            max_cpu_seconds: Some(1),
            threads: Some(16),
            ..FactorOptions::default()
        };
        assert!(collector.check_cpu_budget(&options).is_ok());
        collector.add_cpu_ms(63, 16);
        let error = collector
            .check_cpu_budget(&options)
            .expect_err("63 ms across 16 threads must exceed a 1 s cpu budget");
        match error {
            crate::FactorError::ResourceLimit(detail) => {
                assert!(detail.contains("cpu time budget exceeded"), "{detail}");
                assert!(detail.contains("1008"), "{detail}");
                assert!(detail.contains("1000"), "{detail}");
                assert!(detail.contains("cap 1 s"), "{detail}");
            }
            other => panic!("expected ResourceLimit, got {other:?}"),
        }
    }

    #[test]
    fn collectors_stay_independent_across_threads() {
        let handles: Vec<_> = (0..8_u128)
            .map(|index| {
                std::thread::spawn(move || {
                    let collector = TimingCollector::default();
                    collector.add_cpu_ms(100 + index, 4);
                    collector.cpu_used_ms()
                })
            })
            .collect();
        let mut totals: Vec<u128> = handles
            .into_iter()
            .map(|handle| handle.join().expect("worker joins"))
            .collect();
        totals.sort_unstable();
        let expected: Vec<u128> = (0..8_u128).map(|index| (100 + index) * 4).collect();
        assert_eq!(totals, expected);
    }

    #[test]
    fn concurrent_factor_calls_keep_independent_budgets() {
        let plain = std::thread::spawn(|| {
            let n = parse_decimal("15").expect("valid decimal");
            super::factor(&n, &FactorOptions::default()).map(|report| report.factors)
        });
        let capped = std::thread::spawn(|| {
            let n = parse_decimal("15").expect("valid decimal");
            let options = FactorOptions {
                max_cpu_seconds: Some(0),
                ..FactorOptions::default()
            };
            super::factor(&n, &options).map(|report| report.factors)
        });
        let plain_factors = plain
            .join()
            .expect("plain call joins")
            .expect("plain factors");
        assert_eq!(plain_factors, vec!["3".to_string(), "5".to_string()]);
        let capped_outcome = capped.join().expect("capped call joins");
        match &capped_outcome {
            Err(crate::FactorError::ResourceLimit(detail))
                if detail.contains("cpu time budget") => {}
            other => panic!("zero cpu budget must trip immediately, got {other:?}"),
        }
    }
}
