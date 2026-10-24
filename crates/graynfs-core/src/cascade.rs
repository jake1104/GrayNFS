use crate::integer::Integer;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CascadeStage {
    Fermat,
    Squfof,
    Pm1,
    Pp1,
    Ecm,
    Mpqs,
    Gnfs,
    Rho,
}

impl CascadeStage {
    pub fn name(self) -> &'static str {
        match self {
            CascadeStage::Fermat => "fermat",
            CascadeStage::Squfof => "squfof",
            CascadeStage::Pm1 => "pm1",
            CascadeStage::Pp1 => "pp1",
            CascadeStage::Ecm => "ecm",
            CascadeStage::Mpqs => "mpqs",
            CascadeStage::Gnfs => "gnfs",
            CascadeStage::Rho => "rho",
        }
    }
}

#[derive(Clone, Debug)]
pub struct CascadeSignals {
    pub bits: u64,
    pub target_digits: u32,
    pub failed_curves: u64,
    pub measured_ms_per_curve: Option<u128>,
    pub pilot_relations_per_ms: Option<f64>,
    pub trial_division_primes_checked: Option<u64>,
    pub workers: usize,
    pub pilot_relations_per_sec: Option<f64>,
    pub remaining_ms_at_collect: Option<u128>,
}

impl CascadeSignals {
    pub fn collect(n: &Integer, deadline: Option<Instant>) -> Self {
        let bits = n.bits();
        let digits = (bits * 301 / 1000).min(u64::from(u32::MAX)).max(1) as u32;
        let workers = std::thread::available_parallelism()
            .map(|value| value.get())
            .unwrap_or(1);
        let remaining_ms_at_collect =
            deadline.map(|limit| limit.saturating_duration_since(Instant::now()).as_millis());
        Self {
            bits,
            target_digits: digits,
            failed_curves: 0,
            measured_ms_per_curve: None,
            pilot_relations_per_ms: None,
            trial_division_primes_checked: None,
            workers,
            pilot_relations_per_sec: None,
            remaining_ms_at_collect,
        }
    }

    pub fn observe_workers(&mut self, count: usize) {
        self.workers = count.max(1);
    }

    pub fn observe_options(
        &mut self,
        options: &crate::options::FactorOptions,
        deadline: Option<Instant>,
    ) {
        self.workers = options.effective_workers();
        if let Some(limit) = deadline {
            self.remaining_ms_at_collect =
                Some(limit.saturating_duration_since(Instant::now()).as_millis());
        } else if let Some(span) = options.max_duration {
            self.remaining_ms_at_collect = Some(span.as_millis());
        }
    }

    pub fn collect_with_options(
        n: &Integer,
        options: &crate::options::FactorOptions,
        deadline: Option<Instant>,
    ) -> Self {
        let resolved = deadline.or_else(|| options.effective_deadline(Instant::now()));
        let mut signals = Self::collect(n, resolved);
        signals.observe_options(options, resolved);
        signals
    }

    pub fn observe_trial_division(&mut self, checked: u64) {
        self.trial_division_primes_checked = Some(checked);
    }

    pub fn report_pilot(&mut self, relations: usize, elapsed_ms: u128) {
        let rate = relations as f64 / elapsed_ms.max(1) as f64;
        if rate.is_finite() && rate > 0.0 {
            self.pilot_relations_per_ms = Some(rate);
            self.pilot_relations_per_sec = Some(rate * 1000.0);
        }
    }

    pub fn measured_relations_per_sec(&self) -> Option<f64> {
        if let Some(rate) = self.pilot_relations_per_sec {
            return Some(rate);
        }
        self.pilot_relations_per_ms
            .filter(|rate| rate.is_finite() && *rate > 0.0)
            .map(|rate| rate * 1000.0)
    }

    pub fn probe_budget_ms(n: &Integer) -> u128 {
        crate::options::probe_budget_ms_for_bits(n.bits())
    }

    pub fn effective_remaining_ms(&self, n: &Integer, deadline: Option<Instant>) -> u128 {
        self.remaining_ms(deadline)
            .unwrap_or_else(|| Self::probe_budget_ms(n))
    }

    pub fn report_ecm(&mut self, curves_run: u64, elapsed_ms: u128, found: bool) {
        if curves_run > 0 {
            self.measured_ms_per_curve = Some(elapsed_ms / u128::from(curves_run));
        }
        if !found {
            self.failed_curves = self.failed_curves.saturating_add(curves_run);
        }
    }

    pub fn remaining_ms(&self, deadline: Option<Instant>) -> Option<u128> {
        deadline.map(|limit| limit.saturating_duration_since(Instant::now()).as_millis())
    }
}

pub fn ecm_plan(target_digits: u32) -> (u64, u64, u64) {
    if target_digits <= 20 {
        return (11_000, 1_900_000, 74);
    }
    if target_digits <= 25 {
        return (50_000, 13_000_000, 214);
    }
    let mut steps = u64::from(target_digits).saturating_sub(25).div_ceil(5);
    let mut bound_one = 50_000_u64;
    let mut bound_two = 13_000_000_u64;
    let mut curves = 214_u64;
    while steps > 0 {
        bound_one = bound_one.saturating_mul(5);
        bound_two = bound_two.saturating_mul(7);
        curves = curves.saturating_mul(3);
        steps -= 1;
    }
    (bound_one, bound_two, curves)
}

pub fn project_ecm_cost_ms(signals: &CascadeSignals) -> Option<u128> {
    project_ecm_cost_ms_for_digits(signals, signals.target_digits)
}

pub fn project_ecm_cost_ms_for_digits(signals: &CascadeSignals, digits: u32) -> Option<u128> {
    let per_curve = signals.measured_ms_per_curve?;
    let (_, _, curves) = ecm_plan(digits);
    u128::from(curves).checked_mul(per_curve.max(1))
}

pub fn project_siqs_cost_ms(relations_needed: usize, relations_per_ms: f64) -> Option<u128> {
    if relations_needed == 0 || !relations_per_ms.is_finite() || relations_per_ms <= 0.0 {
        return None;
    }
    let cost_ms = relations_needed as f64 / relations_per_ms;
    if cost_ms.is_finite() {
        Some(cost_ms.ceil() as u128)
    } else {
        None
    }
}

pub fn projected_lower_cost(
    cost_ecm_ms: Option<u128>,
    cost_siqs_ms: Option<u128>,
) -> Option<&'static str> {
    match (cost_ecm_ms, cost_siqs_ms) {
        (Some(ecm), Some(siqs)) if ecm > siqs => Some("siqs"),
        (Some(ecm), Some(siqs)) if ecm < siqs => Some("ecm"),
        (Some(_), Some(_)) => Some("tie"),
        _ => None,
    }
}

pub fn project_ecm_remaining_ms(signals: &CascadeSignals, digits: u32) -> Option<u128> {
    let per_curve = match signals.measured_ms_per_curve {
        Some(cost) if cost > 0 => cost,
        _ => return None,
    };
    let mut tier = digits.min(1000);
    for _ in 0..64 {
        let (_, _, curves) = ecm_plan(tier);
        if u128::from(curves) > u128::from(signals.failed_curves) {
            return u128::from(curves)
                .checked_sub(u128::from(signals.failed_curves))
                .and_then(|remaining| remaining.checked_mul(per_curve));
        }
        tier = tier.saturating_add(5);
    }
    None
}

pub fn project_siqs_remaining_ms(
    signals: &CascadeSignals,
    relations_needed: usize,
) -> Option<u128> {
    if relations_needed == 0 {
        return None;
    }
    let rate = signals.measured_relations_per_sec()?;
    if !rate.is_finite() || rate <= 0.0 {
        return None;
    }
    let cost_ms = relations_needed as f64 / rate * 1000.0;
    if cost_ms.is_finite() {
        Some(cost_ms.ceil() as u128)
    } else {
        None
    }
}

pub fn cost_shares(cost_ecm_ms: Option<u128>, cost_siqs_ms: Option<u128>) -> Option<(f64, f64)> {
    let cost_ecm = cost_ecm_ms? as f64;
    let cost_siqs = cost_siqs_ms? as f64;
    let total = cost_ecm + cost_siqs;
    if total <= 0.0 {
        return Some((0.5, 0.5));
    }
    let ecm_share = cost_siqs / total;
    let siqs_share = cost_ecm / total;
    Some((ecm_share, siqs_share))
}

pub fn should_switch_to_siqs(signals: &CascadeSignals, relations_needed: usize) -> Option<bool> {
    let ecm_next = project_ecm_remaining_ms(signals, escalate_digits(signals.target_digits))?;
    let siqs = project_siqs_remaining_ms(signals, relations_needed)?;
    Some(ecm_next > siqs)
}

pub fn governor_slices(
    signals: &CascadeSignals,
    digits: u32,
    relations_needed: usize,
) -> Option<(f64, f64)> {
    let ecm = project_ecm_remaining_ms(signals, digits);
    let siqs = project_siqs_remaining_ms(signals, relations_needed);
    cost_shares(ecm, siqs)
}

pub fn allocate_slices(
    signals: &CascadeSignals,
    remaining_ms: Option<u128>,
) -> (Option<u128>, Option<u128>, Option<u128>, Option<u128>) {
    let remaining = match remaining_ms {
        Some(value) => value,
        None => return (None, None, None, None),
    };
    let fermat_divisor = if signals.target_digits <= 25 { 32 } else { 64 };
    let fermat = remaining / fermat_divisor;
    let pm1 = remaining / 16;
    let pp1 = remaining / 16;
    let ecm = remaining * 3 / 4;
    if let Some(projected) = project_ecm_cost_ms(signals) {
        let shrunk = fermat.min((projected / 8).max(2000));
        let savings = fermat.saturating_sub(shrunk);
        return (
            Some(shrunk),
            Some(pm1),
            Some(pp1),
            Some(ecm.saturating_add(savings)),
        );
    }
    (Some(fermat), Some(pm1), Some(pp1), Some(ecm))
}

pub fn cap_ecm_slice(
    signals: &CascadeSignals,
    live_remaining_ms: Option<u128>,
    planned_ecm_ms: Option<u128>,
) -> Option<u128> {
    const NEGLIGIBLE_BURN_DIVISOR: u128 = 32;
    let live = live_remaining_ms?;
    let planned = planned_ecm_ms?;
    if signals.failed_curves == 0 {
        return Some(planned);
    }
    let per_curve = match signals.measured_ms_per_curve {
        Some(cost) if cost > 0 => cost,
        _ => return Some(planned),
    };
    let burn = per_curve.saturating_mul(u128::from(signals.failed_curves));
    if burn.saturating_mul(NEGLIGIBLE_BURN_DIVISOR) < live {
        return Some(planned);
    }
    Some(planned.min(live / 2))
}

pub(crate) fn tranche_first_box(
    ecm_budget: Option<u128>,
    live_remaining_ms: Option<u128>,
) -> Option<u128> {
    let box_ms = ecm_budget?;
    let half = box_ms / 2;
    match live_remaining_ms {
        Some(live) => Some(half.min(live / 2)),
        None => Some(half),
    }
}

pub fn cheap_slice_ms(remaining_ms: Option<u128>) -> Option<u128> {
    remaining_ms.map(|remaining| remaining / 16)
}

pub fn cheap_ecm_plan(remaining_ms: Option<u128>) -> (u64, u64, u64) {
    const BOUND_ONE: u64 = 2_000;
    const BOUND_TWO: u64 = 200_000;
    const MIN_CURVES: u64 = 4;
    const MAX_CURVES: u64 = 32;
    const MS_PER_CURVE: u128 = 64;
    let curves = match remaining_ms {
        Some(remaining) => {
            (remaining / MS_PER_CURVE).clamp(MIN_CURVES as u128, MAX_CURVES as u128) as u64
        }
        None => MAX_CURVES,
    };
    (BOUND_ONE, BOUND_TWO, curves)
}

pub fn fermat_cap_for(n: &Integer) -> u64 {
    const CEILING: u64 = 200_000;
    const FLOOR: u64 = 20_000;
    const REACH: f64 = 65_536.0;
    const SCALE: f64 = 85.0;
    let bits = n.bits();
    if bits == 0 {
        return CEILING;
    }
    let sqrt_n = 2_f64.powf(bits as f64 / 2.0);
    let estimated = if sqrt_n.is_finite() && sqrt_n > 0.0 {
        REACH * REACH / (8.0 * sqrt_n)
    } else {
        0.0
    };
    let scaled = (SCALE * estimated).clamp(0.0, CEILING as f64);
    (FLOOR as f64).max(scaled).min(CEILING as f64) as u64
}

pub fn escalate_digits(target_digits: u32) -> u32 {
    target_digits.saturating_add(5)
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ReclaimStep {
    None,
    SiqsWidth,
    EcmCap,
    Threads,
}

impl ReclaimStep {
    pub fn name(self) -> &'static str {
        match self {
            ReclaimStep::None => "none",
            ReclaimStep::SiqsWidth => "siqs_width",
            ReclaimStep::EcmCap => "ecm_cap",
            ReclaimStep::Threads => "threads",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct ReclaimState {
    pub width_done: bool,
    pub cap_done: bool,
    pub threads_done: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct ReclaimPlan {
    pub step: ReclaimStep,
    pub siqs_width: usize,
    pub ecm_digits: u32,
    pub workers: usize,
    pub live_bytes: u64,
    pub budget_bytes: u64,
}

impl ReclaimPlan {
    pub fn width_applied(self) -> bool {
        self.step != ReclaimStep::None
    }

    pub fn cap_applied(self) -> bool {
        self.step == ReclaimStep::EcmCap || self.step == ReclaimStep::Threads
    }

    pub fn threads_applied(self) -> bool {
        self.step == ReclaimStep::Threads
    }
}

pub fn next_reclaim_step(live_bytes: u64, budget_bytes: u64, state: ReclaimState) -> ReclaimStep {
    if !crate::options::memory_under_pressure(live_bytes, budget_bytes) {
        return ReclaimStep::None;
    }
    if !state.width_done {
        return ReclaimStep::SiqsWidth;
    }
    if !state.cap_done {
        return ReclaimStep::EcmCap;
    }
    if !state.threads_done {
        return ReclaimStep::Threads;
    }
    ReclaimStep::None
}

pub fn reclaim_siqs_width(current_width: usize, base_width: usize) -> usize {
    let base = base_width.max(1);
    let mut ladder = vec![base];
    for step in crate::siqs::SIQS_WIDTH_ESCALATION_STEPS {
        ladder.push(base.saturating_mul(step));
    }
    ladder.sort_unstable();
    ladder.dedup();
    let mut out = base;
    for entry in ladder {
        if entry < current_width {
            out = entry;
        }
    }
    out
}

pub fn reclaim_ecm_digits(digits: u32) -> u32 {
    digits.saturating_sub(5)
}

pub fn reclaim_thread_count(workers: usize) -> usize {
    (workers / 2).max(1)
}

pub fn plan_joint_reclaim(
    live_bytes: u64,
    budget_bytes: u64,
    base_width: usize,
    current_width: usize,
    digits: u32,
    workers: usize,
) -> ReclaimPlan {
    let mut state = ReclaimState::default();
    let mut width = current_width.max(1);
    let mut ecm_digits = digits;
    let mut thread_count = workers.max(1);
    let mut live = live_bytes;
    let mut step = ReclaimStep::None;
    loop {
        let next = next_reclaim_step(live, budget_bytes, state);
        match next {
            ReclaimStep::None => break,
            ReclaimStep::SiqsWidth => {
                let narrowed = reclaim_siqs_width(width, base_width);
                if narrowed < width {
                    live = live.saturating_mul(narrowed as u64) / (width.max(1) as u64);
                    width = narrowed;
                }
                state.width_done = true;
                step = ReclaimStep::SiqsWidth;
            }
            ReclaimStep::EcmCap => {
                let lowered = reclaim_ecm_digits(ecm_digits);
                let (_, old_two, _) = ecm_plan(ecm_digits);
                let (_, new_two, _) = ecm_plan(lowered);
                if new_two < old_two && old_two > 0 {
                    live = ((live as u128 * new_two as u128) / old_two as u128)
                        .min(u64::MAX as u128) as u64;
                }
                ecm_digits = lowered;
                state.cap_done = true;
                step = ReclaimStep::EcmCap;
            }
            ReclaimStep::Threads => {
                let cut = reclaim_thread_count(thread_count);
                if cut < thread_count {
                    live = live.saturating_mul(cut as u64) / (thread_count.max(1) as u64);
                    thread_count = cut;
                }
                state.threads_done = true;
                step = ReclaimStep::Threads;
            }
        }
    }
    ReclaimPlan {
        step,
        siqs_width: width,
        ecm_digits,
        workers: thread_count,
        live_bytes: live,
        budget_bytes,
    }
}

pub fn reclaim_note(plan: &ReclaimPlan) -> String {
    format!(
        "reclaim={} siqs_width={} ecm_digits={} workers={} live_bytes={} budget_bytes={}",
        plan.step.name(),
        plan.siqs_width,
        plan.ecm_digits,
        plan.workers,
        plan.live_bytes,
        plan.budget_bytes,
    )
}

pub(crate) fn control_joint_memory(
    n: &Integer,
    options: &crate::options::FactorOptions,
    timings: Option<&crate::factor::TimingCollector>,
) -> ReclaimPlan {
    let needed = estimate_siqs_relations_needed(n, options);
    let columns = needed.saturating_sub(1);
    let live = crate::options::estimate_live_bytes(columns, needed, 0, 0);
    let workers = options.effective_workers();
    let width = options.sieve_a_bound;
    let collect = crate::siqs::siqs_collect_estimate_bytes(columns, width as i64, workers);
    let pressure = collect.saturating_add(live);
    let budget = options.memory_budget_bytes();
    crate::factor::record_memory_shadow(
        timings,
        "cascade-reclaim-plan",
        None,
        columns,
        needed,
        0,
        0,
        width as i64,
        workers,
        budget,
        true,
    );
    let base = crate::options::FactorOptions::default().sieve_a_bound;
    plan_joint_reclaim(
        pressure,
        budget,
        base,
        width,
        JOINT_ECM_TARGET_DIGITS,
        workers,
    )
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct JointResume {
    pub ecm: Option<(u32, u64)>,
    pub siqs_next_poly: usize,
    pub siqs_banked: usize,
}

impl JointResume {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn apply_ecm_progress(&mut self, end: (u32, u64)) {
        self.ecm = Some(end);
    }

    pub fn apply_siqs_progress(&mut self, next_poly: usize, new_relations: usize) {
        self.siqs_next_poly = next_poly.max(self.siqs_next_poly);
        self.siqs_banked = self.siqs_banked.saturating_add(new_relations);
    }
}

pub fn joint_worker_split(total_workers: usize, shares: Option<(f64, f64)>) -> (usize, usize) {
    let total = total_workers.max(1);
    let Some((ecm_share, _)) = shares else {
        return (total, total);
    };
    if !ecm_share.is_finite() {
        return (total, total);
    }
    if total == 1 {
        return (1, 1);
    }
    let clamped = ecm_share.clamp(0.0, 1.0);
    let ecm = ((total as f64 * clamped).round() as usize).clamp(1, total - 1);
    (ecm, (total - ecm).max(1))
}

pub(crate) fn joint_budget_to_deadline(budget_ms: Option<u128>) -> Option<Instant> {
    let ms = budget_ms?;
    let duration = Duration::from_millis(u64::try_from(ms).ok()?);
    Instant::now().checked_add(duration)
}

pub fn estimate_siqs_relations_needed(
    n: &Integer,
    options: &crate::options::FactorOptions,
) -> usize {
    let base = crate::base::FactorBase::for_siqs(
        n,
        options.rational_base_bound,
        options.algebraic_base_bound,
    );
    base.column_count().saturating_add(1)
}

pub(crate) struct JointResult {
    pub divisor: Option<Integer>,
    pub winner: Option<&'static str>,
    pub ecm_resume: Option<(u32, u64)>,
    pub ecm_curves: u64,
    pub ecm_stage_one_ms: u128,
    pub ecm_stage_two_ms: u128,
    pub ecm_ms: u128,
    pub ecm_t1: Option<u128>,
    pub ecm_cap: Option<u128>,
    pub ecm_box_ms: Option<u128>,
    pub siqs_bank: crate::siqs::SiqsBank,
    pub ecm_workers: usize,
    pub siqs_workers: usize,
    pub used_fallback: bool,
}

fn verify_joint_divisor(n: &Integer, candidate: &Integer) -> bool {
    candidate > &Integer::from(1) && candidate < n && n % candidate == Integer::from(0)
}

const JOINT_ECM_TARGET_DIGITS: u32 = 15;

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_joint_ecm_siqs(
    n: &Integer,
    signals: &mut CascadeSignals,
    options: &crate::options::FactorOptions,
    cancelled: &std::sync::atomic::AtomicBool,
    deadline: Option<Instant>,
    timings: Option<&crate::factor::TimingCollector>,
    ecm_resume: Option<(u32, u64)>,
    siqs_bank: crate::siqs::SiqsBank,
) -> JointResult {
    let first = run_single_slice(
        n, signals, options, cancelled, deadline, timings, ecm_resume, siqs_bank,
    );
    if first.divisor.is_some() {
        return first;
    }
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return first;
    }
    if cancelled.load(std::sync::atomic::Ordering::Relaxed) {
        return first;
    }
    let relations_needed = estimate_siqs_relations_needed(n, options);
    let steered = governor_slices(signals, signals.target_digits, relations_needed);
    if steered.is_none() || !first.used_fallback {
        return first;
    }
    if first.ecm_resume.is_none() && first.siqs_bank.len() == 0 {
        return first;
    }
    let second = run_single_slice(
        n,
        signals,
        options,
        cancelled,
        deadline,
        timings,
        first.ecm_resume,
        first.siqs_bank.clone(),
    );
    if second.divisor.is_some() {
        return second;
    }
    if second.ecm_curves == 0 && second.siqs_bank.len() == first.siqs_bank.len() {
        return first;
    }
    second
}

#[allow(clippy::too_many_arguments)]
fn run_single_slice(
    n: &Integer,
    signals: &mut CascadeSignals,
    options: &crate::options::FactorOptions,
    cancelled: &std::sync::atomic::AtomicBool,
    deadline: Option<Instant>,
    timings: Option<&crate::factor::TimingCollector>,
    ecm_resume: Option<(u32, u64)>,
    siqs_bank: crate::siqs::SiqsBank,
) -> JointResult {
    let reclaim = control_joint_memory(n, options, timings);
    if let Some(collector) = timings {
        collector.record_elapsed("reclaim", None, 0, reclaim_note(&reclaim));
    }
    let workers = reclaim.workers;
    let relations_needed = estimate_siqs_relations_needed(n, options);
    let shares = governor_slices(signals, signals.target_digits, relations_needed);
    let (ecm_workers, siqs_workers) = joint_worker_split(workers, shares);
    let used_fallback = shares.is_none();
    let live_remaining = signals.remaining_ms(deadline);
    let legacy_ecm_budget = allocate_slices(signals, live_remaining).3;
    let incoming_banked = siqs_bank.len();
    let ecm_box_opt = if used_fallback {
        legacy_ecm_budget
    } else {
        live_remaining
    };
    let ecm_box_opt = match ecm_box_opt {
        Some(boxed) if reclaim.cap_applied() => Some(boxed / 2),
        other => other,
    };
    let ecm_digits = reclaim.ecm_digits;
    let ecm_deadline = if used_fallback {
        joint_budget_to_deadline(legacy_ecm_budget)
    } else {
        deadline
    };
    if used_fallback {
        let (ecm_report, ecm_t1, ecm_cap, ecm_box_ms) = ecm_tranche_with_resume(
            n,
            ecm_deadline,
            ecm_box_opt,
            live_remaining,
            ecm_workers,
            ecm_resume,
            ecm_digits,
        );
        return finish_joint_arms(
            n,
            signals,
            ecm_report,
            ecm_t1,
            ecm_cap,
            ecm_box_ms,
            None,
            siqs_bank,
            incoming_banked,
            ecm_workers,
            siqs_workers,
            used_fallback,
        );
    }
    let siqs_deadline = deadline;
    let mut siqs_options = options.clone();
    siqs_options.threads = Some(siqs_workers);
    siqs_options.sieve_a_bound = reclaim.siqs_width;
    let run_siqs = options.use_mpqs;
    let ecm_job_deadline = ecm_deadline;
    let ecm_job_box = ecm_box_opt;
    let ecm_job_live = live_remaining;
    let (ecm_report, ecm_t1, ecm_cap, ecm_box_ms, siqs_hit, siqs_bank_out) =
        std::thread::scope(|scope| {
            let ecm_handle = scope.spawn(move || {
                ecm_tranche_with_resume(
                    n,
                    ecm_job_deadline,
                    ecm_job_box,
                    ecm_job_live,
                    ecm_workers,
                    ecm_resume,
                    ecm_digits,
                )
            });
            let siqs_handle = scope.spawn(move || {
                if !run_siqs {
                    return (None, siqs_bank);
                }
                let mut nop = |_: &crate::factor::ProgressEvent| Ok(());
                crate::siqs::siqs_pilot_with_bank(
                    n,
                    &siqs_options,
                    cancelled,
                    siqs_deadline,
                    &mut nop,
                    timings,
                    siqs_bank,
                )
            });
            let (ecm_report, ecm_t1, ecm_cap, ecm_box_ms) =
                ecm_handle.join().unwrap_or_else(|_| {
                    (
                        crate::ecm::ecm_factor_report(
                            n,
                            JOINT_ECM_TARGET_DIGITS,
                            Some(0),
                            ecm_job_deadline,
                            ecm_workers,
                            ecm_resume,
                        ),
                        None,
                        None,
                        None,
                    )
                });
            let (siqs_hit, siqs_bank_out) = siqs_handle
                .join()
                .unwrap_or((None, crate::siqs::SiqsBank::new()));
            (
                ecm_report,
                ecm_t1,
                ecm_cap,
                ecm_box_ms,
                siqs_hit,
                siqs_bank_out,
            )
        });
    finish_joint_arms(
        n,
        signals,
        ecm_report,
        ecm_t1,
        ecm_cap,
        ecm_box_ms,
        siqs_hit,
        siqs_bank_out,
        incoming_banked,
        ecm_workers,
        siqs_workers,
        used_fallback,
    )
}

#[allow(clippy::too_many_arguments)]
fn finish_joint_arms(
    n: &Integer,
    signals: &mut CascadeSignals,
    ecm_report: crate::ecm::EcmReport,
    ecm_t1: Option<u128>,
    ecm_cap: Option<u128>,
    ecm_box_ms: Option<u128>,
    siqs_hit: Option<(Integer, crate::siqs::SiqsOutcome)>,
    siqs_bank: crate::siqs::SiqsBank,
    incoming_banked: usize,
    ecm_workers: usize,
    siqs_workers: usize,
    used_fallback: bool,
) -> JointResult {
    signals.report_ecm(
        ecm_report.curves_run,
        ecm_report.elapsed_ms,
        ecm_report.divisor.is_some(),
    );
    let fresh_relations = siqs_bank.len().saturating_sub(incoming_banked);
    if fresh_relations > 0 {
        signals.report_pilot(siqs_bank.len(), ecm_report.elapsed_ms.max(1));
    }
    let ecm_divisor = ecm_report
        .divisor
        .clone()
        .filter(|candidate| verify_joint_divisor(n, candidate));
    let siqs_divisor = siqs_hit.and_then(|(candidate, outcome)| {
        if verify_joint_divisor(n, &candidate) {
            signals.report_pilot(outcome.relations.max(1), 1);
            Some(candidate)
        } else {
            None
        }
    });
    let (divisor, winner) = match (ecm_divisor, siqs_divisor) {
        (Some(ecm), _) => (Some(ecm), Some("ecm")),
        (None, Some(siqs)) => (Some(siqs), Some("siqs")),
        (None, None) => (None, None),
    };
    let ecm_resume = Some(crate::ecm::ecm_checkpoint(&ecm_report));
    let ecm_curves = ecm_report.curves_run;
    JointResult {
        divisor,
        winner,
        ecm_resume,
        ecm_curves,
        ecm_stage_one_ms: ecm_report.stage_one_ms,
        ecm_stage_two_ms: ecm_report.stage_two_ms,
        ecm_ms: ecm_report.elapsed_ms,
        ecm_t1,
        ecm_cap,
        ecm_box_ms,
        siqs_bank,
        ecm_workers,
        siqs_workers,
        used_fallback,
    }
}

fn ecm_tranche_with_resume(
    n: &Integer,
    ecm_deadline: Option<Instant>,
    ecm_box_opt: Option<u128>,
    live_remaining: Option<u128>,
    ecm_workers: usize,
    ecm_resume: Option<(u32, u64)>,
    target_digits: u32,
) -> (
    crate::ecm::EcmReport,
    Option<u128>,
    Option<u128>,
    Option<u128>,
) {
    let ecm_start = Instant::now();
    let ecm_box_ms =
        ecm_deadline.map(|limit| limit.saturating_duration_since(ecm_start).as_millis());
    let first_box = tranche_first_box(ecm_box_opt, live_remaining);
    let first_slice =
        match (ecm_deadline, first_box) {
            (Some(limit), Some(first)) => Some(limit.min(
                Instant::now() + Duration::from_millis(u64::try_from(first).unwrap_or(u64::MAX)),
            )),
            _ => ecm_deadline,
        };
    let first =
        crate::ecm::ecm_factor_report(n, target_digits, None, first_slice, ecm_workers, ecm_resume);
    let first_ms = first.elapsed_ms;
    if first.divisor.is_some() {
        return (first, first_box, None, ecm_box_ms);
    }
    if first_box.is_none_or(|boxed| first_ms < boxed) {
        return (first, first_box, None, ecm_box_ms);
    }
    let live_two =
        ecm_deadline.map(|limit| limit.saturating_duration_since(Instant::now()).as_millis());
    let spent = ecm_start.elapsed().as_millis();
    let remain_box = match (ecm_box_opt, live_two) {
        (Some(box_ms), Some(live)) => Some(box_ms.saturating_sub(spent).min(live)),
        (Some(box_ms), None) => Some(box_ms.saturating_sub(spent)),
        _ => None,
    };
    let mut local = CascadeSignals::collect(n, None);
    if first.curves_run > 0 {
        local.failed_curves = first.curves_run;
    }
    if first.curves_run > 0 {
        local.measured_ms_per_curve = Some(first.elapsed_ms / u128::from(first.curves_run).max(1));
    }
    let capped_two = cap_ecm_slice(&local, live_two, remain_box);
    let cap_note = match (capped_two, remain_box) {
        (Some(cap), Some(remain)) if cap < remain => Some(cap),
        _ => None,
    };
    let second_slice =
        match (ecm_deadline, capped_two) {
            (Some(limit), Some(cap)) => Some(limit.min(
                Instant::now() + Duration::from_millis(u64::try_from(cap).unwrap_or(u64::MAX)),
            )),
            _ => ecm_deadline,
        };
    let second = crate::ecm::ecm_factor_report(
        n,
        target_digits,
        None,
        second_slice,
        ecm_workers,
        Some(crate::ecm::ecm_checkpoint(&first)),
    );
    let combined = combine_ecm_reports(first, second);
    (combined, first_box, cap_note, ecm_box_ms)
}

fn combine_ecm_reports(
    first: crate::ecm::EcmReport,
    second: crate::ecm::EcmReport,
) -> crate::ecm::EcmReport {
    let divisor = second.divisor.or(first.divisor);
    crate::ecm::EcmReport {
        divisor,
        curves_run: first.curves_run.saturating_add(second.curves_run),
        elapsed_ms: first.elapsed_ms.saturating_add(second.elapsed_ms),
        stage_one_ms: first.stage_one_ms.saturating_add(second.stage_one_ms),
        stage_two_ms: first.stage_two_ms.saturating_add(second.stage_two_ms),
        end_digits: second.end_digits,
        end_curve: second.end_curve,
    }
}

pub(crate) fn fermat_iter_budget(ceiling: u64, slice_ms: Option<u128>, ns_per_iter: u128) -> u64 {
    let Some(slice) = slice_ms else {
        return ceiling;
    };
    let affordable = slice.saturating_mul(1_000_000) / ns_per_iter.max(1);
    ceiling.min(affordable.min(u64::MAX as u128) as u64)
}

#[cfg(test)]
mod cascade_tests {
    use super::fermat_iter_budget;
    use crate::cascade::{
        CascadeSignals, allocate_slices, cheap_ecm_plan, cheap_slice_ms, ecm_plan, escalate_digits,
        project_ecm_cost_ms, project_ecm_cost_ms_for_digits, project_siqs_cost_ms,
        projected_lower_cost, tranche_first_box,
    };
    use crate::integer::Integer;

    #[test]
    fn tight_slices_fund_ecm_majority_with_bounded_cheap_tier() {
        let signals = CascadeSignals::collect(&(Integer::from(1_u8) << 267), None);
        let (fermat, pm1, pp1, ecm) = allocate_slices(&signals, Some(1000));
        assert!(ecm.expect("ecm slice") >= 700);
        assert!(fermat.expect("fermat slice") <= 40);
        assert!(pm1.expect("pm1 slice") <= 70);
        assert!(pp1.expect("pp1 slice") <= 70);
        let cheap = cheap_slice_ms(Some(1000)).expect("cheap slice");
        assert!(cheap <= 70);
        let (bound_one, bound_two, curves) = cheap_ecm_plan(Some(1000));
        assert!(bound_one <= 11_000);
        assert!(bound_two <= 1_900_000);
        assert!((4..=32).contains(&curves));
    }

    #[test]
    fn cheap_tier_curves_derive_from_remaining_deadline() {
        let (_, _, tiny) = cheap_ecm_plan(Some(100));
        let (_, _, small) = cheap_ecm_plan(Some(1000));
        let (_, _, ample) = cheap_ecm_plan(Some(60_000));
        let (_, _, unbounded) = cheap_ecm_plan(None);
        assert_eq!(tiny, 4);
        assert!(small > tiny);
        assert!(small <= ample);
        assert_eq!(ample, 32);
        assert_eq!(unbounded, 32);
    }

    #[test]
    fn iter_budget_respects_ceiling_and_slice() {
        assert_eq!(fermat_iter_budget(200_000, None, 1_000), 200_000);
        assert_eq!(fermat_iter_budget(200_000, Some(60_000), 1), 200_000);
        assert_eq!(fermat_iter_budget(200_000, Some(10), 1_000), 10_000);
        assert_eq!(fermat_iter_budget(200_000, Some(0), 1_000), 0);
    }

    #[test]
    fn siqs_projection_requires_measurable_inputs() {
        assert_eq!(project_siqs_cost_ms(0, 1.0), None);
        assert_eq!(project_siqs_cost_ms(100, 0.0), None);
        assert_eq!(project_siqs_cost_ms(100, f64::NAN), None);
        assert_eq!(project_siqs_cost_ms(100, f64::NEG_INFINITY), None);
        assert_eq!(project_siqs_cost_ms(3000, 3.0), Some(1000));
        assert_eq!(project_siqs_cost_ms(3001, 3.0), Some(1001));
    }

    #[test]
    fn lower_cost_tags_projection_winner() {
        assert_eq!(projected_lower_cost(Some(10), Some(20)), Some("ecm"));
        assert_eq!(projected_lower_cost(Some(20), Some(10)), Some("siqs"));
        assert_eq!(projected_lower_cost(Some(10), Some(10)), Some("tie"));
        assert_eq!(projected_lower_cost(None, Some(10)), None);
        assert_eq!(projected_lower_cost(Some(10), None), None);
        assert_eq!(projected_lower_cost(None, None), None);
    }

    #[test]
    fn ecm_projection_scales_with_digits() {
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        assert_eq!(signals.target_digits, 30);
        signals.measured_ms_per_curve = Some(10);
        let at_target = project_ecm_cost_ms(&signals).expect("calibrated projection");
        let escalated =
            project_ecm_cost_ms_for_digits(&signals, escalate_digits(signals.target_digits))
                .expect("escalated projection");
        assert_eq!(
            at_target,
            u128::from(ecm_plan(signals.target_digits).2) * 10
        );
        assert!(escalated >= at_target);
        assert_eq!(
            escalated,
            u128::from(ecm_plan(escalate_digits(signals.target_digits)).2) * 10
        );
    }

    #[test]
    fn tranche_box_clamps_to_live_remaining() {
        assert_eq!(tranche_first_box(Some(700), Some(800)), Some(350));
        assert_eq!(tranche_first_box(Some(700), Some(400)), Some(200));
        assert_eq!(tranche_first_box(Some(700), None), Some(350));
        assert_eq!(tranche_first_box(None, Some(400)), None);
    }

    #[test]
    fn profiler_carries_live_device_signals() {
        use std::time::{Duration, Instant};
        let n = Integer::from(1_u8) << 100;
        let deadline = Some(Instant::now() + Duration::from_secs(60));
        let mut signals = CascadeSignals::collect(&n, deadline);
        assert_eq!(signals.bits, n.bits());
        assert!(signals.workers >= 1);
        assert!(signals.measured_ms_per_curve.is_none());
        assert!(signals.pilot_relations_per_ms.is_none());
        assert!(signals.pilot_relations_per_sec.is_none());
        signals.report_ecm(4, 80, false);
        assert_eq!(signals.measured_ms_per_curve, Some(20));
        signals.report_pilot(300, 100);
        assert_eq!(signals.pilot_relations_per_ms, Some(3.0));
        assert_eq!(signals.pilot_relations_per_sec, Some(3000.0));
        assert_eq!(signals.measured_relations_per_sec(), Some(3000.0));
        signals.observe_trial_division(168);
        assert_eq!(signals.trial_division_primes_checked, Some(168));
        signals.observe_workers(8);
        assert_eq!(signals.workers, 8);
        assert!(signals.remaining_ms_at_collect.is_some());
        assert!(signals.remaining_ms(deadline).is_some());
        let probe = CascadeSignals::probe_budget_ms(&n);
        assert!((1000..=3_600_000).contains(&probe));
        let effective = signals.effective_remaining_ms(&n, deadline);
        assert!(effective > 0);
        let mut legacy = CascadeSignals::collect(&n, None);
        legacy.pilot_relations_per_ms = Some(2.5);
        assert_eq!(legacy.pilot_relations_per_sec, None);
        assert_eq!(legacy.measured_relations_per_sec(), Some(2500.0));
        legacy.observe_workers(0);
        assert_eq!(legacy.workers, 1);
    }

    #[test]
    fn profiler_honors_options_budget() {
        use crate::options::FactorOptions;
        use std::time::{Duration, Instant};
        let n = Integer::from(1_u8) << 100;
        let options = FactorOptions {
            threads: Some(4),
            max_duration: Some(Duration::from_secs(30)),
            ..FactorOptions::default()
        };
        assert_eq!(options.effective_workers(), 4);
        let from = Instant::now();
        assert!(options.effective_deadline(from).is_some());
        assert!(FactorOptions::default().effective_deadline(from).is_none());
        assert!(FactorOptions::default().effective_workers() >= 1);
        let budgeted = CascadeSignals::collect_with_options(&n, &options, None);
        assert_eq!(budgeted.workers, 4);
        let snapshot = budgeted.remaining_ms_at_collect.expect("budget snapshot");
        assert!((29_000..=30_000).contains(&snapshot));
        let explicit = Some(Instant::now() + Duration::from_secs(60));
        let mut overridden = CascadeSignals::collect(&n, explicit);
        overridden.observe_options(&options, explicit);
        assert_eq!(overridden.workers, 4);
        assert!(overridden.remaining_ms_at_collect.is_some());
    }

    #[test]
    fn cost_model_prefers_ecm_when_curve_rate_is_fast() {
        use crate::cascade::{
            cost_shares, project_ecm_remaining_ms, project_siqs_remaining_ms, projected_lower_cost,
            should_switch_to_siqs,
        };
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        signals.measured_ms_per_curve = Some(2);
        signals.report_pilot(100, 1000);
        let relations_needed = 100_000_usize;
        let ecm = project_ecm_remaining_ms(&signals, signals.target_digits).expect("ecm known");
        let siqs = project_siqs_remaining_ms(&signals, relations_needed).expect("siqs known");
        assert!(ecm < siqs);
        assert_eq!(projected_lower_cost(Some(ecm), Some(siqs)), Some("ecm"));
        let shares = cost_shares(Some(ecm), Some(siqs)).expect("shares known");
        assert!(shares.0 > shares.1);
        assert!((shares.0 + shares.1 - 1.0).abs() < 1e-9);
        assert_eq!(
            should_switch_to_siqs(&signals, relations_needed),
            Some(false)
        );
    }

    #[test]
    fn cost_model_prefers_siqs_when_pilot_rate_is_fast() {
        use crate::cascade::{
            cost_shares, project_ecm_remaining_ms, project_siqs_remaining_ms, projected_lower_cost,
            should_switch_to_siqs,
        };
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        signals.measured_ms_per_curve = Some(500);
        signals.report_pilot(3000, 100);
        let relations_needed = 3000_usize;
        let ecm = project_ecm_remaining_ms(&signals, signals.target_digits).expect("ecm known");
        let siqs = project_siqs_remaining_ms(&signals, relations_needed).expect("siqs known");
        assert!(siqs < ecm);
        assert_eq!(projected_lower_cost(Some(ecm), Some(siqs)), Some("siqs"));
        let shares = cost_shares(Some(ecm), Some(siqs)).expect("shares known");
        assert!(shares.1 > shares.0);
        assert!((shares.0 + shares.1 - 1.0).abs() < 1e-9);
        assert_eq!(
            should_switch_to_siqs(&signals, relations_needed),
            Some(true)
        );
    }

    #[test]
    fn cost_model_unknown_without_measurements() {
        use crate::cascade::{
            cost_shares, project_ecm_remaining_ms, project_siqs_remaining_ms, should_switch_to_siqs,
        };
        let signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        assert_eq!(
            project_ecm_remaining_ms(&signals, signals.target_digits),
            None
        );
        assert_eq!(project_siqs_remaining_ms(&signals, 3000), None);
        assert_eq!(cost_shares(None, Some(100)), None);
        assert_eq!(cost_shares(Some(100), None), None);
        assert_eq!(should_switch_to_siqs(&signals, 3000), None);
        let mut zero_curve = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        zero_curve.measured_ms_per_curve = Some(0);
        assert_eq!(
            project_ecm_remaining_ms(&zero_curve, zero_curve.target_digits),
            None
        );
        let mut rated = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        rated.report_pilot(3000, 100);
        assert_eq!(project_siqs_remaining_ms(&rated, 0), None);
    }

    #[test]
    fn ecm_posterior_escalates_past_failed_curves() {
        use crate::cascade::{escalate_digits, project_ecm_remaining_ms};
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        signals.target_digits = 25;
        signals.measured_ms_per_curve = Some(10);
        signals.failed_curves = 0;
        assert_eq!(project_ecm_remaining_ms(&signals, 25), Some(2140));
        signals.failed_curves = 214;
        assert_eq!(project_ecm_remaining_ms(&signals, 25), Some(4280));
        signals.failed_curves = 642;
        assert_eq!(project_ecm_remaining_ms(&signals, 25), Some(12_840));
        assert_eq!(escalate_digits(25), 30);
    }

    #[test]
    fn cost_shares_split_ties_evenly() {
        use crate::cascade::cost_shares;
        assert_eq!(cost_shares(Some(10), Some(10)), Some((0.5, 0.5)));
        assert_eq!(cost_shares(Some(0), Some(0)), Some((0.5, 0.5)));
    }

    #[test]
    fn governor_slices_follow_measured_cost_shares() {
        use crate::cascade::{
            cost_shares, governor_slices, project_ecm_remaining_ms, project_siqs_remaining_ms,
        };
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        signals.measured_ms_per_curve = Some(10);
        signals.report_pilot(3000, 100);
        let digits = signals.target_digits;
        let relations_needed = 3000_usize;
        let ecm = project_ecm_remaining_ms(&signals, digits).expect("ecm known");
        let siqs = project_siqs_remaining_ms(&signals, relations_needed).expect("siqs known");
        let expected = cost_shares(Some(ecm), Some(siqs)).expect("shares known");
        let got = governor_slices(&signals, digits, relations_needed).expect("governor known");
        assert!((got.0 - expected.0).abs() < 1e-12);
        assert!((got.1 - expected.1).abs() < 1e-12);
        assert!((got.0 + got.1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn governor_slices_fall_back_to_legacy_when_unknown() {
        use crate::cascade::{allocate_slices, governor_slices};
        let fresh = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        assert_eq!(governor_slices(&fresh, fresh.target_digits, 3000), None);
        let mut ecm_only = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        ecm_only.measured_ms_per_curve = Some(10);
        assert_eq!(
            governor_slices(&ecm_only, ecm_only.target_digits, 3000),
            None
        );
        let mut siqs_only = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        siqs_only.report_pilot(3000, 100);
        assert_eq!(
            governor_slices(&siqs_only, siqs_only.target_digits, 3000),
            None
        );
        let remaining = Some(96_000_u128);
        assert_eq!(
            allocate_slices(&fresh, remaining),
            (
                Some(96_000 / 64),
                Some(96_000 / 16),
                Some(96_000 / 16),
                Some(96_000 * 3 / 4),
            )
        );
        let fallback = allocate_slices(&ecm_only, remaining);
        assert!(fallback.0.is_some());
        assert!(fallback.3.is_some());
    }

    #[test]
    fn joint_split_follows_governor_shares() {
        use crate::cascade::{governor_slices, joint_worker_split};
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 100), None);
        signals.measured_ms_per_curve = Some(10);
        signals.report_pilot(3000, 100);
        let digits = signals.target_digits;
        let shares = governor_slices(&signals, digits, 3000).expect("governor known");
        let (ecm_workers, siqs_workers) = joint_worker_split(8, Some(shares));
        let expect_ecm = (8.0 * shares.0).round() as usize;
        assert_eq!(ecm_workers, expect_ecm.clamp(1, 7));
        assert_eq!(siqs_workers, 8 - ecm_workers);
        assert_eq!(ecm_workers + siqs_workers, 8);
        let (fallback_ecm, fallback_siqs) = joint_worker_split(8, None);
        assert_eq!((fallback_ecm, fallback_siqs), (8, 8));
    }

    #[test]
    fn joint_resume_preserves_progress() {
        use crate::cascade::JointResume;
        let mut resume = JointResume::new();
        resume.apply_ecm_progress((12, 16));
        resume.apply_siqs_progress(9, 5);
        assert_eq!(resume.ecm, Some((12, 16)));
        assert_eq!(resume.siqs_next_poly, 9);
        assert_eq!(resume.siqs_banked, 5);
        resume.apply_siqs_progress(12, 3);
        assert_eq!(resume.siqs_banked, 8);
        assert_eq!(resume.siqs_next_poly, 12);
        assert_eq!(resume.ecm, Some((12, 16)));
        let next = resume.ecm.expect("ecm resume");
        assert_eq!(next.1, 16);
    }

    #[test]
    fn fallback_slice_runs_ecm_only_without_siqs_work() {
        use crate::cascade::run_joint_ecm_siqs;
        use crate::factor::TimingCollector;
        use crate::options::FactorOptions;
        use crate::siqs::SiqsBank;
        use std::sync::atomic::AtomicBool;
        use std::time::Duration;
        let n = Integer::from(1000003) * Integer::from(1000033);
        let options = FactorOptions {
            max_duration: Some(Duration::from_millis(50)),
            ..FactorOptions::default()
        };
        let deadline = options.effective_deadline(std::time::Instant::now());
        let mut signals = CascadeSignals::collect_with_options(&n, &options, deadline);
        let timings = TimingCollector::default();
        let cancelled = AtomicBool::new(false);
        let joint = run_joint_ecm_siqs(
            &n,
            &mut signals,
            &options,
            &cancelled,
            deadline,
            Some(&timings),
            None,
            SiqsBank::new(),
        );
        assert!(joint.used_fallback);
        assert_eq!(joint.siqs_bank.len(), 0);
        let stages = timings.snapshot();
        assert!(
            stages.iter().all(|entry| !entry.stage.starts_with("siqs")),
            "fallback must not run siqs work"
        );
        assert!(
            stages.iter().any(|entry| entry.stage == "reclaim"),
            "reclaim decision must be recorded"
        );
        assert!(joint.winner.is_none_or(|winner| winner == "ecm"));
    }

    #[test]
    fn forced_pressure_reclaims_through_threads() {
        use crate::cascade::run_joint_ecm_siqs;
        use crate::factor::TimingCollector;
        use crate::options::FactorOptions;
        use crate::siqs::SiqsBank;
        use std::sync::atomic::AtomicBool;
        use std::time::Duration;
        let n = Integer::from(1000003) * Integer::from(1000033);
        let options = FactorOptions {
            max_duration: Some(Duration::from_millis(50)),
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        let deadline = options.effective_deadline(std::time::Instant::now());
        let mut signals = CascadeSignals::collect_with_options(&n, &options, deadline);
        let full_workers = options.effective_workers();
        let timings = TimingCollector::default();
        let cancelled = AtomicBool::new(false);
        let joint = run_joint_ecm_siqs(
            &n,
            &mut signals,
            &options,
            &cancelled,
            deadline,
            Some(&timings),
            None,
            SiqsBank::new(),
        );
        assert!(joint.used_fallback);
        assert_eq!(joint.ecm_workers, (full_workers / 2).max(1));
        assert_eq!(joint.siqs_workers, (full_workers / 2).max(1));
        let notes: Vec<String> = timings
            .snapshot()
            .iter()
            .filter(|entry| entry.stage == "reclaim")
            .map(|entry| entry.detail.clone())
            .collect();
        assert_eq!(notes.len(), 1);
        assert!(
            notes[0].contains("reclaim=threads"),
            "forced pressure must reclaim through threads, got {}",
            notes[0]
        );
        assert!(joint.winner.is_none_or(|winner| winner == "ecm"));
    }

    #[test]
    fn ecm_miss_cap_leaves_sieve_spillover() {
        use crate::cascade::cap_ecm_slice;
        let mut signals = CascadeSignals::collect(&(Integer::from(1_u8) << 90), None);
        signals.failed_curves = 32;
        signals.measured_ms_per_curve = Some(2);
        let capped = cap_ecm_slice(&signals, Some(800), Some(750));
        assert_eq!(capped, Some(400));
        let live = 800_u128;
        let kept = capped.expect("capped ecm slice");
        assert!(kept < 750, "miss burn must shrink the ecm box");
        assert!(
            live - kept >= live / 2,
            "spillover must keep a live-relative sieve floor"
        );
    }

    #[test]
    fn ecm_win_path_stays_uncapped() {
        use crate::cascade::cap_ecm_slice;
        let mut demo = CascadeSignals::collect(&(Integer::from(1_u8) << 267), None);
        demo.failed_curves = 14;
        demo.measured_ms_per_curve = Some(0);
        assert_eq!(cap_ecm_slice(&demo, Some(900), Some(750)), Some(750));
        let mut fresh = CascadeSignals::collect(&(Integer::from(1_u8) << 267), None);
        fresh.measured_ms_per_curve = None;
        assert_eq!(cap_ecm_slice(&fresh, Some(900), Some(750)), Some(750));
        let mut trickle = CascadeSignals::collect(&(Integer::from(1_u8) << 90), None);
        trickle.failed_curves = 2;
        trickle.measured_ms_per_curve = Some(1);
        assert_eq!(cap_ecm_slice(&trickle, Some(800), Some(750)), Some(750));
        assert_eq!(cap_ecm_slice(&demo, None, Some(750)), None);
    }

    #[test]
    fn reclaim_order_is_width_then_cap_then_threads() {
        use crate::cascade::{ReclaimState, ReclaimStep, next_reclaim_step};
        let budget = 1_000_u64;
        let over = 2_000_u64;
        let fresh = ReclaimState::default();
        assert_eq!(next_reclaim_step(budget, budget, fresh), ReclaimStep::None);
        assert_eq!(next_reclaim_step(500, budget, fresh), ReclaimStep::None);
        assert_eq!(
            next_reclaim_step(over, budget, fresh),
            ReclaimStep::SiqsWidth
        );
        let after_width = ReclaimState {
            width_done: true,
            ..fresh
        };
        assert_eq!(
            next_reclaim_step(over, budget, after_width),
            ReclaimStep::EcmCap
        );
        let after_cap = ReclaimState {
            cap_done: true,
            ..after_width
        };
        assert_eq!(
            next_reclaim_step(over, budget, after_cap),
            ReclaimStep::Threads
        );
        let exhausted = ReclaimState {
            threads_done: true,
            ..after_cap
        };
        assert_eq!(
            next_reclaim_step(over, budget, exhausted),
            ReclaimStep::None
        );
        assert_eq!(next_reclaim_step(500, budget, after_cap), ReclaimStep::None);
    }

    #[test]
    fn reclaim_plan_stages_match_live_pressure() {
        use crate::cascade::{ReclaimStep, plan_joint_reclaim};
        let budget = 1_000_u64;
        let calm = plan_joint_reclaim(500, budget, 700, 700, 25, 8);
        assert_eq!(calm.step, ReclaimStep::None);
        assert_eq!(calm.siqs_width, 700);
        assert_eq!(calm.ecm_digits, 25);
        assert_eq!(calm.workers, 8);
        let storm = plan_joint_reclaim(100_000, budget, 700, 11_200, 30, 8);
        assert_eq!(storm.step, ReclaimStep::Threads);
        assert_eq!(storm.siqs_width, 2_800);
        assert_eq!(storm.ecm_digits, 25);
        assert_eq!(storm.workers, 4);
        let narrow = plan_joint_reclaim(1_100, budget, 700, 2_800, 25, 8);
        assert_eq!(narrow.step, ReclaimStep::SiqsWidth);
        assert_eq!(narrow.siqs_width, 700);
        assert_eq!(narrow.ecm_digits, 25);
        assert_eq!(narrow.workers, 8);
    }

    #[test]
    fn reclaim_steps_reuse_existing_tables() {
        use crate::cascade::{escalate_digits, reclaim_ecm_digits, reclaim_siqs_width};
        use crate::siqs::SIQS_WIDTH_ESCALATION_STEPS;
        assert_eq!(
            SIQS_WIDTH_ESCALATION_STEPS
                .iter()
                .map(|step| 700_usize.saturating_mul(*step))
                .collect::<Vec<_>>(),
            vec![2_800, 11_200]
        );
        assert_eq!(reclaim_siqs_width(11_200, 700), 2_800);
        assert_eq!(reclaim_siqs_width(2_800, 700), 700);
        assert_eq!(reclaim_siqs_width(700, 700), 700);
        assert_eq!(reclaim_ecm_digits(escalate_digits(25)), 25);
        assert_eq!(reclaim_ecm_digits(15), 10);
    }
}
