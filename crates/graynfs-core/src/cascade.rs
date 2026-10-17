use crate::integer::Integer;
use std::time::Instant;

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
}

impl CascadeSignals {
    pub fn collect(n: &Integer, deadline: Option<Instant>) -> Self {
        let _ = deadline;
        let bits = n.bits();
        let digits = (bits * 301 / 1000).min(u64::from(u32::MAX)).max(1) as u32;
        Self {
            bits,
            target_digits: digits,
            failed_curves: 0,
            measured_ms_per_curve: None,
            pilot_relations_per_ms: None,
            trial_division_primes_checked: None,
        }
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
}
