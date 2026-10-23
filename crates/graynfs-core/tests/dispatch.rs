use graynfs_core::cascade::{CascadeSignals, allocate_slices, fermat_cap_for};
use graynfs_core::{FactorOptions, factor, parse_decimal};
use std::time::Duration;

fn uncalibrated() -> CascadeSignals {
    CascadeSignals {
        bits: 60,
        target_digits: 18,
        failed_curves: 0,
        measured_ms_per_curve: None,
        pilot_relations_per_ms: None,
        trial_division_primes_checked: None,
        workers: 1,
        pilot_relations_per_sec: None,
        remaining_ms_at_collect: None,
    }
}

fn calibrated(target_digits: u32, ms_per_curve: u128) -> CascadeSignals {
    CascadeSignals {
        bits: 60,
        target_digits,
        failed_curves: 0,
        measured_ms_per_curve: Some(ms_per_curve),
        pilot_relations_per_ms: None,
        trial_division_primes_checked: None,
        workers: 1,
        pilot_relations_per_sec: None,
        remaining_ms_at_collect: None,
    }
}

#[test]
fn test_uncalibrated_matches_legacy_fractions() {
    let signals = uncalibrated();
    let remaining = 96_000_u128;
    let slices = allocate_slices(&signals, Some(remaining));
    assert_eq!(
        slices,
        (
            Some(remaining / 32),
            Some(remaining / 16),
            Some(remaining / 16),
            Some(remaining * 3 / 4),
        )
    );
    assert_eq!(allocate_slices(&signals, None), (None, None, None, None));
}

#[test]
fn test_slice_shrinks_when_ecm_calibrated() {
    let remaining = 96_000_u128;
    let legacy = allocate_slices(&uncalibrated(), Some(remaining));
    let tuned = allocate_slices(&calibrated(20, 10), Some(remaining));
    let (legacy_fermat, _, _, legacy_ecm) = legacy;
    let (tuned_fermat, _, _, tuned_ecm) = tuned;
    let legacy_fermat = legacy_fermat.expect("legacy fermat budget");
    let legacy_ecm = legacy_ecm.expect("legacy ecm budget");
    let tuned_fermat = tuned_fermat.expect("tuned fermat budget");
    let tuned_ecm = tuned_ecm.expect("tuned ecm budget");
    assert!(tuned_fermat < legacy_fermat);
    assert!(tuned_ecm > legacy_ecm);
    assert_eq!(tuned_fermat + tuned_ecm, legacy_fermat + legacy_ecm);
    assert_eq!(
        allocate_slices(&calibrated(20, 10), None),
        (None, None, None, None)
    );
}

#[test]
fn test_fermat_cap_preserves_legit_wins() {
    let deg3 = parse_decimal("1155285037654380637").expect("valid deg3 input");
    let deg5 = parse_decimal("4606404410723").expect("valid deg5 input");
    let tiny = parse_decimal("15").expect("valid tiny input");
    for value in [&deg3, &deg5, &tiny] {
        let cap = fermat_cap_for(value);
        assert!(cap >= 2330, "cap {cap} must preserve legit wins");
        assert!(cap <= 200_000, "cap {cap} must keep absolute ceiling");
    }
}

fn pilot_rate(report: &graynfs_core::FactorReport) -> Option<f64> {
    report.timings.iter().find_map(|timing| {
        if timing.stage != "pilot_rate" {
            return None;
        }
        timing
            .detail
            .split("rate=")
            .nth(1)
            .and_then(|tail| tail.split_whitespace().next())
            .and_then(|value| value.parse::<f64>().ok())
    })
}

#[test]
fn test_fermat_slice_shrinks_with_target_digits() {
    let remaining = 96_000_u128;
    let mut small = uncalibrated();
    small.target_digits = 18;
    assert_eq!(
        allocate_slices(&small, Some(remaining)).0,
        Some(remaining / 32)
    );
    let mut boundary = uncalibrated();
    boundary.target_digits = 25;
    assert_eq!(
        allocate_slices(&boundary, Some(remaining)).0,
        Some(remaining / 32)
    );
    let mut large = uncalibrated();
    large.target_digits = 26;
    assert_eq!(
        allocate_slices(&large, Some(remaining)).0,
        Some(remaining / 64)
    );
    let mut demo = uncalibrated();
    demo.target_digits = 80;
    assert_eq!(
        allocate_slices(&demo, Some(remaining)).0,
        Some(remaining / 64)
    );
}

#[test]
fn test_pilot_rate_recorded_in_timings() {
    let n = parse_decimal("11783510925019").expect("valid 44-bit input");
    let options = FactorOptions {
        skip_preprocessing: true,
        max_duration: Some(Duration::from_secs(30)),
        ..FactorOptions::default()
    };
    let report = factor(&n, &options).expect("factors");
    assert_eq!(report.factors.len(), 2);
    let rate = pilot_rate(&report).expect("pilot_rate timing entry");
    assert!(rate > 0.0, "expected a positive pilot rate, got {rate}");
}

#[test]
fn test_race_pilot_carries_relations_and_error_kind() {
    let n = parse_decimal("11783510925019").expect("valid 44-bit input");
    let options = FactorOptions {
        skip_preprocessing: true,
        max_duration: Some(Duration::from_secs(30)),
        ..FactorOptions::default()
    };
    let report = factor(&n, &options).expect("factors");
    let pilots: Vec<_> = report
        .timings
        .iter()
        .filter(|timing| timing.stage == "race_pilot")
        .collect();
    assert!(!pilots.is_empty(), "expected race_pilot entries");
    for pilot in pilots {
        assert!(
            pilot.detail.contains("relations="),
            "race_pilot detail lacks relations: {}",
            pilot.detail
        );
        assert!(
            pilot.detail.contains("err="),
            "race_pilot detail lacks err kind: {}",
            pilot.detail
        );
    }
}
