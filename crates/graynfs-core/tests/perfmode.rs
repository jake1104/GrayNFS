use graynfs_core::factor::factor;
use graynfs_core::options::FactorOptions;
use graynfs_core::parse_decimal;

fn cascade_options(perf_mode: bool) -> FactorOptions {
    FactorOptions {
        use_fermat: false,
        use_pm1: false,
        use_pp1: false,
        use_ecm: false,
        use_gnfs: false,
        perf_mode,
        ..FactorOptions::default()
    }
}

fn cascade_detail(report: &graynfs_core::factor::FactorReport) -> String {
    report
        .timings
        .iter()
        .find(|entry| entry.stage == "cascade")
        .unwrap_or_else(|| panic!("expected a cascade timing entry"))
        .detail
        .clone()
}

#[test]
fn test_perf_mode_default_off() {
    assert!(!FactorOptions::default().perf_mode);
}

#[test]
fn test_perf_mode_identical_factors() {
    let number = parse_decimal("100160063").expect("valid input");
    let plain = factor(&number, &cascade_options(false)).expect("factorization succeeds");
    let perf = factor(&number, &cascade_options(true)).expect("factorization succeeds");
    assert_eq!(plain.factors, vec!["10007", "10009"]);
    assert_eq!(perf.factors, plain.factors);
    assert_eq!(perf.method, plain.method);
    assert_eq!(perf.relations, plain.relations);
    assert_eq!(perf.rank, plain.rank);
    let plain_detail = cascade_detail(&plain);
    let perf_detail = cascade_detail(&perf);
    assert!(!plain_detail.contains("perf_mode=true"));
    assert!(perf_detail.contains("perf_mode=true"));
    assert_eq!(perf_detail.replace(" perf_mode=true", ""), plain_detail);
}
