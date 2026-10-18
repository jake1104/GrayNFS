use std::time::Duration;

use graynfs_core::algebraic::{QuadraticElement, QuadraticRing};
use graynfs_core::error::FactorError;
use graynfs_core::factor::factor;
use graynfs_core::options::FactorOptions;
use graynfs_core::polynomial::select_polynomial;
use graynfs_core::{Integer, integer_gcd, parse_decimal};

fn product(factors: &[String]) -> Integer {
    factors
        .iter()
        .map(|value| parse_decimal(value).expect("factor decimal"))
        .fold(Integer::from(1), |accumulator, value| accumulator * value)
}

#[test]
fn factors_small_integers() {
    let options = FactorOptions::default();
    for (input, mut expected) in [
        ("15", vec!["3", "5"]),
        ("77", vec!["7", "11"]),
        ("97", vec!["97"]),
    ] {
        let number = parse_decimal(input).expect("valid input");
        let report = factor(&number, &options).expect("factorization succeeds");
        let mut actual = report.factors.clone();
        actual.sort();
        expected.sort();
        assert_eq!(actual, expected);
        assert_eq!(product(&report.factors), number);
    }
}

#[test]
fn fermat_win_reports_fermat_method() {
    let number = parse_decimal("1000036000099").expect("valid input");
    let report = factor(&number, &FactorOptions::default()).expect("factorization succeeds");
    assert_eq!(report.method, "fermat");
    let mut actual = report.factors.clone();
    actual.sort();
    assert_eq!(actual, vec!["1000003", "1000033"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn pm1_win_reports_pm1_method() {
    let number = parse_decimal("17497577401").expect("valid input");
    let options = FactorOptions {
        use_fermat: false,
        use_pp1: false,
        use_ecm: false,
        use_gnfs: false,
        allow_fallback: false,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("factorization succeeds");
    assert_eq!(report.method, "pm1");
    assert_eq!(product(&report.factors), number);
}

#[test]
fn unsupported_degree_fails_before_fallback() {
    let number = parse_decimal("100160063").expect("valid input");
    let options = FactorOptions {
        degree: Some(6),
        allow_fallback: true,
        ..FactorOptions::default()
    };
    assert!(matches!(
        factor(&number, &options),
        Err(FactorError::InvalidInput(_))
    ));
}

#[test]
fn preserves_duplicate_factors_from_squares() {
    let number = parse_decimal("289").expect("valid input");
    let report = factor(&number, &FactorOptions::default()).expect("factorization succeeds");
    assert_eq!(report.factors, vec!["17", "17"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn selected_quadratic_has_common_root_modulo_n() {
    let number = parse_decimal("1000036000099").expect("valid input");
    let options = FactorOptions {
        degree: Some(2),
        ..FactorOptions::default()
    };
    let polynomial = select_polynomial(&number, &options).expect("polynomial");
    assert_eq!(polynomial.degree, 2);
    assert_eq!(polynomial.f_at_m(), number);
    assert_eq!(
        polynomial.homogeneous(&polynomial.m, &Integer::from(1)),
        number
    );
    assert!(polynomial.constant() > Integer::from(0));
    assert!(polynomial.m > Integer::from(0));
}

#[test]
fn quadratic_ring_finds_exact_algebraic_square_root() {
    let ring = QuadraticRing::new(Integer::from(3));
    let root = QuadraticElement {
        x: Integer::from(2),
        y: Integer::from(1),
    };
    let value = ring.multiply(&root, &root);
    let computed = ring.square_root(&value).expect("exact root");
    let square = ring.multiply(&computed, &computed);
    assert_eq!(square.x, value.x);
    assert_eq!(square.y, value.y);
}

#[test]
fn quadratic_relation_collector_produces_relations() {
    let number = parse_decimal("1000036000099").expect("valid input");
    let options = FactorOptions {
        degree: Some(2),
        rational_base_bound: 2_000,
        algebraic_base_bound: 2_000,
        sieve_a_bound: 700,
        sieve_b_bound: 120,
        ..FactorOptions::default()
    };
    let polynomial = select_polynomial(&number, &options).expect("polynomial");
    let base = graynfs_core::base::FactorBase::build(
        &polynomial,
        &number,
        options.rational_base_bound,
        options.algebraic_base_bound,
    )
    .expect("factor base");
    let mut progress = |_: &graynfs_core::ProgressEvent| Ok(());
    let relations = graynfs_core::relations::collect_relations(
        &number,
        &polynomial,
        &base,
        &options,
        &std::sync::atomic::AtomicBool::new(false),
        None,
        &mut progress,
    )
    .expect("relations");
    assert!(!relations.is_empty());
    assert!(base.column_count() > 0);
}

#[test]
fn rejects_invalid_inputs() {
    let options = FactorOptions::default();
    assert!(matches!(
        factor(&Integer::from(0), &options),
        Err(FactorError::InvalidInput(_))
    ));
    assert!(matches!(
        factor(&Integer::from(-1), &options),
        Err(FactorError::InvalidInput(_))
    ));
    assert!(parse_decimal("not-a-number").is_err());
    assert!(integer_gcd(&Integer::from(2), &Integer::from(4)) == Integer::from(2));
}

#[test]
#[ignore = "requires release profile"]
fn cubic_gnfs_factors_semiprime() {
    let number = parse_decimal("100160063").expect("valid input");
    let options = FactorOptions {
        degree: Some(3),
        rational_base_bound: 500,
        algebraic_base_bound: 500,
        sieve_a_bound: 1_600,
        sieve_b_bound: 700,
        max_duration: Some(Duration::from_secs(90)),
        allow_fallback: false,
        skip_preprocessing: true,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("cubic GNFS factorization succeeds");
    assert_eq!(report.method, "gnfs");
    assert_eq!(product(&report.factors), number);
    assert!(report.factors.iter().any(|value| value == "10007"));
    assert!(report.factors.iter().any(|value| value == "10009"));
}

#[test]
#[ignore = "requires release profile"]
fn cubic_gnfs_handles_denser_factor_base() {
    let number = parse_decimal("100460333").expect("valid input");
    let options = FactorOptions {
        degree: Some(3),
        rational_base_bound: 700,
        algebraic_base_bound: 700,
        sieve_a_bound: 1_800,
        sieve_b_bound: 800,
        max_duration: Some(Duration::from_secs(90)),
        allow_fallback: false,
        skip_preprocessing: true,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("cubic GNFS factorization succeeds");
    assert_eq!(report.method, "gnfs");
    assert_eq!(product(&report.factors), number);
    assert!(report.factors.iter().any(|value| value == "10009"));
    assert!(report.factors.iter().any(|value| value == "10037"));
}

#[test]
#[ignore = "requires release profile"]
fn quartic_gnfs_factors_semiprime() {
    let number = parse_decimal("100160063").expect("valid input");
    let options = FactorOptions {
        degree: Some(4),
        rational_base_bound: 2_000,
        algebraic_base_bound: 2_000,
        sieve_a_bound: 6_000,
        sieve_b_bound: 1_500,
        max_duration: Some(Duration::from_secs(1_200)),
        allow_fallback: false,
        skip_preprocessing: true,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("quartic GNFS factorization succeeds");
    assert_eq!(report.method, "gnfs");
    assert_eq!(product(&report.factors), number);
    assert!(report.factors.iter().any(|value| value == "10007"));
    assert!(report.factors.iter().any(|value| value == "10009"));
}

#[test]
#[ignore = "wide GNFS benchmark"]
fn wide_gnfs_regression() {
    let number = parse_decimal("1000036000099").expect("valid input");
    let options = FactorOptions {
        rational_base_bound: 5_000,
        algebraic_base_bound: 5_000,
        sieve_a_bound: 2_500,
        sieve_b_bound: 300,
        max_duration: Some(Duration::from_secs(120)),
        allow_fallback: false,
        skip_preprocessing: true,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("GNFS factorization succeeds");
    assert_eq!(report.method, "siqs");
    assert_eq!(product(&report.factors), number);
    assert!(report.factors.iter().any(|value| value == "1000003"));
    assert!(report.factors.iter().any(|value| value == "1000033"));
}

#[test]
fn siqs_parallel_factors_small_semiprime() {
    let number = parse_decimal("100160063").expect("valid input");
    let options = FactorOptions {
        degree: Some(2),
        skip_preprocessing: true,
        threads: Some(2),
        max_duration: Some(Duration::from_secs(300)),
        allow_fallback: false,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("SIQS factorization succeeds");
    assert_eq!(report.method, "siqs");
    let mut actual = report.factors.clone();
    actual.sort();
    assert_eq!(actual, vec!["10007", "10009"]);
    assert_eq!(product(&report.factors), number);
    assert!(report.relations > 0, "sieve path must collect relations");
}

#[test]
fn siqs_win_reports_siqs_method() {
    let number = parse_decimal("100460333").expect("valid input");
    let options = FactorOptions {
        degree: Some(2),
        skip_preprocessing: true,
        threads: Some(2),
        max_duration: Some(Duration::from_secs(300)),
        allow_fallback: false,
        ..FactorOptions::default()
    };
    let report = factor(&number, &options).expect("SIQS factorization succeeds");
    assert_eq!(report.method, "siqs");
    let mut actual = report.factors.clone();
    actual.sort();
    assert_eq!(actual, vec!["10009", "10037"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn trial_division_win_reports_trial_division_method() {
    let number = parse_decimal("3000099").expect("valid input");
    let report = factor(&number, &FactorOptions::default()).expect("factorization succeeds");
    assert_eq!(report.method, "trial_division");
    let mut actual = report.factors.clone();
    actual.sort();
    assert_eq!(actual, vec!["1000033", "3"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn prime_input_reports_prime_method() {
    let number = parse_decimal("97").expect("valid input");
    let report = factor(&number, &FactorOptions::default()).expect("factorization succeeds");
    assert_eq!(report.method, "prime");
    assert_eq!(report.factors, vec!["97"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn perfect_square_reports_square_path() {
    let number = parse_decimal("289").expect("valid input");
    let report = factor(&number, &FactorOptions::default()).expect("factorization succeeds");
    assert!(
        report.method.contains("perfect_square"),
        "method was {}",
        report.method
    );
    assert_eq!(report.factors, vec!["17", "17"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn perfect_power_reports_power_path() {
    let number = parse_decimal("1000").expect("valid input");
    let report = factor(&number, &FactorOptions::default()).expect("factorization succeeds");
    assert!(
        report.method.contains("perfect_power"),
        "method was {}",
        report.method
    );
    let mut actual = report.factors.clone();
    actual.sort();
    assert_eq!(actual, vec!["2", "2", "2", "5", "5", "5"]);
    assert_eq!(product(&report.factors), number);
}

#[test]
fn public_parse_decimal_maps_bad_inputs_to_invalid_input() {
    let bad_inputs = [
        "",
        "-",
        "+",
        "1e3",
        "0x10",
        " 12",
        "12 ",
        "1\u{0}",
        "\u{660}\u{661}\u{662}",
        "\u{ff11}\u{ff12}",
        "--1",
        "1,000",
    ];
    for candidate in bad_inputs {
        match parse_decimal(candidate) {
            Ok(value) => panic!("unexpected ok for candidate {:?}: {value}", candidate),
            Err(FactorError::InvalidInput(_)) => {}
            Err(other) => panic!("unexpected error variant for candidate {candidate:?}: {other:?}"),
        }
    }
    assert_eq!(
        parse_decimal("+0").expect("plus zero parses"),
        Integer::from(0)
    );
    assert_eq!(
        parse_decimal("007").expect("leading zeros parse"),
        Integer::from(7)
    );
    assert_eq!(
        parse_decimal("1_000").expect("underscore separator parses"),
        Integer::from(1000)
    );
}
