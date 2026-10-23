use graynfs_core::prime::{
    AKS_PROOF_MAX_BITS, is_baillie_psw_prime, is_frobenius_probable_prime, is_prime_aks_proof,
    is_probable_prime, is_probable_prime_with_options, passes_miller_rabin_base2,
};
use graynfs_core::{FactorOptions, Integer, factor, parse_decimal};

fn carmichael_corpus() -> Vec<(&'static str, Vec<&'static str>)> {
    vec![
        ("561", vec!["3", "11", "17"]),
        ("1105", vec!["5", "13", "17"]),
        ("1729", vec!["7", "13", "19"]),
        ("2465", vec!["5", "17", "29"]),
        ("2821", vec!["7", "13", "31"]),
        ("6601", vec!["7", "23", "41"]),
        ("8911", vec!["7", "19", "67"]),
        ("10585", vec!["5", "29", "73"]),
        ("15841", vec!["7", "31", "73"]),
        ("29341", vec!["13", "37", "61"]),
        ("41041", vec!["7", "11", "13", "41"]),
        ("46657", vec!["13", "37", "97"]),
        ("52633", vec!["7", "73", "103"]),
        ("62745", vec!["3", "5", "47", "89"]),
        ("63973", vec!["7", "13", "19", "37"]),
        ("75361", vec!["11", "13", "17", "31"]),
    ]
}

fn pseudoprime_corpus() -> Vec<(&'static str, Vec<&'static str>)> {
    vec![
        ("2047", vec!["23", "89"]),
        ("3215031751", vec!["151", "751", "28351"]),
        ("2152302898747", vec!["6763", "10627", "29947"]),
        ("3474749660383", vec!["1303", "16927", "157543"]),
        ("341550071728321", vec!["10670053", "32010157"]),
        ("3825123056546413051", vec!["149491", "747451", "34233211"]),
    ]
}

fn prime_corpus() -> Vec<&'static str> {
    vec![
        "2",
        "3",
        "5",
        "7",
        "11",
        "13",
        "17",
        "19",
        "23",
        "29",
        "31",
        "37",
        "41",
        "43",
        "47",
        "101",
        "997",
        "1000003",
        "2147483647",
        "170141183460469231731687303715884105727",
    ]
}

fn frobenius_on() -> FactorOptions {
    FactorOptions {
        use_frobenius: true,
        ..FactorOptions::default()
    }
}

fn aks_feasible(value: &Integer) -> bool {
    value.bits() <= 43
}

fn aks_prime_feasible(value: &Integer) -> bool {
    value.bits() <= 17
}

fn product_of(decimals: &[String]) -> Integer {
    decimals
        .iter()
        .map(|entry| parse_decimal(entry).expect("factor decimal parses"))
        .fold(Integer::from(1), |accumulator, factor| accumulator * factor)
}

fn assert_composite_verdict(value: &Integer) {
    assert!(!is_probable_prime(value));
    assert!(!is_baillie_psw_prime(value));
    assert!(!is_probable_prime_with_options(
        value,
        &FactorOptions::default()
    ));
    assert!(!is_probable_prime_with_options(value, &frobenius_on()));
    assert!(!is_frobenius_probable_prime(value));
    if aks_feasible(value) {
        assert!(!is_prime_aks_proof(value));
    }
}

fn assert_prime_verdict(value: &Integer) {
    assert!(is_probable_prime(value));
    assert!(is_baillie_psw_prime(value));
    assert!(is_probable_prime_with_options(
        value,
        &FactorOptions::default()
    ));
    assert!(is_probable_prime_with_options(value, &frobenius_on()));
    assert!(is_frobenius_probable_prime(value));
    if aks_prime_feasible(value) {
        assert!(is_prime_aks_proof(value));
    }
}

#[test]
fn carmichael_gates_declare_composite() {
    assert_eq!(carmichael_corpus().len(), 16);
    for (digits, expected) in carmichael_corpus() {
        let value = parse_decimal(digits).expect("carmichael parses");
        let check: Integer = expected
            .iter()
            .map(|factor| parse_decimal(factor).expect("factor parses"))
            .fold(Integer::from(1), |accumulator, factor| accumulator * factor);
        assert_eq!(check, value);
        assert_composite_verdict(&value);
        let report = factor(&value, &FactorOptions::default()).expect("factor succeeds");
        assert_eq!(product_of(&report.factors), value);
        let mut actual = report.factors.clone();
        actual.sort();
        let mut wanted: Vec<String> = expected.iter().map(|entry| entry.to_string()).collect();
        wanted.sort();
        assert_eq!(actual, wanted);
    }
}

#[test]
fn pseudoprime_gates_declare_composite() {
    for (digits, expected) in pseudoprime_corpus() {
        let value = parse_decimal(digits).expect("pseudoprime parses");
        let check: Integer = expected
            .iter()
            .map(|factor| parse_decimal(factor).expect("factor parses"))
            .fold(Integer::from(1), |accumulator, factor| accumulator * factor);
        assert_eq!(check, value);
        assert_composite_verdict(&value);
    }
    let value = parse_decimal("2047").expect("2047 parses");
    let report = factor(&value, &FactorOptions::default()).expect("factor succeeds");
    assert_eq!(product_of(&report.factors), value);
}

#[test]
fn prime_gates_declare_prime() {
    for digits in prime_corpus() {
        let value = parse_decimal(digits).expect("prime parses");
        if digits == "170141183460469231731687303715884105727" {
            let wide = (Integer::from(1) << 127) - Integer::from(1);
            assert_eq!(value, wide);
            assert!(value.bits() > AKS_PROOF_MAX_BITS);
            assert!(!is_prime_aks_proof(&value));
            assert!(is_probable_prime(&value));
            assert!(is_baillie_psw_prime(&value));
            assert!(is_probable_prime_with_options(
                &value,
                &FactorOptions::default()
            ));
            assert!(is_probable_prime_with_options(&value, &frobenius_on()));
            assert!(is_frobenius_probable_prime(&value));
            continue;
        }
        assert_prime_verdict(&value);
    }
    let small = parse_decimal("101").expect("101 parses");
    let report = factor(&small, &FactorOptions::default()).expect("factor succeeds");
    assert_eq!(product_of(&report.factors), small);
    assert_eq!(report.factors, vec!["101".to_string()]);
}

#[test]
fn mr_only_is_caught_by_suite() {
    for digits in ["15841", "29341", "52633"] {
        let value = parse_decimal(digits).expect("carmichael parses");
        assert!(passes_miller_rabin_base2(&value));
        assert!(!is_probable_prime(&value));
        assert!(!is_baillie_psw_prime(&value));
        assert!(!is_probable_prime_with_options(&value, &frobenius_on()));
    }
    for digits in ["2047", "3215031751", "2152302898747", "3474749660383"] {
        let value = parse_decimal(digits).expect("pseudoprime parses");
        assert!(passes_miller_rabin_base2(&value));
        assert!(!is_probable_prime(&value));
        assert!(!is_baillie_psw_prime(&value));
        assert!(!is_probable_prime_with_options(&value, &frobenius_on()));
    }
}
