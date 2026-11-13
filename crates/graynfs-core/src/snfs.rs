#[cfg(test)]
use crate::factor::SMALL_COEFFICIENT_SCAN_MIN_BITS;
#[cfg(test)]
use crate::factor::SpecialForm;
use crate::factor::detect_special_form;
use crate::general::integer_root;
#[cfg(test)]
use crate::general::{IrreducibilityOutcome, classify_irreducibility, has_unit_index};
use crate::integer::{Integer, IntegerExt};
use crate::number_field::NumberField;
use crate::polynomial::PolynomialPair;
use num_traits::Zero;

pub struct SnfsForm {
    pub a: u64,
    pub b: u64,
    pub k: u32,
    pub c: u64,
    pub sign: i8,
}

impl SnfsForm {
    pub fn new(a: u64, b: u64, k: u32, c: u64, sign: i8) -> Self {
        Self { a, b, k, c, sign }
    }
}

pub fn is_snfs_eligible(form: &SnfsForm) -> bool {
    if form.c != 1 {
        return false;
    }
    if !(3..=5).contains(&form.k) {
        return false;
    }
    if !(1..=16).contains(&form.a) {
        return false;
    }
    if !(2..=1_000_000).contains(&form.b) {
        return false;
    }
    form.sign == 1 || form.sign == -1
}

pub struct SnfsPolynomial {
    pub field: NumberField,
    pub root: Integer,
}

pub fn snfs_integer(form: &SnfsForm) -> Integer {
    let magnitude = Integer::from(form.a) * Integer::from(form.b).pow(form.k);
    magnitude + Integer::from(form.sign)
}

pub fn snfs_polynomial(form: &SnfsForm) -> Option<SnfsPolynomial> {
    if !is_snfs_eligible(form) {
        return None;
    }
    let degree = form.k as usize;
    let scale = Integer::from(form.a).pow(form.k - 1);
    let constant = if form.sign == 1 { scale } else { -scale };
    let mut coefficients = vec![Integer::from(0); degree + 1];
    coefficients[0] = constant;
    coefficients[degree] = Integer::from(1);
    let field = NumberField::new(degree, coefficients)?;
    let root = Integer::from(form.a) * Integer::from(form.b);
    Some(SnfsPolynomial { field, root })
}

pub(crate) struct SnfsPreset {
    pub polynomial: PolynomialPair,
    pub detail: String,
}

pub(crate) fn snfs_binomial_irreducible(degree: usize, constant: &Integer) -> bool {
    let shifted = -constant.clone();
    let magnitude = crate::integer::integer_abs(&shifted);
    for prime in [2_u32, 3, 5] {
        if !(degree as u32).is_multiple_of(prime) {
            continue;
        }
        if prime % 2 == 0 && shifted < Integer::from(0) {
            continue;
        }
        let root = integer_root(&magnitude, prime as usize);
        if root.pow(prime) == magnitude {
            return false;
        }
    }
    if degree.is_multiple_of(4) && shifted < Integer::from(0) {
        let four = Integer::from(4);
        if (&magnitude % &four).is_zero() {
            let quarter = magnitude / four;
            let root = integer_root(&quarter, 4);
            if root.pow(4) == quarter {
                return false;
            }
        }
    }
    true
}

pub(crate) fn snfs_preset_for(n: &Integer, degree: usize) -> Option<SnfsPreset> {
    if !(3..=5).contains(&degree) {
        return None;
    }
    let detected = detect_special_form(n);
    let form = match detected {
        Some(form) if form.k != 0 && form.k as usize == degree && form.c == 1 => {
            SnfsForm::new(form.a, form.b, form.k, form.c, form.sign)
        }
        _ => recover_exact_form(n, degree)?,
    };
    if !is_snfs_eligible(&form) {
        return None;
    }
    let built = snfs_polynomial(&form)?;
    let scale = Integer::from(form.a).pow(form.k - 1);
    let polynomial = PolynomialPair {
        degree,
        coefficients: built.field.polynomial.clone(),
        m: built.root.clone(),
    };
    if polynomial.evaluate(&polynomial.m) != &scale * n {
        return None;
    }
    if !snfs_binomial_irreducible(degree, &polynomial.coefficients[0]) {
        return None;
    }
    Some(SnfsPreset {
        polynomial,
        detail: format!(
            "snfs a={} b={} k={} c={} sign={}",
            form.a, form.b, form.k, form.c, form.sign
        ),
    })
}

fn recover_exact_form(n: &Integer, degree: usize) -> Option<SnfsForm> {
    let one = Integer::from(1);
    for coefficient in 1..=16_u64 {
        let magnitude = Integer::from(coefficient);
        for sign in [1_i8, -1] {
            let shifted = if sign == 1 { n - &one } else { n + &one };
            if shifted <= Integer::from(0)
                || shifted.clone() % magnitude.clone() != Integer::from(0)
            {
                continue;
            }
            let quotient = shifted / magnitude.clone();
            let root = integer_root(&quotient, degree);
            if root.pow(degree as u32) != quotient {
                continue;
            }
            let base = root.to_u64_saturating();
            if !(2..=1_000_000).contains(&base) {
                continue;
            }
            return Some(SnfsForm::new(coefficient, base, degree as u32, 1, sign));
        }
    }
    None
}

#[cfg(test)]
pub(crate) fn snfs_e2e_corpus() -> Vec<SnfsCorpusEntry> {
    vec![
        SnfsCorpusEntry::new(2, 10, 3, 1),
        SnfsCorpusEntry::new(2, 12, 3, -1),
        SnfsCorpusEntry::new(1, 10, 4, 1),
        SnfsCorpusEntry::new(3, 6, 4, -1),
        SnfsCorpusEntry::new(2, 6, 5, 1),
        SnfsCorpusEntry::new(2, 4, 5, -1),
    ]
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SnfsDecision {
    Promote,
    Fallback(&'static str),
}

#[cfg(test)]
pub(crate) fn gate_snfs_form(form: &SnfsForm) -> SnfsDecision {
    if !is_snfs_eligible(form) {
        return SnfsDecision::Fallback("shape");
    }
    let built = match snfs_polynomial(form) {
        Some(built) => built,
        None => return SnfsDecision::Fallback("shape"),
    };
    let candidate = PolynomialPair {
        degree: built.field.degree,
        coefficients: built.field.polynomial.clone(),
        m: built.root.clone(),
    };
    if matches!(
        classify_irreducibility(&candidate),
        IrreducibilityOutcome::Rejected
    ) {
        return SnfsDecision::Fallback("reducible");
    }
    if !has_unit_index(&candidate) {
        return SnfsDecision::Fallback("unit-index");
    }
    SnfsDecision::Promote
}

#[cfg(test)]
pub(crate) fn gate_special_form(n: &Integer, form: &SpecialForm) -> SnfsDecision {
    if form.k == 0 {
        return SnfsDecision::Fallback("low-weight-only");
    }
    if form.a != 1 && n.bits() < SMALL_COEFFICIENT_SCAN_MIN_BITS {
        return SnfsDecision::Fallback("sub-128-small-coefficient");
    }
    gate_snfs_form(&SnfsForm::new(form.a, form.b, form.k, form.c, form.sign))
}

#[cfg(test)]
pub(crate) fn gate_detected(n: &Integer) -> SnfsDecision {
    match detect_special_form(n) {
        Some(form) => gate_special_form(n, &form),
        None => SnfsDecision::Fallback("undetected"),
    }
}

#[cfg(test)]
pub(crate) fn decide_on_scan(scan: Option<Integer>) -> SnfsDecision {
    match scan {
        Some(_) => SnfsDecision::Promote,
        None => SnfsDecision::Fallback("scan-miss"),
    }
}

pub struct SnfsCorpusEntry {
    pub form: SnfsForm,
    pub n: Integer,
}

impl SnfsCorpusEntry {
    pub fn new(a: u64, b: u64, k: u32, sign: i8) -> Self {
        let form = SnfsForm::new(a, b, k, 1, sign);
        let n = snfs_integer(&form);
        Self { form, n }
    }
}

pub fn snfs_corpus() -> Vec<SnfsCorpusEntry> {
    vec![
        SnfsCorpusEntry::new(1, 10, 3, 1),
        SnfsCorpusEntry::new(2, 10, 3, 1),
        SnfsCorpusEntry::new(2, 10, 3, -1),
        SnfsCorpusEntry::new(1, 10, 4, 1),
        SnfsCorpusEntry::new(3, 6, 4, -1),
        SnfsCorpusEntry::new(1, 10, 5, -1),
        SnfsCorpusEntry::new(2, 6, 5, 1),
        SnfsCorpusEntry::new(16, 2, 5, 1),
    ]
}

#[cfg(test)]
mod tests {
    use super::{Integer, SnfsForm};
    use super::{
        SnfsDecision, decide_on_scan, gate_detected, gate_snfs_form, gate_special_form,
        is_snfs_eligible, snfs_binomial_irreducible, snfs_corpus, snfs_e2e_corpus, snfs_polynomial,
        snfs_preset_for,
    };
    use crate::factor::SpecialForm;
    use crate::general::has_unit_index;
    use crate::polynomial::PolynomialPair;
    use num_traits::Zero;

    fn evaluate(coefficients: &[Integer], point: &Integer) -> Integer {
        let mut result = Integer::from(0);
        for coefficient in coefficients.iter().rev() {
            result = result * point + coefficient;
        }
        result
    }

    #[test]
    fn corpus_covers_required_shapes() {
        let corpus = snfs_corpus();
        assert!(corpus.len() >= 6);
        for degree in [3_u32, 4, 5] {
            assert!(
                corpus.iter().any(|entry| entry.form.k == degree),
                "missing degree {degree}"
            );
        }
        assert!(corpus.iter().any(|entry| entry.form.sign == 1));
        assert!(corpus.iter().any(|entry| entry.form.sign == -1));
        assert!(corpus.iter().any(|entry| entry.form.a == 1));
        assert!(corpus.iter().any(|entry| entry.form.a > 1));
        for entry in &corpus {
            assert!(is_snfs_eligible(&entry.form));
        }
    }

    #[test]
    fn corpus_polynomials_hit_their_integers() {
        for entry in snfs_corpus() {
            let built = snfs_polynomial(&entry.form).expect("eligible form builds");
            assert_eq!(built.field.degree, entry.form.k as usize);
            assert_eq!(built.field.polynomial.len(), entry.form.k as usize + 1);
            assert_eq!(
                built.field.polynomial[entry.form.k as usize],
                Integer::from(1)
            );
            for index in 1..entry.form.k as usize {
                assert!(built.field.polynomial[index].is_zero());
            }
            let scale = Integer::from(entry.form.a).pow(entry.form.k - 1);
            let expected_constant = if entry.form.sign == 1 {
                scale.clone()
            } else {
                -scale.clone()
            };
            assert_eq!(built.field.polynomial[0], expected_constant);
            assert_eq!(
                built.root,
                Integer::from(entry.form.a) * Integer::from(entry.form.b)
            );
            let value = evaluate(&built.field.polynomial, &built.root);
            assert_eq!(value, &scale * &entry.n);
            assert!((value.clone() % &entry.n).is_zero());
            if entry.form.a == 1 {
                assert_eq!(value, entry.n);
            }
        }
    }

    #[test]
    fn eligibility_boundary_values_hold() {
        for sign in [1_i8, -1] {
            assert!(is_snfs_eligible(&SnfsForm::new(1, 2, 3, 1, sign)));
            assert!(is_snfs_eligible(&SnfsForm::new(16, 1_000_000, 5, 1, sign)));
        }
    }

    #[test]
    fn corpus_presets_route_through_irreducibility_and_unit_gates() {
        let expected = [
            "reducible",
            "unit-index",
            "unit-index",
            "reducible",
            "unit-index",
            "reducible",
            "unit-index",
            "unit-index",
        ];
        for (entry, reason) in snfs_corpus().iter().zip(expected) {
            assert_eq!(gate_snfs_form(&entry.form), SnfsDecision::Fallback(reason));
        }
    }

    #[test]
    fn low_weight_only_detections_never_promote() {
        let heleven = Integer::from(1_000_000_192_513_u64);
        let detected = crate::factor::detect_special_form(&heleven).expect("heuristic detection");
        assert_eq!(detected.k, 0);
        assert_eq!(
            gate_detected(&heleven),
            SnfsDecision::Fallback("low-weight-only")
        );
        let synthetic = SpecialForm {
            a: 0,
            b: 0,
            k: 0,
            c: 0,
            sign: 0,
            low_weight: true,
        };
        assert_eq!(
            gate_special_form(&heleven, &synthetic),
            SnfsDecision::Fallback("low-weight-only")
        );
    }

    #[test]
    fn sub_128_small_coefficient_forms_never_promote() {
        let n = Integer::from(2_000_000_001_u64);
        assert!(n.bits() < 128);
        assert!(crate::factor::detect_special_form(&n).is_none());
        assert_eq!(gate_detected(&n), SnfsDecision::Fallback("undetected"));
        let synthetic = SpecialForm {
            a: 2,
            b: 1000,
            k: 3,
            c: 1,
            sign: 1,
            low_weight: false,
        };
        assert_eq!(
            gate_special_form(&n, &synthetic),
            SnfsDecision::Fallback("sub-128-small-coefficient")
        );
    }

    #[test]
    fn above_128_small_coefficient_clears_the_bit_gate() {
        let base = Integer::from(1_000_000_u64);
        let n = Integer::from(2) * base.pow(7_u32) + Integer::from(1);
        assert!(n.bits() >= 128);
        let detected =
            crate::factor::detect_special_form(&n).expect("exact small-coefficient detection");
        assert_eq!(detected.a, 2);
        assert_eq!(
            gate_special_form(&n, &detected),
            SnfsDecision::Fallback("shape")
        );
    }

    #[test]
    fn scan_miss_falls_back_to_generic() {
        assert_eq!(decide_on_scan(None), SnfsDecision::Fallback("scan-miss"));
        assert_eq!(
            decide_on_scan(Some(Integer::from(7))),
            SnfsDecision::Promote
        );
    }

    #[test]
    fn ineligible_shapes_stay_out() {
        let cases = [
            SnfsForm::new(1, 10, 0, 1, 1),
            SnfsForm::new(1, 10, 2, 1, 1),
            SnfsForm::new(1, 10, 6, 1, 1),
            SnfsForm::new(1, 10, 3, 2, 1),
            SnfsForm::new(0, 10, 3, 1, 1),
            SnfsForm::new(17, 10, 3, 1, 1),
            SnfsForm::new(1, 1, 3, 1, 1),
            SnfsForm::new(1, 1_000_001, 3, 1, 1),
            SnfsForm::new(1, 10, 3, 1, 0),
        ];
        for form in &cases {
            assert!(!is_snfs_eligible(form));
            assert!(snfs_polynomial(form).is_none());
        }
    }

    #[test]
    fn shape_ineligible_forms_fall_back_to_generic() {
        for form in [
            SnfsForm::new(1, 10, 2, 1, 1),
            SnfsForm::new(1, 10, 6, 1, 1),
            SnfsForm::new(0, 10, 3, 1, 1),
        ] {
            assert_eq!(gate_snfs_form(&form), SnfsDecision::Fallback("shape"));
        }
    }

    #[test]
    fn scoped_capelli_accepts_e2e_binomials() {
        for entry in snfs_e2e_corpus() {
            let built = snfs_polynomial(&entry.form).expect("eligible form builds");
            assert!(
                snfs_binomial_irreducible(entry.form.k as usize, &built.field.polynomial[0]),
                "a={} b={} k={} sign={}",
                entry.form.a,
                entry.form.b,
                entry.form.k,
                entry.form.sign
            );
        }
    }

    #[test]
    fn scoped_capelli_rejects_genuinely_reducible_binomials() {
        for form in [
            SnfsForm::new(1, 10, 3, 1, 1),
            SnfsForm::new(1, 10, 5, 1, -1),
        ] {
            let built = snfs_polynomial(&form).expect("eligible form builds");
            assert!(!snfs_binomial_irreducible(
                form.k as usize,
                &built.field.polynomial[0]
            ));
        }
        let quartic = SnfsForm::new(1, 10, 4, 1, 1);
        let built = snfs_polynomial(&quartic).expect("eligible form builds");
        assert!(snfs_binomial_irreducible(4, &built.field.polynomial[0]));
    }

    #[test]
    fn scoped_preset_promotes_all_six_despite_unit_index() {
        for entry in snfs_e2e_corpus() {
            let preset =
                snfs_preset_for(&entry.n, entry.form.k as usize).expect("scoped preset promotes");
            assert!(preset.detail.contains("snfs a="));
            let candidate = PolynomialPair {
                degree: entry.form.k as usize,
                coefficients: preset.polynomial.coefficients.clone(),
                m: preset.polynomial.m.clone(),
            };
            assert!(!has_unit_index(&candidate));
        }
    }

    #[test]
    fn large_corpus_entries_pass_production_gates() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("bench")
            .join("snfs")
            .join("corpus-large.json");
        let text = std::fs::read_to_string(&path).expect("bench snfs corpus-large.json readable");
        let document: serde_json::Value =
            serde_json::from_str(&text).expect("corpus-large.json parses");
        let entries = document["entries"].as_array().expect("entries array");
        assert!(entries.len() >= 6);
        let mut degrees = std::collections::HashSet::new();
        let mut signs = std::collections::HashSet::new();
        let mut unit_coefficient = false;
        let mut multi_coefficient = false;
        for entry in entries {
            let a = entry["a"].as_u64().expect("a");
            let b = entry["b"].as_u64().expect("b");
            let k = entry["k"].as_u64().expect("k") as u32;
            let sign = entry["sign"].as_i64().expect("sign") as i8;
            let bits = entry["bits"].as_u64().expect("bits");
            assert!((60..=100).contains(&bits));
            let n = crate::factor::parse_decimal(entry["n"].as_str().expect("n")).expect("valid n");
            assert_eq!(n.bits(), bits);
            let form = SnfsForm::new(a, b, k, 1, sign);
            assert!(is_snfs_eligible(&form));
            let built = snfs_polynomial(&form).expect("eligible form builds");
            assert_eq!(built.field.degree, k as usize);
            assert_eq!(
                built.root,
                Integer::from(a) * Integer::from(b),
                "root m == a*b"
            );
            let scale = Integer::from(a).pow(k - 1);
            assert_eq!(built.field.polynomial[k as usize], Integer::from(1));
            for index in 1..k as usize {
                assert!(built.field.polynomial[index].is_zero());
            }
            let expected_constant = if sign == 1 {
                scale.clone()
            } else {
                -scale.clone()
            };
            assert_eq!(built.field.polynomial[0], expected_constant);
            assert!(snfs_binomial_irreducible(
                k as usize,
                &built.field.polynomial[0]
            ));
            let candidate = PolynomialPair {
                degree: k as usize,
                coefficients: built.field.polynomial.clone(),
                m: built.root.clone(),
            };
            let value = evaluate(&candidate.coefficients, &candidate.m);
            assert_eq!(value, &scale * &n, "scaled identity f(m) == a^(k-1)*N");
            assert!((value % &n).is_zero(), "f(m) mod N == 0");
            let preset = snfs_preset_for(&n, k as usize).expect("production preset promotes");
            assert!(crate::general::snfs_binomial_shape(&preset.polynomial));
            for (key, expected) in [
                ("eligible", true),
                ("scaled_identity", true),
                ("mod_zero", true),
                ("binomial_irreducible", true),
                ("preset_promotes", true),
                ("binomial_shape", true),
            ] {
                assert_eq!(entry["gates"][key], expected, "recorded gate {key}");
            }
            let mut product = Integer::from(1);
            for factor in entry["factors"].as_array().expect("factors array") {
                let prime = crate::factor::parse_decimal(factor.as_str().expect("factor"))
                    .expect("valid factor");
                product *= &prime;
            }
            assert_eq!(product, n, "recorded factors multiply to N");
            degrees.insert(k);
            signs.insert(sign);
            unit_coefficient |= a == 1;
            multi_coefficient |= a > 1;
        }
        for degree in [3_u32, 4, 5] {
            assert!(degrees.contains(&degree), "missing degree {degree}");
        }
        assert!(signs.contains(&1_i8) && signs.contains(&-1_i8));
        assert!(unit_coefficient && multi_coefficient);
    }

    #[test]
    fn e2e_snfs_true_path_and_false_parity() {
        let expected = [
            vec!["3", "23", "29"],
            vec!["5", "691"],
            vec!["73", "137"],
            vec!["13", "13", "23"],
            vec!["103", "151"],
            vec!["23", "89"],
        ];
        for (entry, factors) in snfs_e2e_corpus().iter().zip(expected) {
            let degree = entry.form.k as usize;
            let enabled = crate::options::FactorOptions {
                degree: Some(degree),
                skip_preprocessing: true,
                threads: Some(1),
                ..crate::options::FactorOptions::default()
            };
            let (enabled_result, enabled_timings) =
                crate::factor::factor_with_timings(&entry.n, &enabled, |_| Ok(()));
            let enabled_report = enabled_result.expect("snfs run factors");
            assert!(
                enabled_report.method.starts_with("gnfs"),
                "method stays gnfs for {}, got {}",
                entry.n,
                enabled_report.method
            );
            assert_eq!(enabled_report.factors, factors, "factors for {}", entry.n);
            assert!(
                enabled_timings
                    .iter()
                    .any(|entry| entry.stage == "select_polynomial"
                        && entry.detail.contains("snfs a=")),
                "snfs preset must drive collection for {}: {}",
                entry.n,
                enabled_timings
                    .iter()
                    .map(|entry| format!("{}:{}", entry.stage, entry.detail))
                    .collect::<Vec<_>>()
                    .join(" | ")
            );
            let disabled = crate::options::FactorOptions {
                use_snfs: false,
                ..enabled.clone()
            };
            let (disabled_result, disabled_timings) =
                crate::factor::factor_with_timings(&entry.n, &disabled, |_| Ok(()));
            let disabled_report = disabled_result.expect("generic run factors");
            assert!(
                disabled_report.method.starts_with("gnfs"),
                "method stays gnfs, got {}",
                disabled_report.method
            );
            assert_eq!(disabled_report.factors, factors);
            assert!(
                !disabled_timings
                    .iter()
                    .any(|entry| entry.detail.contains("snfs a=")),
                "disabled run must stay on the generic path"
            );
        }
    }

    #[test]
    fn e2e_snfs_hostile_yield_shape_falls_back_soundly() {
        let n = crate::factor::parse_decimal("513").expect("valid");
        assert!(snfs_preset_for(&n, 5).is_some());
        let enabled = crate::options::FactorOptions {
            degree: Some(5),
            skip_preprocessing: true,
            threads: Some(1),
            ..crate::options::FactorOptions::default()
        };
        let (result, timings) = crate::factor::factor_with_timings(&n, &enabled, |_| Ok(()));
        let report = result.expect("fallback factors");
        assert_eq!(report.factors, vec!["3", "3", "3", "19"]);
        assert!(report.method.starts_with("gnfs"));
        assert!(
            !timings.iter().any(|entry| entry.detail.contains("snfs a=")),
            "hostile shape must fall back to the generic path"
        );
    }

    #[test]
    fn snfs_binomial_shape_gate_matches_only_binomials() {
        for entry in snfs_e2e_corpus() {
            let preset = snfs_preset_for(&entry.n, entry.form.k as usize).expect("preset");
            assert!(crate::general::snfs_binomial_shape(&preset.polynomial));
        }
        let generic = PolynomialPair {
            degree: 3,
            coefficients: vec![
                Integer::from(7),
                Integer::from(2),
                Integer::from(5),
                Integer::from(1),
            ],
            m: Integer::from(11),
        };
        assert!(!crate::general::snfs_binomial_shape(&generic));
    }

    #[test]
    fn e2e_corpus_matches_bench_inputs() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("bench")
            .join("snfs")
            .join("corpus.json");
        let text = std::fs::read_to_string(&path).expect("bench snfs corpus.json readable");
        assert_eq!(snfs_e2e_corpus().len(), 6);
        for entry in snfs_e2e_corpus() {
            assert!(text.contains(&format!("\"n\": \"{}\"", entry.n)));
            assert!(text.contains(&format!("\"degree\": {}", entry.form.k)));
        }
    }

    #[test]
    fn standing_pins_route_around_the_snfs_path() {
        for decimal in [
            "770817968652704265841785071059",
            "35658606617546370626847644936101",
            "1155285037654380637",
        ] {
            let n = crate::factor::parse_decimal(decimal).expect("valid pin");
            for degree in [3_usize, 4, 5] {
                assert!(snfs_preset_for(&n, degree).is_none());
            }
        }
    }
}
