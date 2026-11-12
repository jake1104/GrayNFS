#[cfg(test)]
use crate::factor::{SMALL_COEFFICIENT_SCAN_MIN_BITS, SpecialForm, detect_special_form};
#[cfg(test)]
use crate::general::{IrreducibilityOutcome, classify_irreducibility, has_unit_index};
use crate::integer::Integer;
use crate::number_field::NumberField;
#[cfg(test)]
use crate::polynomial::PolynomialPair;

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
        is_snfs_eligible, snfs_corpus, snfs_polynomial,
    };
    use crate::factor::SpecialForm;
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
}
