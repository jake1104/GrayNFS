use crate::integer::Integer;
use crate::number_field::NumberField;

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
    use super::{is_snfs_eligible, snfs_corpus, snfs_polynomial};
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
}
