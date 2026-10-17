use graynfs_core::number_field::{Element, NumberField};
use graynfs_core::poly::{ModPolynomial, factor_squarefree};
use graynfs_core::{FactorOptions, Integer, factor, parse_decimal};

#[test]
fn modular_polynomial_finds_all_linear_roots() {
    let polynomial = ModPolynomial::new(vec![1, 0, 1], 5);
    let roots = polynomial.roots();
    assert_eq!(roots, vec![2, 3]);
}

fn brute_force_roots(polynomial: &ModPolynomial) -> Vec<u64> {
    (0..polynomial.modulus)
        .filter(|point| polynomial.eval(*point) == 0)
        .collect()
}

#[test]
fn modular_polynomial_roots_match_brute_force() {
    let mut state = 0x2024_0519_7f4a_7c15_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for prime in [
        7_u64, 11, 13, 17, 19, 23, 29, 31, 101, 257, 1009, 2003, 5003,
    ] {
        for degree in 1_usize..=5 {
            for _ in 0..6 {
                let mut coefficients: Vec<u64> = (0..=degree).map(|_| next() % prime).collect();
                if coefficients[degree] == 0 {
                    coefficients[degree] = 1;
                }
                let polynomial = ModPolynomial::new(coefficients, prime);
                assert_eq!(
                    polynomial.roots(),
                    brute_force_roots(&polynomial),
                    "prime={prime} coefficients={:?}",
                    polynomial.coefficients
                );
            }
        }
    }
}

#[test]
fn modular_polynomial_roots_handle_repeated_factors() {
    let polynomial = ModPolynomial::new(vec![0, 0, 0, 1], 5);
    assert_eq!(polynomial.roots(), brute_force_roots(&polynomial));
    let repeated = ModPolynomial::new(vec![0, 0, 4, 1], 5);
    assert_eq!(repeated.roots(), brute_force_roots(&repeated));
    let constant_zero = ModPolynomial::new(vec![0, 0], 7);
    assert!(constant_zero.roots().is_empty());
    let constant = ModPolynomial::new(vec![3], 7);
    assert!(constant.roots().is_empty());
}

#[test]
fn squarefree_factorization_reconstructs_the_input() {
    let mut state = 0x51ed_270b_9f21_c3d5_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for prime in [7_u64, 11, 13, 17, 19, 23, 31, 101, 257, 1009] {
        for degree in 1_usize..=5 {
            for _ in 0..6 {
                let mut coefficients: Vec<u64> = (0..=degree).map(|_| next() % prime).collect();
                if coefficients[degree] == 0 {
                    coefficients[degree] = 1;
                }
                let polynomial = ModPolynomial::new(coefficients, prime);
                if !polynomial.derivative().gcd(&polynomial).is_constant() {
                    assert!(factor_squarefree(&polynomial).is_empty());
                    continue;
                }
                let factors = factor_squarefree(&polynomial);
                assert!(
                    factors.iter().all(|factor| factor.degree() > 0),
                    "prime={prime} coefficients={:?} factors={:?}",
                    polynomial.coefficients,
                    factors
                );
                let mut product = ModPolynomial::one(prime);
                for factor in &factors {
                    product = product.mul(factor).rem(&polynomial);
                }
                assert!(
                    product.is_zero(),
                    "prime={prime} coefficients={:?} factors={:?} product={:?}",
                    polynomial.coefficients,
                    factors,
                    product.coefficients
                );
                let total: usize = factors.iter().map(ModPolynomial::degree).sum();
                assert_eq!(total, polynomial.degree());
            }
        }
    }
}

#[test]
fn modular_polynomial_gcd_and_multiplication_are_consistent() {
    let left = ModPolynomial::new(vec![1, 1], 101);
    let right = ModPolynomial::new(vec![1, 1], 101);
    let product = left.mul(&right);
    assert_eq!(product.coefficients, vec![1, 2, 1]);
    let gcd = left.gcd(&right);
    assert_eq!(gcd.coefficients, vec![1, 1]);
}

#[test]
fn number_field_reduces_using_monic_polynomial() {
    let field = NumberField::new(
        3,
        vec![
            Integer::from(2),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ],
    )
    .expect("cubic field");
    let theta =
        Element::from_coefficients(vec![Integer::from(0), Integer::from(1), Integer::from(0)]);
    let cube = field.multiply(&theta, &field.multiply(&theta, &theta));
    assert_eq!(cube.coefficients[0], Integer::from(-2));
    assert!(cube.coefficients[1].is_zero());
    assert!(cube.coefficients[2].is_zero());
}

#[test]
fn number_field_homogeneous_matches_evaluation() {
    let number = parse_decimal("1000036000099").expect("valid number");
    let m = parse_decimal("10000").expect("valid m");
    let field = NumberField::new(
        3,
        vec![
            Integer::from(36000099),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ],
    )
    .expect("cubic field");
    assert_eq!(field.homogeneous(&m, &Integer::from(1)), number);
}

#[test]
fn lift_budget_keeps_soundness() {
    graynfs_core::hensel::clear_solver_abort();
    let n = parse_decimal("1155285037654380637").expect("valid number");
    let options = FactorOptions {
        degree: Some(3),
        skip_preprocessing: true,
        threads: Some(1),
        ..FactorOptions::default()
    };
    let started = std::time::Instant::now();
    let report = factor(&n, &options).expect("deg3-skip resolves");
    let wall_ms = started.elapsed().as_millis();
    assert_eq!(
        report.factors,
        vec!["1074741869".to_string(), "1074941873".to_string(),]
    );
    let first = parse_decimal("1074741869").expect("valid factor");
    let second = parse_decimal("1074941873").expect("valid factor");
    let product: Integer = first * second;
    assert_eq!(product, n);
    let bound_ms = if cfg!(debug_assertions) {
        120_000
    } else {
        32_500
    };
    assert!(
        wall_ms < bound_ms,
        "dep0 tail must be capped, wall={wall_ms}ms"
    );
}

use num_traits::Zero;
