use num_traits::{Euclid, ToPrimitive, Zero};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use crate::field::{ExtensionField, is_square, square_root};
use crate::integer::{Integer, integer_abs};
use crate::number_field::{Element, NumberField};
use crate::poly::{ModPolynomial, crt_reconstruct, factor_squarefree};

static SOLVER_ABORT_NANOS: AtomicU64 = AtomicU64::new(0);

fn now_nanos() -> u64 {
    let elapsed = Instant::now().elapsed();
    u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX)
}

pub fn clear_solver_abort() {
    SOLVER_ABORT_NANOS.store(0, Ordering::Relaxed);
}

pub fn solver_abort_now() {
    SOLVER_ABORT_NANOS.store(now_nanos(), Ordering::Relaxed);
}

pub fn solver_aborted() -> bool {
    let limit = SOLVER_ABORT_NANOS.load(Ordering::Relaxed);
    limit != 0 && now_nanos() >= limit
}

pub fn budget_expired(deadline: Option<Instant>) -> bool {
    if solver_aborted() {
        return true;
    }
    deadline.is_some_and(|limit| Instant::now() >= limit)
}

const SOLVER_PRIMES: [u64; 24] = [
    7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
];

pub struct PolynomialForm {
    degree: usize,
    polynomial: Vec<Integer>,
    reduction: Vec<Vec<Integer>>,
    squarefree: std::cell::RefCell<std::collections::HashMap<u64, Option<Vec<ModPolynomial>>>>,
}

impl PolynomialForm {
    fn from_field(field: &NumberField) -> Option<Self> {
        let degree = field.degree;
        if degree < 2 || field.polynomial.len() != degree + 1 {
            return None;
        }
        let span = 2 * degree - 1;
        let mut reduction: Vec<Vec<Integer>> = (0..span)
            .map(|index| {
                let mut row = vec![Integer::from(0); degree];
                if index < degree {
                    row[index] = Integer::from(1);
                }
                row
            })
            .collect();
        for index in degree..span {
            for offset in 0..degree {
                let target = index - degree + offset;
                let factor = -field.polynomial[offset].clone();
                if target < degree {
                    reduction[index][target] += &factor;
                } else {
                    let source = reduction[target].clone();
                    for column in 0..degree {
                        reduction[index][column] += &factor * &source[column];
                    }
                }
            }
        }
        Some(Self {
            degree,
            polynomial: field.polynomial.clone(),
            reduction,
            squarefree: std::cell::RefCell::new(std::collections::HashMap::new()),
        })
    }

    fn squarefree_factors(&self, prime: u64) -> Option<Vec<ModPolynomial>> {
        if let Some(cached) = self.squarefree.borrow().get(&prime) {
            return cached.clone();
        }
        let computed = reduced_polynomial(self, prime).and_then(|reduced| {
            reduced
                .derivative()
                .gcd(&reduced)
                .is_constant()
                .then(|| factor_squarefree(&reduced))
        });
        self.squarefree.borrow_mut().insert(prime, computed.clone());
        computed
    }

    fn residuals(&self, value: &[Integer], alpha: &[Integer]) -> Vec<Integer> {
        let degree = self.degree;
        let mut convolution = vec![Integer::from(0); 2 * degree - 1];
        for i in 0..degree {
            for j in 0..degree {
                convolution[i + j] += &value[i] * &value[j];
            }
        }
        let mut out = vec![Integer::from(0); degree];
        for (index, row) in self.reduction.iter().enumerate() {
            if convolution[index].is_zero() {
                continue;
            }
            let weight = &convolution[index];
            for (column, factor) in row.iter().enumerate() {
                if factor.is_zero() {
                    continue;
                }
                out[column] += factor * weight;
            }
        }
        for index in 0..degree {
            out[index] -= &alpha[index];
        }
        out
    }

    fn holds_exact(&self, value: &[Integer], alpha: &[Integer]) -> bool {
        self.residuals(value, alpha)
            .iter()
            .all(|residual| residual.is_zero())
    }
}

const FIELD_ROOT_MISSES: usize = 3;
const SEED_LIMIT: usize = 4;

pub fn find_square_root(
    field: &NumberField,
    value: &Element,
    deadline: Option<Instant>,
) -> Option<Element> {
    if value.coefficients.len() != field.degree {
        return None;
    }
    let form = PolynomialForm::from_field(field)?;
    let alpha = &value.coefficients;
    let alpha_bits = alpha
        .iter()
        .map(integer_abs)
        .map(|coefficient| coefficient.bits())
        .max()
        .unwrap_or(0);
    let mut open: Vec<u64> = Vec::new();
    if let Some(root) = field_root(&form, alpha, alpha_bits / 2 + 32, &mut open, deadline) {
        return Some(root);
    }
    if let Some(root) = field_root(&form, alpha, alpha_bits + 16, &mut open, deadline) {
        return Some(root);
    }
    let capped: Vec<u64> = open
        .iter()
        .copied()
        .take(prime_budget(form.degree))
        .collect();
    if capped.is_empty() {
        return None;
    }
    if let Some(root) = hensel_root(&form, alpha, alpha_bits / 2 + 32, &capped, deadline) {
        return Some(root);
    }
    hensel_root(&form, alpha, alpha_bits + 16, &capped, deadline)
}

fn reduced_polynomial(form: &PolynomialForm, prime: u64) -> Option<ModPolynomial> {
    let coefficients: Vec<u64> = form
        .polynomial
        .iter()
        .map(|value| signed_residue(value, prime) as u64)
        .collect();
    let reduced = ModPolynomial::new(coefficients, prime).monic();
    (reduced.degree() == form.degree).then_some(reduced)
}

fn field_seed(form: &PolynomialForm, alpha: &[Integer], prime: u64) -> (Vec<Vec<i128>>, bool) {
    let Some(factors) = form.squarefree_factors(prime) else {
        return (Vec::new(), false);
    };
    if factors.is_empty() {
        return (Vec::new(), true);
    }
    let alpha_coefficients: Vec<u64> = alpha
        .iter()
        .map(|value| signed_residue(value, prime) as u64)
        .collect();
    let alpha_mod = ModPolynomial::new(alpha_coefficients, prime);
    let mut choices: Vec<Vec<ModPolynomial>> = Vec::with_capacity(factors.len());
    for factor in &factors {
        let Some(field) = ExtensionField::new(prime, factor.clone()) else {
            return (Vec::new(), true);
        };
        let value = field.normalize(&alpha_mod.rem(factor).coefficients);
        if !is_square(&field, &value) {
            return (Vec::new(), false);
        }
        let Some(root) = square_root(&field, &value) else {
            return (Vec::new(), false);
        };
        let negated: Vec<u64> = root
            .iter()
            .map(|value| (prime - value % prime) % prime)
            .collect();
        choices.push(vec![
            ModPolynomial::new(root, prime),
            ModPolynomial::new(negated, prime),
        ]);
    }
    let total = 1_usize << choices.len();
    let mut seeds = Vec::with_capacity(total);
    for index in 0..total {
        let residues: Vec<ModPolynomial> = choices
            .iter()
            .enumerate()
            .map(|(bit, options)| options[(index >> bit) & 1].clone())
            .collect();
        let Some(seed) = crt_reconstruct(&factors, &residues) else {
            return (Vec::new(), true);
        };
        let mut coefficients = seed.coefficients.clone();
        coefficients.resize(form.degree, 0);
        seeds.push(
            coefficients
                .iter()
                .map(|value| i128::from(*value % prime))
                .collect(),
        );
    }
    (seeds, false)
}

fn field_root(
    form: &PolynomialForm,
    alpha: &[Integer],
    bound_bits: u64,
    open: &mut Vec<u64>,
    deadline: Option<Instant>,
) -> Option<Element> {
    let mut misses = 0_usize;
    for &prime in SOLVER_PRIMES.iter() {
        if budget_expired(deadline) {
            return None;
        }
        if misses >= FIELD_ROOT_MISSES {
            return None;
        }
        let (seeds, worth) = field_seed(form, alpha, prime);
        if worth && !open.contains(&prime) {
            open.push(prime);
        }
        if seeds.is_empty() {
            continue;
        }
        let kernel = reduction_mod(form, prime);
        let mut modulus = Integer::from(prime);
        while modulus.bits() <= bound_bits {
            modulus *= prime;
        }
        let mut lifted = false;
        for seed in seeds.iter().take(SEED_LIMIT) {
            if let Some(candidate) =
                lift_seed(form, &kernel, seed, alpha, prime, &modulus, deadline)
            {
                return Some(candidate);
            }
            lifted = true;
        }
        if lifted {
            misses += 1;
        }
    }
    None
}

fn prime_budget(degree: usize) -> usize {
    const BUDGET: u128 = 20_000_000;
    let mut total = 0_u128;
    let mut count = 0_usize;
    for &prime in SOLVER_PRIMES.iter() {
        let cost = (prime as u128).pow(degree as u32);
        if total + cost > BUDGET {
            break;
        }
        total += cost;
        count += 1;
    }
    count.max(1)
}

fn hensel_root(
    form: &PolynomialForm,
    alpha: &[Integer],
    bound_bits: u64,
    primes: &[u64],
    deadline: Option<Instant>,
) -> Option<Element> {
    for &prime in primes.iter() {
        if budget_expired(deadline) {
            return None;
        }
        let mut modulus = Integer::from(prime);
        while modulus.bits() <= bound_bits {
            modulus *= prime;
        }
        if let Some(root) = hensel_for_prime(form, alpha, prime, &modulus, deadline) {
            return Some(root);
        }
    }
    None
}

fn hensel_for_prime(
    form: &PolynomialForm,
    alpha: &[Integer],
    prime: u64,
    modulus: &Integer,
    deadline: Option<Instant>,
) -> Option<Element> {
    let degree = form.degree;
    let p = prime as i128;
    let kernel = reduction_mod(form, prime);
    let residues: Vec<i128> = alpha.iter().map(|value| residue(value, prime)).collect();
    let mut seed = vec![0_i128; degree];
    loop {
        if budget_expired(deadline) {
            return None;
        }
        if system_holds_mod(form, &kernel, &seed, &residues, p) {
            let jacobian = jacobian_matrix(form, &kernel, &seed, p);
            if determinant(jacobian, p) != 0
                && let Some(candidate) =
                    lift_seed(form, &kernel, &seed, alpha, prime, modulus, deadline)
            {
                return Some(candidate);
            }
        }
        if !next_seed(&mut seed, p) {
            return None;
        }
    }
}

fn next_seed(seed: &mut [i128], p: i128) -> bool {
    for coordinate in seed.iter_mut() {
        *coordinate += 1;
        if *coordinate < p {
            return true;
        }
        *coordinate = 0;
    }
    false
}

fn lift_seed(
    form: &PolynomialForm,
    kernel: &[Vec<i128>],
    seed: &[i128],
    alpha: &[Integer],
    prime: u64,
    modulus: &Integer,
    deadline: Option<Instant>,
) -> Option<Element> {
    let degree = form.degree;
    let p = prime as i128;
    let mut value: Vec<Integer> = seed.iter().map(|item| Integer::from(*item)).collect();
    let mut place = Integer::from(prime);
    while &place < modulus {
        if budget_expired(deadline) {
            return None;
        }
        let residuals = form.residuals(&value, alpha);
        let mut rhs = vec![0_i128; degree];
        for index in 0..degree {
            let quotient = &residuals[index] / &place;
            let remainder = &residuals[index] % &place;
            if !remainder.is_zero() {
                return None;
            }
            rhs[index] = reduce(-signed_residue(&quotient, prime), p);
        }
        let current: Vec<i128> = value
            .iter()
            .map(|item| signed_residue(item, prime))
            .collect();
        let jacobian = jacobian_matrix(form, kernel, &current, p);
        let step = linear_solve(jacobian, rhs, p)?;
        for index in 0..degree {
            value[index] += &place * Integer::from(step[index]);
        }
        place *= prime;
    }
    let half = modulus / Integer::from(2);
    let negative_half = -half.clone();
    for coordinate in &mut value {
        if *coordinate > half {
            *coordinate -= modulus;
        } else if *coordinate < negative_half {
            *coordinate += modulus;
        }
    }
    if form.holds_exact(&value, alpha) {
        return Some(Element::from_coefficients(value));
    }
    None
}

fn reduction_mod(form: &PolynomialForm, prime: u64) -> Vec<Vec<i128>> {
    form.reduction
        .iter()
        .map(|row| row.iter().map(|value| residue(value, prime)).collect())
        .collect()
}

fn convolution_mod(seed: &[i128], degree: usize, p: i128) -> Vec<i128> {
    let mut convolution = vec![0_i128; 2 * degree - 1];
    for i in 0..degree {
        for j in 0..degree {
            convolution[i + j] = reduce(convolution[i + j] + seed[i] * seed[j], p);
        }
    }
    convolution
}

fn system_holds_mod(
    form: &PolynomialForm,
    kernel: &[Vec<i128>],
    seed: &[i128],
    residues: &[i128],
    p: i128,
) -> bool {
    let degree = form.degree;
    let convolution = convolution_mod(seed, degree, p);
    for column in 0..degree {
        let mut total = 0_i128;
        for index in 0..2 * degree - 1 {
            if kernel[index][column] == 0 || convolution[index] == 0 {
                continue;
            }
            total += kernel[index][column] * convolution[index];
        }
        if reduce(total - residues[column], p) != 0 {
            return false;
        }
    }
    true
}

#[allow(clippy::needless_range_loop)]
fn jacobian_matrix(
    form: &PolynomialForm,
    kernel: &[Vec<i128>],
    seed: &[i128],
    p: i128,
) -> Vec<Vec<i128>> {
    let degree = form.degree;
    let mut matrix = vec![vec![0_i128; degree]; degree];
    for column in 0..degree {
        for row in 0..degree {
            let mut total = 0_i128;
            for index in 0..2 * degree - 1 {
                let weight = kernel[index][row];
                if weight == 0 {
                    continue;
                }
                for k in 0..degree {
                    let other = index as isize - k as isize;
                    if other < 0 || other >= degree as isize {
                        continue;
                    }
                    let other = other as usize;
                    if k == column {
                        total += weight * seed[other];
                    }
                    if other == column {
                        total += weight * seed[k];
                    }
                }
            }
            matrix[row][column] = reduce(total, p);
        }
    }
    matrix
}
#[allow(clippy::needless_range_loop)]
fn determinant(matrix: Vec<Vec<i128>>, p: i128) -> i128 {
    let size = matrix.len();
    let mut work: Vec<Vec<i128>> = matrix
        .iter()
        .map(|row| row.iter().map(|value| reduce(*value, p)).collect())
        .collect();
    let mut result = 1_i128;
    for column in 0..size {
        let mut pivot = None;
        for row in column..size {
            if work[row][column] != 0 {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else {
            return 0;
        };
        if pivot != column {
            work.swap(pivot, column);
            result = -result;
        }
        result *= work[column][column];
        let Some(inverse) = modular_inverse(work[column][column], p) else {
            return 0;
        };
        for row in (column + 1)..size {
            if work[row][column] == 0 {
                continue;
            }
            let factor = reduce(work[row][column] * inverse, p);
            for index in column..size {
                let update = reduce(factor * work[column][index], p);
                work[row][index] = reduce(work[row][index] - update, p);
            }
        }
    }
    reduce(result, p)
}

#[allow(clippy::needless_range_loop)]
fn linear_solve(mut matrix: Vec<Vec<i128>>, mut rhs: Vec<i128>, p: i128) -> Option<Vec<i128>> {
    let size = matrix.len();
    for column in 0..size {
        let mut pivot = None;
        for row in column..size {
            if matrix[row][column] != 0 {
                pivot = Some(row);
                break;
            }
        }
        let pivot = pivot?;
        if pivot != column {
            matrix.swap(pivot, column);
            rhs.swap(pivot, column);
        }
        let inverse = modular_inverse(matrix[column][column], p)?;
        for index in column..size {
            matrix[column][index] = reduce(matrix[column][index] * inverse, p);
        }
        rhs[column] = reduce(rhs[column] * inverse, p);
        for row in 0..size {
            if row == column || matrix[row][column] == 0 {
                continue;
            }
            let factor = matrix[row][column];
            for index in column..size {
                let update = reduce(factor * matrix[column][index], p);
                matrix[row][index] = reduce(matrix[row][index] - update, p);
            }
            let update = reduce(factor * rhs[column], p);
            rhs[row] = reduce(rhs[row] - update, p);
        }
    }
    Some(rhs)
}

fn modular_inverse(value: i128, p: i128) -> Option<i128> {
    let (mut old_r, mut r) = (reduce(value, p), p);
    let (mut old_t, mut t) = (1_i128, 0_i128);
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_t, t) = (t, old_t - quotient * t);
    }
    if reduce(old_r, p) != 1 {
        return None;
    }
    Some(reduce(old_t, p))
}

fn residue(value: &Integer, prime: u64) -> i128 {
    signed_residue(value, prime)
}

fn signed_residue(value: &Integer, prime: u64) -> i128 {
    let reduced = value.rem_euclid(&Integer::from(prime));
    reduced.to_i128().unwrap_or(0)
}

fn reduce(value: i128, p: i128) -> i128 {
    value.rem_euclid(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pure_field(degree: usize, constant: i64) -> NumberField {
        let mut coefficients = vec![Integer::from(0); degree + 1];
        coefficients[0] = Integer::from(constant);
        coefficients[degree] = Integer::from(1);
        NumberField::new(degree, coefficients).expect("pure field")
    }

    fn mixed_field(degree: usize, values: &[i64]) -> NumberField {
        let mut coefficients: Vec<Integer> =
            values.iter().map(|value| Integer::from(*value)).collect();
        coefficients.resize(degree + 1, Integer::from(0));
        coefficients[degree] = Integer::from(1);
        NumberField::new(degree, coefficients).expect("mixed field")
    }

    #[test]
    fn hensel_solver_matches_cubic_on_known_root() {
        let field = pure_field(3, 2);
        let beta =
            Element::from_coefficients(vec![Integer::from(1), Integer::from(2), Integer::from(3)]);
        let alpha = field.multiply(&beta, &beta);
        let root = find_square_root(&field, &alpha, None).expect("hensel cubic root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }

    #[test]
    fn hensel_solver_finds_quartic_root() {
        let field = pure_field(4, 97);
        let beta = Element::from_coefficients(vec![
            Integer::from(3),
            Integer::from(-5),
            Integer::from(7),
            Integer::from(11),
        ]);
        let alpha = field.multiply(&beta, &beta);
        let root = find_square_root(&field, &alpha, None).expect("hensel quartic root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }

    #[test]
    fn hensel_solver_handles_large_quartic_coefficients() {
        let field = pure_field(4, 160063);
        let two = Integer::from(2);
        let beta = Element::from_coefficients(vec![
            Integer::from(5) * &two.pow(200_u32),
            Integer::from(3) * &two.pow(150_u32),
            Integer::from(7) * &two.pow(100_u32),
            Integer::from(11) * &two.pow(50_u32),
        ]);
        let alpha = field.multiply(&beta, &beta);
        let root = find_square_root(&field, &alpha, None).expect("large hensel quartic root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }

    #[test]
    fn hensel_solver_finds_quintic_and_sextic_roots() {
        for (degree, constant) in [(5_usize, 3_i64), (6_usize, 5_i64)] {
            let field = pure_field(degree, constant);
            let beta = Element::from_coefficients((1..=degree as i64).map(Integer::from).collect());
            let alpha = field.multiply(&beta, &beta);
            let root = find_square_root(&field, &alpha, None).expect("hensel higher root");
            assert_eq!(field.multiply(&root, &root), alpha);
        }
    }

    #[test]
    fn hensel_solver_handles_non_pure_quartic() {
        let field = mixed_field(4, &[-13, 2, -3, 1, 5]);
        let beta = Element::from_coefficients(vec![
            Integer::from(7),
            Integer::from(-4),
            Integer::from(9),
            Integer::from(2),
        ]);
        let alpha = field.multiply(&beta, &beta);
        let root = find_square_root(&field, &alpha, None).expect("non pure quartic root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }

    #[test]
    fn hensel_solver_handles_non_pure_cubic() {
        let field = mixed_field(3, &[-5, 2, -3, 1]);
        let beta = Element::from_coefficients(vec![
            Integer::from(11),
            Integer::from(-7),
            Integer::from(5),
        ]);
        let alpha = field.multiply(&beta, &beta);
        let root = find_square_root(&field, &alpha, None).expect("non pure cubic root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }

    #[test]
    fn hensel_solver_handles_large_non_pure_quartic_coefficients() {
        let field = mixed_field(4, &[-7, 3, -2, 1, 2]);
        let two = Integer::from(2);
        let beta = Element::from_coefficients(vec![
            Integer::from(5) * &two.pow(200_u32),
            Integer::from(9) * &two.pow(150_u32),
            Integer::from(13) * &two.pow(100_u32),
            Integer::from(3) * &two.pow(50_u32),
        ]);
        let alpha = field.multiply(&beta, &beta);
        let root = find_square_root(&field, &alpha, None).expect("large non pure quartic root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }
}
