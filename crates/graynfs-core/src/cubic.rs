use num_traits::{Euclid, ToPrimitive, Zero};
use std::time::Instant;

use crate::hensel::budget_expired;
use crate::integer::{Integer, integer_abs};
use crate::number_field::{Element, NumberField};

const CUBIC_SOLVER_PRIMES: [u64; 24] = [
    7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
];

#[derive(Clone)]
pub struct CubicForm {
    constant: Integer,
    linear: Integer,
    quadratic: Integer,
}

#[derive(Clone, Copy)]
struct CubicMod {
    two_constant: i128,
    two_linear: i128,
    two_quadratic: i128,
    constant_quadratic: i128,
    quadratic_linear: i128,
    quadratic_squared: i128,
}

impl CubicForm {
    fn from_field(field: &NumberField) -> Option<Self> {
        if field.degree != 3 || field.polynomial.len() < 4 {
            return None;
        }
        let constant = field.polynomial[0].clone();
        if constant.is_zero() {
            return None;
        }
        Some(Self {
            constant,
            linear: field.polynomial[1].clone(),
            quadratic: field.polynomial[2].clone(),
        })
    }

    fn to_mod(&self, prime: u64) -> CubicMod {
        let p = prime as i128;
        let a0 = residue(&self.constant, prime);
        let a1 = residue(&self.linear, prime);
        let a2 = residue(&self.quadratic, prime);
        CubicMod {
            two_constant: reduce(2 * a0, p),
            two_linear: reduce(2 * a1, p),
            two_quadratic: reduce(2 * a2, p),
            constant_quadratic: reduce(a0 * a2, p),
            quadratic_linear: reduce(a2 * a1 - a0, p),
            quadratic_squared: reduce(a2 * a2 - a1, p),
        }
    }
}

pub fn find_square_root(
    field: &NumberField,
    value: &Element,
    deadline: Option<Instant>,
) -> Option<Element> {
    if value.coefficients.len() != 3 {
        return None;
    }
    let form = CubicForm::from_field(field)?;
    let alpha = &value.coefficients;
    let alpha_bits = alpha
        .iter()
        .map(integer_abs)
        .map(|coefficient| coefficient.bits())
        .max()
        .unwrap_or(0);
    if let Some(root) = hensel_cubic_root(&form, alpha, (alpha_bits / 2 + 32) as usize, deadline) {
        return Some(root);
    }
    hensel_cubic_root(&form, alpha, (alpha_bits + 16) as usize, deadline)
}

fn hensel_cubic_root(
    form: &CubicForm,
    alpha: &[Integer],
    bound_bits: usize,
    deadline: Option<Instant>,
) -> Option<Element> {
    for &prime in CUBIC_SOLVER_PRIMES.iter() {
        if budget_expired(deadline) {
            return None;
        }
        let mut modulus = Integer::from(prime);
        while modulus.bits() <= bound_bits as u64 {
            modulus *= prime;
        }
        if let Some(root) = hensel_for_prime(form, alpha, prime, &modulus, deadline) {
            return Some(root);
        }
    }
    None
}

fn hensel_for_prime(
    form: &CubicForm,
    alpha: &[Integer],
    prime: u64,
    modulus: &Integer,
    deadline: Option<Instant>,
) -> Option<Element> {
    let p = prime as i128;
    let coefficients = form.to_mod(prime);
    let target = [
        residue(&alpha[0], prime),
        residue(&alpha[1], prime),
        residue(&alpha[2], prime),
    ];
    let current = Integer::from(prime);
    for z in 0..p {
        if budget_expired(deadline) {
            return None;
        }
        for x in 0..p {
            for y in 0..p {
                if !system_holds_mod(x, y, z, coefficients, target, p) {
                    continue;
                }
                let jacobian = jacobian_matrix(x, y, z, coefficients, p);
                if determinant(jacobian) % p == 0 {
                    continue;
                }
                if let Some(candidate) = lift_seed(
                    [x as u64, y as u64, z as u64],
                    form,
                    alpha,
                    prime,
                    &current,
                    modulus,
                ) {
                    return Some(candidate);
                }
            }
        }
    }
    None
}

fn lift_seed(
    seed: [u64; 3],
    form: &CubicForm,
    alpha: &[Integer],
    prime: u64,
    current: &Integer,
    modulus: &Integer,
) -> Option<Element> {
    let p = prime as i128;
    let coefficients = form.to_mod(prime);
    let mut value = [
        Integer::from(seed[0]),
        Integer::from(seed[1]),
        Integer::from(seed[2]),
    ];
    let mut place = current.clone();
    let zero = [Integer::from(0), Integer::from(0), Integer::from(0)];
    let initial_residuals = exact_residuals(form, alpha, &value);
    let mut quotients = [Integer::from(0), Integer::from(0), Integer::from(0)];
    for index in 0..3 {
        quotients[index] = &initial_residuals[index] / &place;
        let remainder = &initial_residuals[index] - &quotients[index] * &place;
        if !remainder.is_zero() {
            return None;
        }
    }
    let jacobian = jacobian_matrix(
        signed_residue(&value[0], prime),
        signed_residue(&value[1], prime),
        signed_residue(&value[2], prime),
        coefficients,
        p,
    );
    let prime_inverse = inverse_mod_prime(jacobian, p)?;
    let mut inverse = [
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
    ];
    for row in 0..3 {
        for column in 0..3 {
            inverse[row * 3 + column] = Integer::from(prime_inverse[row][column]);
        }
    }
    while &place < modulus {
        let next_place = &place * &place;
        let mut residual = [Integer::from(0), Integer::from(0), Integer::from(0)];
        for index in 0..3 {
            residual[index] = quotients[index].rem_euclid(&place);
        }
        let product = matrix_vector_mod(&inverse, &residual, &place);
        let mut step = [Integer::from(0), Integer::from(0), Integer::from(0)];
        for index in 0..3 {
            step[index] = (&place - &product[index]).rem_euclid(&place);
        }
        let mut next_value = [Integer::from(0), Integer::from(0), Integer::from(0)];
        for index in 0..3 {
            next_value[index] = &value[index] + &place * &step[index];
        }
        let linear = jacobian_step(form, &value, &step);
        let quadratic = exact_residuals(form, &zero, &step);
        let mut next_quotients = [Integer::from(0), Integer::from(0), Integer::from(0)];
        for index in 0..3 {
            let numerator = &quotients[index] + &linear[index];
            let divided = &numerator / &place;
            let remainder = &numerator - &divided * &place;
            if !remainder.is_zero() {
                return None;
            }
            next_quotients[index] = divided + &quadratic[index];
        }
        let exact = jacobian_exact(form, &next_value);
        let mut jacobian_wide = [
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ];
        for index in 0..9 {
            jacobian_wide[index] = exact[index].rem_euclid(&next_place);
        }
        let product_matrix = matrix_multiply_mod(&jacobian_wide, &inverse, &next_place);
        let identity = matrix_identity();
        let mut error = [
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ];
        for index in 0..9 {
            error[index] = (&identity[index] - &product_matrix[index]).rem_euclid(&next_place);
        }
        let correction = matrix_multiply_mod(&inverse, &error, &next_place);
        let mut next_inverse = [
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ];
        for index in 0..9 {
            next_inverse[index] = (&inverse[index] + &correction[index]).rem_euclid(&next_place);
        }
        value = next_value;
        quotients = next_quotients;
        inverse = next_inverse;
        place = next_place;
    }
    let half = modulus / Integer::from(2);
    let negative_half = -half.clone();
    for coordinate in &mut value {
        let reduced = coordinate.rem_euclid(modulus);
        *coordinate = reduced;
        if *coordinate > half {
            *coordinate -= modulus;
        } else if *coordinate < negative_half {
            *coordinate += modulus;
        }
    }
    if system_holds_exact(form, alpha, &value) {
        return Some(Element::from_coefficients(value.to_vec()));
    }
    None
}

fn exact_residuals(form: &CubicForm, alpha: &[Integer], value: &[Integer; 3]) -> [Integer; 3] {
    let two = Integer::from(2);
    let x = &value[0];
    let y = &value[1];
    let z = &value[2];
    let yz = y * z;
    let z_squared = z * z;
    let first = x * x - &two * &form.constant * &yz + &form.constant * &form.quadratic * &z_squared
        - &alpha[0];
    let second = &two * x * y - &two * &form.linear * &yz
        + (&form.quadratic * &form.linear - &form.constant) * &z_squared
        - &alpha[1];
    let third = y * y + &two * x * z - &two * &form.quadratic * &yz
        + (&form.quadratic * &form.quadratic - &form.linear) * &z_squared
        - &alpha[2];
    [first, second, third]
}

fn jacobian_step(form: &CubicForm, value: &[Integer; 3], step: &[Integer; 3]) -> [Integer; 3] {
    let two = Integer::from(2);
    let x = &value[0];
    let y = &value[1];
    let z = &value[2];
    let dx = &step[0];
    let dy = &step[1];
    let dz = &step[2];
    let first = &two * x * dx - &two * &form.constant * z * dy
        + (&two * &form.constant * &form.quadratic * z - &two * &form.constant * y) * dz;
    let second = &two * y * dx
        + (&two * x - &two * &form.linear * z) * dy
        + (&two * &form.quadratic * &form.linear * z
            - &two * &form.linear * y
            - &two * &form.constant * z)
            * dz;
    let third = &two * z * dx
        + (&two * y - &two * &form.quadratic * z) * dy
        + (&two * x - &two * &form.quadratic * y + &two * &form.quadratic * &form.quadratic * z
            - &two * &form.linear * z)
            * dz;
    [first, second, third]
}

fn system_holds_mod(
    x: i128,
    y: i128,
    z: i128,
    coefficients: CubicMod,
    target: [i128; 3],
    p: i128,
) -> bool {
    let z_squared = z * z;
    let yz = y * z;
    let first = reduce(
        x * x - coefficients.two_constant * yz + coefficients.constant_quadratic * z_squared
            - target[0],
        p,
    );
    let second = reduce(
        2 * x * y - coefficients.two_linear * yz + coefficients.quadratic_linear * z_squared
            - target[1],
        p,
    );
    let third = reduce(
        y * y + 2 * x * z - coefficients.two_quadratic * yz
            + coefficients.quadratic_squared * z_squared
            - target[2],
        p,
    );
    first == 0 && second == 0 && third == 0
}

fn system_holds_exact(form: &CubicForm, alpha: &[Integer], value: &[Integer; 3]) -> bool {
    exact_residuals(form, alpha, value)
        .iter()
        .all(|residual| residual.is_zero())
}

fn jacobian_matrix(x: i128, y: i128, z: i128, coefficients: CubicMod, p: i128) -> [[i128; 3]; 3] {
    [
        [
            reduce(2 * x, p),
            reduce(-coefficients.two_constant * z, p),
            reduce(
                -coefficients.two_constant * y + 2 * coefficients.constant_quadratic * z,
                p,
            ),
        ],
        [
            reduce(2 * y, p),
            reduce(2 * x - coefficients.two_linear * z, p),
            reduce(
                -coefficients.two_linear * y + 2 * coefficients.quadratic_linear * z,
                p,
            ),
        ],
        [
            reduce(2 * z, p),
            reduce(2 * y - coefficients.two_quadratic * z, p),
            reduce(
                2 * x - coefficients.two_quadratic * y + 2 * coefficients.quadratic_squared * z,
                p,
            ),
        ],
    ]
}

fn determinant(matrix: [[i128; 3]; 3]) -> i128 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
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

fn jacobian_exact(form: &CubicForm, value: &[Integer; 3]) -> [Integer; 9] {
    let two = Integer::from(2);
    let x = &value[0];
    let y = &value[1];
    let z = &value[2];
    let twice_constant = &two * &form.constant;
    let twice_linear = &two * &form.linear;
    let twice_quadratic = &two * &form.quadratic;
    [
        &two * x,
        -&twice_constant * z,
        -&twice_constant * y + &two * &form.constant * &form.quadratic * z,
        &two * y,
        &two * x - &twice_linear * z,
        -&twice_linear * y + &two * (&form.quadratic * &form.linear - &form.constant) * z,
        &two * z,
        &two * y - &twice_quadratic * z,
        &two * x - &twice_quadratic * y
            + &two * (&form.quadratic * &form.quadratic - &form.linear) * z,
    ]
}

fn matrix_identity() -> [Integer; 9] {
    let mut result = [
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
    ];
    result[0] = Integer::from(1);
    result[4] = Integer::from(1);
    result[8] = Integer::from(1);
    result
}

fn matrix_multiply_mod(
    left: &[Integer; 9],
    right: &[Integer; 9],
    modulus: &Integer,
) -> [Integer; 9] {
    let mut result = [
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
    ];
    for row in 0..3 {
        for column in 0..3 {
            let mut sum = Integer::from(0);
            for inner in 0..3 {
                sum += &left[row * 3 + inner] * &right[inner * 3 + column];
            }
            result[row * 3 + column] = sum.rem_euclid(modulus);
        }
    }
    result
}

fn matrix_vector_mod(
    matrix: &[Integer; 9],
    vector: &[Integer; 3],
    modulus: &Integer,
) -> [Integer; 3] {
    let mut result = [Integer::from(0), Integer::from(0), Integer::from(0)];
    for row in 0..3 {
        let mut sum = Integer::from(0);
        for column in 0..3 {
            sum += &matrix[row * 3 + column] * &vector[column];
        }
        result[row] = sum.rem_euclid(modulus);
    }
    result
}

fn inverse_mod_prime(matrix: [[i128; 3]; 3], p: i128) -> Option<[[i128; 3]; 3]> {
    let det = reduce(determinant(matrix), p);
    if det == 0 {
        return None;
    }
    let inverse = modular_inverse(det, p)?;
    let others = [[1, 2], [0, 2], [0, 1]];
    let mut cofactor = [[0_i128; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let minor = matrix[others[row][0]][others[column][0]]
                * matrix[others[row][1]][others[column][1]]
                - matrix[others[row][0]][others[column][1]]
                    * matrix[others[row][1]][others[column][0]];
            cofactor[row][column] = if (row + column) % 2 == 0 {
                minor
            } else {
                -minor
            };
        }
    }
    let mut result = [[0_i128; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            result[row][column] = reduce(cofactor[column][row] * inverse, p);
        }
    }
    Some(result)
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

    fn assert_root(field: &NumberField, beta: &Element) {
        let alpha = field.multiply(beta, beta);
        let root = find_square_root(field, &alpha, None).expect("cubic square root");
        assert_eq!(field.multiply(&root, &root), alpha);
    }

    #[test]
    fn cubic_solver_finds_known_root() {
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
        assert_root(
            &field,
            &Element::from_coefficients(vec![Integer::from(1), Integer::from(2), Integer::from(3)]),
        );
    }

    #[test]
    fn cubic_solver_handles_large_coefficients() {
        let field = NumberField::new(
            3,
            vec![
                Integer::from(3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ],
        )
        .expect("cubic field");
        let two = Integer::from(2);
        let shift = 400_u32;
        assert_root(
            &field,
            &Element::from_coefficients(vec![
                Integer::from(1) * &two.pow(shift),
                Integer::from(3) * &two.pow(shift / 2),
                Integer::from(7) * &two.pow(shift / 3),
            ]),
        );
    }

    #[test]
    fn cubic_solver_handles_non_pure_polynomial() {
        let field = NumberField::new(
            3,
            vec![
                Integer::from(-5),
                Integer::from(2),
                Integer::from(-3),
                Integer::from(1),
            ],
        )
        .expect("cubic field");
        assert_root(
            &field,
            &Element::from_coefficients(vec![
                Integer::from(11),
                Integer::from(-7),
                Integer::from(5),
            ]),
        );
    }

    #[test]
    fn cubic_solver_handles_large_non_pure_coefficients() {
        let field = NumberField::new(
            3,
            vec![
                Integer::from(-7),
                Integer::from(3),
                Integer::from(-2),
                Integer::from(1),
            ],
        )
        .expect("cubic field");
        let two = Integer::from(2);
        let shift = 300_u32;
        assert_root(
            &field,
            &Element::from_coefficients(vec![
                Integer::from(5) * &two.pow(shift),
                Integer::from(9) * &two.pow(shift / 2),
                Integer::from(13) * &two.pow(shift / 3),
            ]),
        );
    }
}
