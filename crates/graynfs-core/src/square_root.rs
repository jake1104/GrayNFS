use num_traits::Zero;

use crate::integer::{Integer, integer_abs, integer_sqrt};

use crate::algebraic::{QuadraticElement, QuadraticRing};
use crate::error::FactorError;
use crate::relations::Relation;

pub fn extract_factor(
    n: &Integer,
    ring: &QuadraticRing,
    evaluation_point: &Integer,
    relations: &[Relation],
    dependency: &[usize],
) -> Result<Integer, FactorError> {
    if dependency.is_empty() {
        return Err(FactorError::LinearAlgebra(
            "empty relation dependency".to_string(),
        ));
    }
    let mut rational_product = Integer::from(1);
    let mut algebraic_product = QuadraticElement {
        x: Integer::from(1),
        y: Integer::from(0),
    };
    for index in dependency {
        let relation = relations
            .get(*index)
            .ok_or_else(|| FactorError::LinearAlgebra("invalid relation index".to_string()))?;
        rational_product *= &relation.rational_value;
        algebraic_product = ring.multiply(&algebraic_product, &relation.algebraic_value);
    }
    let rational_abs = integer_abs(&rational_product);
    if rational_abs.is_zero() {
        return Err(FactorError::AlgebraicSquareRootFailed);
    }
    let rational_root = integer_sqrt(&rational_abs);
    if &rational_root * &rational_root != rational_abs {
        return Err(FactorError::AlgebraicSquareRootFailed);
    }
    let algebraic_root = ring
        .square_root(&algebraic_product)
        .ok_or(FactorError::AlgebraicSquareRootFailed)?;
    let root_variants = [
        algebraic_root.clone(),
        QuadraticElement {
            x: -algebraic_root.x.clone(),
            y: algebraic_root.y.clone(),
        },
        QuadraticElement {
            x: algebraic_root.x.clone(),
            y: -algebraic_root.y.clone(),
        },
        QuadraticElement {
            x: -algebraic_root.x.clone(),
            y: -algebraic_root.y.clone(),
        },
    ];
    for root in root_variants {
        let algebraic_value = ring.evaluate(&root, evaluation_point);
        for rational_value in [rational_root.clone(), -rational_root.clone()] {
            let x = algebraic_value.clone() % n;
            let y = rational_value % n;
            if (&x * &x - &y * &y) % n != Integer::from(0) {
                continue;
            }
            let difference = crate::integer::integer_gcd(&(&x - &y), n);
            if verified_divisor(n, &difference) {
                return Ok(difference);
            }
            let sum = crate::integer::integer_gcd(&(&x + &y), n);
            if verified_divisor(n, &sum) {
                return Ok(sum);
            }
        }
    }
    Err(FactorError::AlgebraicSquareRootFailed)
}

fn verified_divisor(n: &Integer, divisor: &Integer) -> bool {
    if divisor <= &Integer::from(1) || divisor >= n || n % divisor != Integer::from(0) {
        return false;
    }
    let cofactor = n / divisor;
    &cofactor * divisor == *n
}
