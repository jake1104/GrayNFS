use num_traits::{ToPrimitive, Zero};
use std::sync::atomic::AtomicBool;
use std::time::Instant;

use crate::hensel::budget_expired;
use crate::integer::Integer;
use crate::number_field::NumberField;
use crate::poly::mod_inverse;

pub fn lll_reduce(
    mut basis: Vec<Vec<Integer>>,
    delta: f64,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> Vec<Vec<Integer>> {
    if basis.len() < 2 {
        return basis;
    }
    let delta = if delta > 0.25 { delta } else { 0.75 };
    let max_iterations = 64 * basis.len() * basis.len() + 64;
    for _ in 0..max_iterations {
        if budget_expired(deadline, stop) {
            break;
        }
        let mut changed = false;
        for k in 0..basis.len().saturating_sub(1) {
            let norm = vector_norm_f64(&basis[k]);
            if norm <= 0.0 {
                continue;
            }
            let pivot = basis[k].clone();
            for j in (k + 1)..basis.len() {
                let denominator = vector_norm_f64(&pivot);
                if denominator <= 0.0 {
                    continue;
                }
                let numerator = dot(&pivot, &basis[j]).to_f64().unwrap_or(0.0);
                let mu = (numerator / denominator).round() as i64;
                if mu != 0 {
                    let before = basis[j].clone();
                    for index in 0..basis[k].len() {
                        basis[j][index] -= &pivot[index] * Integer::from(mu);
                    }
                    if before != basis[j] {
                        changed = true;
                    }
                }
            }
            if k > 0 {
                let previous = vector_norm_f64(&basis[k - 1]);
                let mu = if previous > 0.0 {
                    dot(&basis[k - 1], &basis[k]).to_f64().unwrap_or(0.0) / previous
                } else {
                    0.0
                };
                if vector_norm_f64(&basis[k]) < (delta - mu * mu) * previous {
                    basis.swap(k, k + 1);
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            break;
        }
    }
    basis
}

pub fn ideal_basis(
    degree: usize,
    prime: u64,
    root: u64,
    field: &NumberField,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> Vec<Vec<Integer>> {
    let mut basis = Vec::with_capacity(degree);
    let mut first = vec![Integer::from(0); degree];
    if root == 0 {
        first[1.min(degree - 1)] = Integer::from(1);
    } else {
        first[0] = Integer::from(1);
        first[1] = -Integer::from(mod_inverse(root % prime, prime));
    }
    basis.push(first);
    for index in 1..degree {
        let mut vector = vec![Integer::from(0); degree];
        vector[index] = Integer::from(prime);
        basis.push(vector);
    }
    let _ = field;
    lll_reduce(basis, 0.75, deadline, stop)
}

pub fn shortest_vector(
    basis: Vec<Vec<Integer>>,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> Option<Vec<Integer>> {
    lll_reduce(basis, 0.75, deadline, stop)
        .into_iter()
        .filter(|vector| vector.iter().any(|value| !value.is_zero()))
        .min_by_key(|vector| vector_norm_squared(vector))
}

fn dot(left: &[Integer], right: &[Integer]) -> Integer {
    left.iter()
        .zip(right.iter())
        .map(|(left, right)| left * right)
        .sum()
}

fn vector_norm_squared(vector: &[Integer]) -> Integer {
    dot(vector, vector)
}

fn vector_norm_f64(vector: &[Integer]) -> f64 {
    vector_norm_squared(vector)
        .to_f64()
        .unwrap_or(f64::INFINITY)
        .sqrt()
}
