use num_traits::Zero;
use std::sync::atomic::AtomicBool;
use std::time::Instant;

use crate::hensel::budget_expired;

const I64_LIMIT: u128 = i64::MAX as u128;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LatticeBasis {
    pub first: (i64, i64),
    pub second: (i64, i64),
}

impl LatticeBasis {
    pub fn point(&self, u: i64, v: i64) -> (i128, i128) {
        (
            i128::from(u) * i128::from(self.first.0) + i128::from(v) * i128::from(self.second.0),
            i128::from(u) * i128::from(self.first.1) + i128::from(v) * i128::from(self.second.1),
        )
    }
}

fn form(a: i128, b: i128, coefficient: i128) -> i128 {
    a * a + coefficient * b * b
}

fn inner(a: (i128, i128), b: (i128, i128), coefficient: i128) -> i128 {
    a.0 * b.0 + coefficient * a.1 * b.1
}

pub fn reduce_basis(
    mut first: (i64, i64),
    mut second: (i64, i64),
    coefficient: i128,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> LatticeBasis {
    if coefficient <= 0 {
        return LatticeBasis { first, second };
    }
    let iterations = 256;
    for _ in 0..iterations {
        if budget_expired(deadline, stop) {
            break;
        }
        let wide_first = (i128::from(first.0), i128::from(first.1));
        let wide_second = (i128::from(second.0), i128::from(second.1));
        let first_norm = form(wide_first.0, wide_first.1, coefficient);
        let second_norm = form(wide_second.0, wide_second.1, coefficient);
        if first_norm <= 0 || second_norm <= 0 {
            break;
        }
        if second_norm < first_norm {
            std::mem::swap(&mut first, &mut second);
            continue;
        }
        let current: (i128, i128) = (i128::from(first.0), i128::from(first.1));
        let other: (i128, i128) = (i128::from(second.0), i128::from(second.1));
        let pivot = form(current.0, current.1, coefficient);
        let projection = inner(current, other, coefficient);
        let factor = round_div(projection, pivot);
        if factor == 0 {
            break;
        }
        let new_a = other.0 - factor * current.0;
        let new_b = other.1 - factor * current.1;
        if new_a.unsigned_abs() > I64_LIMIT || new_b.unsigned_abs() > I64_LIMIT {
            break;
        }
        second = (new_a as i64, new_b as i64);
    }
    LatticeBasis { first, second }
}

fn round_div(numerator: i128, denominator: i128) -> i128 {
    if denominator == 0 {
        return 0;
    }
    let doubled = 2 * numerator.abs();
    if doubled >= denominator.abs() {
        let sign = if numerator < 0 { -1 } else { 1 };
        sign * ((doubled / denominator.abs() + 1) / 2)
    } else {
        0
    }
}

pub fn lattice_for_prime(
    prime: u64,
    root: u64,
    coefficient: i128,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> LatticeBasis {
    let prime_value = i64::try_from(prime).unwrap_or(i64::MAX);
    let root_value = i64::try_from(root).unwrap_or(0);
    reduce_basis(
        (prime_value, 0),
        (root_value, 1),
        coefficient,
        deadline,
        stop,
    )
}

pub fn lattice_size(basis: &LatticeBasis, coefficient: i128, prime: u64) -> (usize, usize) {
    let first_norm = form(
        i128::from(basis.first.0),
        i128::from(basis.first.1),
        coefficient,
    );
    let second_norm = form(
        i128::from(basis.second.0),
        i128::from(basis.second.1),
        coefficient,
    );
    let prime_value = i128::from(prime.max(1));
    let first_extent = (prime_value / first_norm.max(1)).max(1);
    let second_extent = (prime_value / second_norm.max(1)).max(1);
    (
        usize::try_from(first_extent).unwrap_or(1).clamp(1, 4096),
        usize::try_from(second_extent).unwrap_or(1).clamp(1, 4096),
    )
}

pub fn rescale_bounds(
    basis: &LatticeBasis,
    coefficient: i128,
    prime: u64,
    scale: i128,
) -> (usize, usize) {
    let (first, second) = lattice_size(basis, coefficient, prime);
    let factor = usize::try_from(scale.clamp(1, 64)).unwrap_or(1);
    (
        first.saturating_mul(factor).clamp(1, 8192),
        second.saturating_mul(factor).clamp(1, 8192),
    )
}

pub fn zero_point() -> (i128, i128) {
    (0, 0)
}

pub fn point_is_zero(point: (i128, i128)) -> bool {
    point.0.is_zero() && point.1.is_zero()
}
