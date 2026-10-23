use crate::cascade::{ecm_plan, escalate_digits};
use crate::integer::{Integer, integer_gcd};
use crate::montgomery::Montgomery;
use num_traits::Euclid;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

const STAGE_TWO_CAP: u64 = 16_000_000;
const MAX_DIGITS: u32 = 60;
const TORSION_COMPENSATION: u64 = 3;
const STAGE_TWO_TABLE_BYTES: u64 = 64_000_000;

fn prime_table(limit: u64) -> Vec<bool> {
    let size = limit.saturating_add(1) as usize;
    if size < 3 {
        return vec![false; size];
    }
    let mut sieve = vec![true; size];
    sieve[0] = false;
    sieve[1] = false;
    let mut value = 2_usize;
    while value * value < size {
        if sieve[value] {
            let mut multiple = value * value;
            while multiple < size {
                sieve[multiple] = false;
                multiple += value;
            }
        }
        value += 1;
    }
    sieve
}

fn stage_two_prime_table(bound_one: u64, bound_two: u64) -> Option<Vec<bool>> {
    let capped = bound_two.min(bound_one.saturating_mul(16));
    let width = (capped as f64).sqrt() as u64 + 1;
    let entries = capped
        .saturating_add(width.saturating_mul(2))
        .saturating_add(1);
    if entries > STAGE_TWO_TABLE_BYTES {
        return None;
    }
    Some(prime_table(entries.saturating_sub(1)))
}

fn extended_gcd(first: &Integer, second: &Integer) -> (Integer, Integer, Integer) {
    let mut old_remainder = first.clone();
    let mut remainder = second.clone();
    let mut old_coefficient = Integer::from(1);
    let mut coefficient = Integer::from(0);
    while remainder != Integer::from(0) {
        let quotient = &old_remainder / &remainder;
        let next_remainder = &old_remainder - &quotient * &remainder;
        old_remainder = remainder;
        remainder = next_remainder;
        let next_coefficient = &old_coefficient - &quotient * &coefficient;
        old_coefficient = coefficient;
        coefficient = next_coefficient;
    }
    (old_remainder, old_coefficient, coefficient)
}

fn mod_inverse(value: &Integer, modulus: &Integer) -> Option<Integer> {
    let (divisor, coefficient, _) = extended_gcd(value, modulus);
    if divisor == Integer::from(1) || divisor == Integer::from(-1) {
        Some(coefficient.rem_euclid(modulus))
    } else {
        None
    }
}

#[cfg(test)]
fn double_point(
    field: &Montgomery,
    point_x: &[u64],
    point_z: &[u64],
    curve_a24: &[u64],
) -> (Vec<u64>, Vec<u64>) {
    let mut workspace = ecm_workspace(field);
    let mut next_x = vec![0_u64; field.limbs()];
    let mut next_z = vec![0_u64; field.limbs()];
    double_point_into(
        field,
        point_x,
        point_z,
        curve_a24,
        &mut next_x,
        &mut next_z,
        &mut workspace,
    );
    (next_x, next_z)
}

#[cfg(test)]
fn add_points(
    field: &Montgomery,
    first_x: &[u64],
    first_z: &[u64],
    second_x: &[u64],
    second_z: &[u64],
    diff_x: &[u64],
    diff_z: &[u64],
) -> (Vec<u64>, Vec<u64>) {
    let mut workspace = ecm_workspace(field);
    let mut next_x = vec![0_u64; field.limbs()];
    let mut next_z = vec![0_u64; field.limbs()];
    add_points_into(
        field,
        first_x,
        first_z,
        second_x,
        second_z,
        diff_x,
        diff_z,
        &mut next_x,
        &mut next_z,
        &mut workspace,
    );
    (next_x, next_z)
}

pub(crate) struct EcmWorkspace {
    middle: Vec<u64>,
    temps: [Vec<u64>; 6],
}

pub(crate) fn ecm_workspace(field: &Montgomery) -> EcmWorkspace {
    let limbs = field.limbs();
    EcmWorkspace {
        middle: vec![0_u64; field.multiply_scratch_len()],
        temps: [(); 6].map(|_| vec![0_u64; limbs]),
    }
}

pub(crate) fn double_point_into(
    field: &Montgomery,
    point_x: &[u64],
    point_z: &[u64],
    curve_a24: &[u64],
    out_x: &mut [u64],
    out_z: &mut [u64],
    workspace: &mut EcmWorkspace,
) {
    let EcmWorkspace { middle, temps } = workspace;
    let [t0, t1, t2, t3, t4, _t5] = temps;
    let scratch = middle.as_mut_slice();
    field.add_into(t0, point_x, point_z);
    field.subtract_into(t1, point_x, point_z);
    field.multiply_into(t2, t0, t0, scratch);
    field.multiply_into(t3, t1, t1, scratch);
    field.subtract_into(t4, t2, t3);
    field.multiply_into(out_x, t2, t3, scratch);
    field.multiply_into(t0, curve_a24, t4, scratch);
    field.add_into(t1, t3, t0);
    field.multiply_into(out_z, t4, t1, scratch);
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn add_points_into(
    field: &Montgomery,
    first_x: &[u64],
    first_z: &[u64],
    second_x: &[u64],
    second_z: &[u64],
    diff_x: &[u64],
    diff_z: &[u64],
    out_x: &mut [u64],
    out_z: &mut [u64],
    workspace: &mut EcmWorkspace,
) {
    let EcmWorkspace { middle, temps } = workspace;
    let [t0, t1, t2, t3, _t4, _t5] = temps;
    let scratch = middle.as_mut_slice();
    field.subtract_into(t0, first_x, first_z);
    field.add_into(t1, second_x, second_z);
    field.multiply_into(t2, t0, t1, scratch);
    field.add_into(t0, first_x, first_z);
    field.subtract_into(t1, second_x, second_z);
    field.multiply_into(t3, t0, t1, scratch);
    field.add_into(t0, t2, t3);
    field.subtract_into(t1, t2, t3);
    field.multiply_into(t2, t0, t0, scratch);
    field.multiply_into(out_x, diff_z, t2, scratch);
    field.multiply_into(t2, t1, t1, scratch);
    field.multiply_into(out_z, diff_x, t2, scratch);
}

pub(crate) struct LadderScratch {
    low_x: Vec<u64>,
    low_z: Vec<u64>,
    high_x: Vec<u64>,
    high_z: Vec<u64>,
    doubled_x: Vec<u64>,
    doubled_z: Vec<u64>,
    added_x: Vec<u64>,
    added_z: Vec<u64>,
    next_x: Vec<u64>,
    next_z: Vec<u64>,
    spare_x: Vec<u64>,
    spare_z: Vec<u64>,
}

pub(crate) fn ladder_scratch(field: &Montgomery) -> LadderScratch {
    let buffers = [(); 12].map(|_| vec![0_u64; field.limbs()]);
    let [
        low_x,
        low_z,
        high_x,
        high_z,
        doubled_x,
        doubled_z,
        added_x,
        added_z,
        next_x,
        next_z,
        spare_x,
        spare_z,
    ] = buffers;
    LadderScratch {
        low_x,
        low_z,
        high_x,
        high_z,
        doubled_x,
        doubled_z,
        added_x,
        added_z,
        next_x,
        next_z,
        spare_x,
        spare_z,
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn ladder_mult_into(
    field: &Montgomery,
    scratch: &mut LadderScratch,
    base_x: &[u64],
    base_z: &[u64],
    scalar: u64,
    curve_a24: &[u64],
    one_repr: &[u64],
    out_x: &mut [u64],
    out_z: &mut [u64],
    workspace: &mut EcmWorkspace,
) {
    let LadderScratch {
        low_x,
        low_z,
        high_x,
        high_z,
        doubled_x,
        doubled_z,
        added_x,
        added_z,
        next_x,
        next_z,
        spare_x: _,
        spare_z: _,
    } = scratch;
    low_x.copy_from_slice(one_repr);
    low_z.fill(0);
    high_x.copy_from_slice(base_x);
    high_z.copy_from_slice(base_z);
    let width = 64 - scalar.leading_zeros();
    for position in (0..width).rev() {
        double_point_into(
            field, low_x, low_z, curve_a24, doubled_x, doubled_z, workspace,
        );
        add_points_into(
            field, low_x, low_z, high_x, high_z, base_x, base_z, added_x, added_z, workspace,
        );
        if (scalar >> position) & 1 == 0 {
            std::mem::swap(low_x, doubled_x);
            std::mem::swap(low_z, doubled_z);
            std::mem::swap(high_x, added_x);
            std::mem::swap(high_z, added_z);
        } else {
            double_point_into(field, high_x, high_z, curve_a24, next_x, next_z, workspace);
            std::mem::swap(low_x, added_x);
            std::mem::swap(low_z, added_z);
            std::mem::swap(high_x, next_x);
            std::mem::swap(high_z, next_z);
        }
    }
    out_x.copy_from_slice(low_x);
    out_z.copy_from_slice(low_z);
}

#[cfg(test)]
fn ladder_mult(
    field: &Montgomery,
    base_x: &[u64],
    base_z: &[u64],
    scalar: u64,
    curve_a24: &[u64],
) -> (Vec<u64>, Vec<u64>) {
    let mut workspace = ecm_workspace(field);
    let mut scratch = ladder_scratch(field);
    let one_repr = field.convert_in(&Integer::from(1));
    let mut out_x = vec![0_u64; field.limbs()];
    let mut out_z = vec![0_u64; field.limbs()];
    ladder_mult_into(
        field,
        &mut scratch,
        base_x,
        base_z,
        scalar,
        curve_a24,
        &one_repr,
        &mut out_x,
        &mut out_z,
        &mut workspace,
    );
    (out_x, out_z)
}

#[allow(clippy::excessive_precision)]
const PRAC_RATIOS: [f64; 10] = [
    0.61803398874989485,
    0.72360679774997897,
    0.58017872829546410,
    0.63283980608870629,
    0.61242994950949500,
    0.62018198080741576,
    0.61721461653440386,
    0.61834711965622806,
    0.61791440652881789,
    0.61807966846989581,
];

const PRAC_ADD_COST: f64 = 6.0;
const PRAC_DUP_COST: f64 = 5.0;

fn prac_lucas_cost(scalar: u64, ratio: f64) -> f64 {
    let mut rest = (scalar as f64 * ratio + 0.5) as u64;
    if rest >= scalar {
        return PRAC_ADD_COST * scalar as f64;
    }
    let mut first = scalar - rest;
    let mut second = 2 * rest - scalar;
    let mut cost = PRAC_DUP_COST + PRAC_ADD_COST;
    let mut guard = 0_u64;
    while first != second {
        guard += 1;
        assert!(guard < 1_000_000);
        if first < second {
            rest = first;
            first = second;
            second = rest;
        }
        if first - second <= second / 4 && (first + second).is_multiple_of(3) {
            first = (2 * first - second) / 3;
            second = (second - first) / 2;
            cost += 3.0 * PRAC_ADD_COST;
        } else if first - second <= second / 4 && (first - second).is_multiple_of(6) {
            first = (first - second) / 2;
            cost += PRAC_ADD_COST + PRAC_DUP_COST;
        } else if first.div_ceil(4) <= second {
            first -= second;
            cost += PRAC_ADD_COST;
        } else if (first + second).is_multiple_of(2) {
            first = (first - second) / 2;
            cost += PRAC_ADD_COST + PRAC_DUP_COST;
        } else if first.is_multiple_of(2) {
            first /= 2;
            cost += PRAC_ADD_COST + PRAC_DUP_COST;
        } else if first.is_multiple_of(3) {
            first = first / 3 - second;
            cost += 3.0 * PRAC_ADD_COST + PRAC_DUP_COST;
        } else if (first + second).is_multiple_of(3) {
            first = (first - 2 * second) / 3;
            cost += 3.0 * PRAC_ADD_COST + PRAC_DUP_COST;
        } else if (first - second).is_multiple_of(3) {
            first = (first - second) / 3;
            cost += 3.0 * PRAC_ADD_COST + PRAC_DUP_COST;
        } else {
            second /= 2;
            cost += PRAC_ADD_COST + PRAC_DUP_COST;
        }
    }
    cost
}

#[allow(clippy::too_many_arguments)]
fn prac_mult_into(
    field: &Montgomery,
    scratch: &mut LadderScratch,
    base_x: &[u64],
    base_z: &[u64],
    scalar: u64,
    curve_a24: &[u64],
    out_x: &mut [u64],
    out_z: &mut [u64],
    workspace: &mut EcmWorkspace,
) {
    debug_assert!(scalar > 2);
    let LadderScratch {
        low_x,
        low_z,
        high_x,
        high_z,
        doubled_x,
        doubled_z,
        added_x,
        added_z,
        next_x,
        next_z,
        spare_x,
        spare_z,
    } = scratch;
    let mut choice = 0_usize;
    let trials = field.limbs().min(PRAC_RATIOS.len());
    if trials > 1 {
        let mut best = PRAC_ADD_COST * scalar as f64;
        for (candidate, ratio) in PRAC_RATIOS.iter().take(trials).enumerate() {
            let cost = prac_lucas_cost(scalar, *ratio);
            if cost < best {
                best = cost;
                choice = candidate;
            }
        }
    }
    let split = (scalar as f64 * PRAC_RATIOS[choice] + 0.5) as u64;
    let mut first = scalar - split;
    let mut second = 2 * split - scalar;
    high_x.copy_from_slice(base_x);
    high_z.copy_from_slice(base_z);
    doubled_x.copy_from_slice(base_x);
    doubled_z.copy_from_slice(base_z);
    double_point_into(field, base_x, base_z, curve_a24, low_x, low_z, workspace);
    let mut guard = 0_u64;
    while first != second {
        guard += 1;
        assert!(guard < 1_000_000);
        if first < second {
            std::mem::swap(&mut first, &mut second);
            std::mem::swap(low_x, high_x);
            std::mem::swap(low_z, high_z);
        }
        if first - second <= second / 4 && (first + second).is_multiple_of(3) {
            first = (2 * first - second) / 3;
            second = (second - first) / 2;
            add_points_into(
                field, low_x, low_z, high_x, high_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            add_points_into(
                field, added_x, added_z, low_x, low_z, high_x, high_z, next_x, next_z, workspace,
            );
            add_points_into(
                field, high_x, high_z, added_x, added_z, low_x, low_z, spare_x, spare_z, workspace,
            );
            std::mem::swap(spare_x, high_x);
            std::mem::swap(spare_z, high_z);
            std::mem::swap(low_x, next_x);
            std::mem::swap(low_z, next_z);
        } else if first - second <= second / 4 && (first - second).is_multiple_of(6) {
            first = (first - second) / 2;
            add_points_into(
                field, low_x, low_z, high_x, high_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, high_x);
            std::mem::swap(spare_z, high_z);
            double_point_into(field, low_x, low_z, curve_a24, spare_x, spare_z, workspace);
            std::mem::swap(spare_x, low_x);
            std::mem::swap(spare_z, low_z);
        } else if first.div_ceil(4) <= second {
            first -= second;
            add_points_into(
                field, high_x, high_z, low_x, low_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            std::mem::swap(high_x, added_x);
            std::mem::swap(high_z, added_z);
            std::mem::swap(added_x, doubled_x);
            std::mem::swap(added_z, doubled_z);
        } else if (first + second).is_multiple_of(2) {
            first = (first - second) / 2;
            add_points_into(
                field, high_x, high_z, low_x, low_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, high_x);
            std::mem::swap(spare_z, high_z);
            double_point_into(field, low_x, low_z, curve_a24, spare_x, spare_z, workspace);
            std::mem::swap(spare_x, low_x);
            std::mem::swap(spare_z, low_z);
        } else if first.is_multiple_of(2) {
            first /= 2;
            add_points_into(
                field, doubled_x, doubled_z, low_x, low_z, high_x, high_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, doubled_x);
            std::mem::swap(spare_z, doubled_z);
            double_point_into(field, low_x, low_z, curve_a24, spare_x, spare_z, workspace);
            std::mem::swap(spare_x, low_x);
            std::mem::swap(spare_z, low_z);
        } else if first.is_multiple_of(3) {
            first = first / 3 - second;
            double_point_into(field, low_x, low_z, curve_a24, spare_x, spare_z, workspace);
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            add_points_into(
                field, low_x, low_z, high_x, high_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, next_x);
            std::mem::swap(spare_z, next_z);
            add_points_into(
                field, added_x, added_z, low_x, low_z, low_x, low_z, spare_x, spare_z, workspace,
            );
            std::mem::swap(spare_x, low_x);
            std::mem::swap(spare_z, low_z);
            add_points_into(
                field, added_x, added_z, next_x, next_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            std::mem::swap(doubled_x, high_x);
            std::mem::swap(doubled_z, high_z);
            std::mem::swap(high_x, added_x);
            std::mem::swap(high_z, added_z);
        } else if (first + second).is_multiple_of(3) {
            first = (first - 2 * second) / 3;
            add_points_into(
                field, low_x, low_z, high_x, high_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            add_points_into(
                field, added_x, added_z, low_x, low_z, high_x, high_z, spare_x, spare_z, workspace,
            );
            std::mem::swap(spare_x, high_x);
            std::mem::swap(spare_z, high_z);
            double_point_into(field, low_x, low_z, curve_a24, spare_x, spare_z, workspace);
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            add_points_into(
                field, low_x, low_z, added_x, added_z, low_x, low_z, spare_x, spare_z, workspace,
            );
            std::mem::swap(spare_x, low_x);
            std::mem::swap(spare_z, low_z);
        } else if (first - second).is_multiple_of(3) {
            first = (first - second) / 3;
            add_points_into(
                field, low_x, low_z, high_x, high_z, doubled_x, doubled_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            add_points_into(
                field, low_x, low_z, doubled_x, doubled_z, high_x, high_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, doubled_x);
            std::mem::swap(spare_z, doubled_z);
            std::mem::swap(high_x, added_x);
            std::mem::swap(high_z, added_z);
            double_point_into(field, low_x, low_z, curve_a24, spare_x, spare_z, workspace);
            std::mem::swap(spare_x, added_x);
            std::mem::swap(spare_z, added_z);
            add_points_into(
                field, low_x, low_z, added_x, added_z, low_x, low_z, spare_x, spare_z, workspace,
            );
            std::mem::swap(spare_x, low_x);
            std::mem::swap(spare_z, low_z);
        } else {
            second /= 2;
            add_points_into(
                field, doubled_x, doubled_z, high_x, high_z, low_x, low_z, spare_x, spare_z,
                workspace,
            );
            std::mem::swap(spare_x, doubled_x);
            std::mem::swap(spare_z, doubled_z);
            double_point_into(
                field, high_x, high_z, curve_a24, spare_x, spare_z, workspace,
            );
            std::mem::swap(spare_x, high_x);
            std::mem::swap(spare_z, high_z);
        }
    }
    add_points_into(
        field, low_x, low_z, high_x, high_z, doubled_x, doubled_z, spare_x, spare_z, workspace,
    );
    std::mem::swap(spare_x, low_x);
    std::mem::swap(spare_z, low_z);
    debug_assert_eq!(first, 1);
    out_x.copy_from_slice(low_x);
    out_z.copy_from_slice(low_z);
}

enum CurveSetup {
    Curve(Integer, Integer, Integer),
    Divisor(Integer),
    Next,
}

fn random_curve_setup(n: &Integer, counter: u64) -> CurveSetup {
    let point_x = Integer::from(counter + 3);
    let point_y = Integer::from(counter + 7);
    let square = (&point_x * &point_x).rem_euclid(n);
    let divisor = integer_gcd(&square, n);
    if divisor > Integer::from(1) && divisor < *n {
        return CurveSetup::Divisor(divisor);
    }
    if divisor == *n {
        return CurveSetup::Next;
    }
    let target = (&point_y * &point_y - &point_x * &point_x * &point_x - &point_x).rem_euclid(n);
    let inverse = match mod_inverse(&square, n) {
        Some(found) => found,
        None => return CurveSetup::Next,
    };
    let curve_a = (target * &inverse).rem_euclid(n);
    let inverse_four = match mod_inverse(&Integer::from(4), n) {
        Some(found) => found,
        None => return CurveSetup::Next,
    };
    let curve_a24 = (((curve_a + Integer::from(2)).rem_euclid(n)) * inverse_four).rem_euclid(n);
    CurveSetup::Curve(curve_a24, point_x, Integer::from(1))
}

fn scan_snapshots(
    field: &Montgomery,
    product: &[u64],
    snapshots: &[Vec<u64>],
    n: &Integer,
) -> Option<Integer> {
    if snapshots.is_empty() {
        return None;
    }
    let total = field.convert_out(product);
    if integer_gcd(&total, n) > Integer::from(1) {
        for snapshot in snapshots {
            let coordinate = field.convert_out(snapshot);
            if coordinate == Integer::from(0) {
                continue;
            }
            let divisor = integer_gcd(&coordinate, n);
            if divisor > Integer::from(1) && divisor < *n {
                return Some(divisor);
            }
        }
    }
    None
}

fn ecm_stage_one(
    field: &Montgomery,
    point_x: &[u64],
    point_z: &[u64],
    curve_a24: &[u64],
    n: &Integer,
    bound_one: u64,
    deadline: Option<Instant>,
) -> (Option<Integer>, Vec<u64>, Vec<u64>) {
    let mut current_x = point_x.to_vec();
    let mut current_z = point_z.to_vec();
    let mut workspace = ecm_workspace(field);
    let mut scratch = ladder_scratch(field);
    let one_repr = field.convert_in(&Integer::from(1));
    let mut product = one_repr.clone();
    let mut swap = vec![0_u64; field.limbs()];
    let mut next_x = vec![0_u64; field.limbs()];
    let mut next_z = vec![0_u64; field.limbs()];
    let mut snapshots: Vec<Vec<u64>> = Vec::new();
    let mut since_poll = 0_u32;
    for prime in crate::prime::sieve_primes(bound_one as usize) {
        let mut power = prime;
        let mut exponent = 1_u32;
        while power <= bound_one / prime {
            power *= prime;
            exponent += 1;
        }
        for _ in 0..exponent {
            if prime == 2 {
                double_point_into(
                    field,
                    &current_x,
                    &current_z,
                    curve_a24,
                    &mut next_x,
                    &mut next_z,
                    &mut workspace,
                );
            } else if prime == 3 {
                double_point_into(
                    field,
                    &current_x,
                    &current_z,
                    curve_a24,
                    &mut scratch.added_x,
                    &mut scratch.added_z,
                    &mut workspace,
                );
                let (left, right) = (&current_x as &[u64], &current_z as &[u64]);
                add_points_into(
                    field,
                    left,
                    right,
                    &scratch.added_x,
                    &scratch.added_z,
                    left,
                    right,
                    &mut next_x,
                    &mut next_z,
                    &mut workspace,
                );
            } else {
                prac_mult_into(
                    field,
                    &mut scratch,
                    &current_x,
                    &current_z,
                    prime,
                    curve_a24,
                    &mut next_x,
                    &mut next_z,
                    &mut workspace,
                );
            }
            std::mem::swap(&mut current_x, &mut next_x);
            std::mem::swap(&mut current_z, &mut next_z);
            if current_z.iter().all(|&word| word == 0) {
                break;
            }
        }
        if current_z.iter().all(|&word| word == 0) {
            break;
        }
        snapshots.push(current_z.clone());
        field.multiply_into(&mut swap, &product, &current_z, &mut workspace.middle);
        std::mem::swap(&mut product, &mut swap);
        since_poll += 1;
        if since_poll >= 128 {
            since_poll = 0;
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                break;
            }
            if let Some(divisor) = scan_snapshots(field, &product, &snapshots, n) {
                return (Some(divisor), current_x, current_z);
            }
            product.copy_from_slice(&one_repr);
            snapshots.clear();
        }
    }
    if let Some(divisor) = scan_snapshots(field, &product, &snapshots, n) {
        return (Some(divisor), current_x, current_z);
    }
    (None, current_x, current_z)
}

#[allow(clippy::too_many_arguments)]
fn ecm_stage_two(
    field: &Montgomery,
    point_x: &[u64],
    point_z: &[u64],
    curve_a24: &[u64],
    n: &Integer,
    bound_one: u64,
    bound_two: u64,
    prime_table: Option<&[bool]>,
    deadline: Option<Instant>,
) -> Option<Integer> {
    let bound_two = bound_two.min(bound_one.saturating_mul(16));
    let width = (bound_two as f64).sqrt() as u64 + 1;
    if bound_one < width || width < 2 {
        return None;
    }
    let mut workspace = ecm_workspace(field);
    let mut ladder = ladder_scratch(field);
    let limbs = field.limbs();
    let one_repr = field.convert_in(&Integer::from(1));
    let zero_repr = field.convert_in(&Integer::from(0));
    let mut babies_x = Vec::new();
    let mut babies_z = Vec::new();
    babies_x.push(one_repr.clone());
    babies_z.push(zero_repr);
    babies_x.push(point_x.to_vec());
    babies_z.push(point_z.to_vec());
    let mut second_x = vec![0_u64; limbs];
    let mut second_z = vec![0_u64; limbs];
    double_point_into(
        field,
        point_x,
        point_z,
        curve_a24,
        &mut second_x,
        &mut second_z,
        &mut workspace,
    );
    babies_x.push(second_x.clone());
    babies_z.push(second_z.clone());
    let mut previous_x = point_x.to_vec();
    let mut previous_z = point_z.to_vec();
    let mut current_x = second_x;
    let mut current_z = second_z;
    let mut step_x = vec![0_u64; limbs];
    let mut step_z = vec![0_u64; limbs];
    for _ in 3..width {
        add_points_into(
            field,
            &current_x,
            &current_z,
            point_x,
            point_z,
            &previous_x,
            &previous_z,
            &mut step_x,
            &mut step_z,
            &mut workspace,
        );
        babies_x.push(step_x.clone());
        babies_z.push(step_z.clone());
        std::mem::swap(&mut previous_x, &mut current_x);
        std::mem::swap(&mut current_x, &mut step_x);
        std::mem::swap(&mut previous_z, &mut current_z);
        std::mem::swap(&mut current_z, &mut step_z);
    }
    let mut stride_x = vec![0_u64; limbs];
    let mut stride_z = vec![0_u64; limbs];
    ladder_mult_into(
        field,
        &mut ladder,
        point_x,
        point_z,
        width,
        curve_a24,
        &one_repr,
        &mut stride_x,
        &mut stride_z,
        &mut workspace,
    );
    let mut previous_giant_x = vec![0_u64; limbs];
    let mut previous_giant_z = vec![0_u64; limbs];
    ladder_mult_into(
        field,
        &mut ladder,
        point_x,
        point_z,
        bound_one - width,
        curve_a24,
        &one_repr,
        &mut previous_giant_x,
        &mut previous_giant_z,
        &mut workspace,
    );
    let mut giant_x = vec![0_u64; limbs];
    let mut giant_z = vec![0_u64; limbs];
    ladder_mult_into(
        field,
        &mut ladder,
        point_x,
        point_z,
        bound_one,
        curve_a24,
        &one_repr,
        &mut giant_x,
        &mut giant_z,
        &mut workspace,
    );
    let mut giant_products = Vec::new();
    let mut accumulator = one_repr.clone();
    let mut batch = vec![0_u64; limbs];
    let mut tmp: [Vec<u64>; 4] = [(); 4].map(|_| vec![0_u64; limbs]);
    let [b0, b1, b2, b3] = &mut tmp;
    let mut scratch = vec![0_u64; field.multiply_scratch_len()];
    let mut index = bound_one;
    while index <= bound_two + width {
        batch.copy_from_slice(&one_repr);
        let mut position = 0_usize;
        while position < babies_x.len() {
            let baby = position as u64;
            let keep = match prime_table {
                None => true,
                Some(table) => {
                    let low = index.saturating_sub(baby);
                    let high = index.saturating_add(baby);
                    let low_prime = (low as usize) < table.len() && table[low as usize];
                    let high_prime = (high as usize) < table.len() && table[high as usize];
                    low_prime || high_prime
                }
            };
            if keep {
                field.multiply_into(b0, &giant_x, &babies_z[position], &mut scratch);
                field.multiply_into(b1, &babies_x[position], &giant_z, &mut scratch);
                field.subtract_into(b2, b0, b1);
                field.multiply_into(b3, &batch, b2, &mut scratch);
                std::mem::swap(&mut batch, b3);
            }
            position += 1;
        }
        field.multiply_into(b3, &accumulator, &batch, &mut scratch);
        std::mem::swap(&mut accumulator, b3);
        giant_products.push(accumulator.clone());
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            break;
        }
        add_points_into(
            field,
            &giant_x,
            &giant_z,
            &stride_x,
            &stride_z,
            &previous_giant_x,
            &previous_giant_z,
            &mut step_x,
            &mut step_z,
            &mut workspace,
        );
        std::mem::swap(&mut previous_giant_x, &mut giant_x);
        std::mem::swap(&mut giant_x, &mut step_x);
        std::mem::swap(&mut previous_giant_z, &mut giant_z);
        std::mem::swap(&mut giant_z, &mut step_z);
        index += width;
    }
    let divisor = integer_gcd(&field.convert_out(&accumulator), n);
    if divisor > Integer::from(1) && divisor < *n {
        return Some(divisor);
    }
    if divisor == *n && !giant_products.is_empty() {
        let mut low = 0_usize;
        let mut high = giant_products.len();
        while low < high {
            let middle = (low + high) / 2;
            if integer_gcd(&field.convert_out(&giant_products[middle]), n) == Integer::from(1) {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        if low < giant_products.len() {
            let scalar = bound_one + low as u64 * width;
            let mut guilty_x = vec![0_u64; limbs];
            let mut guilty_z = vec![0_u64; limbs];
            ladder_mult_into(
                field,
                &mut ladder,
                point_x,
                point_z,
                scalar,
                curve_a24,
                &one_repr,
                &mut guilty_x,
                &mut guilty_z,
                &mut workspace,
            );
            let mut position = 0_usize;
            while position < babies_x.len() {
                field.multiply_into(b0, &guilty_x, &babies_z[position], &mut scratch);
                field.multiply_into(b1, &babies_x[position], &guilty_z, &mut scratch);
                field.subtract_into(b2, b0, b1);
                let factor = integer_gcd(&field.convert_out(b2), n);
                if factor > Integer::from(1) && factor < *n {
                    return Some(factor);
                }
                position += 1;
            }
        }
    }
    None
}

fn ecm_single_curve(
    field: &Montgomery,
    n: &Integer,
    bound_one: u64,
    bound_two: u64,
    prime_table: Option<&[bool]>,
    sigma: u64,
    deadline: Option<Instant>,
) -> (Option<Integer>, u64, u64) {
    let no_timing = (None, 0, 0);
    let (curve_a24, point_x, point_z) = match random_curve_setup(n, sigma) {
        CurveSetup::Curve(a24, x, z) => (a24, x, z),
        CurveSetup::Divisor(divisor) => return (Some(divisor), 0, 0),
        CurveSetup::Next => return no_timing,
    };
    let base_a24 = field.convert_in(&curve_a24);
    let base_x = field.convert_in(&point_x);
    let base_z = field.convert_in(&point_z);
    let stage_one_start = Instant::now();
    let (found, end_x, end_z) =
        ecm_stage_one(field, &base_x, &base_z, &base_a24, n, bound_one, deadline);
    let stage_one_nanos = stage_one_start.elapsed().as_nanos() as u64;
    if found.is_some() {
        return (found, stage_one_nanos, 0);
    }
    if end_z.iter().all(|&word| word == 0) {
        return (None, stage_one_nanos, 0);
    }
    if bound_two <= bound_one || bound_two > STAGE_TWO_CAP {
        return (None, stage_one_nanos, 0);
    }
    let stage_two_start = Instant::now();
    let stage_two_found = ecm_stage_two(
        field,
        &end_x,
        &end_z,
        &base_a24,
        n,
        bound_one,
        bound_two,
        prime_table,
        deadline,
    );
    let stage_two_nanos = stage_two_start.elapsed().as_nanos() as u64;
    (stage_two_found, stage_one_nanos, stage_two_nanos)
}

pub(crate) struct CheapEcmReport {
    pub divisor: Option<Integer>,
    pub curves_run: u64,
    pub elapsed_ms: u128,
}

pub(crate) fn ecm_cheap_report(
    n: &Integer,
    bound_one: u64,
    bound_two: u64,
    max_curves: u64,
    deadline: Option<Instant>,
    workers: usize,
) -> CheapEcmReport {
    let started = Instant::now();
    let finish = |divisor: Option<Integer>, curves_run: u64| CheapEcmReport {
        divisor,
        curves_run,
        elapsed_ms: started.elapsed().as_millis(),
    };
    if n < &Integer::from(3) {
        return finish(None, 0);
    }
    if n % Integer::from(2) != Integer::from(1) {
        return finish(None, 0);
    }
    let field = match Montgomery::new(n) {
        Some(field) => field,
        None => return finish(None, 0),
    };
    let prime_table = stage_two_prime_table(bound_one, bound_two);
    let workers = workers.max(1);
    let mut curves_used = 0_u64;
    let mut next_curve = 1_u64;
    while next_curve <= max_curves {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return finish(None, curves_used);
        }
        let batch_end = (next_curve + workers as u64).min(max_curves + 1);
        let stop = AtomicBool::new(false);
        let found = Mutex::new(None);
        std::thread::scope(|scope| {
            let stop_flag = &stop;
            let found_slot = &found;
            let table_ref = prime_table.as_deref();
            let field_ref = &field;
            for sigma in next_curve..batch_end {
                let curve = sigma;
                scope.spawn(move || {
                    if stop_flag.load(Ordering::Relaxed) {
                        return;
                    }
                    let (hit, _, _) = ecm_single_curve(
                        field_ref, n, bound_one, bound_two, table_ref, curve, deadline,
                    );
                    if let Some(divisor) = hit {
                        if let Ok(mut slot) = found_slot.lock()
                            && slot.is_none()
                        {
                            *slot = Some(divisor);
                        }
                        stop_flag.store(true, Ordering::Relaxed);
                    }
                });
            }
        });
        curves_used += batch_end - next_curve;
        next_curve = batch_end;
        if let Some(divisor) = found
            .into_inner()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
        {
            return finish(Some(divisor), curves_used);
        }
    }
    finish(None, curves_used)
}

pub(crate) struct EcmReport {
    pub divisor: Option<Integer>,
    pub curves_run: u64,
    pub elapsed_ms: u128,
    pub stage_one_ms: u128,
    pub stage_two_ms: u128,
    pub end_digits: u32,
    pub end_curve: u64,
}

impl EcmReport {
    fn empty(
        divisor: Option<Integer>,
        curves_run: u64,
        started: Instant,
        end_digits: u32,
        end_curve: u64,
    ) -> Self {
        Self {
            divisor,
            curves_run,
            elapsed_ms: started.elapsed().as_millis(),
            stage_one_ms: 0,
            stage_two_ms: 0,
            end_digits,
            end_curve,
        }
    }
}

pub(crate) fn ecm_checkpoint(report: &EcmReport) -> (u32, u64) {
    (report.end_digits, report.end_curve)
}

pub(crate) fn ecm_factor_report(
    n: &Integer,
    target_digits: u32,
    max_curves: Option<u64>,
    deadline: Option<Instant>,
    workers: usize,
    resume: Option<(u32, u64)>,
) -> EcmReport {
    let started = Instant::now();
    let floor_digits = target_digits.max(1);
    if n < &Integer::from(3) {
        return EcmReport::empty(None, 0, started, floor_digits, 1);
    }
    if n % Integer::from(2) != Integer::from(1) {
        return EcmReport::empty(None, 0, started, floor_digits, 1);
    }
    let field = match Montgomery::new(n) {
        Some(field) => field,
        None => return EcmReport::empty(None, 0, started, floor_digits, 1),
    };
    let mut digits = match resume {
        Some((resume_digits, _)) => resume_digits.max(floor_digits),
        None => floor_digits,
    };
    let mut resume_curve = resume.map(|(_, curve)| curve.max(1));
    let mut curves_used = 0_u64;
    let stage_one_nanos = std::sync::atomic::AtomicU64::new(0);
    let stage_two_nanos = std::sync::atomic::AtomicU64::new(0);
    let workers = workers.max(1);
    loop {
        let (bound_one, bound_two, expected) = ecm_plan(digits);
        let budgeted = expected.saturating_mul(TORSION_COMPENSATION);
        let prime_table = stage_two_prime_table(bound_one, bound_two);
        let mut next_curve = resume_curve.take().unwrap_or(1);
        while next_curve <= budgeted {
            let allowance = match max_curves {
                Some(limit) => limit.saturating_sub(curves_used),
                None => u64::MAX,
            };
            if allowance == 0 {
                return EcmReport {
                    divisor: None,
                    curves_run: curves_used,
                    elapsed_ms: started.elapsed().as_millis(),
                    stage_one_ms: u128::from(stage_one_nanos.load(Ordering::Relaxed)) / 1_000_000,
                    stage_two_ms: u128::from(stage_two_nanos.load(Ordering::Relaxed)) / 1_000_000,
                    end_digits: digits,
                    end_curve: next_curve,
                };
            }
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                return EcmReport {
                    divisor: None,
                    curves_run: curves_used,
                    elapsed_ms: started.elapsed().as_millis(),
                    stage_one_ms: u128::from(stage_one_nanos.load(Ordering::Relaxed)) / 1_000_000,
                    stage_two_ms: u128::from(stage_two_nanos.load(Ordering::Relaxed)) / 1_000_000,
                    end_digits: digits,
                    end_curve: next_curve,
                };
            }
            let batch_end = (next_curve + (workers as u64).min(allowance)).min(budgeted + 1);
            let stop = AtomicBool::new(false);
            let found = Mutex::new(None);
            std::thread::scope(|scope| {
                let stop_flag = &stop;
                let found_slot = &found;
                let table_ref = prime_table.as_deref();
                let stage_one_acc = &stage_one_nanos;
                let stage_two_acc = &stage_two_nanos;
                let field_ref = &field;
                for sigma in next_curve..batch_end {
                    let curve = sigma;
                    scope.spawn(move || {
                        if stop_flag.load(Ordering::Relaxed) {
                            return;
                        }
                        let (hit, stage_one, stage_two) = ecm_single_curve(
                            field_ref, n, bound_one, bound_two, table_ref, curve, deadline,
                        );
                        stage_one_acc.fetch_add(stage_one, Ordering::Relaxed);
                        stage_two_acc.fetch_add(stage_two, Ordering::Relaxed);
                        if let Some(divisor) = hit {
                            if let Ok(mut slot) = found_slot.lock()
                                && slot.is_none()
                            {
                                *slot = Some(divisor);
                            }
                            stop_flag.store(true, Ordering::Relaxed);
                        }
                    });
                }
            });
            curves_used += batch_end - next_curve;
            next_curve = batch_end;
            if let Some(divisor) = found
                .into_inner()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
            {
                return EcmReport {
                    divisor: Some(divisor),
                    curves_run: curves_used,
                    elapsed_ms: started.elapsed().as_millis(),
                    stage_one_ms: u128::from(stage_one_nanos.load(Ordering::Relaxed)) / 1_000_000,
                    stage_two_ms: u128::from(stage_two_nanos.load(Ordering::Relaxed)) / 1_000_000,
                    end_digits: digits,
                    end_curve: next_curve,
                };
            }
        }
        if deadline.is_none() || digits >= MAX_DIGITS {
            return EcmReport {
                divisor: None,
                curves_run: curves_used,
                elapsed_ms: started.elapsed().as_millis(),
                stage_one_ms: u128::from(stage_one_nanos.load(Ordering::Relaxed)) / 1_000_000,
                stage_two_ms: u128::from(stage_two_nanos.load(Ordering::Relaxed)) / 1_000_000,
                end_digits: digits,
                end_curve: next_curve,
            };
        }
        digits = escalate_digits(digits);
    }
}

#[cfg(test)]
mod ecm_tests {
    use super::*;

    #[test]
    fn prime_table_marks_primes_only() {
        let table = prime_table(30);
        assert_eq!(table.len(), 31);
        for (value, &flag) in table.iter().enumerate() {
            let expected = [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29].contains(&(value as u64));
            assert_eq!(flag, expected, "value={value}");
        }
    }

    #[test]
    fn inverse_known_values() {
        let modulus = Integer::from(11);
        assert_eq!(
            mod_inverse(&Integer::from(3), &modulus),
            Some(Integer::from(4))
        );
        assert_eq!(
            mod_inverse(&Integer::from(10), &Integer::from(17)),
            Some(Integer::from(12))
        );
        assert_eq!(mod_inverse(&Integer::from(4), &Integer::from(6)), None);
    }

    #[test]
    fn ladder_matches_formulas() {
        let modulus = Integer::from(101);
        let field = Montgomery::new(&modulus).expect("odd modulus");
        let two = field.convert_in(&Integer::from(2));
        let one = field.convert_in(&Integer::from(1));
        let three = field.convert_in(&Integer::from(3));
        let doubled = double_point(&field, &two, &one, &two);
        assert_eq!(field.convert_out(&doubled.0), Integer::from(9));
        assert_eq!(field.convert_out(&doubled.1), Integer::from(35));
        let added = add_points(&field, &two, &one, &three, &one, &one, &one);
        assert_eq!(field.convert_out(&added.0), Integer::from(100));
        assert_eq!(field.convert_out(&added.1), Integer::from(4));
        let direct = ladder_mult(&field, &two, &one, 6, &two);
        let inner = ladder_mult(&field, &two, &one, 3, &two);
        let composed = ladder_mult(&field, &inner.0, &inner.1, 2, &two);
        let left = field.multiply(&direct.0, &composed.1);
        let right = field.multiply(&direct.1, &composed.0);
        assert_eq!(left, right);
    }

    #[test]
    fn prac_matches_ladder_small_field() {
        let modulus = Integer::from(101);
        let field = Montgomery::new(&modulus).expect("odd modulus");
        let two = field.convert_in(&Integer::from(2));
        let one = field.convert_in(&Integer::from(1));
        for prime in crate::prime::sieve_primes(50_000) {
            if prime < 5 {
                continue;
            }
            let mut workspace = ecm_workspace(&field);
            let mut scratch = ladder_scratch(&field);
            let mut prac_x = vec![0_u64; field.limbs()];
            let mut prac_z = vec![0_u64; field.limbs()];
            prac_mult_into(
                &field,
                &mut scratch,
                &two,
                &one,
                prime,
                &two,
                &mut prac_x,
                &mut prac_z,
                &mut workspace,
            );
            let (ladder_x, ladder_z) = ladder_mult(&field, &two, &one, prime, &two);
            let left = field.multiply(&prac_x, &ladder_z);
            let right = field.multiply(&ladder_x, &prac_z);
            assert_eq!(left, right, "prime={prime}");
        }
    }

    #[test]
    fn prac_matches_ladder_large_field() {
        let big = (Integer::from(1) << 267_u32) - Integer::from(441);
        let big = if big.clone() % Integer::from(2) == Integer::from(1) {
            big
        } else {
            big + Integer::from(1)
        };
        let field = Montgomery::new(&big).expect("odd modulus");
        let base = field.convert_in(&Integer::from(5));
        let one = field.convert_in(&Integer::from(1));
        let curve = field.convert_in(&Integer::from(7));
        let mut scalars: Vec<u64> = crate::prime::sieve_primes(50_000)
            .into_iter()
            .filter(|prime| *prime >= 5)
            .collect();
        scalars.extend([
            1_000_003,
            1_000_033,
            1_000_037,
            2_147_483_647,
            4_294_967_291,
        ]);
        for scalar in scalars {
            let mut workspace = ecm_workspace(&field);
            let mut scratch = ladder_scratch(&field);
            let mut prac_x = vec![0_u64; field.limbs()];
            let mut prac_z = vec![0_u64; field.limbs()];
            prac_mult_into(
                &field,
                &mut scratch,
                &base,
                &one,
                scalar,
                &curve,
                &mut prac_x,
                &mut prac_z,
                &mut workspace,
            );
            let (ladder_x, ladder_z) = ladder_mult(&field, &base, &one, scalar, &curve);
            let left = field.multiply(&prac_x, &ladder_z);
            let right = field.multiply(&ladder_x, &prac_z);
            assert_eq!(left, right, "scalar={scalar}");
        }
    }

    #[test]
    fn prac_chain_cost_below_ladder() {
        let mut ladder_total = 0.0_f64;
        let mut prac_total = 0.0_f64;
        for prime in crate::prime::sieve_primes(11000) {
            let mut power = prime;
            let mut exponent = 1_u32;
            while power <= 11000 / prime {
                power *= prime;
                exponent += 1;
            }
            if prime == 2 {
                ladder_total += 5.0 * f64::from(exponent);
                prac_total += 5.0 * f64::from(exponent);
            } else if prime == 3 {
                ladder_total += 11.0 * f64::from(64 - power.leading_zeros());
                prac_total += 11.0 * f64::from(exponent);
            } else {
                let width = 64 - power.leading_zeros();
                ladder_total += 11.0 * width as f64;
                let mut best = 6.0 * prime as f64;
                for ratio in PRAC_RATIOS {
                    best = best.min(prac_lucas_cost(prime, ratio));
                }
                prac_total += best * f64::from(exponent);
            }
        }
        assert!(
            prac_total < 0.8 * ladder_total,
            "prac={prac_total} ladder={ladder_total}"
        );
    }

    #[test]
    fn stage_two_runs_at_raised_cap() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(300));
        let target: Integer = "9671406556917033397649483".parse().expect("prime decimal");
        let field = Montgomery::new(&target).expect("odd modulus");
        let table = stage_two_prime_table(50_000, 13_000_000).expect("d25 table");
        let (found, _, stage_two) = ecm_single_curve(
            &field,
            &target,
            50_000,
            13_000_000,
            Some(&table),
            1,
            deadline,
        );
        assert!(found.is_none());
        assert!(
            stage_two > 0,
            "stage two must execute at B2=13M, skip path taken"
        );
    }

    #[test]
    fn ecm_finds_small_factor() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(300));
        let target = Integer::from(1000003) * Integer::from(1000033);
        let found = ecm_factor_report(&target, 12, Some(30), deadline, 2, None).divisor;
        let divisor = found.expect("ecm finds a factor");
        assert!(divisor > Integer::from(1));
        assert!(divisor < target);
        assert!(target % &divisor == Integer::from(0));
    }

    #[test]
    fn resume_continues_past_tranche_boundary() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(300));
        let target = crate::factor::parse_decimal("3994424150216394531668178182146365185669")
            .expect("battery input");
        let full = ecm_factor_report(&target, 12, Some(222), deadline, 1, None);
        let divisor = full.divisor.clone().expect("fresh run wins");
        assert_eq!(full.curves_run, 16);
        assert!(&target % &divisor == Integer::from(0));
        let partial = ecm_factor_report(&target, 12, Some(15), deadline, 1, None);
        assert!(partial.divisor.is_none());
        assert_eq!((partial.end_digits, partial.end_curve), (12, 16));
        let resumed = ecm_factor_report(
            &target,
            12,
            None,
            deadline,
            1,
            Some((partial.end_digits, partial.end_curve)),
        );
        assert_eq!(resumed.divisor, Some(divisor));
        assert_eq!(resumed.curves_run, 1);
        assert_eq!(full.curves_run, partial.curves_run + resumed.curves_run);
    }

    #[test]
    fn resume_extends_miss_accounting() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(300));
        let prime = Integer::from(10007);
        let first = ecm_factor_report(&prime, 12, Some(10), deadline, 1, None);
        assert!(first.divisor.is_none());
        assert_eq!(first.curves_run, 10);
        let second = ecm_factor_report(
            &prime,
            12,
            Some(10),
            deadline,
            1,
            Some((first.end_digits, first.end_curve)),
        );
        assert!(second.divisor.is_none());
        assert_eq!(second.curves_run, 10);
        let fresh = ecm_factor_report(&prime, 12, Some(20), deadline, 1, None);
        assert_eq!(fresh.curves_run, 20);
    }

    #[test]
    fn ecm_multi_curve_finds_factor() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(300));
        let target = Integer::from(856679421607u64) * Integer::from(114019471u64);
        let field = Montgomery::new(&target).expect("odd modulus");
        let table = stage_two_prime_table(1000, 16000).expect("small table");
        let mut hits = 0_u32;
        for sigma in 1_u64..=8 {
            if ecm_single_curve(&field, &target, 1000, 16000, Some(&table), sigma, deadline)
                .0
                .is_some()
            {
                hits += 1;
            }
        }
        assert!(hits >= 1);
    }

    #[test]
    fn cheap_tier_finds_small_factor() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(300));
        let target = Integer::from(1000003) * Integer::from(1000033);
        let report = ecm_cheap_report(&target, 2_000, 200_000, 32, deadline, 2);
        let divisor = report.divisor.expect("cheap tier finds a factor");
        assert!(divisor > Integer::from(1));
        assert!(divisor < target);
        assert!(target % &divisor == Integer::from(0));
        assert!(report.curves_run > 0);
        assert!(report.curves_run <= 32);
    }

    #[test]
    fn ecm_prime_returns_none() {
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(120));
        assert_eq!(
            ecm_factor_report(&Integer::from(10007), 12, Some(10), deadline, 2, None).divisor,
            None
        );
    }

    #[test]
    fn random_curve_point_on_curve() {
        use num_traits::ToPrimitive;
        let field = 101_u64;
        for counter in 0_u64..8 {
            let modulus = Integer::from(field);
            let (curve_a24, point_x) = match random_curve_setup(&modulus, counter) {
                CurveSetup::Curve(a24, x, _) => (a24, x),
                _ => panic!("counter={counter} degenerate"),
            };
            let a = (4 * (curve_a24.to_u64().expect("small") % field) + field - 2) % field;
            let x = point_x.to_u64().expect("small") % field;
            let rhs = (x * x % field * x % field + a * x % field * x % field + x) % field;
            let symbol = crate::prime::legendre(rhs, field);
            assert!(symbol >= 0, "counter={counter} point off curve");
        }
    }
}
