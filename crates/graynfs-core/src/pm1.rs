use crate::integer::{Integer, integer_gcd};
use crate::montgomery::Montgomery;
use num_traits::Euclid;
use std::time::Instant;

pub(crate) const STAGE_ONE_BOUND: u64 = 100_000;

pub(crate) fn pm1_bound_for(target_digits: u32) -> u64 {
    if target_digits <= 25 {
        return STAGE_ONE_BOUND;
    }
    let mut steps = target_digits.saturating_sub(25).div_ceil(5);
    let mut bound = STAGE_ONE_BOUND;
    while steps > 0 {
        bound = bound.saturating_mul(2);
        steps -= 1;
    }
    bound
}

pub(crate) const PM1_PRIMES_PER_SEC: u64 = 200_000;

pub(crate) fn bound_for_count(count: u64) -> u64 {
    if count < 6 {
        return [2, 3, 5, 7, 11, 13][count as usize];
    }
    let target = count as f64;
    let mut bound = target * target.ln().max(1.0);
    for _ in 0..8 {
        bound = target * bound.ln().max(1.0);
    }
    bound.max(2.0) as u64
}

pub(crate) fn pm1_bound_for_slice(
    target_digits: u32,
    slice_ms: Option<u128>,
    explore: bool,
    trial_full_scan: bool,
) -> u64 {
    let scaled = pm1_bound_for(target_digits);
    if !explore || !trial_full_scan {
        return STAGE_ONE_BOUND;
    }
    let Some(slice) = slice_ms else {
        return STAGE_ONE_BOUND;
    };
    let affordable = slice.saturating_mul(u128::from(PM1_PRIMES_PER_SEC)) / 1000 / 2;
    let capped = bound_for_count(affordable.min(u64::MAX as u128) as u64);
    scaled.min(capped.max(STAGE_ONE_BOUND))
}

fn smooth_ladder(bound: u64) -> Integer {
    let mut product = Integer::from(1);
    for prime in crate::prime::sieve_primes(bound as usize) {
        let mut power = prime;
        while power <= bound / prime {
            power *= prime;
        }
        product *= Integer::from(power);
    }
    product
}

fn jacobi_symbol(numerator: &Integer, denominator: &Integer) -> i32 {
    if denominator < &Integer::from(1) {
        return 0;
    }
    let mut top = numerator.rem_euclid(denominator);
    let mut bottom = denominator.clone();
    let mut sign = 1_i32;
    let zero = Integer::from(0);
    let one = Integer::from(1);
    let two = Integer::from(2);
    loop {
        if bottom == one {
            return sign;
        }
        if top == zero {
            return 0;
        }
        while {
            let half = &top / &two;
            &half * &two == top
        } {
            top = &top / &two;
            let residue = &bottom % Integer::from(8);
            if residue == Integer::from(3) || residue == Integer::from(5) {
                sign = -sign;
            }
        }
        std::mem::swap(&mut top, &mut bottom);
        let top_mod = &top % Integer::from(4);
        let bottom_mod = &bottom % Integer::from(4);
        if top_mod == Integer::from(3) && bottom_mod == Integer::from(3) {
            sign = -sign;
        }
        top = top.rem_euclid(&bottom);
    }
}

fn lucas_chain_mont(
    field: &Montgomery,
    bytes: &[u8],
    width: usize,
    base_repr: &[u64],
    two_repr: &[u64],
) -> Vec<u64> {
    let limbs = field.limbs();
    let mut low = two_repr.to_vec();
    let mut high = base_repr.to_vec();
    let mut doubled_low = vec![0_u64; limbs];
    let mut doubled_high = vec![0_u64; limbs];
    let mut advanced = vec![0_u64; limbs];
    let mut product = vec![0_u64; limbs];
    let mut middle = vec![0_u64; field.multiply_scratch_len()];
    for position in (0..width).rev() {
        let byte = bytes[bytes.len() - 1 - position / 8];
        field.multiply_into(&mut product, &low, &low, &mut middle);
        field.subtract_into(&mut doubled_low, &product, two_repr);
        field.multiply_into(&mut product, &low, &high, &mut middle);
        field.subtract_into(&mut doubled_high, &product, base_repr);
        if (byte >> (position % 8)) & 1 == 0 {
            std::mem::swap(&mut low, &mut doubled_low);
            std::mem::swap(&mut high, &mut doubled_high);
        } else {
            field.multiply_into(&mut product, base_repr, &doubled_high, &mut middle);
            field.subtract_into(&mut advanced, &product, &doubled_low);
            std::mem::swap(&mut low, &mut doubled_high);
            std::mem::swap(&mut high, &mut advanced);
        }
    }
    low
}

fn scan_snapshots_pm(
    field: &Montgomery,
    product: &[u64],
    snapshots: &[Vec<u64>],
    one_repr: &[u64],
    offset: &mut [u64],
    n: &Integer,
) -> Option<Integer> {
    if snapshots.is_empty() {
        return None;
    }
    let total = field.convert_out(product);
    if integer_gcd(&total, n) > Integer::from(1) {
        for snapshot in snapshots {
            field.subtract_into(offset, snapshot, one_repr);
            let coordinate = field.convert_out(offset);
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

pub(crate) fn pollard_pm1(n: &Integer, bound: u64, deadline: Option<Instant>) -> Option<Integer> {
    if n < &Integer::from(3) {
        return None;
    }
    if n % Integer::from(2) != Integer::from(1) {
        return None;
    }
    let field = Montgomery::new(n)?;
    let one_repr = field.convert_in(&Integer::from(1));
    let mut acc = field.convert_in(&Integer::from(2));
    let mut product = one_repr.clone();
    let mut snapshots: Vec<Vec<u64>> = Vec::new();
    let mut middle = vec![0_u64; field.multiply_scratch_len()];
    let mut swap = vec![0_u64; field.limbs()];
    let mut result = vec![0_u64; field.limbs()];
    let mut factor = vec![0_u64; field.limbs()];
    let mut offset = vec![0_u64; field.limbs()];
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return None;
    }
    let mut since_poll = 0_u32;
    let mut found: Option<Integer> = None;
    let mut saturated = false;
    crate::prime::for_each_prime(bound, deadline, |prime| {
        let mut power = prime;
        while power <= bound / prime {
            power *= prime;
        }
        result.copy_from_slice(&one_repr);
        factor.copy_from_slice(&acc);
        let mut exponent = power;
        while exponent > 0 {
            if exponent & 1 == 1 {
                field.multiply_into(&mut swap, &result, &factor, &mut middle);
                std::mem::swap(&mut result, &mut swap);
            }
            field.multiply_into(&mut swap, &factor, &factor, &mut middle);
            std::mem::swap(&mut factor, &mut swap);
            exponent >>= 1;
        }
        acc.copy_from_slice(&result);
        if acc == one_repr {
            if let Some(divisor) =
                scan_snapshots_pm(&field, &product, &snapshots, &one_repr, &mut offset, n)
            {
                found = Some(divisor);
            }
            saturated = true;
            return false;
        }
        snapshots.push(acc.clone());
        field.subtract_into(&mut offset, &acc, &one_repr);
        field.multiply_into(&mut swap, &product, &offset, &mut middle);
        std::mem::swap(&mut product, &mut swap);
        since_poll += 1;
        if since_poll >= 128 {
            since_poll = 0;
            if let Some(divisor) =
                scan_snapshots_pm(&field, &product, &snapshots, &one_repr, &mut offset, n)
            {
                found = Some(divisor);
                return false;
            }
            product.copy_from_slice(&one_repr);
            snapshots.clear();
        }
        true
    });
    if saturated || found.is_some() {
        return found;
    }
    scan_snapshots_pm(&field, &product, &snapshots, &one_repr, &mut offset, n)
}

pub(crate) fn williams_pp1(n: &Integer, bound: u64, deadline: Option<Instant>) -> Option<Integer> {
    if n < &Integer::from(3) {
        return None;
    }
    if n % Integer::from(2) != Integer::from(1) {
        return None;
    }
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return None;
    }
    let field = Montgomery::new(n)?;
    let ladder = smooth_ladder(bound);
    let (_, bytes) = ladder.to_bytes_be();
    let width = ladder.bits() as usize;
    let two_repr = field.convert_in(&Integer::from(2));
    let mut offset = vec![0_u64; field.limbs()];
    let mut parameter = 3_u64;
    let mut used = 0_u32;
    let mut attempts = 0_u32;
    while used < 3 && attempts < 50 {
        attempts += 1;
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return None;
        }
        let discriminant = Integer::from(parameter * parameter - 4);
        match jacobi_symbol(&discriminant, n) {
            0 => {
                let divisor = integer_gcd(&discriminant, n);
                if divisor > Integer::from(1) && divisor < *n {
                    return Some(divisor);
                }
            }
            -1 => {
                used += 1;
                let base_repr = field.convert_in(&Integer::from(parameter));
                let root = lucas_chain_mont(&field, &bytes, width, &base_repr, &two_repr);
                field.subtract_into(&mut offset, &root, &two_repr);
                let coordinate = field.convert_out(&offset);
                if coordinate != Integer::from(0) {
                    let divisor = integer_gcd(&coordinate, n);
                    if divisor > Integer::from(1) && divisor < *n {
                        return Some(divisor);
                    }
                }
            }
            _ => {}
        }
        parameter += 1;
    }
    None
}

#[cfg(test)]
mod pm1_tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn jacobi_known_values() {
        assert_eq!(jacobi_symbol(&Integer::from(2), &Integer::from(7)), 1);
        assert_eq!(jacobi_symbol(&Integer::from(2), &Integer::from(3)), -1);
        assert_eq!(jacobi_symbol(&Integer::from(6), &Integer::from(35)), -1);
        assert_eq!(jacobi_symbol(&Integer::from(4), &Integer::from(35)), 1);
        assert_eq!(jacobi_symbol(&Integer::from(5), &Integer::from(5)), 0);
    }

    #[test]
    fn lucas_known_values() {
        let modulus = Integer::from(1_000_000_007);
        let field = Montgomery::new(&modulus).expect("odd modulus");
        let two_repr = field.convert_in(&Integer::from(2));
        let evaluate = |index: u64, base: u64, expected: u64| {
            let index = Integer::from(index);
            let (_, bytes) = index.to_bytes_be();
            let width = index.bits() as usize;
            let base_repr = field.convert_in(&Integer::from(base));
            let low = lucas_chain_mont(&field, &bytes, width, &base_repr, &two_repr);
            assert_eq!(field.convert_out(&low), Integer::from(expected));
        };
        evaluate(5, 3, 123);
        evaluate(0, 3, 2);
        evaluate(1, 7, 7);
        evaluate(10, 3, 15127);
    }

    #[test]
    fn pm1_smooth_factor() {
        let deadline = Some(Instant::now() + Duration::from_secs(120));
        let target = Integer::from(29) * Integer::from(10007);
        assert_eq!(target, Integer::from(290203));
        let found = pollard_pm1(&target, 10, deadline);
        assert_eq!(found, Some(Integer::from(29)));
        assert_eq!(pollard_pm1(&Integer::from(10007), 10, deadline), None);
    }

    #[test]
    fn slice_bound_prefers_scaled_when_affordable() {
        assert_eq!(
            pm1_bound_for_slice(18, Some(60_000), true, true),
            STAGE_ONE_BOUND
        );
        assert_eq!(pm1_bound_for_slice(30, Some(600_000), true, true), 200_000);
        let wide = pm1_bound_for_slice(80, Some(60_000), true, true);
        assert!(wide > STAGE_ONE_BOUND);
        assert!(wide <= pm1_bound_for(80));
    }

    #[test]
    fn slice_bound_falls_back_on_tiny_or_missing_slice() {
        let tiny = pm1_bound_for_slice(80, Some(100), true, true);
        assert!(tiny <= STAGE_ONE_BOUND * 2, "no explosion, got {tiny}");
        assert_eq!(pm1_bound_for_slice(80, None, true, true), STAGE_ONE_BOUND);
    }

    #[test]
    fn slice_bound_stays_default_without_explore_flag() {
        assert_eq!(
            pm1_bound_for_slice(80, Some(600_000), false, true),
            STAGE_ONE_BOUND
        );
        assert_eq!(
            pm1_bound_for_slice(30, Some(600_000), false, true),
            STAGE_ONE_BOUND
        );
    }

    #[test]
    fn slice_bound_needs_full_scan_for_extension() {
        assert_eq!(
            pm1_bound_for_slice(80, Some(600_000), true, false),
            STAGE_ONE_BOUND
        );
        assert_eq!(
            pm1_bound_for_slice(30, Some(600_000), true, false),
            STAGE_ONE_BOUND
        );
    }

    #[test]
    fn bound_for_count_inverts_prime_counts() {
        let near = bound_for_count(9592);
        let ratio = near as f64 / 100_000.0;
        assert!(ratio > 0.8 && ratio < 1.2, "near 100k, got {near}");
    }

    #[test]
    fn expired_deadline_huge_bound_returns_fast() {
        let past = Some(Instant::now());
        std::thread::sleep(Duration::from_millis(5));
        let target = Integer::from(14699) * Integer::from(80831689);
        let started = Instant::now();
        assert_eq!(pollard_pm1(&target, 100_000_000, past), None);
        assert_eq!(williams_pp1(&target, 100_000_000, past), None);
        assert!(
            started.elapsed() < Duration::from_secs(10),
            "expired deadline must short-circuit before sieving"
        );
    }

    #[test]
    fn pp1_smooth_factor() {
        let deadline = Some(Instant::now() + Duration::from_secs(120));
        let target = Integer::from(14699) * Integer::from(80831689);
        assert_eq!(target, Integer::from(1188144996611u64));
        assert_eq!(pollard_pm1(&target, 100, deadline), None);
        let found = williams_pp1(&target, 100, deadline);
        let divisor = found.expect("pp1 finds a divisor");
        assert!(divisor > Integer::from(1));
        assert!(divisor < target);
        assert!(target % &divisor == Integer::from(0));
    }
}
