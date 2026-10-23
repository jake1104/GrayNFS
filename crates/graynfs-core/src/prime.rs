use crate::integer::Integer;
use num_traits::{Euclid, ToPrimitive, Zero};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

pub const DETERMINISTIC_U64_BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

pub fn is_prime_u64(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    for prime in DETERMINISTIC_U64_BASES {
        if value == prime {
            return true;
        }
        if value.is_multiple_of(prime) {
            return false;
        }
    }
    let one = value - 1;
    let shift = one.trailing_zeros();
    let odd = one >> shift;
    'bases: for base in DETERMINISTIC_U64_BASES {
        let mut root = mod_pow(base, odd, value);
        if root == 1 || root == one {
            continue;
        }
        for _ in 1..shift {
            root = mul_mod(root, root, value);
            if root == one {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

pub fn small_primes() -> &'static [u64] {
    static TABLE: OnceLock<Vec<u64>> = OnceLock::new();
    TABLE.get_or_init(|| sieve_primes(10_000))
}

pub fn sieve_primes(limit: usize) -> Vec<u64> {
    if limit < 2 {
        return Vec::new();
    }
    let mut sieve = vec![true; limit + 1];
    sieve[0] = false;
    sieve[1] = false;
    let mut value = 2_usize;
    while value * value <= limit {
        if sieve[value] {
            let mut multiple = value * value;
            while multiple <= limit {
                sieve[multiple] = false;
                multiple += value;
            }
        }
        value += 1;
    }
    sieve
        .iter()
        .enumerate()
        .filter_map(|(prime, present)| present.then_some(prime as u64))
        .collect()
}

pub fn for_each_prime(
    limit: u64,
    deadline: Option<std::time::Instant>,
    mut visit: impl FnMut(u64) -> bool,
) -> bool {
    if limit < 2 {
        return true;
    }
    if deadline.is_some_and(|reached| std::time::Instant::now() >= reached) {
        return false;
    }
    let root = (limit as f64).sqrt() as usize + 1;
    let base = sieve_primes(root);
    const SEGMENT: u64 = 1_000_000;
    let mut low = 2_u64;
    while low <= limit {
        if deadline.is_some_and(|reached| std::time::Instant::now() >= reached) {
            return false;
        }
        let high = (low + SEGMENT - 1).min(limit);
        let mut sieve = vec![true; (high - low + 1) as usize];
        for &prime in &base {
            let square = prime * prime;
            if square > high {
                break;
            }
            let mut mark = low.div_ceil(prime) * prime;
            if mark < square {
                mark = square;
            }
            while mark <= high {
                sieve[(mark - low) as usize] = false;
                mark += prime;
            }
        }
        for (offset, present) in sieve.iter().enumerate() {
            if *present && !visit(low + offset as u64) {
                return true;
            }
        }
        low = high + 1;
    }
    true
}

pub fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

pub fn mod_pow(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    base %= modulus;
    let mut result = 1_u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod(result, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    result
}

pub fn legendre(value: u64, prime: u64) -> i32 {
    if prime == 2 {
        return if value & 1 == 0 { 1 } else { 0 };
    }
    let result = mod_pow(value % prime, (prime - 1) / 2, prime);
    if result == 1 {
        1
    } else if result == prime - 1 {
        -1
    } else {
        0
    }
}

pub fn mod_sqrt(value: u64, prime: u64) -> Option<u64> {
    if prime == 2 {
        return Some(value & 1);
    }
    let value = value % prime;
    if value == 0 {
        return Some(0);
    }
    if legendre(value, prime) != 1 {
        return None;
    }
    if prime % 4 == 3 {
        return Some(mod_pow(value, (prime + 1) / 4, prime));
    }
    let mut q = prime - 1;
    let mut s = 0_u32;
    while q & 1 == 0 {
        q >>= 1;
        s += 1;
    }
    let mut z = 2_u64;
    while legendre(z, prime) != -1 {
        z += 1;
    }
    let mut c = mod_pow(z, q, prime);
    let mut t = mod_pow(value, q, prime);
    let mut r = mod_pow(value, q.div_ceil(2), prime);
    let mut m = s;
    while t != 1 {
        let mut i = 1_u32;
        let mut t2 = mul_mod(t, t, prime);
        while i < m && t2 != 1 {
            t2 = mul_mod(t2, t2, prime);
            i += 1;
        }
        if i == m {
            return None;
        }
        let exponent = 1_u64 << (m - i - 1);
        let b = mod_pow(c, exponent, prime);
        r = mul_mod(r, b, prime);
        c = mul_mod(b, b, prime);
        t = mul_mod(t, c, prime);
        m = i;
    }
    Some(r)
}

pub fn integer_sqrt(value: &Integer) -> Integer {
    value.sqrt()
}

pub fn is_integer_square(value: &Integer) -> bool {
    if value < &Integer::from(0) {
        return false;
    }
    let root = integer_sqrt(value);
    &root * &root == *value
}

pub fn is_probable_prime_fast(value: &Integer, rounds: u64) -> bool {
    crate::integer::probable_prime(value, rounds)
}

/// Returns the BPSW verdict (deterministic u64 gate, Baillie-PSW otherwise); product checks verify composites but this result itself is not a primality proof.
pub fn is_probable_prime(value: &Integer) -> bool {
    if let Some(small) = value.to_u64() {
        return is_prime_u64(small);
    }
    is_baillie_psw_prime(value)
}

pub fn is_baillie_psw_prime(value: &Integer) -> bool {
    if let Some(small) = value.to_u64() {
        return is_prime_u64(small);
    }
    if value < &Integer::from(2) {
        return false;
    }
    if euler_base2_composite_detected(value) {
        return false;
    }
    if !passes_miller_rabin_base2(value) {
        return false;
    }
    is_strong_lucas_probable_prime(value)
}

static FROBENIUS_CALLS: AtomicU64 = AtomicU64::new(0);

pub fn frobenius_calls() -> u64 {
    FROBENIUS_CALLS.load(Ordering::SeqCst)
}

pub fn reset_frobenius_calls() {
    FROBENIUS_CALLS.store(0, Ordering::SeqCst);
}

pub fn is_probable_prime_with_options(
    value: &Integer,
    options: &crate::options::FactorOptions,
) -> bool {
    if !is_probable_prime(value) {
        return false;
    }
    if !options.use_frobenius {
        return true;
    }
    is_frobenius_probable_prime(value)
}

pub fn is_frobenius_probable_prime(value: &Integer) -> bool {
    FROBENIUS_CALLS.fetch_add(1, Ordering::SeqCst);
    if value < &Integer::from(2) {
        return false;
    }
    for prime in DETERMINISTIC_U64_BASES {
        let small = Integer::from(prime);
        if value == &small {
            return true;
        }
        if value % &small == Integer::from(0) {
            return false;
        }
    }
    if value % Integer::from(2) == Integer::from(0) {
        return false;
    }
    if is_integer_square(value) {
        return false;
    }
    let mut discriminant = Integer::from(5);
    let mut selected: Option<(Integer, Integer)> = None;
    for _ in 0..10_000 {
        match jacobi_symbol(&discriminant, value) {
            -1 => {
                let quotient = (Integer::from(1) - &discriminant) / Integer::from(4);
                selected = Some((discriminant.clone(), quotient));
                break;
            }
            0 => return false,
            _ => {}
        }
        if discriminant > Integer::from(0) {
            discriminant = -&discriminant - Integer::from(2);
        } else {
            discriminant = -&discriminant + Integer::from(2);
        }
    }
    let (discriminant, quotient) = match selected {
        Some(pair) => pair,
        None => return false,
    };
    let discriminant_mod = discriminant.rem_euclid(value);
    let quotient_mod = quotient.rem_euclid(value);
    let exponent = value + &Integer::from(1);
    let (lucas_u, lucas_v) =
        frobenius_lucas_pair(value, &discriminant_mod, &quotient_mod, &exponent);
    if !lucas_u.is_zero() {
        return false;
    }
    lucas_v == (&quotient_mod * Integer::from(2)).rem_euclid(value)
}

fn frobenius_lucas_pair(
    value: &Integer,
    discriminant_mod: &Integer,
    quotient_mod: &Integer,
    exponent: &Integer,
) -> (Integer, Integer) {
    let one = Integer::from(1);
    let two = Integer::from(2);
    let half = (value + &one) / &two;
    let mut lucas_u = Integer::from(0);
    let mut lucas_v = Integer::from(2);
    let mut lucas_q = Integer::from(1);
    let width = exponent.bits();
    let mut index = width;
    while index > 0 {
        index -= 1;
        let doubled_u = (&lucas_u * &lucas_v).rem_euclid(value);
        let doubled_v = ((&lucas_v * &lucas_v) - (&lucas_q * &two)).rem_euclid(value);
        let doubled_q = (&lucas_q * &lucas_q).rem_euclid(value);
        lucas_u = doubled_u;
        lucas_v = doubled_v;
        lucas_q = doubled_q;
        if exponent.bit(index) {
            let next_u = ((&lucas_u + &lucas_v) * &half).rem_euclid(value);
            let next_v = ((discriminant_mod * &lucas_u + &lucas_v) * &half).rem_euclid(value);
            let next_q = (quotient_mod * &lucas_q).rem_euclid(value);
            lucas_u = next_u;
            lucas_v = next_v;
            lucas_q = next_q;
        }
    }
    (lucas_u, lucas_v)
}

fn euler_base2_composite_detected(value: &Integer) -> bool {
    euler_composite_detected(value, &Integer::from(2))
}

fn euler_composite_detected(value: &Integer, base: &Integer) -> bool {
    let zero = Integer::from(0);
    let one = Integer::from(1);
    let two = Integer::from(2);
    if value <= &two {
        return false;
    }
    if value % &two == zero {
        return true;
    }
    if base % value == zero {
        return value != base;
    }
    let symbol = jacobi_symbol(base, value);
    if symbol == 0 {
        return true;
    }
    let exponent = (value - &one) / &two;
    let residue = crate::integer::mod_pow(base, &exponent, value);
    if symbol == 1 {
        residue != one
    } else {
        residue != value - &one
    }
}

pub fn passes_miller_rabin_base2(value: &Integer) -> bool {
    let one = Integer::from(1);
    let two = Integer::from(2);
    if value == &two {
        return true;
    }
    if value < &two || value % &two == Integer::from(0) {
        return false;
    }
    let mut odd = value - &one;
    let mut shifts = 0_u32;
    while &odd % &two == Integer::from(0) {
        odd = &odd / &two;
        shifts += 1;
    }
    let mut root = crate::integer::mod_pow(&two, &odd, value);
    if root == one || root == value - &one {
        return true;
    }
    for _ in 1..shifts {
        root = (&root * &root) % value;
        if root == value - &one {
            return true;
        }
    }
    false
}

pub fn jacobi_symbol(numerator: &Integer, denominator: &Integer) -> i32 {
    let one = Integer::from(1);
    let two = Integer::from(2);
    if denominator == &one {
        return 1;
    }
    if denominator < &two || denominator % &two == Integer::from(0) {
        return 0;
    }
    let mut top = numerator.rem_euclid(denominator);
    let mut bottom = denominator.clone();
    let mut sign = 1_i32;
    let zero = Integer::from(0);
    loop {
        if bottom == one {
            return sign;
        }
        if top == zero {
            return 0;
        }
        while &(&top / &two) * &two == top {
            top = &top / &two;
            let residue = &bottom % Integer::from(8);
            if residue == Integer::from(3) || residue == Integer::from(5) {
                sign = -sign;
            }
        }
        std::mem::swap(&mut top, &mut bottom);
        if &top % Integer::from(4) == Integer::from(3)
            && &bottom % Integer::from(4) == Integer::from(3)
        {
            sign = -sign;
        }
        top = top.rem_euclid(&bottom);
    }
}

pub fn is_strong_lucas_probable_prime(value: &Integer) -> bool {
    if value < &Integer::from(2) {
        return false;
    }
    for prime in DETERMINISTIC_U64_BASES {
        let small = Integer::from(prime);
        if value == &small {
            return true;
        }
        if value % &small == Integer::from(0) {
            return false;
        }
    }
    if value % Integer::from(2) == Integer::from(0) {
        return false;
    }
    if is_integer_square(value) {
        return false;
    }
    let mut discriminant = Integer::from(5);
    for _ in 0..10_000 {
        match jacobi_symbol(&discriminant, value) {
            -1 => return strong_lucas_with_discriminant(value, &discriminant),
            0 => return false,
            _ => {}
        }
        if discriminant > Integer::from(0) {
            discriminant = -&discriminant - Integer::from(2);
        } else {
            discriminant = -&discriminant + Integer::from(2);
        }
    }
    false
}

fn strong_lucas_with_discriminant(value: &Integer, discriminant: &Integer) -> bool {
    let one = Integer::from(1);
    let two = Integer::from(2);
    let selector = Integer::from(1);
    let quotient = (Integer::from(1) - discriminant) / Integer::from(4);
    let quotient_mod = quotient.rem_euclid(value);
    let discriminant_mod = discriminant.rem_euclid(value);
    let half = (value + &one) / &two;
    let mut odd = value + &one;
    let mut shifts = 0_u32;
    while &odd % &two == Integer::from(0) {
        odd = &odd / &two;
        shifts += 1;
    }
    let mut lucas_u = Integer::from(1);
    let mut lucas_v = selector.clone();
    let mut lucas_q = quotient_mod.clone();
    let width = odd.bits();
    let mut index = width;
    while index > 1 {
        index -= 1;
        let doubled_u = (&lucas_u * &lucas_v).rem_euclid(value);
        let doubled_v = ((&lucas_v * &lucas_v) - (&lucas_q * &two)).rem_euclid(value);
        let doubled_q = (&lucas_q * &lucas_q).rem_euclid(value);
        lucas_u = doubled_u;
        lucas_v = doubled_v;
        lucas_q = doubled_q;
        if odd.bit(index - 1) {
            let next_u = (((&selector * &lucas_u) + &lucas_v) * &half).rem_euclid(value);
            let next_v = (((&discriminant_mod * &lucas_u) + (&selector * &lucas_v)) * &half)
                .rem_euclid(value);
            let next_q = (&quotient_mod * &lucas_q).rem_euclid(value);
            lucas_u = next_u;
            lucas_v = next_v;
            lucas_q = next_q;
        }
    }
    if lucas_u.is_zero() {
        return true;
    }
    let mut current_v = lucas_v;
    let mut current_q = lucas_q;
    for _ in 0..shifts {
        if current_v.is_zero() {
            return true;
        }
        current_v = ((&current_v * &current_v) - (&current_q * &two)).rem_euclid(value);
        current_q = (&current_q * &current_q).rem_euclid(value);
    }
    false
}

pub fn is_prime_u64_deterministic(value: u64) -> bool {
    is_prime_u64(value)
}

pub const AKS_PROOF_MAX_BITS: u64 = 64;

static AKS_CALLS: AtomicU64 = AtomicU64::new(0);

pub fn aks_calls() -> u64 {
    AKS_CALLS.load(Ordering::SeqCst)
}

pub fn reset_aks_calls() {
    AKS_CALLS.store(0, Ordering::SeqCst);
}

pub fn is_prime_aks_proof(value: &Integer) -> bool {
    AKS_CALLS.fetch_add(1, Ordering::SeqCst);
    aks::is_prime(value)
}

mod aks {
    use super::AKS_PROOF_MAX_BITS;
    use super::Integer;
    use super::mul_mod;
    use num_traits::ToPrimitive;
    use std::cmp::Ordering;

    pub fn is_prime(value: &Integer) -> bool {
        if value < &Integer::from(2) {
            return false;
        }
        if value.bits() > AKS_PROOF_MAX_BITS {
            return false;
        }
        let n = match value.to_u64() {
            Some(n) => n,
            None => return false,
        };
        if n < 2 {
            return false;
        }
        if n == 2 || n == 3 {
            return true;
        }
        if is_perfect_power(n) {
            return false;
        }
        let floor_log = 63 - n.leading_zeros();
        let ceil_log = if n.is_power_of_two() {
            floor_log
        } else {
            floor_log + 1
        };
        let bound_sq = ceil_log as u64 * ceil_log as u64;
        let mut r = 2_u64;
        loop {
            if r >= n {
                return true;
            }
            if gcd_u64(n, r) > 1 {
                return false;
            }
            if order_exceeds(n, r, bound_sq) {
                break;
            }
            r += 1;
        }
        let top = r.min(n - 1);
        let mut witness = 2_u64;
        while witness <= top {
            if gcd_u64(witness, n) > 1 {
                return false;
            }
            witness += 1;
        }
        if n <= r {
            return true;
        }
        let root = isqrt_u64(euler_phi(r));
        let limit = root * ceil_log as u64;
        let mut a = 1_u64;
        while a <= limit {
            if !passes_poly(n, r, a) {
                return false;
            }
            a += 1;
        }
        true
    }

    fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            let next = a % b;
            a = b;
            b = next;
        }
        a
    }

    fn pow_order(base: u64, exp: u32, target: u64) -> Ordering {
        let goal = target as u128;
        let cap = goal + 1;
        let mut result: u128 = 1;
        let mut factor: u128 = base as u128;
        let mut remaining = exp;
        while remaining > 0 {
            if remaining & 1 == 1 {
                result = result.saturating_mul(factor).min(cap);
            }
            remaining >>= 1;
            if remaining > 0 {
                factor = factor.saturating_mul(factor).min(cap);
            }
        }
        if result == cap {
            Ordering::Greater
        } else {
            result.cmp(&goal)
        }
    }

    fn nth_root_exact(n: u64, exp: u32) -> Option<u64> {
        let mut hi = 1_u64 << (64 / exp + 1);
        if hi > n {
            hi = n;
        }
        let mut lo = 2_u64;
        while lo <= hi {
            let mid = lo + (hi - lo) / 2;
            match pow_order(mid, exp, n) {
                Ordering::Equal => return Some(mid),
                Ordering::Less => lo = mid + 1,
                Ordering::Greater => {
                    if mid == 0 {
                        break;
                    }
                    hi = mid - 1;
                }
            }
        }
        None
    }

    fn is_perfect_power(n: u64) -> bool {
        if n < 4 {
            return false;
        }
        let max_exp = 64 - n.leading_zeros();
        let mut exp = 2_u32;
        while exp <= max_exp {
            if nth_root_exact(n, exp).is_some() {
                return true;
            }
            exp += 1;
        }
        false
    }

    fn order_exceeds(n: u64, r: u64, bound: u64) -> bool {
        let residue = n % r;
        let mut current = residue;
        let mut steps = 0_u64;
        while current != 1 {
            if steps >= bound {
                return true;
            }
            current = mul_mod(current, residue, r);
            steps += 1;
            if steps > r {
                return true;
            }
        }
        false
    }

    fn euler_phi(mut n: u64) -> u64 {
        let mut result = n;
        let mut p = 2_u64;
        while p <= n / p {
            if n.is_multiple_of(p) {
                while n.is_multiple_of(p) {
                    n /= p;
                }
                result -= result / p;
            }
            p += 1;
        }
        if n > 1 {
            result -= result / n;
        }
        result
    }

    fn isqrt_u64(n: u64) -> u64 {
        if n < 2 {
            return n;
        }
        let mut x = 1_u64 << (64 - n.leading_zeros()).div_ceil(2);
        loop {
            let next = (x + n / x) / 2;
            if next >= x {
                return x;
            }
            x = next;
        }
    }

    fn poly_mul(left: &[u64], right: &[u64], degree: usize, modulus: u64) -> Vec<u64> {
        let mut out = vec![0_u64; degree];
        for (i, first) in left.iter().enumerate() {
            if *first == 0 {
                continue;
            }
            for (j, second) in right.iter().enumerate() {
                if *second == 0 {
                    continue;
                }
                let pos = (i + j) % degree;
                let product = (*first as u128 * *second as u128) % modulus as u128;
                out[pos] = ((out[pos] as u128 + product) % modulus as u128) as u64;
            }
        }
        out
    }

    fn poly_pow(base: &[u64], mut exp: u64, degree: usize, modulus: u64) -> Vec<u64> {
        let mut result = vec![0_u64; degree];
        result[0] = 1 % modulus;
        let mut factor = base.to_vec();
        while exp > 0 {
            if exp & 1 == 1 {
                result = poly_mul(&result, &factor, degree, modulus);
            }
            factor = poly_mul(&factor, &factor, degree, modulus);
            exp >>= 1;
        }
        result
    }

    fn passes_poly(n: u64, r: u64, a: u64) -> bool {
        let degree = r as usize;
        let mut base = vec![0_u64; degree];
        base[0] = a % n;
        let one = 1 % r;
        base[one as usize] = (base[one as usize] + 1) % n;
        let left = poly_pow(&base, n, degree, n);
        let mut right = vec![0_u64; degree];
        right[0] = a % n;
        let pos = (n % r) as usize;
        right[pos] = (right[pos] + 1) % n;
        left == right
    }
}

pub fn integer_abs(value: &Integer) -> Integer {
    if value < &Integer::from(0) {
        -value.clone()
    } else {
        value.clone()
    }
}

#[cfg(test)]
mod prime_tests {
    use super::*;
    use std::time::Duration;
    #[test]
    fn unified_u64_gate_matches_both_former_paths() {
        assert_eq!(
            DETERMINISTIC_U64_BASES,
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
        );
        let primes: &[u64] = &[
            2,
            3,
            5,
            7,
            11,
            13,
            17,
            19,
            23,
            29,
            31,
            37,
            41,
            43,
            47,
            53,
            59,
            61,
            67,
            71,
            73,
            79,
            83,
            89,
            97,
            1_000_000_007,
            2_147_483_647,
            18_446_744_073_709_551_557,
        ];
        let composites: &[u64] = &[
            0,
            1,
            4,
            6,
            8,
            9,
            561,
            1105,
            1729,
            2465,
            2821,
            6601,
            8911,
            10585,
            15841,
            29341,
            41041,
            46657,
            52633,
            62745,
            63973,
            75361,
            3_215_031_751,
            2_152_302_898_747,
            3_474_749_660_383,
            341_550_071_728_321,
            3_825_123_056_546_413_051,
            u64::MAX - 1,
            u64::MAX,
        ];
        for value in primes {
            assert!(is_prime_u64(*value), "prime missed: {value}");
            assert!(is_prime_u64_deterministic(*value), "alias missed: {value}");
            assert!(
                is_probable_prime(&Integer::from(*value)),
                "gate missed: {value}"
            );
        }
        for value in composites {
            assert!(!is_prime_u64(*value), "false positive: {value}");
            assert!(
                !is_prime_u64_deterministic(*value),
                "alias false positive: {value}"
            );
            assert!(
                !is_probable_prime(&Integer::from(*value)),
                "gate false positive: {value}"
            );
        }
        for value in primes.iter().chain(composites.iter()) {
            assert_eq!(
                is_prime_u64(*value),
                is_prime_u64_deterministic(*value),
                "path split on {value}"
            );
        }
    }
    #[test]
    fn segmented_walk_matches_sieve() {
        let expected = sieve_primes(5_000_003);
        let mut visited = Vec::new();
        assert!(for_each_prime(5_000_003, None, |prime| {
            visited.push(prime);
            true
        }));
        assert_eq!(visited, expected);
    }
    #[test]
    fn walk_stops_on_visitor_signal() {
        let mut count = 0_usize;
        assert!(for_each_prime(1_000_000, None, |_| {
            count += 1;
            count < 25
        }));
        assert_eq!(count, 25);
    }
    #[test]
    fn walk_aborts_on_expired_deadline() {
        let past = Some(std::time::Instant::now() + Duration::from_millis(1));
        std::thread::sleep(Duration::from_millis(5));
        let mut count = 0_usize;
        assert!(!for_each_prime(200_000_000, past, |_| {
            count += 1;
            true
        }));
        assert_eq!(count, 0);
    }
}

#[cfg(test)]
mod bpsw_tests {
    use super::*;

    #[test]
    fn base_two_pseudoprime_falls_at_lucas_step() {
        let value = Integer::from(2047);
        assert!(passes_miller_rabin_base2(&value));
        assert!(!is_strong_lucas_probable_prime(&value));
        assert!(!is_baillie_psw_prime(&value));
        assert!(!is_probable_prime(&value));
    }

    #[test]
    fn selfridge_lucas_pseudoprimes_stay_composite() {
        for small in [323_u64, 377, 1159, 1829, 19093] {
            let value = Integer::from(small);
            assert!(
                !is_strong_lucas_probable_prime(&value),
                "lucas miss: {small}"
            );
            assert!(!is_baillie_psw_prime(&value), "bpsw miss: {small}");
            assert!(!is_probable_prime(&value), "gate miss: {small}");
        }
    }

    #[test]
    fn strong_lucas_pseudoprimes_fall_at_miller_rabin_step() {
        for small in [5459_u64, 5777, 10877, 16109, 18971] {
            let value = Integer::from(small);
            assert!(
                is_strong_lucas_probable_prime(&value),
                "lucas miss: {small}"
            );
            assert!(!passes_miller_rabin_base2(&value), "mr2 miss: {small}");
            assert!(!is_baillie_psw_prime(&value), "bpsw miss: {small}");
            assert!(!is_probable_prime(&value), "gate miss: {small}");
        }
    }

    #[test]
    fn bpsw_matches_deterministic_gate_below_5000() {
        for small in 0_u64..5000 {
            let value = Integer::from(small);
            let expected = is_prime_u64(small);
            let combined =
                passes_miller_rabin_base2(&value) && is_strong_lucas_probable_prime(&value);
            assert_eq!(combined, expected, "lucas-stack split: {small}");
            assert_eq!(
                is_baillie_psw_prime(&value),
                expected,
                "bpsw split: {small}"
            );
            assert_eq!(is_probable_prime(&value), expected, "gate split: {small}");
        }
    }

    #[test]
    fn extra_strong_and_carmichael_cases_stay_composite() {
        for small in [
            561_u64,
            1105,
            1729,
            989,
            3239,
            5777,
            10877,
            3_215_031_751,
            2_152_302_898_747,
        ] {
            let value = Integer::from(small);
            assert!(!is_baillie_psw_prime(&value), "bpsw miss: {small}");
            assert!(!is_probable_prime(&value), "gate miss: {small}");
        }
    }

    #[test]
    fn bpsw_accepts_primes_on_both_sides_of_u64() {
        for small in [2_u64, 3, 5, 4231, 5347, 1_000_000_007, 2_147_483_647] {
            let value = Integer::from(small);
            assert!(is_baillie_psw_prime(&value), "bpsw prime miss: {small}");
        }
        let mersenne127: Integer = (Integer::from(1) << 127) - Integer::from(1);
        assert!(passes_miller_rabin_base2(&mersenne127));
        assert!(is_strong_lucas_probable_prime(&mersenne127));
        assert!(is_baillie_psw_prime(&mersenne127));
        assert!(is_probable_prime(&mersenne127));
        let big_composite = &mersenne127 * &Integer::from(3);
        assert!(!is_baillie_psw_prime(&big_composite));
        assert!(!is_probable_prime(&big_composite));
    }

    #[test]
    fn euler_filter_pass_still_reaches_bpsw_verdict() {
        for small in [561_u64, 1105, 1729] {
            let value = Integer::from(small);
            assert!(
                !euler_base2_composite_detected(&value),
                "euler filter miss: {small}"
            );
            assert!(!is_baillie_psw_prime(&value), "bpsw miss: {small}");
            assert!(!is_probable_prime(&value), "gate miss: {small}");
        }
    }

    #[test]
    fn dual_filter_pass_still_falls_at_lucas_step() {
        for small in [2047_u64, 3277, 8321] {
            let value = Integer::from(small);
            assert!(
                !euler_base2_composite_detected(&value),
                "euler filter miss: {small}"
            );
            assert!(passes_miller_rabin_base2(&value), "mr2 miss: {small}");
            assert!(
                !is_strong_lucas_probable_prime(&value),
                "lucas miss: {small}"
            );
            assert!(!is_baillie_psw_prime(&value), "bpsw miss: {small}");
            assert!(!is_probable_prime(&value), "gate miss: {small}");
        }
    }

    #[test]
    fn euler_filter_covers_big_gate_path() {
        let mersenne127: Integer = (Integer::from(1) << 127) - Integer::from(1);
        assert!(!euler_base2_composite_detected(&mersenne127));
        assert!(is_baillie_psw_prime(&mersenne127));
        let big_composite = &mersenne127 * &Integer::from(3);
        assert!(euler_base2_composite_detected(&big_composite));
        assert!(!is_baillie_psw_prime(&big_composite));
        assert!(!is_probable_prime(&big_composite));
    }

    #[test]
    fn jacobi_symbol_matches_known_values() {
        assert_eq!(jacobi_symbol(&Integer::from(2), &Integer::from(7)), 1);
        assert_eq!(jacobi_symbol(&Integer::from(2), &Integer::from(3)), -1);
        assert_eq!(
            jacobi_symbol(&Integer::from(1001), &Integer::from(9907)),
            -1
        );
        assert_eq!(jacobi_symbol(&Integer::from(19), &Integer::from(45)), 1);
        assert_eq!(jacobi_symbol(&Integer::from(8), &Integer::from(21)), -1);
        assert_eq!(jacobi_symbol(&Integer::from(5), &Integer::from(5)), 0);
        assert_eq!(jacobi_symbol(&Integer::from(30), &Integer::from(59)), -1);
    }
}

#[cfg(test)]
mod frobenius_tests {
    use super::*;
    use crate::options::FactorOptions;

    const SMALL_PRIMES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    const CARMICHAEL: [u64; 3] = [561, 1105, 1729];
    const BASE_TWO_PSEUDOPRIMES: [u64; 3] = [2047, 3277, 8321];
    const STRONG_LUCAS_PSEUDOPRIMES: [u64; 5] = [5459, 5777, 10877, 16109, 18971];

    fn flag_off() -> FactorOptions {
        FactorOptions::default()
    }

    fn flag_on() -> FactorOptions {
        FactorOptions {
            use_frobenius: true,
            ..FactorOptions::default()
        }
    }

    fn mersenne127() -> Integer {
        (Integer::from(1) << 127) - Integer::from(1)
    }

    static FROBENIUS_SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

    fn serial_guard() -> std::sync::MutexGuard<'static, ()> {
        FROBENIUS_SERIAL
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    #[test]
    fn default_paths_never_call_frobenius() {
        let _guard = serial_guard();
        reset_frobenius_calls();
        for small in SMALL_PRIMES {
            assert!(is_probable_prime(&Integer::from(small)));
        }
        for small in CARMICHAEL
            .iter()
            .chain(BASE_TWO_PSEUDOPRIMES.iter())
            .chain(STRONG_LUCAS_PSEUDOPRIMES.iter())
        {
            assert!(!is_probable_prime(&Integer::from(*small)));
        }
        assert!(is_probable_prime(&mersenne127()));
        assert!(is_baillie_psw_prime(&mersenne127()));
        assert!(!is_probable_prime(&(&mersenne127() * &Integer::from(3))));
        assert!(!is_baillie_psw_prime(&(&mersenne127() * &Integer::from(3))));
        let report =
            crate::factor::factor(&Integer::from(10403), &flag_off()).expect("10403 factors");
        assert!(
            report
                .timings
                .iter()
                .all(|entry| entry.stage != "frobenius")
        );
        assert_eq!(frobenius_calls(), 0);
    }

    #[test]
    fn flag_off_gate_matches_plain_gate_without_calling_frobenius() {
        let _guard = serial_guard();
        reset_frobenius_calls();
        let off = flag_off();
        let mut corpus: Vec<Integer> = SMALL_PRIMES
            .iter()
            .map(|small| Integer::from(*small))
            .collect();
        for small in CARMICHAEL
            .iter()
            .chain(BASE_TWO_PSEUDOPRIMES.iter())
            .chain(STRONG_LUCAS_PSEUDOPRIMES.iter())
        {
            corpus.push(Integer::from(*small));
        }
        corpus.push(mersenne127());
        corpus.push(&mersenne127() * &Integer::from(3));
        for value in &corpus {
            assert_eq!(
                is_probable_prime_with_options(value, &off),
                is_probable_prime(value)
            );
        }
        assert_eq!(frobenius_calls(), 0);
    }

    #[test]
    fn flag_on_accepts_known_primes_as_and_composition() {
        let _guard = serial_guard();
        let on = flag_on();
        for small in SMALL_PRIMES {
            let value = Integer::from(small);
            assert!(is_frobenius_probable_prime(&value));
            assert!(is_probable_prime_with_options(&value, &on));
            assert_eq!(
                is_probable_prime_with_options(&value, &on),
                is_baillie_psw_prime(&value) && is_frobenius_probable_prime(&value)
            );
        }
        let big = mersenne127();
        assert!(is_frobenius_probable_prime(&big));
        assert!(is_probable_prime_with_options(&big, &on));
        assert_eq!(
            is_probable_prime_with_options(&big, &on),
            is_baillie_psw_prime(&big) && is_frobenius_probable_prime(&big)
        );
    }

    #[test]
    fn flag_on_rejects_composites_as_and_composition() {
        let _guard = serial_guard();
        let on = flag_on();
        for small in CARMICHAEL
            .iter()
            .chain(BASE_TWO_PSEUDOPRIMES.iter())
            .chain(STRONG_LUCAS_PSEUDOPRIMES.iter())
        {
            let value = Integer::from(*small);
            assert!(!is_probable_prime_with_options(&value, &on));
            assert_eq!(
                is_probable_prime_with_options(&value, &on),
                is_baillie_psw_prime(&value) && is_frobenius_probable_prime(&value)
            );
        }
        let big_composite = &mersenne127() * &Integer::from(3);
        assert!(!is_probable_prime_with_options(&big_composite, &on));
        assert_eq!(
            is_probable_prime_with_options(&big_composite, &on),
            is_baillie_psw_prime(&big_composite) && is_frobenius_probable_prime(&big_composite)
        );
    }

    #[test]
    fn fast_path_with_any_rounds_never_calls_frobenius() {
        let _guard = serial_guard();
        reset_frobenius_calls();
        let prime = Integer::from(101);
        for rounds in [0_u64, 1, 64, 1000] {
            assert!(super::is_probable_prime_fast(&prime, rounds));
        }
        let lucas_pass = Integer::from(5459);
        assert!(!super::is_probable_prime_fast(&lucas_pass, u64::MAX));
        assert_eq!(frobenius_calls(), 0);
    }

    #[test]
    fn frobenius_pins_documented_strong_lucas_cases() {
        let _guard = serial_guard();
        for small in [5459_u64, 16109, 18971] {
            assert!(
                !is_frobenius_probable_prime(&Integer::from(small)),
                "frobenius must reject strong-Lucas pseudoprime {small}"
            );
        }
        for small in [5777_u64, 10877] {
            assert!(
                is_frobenius_probable_prime(&Integer::from(small)),
                "{small} is a Frobenius pseudoprime for Selfridge parameters"
            );
            assert!(
                !is_probable_prime_with_options(&Integer::from(small), &flag_on()),
                "strengthened gate still rejects {small} through BPSW"
            );
        }
    }
}

#[cfg(test)]
mod rounds_proof_tests {
    use super::*;
    use crate::integer::probable_prime as fast_rounds_prime;

    #[test]
    fn max_rounds_never_enables_proof_path() {
        let lucas_pass = Integer::from(5459);
        assert!(is_strong_lucas_probable_prime(&lucas_pass));
        assert!(!is_baillie_psw_prime(&lucas_pass));
        assert!(!is_probable_prime_fast(&lucas_pass, u64::MAX));
        assert!(!fast_rounds_prime(&lucas_pass, u64::MAX));
    }

    #[test]
    fn rounds_tuning_stays_inside_fast_path() {
        let prime = Integer::from(101);
        for rounds in [0_u64, 1, 12, 16, 17, 64, 1000] {
            assert!(is_probable_prime_fast(&prime, rounds));
        }
        let lucas_pass = Integer::from(5459);
        for rounds in [0_u64, 1, 12, 16, 17, 64, 1000, u64::MAX] {
            assert!(!is_probable_prime_fast(&lucas_pass, rounds));
        }
    }

    #[test]
    fn explicit_proof_accepts_prime_without_rounds() {
        let prime = Integer::from(101);
        assert!(is_baillie_psw_prime(&prime));
        assert!(is_probable_prime_with_options(
            &prime,
            &crate::options::FactorOptions::default()
        ));
        assert!(is_probable_prime(&prime));
        let lucas_pass = Integer::from(5459);
        assert!(!is_baillie_psw_prime(&lucas_pass));
        assert!(!is_probable_prime_with_options(
            &lucas_pass,
            &crate::options::FactorOptions::default()
        ));
    }
}

#[cfg(test)]
mod aks_tests {
    use super::*;

    static AKS_SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

    fn aks_guard() -> std::sync::MutexGuard<'static, ()> {
        AKS_SERIAL
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    const AKS_CORPUS_PRIMES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    const AKS_CORPUS_COMPOSITES: [u64; 11] = [4, 6, 8, 9, 561, 1105, 1729, 2047, 3277, 5459, 10403];

    #[test]
    fn hot_path_never_enters_aks() {
        let _guard = aks_guard();
        reset_aks_calls();
        for small in AKS_CORPUS_PRIMES {
            assert!(is_probable_prime(&Integer::from(small)));
        }
        for small in AKS_CORPUS_COMPOSITES {
            assert!(!is_probable_prime(&Integer::from(small)));
        }
        for small in AKS_CORPUS_PRIMES {
            assert!(is_baillie_psw_prime(&Integer::from(small)));
        }
        for small in AKS_CORPUS_COMPOSITES {
            assert!(!is_baillie_psw_prime(&Integer::from(small)));
        }
        let report = crate::factor::factor(
            &Integer::from(10403),
            &crate::options::FactorOptions::default(),
        )
        .expect("10403 factors");
        assert!(report.timings.iter().all(|entry| entry.stage != "aks"));
        assert_eq!(aks_calls(), 0);
    }

    #[test]
    fn fast_path_with_any_rounds_never_calls_aks() {
        let _guard = aks_guard();
        reset_aks_calls();
        let prime = Integer::from(101);
        for rounds in [0_u64, 1, 64, 1000] {
            assert!(super::is_probable_prime_fast(&prime, rounds));
        }
        let lucas_pass = Integer::from(5459);
        assert!(!super::is_probable_prime_fast(&lucas_pass, u64::MAX));
        assert_eq!(aks_calls(), 0);
        assert!(is_prime_aks_proof(&prime));
        assert_eq!(aks_calls(), 1);
    }

    #[test]
    fn over_cap_input_refused() {
        let _guard = aks_guard();
        reset_aks_calls();
        let wide = Integer::from(1) << 70;
        assert!(!is_prime_aks_proof(&wide));
        let wider = Integer::from(1) << 100;
        assert!(!is_prime_aks_proof(&wider));
        let huge = (Integer::from(1) << 127) - Integer::from(1);
        assert!(!is_prime_aks_proof(&huge));
        assert_eq!(aks_calls(), 3);
    }

    #[test]
    fn explicit_proof_accepts_primes_rejects_composites() {
        let _guard = aks_guard();
        reset_aks_calls();
        for small in AKS_CORPUS_PRIMES {
            assert!(
                is_prime_aks_proof(&Integer::from(small)),
                "aks must accept prime {small}"
            );
        }
        for small in AKS_CORPUS_COMPOSITES {
            assert!(
                !is_prime_aks_proof(&Integer::from(small)),
                "aks must reject composite {small}"
            );
        }
        assert!(aks_calls() > 0);
    }

    #[test]
    fn aks_matches_sieve_below_200() {
        let _guard = aks_guard();
        let sieve = sieve_primes(200);
        for small in 0_u64..200 {
            let expected = sieve.contains(&small);
            assert_eq!(
                is_prime_aks_proof(&Integer::from(small)),
                expected,
                "aks split at {small}"
            );
            assert_eq!(is_prime_u64(small), expected, "sieve split at {small}");
        }
    }
}
