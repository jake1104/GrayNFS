use crate::integer::{Integer, IntegerExt};
use std::sync::OnceLock;

pub fn small_primes() -> &'static [u64] {
    static TABLE: OnceLock<Vec<u64>> = OnceLock::new();
    TABLE.get_or_init(|| sieve_primes(10_000))
}

pub fn is_prime_u64(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    for prime in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if value.is_multiple_of(prime) {
            return value == prime;
        }
    }
    let mut d = value - 1;
    let mut s = 0_u32;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }
    let bases = [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for base in bases {
        let mut x = mod_pow(base % value, d, value);
        if x == 1 || x == value - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod(x, x, value);
            if x == value - 1 {
                break;
            }
        }
        if x != value - 1 {
            return false;
        }
    }
    true
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

pub fn is_probable_prime(value: &Integer) -> bool {
    value >= &Integer::from(2) && value.is_probably_prime(32)
}

pub fn is_prime_u64_deterministic(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    const SMALL_PRIMES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if value < 2 {
        return false;
    }
    for prime in SMALL_PRIMES {
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
    'bases: for base in SMALL_PRIMES {
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
