use num_bigint::BigInt;
use num_traits::{Num, One, Signed, ToPrimitive, Zero};

pub type Integer = BigInt;

pub trait IntegerExt {
    fn to_u64_saturating(&self) -> u64;
    fn to_u128_saturating(&self) -> u128;
    fn is_probably_prime(&self, rounds: u64) -> bool;
    fn sqrt(&self) -> Integer;
}

impl IntegerExt for Integer {
    fn to_u64_saturating(&self) -> u64 {
        self.to_u64()
            .unwrap_or_else(|| if self.is_negative() { 0 } else { u64::MAX })
    }

    fn to_u128_saturating(&self) -> u128 {
        self.to_u128()
            .unwrap_or_else(|| if self.is_negative() { 0 } else { u128::MAX })
    }

    fn is_probably_prime(&self, rounds: u64) -> bool {
        probable_prime(self, rounds)
    }

    fn sqrt(&self) -> Integer {
        integer_sqrt(self)
    }
}

pub fn integer_sqrt(value: &Integer) -> Integer {
    if value <= &Integer::from(0) {
        return Integer::from(0);
    }
    let bits = value.bits();
    let mut guess = Integer::from(1) << bits.div_ceil(2);
    loop {
        let next = (&guess + value / &guess) / Integer::from(2);
        if next >= guess {
            return guess;
        }
        guess = next;
    }
}

pub fn integer_abs(value: &Integer) -> Integer {
    if value.is_negative() {
        -value.clone()
    } else {
        value.clone()
    }
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    b >>= b.trailing_zeros();
    loop {
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        b -= a;
        if b == 0 {
            return a << shift;
        }
        b >>= b.trailing_zeros();
    }
}

pub fn integer_gcd(left: &Integer, right: &Integer) -> Integer {
    let a = integer_abs(left);
    let b = integer_abs(right);
    if let (Some(small_a), Some(small_b)) = (a.to_u64(), b.to_u64()) {
        return Integer::from(gcd_u64(small_a, small_b));
    }
    let mut a = a;
    let mut b = b;
    while !b.is_zero() {
        let remainder = &a % &b;
        a = b;
        b = remainder;
    }
    a
}

pub fn probable_prime(value: &Integer, rounds: u64) -> bool {
    if value < &Integer::from(2) {
        return false;
    }
    let small_primes = crate::prime::DETERMINISTIC_U64_BASES;
    for prime in small_primes {
        let prime_value = Integer::from(prime);
        if value == &prime_value {
            return true;
        }
        if (value % &prime_value).is_zero() {
            return false;
        }
    }
    let mut odd = value - Integer::from(1);
    let mut shifts = 0_u32;
    while (&odd % Integer::from(2)).is_zero() {
        odd /= Integer::from(2);
        shifts += 1;
    }
    let bases: &[u64] = if value.bits() <= 64 {
        &crate::prime::DETERMINISTIC_U64_BASES
    } else {
        &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    };
    let actual_rounds = rounds.max(bases.len() as u64);
    for index in 0..actual_rounds {
        let base = Integer::from(bases[(index as usize) % bases.len()]);
        let mut x = mod_pow(&base, &odd, value);
        if x == Integer::one() || x == value - Integer::from(1) {
            continue;
        }
        let mut passed = false;
        for _ in 1..shifts {
            x = (&x * &x) % value;
            if x == value - Integer::from(1) {
                passed = true;
                break;
            }
        }
        if !passed {
            return false;
        }
    }
    true
}

pub fn mod_pow(base: &Integer, exponent: &Integer, modulus: &Integer) -> Integer {
    if modulus <= &Integer::from(1) {
        return Integer::from(0);
    }
    let mut base = base % modulus;
    let mut result = Integer::one();
    let mut exponent = exponent.clone();
    let two = Integer::from(2);
    while !exponent.is_zero() {
        if (&exponent % &two) == Integer::one() {
            result = (&result * &base) % modulus;
        }
        base = (&base * &base) % modulus;
        exponent /= &two;
    }
    result
}

pub fn parse_decimal(value: &str) -> Result<Integer, String> {
    if value.is_empty() {
        return Err("N must not be empty".to_string());
    }
    Integer::from_str_radix(value, 10).map_err(|error| error.to_string())
}

#[cfg(test)]
mod parse_decimal_tests {
    use super::{Integer, parse_decimal};

    const HOSTILE_CANDIDATES: &[&str] = &[
        "",
        "-",
        "+",
        "+0",
        "-0",
        "0",
        "0000",
        "007",
        "1e3",
        "0x10",
        " 12",
        "12 ",
        "\t12",
        "\n12",
        "12\n",
        "\u{0}",
        "1\u{0}",
        "12\u{0}34",
        "\u{660}\u{661}\u{662}",
        "\u{ff11}\u{ff12}",
        "abc",
        "--1",
        "++1",
        "1-",
        "- 1",
        "1_000",
        "1,000",
        ".5",
        "5.",
        "1.0",
    ];

    #[test]
    fn parse_decimal_never_panics_and_only_returns_ok_or_err() {
        for candidate in HOSTILE_CANDIDATES {
            match parse_decimal(candidate) {
                Ok(value) => {
                    let _ = value.to_string();
                }
                Err(message) => {
                    assert!(
                        !message.is_empty(),
                        "empty error message for candidate {:?}",
                        candidate
                    );
                }
            }
        }
    }

    #[test]
    fn parse_decimal_classifies_expected_candidates() {
        assert!(parse_decimal("").is_err());
        assert!(parse_decimal("-").is_err());
        assert!(parse_decimal("+").is_err());
        assert_eq!(parse_decimal("+0").unwrap(), Integer::from(0));
        assert_eq!(parse_decimal("-0").unwrap(), Integer::from(0));
        assert_eq!(parse_decimal("0").unwrap(), Integer::from(0));
        assert_eq!(parse_decimal("0000").unwrap(), Integer::from(0));
        assert_eq!(parse_decimal("007").unwrap(), Integer::from(7));
        assert!(parse_decimal("1e3").is_err());
        assert!(parse_decimal("0x10").is_err());
        assert!(parse_decimal(" 12").is_err());
        assert!(parse_decimal("12 ").is_err());
        assert!(parse_decimal("\t12").is_err());
        assert!(parse_decimal("1\u{0}").is_err());
        assert!(parse_decimal("\u{660}\u{661}\u{662}").is_err());
        assert!(parse_decimal("\u{ff11}\u{ff12}").is_err());
        assert!(parse_decimal("--1").is_err());
        assert!(parse_decimal("1,000").is_err());
        assert_eq!(parse_decimal("1_000").unwrap(), Integer::from(1000));
    }

    #[test]
    fn parse_decimal_accepts_very_long_decimal_strings() {
        let digits = "9".repeat(100_000);
        let value = parse_decimal(&digits).expect("long decimal string parses");
        assert_eq!(value.to_string(), digits);
    }
}
