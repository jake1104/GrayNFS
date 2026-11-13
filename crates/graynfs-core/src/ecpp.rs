use crate::integer::{integer_gcd, integer_sqrt};
use crate::prime::{is_prime_u64, sieve_primes};
use num_bigint::BigInt;
use num_traits::{Euclid, Num, One, Signed, ToPrimitive, Zero};
use std::sync::OnceLock;

pub const ECPP_CERT_FORMAT: &str = "graynfs-ecpp-cert/1";
pub const ECPP_MIN_BITS: u64 = 64;
pub const ECPP_MAX_BITS: u64 = 128;

const HEEGNER_ROWS: [(i64, i64); 9] = [
    (-3, 0),
    (-4, 1728),
    (-7, -3375),
    (-8, 8000),
    (-11, -32768),
    (-19, -884736),
    (-43, -884736000),
    (-67, -147197952000),
    (-163, -262537412640768000),
];

const STRIP_BOUND: usize = 200_000;
const RHO_EVAL_BUDGET: u64 = 40_000;
const RHO_ROUNDS: u32 = 2;
const POINT_X_TRIES: u64 = 12;
const NON_RESIDUE_TRIES: u64 = 1000;
const MAX_CHAIN_DEPTH: u32 = 16;
const LEVEL_CANDIDATE_CAP: usize = 8;
const SMALL_FILTER_BOUND: usize = 1000;
const MILLER_RABIN_BASES: [u64; 8] = [2, 3, 5, 7, 11, 13, 17, 19];

#[derive(Clone, Debug)]
pub struct EcppLevel {
    pub n: BigInt,
    pub a: BigInt,
    pub b: BigInt,
    pub x: BigInt,
    pub y: BigInt,
    pub m: BigInt,
    pub q: BigInt,
    pub discriminant: BigInt,
}

#[derive(Clone, Debug)]
pub struct EcppCertificate {
    pub n: BigInt,
    pub levels: Vec<EcppLevel>,
}

#[derive(Clone, Debug)]
struct ProjectivePoint {
    x: BigInt,
    y: BigInt,
    z: BigInt,
}

fn strip_primes() -> &'static [u64] {
    static TABLE: OnceLock<Vec<u64>> = OnceLock::new();
    TABLE.get_or_init(|| sieve_primes(STRIP_BOUND))
}

fn small_filter_primes() -> &'static [u64] {
    static TABLE: OnceLock<Vec<u64>> = OnceLock::new();
    TABLE.get_or_init(|| sieve_primes(SMALL_FILTER_BOUND))
}

pub fn parse_canonical_decimal(text: &str) -> Option<BigInt> {
    if text.is_empty() {
        return None;
    }
    if text.len() > 1 && text.starts_with('0') {
        return None;
    }
    if !text.bytes().all(|unit| unit.is_ascii_digit()) {
        return None;
    }
    BigInt::from_str_radix(text, 10).ok()
}

fn parse_discriminant_decimal(text: &str) -> Option<BigInt> {
    let digits = text.strip_prefix('-').unwrap_or(text);
    let magnitude = parse_canonical_decimal(digits)?;
    if magnitude.is_zero() {
        return None;
    }
    if text.starts_with('-') {
        Some(-magnitude)
    } else {
        Some(magnitude)
    }
}

fn is_square(value: &BigInt) -> bool {
    if value.is_negative() {
        return false;
    }
    let root = integer_sqrt(value);
    &root * &root == *value
}

fn atkin_morain_edge(value: &BigInt) -> BigInt {
    let fourth = integer_sqrt(&integer_sqrt(value));
    let shifted = fourth + BigInt::one();
    &shifted * &shifted
}

fn hasse_interval_ok(value: &BigInt, order: &BigInt) -> bool {
    if order <= &BigInt::zero() {
        return false;
    }
    let root = integer_sqrt(value);
    let low = &root - BigInt::one();
    let high = &root + BigInt::one();
    *order > &low * &low && *order < &high * &high
}

fn point_double(point: &ProjectivePoint, curve_a: &BigInt, modulus: &BigInt) -> ProjectivePoint {
    if point.z.rem_euclid(modulus).is_zero() {
        return ProjectivePoint {
            x: BigInt::zero(),
            y: BigInt::one(),
            z: BigInt::zero(),
        };
    }
    let xx = (&point.x * &point.x).rem_euclid(modulus);
    let yy = (&point.y * &point.y).rem_euclid(modulus);
    let yyyy = (&yy * &yy).rem_euclid(modulus);
    let zz = (&point.z * &point.z).rem_euclid(modulus);
    let zzzz = (&zz * &zz).rem_euclid(modulus);
    let sum = (&point.x + &yy).rem_euclid(modulus);
    let doubled = (((&sum * &sum).rem_euclid(modulus) - &xx - &yyyy).rem_euclid(modulus)
        * BigInt::from(2))
    .rem_euclid(modulus);
    let pace = ((BigInt::from(3) * &xx).rem_euclid(modulus)
        + (curve_a * &zzzz).rem_euclid(modulus))
    .rem_euclid(modulus);
    let fore = (&pace * &pace - BigInt::from(2) * &doubled).rem_euclid(modulus);
    let next_x = fore.clone();
    let next_y = (&pace * &(doubled - &fore) - BigInt::from(8) * &yyyy).rem_euclid(modulus);
    let next_z = (BigInt::from(2) * &point.y * &point.z).rem_euclid(modulus);
    ProjectivePoint {
        x: next_x,
        y: next_y,
        z: next_z,
    }
}

fn point_add(
    left: &ProjectivePoint,
    right: &ProjectivePoint,
    curve_a: &BigInt,
    modulus: &BigInt,
) -> ProjectivePoint {
    if left.z.rem_euclid(modulus).is_zero() {
        return right.clone();
    }
    if right.z.rem_euclid(modulus).is_zero() {
        return left.clone();
    }
    let z1z1 = (&left.z * &left.z).rem_euclid(modulus);
    let z2z2 = (&right.z * &right.z).rem_euclid(modulus);
    let u1 = (&left.x * &z2z2).rem_euclid(modulus);
    let u2 = (&right.x * &z1z1).rem_euclid(modulus);
    let s1 = ((&left.y * &right.z).rem_euclid(modulus) * &z2z2).rem_euclid(modulus);
    let s2 = ((&right.y * &left.z).rem_euclid(modulus) * &z1z1).rem_euclid(modulus);
    if u1 == u2 {
        if s1 == s2 {
            return point_double(left, curve_a, modulus);
        }
        return ProjectivePoint {
            x: BigInt::zero(),
            y: BigInt::one(),
            z: BigInt::zero(),
        };
    }
    let held = (&u2 - &u1).rem_euclid(modulus);
    let doubled = (&held * BigInt::from(2)).rem_euclid(modulus);
    let phased = (&doubled * &doubled).rem_euclid(modulus);
    let raced = (&held * &phased).rem_euclid(modulus);
    let paced = (&s2 - &s1).rem_euclid(modulus);
    let paced_two = (&paced * BigInt::from(2)).rem_euclid(modulus);
    let vee = (&u1 * &phased).rem_euclid(modulus);
    let fore = (&paced_two * &paced_two - &raced - BigInt::from(2) * &vee).rem_euclid(modulus);
    let next_x = fore.clone();
    let next_y =
        (&paced_two * &(&vee - &fore) - BigInt::from(2) * &s1 * &raced).rem_euclid(modulus);
    let z_sum = (&left.z + &right.z).rem_euclid(modulus);
    let next_z =
        ((&z_sum * &z_sum - &z1z1 - &z2z2).rem_euclid(modulus) * &held).rem_euclid(modulus);
    ProjectivePoint {
        x: next_x,
        y: next_y,
        z: next_z,
    }
}

fn scalar_multiply(
    point: &ProjectivePoint,
    scalar: &BigInt,
    curve_a: &BigInt,
    modulus: &BigInt,
) -> ProjectivePoint {
    let mut total = ProjectivePoint {
        x: BigInt::zero(),
        y: BigInt::one(),
        z: BigInt::zero(),
    };
    let width = scalar.bits();
    let mut index = width;
    while index > 0 {
        index -= 1;
        total = point_double(&total, curve_a, modulus);
        if scalar.bit(index) {
            total = point_add(&total, point, curve_a, modulus);
        }
    }
    total
}

fn verify_level(level: &EcppLevel) -> bool {
    if level.n < BigInt::from(2) {
        return false;
    }
    if level.q < BigInt::from(2) {
        return false;
    }
    if level.a.is_negative()
        || level.b.is_negative()
        || level.x.is_negative()
        || level.y.is_negative()
        || level.m.is_negative()
        || level.q.is_negative()
    {
        return false;
    }
    let modulus = &level.n;
    let delta = ((BigInt::from(4) * &level.a * &level.a * &level.a).rem_euclid(modulus)
        + (BigInt::from(27) * &level.b * &level.b).rem_euclid(modulus))
    .rem_euclid(modulus);
    if integer_gcd(&delta, modulus) != BigInt::one() {
        return false;
    }
    let lhs = (&level.y * &level.y).rem_euclid(modulus);
    let rhs = ((&level.x * &level.x * &level.x).rem_euclid(modulus)
        + (&level.a * &level.x).rem_euclid(modulus)
        + &level.b)
        .rem_euclid(modulus);
    if lhs != rhs {
        return false;
    }
    if !hasse_interval_ok(modulus, &level.m) {
        return false;
    }
    if (&level.m % &level.q) != BigInt::zero() {
        return false;
    }
    if level.q <= atkin_morain_edge(modulus) {
        return false;
    }
    let base = ProjectivePoint {
        x: level.x.rem_euclid(modulus),
        y: level.y.rem_euclid(modulus),
        z: BigInt::one(),
    };
    let full = scalar_multiply(&base, &level.m, &level.a, modulus);
    if !full.z.rem_euclid(modulus).is_zero() {
        return false;
    }
    let cofactor = &level.m / &level.q;
    let part = scalar_multiply(&base, &cofactor, &level.a, modulus);
    let depth = part.z.rem_euclid(modulus);
    if depth.is_zero() {
        return false;
    }
    integer_gcd(&depth, modulus) == BigInt::one()
}

pub fn verify_ecpp_certificate(cert: &EcppCertificate) -> bool {
    if cert.levels.is_empty() {
        return false;
    }
    if cert.n < BigInt::from(2) {
        return false;
    }
    if cert.levels[0].n != cert.n {
        return false;
    }
    for (index, level) in cert.levels.iter().enumerate() {
        if index + 1 < cert.levels.len() && level.q != cert.levels[index + 1].n {
            return false;
        }
        if !verify_level(level) {
            return false;
        }
    }
    match cert.levels.last() {
        Some(final_level) => match final_level.q.to_u64() {
            Some(small) => is_prime_u64(small),
            None => false,
        },
        None => false,
    }
}

pub fn certificate_to_json(cert: &EcppCertificate) -> String {
    let levels: Vec<serde_json::Value> = cert
        .levels
        .iter()
        .map(|level| {
            serde_json::json!({
                "n": level.n.to_string(),
                "a": level.a.to_string(),
                "b": level.b.to_string(),
                "x": level.x.to_string(),
                "y": level.y.to_string(),
                "m": level.m.to_string(),
                "q": level.q.to_string(),
                "discriminant": level.discriminant.to_string(),
            })
        })
        .collect();
    serde_json::json!({
        "format": ECPP_CERT_FORMAT,
        "n": cert.n.to_string(),
        "levels": levels,
    })
    .to_string()
}

fn schema_error(detail: String) -> crate::error::FactorError {
    crate::error::FactorError::InvalidInput(format!("ecpp certificate schema violation: {detail}"))
}

fn read_decimal_field(
    body: &serde_json::Value,
    level_index: Option<usize>,
    name: &str,
) -> Result<BigInt, crate::error::FactorError> {
    let text = body
        .get(name)
        .and_then(|entry| entry.as_str())
        .ok_or_else(|| {
            schema_error(match level_index {
                Some(position) => format!("levels[{position}].{name} missing or not a string"),
                None => format!("top-level {name} missing or not a string"),
            })
        })?;
    parse_canonical_decimal(text).ok_or_else(|| {
        schema_error(match level_index {
            Some(position) => format!("levels[{position}].{name} not canonical decimal"),
            None => format!("top-level {name} not canonical decimal"),
        })
    })
}

pub fn verify_ecpp_json(value: &BigInt, text: &str) -> Result<bool, crate::error::FactorError> {
    let body: serde_json::Value = serde_json::from_str(text)
        .map_err(|failure| schema_error(format!("not JSON: {failure}")))?;
    let format = body
        .get("format")
        .and_then(|entry| entry.as_str())
        .ok_or_else(|| schema_error("format missing or not a string".to_string()))?;
    if format != ECPP_CERT_FORMAT {
        return Err(schema_error(format!("unexpected format tag: {format}")));
    }
    let claimed = read_decimal_field(&body, None, "n")?;
    if claimed != *value {
        return Err(schema_error(
            "top-level n differs from requested n".to_string(),
        ));
    }
    let rows = body
        .get("levels")
        .and_then(|entry| entry.as_array())
        .ok_or_else(|| schema_error("levels missing or not an array".to_string()))?;
    if rows.is_empty() {
        return Err(schema_error("levels array empty".to_string()));
    }
    let mut levels = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        let discriminant_text = row
            .get("discriminant")
            .and_then(|entry| entry.as_str())
            .ok_or_else(|| {
                schema_error(format!(
                    "levels[{index}].discriminant missing or not a string"
                ))
            })?;
        let discriminant = parse_discriminant_decimal(discriminant_text).ok_or_else(|| {
            schema_error(format!(
                "levels[{index}].discriminant not a negative integer"
            ))
        })?;
        if discriminant >= BigInt::zero() {
            return Err(schema_error(format!(
                "levels[{index}].discriminant not a negative integer"
            )));
        }
        levels.push(EcppLevel {
            n: read_decimal_field(row, Some(index), "n")?,
            a: read_decimal_field(row, Some(index), "a")?,
            b: read_decimal_field(row, Some(index), "b")?,
            x: read_decimal_field(row, Some(index), "x")?,
            y: read_decimal_field(row, Some(index), "y")?,
            m: read_decimal_field(row, Some(index), "m")?,
            q: read_decimal_field(row, Some(index), "q")?,
            discriminant,
        });
    }
    if levels[0].n != claimed {
        return Err(schema_error(
            "levels[0].n differs from top-level n".to_string(),
        ));
    }
    for (index, level) in levels.iter().enumerate() {
        if index + 1 < levels.len() && level.q != levels[index + 1].n {
            return Err(schema_error(format!(
                "levels[{index}].q differs from levels[{}].n",
                index + 1
            )));
        }
    }
    Ok(verify_ecpp_certificate(&EcppCertificate {
        n: claimed,
        levels,
    }))
}

fn kronecker_symbol(top: &BigInt, bottom: &BigInt) -> i32 {
    if bottom <= &BigInt::one() {
        return 0;
    }
    if bottom % BigInt::from(2) == BigInt::zero() {
        return 0;
    }
    let mut head = top.rem_euclid(bottom);
    let mut tail = bottom.clone();
    let mut sign = 1_i32;
    let zero = BigInt::zero();
    let one = BigInt::one();
    let two = BigInt::from(2);
    loop {
        if tail == one {
            return sign;
        }
        if head == zero {
            return 0;
        }
        while (&head % &two) == zero {
            head /= &two;
            let residue = &tail % BigInt::from(8);
            if residue == BigInt::from(3) || residue == BigInt::from(5) {
                sign = -sign;
            }
        }
        std::mem::swap(&mut head, &mut tail);
        if (&head % BigInt::from(4)) == BigInt::from(3)
            && (&tail % BigInt::from(4)) == BigInt::from(3)
        {
            sign = -sign;
        }
        head = head.rem_euclid(&tail);
    }
}

fn miller_rabin_filter(value: &BigInt) -> bool {
    let zero = BigInt::zero();
    let one = BigInt::one();
    let two = BigInt::from(2);
    if value < &two {
        return false;
    }
    for prime in small_filter_primes() {
        let small = BigInt::from(*prime);
        if value == &small {
            return true;
        }
        if (value % &small) == zero {
            return false;
        }
    }
    if (value % &two) == zero {
        return false;
    }
    let mut odd = value - &one;
    let mut shifts = 0_u32;
    while (&odd % &two) == zero {
        odd /= &two;
        shifts += 1;
    }
    for base in MILLER_RABIN_BASES {
        let trial = BigInt::from(base) % value;
        if trial == zero {
            continue;
        }
        let mut root = crate::integer::mod_pow(&trial, &odd, value);
        if root == one || root == value - &one {
            continue;
        }
        let mut passed = false;
        for _ in 1..shifts {
            root = (&root * &root) % value;
            if root == value - &one {
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

fn modular_sqrt(value: &BigInt, prime: &BigInt) -> Option<BigInt> {
    let one = BigInt::one();
    let two = BigInt::from(2);
    if prime == &two {
        return Some(value.rem_euclid(prime));
    }
    if prime < &two || (prime % &two) == BigInt::zero() {
        return None;
    }
    let residue = value.rem_euclid(prime);
    if residue.is_zero() {
        return Some(BigInt::zero());
    }
    let half = (prime - &one) / &two;
    if crate::integer::mod_pow(&residue, &half, prime) != one {
        return None;
    }
    if (prime % BigInt::from(4)) == BigInt::from(3) {
        let exponent = (prime + &one) / BigInt::from(4);
        return Some(crate::integer::mod_pow(&residue, &exponent, prime));
    }
    let mut quotient = prime - &one;
    let mut power = 0_u32;
    while (&quotient % &two) == BigInt::zero() {
        quotient /= &two;
        power += 1;
    }
    let mut probe = two.clone();
    let mut guard = 0_u64;
    while crate::integer::mod_pow(&probe, &half, prime) != prime - &one {
        probe += &one;
        guard += 1;
        if guard > NON_RESIDUE_TRIES {
            return None;
        }
    }
    let mut factor = crate::integer::mod_pow(&probe, &quotient, prime);
    let shifted = (&quotient + &one) / &two;
    let mut root = crate::integer::mod_pow(&residue, &shifted, prime);
    let mut widgets = crate::integer::mod_pow(&residue, &quotient, prime);
    let mut depth = power;
    while widgets != one {
        let mut distance = 1_u32;
        let mut squared = (&widgets * &widgets).rem_euclid(prime);
        while distance < depth && squared != one {
            squared = (&squared * &squared).rem_euclid(prime);
            distance += 1;
        }
        if distance == depth {
            return None;
        }
        let exponent = BigInt::one() << (depth - distance - 1);
        let turn = crate::integer::mod_pow(&factor, &exponent, prime);
        root = (&root * &turn).rem_euclid(prime);
        factor = (&turn * &turn).rem_euclid(prime);
        widgets = (&widgets * &factor).rem_euclid(prime);
        depth = distance;
    }
    Some(root)
}

fn cornacchia_traces(value: &BigInt, root: &BigInt, absolute: &BigInt) -> Vec<BigInt> {
    let four_n = BigInt::from(4) * value;
    let mut remainders = Vec::new();
    let mut head = value.clone();
    let mut tail = root.rem_euclid(value);
    while !tail.is_zero() {
        let next = (&head % &tail).rem_euclid(&tail);
        head = tail;
        tail = next;
        if &head * &head < four_n && !head.is_zero() && !remainders.contains(&head) {
            remainders.push(head.clone());
        }
    }
    let mut traces = Vec::new();
    for base in remainders {
        for scaled in [base.clone(), BigInt::from(2) * &base] {
            if &scaled * &scaled > four_n {
                continue;
            }
            let gap = &four_n - &scaled * &scaled;
            if (&gap % absolute) != BigInt::zero() {
                continue;
            }
            let quotient = gap / absolute;
            if is_square(&quotient) && !traces.contains(&scaled) {
                traces.push(scaled);
            }
        }
    }
    traces
}

fn trace_orders(value: &BigInt, trace: &BigInt, root_u: &BigInt, discriminant: i64) -> Vec<BigInt> {
    let one = BigInt::one();
    let two = BigInt::from(2);
    let three = BigInt::from(3);
    let mut orders = Vec::new();
    let mut push_trace = |held: BigInt| {
        let order = value + &one - &held;
        if !orders.contains(&order) {
            orders.push(order);
        }
    };
    match discriminant {
        -3 => {
            push_trace(trace.clone());
            push_trace(-trace.clone());
            let plus = trace + &three * root_u;
            let minus = trace - &three * root_u;
            for signed in [&plus, &minus] {
                for version in [signed.clone(), -signed.clone()] {
                    if version.rem_euclid(&two) == BigInt::zero() {
                        push_trace(&version / &two);
                    }
                }
            }
        }
        -4 => {
            if trace.rem_euclid(&two) == BigInt::zero() {
                let half = trace / &two;
                push_trace(BigInt::from(2) * &half);
                push_trace(-BigInt::from(2) * &half);
                push_trace(BigInt::from(2) * root_u);
                push_trace(-BigInt::from(2) * root_u);
            }
        }
        _ => {
            push_trace(trace.clone());
            push_trace(-trace.clone());
        }
    }
    orders
}

fn strip_small_factors(value: &BigInt) -> (BigInt, BigInt) {
    let mut smooth = BigInt::one();
    let mut rest = value.clone();
    for prime in strip_primes() {
        let small = BigInt::from(*prime);
        if &small * &small > rest {
            break;
        }
        while (&rest % &small) == BigInt::zero() {
            rest /= &small;
            smooth *= &small;
        }
        if rest == BigInt::one() {
            break;
        }
    }
    (smooth, rest)
}

fn rho_divisor(value: &BigInt) -> Option<BigInt> {
    let zero = BigInt::zero();
    let one = BigInt::one();
    let two = BigInt::from(2);
    if value <= &one {
        return None;
    }
    if (value % &two) == zero {
        return Some(two);
    }
    for prime in small_filter_primes() {
        let small = BigInt::from(*prime);
        if value == &small {
            return None;
        }
        if (value % &small) == zero {
            return Some(small);
        }
    }
    for round in 1..=RHO_ROUNDS {
        let constant = BigInt::from(round);
        let mut first = two.clone();
        let mut second = two.clone();
        let mut block = one.clone();
        let mut spent = 0_u64;
        let mut divisor = one.clone();
        while divisor == one && spent < RHO_EVAL_BUDGET {
            first = ((&first * &first).rem_euclid(value) + &constant).rem_euclid(value);
            second = ((&second * &second).rem_euclid(value) + &constant).rem_euclid(value);
            second = ((&second * &second).rem_euclid(value) + &constant).rem_euclid(value);
            let gap = if first > second {
                &first - &second
            } else {
                &second - &first
            };
            if gap.is_zero() {
                break;
            }
            block = (&block * &gap).rem_euclid(value);
            spent += 1;
            if spent.is_multiple_of(256) {
                divisor = integer_gcd(&block, value);
                if divisor == value.clone() {
                    break;
                }
                block = one.clone();
            }
        }
        if divisor == one && block != one {
            divisor = integer_gcd(&block, value);
        }
        if divisor != one && divisor != *value {
            return Some(divisor);
        }
    }
    None
}

fn collect_prime_factors(value: &BigInt, found: &mut Vec<BigInt>, depth: u32) {
    if value <= &BigInt::one() || depth > 16 {
        return;
    }
    let (_, rest) = strip_small_factors(value);
    if rest == BigInt::one() {
        return;
    }
    if miller_rabin_filter(&rest) {
        found.push(rest);
        return;
    }
    if let Some(divisor) = rho_divisor(&rest) {
        let other = &rest / &divisor;
        collect_prime_factors(&divisor, found, depth + 1);
        collect_prime_factors(&other, found, depth + 1);
    }
}

fn hunt_quotient(modulus: &BigInt, order: &BigInt) -> Option<BigInt> {
    if !hasse_interval_ok(modulus, order) {
        return None;
    }
    let edge = atkin_morain_edge(modulus);
    let mut factors = Vec::new();
    collect_prime_factors(order, &mut factors, 0);
    let mut best: Option<BigInt> = None;
    for factor in factors {
        if factor > edge && factor < *modulus {
            let replace = match &best {
                Some(current) => factor < *current,
                None => true,
            };
            if replace {
                best = Some(factor);
            }
        }
    }
    best
}

fn find_non_residue(modulus: &BigInt) -> Option<BigInt> {
    let one = BigInt::one();
    let exponent = (modulus - &one) / BigInt::from(2);
    let target = modulus - &one;
    let mut probe = BigInt::from(2);
    for _ in 0..NON_RESIDUE_TRIES {
        if crate::integer::mod_pow(&probe, &exponent, modulus) == target {
            return Some(probe);
        }
        probe += &one;
    }
    None
}

fn base_curve(joint: &BigInt, modulus: &BigInt) -> (BigInt, BigInt) {
    let residue = joint.rem_euclid(modulus);
    if residue.is_zero() {
        return (BigInt::zero(), BigInt::one());
    }
    if residue == BigInt::from(1728) {
        return (BigInt::one(), BigInt::zero());
    }
    let spread = (BigInt::from(1728) - &residue).rem_euclid(modulus);
    let curve_a = (BigInt::from(3) * &residue * &spread).rem_euclid(modulus);
    let curve_b =
        ((BigInt::from(2) * &residue).rem_euclid(modulus) * &spread * &spread).rem_euclid(modulus);
    (curve_a, curve_b)
}

fn twist_curves(
    curve_a: &BigInt,
    curve_b: &BigInt,
    joint: &BigInt,
    modulus: &BigInt,
    stir: &BigInt,
) -> Vec<(BigInt, BigInt)> {
    let residue = joint.rem_euclid(modulus);
    let mut curves = Vec::new();
    if residue.is_zero() {
        let mut power = BigInt::one();
        for _ in 0..6 {
            curves.push((BigInt::zero(), power.clone()));
            power = (&power * stir).rem_euclid(modulus);
        }
        return curves;
    }
    if residue == BigInt::from(1728) {
        let mut power = BigInt::one();
        for _ in 0..4 {
            curves.push((power.clone(), BigInt::zero()));
            power = (&power * stir).rem_euclid(modulus);
        }
        return curves;
    }
    curves.push((curve_a.clone(), curve_b.clone()));
    let square = (stir * stir).rem_euclid(modulus);
    let cubed = (&square * stir).rem_euclid(modulus);
    curves.push((
        (curve_a * &square).rem_euclid(modulus),
        (curve_b * &cubed).rem_euclid(modulus),
    ));
    curves
}

fn build_level(
    modulus: &BigInt,
    discriminant: i64,
    joint: &BigInt,
    order: &BigInt,
    quotient: &BigInt,
) -> Option<EcppLevel> {
    let stir = find_non_residue(modulus)?;
    let (curve_a, curve_b) = base_curve(joint, modulus);
    let cofactor = order / quotient;
    for (held_a, held_b) in twist_curves(&curve_a, &curve_b, joint, modulus, &stir) {
        let delta = ((BigInt::from(4) * &held_a * &held_a * &held_a).rem_euclid(modulus)
            + (BigInt::from(27) * &held_b * &held_b).rem_euclid(modulus))
        .rem_euclid(modulus);
        if integer_gcd(&delta, modulus) != BigInt::one() {
            continue;
        }
        let mut probe_x = BigInt::zero();
        for _ in 0..POINT_X_TRIES {
            let rhs = ((&probe_x * &probe_x * &probe_x).rem_euclid(modulus)
                + (&held_a * &probe_x).rem_euclid(modulus)
                + &held_b)
                .rem_euclid(modulus);
            probe_x += BigInt::one();
            let probe_y = match modular_sqrt(&rhs, modulus) {
                Some(root) => root,
                None => continue,
            };
            let point = ProjectivePoint {
                x: (&probe_x - BigInt::one()).rem_euclid(modulus),
                y: probe_y,
                z: BigInt::one(),
            };
            let full = scalar_multiply(&point, order, &held_a, modulus);
            if !full.z.rem_euclid(modulus).is_zero() {
                continue;
            }
            let part = scalar_multiply(&point, &cofactor, &held_a, modulus);
            let depth = part.z.rem_euclid(modulus);
            if depth.is_zero() {
                continue;
            }
            if integer_gcd(&depth, modulus) != BigInt::one() {
                continue;
            }
            return Some(EcppLevel {
                n: modulus.clone(),
                a: held_a,
                b: held_b,
                x: point.x,
                y: point.y,
                m: order.clone(),
                q: quotient.clone(),
                discriminant: BigInt::from(discriminant),
            });
        }
    }
    None
}

fn find_levels(modulus: &BigInt) -> Vec<EcppLevel> {
    let mut scored: Vec<(BigInt, EcppLevel)> = Vec::new();
    for (discriminant, joint) in HEEGNER_ROWS {
        let absolute = BigInt::from(discriminant.unsigned_abs());
        if integer_gcd(&absolute, modulus) != BigInt::one() {
            continue;
        }
        if kronecker_symbol(&BigInt::from(discriminant), modulus) != 1 {
            continue;
        }
        let stained = BigInt::from(discriminant).rem_euclid(modulus);
        let root = match modular_sqrt(&stained, modulus) {
            Some(value) => value,
            None => continue,
        };
        for trace in cornacchia_traces(modulus, &root, &absolute) {
            let four_n = BigInt::from(4) * modulus;
            let square_gap = &four_n - &trace * &trace;
            if (&square_gap % &absolute) != BigInt::zero() {
                continue;
            }
            let root_u = integer_sqrt(&(square_gap / &absolute));
            for order in trace_orders(modulus, &trace, &root_u, discriminant) {
                if !hasse_interval_ok(modulus, &order) {
                    continue;
                }
                let quotient = match hunt_quotient(modulus, &order) {
                    Some(value) => value,
                    None => continue,
                };
                let joint_big = BigInt::from(joint);
                match build_level(modulus, discriminant, &joint_big, &order, &quotient) {
                    Some(level) => {
                        scored.push((level.q.clone(), level));
                        scored.sort_by(|left, right| left.0.cmp(&right.0));
                        scored.truncate(LEVEL_CANDIDATE_CAP);
                        if let Some((best, _)) = scored.first()
                            && best.bits() * 2 <= modulus.bits() + 8
                        {
                            return scored.into_iter().map(|(_, level)| level).collect();
                        }
                    }
                    None => continue,
                }
            }
        }
        scored.sort_by(|left, right| left.0.cmp(&right.0));
        scored.truncate(LEVEL_CANDIDATE_CAP);
        if let Some((best, _)) = scored.first()
            && best.bits() * 2 <= modulus.bits() + 8
        {
            break;
        }
    }
    scored.sort_by(|left, right| left.0.cmp(&right.0));
    scored.truncate(LEVEL_CANDIDATE_CAP);
    scored.into_iter().map(|(_, level)| level).collect()
}

fn prove_from(current: &BigInt, levels: &mut Vec<EcppLevel>, depth: u32, top: bool) -> bool {
    if depth > MAX_CHAIN_DEPTH {
        return false;
    }
    if !top && let Some(small) = current.to_u64() {
        return is_prime_u64(small);
    }
    for candidate in find_levels(current) {
        let next = candidate.q.clone();
        levels.push(candidate);
        if prove_from(&next, levels, depth + 1, false) {
            return true;
        }
        levels.pop();
    }
    false
}

pub fn prove_ecpp(value: &BigInt) -> Option<EcppCertificate> {
    if value < &BigInt::from(2) {
        return None;
    }
    if value.bits() < ECPP_MIN_BITS || value.bits() > ECPP_MAX_BITS {
        return None;
    }
    for prime in small_filter_primes() {
        if (value % BigInt::from(*prime)) == BigInt::zero() {
            return None;
        }
    }
    let mut levels = Vec::new();
    if !prove_from(value, &mut levels, 0, true) {
        return None;
    }
    let cert = EcppCertificate {
        n: value.clone(),
        levels,
    };
    if verify_ecpp_certificate(&cert) {
        Some(cert)
    } else {
        None
    }
}

#[cfg(test)]
mod proof_tests {
    use super::*;
    use crate::integer::mod_pow as big_mod_pow;
    use crate::prime::is_prime_u64;

    fn decimal(text: &str) -> BigInt {
        BigInt::from_str_radix(text, 10).unwrap()
    }

    #[test]
    fn canonical_decimal_gate() {
        assert_eq!(parse_canonical_decimal("0"), Some(BigInt::zero()));
        assert_eq!(parse_canonical_decimal("101"), Some(BigInt::from(101)));
        assert_eq!(parse_canonical_decimal(""), None);
        assert_eq!(parse_canonical_decimal("007"), None);
        assert_eq!(parse_canonical_decimal("00"), None);
        assert_eq!(parse_canonical_decimal("-5"), None);
        assert_eq!(parse_canonical_decimal(" 12"), None);
        assert_eq!(parse_canonical_decimal("12 "), None);
        assert_eq!(parse_canonical_decimal("1_000"), None);
        assert_eq!(parse_canonical_decimal("abc"), None);
        assert_eq!(parse_discriminant_decimal("-3"), Some(BigInt::from(-3)));
        assert_eq!(parse_discriminant_decimal("3"), Some(BigInt::from(3)));
        assert_eq!(parse_discriminant_decimal("0"), None);
        assert_eq!(parse_discriminant_decimal("-0"), None);
        assert_eq!(parse_discriminant_decimal(""), None);
    }

    #[test]
    fn kronecker_matches_known_values() {
        assert_eq!(kronecker_symbol(&BigInt::from(-3), &BigInt::from(7)), 1);
        assert_eq!(kronecker_symbol(&BigInt::from(-3), &BigInt::from(5)), -1);
        assert_eq!(kronecker_symbol(&BigInt::from(5), &BigInt::from(5)), 0);
        assert_eq!(kronecker_symbol(&BigInt::from(2), &BigInt::from(7)), 1);
        assert_eq!(kronecker_symbol(&BigInt::from(2), &BigInt::from(3)), -1);
        assert_eq!(kronecker_symbol(&BigInt::from(1), &BigInt::from(101)), 1);
        assert_eq!(kronecker_symbol(&BigInt::from(-4), &BigInt::from(8)), 0);
    }

    #[test]
    fn sqrt_matches_euler_criterion() {
        for (value, prime) in [(2_i64, 7_i64), (3, 11), (5, 101), (10, 1009)] {
            let big_value = BigInt::from(value);
            let big_prime = BigInt::from(prime);
            match modular_sqrt(&big_value, &big_prime) {
                Some(root) => {
                    assert_eq!(
                        (&root * &root).rem_euclid(&big_prime),
                        big_value.rem_euclid(&big_prime)
                    );
                    assert_eq!(
                        big_mod_pow(
                            &big_value,
                            &((&big_prime - BigInt::one()) / BigInt::from(2)),
                            &big_prime
                        ),
                        BigInt::one()
                    );
                }
                None => {
                    assert_ne!(
                        big_mod_pow(
                            &big_value,
                            &((&big_prime - BigInt::one()) / BigInt::from(2)),
                            &big_prime
                        ),
                        BigInt::one()
                    );
                }
            }
        }
        let root_two = modular_sqrt(&BigInt::from(2), &BigInt::from(7)).expect("root");
        assert_eq!(
            (&root_two * &root_two).rem_euclid(&BigInt::from(7)),
            BigInt::from(2)
        );
        assert_eq!(modular_sqrt(&BigInt::from(3), &BigInt::from(7)), None);
    }

    #[test]
    fn miller_rabin_filter_matches_u64_gate_below_5000() {
        for small in 0_u64..5000 {
            assert_eq!(
                miller_rabin_filter(&BigInt::from(small)),
                is_prime_u64(small),
                "filter split at {small}"
            );
        }
    }

    #[test]
    fn projective_matches_affine_on_small_prime() {
        let modulus = BigInt::from(1009);
        let curve_a = BigInt::from(37);
        let curve_b = BigInt::from(41);
        let mut checked = 0_u32;
        let mut probe_x = BigInt::zero();
        while checked < 4 {
            let rhs = ((&probe_x * &probe_x * &probe_x).rem_euclid(&modulus)
                + (&curve_a * &probe_x).rem_euclid(&modulus)
                + &curve_b)
                .rem_euclid(&modulus);
            probe_x += BigInt::one();
            let probe_y = match modular_sqrt(&rhs, &modulus) {
                Some(root) => root,
                None => continue,
            };
            let affine_x = (&probe_x - BigInt::one()).rem_euclid(&modulus);
            for scalar in [1_u64, 2, 7, 123, 1009] {
                let point = ProjectivePoint {
                    x: affine_x.clone(),
                    y: probe_y.clone(),
                    z: BigInt::one(),
                };
                let got = scalar_multiply(&point, &BigInt::from(scalar), &curve_a, &modulus);
                let want = affine_multiply(&affine_x, &probe_y, scalar, &curve_a, &modulus);
                if got.z.rem_euclid(&modulus).is_zero() {
                    assert!(want.is_none(), "projective infinity must match affine");
                } else {
                    let (want_x, want_y) = want.expect("affine point present");
                    let inv_z = big_mod_pow(&got.z, &(modulus.clone() - BigInt::from(2)), &modulus);
                    let inv_z2 = (&inv_z * &inv_z).rem_euclid(&modulus);
                    let got_x = (&got.x * &inv_z2).rem_euclid(&modulus);
                    let got_y = (&got.y * &inv_z * &inv_z2).rem_euclid(&modulus);
                    assert_eq!(got_x, want_x);
                    assert_eq!(got_y, want_y);
                }
            }
            checked += 1;
        }
    }

    fn affine_multiply(
        base_x: &BigInt,
        base_y: &BigInt,
        scalar: u64,
        curve_a: &BigInt,
        modulus: &BigInt,
    ) -> Option<(BigInt, BigInt)> {
        let mut acc: Option<(BigInt, BigInt)> = None;
        let mut head = Some((base_x.clone(), base_y.clone()));
        let mut rest = scalar;
        while rest > 0 {
            if rest & 1 == 1 {
                acc = match (acc, head.clone()) {
                    (Some(pair), Some(top)) => {
                        affine_add(&pair.0, &pair.1, &top.0, &top.1, curve_a, modulus)
                    }
                    (None, top) => top,
                    (pair, None) => pair,
                };
            }
            head = match head {
                Some(top) => affine_add(&top.0, &top.1, &top.0, &top.1, curve_a, modulus),
                None => None,
            };
            rest >>= 1;
        }
        acc
    }

    fn affine_add(
        x1: &BigInt,
        y1: &BigInt,
        x2: &BigInt,
        y2: &BigInt,
        curve_a: &BigInt,
        modulus: &BigInt,
    ) -> Option<(BigInt, BigInt)> {
        if x1 == x2 {
            if ((y1 + y2).rem_euclid(modulus)).is_zero() {
                return None;
            }
            let slope = ((BigInt::from(3) * x1 * x1 + curve_a).rem_euclid(modulus)
                * big_mod_pow(
                    &(BigInt::from(2) * y1).rem_euclid(modulus),
                    &(modulus - BigInt::from(2)),
                    modulus,
                ))
            .rem_euclid(modulus);
            let out_x = (&slope * &slope - x1 - x2).rem_euclid(modulus);
            let out_y = (&slope * &(x1 - &out_x) - y1).rem_euclid(modulus);
            return Some((out_x, out_y));
        }
        let slope = (((y2 - y1).rem_euclid(modulus))
            * big_mod_pow(
                &(x2 - x1).rem_euclid(modulus),
                &(modulus - BigInt::from(2)),
                modulus,
            ))
        .rem_euclid(modulus);
        let out_x = (&slope * &slope - x1 - x2).rem_euclid(modulus);
        let out_y = (&slope * &(x1 - &out_x) - y1).rem_euclid(modulus);
        Some((out_x, out_y))
    }

    fn count_curve_order(curve_a: &BigInt, curve_b: &BigInt, prime: &BigInt) -> BigInt {
        let mut total = BigInt::one();
        let mut probe_x = BigInt::zero();
        while probe_x < *prime {
            let rhs = ((&probe_x * &probe_x * &probe_x).rem_euclid(prime)
                + (curve_a * &probe_x).rem_euclid(prime)
                + curve_b)
                .rem_euclid(prime);
            if rhs.is_zero() {
                total += BigInt::one();
            } else {
                let half = (prime - BigInt::one()) / BigInt::from(2);
                if big_mod_pow(&rhs, &half, prime) == BigInt::one() {
                    total += BigInt::from(2);
                }
            }
            probe_x += BigInt::one();
        }
        total
    }

    #[test]
    fn prover_order_matches_brute_force_count() {
        for text in ["101", "1009"] {
            let prime = decimal(text);
            let mut levels = Vec::new();
            assert!(
                prove_from(&prime, &mut levels, 0, true),
                "small prime must chain: {text}"
            );
            assert!(!levels.is_empty());
            let first = &levels[0];
            let counted = count_curve_order(&first.a, &first.b, &prime);
            assert_eq!(first.m, counted, "order mismatch at {text}");
            let cert = EcppCertificate {
                n: prime.clone(),
                levels,
            };
            assert!(verify_ecpp_certificate(&cert));
        }
    }

    #[test]
    fn tamper_battery_rejects_single_field_mutants() {
        let prime = decimal("9223372036854775837");
        let cert = prove_ecpp(&prime).expect("corpus prime must certify");
        assert!(verify_ecpp_certificate(&cert));
        let base = certificate_to_json(&cert);
        assert!(verify_ecpp_json(&prime, &base).expect("schema must hold"));
        let fields = ["a", "b", "x", "y", "m", "q"];
        for field in fields {
            let mut body: serde_json::Value = serde_json::from_str(&base).expect("parses");
            let bumped =
                decimal(body["levels"][0][field].as_str().expect("field present")) + BigInt::one();
            body["levels"][0][field] = serde_json::Value::String(bumped.to_string());
            let text = body.to_string();
            assert!(
                !verify_ecpp_json(&prime, &text).unwrap_or(false),
                "mutant {field} must not verify"
            );
        }
        let mut format_body: serde_json::Value = serde_json::from_str(&base).expect("parses");
        format_body["format"] = serde_json::Value::String("graynfs-ecpp-cert/0".to_string());
        assert!(verify_ecpp_json(&prime, &format_body.to_string()).is_err());
        let mut link_body: serde_json::Value = serde_json::from_str(&base).expect("parses");
        if link_body["levels"].as_array().expect("array").len() > 1 {
            let bumped =
                decimal(link_body["levels"][1]["n"].as_str().expect("present")) + BigInt::one();
            link_body["levels"][1]["n"] = serde_json::Value::String(bumped.to_string());
            assert!(verify_ecpp_json(&prime, &link_body.to_string()).is_err());
        }
        let other = decimal("604462909807314587353111");
        assert!(verify_ecpp_json(&other, &base).is_err());
        assert!(verify_ecpp_json(&prime, "{not json").is_err());
        assert!(verify_ecpp_json(&prime, "{}").is_err());
        let mut lead_body: serde_json::Value = serde_json::from_str(&base).expect("parses");
        lead_body["levels"][0]["a"] = serde_json::Value::String("007".to_string());
        assert!(verify_ecpp_json(&prime, &lead_body.to_string()).is_err());
        let mut disc_body: serde_json::Value = serde_json::from_str(&base).expect("parses");
        disc_body["levels"][0]["discriminant"] = serde_json::Value::String("3".to_string());
        assert!(verify_ecpp_json(&prime, &disc_body.to_string()).is_err());
    }

    #[test]
    fn verifier_rejects_degenerate_quotients_without_panic() {
        let prime = decimal("9223372036854775837");
        let cert = prove_ecpp(&prime).expect("corpus prime must certify");
        let mut zero_q = cert.levels[0].clone();
        zero_q.q = BigInt::zero();
        assert!(!verify_level(&zero_q));
        let mut one_q = cert.levels[0].clone();
        one_q.q = BigInt::one();
        assert!(!verify_level(&one_q));
        let text = certificate_to_json(&cert);
        let mut body: serde_json::Value = serde_json::from_str(&text).expect("parses");
        body["levels"][0]["q"] = serde_json::Value::String("0".to_string());
        assert!(!verify_ecpp_json(&prime, &body.to_string()).expect("schema holds"));
    }

    #[test]
    fn scope_gates_reject_out_of_range_inputs() {
        assert!(prove_ecpp(&BigInt::zero()).is_none());
        assert!(prove_ecpp(&BigInt::one()).is_none());
        assert!(prove_ecpp(&BigInt::from(101)).is_none());
        assert!(prove_ecpp(&(BigInt::one() << 63)).is_none());
        assert!(prove_ecpp(&(BigInt::one() << 129)).is_none());
        assert!(prove_ecpp(&(BigInt::one() << 128)).is_none());
    }
}
