use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use num_traits::ToPrimitive;

use crate::algebraic::QuadraticElement;
use crate::base::FactorBase;
use crate::error::FactorError;
use crate::factor::{ProgressEvent, Timer, TimingCollector};
use crate::integer::{Integer, IntegerExt, integer_gcd};
use crate::linear_algebra::rank_and_dependencies;
use crate::options::{
    FactorOptions, estimate_collect_bytes, estimate_live_bytes, estimate_total_bytes,
};
use crate::polynomial::PolynomialPair;
use crate::prime::{is_prime_u64_deterministic, legendre, mod_sqrt, sieve_primes};
use crate::relations::Relation;

pub(crate) struct SiqsPolynomial {
    pub(crate) coeff_a: Integer,
    pub(crate) coeff_b: Integer,
    pub(crate) coeff_c: Integer,
    pub(crate) prime_factors: Vec<u64>,
}

impl SiqsPolynomial {
    pub(crate) fn evaluate(&self, x: &Integer) -> Integer {
        &self.coeff_a * x * x + Integer::from(2) * &self.coeff_b * x + &self.coeff_c
    }

    pub(crate) fn x_value(&self, x: &Integer) -> Integer {
        &self.coeff_a * x + &self.coeff_b
    }

    /// Test-only oracle for the cached path: direct per-poly root computation.
    #[cfg(test)]
    pub(crate) fn roots_mod(&self, prime: u64) -> Vec<u64> {
        self.roots_mod_with_gamma(
            prime,
            mod_sqrt(
                n_mod_for_gamma(&self.coeff_a, &self.coeff_b, &self.coeff_c, prime),
                prime,
            ),
        )
    }

    /// Root computation with a precomputed `gamma = sqrt(n) mod prime` (X1a).
    /// Identical to [`roots_mod`](Self::roots_mod): the discriminant
    /// `b*b - a*c` equals `n` exactly as an integer (exact division is
    /// enforced at polynomial construction), so `disc mod prime` is
    /// `n mod prime` and its square root is the cached `gamma`.
    /// The `a_mod == 0` linear path never consults `gamma` and is unchanged.
    pub(crate) fn roots_mod_with_gamma(&self, prime: u64, gamma: Option<u64>) -> Vec<u64> {
        let modulus = Integer::from(prime);
        let a_mod = (&self.coeff_a % &modulus).to_u64_saturating();
        let b_mod = (&self.coeff_b % &modulus).to_u64_saturating();
        if a_mod == 0 {
            let denom = (2 * b_mod) % prime;
            if denom == 0 {
                return Vec::new();
            }
            let c_mod = ((&self.coeff_c % &modulus) + &modulus) % &modulus;
            let inv = mod_inverse_u64(denom, prime);
            return vec![((prime - c_mod.to_u64_saturating()) % prime * inv) % prime];
        }
        let Some(gamma) = gamma else {
            return Vec::new();
        };
        let inv_a = mod_inverse_u64(a_mod, prime);
        let neg_b = (prime - b_mod) % prime;
        let mut roots = vec![((neg_b + gamma) % prime * inv_a) % prime];
        let other = ((neg_b + prime - gamma) % prime * inv_a) % prime;
        if other != roots[0] {
            roots.push(other);
        }
        roots
    }
}

/// Residue whose square root is the SIQS `gamma` for `prime`: `n mod prime`
/// recovered from the polynomial identity `b*b - a*c == n`.
#[cfg(test)]
fn n_mod_for_gamma(coeff_a: &Integer, coeff_b: &Integer, coeff_c: &Integer, prime: u64) -> u64 {
    let modulus = Integer::from(prime);
    let disc = (coeff_b * coeff_b - coeff_a * coeff_c) % &modulus;
    disc.to_u64_saturating()
}

/// Per-attempt table of `gamma = sqrt(n) mod prime` (X1a): one `mod_sqrt`
/// per base prime per attempt instead of one per (polynomial, prime).
/// Primes beyond the table length (cannot occur for base primes by
/// construction) fall back to direct computation via `None`-handling callers.
#[derive(Clone, Debug)]
pub(crate) struct SiqsGammaTable {
    gammas: Vec<Option<u64>>,
}

impl SiqsGammaTable {
    pub(crate) fn build(n: &Integer, primes: &[u64]) -> Self {
        let max = primes.iter().copied().max().unwrap_or(0);
        let mut gammas = vec![None; max as usize + 1];
        let mut done = std::collections::HashSet::new();
        for &prime in primes {
            if prime < 2 || !done.insert(prime) {
                continue;
            }
            let residue = (n % Integer::from(prime)).to_u64_saturating();
            gammas[prime as usize] = mod_sqrt(residue, prime);
        }
        Self { gammas }
    }

    pub(crate) fn get(&self, prime: u64) -> Option<u64> {
        self.gammas.get(prime as usize).copied().flatten()
    }

    /// Direct query used only by tests to distinguish "absent" from "non-residue".
    #[cfg(test)]
    fn raw(&self, prime: u64) -> Option<Option<u64>> {
        self.gammas.get(prime as usize).copied()
    }
}

fn mod_inverse_u64(value: u64, modulus: u64) -> u64 {
    let mut old_r = modulus as i128;
    let mut r = value as i128;
    let mut old_s = 0_i128;
    let mut s = 1_i128;
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
    }
    ((old_s % modulus as i128 + modulus as i128) % modulus as i128) as u64
}

fn egcd_inverse(value: &Integer, modulus: &Integer) -> Option<Integer> {
    let mut old_r = modulus.clone();
    let mut r = value % modulus;
    let mut old_s = Integer::from(0);
    let mut s = Integer::from(1);
    while r != Integer::from(0) {
        let quotient = &old_r / &r;
        let next_r = &old_r - &quotient * &r;
        old_r = r;
        r = next_r;
        let next_s = &old_s - &quotient * &s;
        old_s = s;
        s = next_s;
    }
    if old_r != Integer::from(1) && old_r != Integer::from(-1) {
        return None;
    }
    let mut result = old_s % modulus;
    if result < Integer::from(0) {
        result += modulus;
    }
    Some(result)
}

/// One SIQS-suitable prime with its cached square root (X1d): the root of
/// `n` modulo `prime`, computed once per attempt instead of once per
/// generated polynomial.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SiqsSuitablePrime {
    pub(crate) prime: u64,
    pub(crate) gamma: u64,
}

/// Suitable-prime list built once per attempt (X1d). The prime sequence is
/// the exact prefix the per-index builder would produce for any `skip` cap
/// (sieve order and filters are deterministic), so windowing below yields
/// the identical polynomial sequence.
pub(crate) fn suitable_primes(n: &Integer, algebraic_bound: usize) -> Vec<SiqsSuitablePrime> {
    let mut suitable = Vec::new();
    if n <= &Integer::from(1) {
        return suitable;
    }
    for prime in sieve_primes(algebraic_bound) {
        if prime < 3 || n % prime == Integer::from(0) {
            continue;
        }
        let n_mod = (n % Integer::from(prime)).to_u64_saturating();
        if legendre(n_mod, prime) != 1 {
            continue;
        }
        let Some(gamma) = mod_sqrt(n_mod, prime) else {
            continue;
        };
        suitable.push(SiqsSuitablePrime { prime, gamma });
    }
    suitable
}

/// Legacy per-index constructor, kept as the test oracle for the hoisted
/// sequence: sieve the suitable list fresh, then delegate.
#[cfg(test)]
pub(crate) fn generate_polynomial(
    n: &Integer,
    algebraic_bound: usize,
    target_a: &Integer,
    skip: usize,
) -> Option<SiqsPolynomial> {
    let suitable = suitable_primes(n, algebraic_bound);
    generate_polynomial_from(n, &suitable, target_a, skip, None)
}

/// Polynomial generation from a prebuilt suitable list (X1d): window
/// selection, greedy factor take, and CRT combination are the legacy
/// algorithm verbatim, so outputs are identical for every `skip`.
pub(crate) fn generate_polynomial_from(
    n: &Integer,
    suitable: &[SiqsSuitablePrime],
    target_a: &Integer,
    skip: usize,
    factor_cap: Option<usize>,
) -> Option<SiqsPolynomial> {
    if n <= &Integer::from(1) {
        return None;
    }
    // Exact legacy cap: the per-index builder stopped the sieve at
    // skip + 65 suitable primes, so the greedy window never sees more.
    let capped = suitable.len().min(skip + 65);
    if capped <= skip + 2 {
        return None;
    }
    let window = &suitable[skip..capped];
    let mut coeff_a = Integer::from(1);
    let mut factors = Vec::new();
    for entry in window {
        if &coeff_a >= target_a {
            break;
        }
        if factor_cap.is_some_and(|cap| factors.len() >= cap) {
            break;
        }
        coeff_a *= Integer::from(entry.prime);
        factors.push(entry.prime);
    }
    if factors.len() < 2 || &coeff_a < target_a {
        return None;
    }
    let mut coeff_b = Integer::from(0);
    for prime in &factors {
        let modulus = Integer::from(*prime);
        let gamma = suitable_gamma(suitable, *prime)?;
        let rest = &coeff_a / &modulus;
        let inv = egcd_inverse(&rest, &modulus)?;
        let term = Integer::from(gamma) * rest * inv;
        coeff_b += term;
    }
    coeff_b %= &coeff_a;
    if &coeff_b % Integer::from(2) == Integer::from(0) {
        coeff_b = &coeff_a - &coeff_b;
    }
    let coeff_c = (&coeff_b * &coeff_b - n) / &coeff_a;
    if (&coeff_b * &coeff_b - n) % &coeff_a != Integer::from(0) {
        return None;
    }
    Some(SiqsPolynomial {
        coeff_a,
        coeff_b,
        coeff_c,
        prime_factors: factors,
    })
}

/// Cached `gamma` lookup for a factor prime known to come from `suitable`.
fn suitable_gamma(suitable: &[SiqsSuitablePrime], prime: u64) -> Option<u64> {
    suitable
        .iter()
        .find(|entry| entry.prime == prime)
        .map(|entry| entry.gamma)
}

/// Single-prime-flip incremental polynomial update (X1d Gray-code
/// mechanism): replace `out_prime` by the fresh prime `in_prime` with O(1)
/// BigInt ops. `sign` is the input polynomial's uniform CRT sign (`+1` if
/// `b` takes `+gamma` on the kept factors, `-1` for the conjugate
/// choice); installing `sign * in_gamma` on the new prime keeps the whole
/// solution in the `sign` class, so by CRT uniqueness the reduced
/// representative — and hence the parity-fix branch — matches a full O(k)
/// recomputation exactly (proven by `gray_flip_matches_crt_characterization`).
/// Test-only: Gray-code resequence wiring was reverted (see X1 report), so
/// the flip survives as the proven reference for the O(1) update.
#[cfg(test)]
fn gray_flip_polynomial(
    n: &Integer,
    poly: &SiqsPolynomial,
    out_prime: u64,
    in_prime: u64,
    in_gamma: u64,
    sign: i32,
) -> Option<SiqsPolynomial> {
    if out_prime == in_prime {
        return None;
    }
    if &poly.coeff_a % in_prime == Integer::from(0) {
        return None;
    }
    let half = &poly.coeff_a / Integer::from(out_prime);
    if &half * Integer::from(out_prime) != poly.coeff_a {
        return None;
    }
    let new_a = &half * Integer::from(in_prime);
    let modulus = Integer::from(in_prime);
    let half_mod = &half % &modulus;
    let inv = egcd_inverse(&half_mod, &modulus)?;
    let target = Integer::from(if sign >= 0 {
        in_gamma
    } else {
        (in_prime - in_gamma) % in_prime
    });
    let b_mod = &poly.coeff_b % &modulus;
    let delta = (target - b_mod + &modulus) % &modulus;
    let step = (delta * inv) % &modulus;
    let mut new_b = (&poly.coeff_b + &half * &step) % &new_a;
    if &new_b % Integer::from(2) == Integer::from(0) {
        new_b = &new_a - &new_b;
    }
    let new_c = (&new_b * &new_b - n) / &new_a;
    if (&new_b * &new_b - n) % &new_a != Integer::from(0) {
        return None;
    }
    let mut factors = poly.prime_factors.clone();
    let position = factors.iter().position(|&prime| prime == out_prime)?;
    factors[position] = in_prime;
    Some(SiqsPolynomial {
        coeff_a: new_a,
        coeff_b: new_b,
        coeff_c: new_c,
        prime_factors: factors,
    })
}

pub(crate) fn target_a_bits(n_bits: u64, interval_bits: u32) -> u64 {
    n_bits.div_ceil(2) + 1 - u64::from(interval_bits)
}

pub(crate) fn siqs_extract(
    n: &Integer,
    xs: &[Integer],
    smooth_products: &[Integer],
) -> Result<Integer, FactorError> {
    if xs.is_empty() || xs.len() != smooth_products.len() {
        return Err(FactorError::AlgebraicSquareRootFailed);
    }
    let mut big_x = Integer::from(1);
    let mut big_y = Integer::from(1);
    for (x, product) in xs.iter().zip(smooth_products.iter()) {
        big_x *= x;
        big_y *= product;
    }
    big_x %= n;
    let root = big_y.sqrt();
    if &root * &root != big_y {
        return Err(FactorError::AlgebraicSquareRootFailed);
    }
    if (&big_x * &big_x - &big_y) % n != Integer::from(0) {
        return Err(FactorError::AlgebraicSquareRootFailed);
    }
    for candidate in [
        integer_gcd(&(&big_x - &root), n),
        integer_gcd(&(&big_x + &root), n),
    ] {
        if candidate > Integer::from(1) && &candidate < n && n % &candidate == Integer::from(0) {
            return Ok(candidate);
        }
    }
    Err(FactorError::AlgebraicSquareRootFailed)
}

pub(crate) struct SiqsOutcome {
    pub(crate) display: PolynomialPair,
    pub(crate) relations: usize,
    pub(crate) rank: usize,
}

fn check_cancelled(cancelled: &AtomicBool, deadline: Option<Instant>) -> Result<(), FactorError> {
    if cancelled.load(Ordering::Relaxed) {
        return Err(FactorError::Cancelled);
    }
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return Err(FactorError::ResourceLimit(
            "siqs collection deadline exceeded".to_string(),
        ));
    }
    Ok(())
}

fn toggle_first_column(base: &FactorBase, mask: &mut [u64], prime: u64) {
    let Some(indices) = base.algebraic_by_prime.get(&prime) else {
        return;
    };
    let Some(first) = indices.first() else {
        return;
    };
    let Some(entry) = base.algebraic_primes.get(*first) else {
        return;
    };
    mask[entry.id / 64] ^= 1_u64 << (entry.id % 64);
}

/// Per-attempt dense toggle index (X1b): `prime -> column` arrays replacing
/// per-toggle `HashMap` lookups in the classify hot path.
/// `u32::MAX` marks "absent". Semantics match
/// [`FactorBase::toggle_rational`] + first-column fallback exactly:
/// a rational hit toggles only the rational column, otherwise the first
/// algebraic column for that prime (if any) is toggled.
#[derive(Clone, Debug)]
pub(crate) struct SiqsToggleIndex {
    rational_col: Vec<u32>,
    first_alg_col: Vec<u32>,
}

impl SiqsToggleIndex {
    pub(crate) const ABSENT: u32 = u32::MAX;

    pub(crate) fn build(base: &FactorBase) -> Self {
        let max = base
            .rational_primes
            .iter()
            .chain(base.algebraic_modulus_primes.iter())
            .copied()
            .max()
            .unwrap_or(0);
        let mut index = Self {
            rational_col: vec![Self::ABSENT; max as usize + 1],
            first_alg_col: vec![Self::ABSENT; max as usize + 1],
        };
        for (column, prime) in base.rational_primes.iter().enumerate() {
            index.rational_col[*prime as usize] = column as u32;
        }
        for (prime, slots) in &base.algebraic_by_prime {
            let Some(first) = slots.first() else {
                continue;
            };
            let Some(entry) = base.algebraic_primes.get(*first) else {
                continue;
            };
            let slot = &mut index.first_alg_col[*prime as usize];
            if *slot == Self::ABSENT {
                *slot = entry.id as u32;
            }
        }
        index
    }

    fn toggle(&self, mask: &mut [u64], prime: u64) {
        let Some(&column) = self.rational_col.get(prime as usize) else {
            return;
        };
        if column != Self::ABSENT {
            mask[column as usize / 64] ^= 1_u64 << (column as usize % 64);
            return;
        }
        if let Some(&column) = self.first_alg_col.get(prime as usize)
            && column != Self::ABSENT
        {
            mask[column as usize / 64] ^= 1_u64 << (column as usize % 64);
        }
    }
}

/// Barrett reduction for the u64 trial-division path (X1c): one precomputed
/// `mu = floor(2^64 / prime)` per prime replaces each trial `rest % prime`
/// DIV with a MUL plus at most one correction. The BigInt path is untouched.
/// Correctness: `q_hat = (value * mu) >> 64` satisfies `q - 1 <= q_hat <= q`
/// for the true quotient `q` (since `mu <= 2^64 / prime` gives `q_hat <= q`,
/// and `mu > 2^64 / prime - 1` with `value < 2^64` gives `q - q_hat < 2`),
/// so a single `r >= prime` correction yields the exact remainder.
fn barrett_mu(prime: u64) -> u64 {
    debug_assert!(prime >= 2);
    ((1_u128 << 64) / prime as u128) as u64
}

fn barrett_divmod(value: u64, prime: u64, mu: u64) -> (u64, u64) {
    debug_assert!(prime >= 2);
    let mut quotient = ((value as u128 * mu as u128) >> 64) as u64;
    let mut remainder = value - quotient * prime;
    if remainder >= prime {
        remainder -= prime;
        quotient += 1;
    }
    debug_assert!(remainder < prime);
    debug_assert_eq!(quotient, value / prime);
    debug_assert_eq!(remainder, value % prime);
    (quotient, remainder)
}

/// Per-attempt Barrett table (X1c): `mu` for every prime up to the largest
/// trial prime (base primes and polynomial factors alike). Out-of-range
/// primes fall back to native division.
#[derive(Clone, Debug)]
pub(crate) struct SiqsBarrettTable {
    mu: Vec<u64>,
}

impl SiqsBarrettTable {
    pub(crate) fn build(cover_through: u64) -> Self {
        let mut mu = vec![0_u64; cover_through as usize + 1];
        for prime in 2..=cover_through {
            mu[prime as usize] = barrett_mu(prime);
        }
        Self { mu }
    }

    pub(crate) fn divmod(&self, value: u64, prime: u64) -> (u64, u64) {
        match self.mu.get(prime as usize) {
            Some(&mu) if prime >= 2 => barrett_divmod(value, prime, mu),
            _ => (value / prime, value % prime),
        }
    }
}

fn toggle_parity(base: &FactorBase, mask: &mut [u64], prime: u64) {
    if !base.toggle_rational(mask, prime).unwrap_or(false) {
        toggle_first_column(base, mask, prime);
    }
}

fn toggle_parity_fast(base: &FactorBase, toggles: &SiqsToggleIndex, mask: &mut [u64], prime: u64) {
    if (prime as usize) < toggles.rational_col.len() {
        toggles.toggle(mask, prime);
    } else {
        toggle_parity(base, mask, prime);
    }
}

fn strip_prime_u64(
    rest: &mut u64,
    prime: u64,
    base: &FactorBase,
    toggles: &SiqsToggleIndex,
    barrett: &SiqsBarrettTable,
    mask: &mut [u64],
) {
    if *rest < prime {
        return;
    }
    let mut exponent = 0_u32;
    loop {
        let (quotient, remainder) = barrett.divmod(*rest, prime);
        if remainder != 0 {
            break;
        }
        *rest = quotient;
        exponent += 1;
    }
    if exponent & 1 != 0 {
        toggle_parity_fast(base, toggles, mask, prime);
    }
}

#[derive(Debug)]
enum Classify {
    Smooth,
    Partial(u64),
    Partial2(u64, u64),
    Reject,
}

const SPLITTER_PRIME_CAP: usize = 1_000_000;

pub(crate) fn splitter_primes(base: &FactorBase, lp_bound: u64) -> Vec<u64> {
    let floor = base
        .rational_primes
        .last()
        .copied()
        .max(base.algebraic_modulus_primes.last().copied())
        .unwrap_or(1);
    let cap = (lp_bound as usize).min(SPLITTER_PRIME_CAP);
    if cap <= floor as usize {
        return Vec::new();
    }
    sieve_primes(cap)
        .into_iter()
        .filter(|prime| *prime > floor)
        .collect()
}

/// Per-polynomial SIQS measurement counters (X0: reporting only, no behavior change).
#[derive(Clone, Debug, Default)]
pub(crate) struct SiqsPolyCounters {
    pub(crate) sieve_mark_nanos: u64,
    pub(crate) classify_u64_count: u64,
    pub(crate) classify_u64_nanos: u64,
    pub(crate) classify_bigint_count: u64,
    pub(crate) classify_bigint_nanos: u64,
    pub(crate) splitter_count: u64,
    pub(crate) splitter_nanos: u64,
    pub(crate) merge_count: u64,
    pub(crate) evaluate_calls: u64,
    pub(crate) evaluate_nanos: u64,
    pub(crate) roots_nanos: u64,
    pub(crate) generate_nanos: u64,
    pub(crate) score_max: u64,
    pub(crate) pass_evaluated: [u64; 3],
    pub(crate) pass_relations: [u64; 3],
}

impl SiqsPolyCounters {
    pub(crate) fn accumulate(&mut self, other: &Self) {
        self.sieve_mark_nanos += other.sieve_mark_nanos;
        self.classify_u64_count += other.classify_u64_count;
        self.classify_u64_nanos += other.classify_u64_nanos;
        self.classify_bigint_count += other.classify_bigint_count;
        self.classify_bigint_nanos += other.classify_bigint_nanos;
        self.splitter_count += other.splitter_count;
        self.splitter_nanos += other.splitter_nanos;
        self.merge_count += other.merge_count;
        self.evaluate_calls += other.evaluate_calls;
        self.evaluate_nanos += other.evaluate_nanos;
        self.roots_nanos += other.roots_nanos;
        self.generate_nanos += other.generate_nanos;
        self.score_max = self.score_max.max(other.score_max);
        for index in 0..3 {
            self.pass_evaluated[index] += other.pass_evaluated[index];
            self.pass_relations[index] += other.pass_relations[index];
        }
    }

    pub(crate) fn millis(nanos: u64) -> f64 {
        nanos as f64 / 1_000_000.0
    }
}

pub(crate) fn siqs_rescue_floor(scores: &[u32]) -> u64 {
    if scores.is_empty() {
        return 0;
    }
    let mut ranked = scores.to_vec();
    ranked.sort_unstable();
    let rank = ranked.len().saturating_mul(20) / 100;
    ranked[rank.min(ranked.len() - 1)] as u64
}

pub(crate) fn siqs_trial_redundant(prime: u64, always_trial: &[u64]) -> bool {
    always_trial.contains(&prime)
}

fn split_residual_counted(
    rest: u64,
    splitters: &[u64],
    lp_bound: u64,
    stats: &mut SiqsPolyCounters,
) -> Classify {
    let start = Instant::now();
    let outcome = split_residual_inner(rest, splitters, lp_bound);
    stats.splitter_count += 1;
    stats.splitter_nanos += start.elapsed().as_nanos() as u64;
    outcome
}

fn split_residual_inner(mut rest: u64, splitters: &[u64], lp_bound: u64) -> Classify {
    if rest == 1 {
        return Classify::Smooth;
    }
    if rest <= lp_bound && is_prime_u64_deterministic(rest) {
        return Classify::Partial(rest);
    }
    if (rest as u128) > (lp_bound as u128).saturating_mul(lp_bound as u128) {
        return Classify::Reject;
    }
    let mut lps: Vec<u64> = Vec::new();
    for &prime in splitters {
        if rest == 1 {
            break;
        }
        if (prime as u128) * (prime as u128) > rest as u128 {
            break;
        }
        while rest.is_multiple_of(prime) {
            lps.push(prime);
            rest /= prime;
        }
    }
    if rest > 1 {
        if rest > lp_bound || !is_prime_u64_deterministic(rest) {
            return Classify::Reject;
        }
        lps.push(rest);
    }
    lps.sort_unstable();
    let mut i = 0;
    while i + 1 < lps.len() {
        if lps[i] == lps[i + 1] {
            lps.remove(i);
            lps.remove(i);
        } else {
            i += 1;
        }
    }
    if lps.iter().any(|&prime| prime > lp_bound) {
        return Classify::Reject;
    }
    match lps.len() {
        0 => Classify::Smooth,
        1 => Classify::Partial(lps[0]),
        2 => Classify::Partial2(lps[0], lps[1]),
        _ => Classify::Reject,
    }
}

#[cfg(test)]
fn classify_siqs(
    signed: &Integer,
    base: &FactorBase,
    mask: &mut [u64],
    lp_bound: u64,
    trial_primes: &[u64],
    always_trial: &[u64],
    splitters: &[u64],
) -> Classify {
    let mut discard = SiqsPolyCounters::default();
    let toggles = SiqsToggleIndex::build(base);
    let cover = base
        .rational_primes
        .iter()
        .chain(base.algebraic_modulus_primes.iter())
        .chain(trial_primes.iter())
        .chain(always_trial.iter())
        .copied()
        .max()
        .unwrap_or(2)
        .max(2);
    let barrett = SiqsBarrettTable::build(cover);
    classify_siqs_counted(
        signed,
        base,
        &toggles,
        &barrett,
        mask,
        lp_bound,
        trial_primes,
        always_trial,
        splitters,
        &mut discard,
    )
}

#[allow(clippy::too_many_arguments)]
fn classify_siqs_counted(
    signed: &Integer,
    base: &FactorBase,
    toggles: &SiqsToggleIndex,
    barrett: &SiqsBarrettTable,
    mask: &mut [u64],
    lp_bound: u64,
    trial_primes: &[u64],
    always_trial: &[u64],
    splitters: &[u64],
    stats: &mut SiqsPolyCounters,
) -> Classify {
    if signed == &Integer::from(0) {
        return Classify::Reject;
    }
    let negative = signed < &Integer::from(0);
    if negative {
        mask[base.sign_column / 64] ^= 1_u64 << (base.sign_column % 64);
    }
    let magnitude = if negative {
        -signed.clone()
    } else {
        signed.clone()
    };
    if let Some(mut rest) = magnitude.to_u64() {
        let start = Instant::now();
        for prime in std::iter::once(&2_u64).chain(always_trial.iter()) {
            if rest == 1 {
                break;
            }
            strip_prime_u64(&mut rest, *prime, base, toggles, barrett, mask);
        }
        for prime in trial_primes.iter() {
            if rest == 1 {
                break;
            }
            if siqs_trial_redundant(*prime, always_trial) {
                continue;
            }
            strip_prime_u64(&mut rest, *prime, base, toggles, barrett, mask);
        }
        let outcome = if rest == 1 {
            Classify::Smooth
        } else {
            split_residual_counted(rest, splitters, lp_bound, stats)
        };
        stats.classify_u64_count += 1;
        stats.classify_u64_nanos += start.elapsed().as_nanos() as u64;
        return outcome;
    }
    let start = Instant::now();
    let mut residual = magnitude;
    let prime_at = |index: usize| -> u64 {
        if index == 0 {
            2
        } else if index <= always_trial.len() {
            always_trial[index - 1]
        } else {
            trial_primes[index - 1 - always_trial.len()]
        }
    };
    let total = 1 + always_trial.len() + trial_primes.len();
    let mut demoted: Option<u64> = None;
    for index in 0..total {
        if residual == Integer::from(1) {
            break;
        }
        let prime = prime_at(index);
        if index > always_trial.len() && siqs_trial_redundant(prime, always_trial) {
            continue;
        }
        let prime_value = Integer::from(prime);
        if residual < prime_value {
            continue;
        }
        let mut exponent = 0_u32;
        while (&residual % &prime_value) == Integer::from(0) {
            residual /= &prime_value;
            exponent += 1;
        }
        if exponent & 1 != 0 {
            toggle_parity_fast(base, toggles, mask, prime);
        }
        if let Some(small) = residual.to_u64() {
            let mut rest = small;
            for remaining in (index + 1)..total {
                if rest == 1 {
                    break;
                }
                let trial_prime = prime_at(remaining);
                if remaining > always_trial.len() && siqs_trial_redundant(trial_prime, always_trial)
                {
                    continue;
                }
                strip_prime_u64(&mut rest, trial_prime, base, toggles, barrett, mask);
            }
            demoted = Some(rest);
            break;
        }
    }
    let outcome = match demoted {
        Some(rest) => {
            if rest == 1 {
                Classify::Smooth
            } else {
                split_residual_counted(rest, splitters, lp_bound, stats)
            }
        }
        None => {
            if residual == Integer::from(1) {
                Classify::Smooth
            } else {
                match residual.to_u64() {
                    Some(small) => split_residual_counted(small, splitters, lp_bound, stats),
                    None => Classify::Reject,
                }
            }
        }
    };
    stats.classify_bigint_count += 1;
    stats.classify_bigint_nanos += start.elapsed().as_nanos() as u64;
    outcome
}

fn sieve_roots_for(
    base: &FactorBase,
    poly: &SiqsPolynomial,
    gammas: &SiqsGammaTable,
) -> Vec<(u64, Vec<u64>, u32)> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for prime in base
        .rational_primes
        .iter()
        .chain(base.algebraic_modulus_primes.iter())
    {
        if *prime < 3 || !seen.insert(*prime) {
            continue;
        }
        let roots = poly.roots_mod_with_gamma(*prime, gammas.get(*prime));
        if roots.is_empty() {
            continue;
        }
        let weight = (u64::BITS - prime.leading_zeros()).max(1);
        out.push((*prime, roots, weight));
    }
    out
}

type PolyOutput = (Vec<Relation>, Vec<(Relation, Vec<u64>)>);

struct LpEdge {
    rel: Relation,
}

struct LpGraph {
    adj: std::collections::BTreeMap<u64, Vec<(u64, usize)>>,
    edges: Vec<LpEdge>,
}

impl LpGraph {
    fn new() -> Self {
        Self {
            adj: std::collections::BTreeMap::new(),
            edges: Vec::new(),
        }
    }

    fn ingest(&mut self, rel: Relation, mut lps: Vec<u64>, n: &Integer, out: &mut Vec<Relation>) {
        lps.sort_unstable();
        let ends = match lps.len() {
            1 => [0, lps[0]],
            2 => [lps[0], lps[1]],
            _ => return,
        };
        let mut parent: std::collections::BTreeMap<u64, (u64, usize)> =
            std::collections::BTreeMap::new();
        let mut queue = std::collections::VecDeque::from([ends[0]]);
        parent.insert(ends[0], (ends[0], usize::MAX));
        let mut found = false;
        while let Some(node) = queue.pop_front() {
            if node == ends[1] {
                found = true;
                break;
            }
            if let Some(neighbors) = self.adj.get(&node) {
                for &(next, edge_idx) in neighbors {
                    if let std::collections::btree_map::Entry::Vacant(slot) = parent.entry(next) {
                        slot.insert((node, edge_idx));
                        queue.push_back(next);
                    }
                }
            }
        }
        if !found {
            let edge_idx = self.edges.len();
            self.edges.push(LpEdge { rel });
            self.adj
                .entry(ends[0])
                .or_default()
                .push((ends[1], edge_idx));
            self.adj
                .entry(ends[1])
                .or_default()
                .push((ends[0], edge_idx));
            return;
        }
        let mut combined_a = rel.a.clone();
        let mut combined_value = rel.rational_value.clone();
        let mut combined_mask = rel.mask.clone();
        let mut node = ends[1];
        while node != ends[0] {
            let (prev, edge_idx) = parent[&node];
            let edge = &self.edges[edge_idx];
            combined_a = (&combined_a * &edge.rel.a) % n;
            combined_value = &combined_value * &edge.rel.rational_value;
            for (target, source) in combined_mask.iter_mut().zip(edge.rel.mask.iter()) {
                *target ^= *source;
            }
            node = prev;
        }
        if combined_mask.iter().all(|word| *word == 0) {
            return;
        }
        out.push(Relation {
            a: combined_a.clone(),
            b: rel.b.clone(),
            rational_value: combined_value,
            algebraic_value: QuadraticElement {
                x: combined_a,
                y: Integer::from(0),
            },
            mask: combined_mask,
        });
    }
}

fn siqs_worker_count(options: &FactorOptions) -> usize {
    let available = std::thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    options.threads.unwrap_or(available).max(1)
}

fn push_unique(
    seen: &mut std::collections::HashSet<Vec<u64>>,
    out: &mut Vec<Relation>,
    relation: Relation,
) {
    if seen.insert(relation.mask.clone()) {
        out.push(relation);
    }
}

pub(crate) const SIQS_FILTER_MERGE_CAPS: [usize; 3] = [64, 128, 256];
pub(crate) const SIQS_FILTER_OVERWEIGHT_CAP: usize = 256;
#[cfg(test)]
pub(crate) const SIQS_FILTER_MEMORY_BYTES: u64 = 256 * 1024 * 1024;
pub(crate) const SIQS_WIDTH_ESCALATION_STEPS: [usize; 2] = [4, 16];

pub(crate) const SIQS_TABLE_SMALL_DIGITS: u64 = 36;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SiqsTableEntry {
    pub(crate) rational_bound: usize,
    pub(crate) algebraic_bound: usize,
    pub(crate) half_interval: usize,
    pub(crate) lp_bound: usize,
    pub(crate) a_factors: Option<usize>,
}

pub(crate) fn siqs_decimal_digits(n: &Integer) -> u64 {
    n.bits() * 301 / 1000 + 1
}

pub(crate) fn siqs_table_bounds(digits: u64) -> SiqsTableEntry {
    if digits <= SIQS_TABLE_SMALL_DIGITS {
        let defaults = FactorOptions::default();
        return SiqsTableEntry {
            rational_bound: defaults.rational_base_bound,
            algebraic_bound: defaults.algebraic_base_bound,
            half_interval: defaults.sieve_a_bound,
            lp_bound: defaults.large_prime_bound,
            a_factors: None,
        };
    }
    if digits <= 60 {
        return SiqsTableEntry {
            rational_bound: 16_000,
            algebraic_bound: 64_000,
            half_interval: 65_536,
            lp_bound: 2_000_000,
            a_factors: None,
        };
    }
    if digits <= 74 {
        return SiqsTableEntry {
            rational_bound: 64_000,
            algebraic_bound: 250_000,
            half_interval: 65_536,
            lp_bound: 8_000_000,
            a_factors: Some(9),
        };
    }
    if digits <= 80 {
        return SiqsTableEntry {
            rational_bound: 100_000,
            algebraic_bound: 300_000,
            half_interval: 196_608,
            lp_bound: 67_108_864,
            a_factors: Some(10),
        };
    }
    SiqsTableEntry {
        rational_bound: 200_000,
        algebraic_bound: 600_000,
        half_interval: 393_216,
        lp_bound: 134_217_728,
        a_factors: Some(10),
    }
}

pub(crate) fn siqs_memory_budget_bytes(options: &FactorOptions) -> u64 {
    options.memory_budget_bytes()
}

pub(crate) fn siqs_collect_estimate_bytes(columns: usize, half_width: i64, workers: usize) -> u64 {
    estimate_collect_bytes(columns, half_width, workers)
}

pub(crate) fn siqs_rank_estimate_bytes(columns: usize, relations: usize) -> u64 {
    estimate_total_bytes(columns, relations)
}

fn resolve_siqs_table(n: &Integer, options: &mut FactorOptions) -> Option<usize> {
    let defaults = FactorOptions::default();
    let entry = siqs_table_bounds(siqs_decimal_digits(n));
    if options.rational_base_bound == defaults.rational_base_bound {
        options.rational_base_bound = entry.rational_bound;
    }
    if options.algebraic_base_bound == defaults.algebraic_base_bound {
        options.algebraic_base_bound = entry.algebraic_bound;
    }
    if options.sieve_a_bound == defaults.sieve_a_bound {
        options.sieve_a_bound = entry.half_interval;
    }
    if options.large_prime_bound == defaults.large_prime_bound {
        options.large_prime_bound = entry.lp_bound;
    }
    let table_active = options.rational_base_bound == entry.rational_bound
        || options.algebraic_base_bound == entry.algebraic_bound
        || options.sieve_a_bound == entry.half_interval
        || options.large_prime_bound == entry.lp_bound;
    if table_active { entry.a_factors } else { None }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SiqsFilterStats {
    pub(crate) zero_masks_dropped: usize,
    pub(crate) singletons_removed: usize,
    pub(crate) merges: usize,
    pub(crate) overweight_dropped: usize,
}

fn mask_weight(mask: &[u64]) -> usize {
    mask.iter().map(|word| word.count_ones() as usize).sum()
}

fn mask_columns(mask: &[u64], column_count: usize) -> Vec<usize> {
    let mut columns = Vec::new();
    for (word_index, word) in mask.iter().enumerate() {
        let mut rest = *word;
        while rest != 0 {
            let bit = rest.trailing_zeros() as usize;
            let column = word_index * 64 + bit;
            if column < column_count {
                columns.push(column);
            }
            rest &= rest - 1;
        }
    }
    columns
}

fn merge_relation_pair(n: &Integer, first: &Relation, second: &Relation) -> Relation {
    let mut mask = first.mask.clone();
    for (target, source) in mask.iter_mut().zip(second.mask.iter()) {
        *target ^= *source;
    }
    let combined_a = (&first.a * &second.a) % n;
    let combined_value = &first.rational_value * &second.rational_value;
    Relation {
        a: combined_a.clone(),
        b: first.b.clone(),
        rational_value: combined_value,
        algebraic_value: QuadraticElement {
            x: combined_a,
            y: Integer::from(0),
        },
        mask,
    }
}

#[cfg(test)]
pub(crate) fn filter_relations(
    n: &Integer,
    relations: Vec<Relation>,
    column_count: usize,
    merge_caps: &[usize],
    overweight_cap: usize,
) -> (Vec<Relation>, SiqsFilterStats) {
    filter_relations_with_budget(
        n,
        relations,
        column_count,
        merge_caps,
        overweight_cap,
        SIQS_FILTER_MEMORY_BYTES,
    )
}

pub(crate) fn filter_relations_with_budget(
    n: &Integer,
    relations: Vec<Relation>,
    column_count: usize,
    merge_caps: &[usize],
    overweight_cap: usize,
    budget_bytes: u64,
) -> (Vec<Relation>, SiqsFilterStats) {
    let mut stats = SiqsFilterStats::default();
    let mask_words = column_count.div_ceil(64);
    let estimate = relations.len() as u64 * mask_words as u64 * 8;
    if estimate > budget_bytes {
        return (relations, stats);
    }
    let mut kept: Vec<Relation> = Vec::with_capacity(relations.len());
    for relation in relations {
        if relation.mask.iter().all(|word| *word == 0) {
            stats.zero_masks_dropped += 1;
            continue;
        }
        kept.push(relation);
    }
    let mut alive = vec![true; kept.len()];
    let mut alive_count = kept.len();
    loop {
        let mut work = FilterWork::build(&kept, &alive, column_count);
        let singles = work.remove_singletons(&kept, &mut alive);
        stats.singletons_removed += singles;
        alive_count -= singles;
        let mut round_merges = 0_usize;
        for cap in merge_caps {
            loop {
                let merges = work.merge_pairs(n, &mut kept, &mut alive, *cap);
                if merges == 0 {
                    break;
                }
                round_merges += merges;
                alive_count -= merges;
                if alive_count == 0 {
                    break;
                }
            }
            if alive_count == 0 {
                break;
            }
        }
        stats.merges += round_merges;
        if singles == 0 && round_merges == 0 {
            break;
        }
        if alive_count == 0 {
            break;
        }
    }
    if alive_count >= 3 * column_count.max(1) {
        let mut dropped = 0_usize;
        for (index, relation) in kept.iter().enumerate() {
            if alive[index] && mask_weight(&relation.mask) > overweight_cap {
                alive[index] = false;
                dropped += 1;
            }
        }
        alive_count -= dropped;
        stats.overweight_dropped = dropped;
    }
    let mut out = Vec::with_capacity(alive_count);
    for (index, relation) in kept.into_iter().enumerate() {
        if alive[index] {
            out.push(relation);
        }
    }
    (out, stats)
}

struct FilterWork {
    members: Vec<Vec<usize>>,
    weight: Vec<u32>,
}

impl FilterWork {
    fn build(relations: &[Relation], alive: &[bool], column_count: usize) -> Self {
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); column_count];
        let mut weight = vec![0_u32; column_count];
        for (index, relation) in relations.iter().enumerate() {
            if !alive[index] {
                continue;
            }
            for column in mask_columns(&relation.mask, column_count) {
                members[column].push(index);
                weight[column] += 1;
            }
        }
        Self { members, weight }
    }

    fn alive_members(
        &self,
        relations: &[Relation],
        alive: &[bool],
        column: usize,
        column_count: usize,
    ) -> Vec<usize> {
        let mut out = Vec::new();
        for &index in &self.members[column] {
            if alive[index] && mask_columns(&relations[index].mask, column_count).contains(&column)
            {
                out.push(index);
            }
        }
        out
    }

    fn remove_singletons(&mut self, relations: &[Relation], alive: &mut [bool]) -> usize {
        let column_count = self.weight.len();
        let mut queue = Vec::new();
        for (column, weight) in self.weight.iter().enumerate() {
            if *weight == 1 {
                queue.push(column);
            }
        }
        let mut removed = 0_usize;
        while let Some(column) = queue.pop() {
            if self.weight[column] != 1 {
                continue;
            }
            let owners = self.alive_members(relations, alive, column, column_count);
            if owners.len() != 1 {
                continue;
            }
            let index = owners[0];
            alive[index] = false;
            removed += 1;
            for owned in mask_columns(&relations[index].mask, column_count) {
                self.weight[owned] -= 1;
                if self.weight[owned] == 1 {
                    queue.push(owned);
                }
            }
        }
        removed
    }

    fn merge_pairs(
        &mut self,
        n: &Integer,
        relations: &mut [Relation],
        alive: &mut [bool],
        cap: usize,
    ) -> usize {
        let column_count = self.weight.len();
        let mut consumed = vec![false; relations.len()];
        let mut merges = 0_usize;
        for column in 0..column_count {
            if self.weight[column] != 2 {
                continue;
            }
            let owners = self.alive_members(relations, alive, column, column_count);
            if owners.len() != 2 {
                continue;
            }
            let (first_index, second_index) = (owners[0], owners[1]);
            if consumed[first_index] || consumed[second_index] {
                continue;
            }
            let mut weight = 0_usize;
            for (left, right) in relations[first_index]
                .mask
                .iter()
                .zip(relations[second_index].mask.iter())
            {
                weight += (left ^ right).count_ones() as usize;
            }
            if weight == 0 || weight > cap {
                continue;
            }
            let merged = merge_relation_pair(n, &relations[first_index], &relations[second_index]);
            let old_first = mask_columns(&relations[first_index].mask, column_count);
            let old_second = mask_columns(&relations[second_index].mask, column_count);
            let new_mask = mask_columns(&merged.mask, column_count);
            relations[first_index] = merged;
            alive[second_index] = false;
            for owned in &old_second {
                self.weight[*owned] -= 1;
            }
            for owned in &old_first {
                if !new_mask.contains(owned) {
                    self.weight[*owned] -= 1;
                }
            }
            for owned in &new_mask {
                if !old_first.contains(owned) {
                    self.weight[*owned] += 1;
                    self.members[*owned].push(first_index);
                }
            }
            consumed[first_index] = true;
            consumed[second_index] = true;
            merges += 1;
        }
        merges
    }
}

#[allow(clippy::too_many_arguments)]
fn collect_one_poly<F>(
    n: &Integer,
    poly: &SiqsPolynomial,
    nonce: u64,
    base: &FactorBase,
    gammas: &SiqsGammaTable,
    toggles: &SiqsToggleIndex,
    barrett: &SiqsBarrettTable,
    half_width: i64,
    lp_bound: u64,
    splitters: &[u64],
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
    target: u64,
) -> Result<(PolyOutput, SiqsPolyCounters), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let mut stats = SiqsPolyCounters::default();
    let width = (half_width as usize) * 2 + 1;
    let words = base.column_count().div_ceil(64);
    let roots_start = Instant::now();
    let roots = sieve_roots_for(base, poly, gammas);
    stats.roots_nanos += roots_start.elapsed().as_nanos() as u64;
    // X1e verdict (measured, see X1 report): scores stay tiny in practice
    // (score_max <= 78 across the 82/100-bit canonical corpus vs u16::MAX)
    // and evaluate is ~1-4% of collect CPU, so u32 scores and the plain
    // BigInt evaluate stay: a u16/Horner rewrite would add risk for no
    // measurable win. Revisit only if score_max approaches 65535 or the
    // evaluate share exceeds 5%.
    let mut scores = vec![0_u32; width];
    let mut hits: Vec<Vec<u64>> = vec![Vec::new(); width];
    let mark_start = Instant::now();
    let mut poly_score_max = 0_u64;
    for (prime, hit_roots, weight) in &roots {
        for root in hit_roots {
            let start = ((*root as i64 + half_width).rem_euclid(*prime as i64)) as usize;
            let mut index = start;
            while index < width {
                scores[index] = scores[index].saturating_add(*weight);
                poly_score_max = poly_score_max.max(scores[index] as u64);
                let slot = &mut hits[index];
                if slot.last().copied() != Some(*prime) {
                    slot.push(*prime);
                }
                index += *prime as usize;
            }
        }
    }
    stats.score_max = stats.score_max.max(poly_score_max);
    stats.sieve_mark_nanos += mark_start.elapsed().as_nanos() as u64;
    let ends = [
        poly.evaluate(&Integer::from(-half_width)),
        poly.evaluate(&Integer::from(half_width)),
    ];
    stats.evaluate_calls += 2;
    let mut max_bits = 1_u64;
    for end in &ends {
        let mag = if end < &Integer::from(0) {
            -end.clone()
        } else {
            end.clone()
        };
        max_bits = max_bits.max(mag.bits());
    }
    // X1f verdict (measured, see X1 report): the three-pass ladder stays.
    // Collapsing to a single adaptive threshold was allowed only if passes
    // 1/2 never add relations on the canonical corpus; instead they add
    // 48+78 relations on the decisive 100-bit attempt (and contribute on
    // every other SIQS attempt measured), so the rescue passes earn their
    // keep. The pass histogram in the siqs_collect detail line monitors this.
    let ladder = [max_bits.saturating_sub(8), max_bits.saturating_sub(16), 0];
    let mut rescue_floor = 0_u64;
    let mut rescue_ready = false;
    let mut smooth: Vec<Relation> = Vec::new();
    let mut partials: Vec<(Relation, Vec<u64>)> = Vec::new();
    // X1b: one mask buffer per polynomial, cleared per candidate; only the
    // rare accepted relations pay for a clone.
    let mut scratch_mask = vec![0_u64; words];
    for (pass, base_threshold) in ladder.iter().enumerate() {
        let threshold = if pass == 2 {
            if !rescue_ready {
                rescue_floor = siqs_rescue_floor(&scores);
                rescue_ready = true;
            }
            rescue_floor
        } else {
            *base_threshold
        };
        let pass_slot = pass.min(2);
        for (offset, score) in scores.iter().enumerate() {
            if offset & 1023 == 0 {
                check_cancelled(cancelled, deadline)?;
            }
            if (*score as u64) < threshold {
                continue;
            }
            let x = Integer::from(offset as i64 - half_width);
            let eval_start = Instant::now();
            let value = poly.evaluate(&x);
            stats.evaluate_nanos += eval_start.elapsed().as_nanos() as u64;
            stats.evaluate_calls += 1;
            if value == Integer::from(0) {
                continue;
            }
            stats.pass_evaluated[pass_slot] += 1;
            let signed = &poly.coeff_a * &value;
            scratch_mask.fill(0);
            match classify_siqs_counted(
                &signed,
                base,
                toggles,
                barrett,
                &mut scratch_mask,
                lp_bound,
                &hits[offset],
                &poly.prime_factors,
                splitters,
                &mut stats,
            ) {
                Classify::Reject => continue,
                Classify::Smooth => {}
                Classify::Partial(lp) => {
                    if scratch_mask.iter().all(|word| *word == 0) {
                        continue;
                    }
                    let point = poly.x_value(&x) % n;
                    partials.push((
                        Relation {
                            a: point.clone() % n,
                            b: Integer::from(nonce),
                            rational_value: signed.clone(),
                            algebraic_value: QuadraticElement {
                                x: point,
                                y: Integer::from(0),
                            },
                            mask: scratch_mask.clone(),
                        },
                        vec![lp],
                    ));
                    stats.pass_relations[pass_slot] += 1;
                    continue;
                }
                Classify::Partial2(first, second) => {
                    if scratch_mask.iter().all(|word| *word == 0) {
                        continue;
                    }
                    let point = poly.x_value(&x) % n;
                    partials.push((
                        Relation {
                            a: point.clone() % n,
                            b: Integer::from(nonce),
                            rational_value: signed.clone(),
                            algebraic_value: QuadraticElement {
                                x: point,
                                y: Integer::from(0),
                            },
                            mask: scratch_mask.clone(),
                        },
                        vec![first, second],
                    ));
                    stats.pass_relations[pass_slot] += 1;
                    continue;
                }
            }
            if scratch_mask.iter().all(|word| *word == 0) {
                continue;
            }
            let point = poly.x_value(&x) % n;
            smooth.push(Relation {
                a: point.clone() % n,
                b: Integer::from(nonce),
                rational_value: signed,
                algebraic_value: QuadraticElement {
                    x: point,
                    y: Integer::from(0),
                },
                mask: scratch_mask.clone(),
            });
            stats.pass_relations[pass_slot] += 1;
        }
        if !smooth.is_empty() || pass == ladder.len() - 1 {
            break;
        }
        progress(&ProgressEvent {
            stage: "collect".to_string(),
            current: (smooth.len() + partials.len()) as u64,
            total: target,
            message: "sieving siqs polynomial".to_string(),
        })?;
    }
    Ok(((smooth, partials), stats))
}

fn siqs_extract_relations(
    n: &Integer,
    relations: &[Relation],
    dependency: &[usize],
) -> Result<Integer, FactorError> {
    let mut xs = Vec::with_capacity(dependency.len());
    let mut products = Vec::with_capacity(dependency.len());
    for index in dependency {
        let relation = relations
            .get(*index)
            .ok_or_else(|| FactorError::LinearAlgebra("invalid relation index".to_string()))?;
        xs.push(relation.algebraic_value.x.clone());
        products.push(relation.rational_value.clone());
    }
    siqs_extract(n, &xs, &products)
}

pub(crate) const SIQS_POLY_INDEX_CAP: usize = 4096;

fn window_reaches_target(
    suitable: &[SiqsSuitablePrime],
    skip: usize,
    target_a: &Integer,
    factor_cap: Option<usize>,
) -> bool {
    let capped = suitable.len().min(skip + 65);
    if capped <= skip + 2 {
        return false;
    }
    let mut product = Integer::from(1);
    let mut factors = 0_usize;
    for entry in &suitable[skip..capped] {
        if &product >= target_a {
            break;
        }
        if factor_cap.is_some_and(|cap| factors >= cap) {
            break;
        }
        product *= Integer::from(entry.prime);
        factors += 1;
    }
    factors >= 2 && &product >= target_a
}

pub(crate) fn first_viable_poly_index(
    suitable: &[SiqsSuitablePrime],
    target_a: &Integer,
    factor_cap: Option<usize>,
) -> usize {
    for index in 0..=SIQS_POLY_INDEX_CAP {
        if window_reaches_target(suitable, index * 3, target_a, factor_cap) {
            return index;
        }
    }
    usize::MAX
}

#[derive(Clone, Debug, Default)]
pub(crate) struct SiqsBank {
    pub relations: Vec<Relation>,
    pub next_poly: usize,
}

impl SiqsBank {
    pub(crate) fn new() -> Self {
        Self {
            relations: Vec::new(),
            next_poly: 0,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.relations.len()
    }

    pub(crate) fn ingest_new(&mut self, fresh: Vec<Relation>, next_poly: usize) {
        let mut seen: std::collections::HashSet<Vec<u64>> = self
            .relations
            .iter()
            .map(|item| item.mask.clone())
            .collect();
        for relation in fresh {
            if seen.insert(relation.mask.clone()) {
                self.relations.push(relation);
            }
        }
        self.next_poly = next_poly.max(self.next_poly);
    }
}

pub(crate) fn siqs_pilot_with_bank<F>(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
    timings: Option<&TimingCollector>,
    mut bank: SiqsBank,
) -> (Option<(Integer, SiqsOutcome)>, SiqsBank)
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let mut working = options.clone();
    let factor_cap = resolve_siqs_table(n, &mut working);
    let memory_budget = siqs_memory_budget_bytes(options);
    let base = FactorBase::for_siqs(n, working.rational_base_bound, working.algebraic_base_bound);
    let columns = base.column_count();
    let half_width = match i64::try_from(working.sieve_a_bound) {
        Ok(value) if value > 0 => value,
        _ => return (None, bank),
    };
    let workers = siqs_worker_count(&working);
    let pilot_collect_allow =
        siqs_collect_estimate_bytes(columns, half_width, workers) <= memory_budget;
    crate::factor::record_memory_shadow(
        timings,
        "siqs-pilot-precheck",
        Some(2),
        columns,
        0,
        0,
        0,
        half_width,
        workers,
        memory_budget,
        pilot_collect_allow,
    );
    if !pilot_collect_allow {
        return (None, bank);
    }
    let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
    let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
    let suitable = suitable_primes(n, working.algebraic_base_bound);
    let poly_floor = first_viable_poly_index(&suitable, &target_a, factor_cap);
    if poly_floor == usize::MAX {
        return (None, bank);
    }
    let mut next_poly = bank.next_poly.max(poly_floor);
    let mut relations = std::mem::take(&mut bank.relations);
    let mut seen_masks: std::collections::HashSet<Vec<u64>> =
        relations.iter().map(|item| item.mask.clone()).collect();
    let mut lp_graph = LpGraph::new();
    let display: Option<PolynomialPair> =
        generate_polynomial_from(n, &suitable, &target_a, 0, factor_cap).map(|poly| {
            PolynomialPair {
                degree: 2,
                coefficients: vec![
                    poly.coeff_c.clone(),
                    Integer::from(2) * &poly.coeff_b,
                    poly.coeff_a.clone(),
                ],
                m: poly.coeff_b.clone(),
            }
        });
    let collect_timer = Timer::new(timings, "siqs_collect", Some(2));
    let lp_bound = u64::try_from(working.large_prime_bound).unwrap_or(u64::MAX);
    let splitters = splitter_primes(&base, lp_bound);
    let gamma_primes: Vec<u64> = base
        .rational_primes
        .iter()
        .chain(base.algebraic_modulus_primes.iter())
        .copied()
        .collect();
    let gammas = SiqsGammaTable::build(n, &gamma_primes);
    let toggles = SiqsToggleIndex::build(&base);
    let barrett_cover = gamma_primes
        .iter()
        .copied()
        .max()
        .unwrap_or(2)
        .max(working.algebraic_base_bound as u64)
        .max(2);
    let barrett = SiqsBarrettTable::build(barrett_cover);
    let mut totals = SiqsPolyCounters::default();
    let mut polys_done = 0_usize;
    let mut cpu_last = Instant::now();
    'rounds: loop {
        if cancelled.load(Ordering::Relaxed) {
            break;
        }
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            break;
        }
        if let Some(collector) = timings {
            collector.add_cpu_ms(cpu_last.elapsed().as_millis(), workers);
            cpu_last = Instant::now();
            if collector.check_cpu_budget(options).is_err() {
                break;
            }
        }
        if next_poly > SIQS_POLY_INDEX_CAP {
            break;
        }
        let batch: Vec<usize> = (0..workers)
            .map(|w| next_poly + w)
            .take_while(|&i| i <= SIQS_POLY_INDEX_CAP)
            .collect();
        if batch.is_empty() {
            break;
        }
        next_poly += batch.len();
        let mut round_out: Vec<(usize, (PolyOutput, SiqsPolyCounters))> =
            Vec::with_capacity(batch.len());
        let mut round_generated = 0_usize;
        let mut round_failed = false;
        let target_a_ref: &Integer = &target_a;
        let base_ref: &FactorBase = &base;
        let suitable_worker: &[SiqsSuitablePrime] = &suitable;
        let gammas_ref: &SiqsGammaTable = &gammas;
        let toggles_ref: &SiqsToggleIndex = &toggles;
        let barrett_ref: &SiqsBarrettTable = &barrett;
        let splitters_ref: &[u64] = &splitters;
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(batch.len());
            for &index in &batch {
                handles.push(scope.spawn(move || {
                    let gen_start = Instant::now();
                    let poly = match generate_polynomial_from(
                        n,
                        suitable_worker,
                        target_a_ref,
                        index * 3,
                        factor_cap,
                    ) {
                        Some(poly) => poly,
                        None => {
                            let mut empty = SiqsPolyCounters::default();
                            empty.generate_nanos += gen_start.elapsed().as_nanos() as u64;
                            return (
                                index,
                                false,
                                Ok((
                                    (Vec::<Relation>::new(), Vec::<(Relation, Vec<u64>)>::new()),
                                    empty,
                                )),
                            );
                        }
                    };
                    let gen_nanos = gen_start.elapsed().as_nanos() as u64;
                    let mut nop = |_: &ProgressEvent| Ok(());
                    let out = collect_one_poly(
                        n,
                        &poly,
                        index as u64,
                        base_ref,
                        gammas_ref,
                        toggles_ref,
                        barrett_ref,
                        half_width,
                        lp_bound,
                        splitters_ref,
                        cancelled,
                        deadline,
                        &mut nop,
                        columns as u64,
                    )
                    .map(|(output, mut poly_stats)| {
                        poly_stats.generate_nanos += gen_nanos;
                        (output, poly_stats)
                    });
                    (index, true, out)
                }));
            }
            for handle in handles {
                match handle.join() {
                    Ok((index, generated, Ok(out))) => {
                        if generated {
                            round_generated += 1;
                            polys_done += 1;
                        }
                        round_out.push((index, out))
                    }
                    Ok((_, _, Err(_))) => {
                        round_failed = true;
                    }
                    Err(_) => {
                        round_failed = true;
                    }
                }
            }
        });
        if round_failed {
            break 'rounds;
        }
        if round_generated == 0 {
            break;
        }
        round_out.sort_by_key(|(index, _)| *index);
        for (_, ((smooth, partials), poly_stats)) in round_out {
            totals.accumulate(&poly_stats);
            for relation in smooth {
                push_unique(&mut seen_masks, &mut relations, relation);
            }
            for (relation, lps) in partials {
                let mut combined: Vec<Relation> = Vec::new();
                lp_graph.ingest(relation, lps, n, &mut combined);
                totals.merge_count += combined.len() as u64;
                for relation in combined {
                    push_unique(&mut seen_masks, &mut relations, relation);
                }
            }
        }
        let _ = progress(&ProgressEvent {
            stage: "collect".to_string(),
            current: relations.len() as u64,
            total: columns as u64,
            message: "sieving siqs pilot polynomial".to_string(),
        });
        let live_bytes = estimate_live_bytes(columns, relations.len(), 0, 0);
        crate::factor::record_memory_shadow(
            timings,
            "siqs-pilot-live",
            Some(2),
            columns,
            relations.len(),
            0,
            0,
            half_width,
            workers,
            memory_budget,
            live_bytes <= memory_budget,
        );
        if live_bytes > memory_budget {
            break 'rounds;
        }
        if relations.len() >= columns {
            break;
        }
        if polys_done >= workers {
            break;
        }
    }
    collect_timer.finish(format!(
        "pilot polys={polys_done} relations={} banked={}",
        relations.len(),
        relations.len()
    ));
    bank.ingest_new(relations, next_poly);
    if bank.relations.len() <= 1 {
        return (None, bank);
    }
    let need_bytes = siqs_rank_estimate_bytes(columns, bank.relations.len());
    crate::factor::record_memory_shadow(
        timings,
        "siqs-pilot-rank",
        Some(2),
        columns,
        bank.relations.len(),
        0,
        0,
        half_width,
        workers,
        memory_budget,
        need_bytes <= memory_budget,
    );
    if need_bytes > memory_budget {
        return (None, bank);
    }
    let (rank, dependencies) = {
        let timer = Timer::new(timings, "siqs_rank", Some(2));
        let value = rank_and_dependencies(&bank.relations, columns, 128, 1, workers);
        timer.finish(format!("rank={} of {}", value.0, bank.relations.len()));
        value
    };
    if bank.relations.len() <= rank {
        return (None, bank);
    }
    for (di, index) in crate::linear_algebra::order_dependencies_by_cost(&dependencies)
        .iter()
        .enumerate()
    {
        let timer = Timer::new(timings, "siqs_sqrt", Some(2));
        let attempt_result = siqs_extract_relations(n, &bank.relations, &dependencies[*index]);
        let found = attempt_result.is_ok();
        timer.finish(format!("dependency={di} found={found}"));
        if let Ok(divisor) = attempt_result {
            let outcome = SiqsOutcome {
                display: match display {
                    Some(pair) => pair,
                    None => return (None, bank),
                },
                relations: bank.relations.len(),
                rank,
            };
            let hit = (divisor, outcome);
            return (Some(hit), bank);
        }
    }
    (None, bank)
}

pub(crate) fn siqs_divisor<F>(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
    timings: Option<&TimingCollector>,
) -> Result<(Integer, SiqsOutcome), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let mut working = options.clone();
    let factor_cap = resolve_siqs_table(n, &mut working);
    let memory_budget = siqs_memory_budget_bytes(options);
    let mut last_error = FactorError::SearchExhausted;
    let entry = working.clone();
    for attempt in 0..12_u32 {
        let phase = attempt % 4;
        if phase == 0 && attempt > 0 {
            let starved = matches!(last_error, FactorError::InsufficientRelations { .. });
            let expired = deadline.is_some_and(|limit| Instant::now() >= limit);
            let Some(step) = SIQS_WIDTH_ESCALATION_STEPS.get(attempt as usize / 4 - 1) else {
                break;
            };
            if !starved || expired {
                break;
            }
            working = entry.clone();
            working.sieve_a_bound = entry.sieve_a_bound.saturating_mul(*step);
        }
        if phase > 0 {
            if working.rational_base_bound < 50_000 {
                working.rational_base_bound = (working.rational_base_bound * 2).min(50_000);
            }
            if working.algebraic_base_bound < 250_000 {
                working.algebraic_base_bound = (working.algebraic_base_bound * 2).min(250_000);
            }
            if phase == 2 && working.sieve_a_bound < 500_000 {
                working.sieve_a_bound = (working.sieve_a_bound * 2).min(500_000);
            }
        }
        let base = {
            let timer = Timer::new(timings, "siqs_base", Some(2));
            let value =
                FactorBase::for_siqs(n, working.rational_base_bound, working.algebraic_base_bound);
            timer.finish(format!("columns={}", value.column_count()));
            value
        };
        let columns = base.column_count();
        let attempt_deadline = deadline.map(|limit| {
            let remaining = limit.saturating_duration_since(Instant::now());
            Instant::now() + remaining.div_f32((4_u32 - phase) as f32)
        });
        let half_width = i64::try_from(working.sieve_a_bound).map_err(|_| {
            FactorError::InvalidInput("sieve a bound exceeds the supported range".to_string())
        })?;
        if half_width <= 0 {
            return Err(FactorError::InvalidInput(
                "sieve a bound must be greater than zero".to_string(),
            ));
        }
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let mut relations: Vec<Relation> = Vec::new();
        let mut seen_masks: std::collections::HashSet<Vec<u64>> = std::collections::HashSet::new();
        let mut lp_graph = LpGraph::new();
        // X1d: suitable-prime list (with cached square roots) built once per
        // attempt and shared read-only across workers, instead of one sieve
        // plus one sqrt-chain per generated polynomial.
        let suitable = suitable_primes(n, working.algebraic_base_bound);
        let suitable_ref: &[SiqsSuitablePrime] = &suitable;
        let display: Option<PolynomialPair> =
            generate_polynomial_from(n, suitable_ref, &target_a, 0, factor_cap).map(|poly| {
                PolynomialPair {
                    degree: 2,
                    coefficients: vec![
                        poly.coeff_c.clone(),
                        Integer::from(2) * &poly.coeff_b,
                        poly.coeff_a.clone(),
                    ],
                    m: poly.coeff_b.clone(),
                }
            });
        let collect_timer = Timer::new(timings, "siqs_collect", Some(2));
        let lp_bound = u64::try_from(working.large_prime_bound).unwrap_or(u64::MAX);
        let splitters = splitter_primes(&base, lp_bound);
        // X1a: one sqrt(n) mod p per base prime per attempt; shared read-only
        // across the worker threads instead of one mod_sqrt per (poly, prime).
        let gamma_primes: Vec<u64> = base
            .rational_primes
            .iter()
            .chain(base.algebraic_modulus_primes.iter())
            .copied()
            .collect();
        let gammas = SiqsGammaTable::build(n, &gamma_primes);
        // X1b: dense prime->column table shared read-only across workers.
        let toggles = SiqsToggleIndex::build(&base);
        // X1c: Barrett mu per prime; cover polynomial factors too (they are
        // bounded by the algebraic bound by construction).
        let barrett_cover = gamma_primes
            .iter()
            .copied()
            .max()
            .unwrap_or(2)
            .max(working.algebraic_base_bound as u64)
            .max(2);
        let barrett = SiqsBarrettTable::build(barrett_cover);
        let workers = siqs_worker_count(options);
        let main_collect_allow =
            siqs_collect_estimate_bytes(columns, half_width, workers) <= memory_budget;
        crate::factor::record_memory_shadow(
            timings,
            "siqs-main-precheck",
            Some(2),
            columns,
            0,
            0,
            0,
            half_width,
            workers,
            memory_budget,
            main_collect_allow,
        );
        if !main_collect_allow {
            return Err(FactorError::ResourceLimit(
                "siqs sieve arrays exceed memory budget".to_string(),
            ));
        }
        let poly_floor = first_viable_poly_index(suitable_ref, &target_a, factor_cap);
        let mut next_poly = poly_floor;
        let mut polys_done = 0_usize;
        let mut totals = SiqsPolyCounters::default();
        let mut cpu_last = Instant::now();
        'rounds: loop {
            if cancelled.load(Ordering::Relaxed) {
                return Err(FactorError::Cancelled);
            }
            if attempt_deadline.is_some_and(|limit| Instant::now() >= limit) {
                break;
            }
            if let Some(collector) = timings {
                collector.add_cpu_ms(cpu_last.elapsed().as_millis(), workers);
                cpu_last = Instant::now();
                collector.check_cpu_budget(options)?;
            }
            if next_poly > SIQS_POLY_INDEX_CAP {
                break;
            }
            let batch: Vec<usize> = (0..workers)
                .map(|w| next_poly + w)
                .take_while(|&i| i <= SIQS_POLY_INDEX_CAP)
                .collect();
            if batch.is_empty() {
                break;
            }
            next_poly += batch.len();
            let mut round_out: Vec<(usize, (PolyOutput, SiqsPolyCounters))> =
                Vec::with_capacity(batch.len());
            let mut round_generated = 0_usize;
            let mut round_stopped = false;
            let mut round_failed: Option<FactorError> = None;
            let target_a_ref: &Integer = &target_a;
            let base_ref: &FactorBase = &base;
            let suitable_worker: &[SiqsSuitablePrime] = suitable_ref;
            let gammas_ref: &SiqsGammaTable = &gammas;
            let toggles_ref: &SiqsToggleIndex = &toggles;
            let barrett_ref: &SiqsBarrettTable = &barrett;
            let splitters_ref: &[u64] = &splitters;
            std::thread::scope(|scope| {
                let mut handles = Vec::with_capacity(batch.len());
                for &index in &batch {
                    handles.push(scope.spawn(move || {
                        let gen_start = Instant::now();
                        let poly = match generate_polynomial_from(
                            n,
                            suitable_worker,
                            target_a_ref,
                            index * 3,
                            factor_cap,
                        ) {
                            Some(poly) => poly,
                            None => {
                                let mut empty = SiqsPolyCounters::default();
                                empty.generate_nanos += gen_start.elapsed().as_nanos() as u64;
                                return (
                                    index,
                                    false,
                                    Ok((
                                        (
                                            Vec::<Relation>::new(),
                                            Vec::<(Relation, Vec<u64>)>::new(),
                                        ),
                                        empty,
                                    )),
                                );
                            }
                        };
                        let gen_nanos = gen_start.elapsed().as_nanos() as u64;
                        let mut nop = |_: &ProgressEvent| Ok(());
                        let out = collect_one_poly(
                            n,
                            &poly,
                            index as u64,
                            base_ref,
                            gammas_ref,
                            toggles_ref,
                            barrett_ref,
                            half_width,
                            lp_bound,
                            splitters_ref,
                            cancelled,
                            attempt_deadline,
                            &mut nop,
                            columns as u64,
                        )
                        .map(|(output, mut poly_stats)| {
                            poly_stats.generate_nanos += gen_nanos;
                            (output, poly_stats)
                        });
                        (index, true, out)
                    }));
                }
                for handle in handles {
                    match handle.join() {
                        Ok((index, generated, Ok(out))) => {
                            if generated {
                                round_generated += 1;
                                polys_done += 1;
                            }
                            round_out.push((index, out))
                        }
                        Ok((_, _, Err(FactorError::ResourceLimit(_)))) => {
                            round_stopped = true;
                        }
                        Ok((_, _, Err(other))) => {
                            round_failed = Some(other);
                        }
                        Err(_) => {
                            round_failed = Some(FactorError::Cancelled);
                        }
                    }
                }
            });
            if let Some(error) = round_failed {
                return Err(error);
            }
            if round_generated == 0 {
                break;
            }
            round_out.sort_by_key(|(index, _)| *index);
            for (_, ((smooth, partials), poly_stats)) in round_out {
                totals.accumulate(&poly_stats);
                for relation in smooth {
                    push_unique(&mut seen_masks, &mut relations, relation);
                }
                for (relation, lps) in partials {
                    let mut combined: Vec<Relation> = Vec::new();
                    lp_graph.ingest(relation, lps, n, &mut combined);
                    totals.merge_count += combined.len() as u64;
                    for relation in combined {
                        push_unique(&mut seen_masks, &mut relations, relation);
                    }
                }
            }
            progress(&ProgressEvent {
                stage: "collect".to_string(),
                current: relations.len() as u64,
                total: columns as u64,
                message: "sieving siqs polynomials".to_string(),
            })?;
            let live_bytes = estimate_live_bytes(columns, relations.len(), 0, 0);
            let main_live_allow = live_bytes <= memory_budget;
            crate::factor::record_memory_shadow(
                timings,
                "siqs-main-live",
                Some(2),
                columns,
                relations.len(),
                0,
                0,
                half_width,
                workers,
                memory_budget,
                main_live_allow,
            );
            if !main_live_allow {
                return Err(FactorError::ResourceLimit(format!(
                    "relation collection exceeded the memory budget of {} MiB (live estimate {} MiB)",
                    memory_budget.div_ceil(1024 * 1024),
                    live_bytes.div_ceil(1024 * 1024)
                )));
            }
            if relations.len() >= columns {
                break;
            }
            if round_stopped {
                break 'rounds;
            }
        }
        relations = {
            let timer = Timer::new(timings, "siqs_filter", Some(2));
            let (filtered, filter_stats) = if relations.len() > columns {
                filter_relations_with_budget(
                    n,
                    relations,
                    columns,
                    &SIQS_FILTER_MERGE_CAPS,
                    SIQS_FILTER_OVERWEIGHT_CAP,
                    memory_budget,
                )
            } else {
                (relations, SiqsFilterStats::default())
            };
            timer.finish(format!(
                "rows={} zero={} singletons={} merges={} overweight={}",
                filtered.len(),
                filter_stats.zero_masks_dropped,
                filter_stats.singletons_removed,
                filter_stats.merges,
                filter_stats.overweight_dropped,
            ));
            filtered
        };
        collect_timer.finish(format!(
            "attempt={attempt} polys={polys_done} relations={} sieve_mark_ms={:.2} roots_ms={:.2} gen_ms={:.2} classify_u64={}({:.2}ms) classify_bigint={}({:.2}ms) lp_splitter={}({:.2}ms) merge={} evals={}({:.2}ms) score_max={} pass_eval={}/{}/{} pass_rel={}/{}/{}",
            relations.len(),
            SiqsPolyCounters::millis(totals.sieve_mark_nanos),
            SiqsPolyCounters::millis(totals.roots_nanos),
            SiqsPolyCounters::millis(totals.generate_nanos),
            totals.classify_u64_count,
            SiqsPolyCounters::millis(totals.classify_u64_nanos),
            totals.classify_bigint_count,
            SiqsPolyCounters::millis(totals.classify_bigint_nanos),
            totals.splitter_count,
            SiqsPolyCounters::millis(totals.splitter_nanos),
            totals.merge_count,
            totals.evaluate_calls,
            SiqsPolyCounters::millis(totals.evaluate_nanos),
            totals.score_max,
            totals.pass_evaluated[0],
            totals.pass_evaluated[1],
            totals.pass_evaluated[2],
            totals.pass_relations[0],
            totals.pass_relations[1],
            totals.pass_relations[2],
        ));
        if relations.len() <= 1 {
            last_error = FactorError::InsufficientRelations {
                found: relations.len(),
                required: 2,
            };
            continue;
        }
        let need_bytes = siqs_rank_estimate_bytes(columns, relations.len());
        crate::factor::record_memory_shadow(
            timings,
            "siqs-main-rank",
            Some(2),
            columns,
            relations.len(),
            0,
            0,
            half_width,
            workers,
            memory_budget,
            need_bytes <= memory_budget,
        );
        if need_bytes > memory_budget {
            return Err(FactorError::ResourceLimit(format!(
                "the relation matrix total would need about {} MiB, above the memory budget of {} MiB",
                need_bytes.div_ceil(1024 * 1024),
                memory_budget.div_ceil(1024 * 1024)
            )));
        }
        let (rank, dependencies) = {
            let timer = Timer::new(timings, "siqs_rank", Some(2));
            let value = rank_and_dependencies(&relations, columns, 128, 1, workers);
            timer.finish(format!("rank={} of {}", value.0, relations.len()));
            value
        };
        if relations.len() <= rank {
            last_error = FactorError::InsufficientRelations {
                found: relations.len(),
                required: rank + 1,
            };
            continue;
        }
        for (di, index) in crate::linear_algebra::order_dependencies_by_cost(&dependencies)
            .iter()
            .enumerate()
        {
            let timer = Timer::new(timings, "siqs_sqrt", Some(2));
            let attempt_result = siqs_extract_relations(n, &relations, &dependencies[*index]);
            let found = attempt_result.is_ok();
            timer.finish(format!("dependency={di} found={found}"));
            if let Ok(divisor) = attempt_result {
                return Ok((
                    divisor,
                    SiqsOutcome {
                        display: display.ok_or(FactorError::NoPolynomial)?,
                        relations: relations.len(),
                        rank,
                    },
                ));
            }
        }
        last_error = FactorError::AlgebraicSquareRootFailed;
    }
    Err(last_error)
}

#[cfg(test)]
mod siqs_tests {
    use super::*;

    #[test]
    fn polynomial_satisfies_congruence() {
        let n = Integer::from(100160063);
        let target = Integer::from(1) << target_a_bits(27, 11);
        let poly = generate_polynomial(&n, 2000, &target, 0).expect("poly");
        assert!((&poly.coeff_b * &poly.coeff_b - &n) % &poly.coeff_a == Integer::from(0));
        assert_eq!(
            (&poly.coeff_b * &poly.coeff_b - &n) / &poly.coeff_a,
            poly.coeff_c
        );
        assert!(poly.coeff_c <= Integer::from(0));
        assert!(poly.coeff_a.bits() >= 2);
        assert!(!poly.prime_factors.is_empty());
    }

    #[test]
    fn x_squared_matches_a_times_value_mod_n() {
        let n = Integer::from(100160063);
        let target = Integer::from(1) << target_a_bits(27, 11);
        let poly = generate_polynomial(&n, 2000, &target, 0).expect("poly");
        for raw in [-3_i64, -1, 0, 1, 2, 7] {
            let x = Integer::from(raw);
            let left = poly.x_value(&x) * poly.x_value(&x) % &n;
            let right = (&poly.coeff_a * poly.evaluate(&x)) % &n;
            let right = (right + &n) % &n;
            assert_eq!(left, right, "x={raw}");
        }
    }

    #[test]
    fn sieve_roots_hit_the_polynomial() {
        let n = Integer::from(100160063);
        let target = Integer::from(1) << target_a_bits(27, 11);
        let poly = generate_polynomial(&n, 2000, &target, 0).expect("poly");
        let mut checked = 0;
        for prime in sieve_primes(200) {
            if prime < 3 {
                continue;
            }
            for root in poly.roots_mod(prime) {
                let value = poly.evaluate(&Integer::from(root));
                assert_eq!(value % Integer::from(prime), Integer::from(0));
                checked += 1;
            }
        }
        assert!(checked > 0);
        let _ = checked;
        for prime in poly.prime_factors.iter().copied() {
            if prime == 0 {
                continue;
            }
            assert!(!poly.roots_mod(prime).is_empty());
        }
    }

    #[test]
    fn cached_gamma_roots_match_direct_roots() {
        // X1a proof: the per-attempt gamma table reproduces roots_mod
        // bit-identically on every (poly, prime) pair, including the
        // a_mod == 0 linear path (whose factor primes divide coeff_a).
        let n = Integer::from(11_783_510_925_019_u64);
        let half_width = 700_i64;
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let gamma_primes: Vec<u64> = base
            .rational_primes
            .iter()
            .chain(base.algebraic_modulus_primes.iter())
            .copied()
            .collect();
        let gammas = SiqsGammaTable::build(&n, &gamma_primes);
        let mut checked = 0_usize;
        for skip in [0_usize, 3, 6, 9, 12] {
            let Some(poly) = generate_polynomial(&n, 2000, &target_a, skip) else {
                continue;
            };
            let mut primes: Vec<u64> = gamma_primes.clone();
            primes.extend(poly.prime_factors.iter().copied());
            primes.sort_unstable();
            primes.dedup();
            for prime in primes {
                if prime < 3 {
                    continue;
                }
                let direct = poly.roots_mod(prime);
                let cached = poly.roots_mod_with_gamma(prime, gammas.get(prime));
                assert_eq!(direct, cached, "skip={skip} prime={prime}");
                // Every table prime resolves identically through both tables.
                assert_eq!(
                    gammas.raw(prime),
                    Some(mod_sqrt(
                        (&n % Integer::from(prime)).to_u64_saturating(),
                        prime
                    ))
                );
                checked += 1;
            }
        }
        assert!(checked > 100, "need coverage");
    }

    #[test]
    fn barrett_divmod_matches_native_division() {
        // X1c proof: Barrett divmod is exactly native division, over small
        // primes x edge values plus 50k pseudo-random (value, prime) pairs
        // (debug builds additionally re-check every call inside
        // barrett_divmod itself via debug_assert_eq).
        let primes = sieve_primes(2000);
        let table = SiqsBarrettTable::build(2000);
        let mut checked = 0_u64;
        for &prime in &primes {
            if prime < 2 {
                continue;
            }
            let mut values = vec![
                0_u64,
                1,
                prime - 1,
                prime,
                prime + 1,
                u64::MAX,
                u64::MAX - 1,
                u64::MAX / prime,
            ];
            let mut state = prime.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
            for _ in 0..16 {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                values.push(state);
            }
            for value in values {
                assert_eq!(table.divmod(value, prime), (value / prime, value % prime));
                checked += 1;
            }
        }
        assert!(checked > 5_000);
        // Out-of-range primes fall back to native division exactly.
        for (value, prime) in [(u64::MAX, 2003), (123456789, 5003), (0, 999983)] {
            assert_eq!(table.divmod(value, prime), (value / prime, value % prime));
        }
    }

    #[test]
    fn barrett_strip_loop_reports_timing() {
        // X1c supporting evidence: time a native strip loop vs the Barrett
        // strip loop over identical pseudo-random workloads. Informational
        // only (no time assertion; CI noise); the verdict proof is the
        // exactness test above plus unchanged classify outcomes.
        let primes = sieve_primes(50_000);
        let table = SiqsBarrettTable::build(50_000);
        let mut state = 0x123456789ABCDEF_u64;
        let mut workloads = Vec::with_capacity(200_000);
        for _ in 0..200_000 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let prime = primes[(state as usize) % primes.len()].max(2);
            workloads.push((state | 1, prime));
        }
        let start = Instant::now();
        let mut native_stripped = 0_u64;
        for &(value, prime) in &workloads {
            let mut rest = value;
            while rest.is_multiple_of(prime) {
                rest /= prime;
                native_stripped += 1;
            }
        }
        let native_ms = start.elapsed().as_secs_f64() * 1000.0;
        let start = Instant::now();
        let mut barrett_stripped = 0_u64;
        for &(value, prime) in &workloads {
            let mut rest = value;
            loop {
                let (quotient, remainder) = table.divmod(rest, prime);
                if remainder != 0 {
                    break;
                }
                rest = quotient;
                barrett_stripped += 1;
            }
        }
        let barrett_ms = start.elapsed().as_secs_f64() * 1000.0;
        assert_eq!(native_stripped, barrett_stripped);
        if std::env::var("GRAYNFS_VERBOSE").is_ok() {
            eprintln!("X1c strip-bench native_ms={native_ms:.2} barrett_ms={barrett_ms:.2}");
        }
    }

    #[test]
    fn dense_toggle_matches_hashmap_toggle() {
        // X1b proof: the dense prime->column table toggles exactly the bits
        // the HashMap path toggles, for rational primes, algebraic-only
        // primes, unknown primes, and out-of-table fallback primes.
        let n = Integer::from(11_783_510_925_019_u64);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let toggles = SiqsToggleIndex::build(&base);
        let words = base.column_count().div_ceil(64);
        let mut primes: Vec<u64> = base
            .rational_primes
            .iter()
            .chain(base.algebraic_modulus_primes.iter())
            .copied()
            .collect();
        primes.extend([2, 3, 5, 999_983, 1_000_003, 50_000_007]);
        primes.sort_unstable();
        primes.dedup();
        for prime in primes {
            let mut slow = vec![0_u64; words];
            toggle_parity(&base, &mut slow, prime);
            let mut fast = vec![0_u64; words];
            toggle_parity_fast(&base, &toggles, &mut fast, prime);
            assert_eq!(slow, fast, "prime={prime}");
        }
    }

    #[test]
    fn hoisted_generation_matches_legacy_sequence() {
        // X1d proof: building the suitable list once per attempt reproduces
        // the legacy per-index sequence exactly (same a/b/c/factors for the
        // first N indices, same None pattern past the list end), because the
        // sieve order, filters, skip cap, greedy take, and CRT are unchanged.
        let n = Integer::from(11_783_510_925_019_u64);
        let half_width = 700_i64;
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let suitable = suitable_primes(&n, 2000);
        assert!(!suitable.is_empty());
        let mut some_count = 0_usize;
        for skip in 0..64_usize {
            let legacy = generate_polynomial(&n, 2000, &target_a, skip);
            let hoisted = generate_polynomial_from(&n, &suitable, &target_a, skip, None);
            match (legacy, hoisted) {
                (None, None) => {}
                (Some(left), Some(right)) => {
                    assert_eq!(left.coeff_a, right.coeff_a, "skip={skip}");
                    assert_eq!(left.coeff_b, right.coeff_b, "skip={skip}");
                    assert_eq!(left.coeff_c, right.coeff_c, "skip={skip}");
                    assert_eq!(left.prime_factors, right.prime_factors, "skip={skip}");
                    assert!(
                        (&right.coeff_b * &right.coeff_b - &n) % &right.coeff_a == Integer::from(0)
                    );
                    some_count += 1;
                }
                _ => panic!("legacy/hoisted mismatch at skip={skip}"),
            }
        }
        assert!(some_count > 32, "need sequence coverage");
        // Past the suitable-list end both spellings agree on None.
        assert!(generate_polynomial(&n, 2000, &target_a, 10_000_000).is_none());
        assert!(generate_polynomial_from(&n, &suitable, &target_a, 10_000_000, None).is_none());
    }

    /// Independent full-CRT reference (same legacy formula, separate path
    /// from the flip): used to prove the O(1) update is bit-exact.
    #[cfg(test)]
    fn full_crt_reference(
        n: &Integer,
        factors: &[u64],
        gamma_of: &dyn Fn(u64) -> Option<u64>,
        target_a: &Integer,
    ) -> Option<SiqsPolynomial> {
        let mut coeff_a = Integer::from(1);
        for prime in factors {
            coeff_a *= Integer::from(*prime);
        }
        if factors.len() < 2 || &coeff_a < target_a {
            return None;
        }
        let mut coeff_b = Integer::from(0);
        for prime in factors {
            let modulus = Integer::from(*prime);
            let gamma = gamma_of(*prime)?;
            let rest = &coeff_a / &modulus;
            let inv = egcd_inverse(&rest, &modulus)?;
            coeff_b += Integer::from(gamma) * rest * inv;
        }
        coeff_b %= &coeff_a;
        if &coeff_b % Integer::from(2) == Integer::from(0) {
            coeff_b = &coeff_a - &coeff_b;
        }
        let coeff_c = (&coeff_b * &coeff_b - n) / &coeff_a;
        if (&coeff_b * &coeff_b - n) % &coeff_a != Integer::from(0) {
            return None;
        }
        Some(SiqsPolynomial {
            coeff_a,
            coeff_b,
            coeff_c,
            prime_factors: factors.to_vec(),
        })
    }

    #[test]
    fn gray_flip_matches_crt_characterization() {
        // X1d proof: N chained single-prime flips each satisfy the exact
        // CRT characterization (invariant b^2 == n mod A, surviving
        // congruences preserved, new-factor congruence installed, exact
        // a'/c' relations, odd b'), hence equal a full O(k) recomputation.
        let n = Integer::from(11_783_510_925_019_u64);
        let half_width = 700_i64;
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let suitable = suitable_primes(&n, 2000);
        let gamma_of = |prime: u64| {
            suitable
                .iter()
                .find(|entry| entry.prime == prime)
                .map(|entry| entry.gamma)
        };
        let mut poly = generate_polynomial_from(&n, &suitable, &target_a, 0, None).expect("poly 0");
        for (pool_cursor, step) in (0..32_usize).enumerate() {
            let out_prime = poly.prime_factors[step % poly.prime_factors.len()];
            let in_entry = *suitable
                .iter()
                .skip(pool_cursor)
                .find(|entry| !poly.prime_factors.contains(&entry.prime))
                .expect("fresh prime");
            let in_gamma = gamma_of(in_entry.prime).expect("gamma");
            // Detect the chain's uniform CRT sign from a kept factor, so the
            // flip installs the matching conjugate and stays bit-exact.
            let kept = *poly
                .prime_factors
                .iter()
                .find(|&&prime| prime != out_prime)
                .expect("kept factor");
            let kept_gamma = gamma_of(kept).expect("gamma");
            let kept_mod = Integer::from(kept);
            let kept_residue = &poly.coeff_b % &kept_mod;
            let sign = if kept_residue == Integer::from(kept_gamma) {
                1
            } else {
                assert_eq!(
                    kept_residue,
                    Integer::from((kept - kept_gamma) % kept),
                    "step={step}"
                );
                -1
            };
            let next = gray_flip_polynomial(&n, &poly, out_prime, in_entry.prime, in_gamma, sign)
                .expect("flip");
            // Invariant: b^2 == n (mod A) with exact c.
            assert_eq!(
                (&next.coeff_b * &next.coeff_b - &n) % &next.coeff_a,
                Integer::from(0),
                "step={step}"
            );
            assert_eq!(
                (&next.coeff_b * &next.coeff_b - &n) / &next.coeff_a,
                next.coeff_c,
                "step={step}"
            );
            // Exact factor swap.
            assert_eq!(
                next.coeff_a,
                &poly.coeff_a / Integer::from(out_prime) * Integer::from(in_entry.prime),
                "step={step}"
            );
            // Surviving congruences preserved up to the single global parity
            // outcome (shared with the reference by the sign-aware install).
            let half = &poly.coeff_a / Integer::from(out_prime);
            let diff = &next.coeff_b - &poly.coeff_b;
            let sum = &next.coeff_b + &poly.coeff_b;
            let diff_ok = ((&diff % &half) + &half) % &half == Integer::from(0);
            let sum_ok = ((&sum % &half) + &half) % &half == Integer::from(0);
            assert!(diff_ok || sum_ok, "step={step}");
            // Legacy parity convention (odd b).
            assert!(&next.coeff_b % Integer::from(2) != Integer::from(0));
            // SIQS validity: x^2 == A * Q(x) (mod n) at sample positions.
            for raw in [-3_i64, -1, 0, 1, 2, 7] {
                let x = Integer::from(raw);
                let left = next.x_value(&x) * next.x_value(&x) % &n;
                let right = ((&next.coeff_a * next.evaluate(&x)) % &n + &n) % &n;
                assert_eq!(left, right, "step={step} x={raw}");
            }
            // Bit-exact against the independent full-CRT reference.
            let reference = full_crt_reference(&n, &next.prime_factors, &gamma_of, &target_a)
                .expect("reference");
            assert_eq!(next.coeff_a, reference.coeff_a, "step={step}");
            assert_eq!(next.coeff_b, reference.coeff_b, "step={step}");
            assert_eq!(next.coeff_c, reference.coeff_c, "step={step}");
            poly = next;
        }
    }

    #[test]
    fn hit_classify_matches_full_base() {
        let n = Integer::from(11_783_510_925_019_u64);
        let half_width = 700_i64;
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let splitters = splitter_primes(&base, 50_000);
        let width = (half_width as usize) * 2 + 1;
        let mut total_checked = 0_usize;
        let mut total_disagreements = 0_usize;
        for skip in [0_usize, 3, 6, 9, 12] {
            let Some(poly) = generate_polynomial(&n, 2000, &target_a, skip) else {
                continue;
            };
            let gamma_primes: Vec<u64> = base
                .rational_primes
                .iter()
                .chain(base.algebraic_modulus_primes.iter())
                .copied()
                .collect();
            let gammas = SiqsGammaTable::build(&n, &gamma_primes);
            let roots = sieve_roots_for(&base, &poly, &gammas);
            let mut hits: Vec<Vec<u64>> = vec![Vec::new(); width];
            for (prime, hit_roots, _) in &roots {
                for root in hit_roots {
                    let start = ((*root as i64 + half_width).rem_euclid(*prime as i64)) as usize;
                    let mut index = start;
                    while index < width {
                        let slot = &mut hits[index];
                        if slot.last().copied() != Some(*prime) {
                            slot.push(*prime);
                        }
                        index += *prime as usize;
                    }
                }
            }
            let words = base.column_count().div_ceil(64);
            let mut checked = 0_usize;
            let mut disagreements = 0_usize;
            let mut first_mismatches: Vec<String> = Vec::new();
            for offset in (0..width).step_by(7) {
                let x = Integer::from(offset as i64 - half_width);
                let value = poly.evaluate(&x);
                if value == Integer::from(0) {
                    continue;
                }
                let signed = &poly.coeff_a * &value;
                let mut mask_hit = vec![0_u64; words];
                let hit_out = classify_siqs(
                    &signed,
                    &base,
                    &mut mask_hit,
                    50_000,
                    &hits[offset],
                    &poly.prime_factors,
                    &splitters,
                );
                let mut mask_full = vec![0_u64; words];
                let full_primes: Vec<u64> = {
                    let mut seen = std::collections::HashSet::new();
                    let mut list = vec![2_u64];
                    for p in base
                        .rational_primes
                        .iter()
                        .chain(base.algebraic_modulus_primes.iter())
                    {
                        if seen.insert(*p) {
                            list.push(*p);
                        }
                    }
                    list.sort_unstable();
                    list
                };
                let full_out = classify_siqs(
                    &signed,
                    &base,
                    &mut mask_full,
                    50_000,
                    &full_primes,
                    &[],
                    &splitters,
                );
                let same_outcome = std::mem::discriminant(&hit_out)
                    == std::mem::discriminant(&full_out)
                    && match (&hit_out, &full_out) {
                        (Classify::Partial(a), Classify::Partial(b)) => a == b,
                        _ => true,
                    };
                checked += 1;
                let mask_matters = matches!(
                    (&hit_out, &full_out),
                    (Classify::Smooth, Classify::Smooth)
                        | (Classify::Partial(_), Classify::Partial(_))
                );
                if !same_outcome || (mask_matters && mask_hit != mask_full) {
                    disagreements += 1;
                    if first_mismatches.len() < 3 {
                        first_mismatches
                            .push(format!("x={offset} hit={hit_out:?} full={full_out:?}"));
                    }
                }
                if checked >= 400 {
                    break;
                }
            }
            assert!(checked > 100, "need coverage");
            assert_eq!(
                disagreements,
                0,
                "hit vs full-base mismatch: {}",
                first_mismatches.join("; ")
            );
            total_checked += checked;
            total_disagreements += disagreements;
        }
        assert!(total_checked > 500, "need multi-poly coverage");
        assert_eq!(total_disagreements, 0, "hit vs full-base mismatch");
    }

    #[cfg(test)]
    fn classify_siqs_bigint_oracle(
        signed: &Integer,
        base: &FactorBase,
        mask: &mut [u64],
        lp_bound: u64,
        trial_primes: &[u64],
        always_trial: &[u64],
        splitters: &[u64],
    ) -> Classify {
        if signed == &Integer::from(0) {
            return Classify::Reject;
        }
        let toggles = SiqsToggleIndex::build(base);
        let negative = signed < &Integer::from(0);
        if negative {
            mask[base.sign_column / 64] ^= 1_u64 << (base.sign_column % 64);
        }
        let magnitude = if negative {
            -signed.clone()
        } else {
            signed.clone()
        };
        if magnitude.to_u64().is_some() {
            panic!("oracle pair must enter the BigInt arm");
        }
        let mut residual = magnitude;
        for prime in std::iter::once(&2_u64)
            .chain(always_trial.iter())
            .chain(trial_primes.iter())
        {
            if residual == Integer::from(1) {
                break;
            }
            let prime_value = Integer::from(*prime);
            if residual < prime_value {
                continue;
            }
            let mut exponent = 0_u32;
            while (&residual % &prime_value) == Integer::from(0) {
                residual /= &prime_value;
                exponent += 1;
            }
            if exponent & 1 != 0 {
                toggle_parity_fast(base, &toggles, mask, *prime);
            }
        }
        if residual == Integer::from(1) {
            return Classify::Smooth;
        }
        match residual.to_u64() {
            Some(small) => {
                let mut discard = SiqsPolyCounters::default();
                split_residual_counted(small, splitters, lp_bound, &mut discard)
            }
            None => Classify::Reject,
        }
    }

    #[test]
    fn demoted_classify_matches_bigint_oracle() {
        let same_outcome = |left: &Classify, right: &Classify| -> bool {
            match (left, right) {
                (Classify::Smooth, Classify::Smooth) => true,
                (Classify::Partial(a), Classify::Partial(b)) => a == b,
                (Classify::Partial2(a, b), Classify::Partial2(c, d)) => a == c && b == d,
                (Classify::Reject, Classify::Reject) => true,
                _ => false,
            }
        };
        let n = Integer::from(100160063);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let words = base.column_count().div_ceil(64);
        let splitters = splitter_primes(&base, 50_000);
        let big: Integer = Integer::from(1) << 70_u64;
        let over = (Integer::from(1) << 65_u64) + Integer::from(165);
        let pairs: Vec<(Integer, Vec<u64>, Vec<u64>)> = vec![
            (big.clone() * Integer::from(15), vec![3, 5, 999_983], vec![]),
            (
                (Integer::from(1) << 80_u64) * Integer::from(3 * 10_007),
                vec![2, 3, 5],
                vec![],
            ),
            (
                big.clone() * Integer::from(10_007 * 10_009),
                vec![2, 3],
                vec![],
            ),
            (
                big.clone() * Integer::from(10_007 * 10_009) * Integer::from(10_037),
                vec![2, 3],
                vec![],
            ),
            (
                big.clone() * Integer::from(1_000_003) * Integer::from(1_000_033),
                vec![2, 3, 5],
                vec![],
            ),
            (
                -((Integer::from(1) << 80_u64) * Integer::from(3 * 10_007)),
                vec![2, 3, 5],
                vec![],
            ),
            (over * Integer::from(3), vec![3, 5, 7], vec![]),
            (
                (Integer::from(1) << 61_u64) * Integer::from(15),
                vec![3, 5, 999_983],
                vec![],
            ),
        ];
        assert!(!pairs.is_empty());
        for (signed, trial, always) in &pairs {
            let mut mask_got = vec![0_u64; words];
            let got = classify_siqs(
                signed,
                &base,
                &mut mask_got,
                50_000,
                trial,
                always,
                &splitters,
            );
            let mut mask_want = vec![0_u64; words];
            let want = classify_siqs_bigint_oracle(
                signed,
                &base,
                &mut mask_want,
                50_000,
                trial,
                always,
                &splitters,
            );
            assert!(
                same_outcome(&got, &want),
                "outcome drift signed-bits={} got={got:?} want={want:?}",
                signed.bits()
            );
            assert_eq!(mask_got, mask_want, "mask drift");
        }
        let entry = (Integer::from(1) << 61_u64) * Integer::from(15);
        assert!(entry.to_u64().is_none());
        let shrunk = &entry / Integer::from(1_u64 << 61);
        assert_eq!(shrunk, Integer::from(15));
        assert!(shrunk.to_u64().is_some());
    }

    #[test]
    fn mod_sqrt_complete_on_small_primes() {
        use crate::prime::{legendre, mod_sqrt};
        for prime in [3_u64, 5, 7, 11, 13, 17, 19, 101] {
            for residue in 0..prime {
                let is_qr = (0..prime).any(|r| (r * r) % prime == residue);
                let leg = legendre(residue, prime);
                let root = mod_sqrt(residue, prime);
                match (is_qr, residue == 0) {
                    (_, true) => assert_eq!(root, Some(0), "p={prime} r=0"),
                    (true, false) => {
                        assert_eq!(leg, 1, "p={prime} r={residue}");
                        let r = root.unwrap_or_else(|| panic!("p={prime} r={residue} None on QR"));
                        assert_eq!((r * r) % prime, residue, "p={prime} r={residue} wrong root");
                    }
                    (false, false) => {
                        assert_eq!(leg, -1, "p={prime} r={residue}");
                        assert_eq!(root, None, "p={prime} r={residue} root on non-QR");
                    }
                }
            }
        }
    }

    #[test]
    fn siqs_extract_recovers_factor() {
        let n = Integer::from(35);
        let found = siqs_extract(&n, &[Integer::from(6)], &[Integer::from(1)]).expect("factor");
        assert!(&n % &found == Integer::from(0));
        assert!(found > Integer::from(1) && found < n.clone());
    }

    #[test]
    fn siqs_extract_rejects_non_square() {
        let n = Integer::from(35);
        assert!(siqs_extract(&n, &[Integer::from(6)], &[Integer::from(2)]).is_err());
    }

    #[test]
    fn classify_detects_single_large_prime() {
        let n = Integer::from(100160063);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let words = base.column_count().div_ceil(64);
        let lp = 10_007_u64;
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(2 * 3 * lp);
        match classify_siqs(&signed, &base, &mut mask, 50_000, &[2, 3, 5], &[], &[]) {
            Classify::Partial(found) => assert_eq!(found, lp),
            other => panic!("expected partial, got {other:?}"),
        }
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(2 * 3 * 5);
        assert!(matches!(
            classify_siqs(&signed, &base, &mut mask, 50_000, &[2, 3, 5], &[], &[]),
            Classify::Smooth
        ));
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(10_007 * 10_009);
        assert!(matches!(
            classify_siqs(&signed, &base, &mut mask, 50_000, &[2, 3, 5], &[], &[]),
            Classify::Reject
        ));
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(2 * 3 * 10_000_000_019_u64);
        assert!(matches!(
            classify_siqs(&signed, &base, &mut mask, 50_000, &[2, 3, 5], &[], &[]),
            Classify::Reject
        ));
    }

    #[test]
    fn classify_detects_double_large_prime() {
        let n = Integer::from(100160063);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let words = base.column_count().div_ceil(64);
        let splitters = splitter_primes(&base, 50_000);
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(2 * 3 * 10_007 * 10_009);
        match classify_siqs(&signed, &base, &mut mask, 50_000, &[2, 3], &[], &splitters) {
            Classify::Partial2(a, b) => assert_eq!((a, b), (10_007, 10_009)),
            other => panic!("expected partial2, got {other:?}"),
        }
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(2 * 10_007 * 10_007);
        assert!(matches!(
            classify_siqs(&signed, &base, &mut mask, 50_000, &[2], &[], &splitters),
            Classify::Smooth
        ));
        let mut mask = vec![0_u64; words];
        let signed = Integer::from(10_007_u64 * 10_009 * 10_037);
        assert!(matches!(
            classify_siqs(&signed, &base, &mut mask, 50_000, &[2], &[], &splitters),
            Classify::Reject
        ));
    }

    #[test]
    fn lp_cycles_combine_to_valid_relations() {
        let n = Integer::from(1_000_003);
        let mut graph = LpGraph::new();
        let mut out: Vec<Relation> = Vec::new();
        let made = |a: u64, v: u64, word: u64| Relation {
            a: Integer::from(a),
            b: Integer::from(0_u64),
            rational_value: Integer::from(v),
            algebraic_value: QuadraticElement {
                x: Integer::from(a),
                y: Integer::from(0),
            },
            mask: vec![word],
        };
        graph.ingest(made(2, 4, 0b01), vec![7, 11], &n, &mut out);
        assert!(out.is_empty());
        graph.ingest(made(2, 4, 0b01), vec![11, 13], &n, &mut out);
        assert!(out.is_empty());
        graph.ingest(made(2, 4, 0b01), vec![13, 7], &n, &mut out);
        assert_eq!(out.len(), 1);
        let combined = &out[0];
        assert_eq!(combined.mask, vec![0b01]);
        assert_eq!(combined.a, Integer::from(8));
        assert_eq!(combined.rational_value, Integer::from(64));
        assert_eq!(
            (&combined.a * &combined.a) % &n,
            &combined.rational_value % &n
        );
        let mut paired: Vec<Relation> = Vec::new();
        graph.ingest(made(2, 4, 0b11), vec![17, 0], &n, &mut paired);
        assert!(paired.is_empty());
        graph.ingest(made(2, 4, 0b01), vec![17, 0], &n, &mut paired);
        assert_eq!(paired.len(), 1);
        assert_eq!(paired[0].mask, vec![0b10]);
    }

    fn filter_made(a: u64, v: u64, word: u64) -> Relation {
        Relation {
            a: Integer::from(a),
            b: Integer::from(0_u64),
            rational_value: Integer::from(v),
            algebraic_value: QuadraticElement {
                x: Integer::from(a),
                y: Integer::from(0),
            },
            mask: vec![word],
        }
    }

    #[test]
    fn filter_removes_singleton_rows_and_keeps_dependencies() {
        let n = Integer::from(1_000_003);
        let rows = vec![
            filter_made(2, 4, 0b0011),
            filter_made(3, 9, 0b0110),
            filter_made(5, 25, 0b0101),
            filter_made(7, 49, 0b1000),
        ];
        let (_, before) = rank_and_dependencies(&rows, 4, 128, 1, 1);
        let (filtered, stats) = filter_relations(&n, rows, 4, &[], 256);
        assert_eq!(stats.singletons_removed, 1);
        assert_eq!(stats.merges, 0);
        assert_eq!(filtered.len(), 3);
        let (_, after) = rank_and_dependencies(&filtered, 4, 128, 1, 1);
        assert_eq!(before.len(), after.len());
        assert_eq!(after, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn filter_merges_weight_two_columns() {
        let n = Integer::from(1_000_003);
        let rows = vec![
            filter_made(2, 4, 0b0011),
            filter_made(3, 9, 0b0110),
            filter_made(5, 25, 0b0101),
            filter_made(7, 49, 0b1001),
        ];
        let (_, before) = rank_and_dependencies(&rows, 4, 128, 1, 1);
        let (filtered, stats) = filter_relations(&n, rows, 4, &[64], 256);
        assert_eq!(stats.singletons_removed, 1);
        assert_eq!(stats.merges, 1);
        assert_eq!(filtered.len(), 2);
        for relation in &filtered {
            assert_eq!(relation.mask, vec![0b0110]);
        }
        let merged = &filtered[0];
        assert_eq!(merged.a, Integer::from(10));
        assert_eq!(merged.rational_value, Integer::from(100));
        let (_, after) = rank_and_dependencies(&filtered, 4, 128, 1, 1);
        assert_eq!(before.len(), after.len());
    }

    #[test]
    fn filter_respects_merge_weight_cap() {
        let n = Integer::from(1_000_003);
        let rows = vec![
            filter_made(2, 4, 0b0011),
            filter_made(3, 9, 0b0110),
            filter_made(5, 25, 0b0101),
            filter_made(7, 49, 0b1001),
        ];
        let (filtered, stats) = filter_relations(&n, rows, 4, &[0], 256);
        assert_eq!(stats.merges, 0);
        assert_eq!(filtered.len(), 3);
    }

    #[test]
    fn filter_drops_zero_masks() {
        let n = Integer::from(1_000_003);
        let rows = vec![
            filter_made(2, 4, 0b0000),
            filter_made(3, 9, 0b0011),
            filter_made(5, 25, 0b0110),
            filter_made(7, 49, 0b0101),
        ];
        let (filtered, stats) = filter_relations(&n, rows, 4, &[], 256);
        assert_eq!(stats.zero_masks_dropped, 1);
        assert_eq!(filtered.len(), 3);
    }

    #[test]
    fn filter_empty_input_stays_empty() {
        let n = Integer::from(1_000_003);
        let (filtered, stats) = filter_relations(&n, Vec::new(), 10, &[64], 256);
        assert!(filtered.is_empty());
        assert_eq!(stats, SiqsFilterStats::default());
    }

    #[test]
    fn filter_identical_pair_survives_for_rank() {
        let n = Integer::from(1_000_003);
        let rows = vec![filter_made(2, 4, 0b0011), filter_made(3, 9, 0b0011)];
        let (filtered, stats) = filter_relations(&n, rows, 4, &[64], 256);
        assert_eq!(stats.merges, 0);
        assert_eq!(filtered.len(), 2);
        let (_, dependencies) = rank_and_dependencies(&filtered, 4, 128, 1, 1);
        assert_eq!(dependencies, vec![vec![0, 1]]);
    }

    #[test]
    fn filter_drops_overweight_row_and_keeps_dependencies() {
        let n = Integer::from(1_000_003);
        let mut rows = vec![filter_made(2, 4, 0b1111)];
        for a in 3_u64..14 {
            rows.push(filter_made(a, a * a, 0b0110));
        }
        rows.push(filter_made(20, 400, 0b1001));
        let (filtered, stats) = filter_relations(&n, rows, 4, &[], 3);
        assert_eq!(stats.singletons_removed, 0);
        assert_eq!(stats.merges, 0);
        assert_eq!(stats.overweight_dropped, 1);
        assert_eq!(filtered.len(), 12);
        assert!(
            filtered
                .iter()
                .all(|relation| relation.mask != vec![0b1111])
        );
        let (after_rank, after) = rank_and_dependencies(&filtered, 4, 128, 1, 1);
        assert_eq!(after_rank, 2);
        assert_eq!(after.len(), 10);
        for dependency in &after {
            let mut word = 0_u64;
            for &index in dependency {
                word ^= filtered[index].mask[0];
            }
            assert_eq!(word, 0);
        }
    }

    #[test]
    fn table_reproduces_defaults_at_small_digits() {
        let entry = siqs_table_bounds(30);
        let defaults = FactorOptions::default();
        assert_eq!(entry.rational_bound, defaults.rational_base_bound);
        assert_eq!(entry.algebraic_bound, defaults.algebraic_base_bound);
        assert_eq!(entry.half_interval, defaults.sieve_a_bound);
        assert_eq!(entry.lp_bound, defaults.large_prime_bound);
        assert_eq!(entry.a_factors, None);
        assert_eq!(siqs_table_bounds(SIQS_TABLE_SMALL_DIGITS).a_factors, None);
    }

    #[test]
    fn table_matches_msieve_reference_at_80_digits() {
        let entry = siqs_table_bounds(80);
        assert_eq!(entry.half_interval * 2, 393_216);
        assert_eq!(entry.lp_bound, 67_108_864);
        assert_eq!(entry.a_factors, Some(10));
        assert_eq!(entry.algebraic_bound, 300_000);
    }

    #[test]
    fn table_ninety_digit_row_is_memory_bounded() {
        let entry = siqs_table_bounds(90);
        assert_eq!(entry.a_factors, Some(10));
        let options = FactorOptions::default();
        let budget = siqs_memory_budget_bytes(&options);
        let columns = 70_000_usize;
        let estimate = siqs_collect_estimate_bytes(columns, entry.half_interval as i64, 16);
        assert!(estimate < budget);
        let overflow = siqs_rank_estimate_bytes(200_000, 200_000);
        assert!(overflow > budget);
    }

    #[test]
    fn explicit_memory_budget_overrides_default() {
        let options = FactorOptions {
            max_memory_mb: Some(512),
            ..FactorOptions::default()
        };
        assert_eq!(siqs_memory_budget_bytes(&options), 512 * 1024 * 1024);
        assert_eq!(
            siqs_memory_budget_bytes(&FactorOptions::default()),
            crate::options::DEFAULT_MEMORY_BUDGET_BYTES
        );
    }

    #[test]
    fn zero_memory_budget_rejects_siqs_collect() {
        let n = Integer::from(100160063);
        let options = FactorOptions {
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        let cancelled = AtomicBool::new(false);
        let mut progress = |_: &ProgressEvent| Ok(());
        let outcome = siqs_divisor(&n, &options, &cancelled, None, &mut progress, None);
        assert!(
            matches!(outcome, Err(FactorError::ResourceLimit(_))),
            "zero budget must fail before sieving"
        );
    }

    #[test]
    fn table_respects_explicit_bounds() {
        let n = Integer::from(770817968652704265841785071059_u128);
        let mut options = FactorOptions {
            rational_base_bound: 40_000,
            ..FactorOptions::default()
        };
        let cap = resolve_siqs_table(&n, &mut options);
        assert_eq!(options.rational_base_bound, 40_000);
        assert_eq!(options.algebraic_base_bound, 2_000);
        assert_eq!(cap, None);
    }

    #[test]
    fn tiny_filter_budget_skips_merge_pass() {
        let n = Integer::from(1_000_003);
        let rows = vec![
            filter_made(2, 4, 0b0011),
            filter_made(3, 9, 0b0110),
            filter_made(5, 25, 0b0101),
            filter_made(7, 49, 0b1001),
        ];
        let (merged, _) = filter_relations(&n, rows.clone(), 4, &[64], 256);
        assert_eq!(merged.len(), 2);
        let (skipped, stats) = filter_relations_with_budget(&n, rows, 4, &[64], 256, 1);
        assert_eq!(skipped.len(), 4);
        assert_eq!(stats, SiqsFilterStats::default());
    }

    #[test]
    fn poly_floor_is_zero_at_canonical_scale() {
        let n = Integer::from(11_783_510_925_019_u64);
        let half_width = 700_i64;
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let suitable = suitable_primes(&n, 2000);
        assert!(!suitable.is_empty());
        assert_eq!(first_viable_poly_index(&suitable, &target_a, None), 0);
        assert!(window_reaches_target(&suitable, 0, &target_a, None));
    }

    #[test]
    fn poly_floor_skips_dead_windows_under_factor_cap() {
        let n = Integer::from(11_783_510_925_019_u64);
        let suitable = suitable_primes(&n, 2000);
        let scaled_target = Integer::from(1) << 51;
        let floor = first_viable_poly_index(&suitable, &scaled_target, Some(8));
        assert!(floor > 0 && floor <= SIQS_POLY_INDEX_CAP);
        assert!(window_reaches_target(
            &suitable,
            floor * 3,
            &scaled_target,
            Some(8)
        ));
        assert!(!window_reaches_target(
            &suitable,
            (floor - 1) * 3,
            &scaled_target,
            Some(8)
        ));
        let viable = generate_polynomial_from(&n, &suitable, &scaled_target, floor * 3, Some(8));
        assert!(viable.is_some());
        let dead = generate_polynomial_from(&n, &suitable, &scaled_target, 0, Some(8));
        assert!(dead.is_none());
    }

    #[test]
    fn poly_floor_reports_unreachable_target() {
        let n = Integer::from(11_783_510_925_019_u64);
        let suitable = suitable_primes(&n, 2000);
        let unreachable = Integer::from(1) << 400;
        assert_eq!(
            first_viable_poly_index(&suitable, &unreachable, Some(10)),
            usize::MAX
        );
        assert!(!window_reaches_target(&suitable, 0, &unreachable, Some(10)));
    }

    #[test]
    fn banked_collection_resumes_past_previous_poly() {
        use crate::algebraic::QuadraticElement;
        use crate::relations::Relation;
        fn made(a: u64, word: u64) -> Relation {
            Relation {
                a: Integer::from(a),
                b: Integer::from(0_u64),
                rational_value: Integer::from(a),
                algebraic_value: QuadraticElement {
                    x: Integer::from(a),
                    y: Integer::from(0),
                },
                mask: vec![word],
            }
        }
        let mut bank = SiqsBank::new();
        bank.ingest_new(vec![made(2, 0b01), made(3, 0b10)], 5);
        assert_eq!(bank.len(), 2);
        assert_eq!(bank.next_poly, 5);
        bank.ingest_new(vec![made(5, 0b100), made(2, 0b01)], 9);
        assert_eq!(bank.next_poly, 9);
        assert_eq!(bank.len(), 3);
    }

    #[test]
    #[ignore = "2026.10.17.8.2 starvation boundary; debug slow but 41.70s green observed (B2 evidence, review 2026-10-20); run under release"]
    fn factors_105bit_by_default_path() {
        let n = crate::parse_decimal("35658606617546370626847644936101").expect("valid input");
        let options = FactorOptions {
            threads: Some(4),
            ..FactorOptions::default()
        };
        let cancelled = AtomicBool::new(false);
        let deadline = Some(Instant::now() + std::time::Duration::from_secs(110));
        let mut progress = |_: &ProgressEvent| Ok(());
        let (divisor, _) = siqs_divisor(&n, &options, &cancelled, deadline, &mut progress, None)
            .expect("105-bit factors by default path");
        assert!(divisor > Integer::from(1) && divisor < n);
        assert_eq!(&n % &divisor, Integer::from(0));
        let mut found = vec![divisor.clone(), &n / &divisor];
        found.sort();
        assert_eq!(
            found,
            vec![
                Integer::from(4373243081080303_u64),
                Integer::from(8153813075658667_u64)
            ]
        );
    }

    #[test]
    fn width_escalation_replays_sweep_grid() {
        let entry = FactorOptions::default().sieve_a_bound;
        let grown: Vec<usize> = SIQS_WIDTH_ESCALATION_STEPS
            .iter()
            .map(|step| entry.saturating_mul(*step))
            .collect();
        assert_eq!(grown, vec![2800, 11200]);
    }

    #[test]
    fn rescue_floor_tracks_measured_score_rank() {
        assert_eq!(siqs_rescue_floor(&[]), 0);
        assert_eq!(siqs_rescue_floor(&[41]), 41);
        assert_eq!(siqs_rescue_floor(&[7, 7, 7, 7]), 7);
        let scores: Vec<u32> = (0..100).collect();
        assert_eq!(siqs_rescue_floor(&scores), 20);
        let mut wide: Vec<u32> = (0..1401).collect();
        wide.reverse();
        assert_eq!(siqs_rescue_floor(&wide), 280);
    }

    #[test]
    fn trial_redundant_flags_always_overlap_only() {
        let always = vec![101_u64, 103, 107];
        assert!(siqs_trial_redundant(101, &always));
        assert!(siqs_trial_redundant(107, &always));
        assert!(!siqs_trial_redundant(2, &always));
        assert!(!siqs_trial_redundant(109, &always));
        assert!(!siqs_trial_redundant(101, &[]));
    }

    #[test]
    fn redundant_trial_skip_keeps_classify_bit_identical() {
        let n = Integer::from(11_783_510_925_019_u64);
        let base = FactorBase::for_siqs(&n, 2000, 2000);
        let words = base.column_count().div_ceil(64);
        let splitters = splitter_primes(&base, 50_000);
        let always = vec![3_u64, 5, 101];
        let with_dup = vec![3_u64, 5, 7, 101, 103];
        let without_dup = vec![7_u64, 103];
        let signed = Integer::from(2_u64)
            * Integer::from(3_u64)
            * Integer::from(5_u64)
            * Integer::from(7_u64)
            * Integer::from(101_u64)
            * Integer::from(103_u64)
            * Integer::from(10_007_u64);
        let mut mask_dup = vec![0_u64; words];
        let dup_out = classify_siqs(
            &signed,
            &base,
            &mut mask_dup,
            50_000,
            &with_dup,
            &always,
            &splitters,
        );
        let mut mask_nodup = vec![0_u64; words];
        let nodup_out = classify_siqs(
            &signed,
            &base,
            &mut mask_nodup,
            50_000,
            &without_dup,
            &always,
            &splitters,
        );
        assert_eq!(mask_dup, mask_nodup);
        let same = std::mem::discriminant(&dup_out) == std::mem::discriminant(&nodup_out);
        assert!(same, "dup={dup_out:?} nodup={nodup_out:?}");
    }

    #[test]
    fn factor_cap_limits_a_factors_without_breaking_congruence() {
        let n = Integer::from(11_783_510_925_019_u64);
        let half_width = 700_i64;
        let interval_bits = (64 - (half_width as u64).leading_zeros()).max(1);
        let target_a = Integer::from(1) << target_a_bits(n.bits(), interval_bits);
        let suitable = suitable_primes(&n, 2000);
        let uncapped =
            generate_polynomial_from(&n, &suitable, &target_a, 0, None).expect("uncapped poly");
        let natural = uncapped.prime_factors.len();
        assert!(natural >= 2);
        let same = generate_polynomial_from(&n, &suitable, &target_a, 0, Some(natural))
            .expect("cap at natural count keeps poly");
        assert_eq!(same.coeff_a, uncapped.coeff_a);
        assert_eq!(same.coeff_b, uncapped.coeff_b);
        if let Some(tight) = generate_polynomial_from(&n, &suitable, &target_a, 0, Some(2)) {
            assert!(tight.prime_factors.len() <= 2);
            assert!((&tight.coeff_b * &tight.coeff_b - &n) % &tight.coeff_a == Integer::from(0));
        }
    }
}
