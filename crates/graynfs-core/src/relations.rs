use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::algebraic::{QuadraticElement, QuadraticRing};
use crate::base::FactorBase;
use crate::error::FactorError;
use crate::factor::ProgressEvent;
use crate::integer::{Integer, IntegerExt, integer_gcd};
use crate::options::FactorOptions;
use crate::polynomial::PolynomialPair;
use crate::prime::is_probable_prime;

const PROGRESS_BATCHES: usize = 8;
const SIEVE_ARRAY_BUDGET_BYTES: u64 = 1 << 30;

pub(crate) fn sieve_array_bytes(blocks: usize, width: usize) -> u64 {
    blocks as u64 * width as u64 * 9
}

#[derive(Clone, Debug)]
pub struct Relation {
    pub a: Integer,
    pub b: Integer,
    pub rational_value: Integer,
    pub algebraic_value: QuadraticElement,
    pub mask: Vec<u64>,
}

#[derive(Clone, Copy)]
struct SieveRoot {
    prime: u64,
    delta: u64,
    weight: u32,
}

enum AlgebraicOutcome {
    Reject,
    Smooth,
    Partial(u64),
}

struct PartialHit {
    relation: Relation,
    large_prime: u64,
}

struct RangeOutput {
    smooth: Vec<Relation>,
    partials: Vec<PartialHit>,
}

impl RangeOutput {
    fn len(&self) -> usize {
        self.smooth.len() + self.partials.len()
    }
}

#[derive(Clone, Copy)]
struct PassContext<'a> {
    polynomial: &'a PolynomialPair,
    base: &'a FactorBase,
    cancelled: &'a AtomicBool,
    deadline: Option<Instant>,
    rational_cache: &'a [Option<Vec<u64>>],
    sieve_roots: &'a [SieveRoot],
    options: &'a FactorOptions,
}

fn pass_workers(rows: i64, options: &FactorOptions) -> usize {
    let available = std::thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    let requested = usize::try_from(rows.max(1)).unwrap_or(1);
    let allowed = options.threads.unwrap_or(available).max(1);
    allowed.min(requested).max(1)
}

pub fn collect_relations<F>(
    _n: &Integer,
    polynomial: &PolynomialPair,
    base: &FactorBase,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    progress: &mut F,
) -> Result<Vec<Relation>, FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    if polynomial.degree != 2 {
        return Err(FactorError::UnsupportedDegree(polynomial.degree));
    }
    let max_a = i64::try_from(options.sieve_a_bound).map_err(|_| {
        FactorError::InvalidInput("sieve a bound exceeds the supported range".to_string())
    })?;
    let max_b = i64::try_from(options.sieve_b_bound).map_err(|_| {
        FactorError::InvalidInput("sieve b bound exceeds the supported range".to_string())
    })?;
    let rational_cache = build_rational_cache(max_a, base);
    let sieve_roots = build_sieve_roots(polynomial, base);
    let initial_threshold = match base.algebraic_modulus_primes.len() {
        0..=200 => 4,
        201..=600 => 8,
        601..=1500 => 12,
        1501..=5000 => 20,
        _ => 28,
    };
    let thresholds = if sieve_roots.is_empty() {
        vec![0]
    } else {
        let mut values = vec![initial_threshold];
        if initial_threshold > 4 {
            values.push(4);
        }
        if initial_threshold > 1 {
            values.push(1);
        }
        values.push(0);
        values.sort_unstable_by(|left, right| right.cmp(left));
        values.dedup();
        values
    };
    let mut best = Vec::new();
    let pairs = (max_a as u64)
        .saturating_mul(2)
        .saturating_add(1)
        .saturating_mul(max_b as u64);
    for threshold in thresholds {
        if threshold == 0 && pairs > 2_000_000 {
            continue;
        }
        let relations = collect_pass(
            max_a,
            max_b,
            threshold,
            PassContext {
                polynomial,
                base,
                cancelled,
                deadline,
                rational_cache: &rational_cache,
                sieve_roots: &sieve_roots,
                options,
            },
            progress,
        )?;
        if relations.len() > best.len() {
            best = relations;
        }
        if best.len() >= base.column_count() {
            return Ok(best);
        }
        if !best.is_empty() {
            break;
        }
    }
    if best.is_empty() {
        return Err(FactorError::InsufficientRelations {
            found: 0,
            required: base.column_count(),
        });
    }
    Ok(best)
}

fn build_rational_cache(max_a: i64, base: &FactorBase) -> Vec<Option<Vec<u64>>> {
    let mut cache = Vec::with_capacity((max_a as usize) * 2 + 1);
    for r in -max_a..=max_a {
        cache.push(rational_mask(r, base));
    }
    cache
}

fn build_sieve_roots(polynomial: &PolynomialPair, base: &FactorBase) -> Vec<SieveRoot> {
    let mut roots = Vec::with_capacity(base.algebraic_primes.len());
    for algebraic in &base.algebraic_primes {
        let prime = algebraic.p;
        let m_mod = (&polynomial.m % Integer::from(prime)).to_u64_saturating();
        let delta = ((u128::from(algebraic.root) + u128::from(prime) - u128::from(m_mod))
            % u128::from(prime)) as u64;
        let weight = (u64::BITS - prime.leading_zeros()).max(1);
        roots.push(SieveRoot {
            prime,
            delta,
            weight,
        });
    }
    roots
}

#[allow(clippy::too_many_arguments)]
fn collect_pass<F>(
    max_a: i64,
    max_b: i64,
    threshold: u32,
    context: PassContext<'_>,
    progress: &mut F,
) -> Result<Vec<Relation>, FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let width = (max_a as usize) * 2 + 1;
    let target = context.base.column_count() as u64;
    let workers = pass_workers(max_b, context.options);
    if workers <= 1 {
        let mut smooth = Vec::new();
        let mut partials = Vec::new();
        let span = (max_b / PROGRESS_BATCHES as i64).max(1);
        let mut start = 1_i64;
        while start <= max_b {
            let end = start.saturating_add(span).min(max_b + 1);
            let output = collect_b_range(start, end, max_a, width, threshold, &context)?;
            smooth.extend(output.smooth);
            partials.extend(output.partials);
            start = end;
            report(progress, smooth.len() + partials.len(), target)?;
        }
        return Ok(finish_pass(&context, smooth, partials));
    }
    let chunk = (max_b + workers as i64 - 1) / workers as i64;
    let ranges: Vec<(i64, i64)> = (0..workers as i64)
        .map(|worker| {
            let start = 1 + worker * chunk;
            let end = start.saturating_add(chunk).min(max_b + 1);
            (start, end)
        })
        .filter(|&(start, end)| start < end)
        .collect();
    let (sender, receiver) = std::sync::mpsc::channel();
    let mut smooth: Vec<Relation> = Vec::new();
    let mut partials: Vec<PartialHit> = Vec::new();
    let mut failure: Option<FactorError> = None;
    let mut reported = 0_usize;
    std::thread::scope(|scope| {
        let handles: Vec<_> = ranges
            .iter()
            .map(|&(start, stop)| {
                let sender = sender.clone();
                scope.spawn(move || {
                    let outcome = collect_b_range(start, stop, max_a, width, threshold, &context);
                    let found = outcome.as_ref().map(RangeOutput::len).unwrap_or(0);
                    let _ = sender.send(found);
                    outcome
                })
            })
            .collect();
        drop(sender);
        while let Ok(found) = receiver.recv() {
            reported += found;
            if failure.is_none()
                && let Err(error) = report(progress, reported, target)
            {
                failure.get_or_insert(error);
            }
        }
        for handle in handles {
            match handle.join() {
                Ok(Ok(output)) => {
                    smooth.extend(output.smooth);
                    partials.extend(output.partials);
                }
                Ok(Err(error)) => {
                    failure.get_or_insert(error);
                }
                Err(_) => {
                    failure.get_or_insert(FactorError::Cancelled);
                }
            };
        }
    });
    if let Some(error) = failure {
        return Err(error);
    }
    Ok(finish_pass(&context, smooth, partials))
}

fn mask_is_empty(mask: &[u64]) -> bool {
    mask.iter().all(|word| *word == 0)
}

fn merge_pair(ring: &QuadraticRing, left: &Relation, right: &Relation) -> Relation {
    let mut mask = left.mask.clone();
    for (target, source) in mask.iter_mut().zip(right.mask.iter()) {
        *target ^= *source;
    }
    Relation {
        a: &left.a * &right.a,
        b: &left.b * &right.b,
        rational_value: &left.rational_value * &right.rational_value,
        algebraic_value: ring.multiply(&left.algebraic_value, &right.algebraic_value),
        mask,
    }
}

fn fold_group(ring: &QuadraticRing, mut group: Vec<Relation>) -> Vec<Relation> {
    let mut folded = Vec::with_capacity(group.len() / 2);
    while let (Some(left), Some(right)) = (group.pop(), group.pop()) {
        folded.push(merge_pair(ring, &left, &right));
    }
    folded
}

fn finish_pass(
    context: &PassContext<'_>,
    mut smooth: Vec<Relation>,
    partials: Vec<PartialHit>,
) -> Vec<Relation> {
    if !partials.is_empty() {
        let ring = QuadraticRing::new(context.polynomial.constant());
        let mut groups: HashMap<u64, Vec<Relation>> = HashMap::new();
        for hit in partials {
            groups
                .entry(hit.large_prime)
                .or_default()
                .push(hit.relation);
        }
        let mut keys: Vec<u64> = groups.keys().copied().collect();
        keys.sort_unstable();
        for key in keys {
            smooth.extend(fold_group(&ring, groups.remove(&key).unwrap_or_default()));
        }
    }
    smooth.retain(|relation| !mask_is_empty(&relation.mask));
    smooth
}

fn report<F>(progress: &mut F, found: usize, target: u64) -> Result<(), FactorError>
where
    F: FnMut(&ProgressEvent) -> Result<(), FactorError>,
{
    let event = ProgressEvent {
        stage: "collect".to_string(),
        current: found as u64,
        total: target,
        message: "sieving the quadratic base".to_string(),
    };
    progress(&event)
}

fn divisor_primes(
    divisors: &[u64],
    epochs: &[u8],
    epoch: u8,
    width: usize,
    offset: usize,
    sieve_roots: &[SieveRoot],
    out: &mut Vec<u64>,
) {
    out.clear();
    for block_index in 0..divisors.len() / width {
        let flat = block_index * width + offset;
        if epochs[flat] != epoch {
            continue;
        }
        let mut word = divisors[flat];
        while word != 0 {
            let bit = word.trailing_zeros() as usize;
            word &= word - 1;
            let root_index = block_index * 64 + bit;
            if let Some(root) = sieve_roots.get(root_index) {
                out.push(root.prime);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn classify_algebraic_known(
    norm: u128,
    a: u128,
    b: u128,
    base: &FactorBase,
    mask: &mut [u64],
    base_bound: u64,
    lp_bound: u64,
    primes: &[u64],
) -> AlgebraicOutcome {
    let words_fit = a <= u128::from(u64::MAX) && b <= u128::from(u64::MAX);
    let mut residual = norm;
    for prime in primes.iter().copied() {
        if residual == 1 {
            break;
        }
        let prime_value = u128::from(prime);
        if residual < prime_value {
            break;
        }
        if residual <= u128::from(u64::MAX) {
            let prime_small = prime;
            let mut small = residual as u64;
            if small < prime_small {
                break;
            }
            let mut exponent = 0_u32;
            while small.is_multiple_of(prime_small) {
                small /= prime_small;
                exponent += 1;
            }
            if exponent == 0 {
                continue;
            }
            residual = u128::from(small);
            if exponent & 1 != 0 {
                if words_fit {
                    base.toggle_algebraic_u64(
                        mask,
                        prime,
                        (a as u64) % prime_small,
                        (b as u64) % prime_small,
                    );
                } else {
                    base.toggle_algebraic_u64(
                        mask,
                        prime,
                        (a % prime_value) as u64,
                        (b % prime_value) as u64,
                    );
                }
            }
            continue;
        }
        let mut exponent = 0_u32;
        while residual.is_multiple_of(prime_value) {
            residual /= prime_value;
            exponent += 1;
        }
        if exponent == 0 {
            continue;
        }
        if exponent & 1 != 0 {
            base.toggle_algebraic_u64(
                mask,
                prime,
                (a % prime_value) as u64,
                (b % prime_value) as u64,
            );
        }
    }
    if residual == 1 {
        AlgebraicOutcome::Smooth
    } else {
        large_prime_part(residual, base_bound, lp_bound)
    }
}

fn collect_b_range(
    b_start: i64,
    b_end: i64,
    max_a: i64,
    width: usize,
    threshold: u32,
    context: &PassContext<'_>,
) -> Result<RangeOutput, FactorError> {
    let PassContext {
        polynomial,
        base,
        cancelled,
        deadline,
        rational_cache,
        sieve_roots,
        options,
    } = *context;
    let base_bound = base.algebraic_modulus_primes.last().copied().unwrap_or(1);
    let lp_bound = u64::try_from(options.large_prime_bound).unwrap_or(u64::MAX);
    let mut smooth = Vec::new();
    let mut partials = Vec::new();
    let root_count = sieve_roots.len();
    let block_count = root_count.div_ceil(64);
    if sieve_array_bytes(block_count, width) > SIEVE_ARRAY_BUDGET_BYTES {
        return Err(FactorError::ResourceLimit(
            "sieve arrays exceed memory budget".to_string(),
        ));
    }
    let mut divisors = vec![0_u64; block_count * width];
    let mut epochs = vec![0_u8; block_count * width];
    let mut scores = vec![0_u32; width];
    let mut candidates = Vec::new();
    let mut epoch = 1_u8;
    for b in b_start..b_end {
        check_cancelled(cancelled, deadline)?;
        epoch = epoch.wrapping_add(1);
        if epoch == 0 {
            divisors.fill(0);
            epochs.fill(0);
            epoch = 1;
        }
        let b_value = Integer::from(b);
        scores.fill(0);
        if root_count > 0 {
            for block_index in 0..block_count {
                let base_index = block_index * width;
                let low = block_index * 64;
                let high = root_count.min(low + 64);
                for (offset, root) in sieve_roots[low..high].iter().enumerate() {
                    if (low + offset) & 255 == 0 {
                        check_cancelled(cancelled, deadline)?;
                    }
                    let residue = ((u128::from(b as u64) * u128::from(root.delta))
                        % u128::from(root.prime)) as u64;
                    let mut board = DivisorBoard {
                        words: &mut divisors[base_index..base_index + width],
                        epochs: &mut epochs[base_index..base_index + width],
                        epoch,
                    };
                    mark_progression(
                        &mut scores,
                        &mut board,
                        1_u64 << offset,
                        max_a as usize,
                        root,
                        residue,
                    );
                }
            }
        }
        for (offset, cached_mask) in rational_cache.iter().enumerate() {
            let r = offset as i64 - max_a;
            if r == 0 {
                continue;
            }
            if scores[offset] < threshold {
                continue;
            }
            let Some(rational_mask) = cached_mask else {
                continue;
            };
            let a_value = &polynomial.m * &b_value + Integer::from(r);
            if a_value <= Integer::from(0) || !is_primitive(&a_value, &b_value) {
                continue;
            }
            let mut mask = rational_mask.clone();
            let outcome = match norm_u128(&a_value, &b_value, polynomial) {
                Some(norm) => {
                    divisor_primes(
                        &divisors,
                        &epochs,
                        epoch,
                        width,
                        offset,
                        sieve_roots,
                        &mut candidates,
                    );
                    classify_algebraic_known(
                        norm,
                        a_value.to_u128_saturating(),
                        b_value.to_u128_saturating(),
                        base,
                        &mut mask,
                        base_bound,
                        lp_bound,
                        &candidates,
                    )
                }
                None => classify_algebraic(
                    &a_value, &b_value, polynomial, base, &mut mask, base_bound, lp_bound,
                ),
            };
            if matches!(outcome, AlgebraicOutcome::Reject) {
                continue;
            }
            let relation = Relation {
                a: a_value.clone(),
                b: b_value.clone(),
                rational_value: Integer::from(r),
                algebraic_value: QuadraticElement {
                    x: a_value,
                    y: -b_value.clone(),
                },
                mask,
            };
            match outcome {
                AlgebraicOutcome::Smooth => smooth.push(relation),
                AlgebraicOutcome::Partial(large_prime) => partials.push(PartialHit {
                    relation,
                    large_prime,
                }),
                AlgebraicOutcome::Reject => continue,
            }
        }
    }
    Ok(RangeOutput { smooth, partials })
}

fn check_cancelled(cancelled: &AtomicBool, deadline: Option<Instant>) -> Result<(), FactorError> {
    if cancelled.load(Ordering::Relaxed) {
        return Err(FactorError::Cancelled);
    }
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return Err(FactorError::ResourceLimit(
            "relation collection deadline exceeded".to_string(),
        ));
    }
    Ok(())
}

fn norm_u128(a: &Integer, b: &Integer, polynomial: &PolynomialPair) -> Option<u128> {
    if a.is_negative() || b.is_negative() {
        return None;
    }
    let a_value = a.to_u128_saturating();
    let b_value = b.to_u128_saturating();
    let square = a_value.checked_mul(a_value)?;
    let b_square = b_value.checked_mul(b_value)?;
    let constant = polynomial.constant().to_u128_saturating();
    let mixed = polynomial
        .coefficients
        .get(1)
        .map(|value| value.to_u128_saturating())
        .unwrap_or(0);
    let product = constant.checked_mul(b_square)?;
    let cross = mixed.checked_mul(a_value.checked_mul(b_value)?)?;
    square.checked_sub(cross)?.checked_add(product)
}

struct DivisorBoard<'a> {
    words: &'a mut [u64],
    epochs: &'a mut [u8],
    epoch: u8,
}

fn mark_progression(
    scores: &mut [u32],
    board: &mut DivisorBoard,
    bit: u64,
    max_a: usize,
    root: &SieveRoot,
    residue: u64,
) {
    let start = ((u128::from(residue) + max_a as u128) % u128::from(root.prime)) as usize;
    let step = match usize::try_from(root.prime) {
        Ok(step) if step > 0 => step,
        _ => scores.len(),
    };
    let mut index = start;
    while index < scores.len() {
        scores[index] = scores[index].saturating_add(root.weight);
        if board.epochs[index] == board.epoch {
            board.words[index] |= bit;
        } else {
            board.words[index] = bit;
            board.epochs[index] = board.epoch;
        }
        if step >= scores.len() - index {
            break;
        }
        index += step;
    }
}

fn is_primitive(a: &Integer, b: &Integer) -> bool {
    if !a.is_negative()
        && !b.is_negative()
        && let (Some(a_value), Some(b_value)) = (a.to_u128(), b.to_u128())
    {
        return gcd_u128(a_value, b_value) == 1;
    }
    integer_gcd(a, b) == Integer::from(1)
}

fn gcd_u128(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn rational_mask(r: i64, base: &FactorBase) -> Option<Vec<u64>> {
    if r == 0 {
        return None;
    }
    let mut value = r.unsigned_abs();
    let mut mask = vec![0_u64; base.column_count().div_ceil(64)];
    if r < 0 {
        mask[base.sign_column / 64] ^= 1_u64 << (base.sign_column % 64);
    }
    for prime in &base.rational_primes {
        if value == 1 {
            break;
        }
        if value < *prime {
            break;
        }
        while value.is_multiple_of(*prime) {
            let index = base.rational_indices[prime];
            mask[index / 64] ^= 1_u64 << (index % 64);
            value /= *prime;
        }
    }
    (value == 1).then_some(mask)
}

fn classify_algebraic(
    a: &Integer,
    b: &Integer,
    polynomial: &PolynomialPair,
    base: &FactorBase,
    mask: &mut [u64],
    base_bound: u64,
    lp_bound: u64,
) -> AlgebraicOutcome {
    let constant = polynomial.constant();
    if !a.is_negative() && !b.is_negative() {
        let a_value = a.to_u128_saturating();
        let b_value = b.to_u128_saturating();
        let c_value = constant.to_u128_saturating();
        if let Some(norm) = a_value.checked_mul(a_value).and_then(|value| {
            value.checked_add(c_value.checked_mul(b_value.checked_mul(b_value)?)?)
        }) {
            return classify_algebraic_u128(
                norm, a_value, b_value, base, mask, base_bound, lp_bound,
            );
        }
    }
    let norm = a * a + constant * (b * b);
    classify_algebraic_bigint(norm, a, b, base, mask, base_bound, lp_bound)
}

fn classify_algebraic_u128(
    mut norm: u128,
    a: u128,
    b: u128,
    base: &FactorBase,
    mask: &mut [u64],
    base_bound: u64,
    lp_bound: u64,
) -> AlgebraicOutcome {
    for prime in &base.algebraic_modulus_primes {
        if norm == 1 {
            break;
        }
        let prime_value = *prime as u128;
        if norm < prime_value {
            break;
        }
        let mut exponent = 0_u32;
        while norm.is_multiple_of(prime_value) {
            norm /= prime_value;
            exponent += 1;
        }
        if exponent == 0 {
            continue;
        }
        if !base.has_algebraic_root(*prime) {
            if exponent & 1 != 0 {
                return AlgebraicOutcome::Reject;
            }
            continue;
        }
        if exponent & 1 != 0 {
            base.toggle_algebraic_u64(
                mask,
                *prime,
                (a % prime_value) as u64,
                (b % prime_value) as u64,
            );
        }
    }
    if norm == 1 {
        return AlgebraicOutcome::Smooth;
    }
    large_prime_part(norm, base_bound, lp_bound)
}

fn classify_algebraic_bigint(
    mut norm: Integer,
    a: &Integer,
    b: &Integer,
    base: &FactorBase,
    mask: &mut [u64],
    base_bound: u64,
    lp_bound: u64,
) -> AlgebraicOutcome {
    for prime in &base.algebraic_modulus_primes {
        if norm.is_one() {
            break;
        }
        let prime_value = Integer::from(*prime);
        if norm < prime_value {
            break;
        }
        let mut exponent = 0_u32;
        while (&norm % &prime_value).is_zero() {
            norm /= &prime_value;
            exponent += 1;
        }
        if exponent == 0 {
            continue;
        }
        if !base.has_algebraic_root(*prime) {
            if exponent & 1 != 0 {
                return AlgebraicOutcome::Reject;
            }
            continue;
        }
        if exponent & 1 != 0 {
            base.toggle_algebraic_bigint(mask, *prime, a, b);
        }
    }
    if norm.is_one() {
        return AlgebraicOutcome::Smooth;
    }
    match norm.to_u128() {
        Some(value) => large_prime_part(value, base_bound, lp_bound),
        None => AlgebraicOutcome::Reject,
    }
}

fn large_prime_part(residual: u128, base_bound: u64, lp_bound: u64) -> AlgebraicOutcome {
    let Ok(small) = u64::try_from(residual) else {
        return AlgebraicOutcome::Reject;
    };
    if small <= base_bound || small > lp_bound {
        return AlgebraicOutcome::Reject;
    }
    if is_probable_prime(&Integer::from(small)) {
        AlgebraicOutcome::Partial(small)
    } else {
        AlgebraicOutcome::Reject
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sieve_budget_rejects_terabyte_layout() {
        assert!(sieve_array_bytes(339, 1_000_001) > SIEVE_ARRAY_BUDGET_BYTES);
        assert!(sieve_array_bytes(120, 32_001) < SIEVE_ARRAY_BUDGET_BYTES);
        assert_eq!(sieve_array_bytes(0, 1_000_001), 0);
    }

    fn relation_with_mask(word: u64) -> Relation {
        Relation {
            a: Integer::from(1_u64),
            b: Integer::from(1_u64),
            rational_value: Integer::from(1_u64),
            algebraic_value: QuadraticElement {
                x: Integer::from(1_u64),
                y: -Integer::from(1_u64),
            },
            mask: vec![word],
        }
    }

    #[test]
    fn fold_group_pairs_only() {
        let ring = QuadraticRing::new(Integer::from(2_u64));
        for (size, expected) in [(0_usize, 0_usize), (1, 0), (2, 1), (3, 1), (4, 2), (5, 2)] {
            let group: Vec<Relation> = (0..size)
                .map(|index| relation_with_mask(1_u64 << index))
                .collect();
            let folded = fold_group(&ring, group);
            assert_eq!(folded.len(), expected, "group size {size}");
        }
    }
}
