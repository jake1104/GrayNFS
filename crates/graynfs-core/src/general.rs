use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use num_traits::{Euclid, One, ToPrimitive, Zero};

use crate::error::FactorError;
use crate::factor::{StageTiming, Timer, TimingCollector};
use crate::integer::{Integer, IntegerExt, integer_abs, integer_gcd, integer_sqrt};
use crate::lattice::{ideal_basis, shortest_vector};
use crate::linear_algebra::{Masked, rank_and_dependencies};
use crate::number_field::{Element, NumberField, polynomial_discriminant};
use crate::options::{
    FactorOptions, estimate_collect_bytes, estimate_live_bytes, estimate_total_bytes,
};
use crate::poly::ModPolynomial;
use crate::poly::is_irreducible;
use crate::poly::mod_inverse;
use crate::polynomial::PolynomialPair;
use crate::prime::{is_prime_u64_deterministic, mul_mod, sieve_primes};

const UNIT_POWER_LIMIT: usize = 8;
const UNIT_TARGET_LIMIT: usize = 8;
const UNIT_SEARCH_LIMIT: usize = 4;

fn residue_mod(value: &Integer, prime: u64) -> u64 {
    value.rem_euclid(&Integer::from(prime)).to_u64_saturating()
}

#[derive(Clone, Debug)]
struct GeneralRelation {
    rational_value: Integer,
    algebraic_value: Element,
    mask: Vec<u64>,
    large_prime: Option<(u64, u64)>,
    double_edge: Option<[(u64, u64); 2]>,
}

impl Masked for GeneralRelation {
    fn mask(&self) -> &[u64] {
        &self.mask
    }
}

struct ClassifiedNorm {
    mask: Vec<u64>,
    large_primes: Vec<(u64, u64)>,
}

struct LargePrimeCollector {
    pending: HashMap<(u64, u64), usize>,
}

impl LargePrimeCollector {
    fn new() -> Self {
        Self {
            pending: HashMap::new(),
        }
    }

    fn pair(&mut self, key: (u64, u64), index: usize) -> Option<usize> {
        match self.pending.remove(&key) {
            Some(previous) => Some(previous),
            None => {
                self.pending.insert(key, index);
                None
            }
        }
    }
}

struct PartialBank {
    singles: Vec<GeneralRelation>,
    single_sources: Vec<(i64, i64)>,
    single_aux: Vec<u64>,
    single_lps: Vec<Vec<(u64, u64)>>,
    collector: LargePrimeCollector,
    rat_pending: HashMap<u64, usize>,
    mixed_pending: HashMap<((u64, u64), u64), usize>,
    triad_pending: HashMap<(Vec<(u64, u64)>, u64), usize>,
    doubles: Vec<GeneralRelation>,
    double_sources: Vec<(i64, i64)>,
    double_aux: Vec<u64>,
    double_parent: HashMap<(u64, u64), (u64, u64)>,
    double_members: HashMap<(u64, u64), Vec<usize>>,
}

impl PartialBank {
    fn new() -> Self {
        Self {
            singles: Vec::new(),
            single_sources: Vec::new(),
            single_aux: Vec::new(),
            single_lps: Vec::new(),
            collector: LargePrimeCollector::new(),
            rat_pending: HashMap::new(),
            mixed_pending: HashMap::new(),
            triad_pending: HashMap::new(),
            doubles: Vec::new(),
            double_sources: Vec::new(),
            double_aux: Vec::new(),
            double_parent: HashMap::new(),
            double_members: HashMap::new(),
        }
    }

    fn find_vertex(&mut self, vertex: (u64, u64)) -> (u64, u64) {
        let mut root = vertex;
        while let Some(next) = self.double_parent.get(&root) {
            if *next == root {
                break;
            }
            root = *next;
        }
        let mut current = vertex;
        while current != root {
            let next = self.double_parent[&current];
            self.double_parent.insert(current, root);
            current = next;
        }
        root
    }

    fn add_double_edge(&mut self, edge: [(u64, u64); 2], index: usize) -> Option<Vec<usize>> {
        let left = self.find_vertex(edge[0]);
        let right = self.find_vertex(edge[1]);
        if left != right {
            self.double_parent.insert(right, left);
            let mut members = self.double_members.remove(&left).unwrap_or_default();
            members.extend(self.double_members.remove(&right).unwrap_or_default());
            members.push(index);
            self.double_members.insert(left, members);
            return None;
        }
        let members = self.double_members.remove(&left).unwrap_or_default();
        let Some(path) = double_cycle_path(&self.doubles, &members, edge[0], edge[1]) else {
            let mut all = members;
            all.push(index);
            self.double_members.insert(left, all);
            return None;
        };
        for member in &path {
            if let Some(endpoints) = self.doubles[*member].double_edge {
                self.double_parent.remove(&endpoints[0]);
                self.double_parent.remove(&endpoints[1]);
            }
        }
        let used: HashSet<usize> = path.iter().copied().collect();
        let leftover: Vec<usize> = members
            .into_iter()
            .filter(|member| !used.contains(member))
            .collect();
        if let Some(&first) = leftover.first()
            && let Some(endpoints) = self.doubles[first].double_edge
        {
            self.double_parent.insert(endpoints[0], endpoints[0]);
            self.double_parent.insert(endpoints[1], endpoints[1]);
            let mut group = vec![first];
            for member in leftover.iter().skip(1) {
                if let Some(ends) = self.doubles[*member].double_edge {
                    let left_root = self.find_vertex(ends[0]);
                    let right_root = self.find_vertex(ends[1]);
                    if left_root != right_root {
                        self.double_parent.insert(right_root, left_root);
                    }
                }
                group.push(*member);
            }
            let root = self.find_vertex(endpoints[0]);
            self.double_members.insert(root, group);
        }
        let mut cycle = path;
        cycle.push(index);
        Some(cycle)
    }
}

fn double_cycle_path(
    doubles: &[GeneralRelation],
    members: &[usize],
    from: (u64, u64),
    to: (u64, u64),
) -> Option<Vec<usize>> {
    if from == to {
        return Some(Vec::new());
    }
    let mut adjacency: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
    for member in members {
        if let Some(endpoints) = doubles[*member].double_edge {
            adjacency.entry(endpoints[0]).or_default().push(*member);
            adjacency.entry(endpoints[1]).or_default().push(*member);
        }
    }
    let mut queue: VecDeque<(u64, u64)> = VecDeque::new();
    let mut came_from: HashMap<(u64, u64), ((u64, u64), usize)> = HashMap::new();
    queue.push_back(from);
    while let Some(vertex) = queue.pop_front() {
        if vertex == to {
            break;
        }
        let Some(edges) = adjacency.get(&vertex) else {
            continue;
        };
        for member in edges {
            let endpoints = doubles[*member].double_edge.unwrap_or_default();
            let next = if endpoints[0] == vertex {
                endpoints[1]
            } else {
                endpoints[0]
            };
            if next == from || came_from.contains_key(&next) {
                continue;
            }
            came_from.insert(next, (vertex, *member));
            queue.push_back(next);
        }
    }
    if !came_from.contains_key(&to) {
        return None;
    }
    let mut path = Vec::new();
    let mut current = to;
    while current != from {
        let (previous, member) = came_from.get(&current).copied()?;
        path.push(member);
        current = previous;
    }
    path.reverse();
    Some(path)
}

#[derive(Clone, Debug)]
struct GeneralBase {
    rational_primes: Vec<u64>,
    rational_indices: HashMap<u64, usize>,
    algebraic_roots: Vec<(u64, u64, usize)>,
    algebraic_index: HashMap<u64, Vec<(u64, usize)>>,
    modulus_primes: Vec<u64>,
    discriminant: Integer,
    algebraic_bound: u64,
    large_prime_bound: u64,
    large_prime_roots: HashMap<u64, Vec<u64>>,
    sign_column: usize,
    lp_columns: usize,
}

impl GeneralBase {
    fn build(
        polynomial: &PolynomialPair,
        n: &Integer,
        rational_bound: usize,
        algebraic_bound: usize,
        large_prime_bound: usize,
    ) -> Result<Self, FactorError> {
        let mut rational_primes = Vec::new();
        let mut rational_indices = HashMap::new();
        for prime in sieve_primes(rational_bound) {
            if n % Integer::from(prime) == Integer::from(0) {
                continue;
            }
            rational_indices.insert(prime, rational_primes.len());
            rational_primes.push(prime);
        }
        let mut algebraic_roots = Vec::new();
        let mut modulus_primes = Vec::new();
        let discriminant = discriminant(polynomial);
        for prime in sieve_primes(algebraic_bound) {
            if n % Integer::from(prime) == Integer::from(0)
                || discriminant.clone() % Integer::from(prime) == Integer::from(0)
            {
                continue;
            }
            modulus_primes.push(prime);
            let coefficients = polynomial
                .coefficients
                .iter()
                .map(|value| residue_mod(value, prime))
                .collect();
            let roots = ModPolynomial::new(coefficients, prime).roots();
            for root in roots {
                algebraic_roots.push((prime, root, 0));
            }
        }
        for (index, entry) in algebraic_roots.iter_mut().enumerate() {
            entry.2 = rational_primes.len() + index;
        }
        let mut algebraic_index: HashMap<u64, Vec<(u64, usize)>> = HashMap::new();
        for (prime, root, column) in algebraic_roots.iter().copied() {
            algebraic_index
                .entry(prime)
                .or_default()
                .push((root, column));
        }
        let sign_column = rational_primes.len() + algebraic_roots.len();
        Ok(Self {
            rational_primes,
            rational_indices,
            algebraic_roots,
            algebraic_index,
            modulus_primes,
            discriminant,
            algebraic_bound: bounded_u64(algebraic_bound),
            large_prime_bound: bounded_u64(large_prime_bound),
            large_prime_roots: HashMap::new(),
            sign_column,
            lp_columns: 0,
        })
    }

    fn column_count(&self) -> usize {
        self.sign_column + 1 + self.lp_columns
    }

    fn mask_words(&self) -> usize {
        self.column_count().div_ceil(64)
    }

    fn sieve_roots(&self, m: &Integer) -> Vec<(u64, u64, u32)> {
        self.algebraic_roots
            .iter()
            .map(|(prime, root, _)| {
                let m_mod = residue_mod(m, *prime);
                let delta = (root + prime - m_mod) % prime;
                let weight = (u64::BITS - prime.leading_zeros()).max(1);
                (*prime, delta, weight)
            })
            .collect()
    }
}

fn bounded_u64(value: usize) -> u64 {
    u64::try_from(value).unwrap_or(u64::MAX)
}

fn scan_failure_taxonomy(units_enabled: bool, found: bool) -> &'static str {
    if found {
        "found"
    } else if units_enabled {
        "unit_correction_failed"
    } else {
        "sqrt_failed"
    }
}

fn dependency_workers(dep_count: usize, options: &FactorOptions) -> usize {
    let threads = options_thread_count(options);
    threads.min(dep_count).max(1)
}

fn options_thread_count(options: &FactorOptions) -> usize {
    let available = std::thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    options.threads.unwrap_or(available).max(1)
}

#[allow(clippy::too_many_arguments)]
fn scan_dependencies(
    n: &Integer,
    field: &NumberField,
    evaluation_point: &Integer,
    relations: &[GeneralRelation],
    dependencies: &[Vec<usize>],
    units: &[Element],
    base: &GeneralBase,
    deadline: Option<Instant>,
    units_enabled: bool,
    limit: usize,
    collector: &TimingCollector,
    degree: usize,
    options: &FactorOptions,
) -> Option<Integer> {
    let attempt_count = limit.min(dependencies.len());
    if attempt_count == 0 {
        return None;
    }
    let smooth_words = (base.sign_column + 1).div_ceil(64);
    let mut head: Vec<usize> = Vec::with_capacity(attempt_count);
    let mut tail: Vec<usize> = Vec::new();
    for (index, dependency) in dependencies.iter().enumerate().take(attempt_count) {
        let anchored = dependency.iter().any(|member| {
            relations.get(*member).is_some_and(|relation| {
                relation
                    .mask
                    .iter()
                    .take(smooth_words)
                    .any(|word| *word != 0)
            })
        });
        if anchored {
            head.push(index);
        } else {
            tail.push(index);
        }
    }
    head.extend(tail);
    let order = head;
    let workers = dependency_workers(attempt_count, options);
    let scan_started = Instant::now();
    let units_tried = units.len();
    if workers <= 1 {
        for index in &order {
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                collector.record_elapsed(
                    "square_root_scan",
                    Some(degree),
                    scan_started.elapsed().as_millis(),
                    format!(
                        "attempts={} workers={} units={} units_tried={} slowest_ms=0 found=false outcome={}",
                        attempt_count,
                        workers,
                        units_enabled,
                        units_tried,
                        scan_failure_taxonomy(units_enabled, false)
                    ),
                );
                return None;
            }
            let attempt = extract_factor(
                n,
                field,
                evaluation_point,
                relations,
                &dependencies[*index],
                units,
                base,
                deadline,
                None,
                units_enabled,
            );
            if let Some(divisor) = attempt {
                collector.record_elapsed(
                    "square_root_scan",
                    Some(degree),
                    scan_started.elapsed().as_millis(),
                    format!(
                        "attempts={} workers={} units={} units_tried={} slowest_ms=0 found=true outcome={}",
                        attempt_count,
                        workers,
                        units_enabled,
                        units_tried,
                        scan_failure_taxonomy(units_enabled, true)
                    ),
                );
                return Some(divisor);
            }
        }
        collector.record_elapsed(
            "square_root_scan",
            Some(degree),
            scan_started.elapsed().as_millis(),
            format!(
                "attempts={} workers={} units={} units_tried={} slowest_ms=0 found=false outcome={}",
                attempt_count,
                workers,
                units_enabled,
                units_tried,
                scan_failure_taxonomy(units_enabled, false)
            ),
        );
        return None;
    }
    let found = Mutex::new(None);
    let stop = AtomicBool::new(false);
    let attempted: Mutex<Vec<(usize, usize)>> = Mutex::new(Vec::new());
    let chunk = attempt_count.div_ceil(workers);
    let queue: Vec<Vec<usize>> = (0..workers)
        .map(|worker| {
            let start = (worker * chunk).min(attempt_count);
            let end = start.saturating_add(chunk).min(attempt_count);
            (start..end).collect()
        })
        .filter(|value: &Vec<usize>| !value.is_empty())
        .collect();
    let found = &found;
    let stop = &stop;
    let attempted = &attempted;
    let order = &order;
    std::thread::scope(|scope| {
        for slice in &queue {
            let slice: &[usize] = slice;
            scope.spawn(move || {
                for index in slice {
                    if stop.load(Ordering::Relaxed) {
                        break;
                    }
                    if deadline.is_some_and(|limit| Instant::now() >= limit) {
                        break;
                    }
                    let true_index = order[*index];
                    let dependency = &dependencies[true_index];
                    let started = Instant::now();
                    let attempt = extract_factor(
                        n,
                        field,
                        evaluation_point,
                        relations,
                        dependency,
                        units,
                        base,
                        deadline,
                        Some(stop),
                        units_enabled,
                    );
                    let elapsed_ms = started.elapsed().as_millis();
                    if let Ok(mut entries) = attempted.lock() {
                        entries.push((true_index, elapsed_ms as usize));
                    }
                    if let Some(divisor) = attempt
                        && let Ok(mut slot) = found.lock()
                        && slot.is_none()
                    {
                        *slot = Some((true_index, divisor));
                        stop.store(true, Ordering::Relaxed);
                    }
                }
            });
        }
    });
    let mut entries = match attempted.lock() {
        Ok(value) => value.clone(),
        Err(_) => Vec::new(),
    };
    entries.sort_by_key(|entry| entry.0);
    let found_index = found
        .lock()
        .map(|value| value.as_ref().map(|entry| entry.0))
        .unwrap_or(None);
    let total_ms: u128 = entries.iter().map(|entry| entry.1 as u128).sum();
    let slowest = entries.iter().map(|entry| entry.1).max().unwrap_or(0);
    collector.record_elapsed(
        "square_root_scan",
        Some(degree),
        total_ms / u128::try_from(workers).unwrap_or(1).max(1),
        format!(
            "attempts={} workers={} units={} units_tried={} slowest_ms={} found={} outcome={}",
            entries.len(),
            workers,
            units_enabled,
            units.len(),
            slowest,
            found_index.is_some(),
            scan_failure_taxonomy(units_enabled, found_index.is_some())
        ),
    );
    for (index, elapsed_ms) in entries
        .iter()
        .filter(|(index, _)| *index < 8 || Some(*index) == found_index)
    {
        collector.record_elapsed(
            "square_root",
            Some(degree),
            *elapsed_ms as u128,
            format!(
                "dependency={index} units={units_enabled} threads={workers} found={}",
                found_index == Some(*index)
            ),
        );
    }
    match found.lock() {
        Ok(slot) => slot.as_ref().map(|entry| entry.1.clone()),
        Err(_) => None,
    }
}

pub fn gnfs_divisor(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
) -> Result<(Integer, PolynomialPair, usize, usize, Vec<StageTiming>), FactorError> {
    let collector = TimingCollector::default();
    match gnfs_divisor_timed(n, options, cancelled, deadline, &collector) {
        Ok((divisor, polynomial, relations, rank)) => {
            Ok((divisor, polynomial, relations, rank, collector.snapshot()))
        }
        Err(error) => Err(error),
    }
}

fn gnfs_divisor_timed(
    n: &Integer,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    collector: &TimingCollector,
) -> Result<(Integer, PolynomialPair, usize, usize), FactorError> {
    let total = Instant::now();
    let degree = options.degree.unwrap_or(3);
    if degree < 3 {
        return Err(FactorError::UnsupportedDegree(degree));
    }
    let polynomial = {
        let timer = Timer::new(Some(collector), "select_polynomial", Some(degree));
        let selection = select_polynomial_with_stats(n, degree, options.polynomial_candidates)?;
        let detail = format!(
            "m={} rejected_small={} rejected_final={}",
            selection.polynomial.m, selection.small_rejections, selection.final_rejections
        );
        timer.finish(detail);
        selection.polynomial
    };
    let field = {
        let timer = Timer::new(Some(collector), "number_field", Some(degree));
        let value = NumberField::new(degree, polynomial.coefficients.clone())
            .ok_or(FactorError::NoPolynomial)?;
        timer.finish(format!("degree={}", value.degree));
        value
    };
    let mut base = {
        let timer = Timer::new(Some(collector), "factor_base", Some(degree));
        let value = GeneralBase::build(
            &polynomial,
            n,
            options.rational_base_bound,
            options.algebraic_base_bound,
            options.large_prime_bound,
        )?;
        timer.finish(format!("columns={}", value.column_count()));
        let columns = value.column_count().max(1);
        let rows = (columns.saturating_mul(3) / 2 + 16).max(16);
        let need_bytes = estimate_total_bytes(columns, rows);
        let budget_bytes = options.memory_budget_bytes();
        let budget_mb = budget_bytes.div_ceil(1024 * 1024);
        let need_mb = need_bytes.div_ceil(1024 * 1024);
        if need_mb as u64 > budget_mb {
            return Err(FactorError::ResourceLimit(format!(
                "the relation matrix would need about {need_mb} MiB, above the memory budget of {budget_mb} MiB"
            )));
        }
        value
    };
    let bound = unit_search_bound(field.degree);
    let mut units: Option<Vec<Element>> = None;
    let (mut relations, mut bank): (Vec<GeneralRelation>, PartialBank) = {
        let timer = Timer::new(Some(collector), "relation_collection", Some(degree));
        let (complete, bank) = collect_relations(
            n,
            &polynomial,
            &field,
            &mut base,
            options,
            cancelled,
            deadline,
            collector,
            degree,
        )?;
        let value: Vec<GeneralRelation> = complete
            .into_iter()
            .filter(|relation| relation.mask.iter().any(|word| *word != 0))
            .collect();
        timer.finish(format!(
            "relations={} columns={}",
            value.len(),
            base.column_count()
        ));
        (value, bank)
    };
    {
        let need_bytes = estimate_total_bytes(base.column_count(), relations.len());
        let budget_bytes = options.memory_budget_bytes();
        if need_bytes > budget_bytes {
            return Err(FactorError::ResourceLimit(format!(
                "the relation matrix would need about {} MiB, above the memory budget of {} MiB",
                need_bytes.div_ceil(1024 * 1024),
                budget_bytes.div_ceil(1024 * 1024)
            )));
        }
        if let Some(error) = live_budget_exceeded(
            options,
            base.column_count(),
            relations.len(),
            bank.singles.len(),
            bank.doubles.len(),
        ) {
            return Err(error);
        }
    }
    let (mut rank, mut dependencies) = {
        let timer = Timer::new(Some(collector), "matrix_rank", Some(degree));
        let value = rank_and_dependencies(
            &relations,
            base.column_count(),
            1024,
            2,
            options_thread_count(options),
        );
        timer.finish(format!("rank={} of {}", value.0, relations.len()));
        value
    };
    if relations.len() <= rank && (!bank.singles.is_empty() || !bank.doubles.is_empty()) {
        if let Some(error) = live_budget_exceeded(
            options,
            base.column_count(),
            relations.len(),
            bank.singles.len(),
            bank.doubles.len(),
        ) {
            return Err(error);
        }
        let drain_started = Instant::now();
        let pre_rows = relations.len();
        let pre_columns = base.column_count();
        drain_partial_bank(&mut base, &mut bank, &mut relations);
        let filtered: Vec<GeneralRelation> = relations
            .into_iter()
            .filter(|relation| relation.mask.iter().any(|word| *word != 0))
            .collect();
        relations = filtered;
        let rerank = rank_and_dependencies(
            &relations,
            base.column_count(),
            1024,
            2,
            options_thread_count(options),
        );
        rank = rerank.0;
        dependencies = rerank.1;
        let detail = format!(
            "rows={}+{} columns={}+{} rank={}",
            pre_rows,
            relations.len().saturating_sub(pre_rows),
            pre_columns,
            base.column_count().saturating_sub(pre_columns),
            rank
        );
        collector.record_elapsed(
            "drain",
            Some(degree),
            drain_started.elapsed().as_millis(),
            detail,
        );
        let need_bytes = estimate_total_bytes(base.column_count(), relations.len());
        let budget_bytes = options.memory_budget_bytes();
        if need_bytes > budget_bytes {
            return Err(FactorError::ResourceLimit(format!(
                "the relation matrix would need about {} MiB, above the memory budget of {} MiB",
                need_bytes.div_ceil(1024 * 1024),
                budget_bytes.div_ceil(1024 * 1024)
            )));
        }
        if let Some(error) = live_budget_exceeded(
            options,
            base.column_count(),
            relations.len(),
            bank.singles.len(),
            bank.doubles.len(),
        ) {
            return Err(error);
        }
    }
    if relations.len() <= rank {
        collector.record_elapsed(
            "matrix_rank",
            Some(degree),
            0,
            format!(
                "rank_deficient relations={} rank={} dependencies={}",
                relations.len(),
                rank,
                dependencies.len()
            ),
        );
        return Err(FactorError::InsufficientRelations {
            found: relations.len(),
            required: rank + 1,
        });
    }
    {
        let timer = Timer::new(Some(collector), "dependency_search", Some(degree));
        timer.finish(format!("dependencies={}", dependencies.len()));
    }
    for (units_enabled, limit) in [(false, 1024usize), (true, 256usize)] {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            break;
        }
        if units_enabled && units.is_none() {
            let timer = Timer::new(Some(collector), "unit_search", Some(degree));
            let value = find_units(&field, bound, UNIT_SEARCH_LIMIT);
            timer.finish(format!("bound={} found={}", bound, value.len()));
            units = Some(value);
        }
        let scan_units: &[Element] = units.as_deref().unwrap_or(&[]);
        let scan = scan_dependencies(
            n,
            &field,
            &polynomial.m,
            &relations,
            &dependencies,
            scan_units,
            &base,
            deadline,
            units_enabled,
            limit,
            collector,
            degree,
            options,
        );
        if let Some(divisor) = scan {
            if units.is_none() {
                collector.record_elapsed(
                    "unit_search",
                    Some(degree),
                    0,
                    "deferred: first scan won".to_string(),
                );
            }
            collector.record(
                "total",
                Some(degree),
                total,
                format!(
                    "relations={} rank={} dependencies={}",
                    relations.len(),
                    rank,
                    dependencies.len()
                ),
            );
            return Ok((divisor, polynomial, relations.len(), rank));
        }
    }
    collector.record(
        "total",
        Some(degree),
        total,
        format!(
            "failed sqrt_failed unit_correction_failed relations={} rank={} dependencies={} units_tried={} bound={}",
            relations.len(),
            rank,
            dependencies.len(),
            units.as_deref().map_or(0, |value| value.len()),
            bound
        ),
    );
    Err(FactorError::AlgebraicSquareRootFailed)
}

const SMALL_INDEX_PRIMES: [u64; 10] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];

const IRREDUCIBILITY_PRIMES: [u64; 10] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];

const IRREDUCIBILITY_RETRY_PRIMES: [u64; 2] = [1009, 1013];

enum IrreducibilityOutcome {
    AcceptedSmall,
    AcceptedRetry,
    Rejected,
}

fn irreducible_mod_any(candidate: &PolynomialPair, primes: &[u64]) -> bool {
    for prime in primes {
        let residues: Vec<u64> = candidate
            .coefficients
            .iter()
            .map(|value| residue_mod(value, *prime))
            .collect();
        let reduced = ModPolynomial::new(residues, *prime);
        if reduced.degree() != candidate.degree {
            continue;
        }
        if is_irreducible(&reduced, *prime) {
            return true;
        }
    }
    false
}

fn classify_irreducibility(candidate: &PolynomialPair) -> IrreducibilityOutcome {
    if irreducible_mod_any(candidate, &IRREDUCIBILITY_PRIMES) {
        IrreducibilityOutcome::AcceptedSmall
    } else if irreducible_mod_any(candidate, &IRREDUCIBILITY_RETRY_PRIMES) {
        IrreducibilityOutcome::AcceptedRetry
    } else {
        IrreducibilityOutcome::Rejected
    }
}

#[cfg(test)]
fn is_irreducible_candidate(candidate: &PolynomialPair) -> bool {
    !matches!(
        classify_irreducibility(candidate),
        IrreducibilityOutcome::Rejected
    )
}

fn has_unit_index(polynomial: &PolynomialPair) -> bool {
    let value = integer_abs(&discriminant(polynomial));
    if value.is_zero() {
        return false;
    }
    SMALL_INDEX_PRIMES
        .iter()
        .all(|prime| value.clone() % Integer::from(*prime) != Integer::from(0))
}

fn polynomial_sample_score(candidate: &PolynomialPair) -> f64 {
    let degree = candidate.degree;
    let mut coefficients = Vec::with_capacity(degree + 1);
    for coefficient in candidate.coefficients.iter() {
        let Some(value) = coefficient.to_f64() else {
            return f64::INFINITY;
        };
        coefficients.push(value);
    }
    let m_value = candidate.m.to_f64().unwrap_or(f64::INFINITY);
    if !m_value.is_finite() {
        return f64::INFINITY;
    }
    let mut total = 0.0;
    let mut count = 0_usize;
    for b in 1..=8_i64 {
        let b_value = b as f64;
        let scale = b_value.powi(degree as i32);
        for r in [1_i64, -1, 2, -2, 3, -3, 5, -5, 7, -7, 11, -11] {
            let x = (m_value * b_value + r as f64) / b_value;
            let mut value = 0.0;
            for coefficient in coefficients.iter().rev() {
                value = value * x + coefficient;
            }
            let norm = (value * scale).abs();
            if norm.is_finite() && norm > 0.0 {
                total += norm.log2();
                count += 1;
            }
        }
    }
    if count == 0 {
        return f64::INFINITY;
    }
    total / count as f64 + polynomial_root_penalty(candidate)
}

fn polynomial_root_penalty(candidate: &PolynomialPair) -> f64 {
    let mut penalty = 0.0;
    for prime in sieve_primes(128) {
        let prime_value = prime as f64;
        let residues: Vec<u64> = candidate
            .coefficients
            .iter()
            .map(|value| residue_mod(value, prime))
            .collect();
        let roots = ModPolynomial::new(residues, prime).roots().len() as f64;
        penalty += prime_value.log2() / (prime_value - 1.0) * (1.0 - roots);
    }
    penalty
}

#[cfg(test)]
pub fn inspect_polynomial(n: &Integer, degree: usize) -> Result<PolynomialPair, FactorError> {
    select_polynomial(n, degree, FactorOptions::default().polynomial_candidates)
}

struct PolynomialSelection {
    polynomial: PolynomialPair,
    small_rejections: usize,
    final_rejections: usize,
}

#[cfg(test)]
fn select_polynomial(
    n: &Integer,
    degree: usize,
    candidates: usize,
) -> Result<PolynomialPair, FactorError> {
    select_polynomial_with_stats(n, degree, candidates).map(|selection| selection.polynomial)
}

fn select_polynomial_with_stats(
    n: &Integer,
    degree: usize,
    candidates: usize,
) -> Result<PolynomialSelection, FactorError> {
    let mut m = integer_root(n, degree);
    let mut best: Option<(PolynomialPair, f64)> = None;
    let mut small_rejections = 0_usize;
    let mut final_rejections = 0_usize;
    for _ in 0..candidates.max(1) {
        let power = integer_pow(&m, degree);
        if m >= Integer::from(1) && n > &power {
            let constant = n - &power;
            if !is_perfect_power(&constant, degree) {
                if degree <= 2 {
                    let mut coefficients = vec![Integer::from(0); degree + 1];
                    coefficients[0] = constant;
                    coefficients[degree] = Integer::from(1);
                    return Ok(PolynomialSelection {
                        polynomial: PolynomialPair {
                            degree,
                            coefficients,
                            m,
                        },
                        small_rejections,
                        final_rejections,
                    });
                }
                for linear in [1_i64, 2, 4, 5, 7, 8, 10, 11] {
                    for quadratic in [1_i64, 2, 3, 5, 7, 11] {
                        let mut coefficients = vec![Integer::from(0); degree + 1];
                        coefficients[degree] = Integer::from(1);
                        let mut offset = Integer::from(0);
                        for (index, coefficient) in
                            coefficients.iter_mut().enumerate().take(degree).skip(1)
                        {
                            *coefficient =
                                Integer::from(if index == 1 { linear } else { quadratic });
                            let value = coefficient.clone();
                            offset += value * integer_pow(&m, index);
                        }
                        let candidate_constant = n - &power - &offset;
                        if candidate_constant <= Integer::from(0)
                            || is_perfect_power(&candidate_constant, degree)
                        {
                            continue;
                        }
                        coefficients[0] = candidate_constant;
                        let candidate = PolynomialPair {
                            degree,
                            coefficients,
                            m: m.clone(),
                        };
                        match classify_irreducibility(&candidate) {
                            IrreducibilityOutcome::Rejected => {
                                small_rejections += 1;
                                final_rejections += 1;
                                continue;
                            }
                            IrreducibilityOutcome::AcceptedRetry => {
                                small_rejections += 1;
                            }
                            IrreducibilityOutcome::AcceptedSmall => {}
                        }
                        if has_unit_index(&candidate) {
                            let score = polynomial_sample_score(&candidate);
                            let replace = match &best {
                                None => true,
                                Some((_, best_score)) => score < *best_score,
                            };
                            if replace {
                                best = Some((candidate, score));
                            }
                        }
                    }
                }
            }
        }
        m -= Integer::from(1);
    }
    best.map(|(candidate, _)| PolynomialSelection {
        polynomial: candidate,
        small_rejections,
        final_rejections,
    })
    .ok_or(FactorError::NoPolynomial)
}

struct ScanContext<'a> {
    n: &'a Integer,
    polynomial: &'a PolynomialPair,
    field: &'a NumberField,
    base: &'a GeneralBase,
    lp_cache: &'a Mutex<HashMap<u64, Vec<u64>>>,
    max_a: i64,
    max_b: i64,
    allow_double: bool,
}

fn check_limits(cancelled: &AtomicBool, deadline: Option<Instant>) -> Result<(), FactorError> {
    if cancelled.load(Ordering::Relaxed) {
        return Err(FactorError::Cancelled);
    }
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return Err(FactorError::ResourceLimit(
            "general relation deadline exceeded".to_string(),
        ));
    }
    Ok(())
}

pub(crate) fn live_budget_exceeded(
    options: &FactorOptions,
    columns: usize,
    relations: usize,
    singles: usize,
    doubles: usize,
) -> Option<FactorError> {
    let need = estimate_live_bytes(columns, relations, singles, doubles);
    let budget = options.memory_budget_bytes();
    if need > budget {
        Some(FactorError::ResourceLimit(format!(
            "relation collection exceeded the memory budget of {} MiB (live estimate {} MiB)",
            budget.div_ceil(1024 * 1024),
            need.div_ceil(1024 * 1024)
        )))
    } else {
        None
    }
}

#[allow(clippy::too_many_arguments)]
fn collect_relations(
    n: &Integer,
    polynomial: &PolynomialPair,
    field: &NumberField,
    base: &mut GeneralBase,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    collector: &TimingCollector,
    degree: usize,
) -> Result<(Vec<GeneralRelation>, PartialBank), FactorError> {
    let max_a = i64::try_from(options.sieve_a_bound)
        .map_err(|_| FactorError::InvalidInput("sieve a bound is too large".to_string()))?;
    let max_b = i64::try_from(options.sieve_b_bound)
        .map_err(|_| FactorError::InvalidInput("sieve b bound is too large".to_string()))?;
    let workers = options_thread_count(options);
    if estimate_collect_bytes(base.column_count(), max_a, workers) > options.memory_budget_bytes() {
        return Err(FactorError::ResourceLimit(format!(
            "general sieve arrays would exceed the memory budget of {} MiB",
            options.memory_budget_bytes().div_ceil(1024 * 1024)
        )));
    }
    let sieve_roots = base.sieve_roots(&polynomial.m);
    let lp_cache = Mutex::new(std::mem::take(&mut base.large_prime_roots));
    let mut context = ScanContext {
        n,
        polynomial,
        field,
        base,
        lp_cache: &lp_cache,
        max_a,
        max_b,
        allow_double: options.enable_double_large_primes,
    };
    let mut consumed: HashSet<(i64, i64)> = HashSet::new();
    let mut bank = PartialBank::new();
    let mut complete: Vec<GeneralRelation> = Vec::new();
    let mut attempt_a = max_a;
    let mut attempt_b = max_b;
    let mut usable_total = 0_usize;
    for _ in 0..5 {
        let pass_started = Instant::now();
        let rational_scores = build_rational_scores(context.base, attempt_a);
        let rational_masks = RationalMasks::build(context.base, attempt_a);
        context.max_a = attempt_a;
        context.max_b = attempt_b;
        complete.extend(scan_pass(
            &mut context,
            &sieve_roots,
            &rational_scores,
            &rational_masks,
            options,
            cancelled,
            deadline,
            true,
            usable_total,
            &mut consumed,
            &mut bank,
        )?);
        usable_total = complete
            .iter()
            .filter(|relation| relation.mask.iter().any(|word| *word != 0))
            .count();
        if let Some(error) = live_budget_exceeded(
            options,
            context.base.column_count(),
            complete.len(),
            bank.singles.len(),
            bank.doubles.len(),
        ) {
            return Err(error);
        }
        collector.record_elapsed(
            "collect",
            Some(degree),
            pass_started.elapsed().as_millis(),
            format!(
                "sieve_pass={} columns={} target={} usable={}",
                complete.len(),
                context.base.column_count(),
                target_relations(context.base.column_count()),
                usable_total
            ),
        );
        if usable_total >= target_relations(context.base.column_count()) {
            break;
        }
        let next_a = (attempt_a as usize)
            .saturating_mul(2)
            .min(COLLECT_SIEVE_A_CAP) as i64;
        let next_b = (attempt_b as usize)
            .saturating_mul(2)
            .min(COLLECT_SIEVE_B_CAP) as i64;
        if next_a == attempt_a && next_b == attempt_b {
            break;
        }
        attempt_a = next_a;
        attempt_b = next_b;
    }
    if complete.len() < 64 {
        let rational_scores = build_rational_scores(context.base, context.max_a);
        let rational_masks = RationalMasks::build(context.base, context.max_a);
        complete.extend(scan_pass(
            &mut context,
            &sieve_roots,
            &rational_scores,
            &rational_masks,
            options,
            cancelled,
            deadline,
            false,
            usable_total,
            &mut consumed,
            &mut bank,
        )?);
    }
    base.large_prime_roots = lp_cache.into_inner().unwrap_or_default();
    let words = base.mask_words();
    for relation in &mut complete {
        grow_mask(&mut relation.mask, words);
    }
    if complete.is_empty() {
        return Err(FactorError::InsufficientRelations {
            found: 0,
            required: base.column_count(),
        });
    }
    Ok((complete, bank))
}

fn drain_partial_bank(
    base: &mut GeneralBase,
    bank: &mut PartialBank,
    complete: &mut Vec<GeneralRelation>,
) {
    if bank.singles.is_empty() && bank.doubles.is_empty() {
        return;
    }
    let mut keys: Vec<(u64, u64)> = Vec::new();
    let mut seen: HashSet<(u64, u64)> = HashSet::new();
    for lps in bank.single_lps.iter() {
        for key in lps.iter() {
            if seen.insert(*key) {
                keys.push(*key);
            }
        }
    }
    for relation in bank.doubles.iter() {
        if let Some(edge) = relation.double_edge {
            for key in edge {
                if seen.insert(key) {
                    keys.push(key);
                }
            }
        }
    }
    for aux in bank
        .single_aux
        .iter()
        .chain(bank.double_aux.iter())
        .copied()
    {
        if aux != 1 {
            let key = (aux, u64::MAX);
            if seen.insert(key) {
                keys.push(key);
            }
        }
    }
    if keys.is_empty() {
        base.lp_columns = 0;
        return;
    }
    let base_columns = base.column_count();
    base.lp_columns = keys.len();
    let mut index: HashMap<(u64, u64), usize> = HashMap::with_capacity(keys.len());
    for (position, key) in keys.iter().enumerate() {
        index.insert(*key, base_columns + position);
    }
    let words = (base_columns + keys.len()).div_ceil(64);
    for relation in complete.iter_mut() {
        grow_mask(&mut relation.mask, words);
    }
    let singles = std::mem::take(&mut bank.singles);
    let lps_list = std::mem::take(&mut bank.single_lps);
    let aux_list = std::mem::take(&mut bank.single_aux);
    for ((mut relation, lps), aux) in singles.into_iter().zip(lps_list).zip(aux_list) {
        grow_mask(&mut relation.mask, words);
        for key in lps.iter() {
            if let Some(column) = index.get(key) {
                relation.mask[column / 64] |= 1_u64 << (column % 64);
            }
        }
        if aux != 1
            && let Some(column) = index.get(&(aux, u64::MAX))
        {
            relation.mask[column / 64] |= 1_u64 << (column % 64);
        }
        relation.large_prime = None;
        complete.push(relation);
    }
    let doubles = std::mem::take(&mut bank.doubles);
    let double_aux = std::mem::take(&mut bank.double_aux);
    for (mut relation, aux) in doubles.into_iter().zip(double_aux) {
        let edge = relation.double_edge.take();
        grow_mask(&mut relation.mask, words);
        if let Some(edge) = edge {
            for key in edge {
                if let Some(column) = index.get(&key) {
                    relation.mask[column / 64] |= 1_u64 << (column % 64);
                }
            }
        }
        if aux != 1
            && let Some(column) = index.get(&(aux, u64::MAX))
        {
            relation.mask[column / 64] |= 1_u64 << (column % 64);
        }
        complete.push(relation);
    }
}

fn target_relations(columns: usize) -> usize {
    columns + columns / 2 + 16
}

fn build_rational_scores(base: &GeneralBase, max_a: i64) -> Vec<u32> {
    let max_a_usize = max_a.max(0) as usize;
    let width = max_a_usize * 2 + 1;
    let mut scores = vec![0_u32; width];
    for prime in &base.rational_primes {
        let weight = log2_weight(*prime);
        let mut power = *prime;
        while power <= width as u64 {
            mark_general_residue(&mut scores, max_a_usize, power, 0, weight);
            match power.checked_mul(*prime) {
                Some(next) => power = next,
                None => break,
            }
        }
    }
    scores
}

fn log2_weight(prime: u64) -> u32 {
    (u64::BITS - prime.leading_zeros()).max(1)
}

const COLLECT_WAVE_B: i64 = 32;
const COLLECT_SIEVE_A_CAP: usize = 500_000;
const COLLECT_SIEVE_B_CAP: usize = 60_000;

struct RawHit {
    relation: GeneralRelation,
    large_primes: Vec<(u64, u64)>,
    source: (i64, i64),
    rat_lp: u64,
}

struct BatchParams<'a> {
    n: &'a Integer,
    polynomial: &'a PolynomialPair,
    field: &'a NumberField,
    base: &'a GeneralBase,
    lp_cache: &'a Mutex<HashMap<u64, Vec<u64>>>,
    sieve_roots: &'a [(u64, u64, u32)],
    rational_scores: &'a [u32],
    rational_masks: &'a RationalMasks,
    rational_thresholds: &'a [u32],
    trial_primes: &'a [u64],
    allow_double: bool,
    degree: u32,
    constant_bits: f64,
    constant_magnitude: f64,
    m_magnitude: f64,
    large_prime_allowance: f64,
    max_a: i64,
    width: usize,
    use_sieve: bool,
}

fn scan_single_b(
    b: i64,
    params: &BatchParams<'_>,
    consumed: &HashSet<(i64, i64)>,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    scratch: &mut Vec<u32>,
) -> Result<Vec<RawHit>, FactorError> {
    check_limits(cancelled, deadline)?;
    let mut hits = Vec::new();
    let b_unsigned = b.unsigned_abs();
    let m_offset = params.polynomial.m.clone() * Integer::from(b);
    let m_small = m_offset.to_i64().unwrap_or(i64::MAX);
    let lower_index = (params.max_a + 1 - m_small).max(0) as usize;
    let scores: Option<&[u32]> = if params.use_sieve {
        scratch.clear();
        scratch.resize(params.width, 0);
        for (prime, delta, weight) in params.sieve_roots {
            let residue = ((b as u128 * *delta as u128) % *prime as u128) as u64;
            mark_general_residue(scratch, params.max_a as usize, *prime, residue, *weight);
        }
        Some(scratch.as_slice())
    } else {
        None
    };
    let b_float = b as f64;
    let b_magnitude = b_float.max(1.0);
    let b_degree_power = b_magnitude.powf(params.degree as f64);
    let m_times_b = params.m_magnitude * b_float;
    let allowance = params.large_prime_allowance * if params.allow_double { 2.0 } else { 1.0 };
    let degree_float = params.degree as f64;
    let estimate_limit = allowance + 8_f64;
    let estimate_residual = params.constant_magnitude * b_degree_power;
    let fast_reject = m_times_b.is_finite() && estimate_residual.is_finite() && params.degree <= 6;
    let a_low = m_times_b - params.max_a as f64;
    let a_high = m_times_b + params.max_a as f64;
    let constant_exponent = if fast_reject && a_low >= 1.0 && a_high <= 1e50 {
        let slack = 8.0 * f64::EPSILON;
        let low = (a_low * (1.0 - slack)).abs().max(1.0).to_bits() >> 52 & 0x7ff;
        let high = (a_high * (1.0 + slack)).abs().max(1.0).to_bits() >> 52 & 0x7ff;
        (low == high).then_some(low as f64)
    } else {
        None
    };
    let constant_bound =
        constant_exponent.map(|exponent| degree_float * (exponent - 1023.0) - 1e-6);
    for index in lower_index..params.width {
        if index % 64 == 0 {
            check_limits(cancelled, deadline)?;
        }
        let r = index as i64 - params.max_a;
        if r == 0 {
            continue;
        }
        if scores.as_ref().is_some_and(|scores| scores[index] == 0) {
            continue;
        }
        if params.use_sieve {
            let score = scores.as_ref().map(|scores| scores[index]).unwrap_or(0);
            if let Some(bound) = constant_bound
                && (score as f64) + estimate_limit < bound
            {
                continue;
            }
            let a_magnitude = (m_times_b + r as f64).abs().max(1.0);
            let a_bits = a_magnitude.to_bits();
            if constant_bound.is_none()
                && fast_reject
                && a_magnitude <= 1e50
                && (score as f64) + estimate_limit
                    < degree_float * (((a_bits >> 52) & 0x7ff) as f64 - 1023.0) - 1e-6
            {
                continue;
            }
            let norm_estimate =
                a_magnitude.powf(params.degree as f64) + params.constant_magnitude * b_degree_power;
            let norm_bits = if norm_estimate.is_finite() && norm_estimate > 0.0 {
                norm_estimate.log2()
            } else {
                (params.degree as f64) * a_magnitude.max(b_magnitude).log2()
                    + params.constant_bits
                    + 1.0
            };
            if (score as f64) + allowance + 8_f64 < norm_bits {
                continue;
            }
        }
        if params.use_sieve
            && r.abs() != 1
            && params
                .rational_scores
                .get(index)
                .is_some_and(|score| *score == 0)
        {
            continue;
        }
        if params.use_sieve
            && params
                .rational_scores
                .get(index)
                .is_some_and(|score| *score < params.rational_thresholds[index])
        {
            continue;
        }
        if consumed.contains(&(b, r)) {
            continue;
        }
        if gcd_u64(r.unsigned_abs(), b_unsigned) != 1 {
            continue;
        }
        let Some((mask, rat_lp)) = params.rational_masks.mask(r) else {
            continue;
        };
        let b_integer = Integer::from(b);
        let a_integer = &m_offset + Integer::from(r);
        let norm = params.field.homogeneous(&a_integer, &b_integer);
        let Some(classified) = classify_norm(
            norm,
            &a_integer,
            &b_integer,
            params.n,
            params.polynomial,
            params.base,
            params.lp_cache,
            mask,
            params.trial_primes,
        ) else {
            continue;
        };
        let ClassifiedNorm { mask, large_primes } = classified;
        let mut algebraic = vec![Integer::from(0); params.field.degree];
        algebraic[0] = a_integer.clone();
        algebraic[1] = -b_integer.clone();
        let relation = GeneralRelation {
            rational_value: Integer::from(r),
            algebraic_value: Element::from_coefficients(algebraic),
            mask,
            large_prime: None,
            double_edge: None,
        };
        hits.push(RawHit {
            relation,
            large_primes,
            source: (b, r),
            rat_lp,
        });
    }
    Ok(hits)
}

fn cycle_aux_square(bank: &PartialBank, cycle: &[usize]) -> bool {
    let Some(first) = cycle
        .first()
        .and_then(|member| bank.double_aux.get(*member))
        .copied()
    else {
        return false;
    };
    if !cycle
        .iter()
        .all(|member| bank.double_aux.get(*member) == Some(&first))
    {
        return false;
    }
    first == 1 || cycle.len().is_multiple_of(2)
}

#[allow(clippy::too_many_arguments)]
fn merge_hit(
    field: &NumberField,
    words: usize,
    allow_double: bool,
    bank: &mut PartialBank,
    consumed: &mut HashSet<(i64, i64)>,
    complete: &mut Vec<GeneralRelation>,
    hit: RawHit,
) {
    let RawHit {
        mut relation,
        large_primes,
        rat_lp,
        source,
    } = hit;
    if large_primes.len() == 1 && rat_lp == 1 {
        let key = large_primes[0];
        relation.large_prime = Some(key);
        let current = bank.singles.len();
        if let Some(previous) = bank.collector.pair(key, current) {
            consumed.insert(bank.single_sources[previous]);
            consumed.insert(source);
            complete.push(combine_relations(
                field,
                words,
                &bank.singles[previous],
                &relation,
            ));
            return;
        }
        bank.single_aux.push(1);
        bank.single_lps.push(large_primes.clone());
        bank.singles.push(relation);
        bank.single_sources.push(source);
        return;
    }
    if large_primes.len() == 1 {
        let key = large_primes[0];
        relation.large_prime = Some(key);
        let current = bank.singles.len();
        if let Some(previous) = bank.mixed_pending.remove(&(key, rat_lp)) {
            consumed.insert(bank.single_sources[previous]);
            consumed.insert(source);
            complete.push(combine_relations(
                field,
                words,
                &bank.singles[previous],
                &relation,
            ));
            return;
        }
        bank.single_aux.push(rat_lp);
        bank.single_lps.push(large_primes.clone());
        bank.singles.push(relation);
        bank.single_sources.push(source);
        bank.mixed_pending.insert((key, rat_lp), current);
        return;
    }
    if large_primes.is_empty() && rat_lp != 1 {
        let current = bank.singles.len();
        if let Some(previous) = bank.rat_pending.remove(&rat_lp) {
            consumed.insert(bank.single_sources[previous]);
            consumed.insert(source);
            complete.push(combine_relations(
                field,
                words,
                &bank.singles[previous],
                &relation,
            ));
            return;
        }
        bank.single_aux.push(rat_lp);
        bank.single_lps.push(large_primes.clone());
        bank.singles.push(relation);
        bank.single_sources.push(source);
        bank.rat_pending.insert(rat_lp, current);
        return;
    }
    if large_primes.len() == 2 {
        if !allow_double {
            return;
        }
        let edge = [large_primes[0], large_primes[1]];
        relation.double_edge = Some(edge);
        let current = bank.doubles.len();
        bank.doubles.push(relation);
        bank.double_aux.push(rat_lp);
        bank.double_sources.push(source);
        if let Some(cycle) = bank.add_double_edge(edge, current)
            && cycle_aux_square(bank, &cycle)
        {
            let mut combined = bank.doubles[cycle[0]].clone();
            for member in cycle.iter().skip(1) {
                combined = combine_relations(field, words, &combined, &bank.doubles[*member]);
                consumed.insert(bank.double_sources[*member]);
            }
            consumed.insert(bank.double_sources[cycle[0]]);
            consumed.insert(source);
            complete.push(combined);
        }
        return;
    }
    if large_primes.len() == 3 {
        if !allow_double {
            return;
        }
        let mut key = large_primes.clone();
        key.sort_unstable();
        let lookup = (key, rat_lp);
        let current = bank.singles.len();
        if let Some(previous) = bank.triad_pending.remove(&lookup) {
            consumed.insert(bank.single_sources[previous]);
            consumed.insert(source);
            complete.push(combine_relations(
                field,
                words,
                &bank.singles[previous],
                &relation,
            ));
            return;
        }
        bank.single_aux.push(rat_lp);
        bank.single_lps.push(large_primes.clone());
        bank.singles.push(relation);
        bank.single_sources.push(source);
        bank.triad_pending.insert(lookup, current);
        return;
    }
    if !large_primes.is_empty() {
        return;
    }
    consumed.insert(source);
    complete.push(relation);
}

#[allow(clippy::too_many_arguments)]
fn scan_pass(
    context: &mut ScanContext<'_>,
    sieve_roots: &[(u64, u64, u32)],
    rational_scores: &[u32],
    rational_masks: &RationalMasks,
    options: &FactorOptions,
    cancelled: &AtomicBool,
    deadline: Option<Instant>,
    use_sieve: bool,
    already_usable: usize,
    consumed: &mut HashSet<(i64, i64)>,
    bank: &mut PartialBank,
) -> Result<Vec<GeneralRelation>, FactorError> {
    let mut complete: Vec<GeneralRelation> = Vec::new();
    let mut usable = 0_usize;
    let max_a = context.max_a;
    let max_b = context.max_b;
    let width = (max_a as usize) * 2 + 1;
    let degree = context.field.degree as u32;
    let constant_bits = integer_abs(&context.polynomial.coefficients[0]).bits() as f64;
    let constant_magnitude = integer_abs(&context.polynomial.coefficients[0])
        .to_f64()
        .unwrap_or(0.0);
    let m_magnitude = context.polynomial.m.to_f64().unwrap_or(0.0);
    let large_prime_allowance =
        (u64::BITS - context.base.large_prime_bound.max(1).leading_zeros()) as f64;
    let trial_limit = (context.base.large_prime_bound.min(1_000_000) as usize).max(2);
    let trial_primes = sieve_primes(trial_limit);
    let rational_thresholds: Vec<u32> = (0..width)
        .map(|index| {
            let r = index as i64 - max_a;
            (r.unsigned_abs().max(1) as f64).log2().ceil() as u32
        })
        .collect();
    let target = target_relations(context.base.column_count());
    let params = BatchParams {
        n: context.n,
        polynomial: context.polynomial,
        field: context.field,
        base: context.base,
        lp_cache: context.lp_cache,
        sieve_roots,
        rational_scores,
        rational_masks,
        rational_thresholds: &rational_thresholds,
        trial_primes: &trial_primes,
        allow_double: context.allow_double,
        degree,
        constant_bits,
        constant_magnitude,
        m_magnitude,
        large_prime_allowance,
        max_a,
        width,
        use_sieve,
    };
    let workers = dependency_workers(max_b.max(1) as usize, options);
    let mut start = 1_i64;
    while start <= max_b {
        check_limits(cancelled, deadline)?;
        if already_usable + usable >= target {
            break;
        }
        let stride = workers as i64 * COLLECT_WAVE_B;
        let end = (start + stride - 1).min(max_b);
        let length = end - start + 1;
        let chunk = (length - 1) / workers as i64 + 1;
        let mut ranges: Vec<(i64, i64)> = Vec::new();
        let mut cursor = start;
        while cursor <= end {
            let stop = (cursor + chunk - 1).min(end);
            ranges.push((cursor, stop));
            cursor = stop + 1;
        }
        let consumed_shared: &HashSet<(i64, i64)> = consumed;
        let params = &params;
        let batches = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (lower, upper) in ranges {
                handles.push(scope.spawn(move || -> Result<Vec<RawHit>, FactorError> {
                    let mut hits: Vec<RawHit> = Vec::new();
                    let mut scores_scratch: Vec<u32> = Vec::new();
                    let mut b = lower;
                    while b <= upper {
                        check_limits(cancelled, deadline)?;
                        hits.extend(scan_single_b(
                            b,
                            params,
                            consumed_shared,
                            cancelled,
                            deadline,
                            &mut scores_scratch,
                        )?);
                        b += 1;
                    }
                    Ok(hits)
                }));
            }
            handles
                .into_iter()
                .map(|handle| handle.join().unwrap_or(Err(FactorError::Cancelled)))
                .collect::<Vec<Result<Vec<RawHit>, FactorError>>>()
        });
        let mut finished = false;
        let mut since_check = 0_usize;
        for batch in batches {
            for hit in batch? {
                let before = complete.len();
                merge_hit(
                    context.field,
                    context.base.mask_words(),
                    context.allow_double,
                    bank,
                    consumed,
                    &mut complete,
                    hit,
                );
                if complete.len() > before
                    && complete
                        .last()
                        .is_some_and(|relation| relation.mask.iter().any(|word| *word != 0))
                {
                    usable += 1;
                }
                since_check += 1;
                if since_check >= 256 {
                    since_check = 0;
                    if let Some(error) = live_budget_exceeded(
                        options,
                        context.base.column_count(),
                        already_usable + complete.len(),
                        bank.singles.len(),
                        bank.doubles.len(),
                    ) {
                        return Err(error);
                    }
                }
                if already_usable + usable >= target {
                    finished = true;
                    break;
                }
            }
            if finished {
                break;
            }
        }
        if finished {
            break;
        }
        start = end + 1;
    }
    Ok(complete)
}

fn combine_relations(
    field: &NumberField,
    words: usize,
    left: &GeneralRelation,
    right: &GeneralRelation,
) -> GeneralRelation {
    let mut left_mask = left.mask.clone();
    let mut right_mask = right.mask.clone();
    grow_mask(&mut left_mask, words);
    grow_mask(&mut right_mask, words);
    let mask = left_mask
        .iter()
        .zip(right_mask.iter())
        .map(|(left_word, right_word)| left_word ^ right_word)
        .collect();
    GeneralRelation {
        rational_value: &left.rational_value * &right.rational_value,
        algebraic_value: field.multiply(&left.algebraic_value, &right.algebraic_value),
        mask,
        large_prime: None,
        double_edge: None,
    }
}

fn grow_mask(mask: &mut Vec<u64>, words: usize) {
    if mask.len() < words {
        mask.resize(words, 0);
    }
}

fn gcd_u64(left: u64, mut right: u64) -> u64 {
    let mut left = left;
    while right != 0 {
        let next = left % right;
        left = right;
        right = next;
    }
    left
}

type RationalMaskEntry = Option<(Vec<u64>, u64)>;
type RationalMaskShard = Mutex<HashMap<i64, RationalMaskEntry>>;
type RationalMaskCache = Vec<RationalMaskShard>;

struct RationalMasks {
    max_a: i64,
    words: usize,
    sign_column: usize,
    rational_primes: Vec<u64>,
    rational_indices: HashMap<u64, usize>,
    large_prime_bound: u64,
    cache: RationalMaskCache,
}

impl RationalMasks {
    fn build(base: &GeneralBase, max_a: i64) -> Self {
        let mut cache = Vec::with_capacity(64);
        for _ in 0..64 {
            cache.push(Mutex::new(HashMap::new()));
        }
        Self {
            max_a,
            words: base.mask_words(),
            sign_column: base.sign_column,
            rational_primes: base.rational_primes.clone(),
            rational_indices: base.rational_indices.clone(),
            large_prime_bound: base.large_prime_bound,
            cache,
        }
    }

    fn shard(&self, value: i64) -> &RationalMaskShard {
        &self.cache[(value as u64 % self.cache.len() as u64) as usize]
    }

    fn compute(&self, value: i64) -> Option<(Vec<u64>, u64)> {
        let mut remaining = value.unsigned_abs();
        let mut mask = vec![0_u64; self.words];
        if value < 0 {
            mask[self.sign_column / 64] ^= 1_u64 << (self.sign_column % 64);
        }
        for prime in &self.rational_primes {
            if remaining == 1 {
                break;
            }
            while remaining.is_multiple_of(*prime) {
                let index = self.rational_indices[prime];
                mask[index / 64] ^= 1_u64 << (index % 64);
                remaining /= *prime;
            }
        }
        if remaining == 1 {
            return Some((mask, 1));
        }
        if remaining <= self.large_prime_bound && is_prime_u64_deterministic(remaining) {
            return Some((mask, remaining));
        }
        None
    }

    fn mask(&self, value: i64) -> Option<(Vec<u64>, u64)> {
        if value == 0 || value < -self.max_a || value > self.max_a {
            return None;
        }
        if let Some(hit) = self
            .shard(value)
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .get(&value)
            .cloned()
        {
            return hit;
        }
        let result = self.compute(value);
        self.shard(value)
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .insert(value, result.clone());
        result
    }
}

#[allow(clippy::too_many_arguments)]
fn classify_norm_small(
    start: u128,
    a: &Integer,
    b: &Integer,
    n: &Integer,
    polynomial: &PolynomialPair,
    base: &GeneralBase,
    lp_cache: &Mutex<HashMap<u64, Vec<u64>>>,
    mut mask: Vec<u64>,
    trial_primes: &[u64],
) -> Option<ClassifiedNorm> {
    let mut remaining = start;
    for prime in &base.modulus_primes {
        let prime_value = u128::from(*prime);
        if remaining == 1 {
            break;
        }
        if remaining < prime_value {
            break;
        }
        let mut exponent = 0_u32;
        while remaining.is_multiple_of(prime_value) {
            remaining /= prime_value;
            exponent += 1;
        }
        if exponent == 0 {
            continue;
        }
        let Some(columns) = base.algebraic_index.get(prime) else {
            if exponent & 1 != 0 {
                return None;
            }
            continue;
        };
        if exponent & 1 != 0 {
            let a_mod = residue_mod(a, *prime);
            let b_mod = residue_mod(b, *prime);
            for (root, column) in columns.iter().copied() {
                if a_mod == (b_mod * root) % prime {
                    mask[column / 64] ^= 1_u64 << (column % 64);
                }
            }
        }
    }
    if remaining == 1 {
        return Some(ClassifiedNorm {
            mask,
            large_primes: Vec::new(),
        });
    }
    classify_large_prime(
        Integer::from(remaining),
        a,
        b,
        n,
        polynomial,
        base,
        lp_cache,
        mask,
        trial_primes,
    )
}

#[allow(clippy::too_many_arguments)]
fn classify_norm(
    mut norm: Integer,
    a: &Integer,
    b: &Integer,
    n: &Integer,
    polynomial: &PolynomialPair,
    base: &GeneralBase,
    lp_cache: &Mutex<HashMap<u64, Vec<u64>>>,
    mut mask: Vec<u64>,
    trial_primes: &[u64],
) -> Option<ClassifiedNorm> {
    if norm >= Integer::from(0)
        && let Some(start) = norm.to_u128()
    {
        return classify_norm_small(
            start,
            a,
            b,
            n,
            polynomial,
            base,
            lp_cache,
            mask,
            trial_primes,
        );
    }
    for prime in &base.modulus_primes {
        let prime_integer = Integer::from(*prime);
        if norm.is_one() {
            break;
        }
        if norm < prime_integer {
            break;
        }
        let mut exponent = 0_u32;
        while (&norm % &prime_integer).is_zero() {
            norm /= &prime_integer;
            exponent += 1;
        }
        if exponent == 0 {
            continue;
        }
        let Some(columns) = base.algebraic_index.get(prime) else {
            if exponent & 1 != 0 {
                return None;
            }
            continue;
        };
        if exponent & 1 != 0 {
            let a_mod = residue_mod(a, *prime);
            let b_mod = residue_mod(b, *prime);
            for (root, column) in columns.iter().copied() {
                if a_mod == (b_mod * root) % prime {
                    mask[column / 64] ^= 1_u64 << (column % 64);
                }
            }
        }
    }
    if norm.is_one() {
        return Some(ClassifiedNorm {
            mask,
            large_primes: Vec::new(),
        });
    }
    classify_large_prime(
        norm,
        a,
        b,
        n,
        polynomial,
        base,
        lp_cache,
        mask,
        trial_primes,
    )
}

#[allow(clippy::too_many_arguments)]
fn match_lp_root(
    prime: u64,
    a: &Integer,
    b: &Integer,
    n: &Integer,
    polynomial: &PolynomialPair,
    base: &GeneralBase,
    lp_cache: &Mutex<HashMap<u64, Vec<u64>>>,
) -> Option<u64> {
    let prime_integer = Integer::from(prime);
    if n % &prime_integer == Integer::from(0)
        || &base.discriminant % &prime_integer == Integer::from(0)
    {
        return None;
    }
    let a_mod = residue_mod(a, prime);
    let b_mod = residue_mod(b, prime);
    if b_mod != 0 {
        let root = mul_mod(a_mod, mod_inverse(b_mod, prime), prime);
        let coefficients = polynomial
            .coefficients
            .iter()
            .map(|value| residue_mod(value, prime))
            .collect();
        if ModPolynomial::new(coefficients, prime).eval(root) == 0 {
            return Some(root);
        }
        return None;
    }
    let roots = lp_cache
        .lock()
        .unwrap_or_else(|error| error.into_inner())
        .get(&prime)
        .cloned()
        .unwrap_or_else(|| {
            let coefficients = polynomial
                .coefficients
                .iter()
                .map(|value| residue_mod(value, prime))
                .collect();
            let roots = ModPolynomial::new(coefficients, prime).roots();
            lp_cache
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .insert(prime, roots.clone());
            roots
        });
    roots
        .iter()
        .copied()
        .find(|root| a_mod == mul_mod(b_mod, *root, prime))
}

#[allow(clippy::too_many_arguments)]
fn classify_large_prime(
    norm: Integer,
    a: &Integer,
    b: &Integer,
    n: &Integer,
    polynomial: &PolynomialPair,
    base: &GeneralBase,
    lp_cache: &Mutex<HashMap<u64, Vec<u64>>>,
    mask: Vec<u64>,
    trial_primes: &[u64],
) -> Option<ClassifiedNorm> {
    let bound = base.large_prime_bound;
    if norm.bits() > 64 {
        return None;
    }
    let cofactor = norm.to_u64_saturating();
    let span = bound as u128;
    if (cofactor as u128) > span * span * span {
        return None;
    }
    let mut remaining = cofactor;
    let mut kernel: Vec<u64> = Vec::new();
    if is_prime_u64_deterministic(cofactor) {
        if cofactor > base.algebraic_bound
            && let Some(root) = match_lp_root(cofactor, a, b, n, polynomial, base, lp_cache)
        {
            return Some(ClassifiedNorm {
                mask,
                large_primes: vec![(cofactor, root)],
            });
        }
        return None;
    }
    for prime in trial_primes {
        let prime = *prime;
        if (prime as u128) * (prime as u128) > remaining as u128 {
            break;
        }
        if !remaining.is_multiple_of(prime) {
            continue;
        }
        let mut exponent = 0_u32;
        while remaining.is_multiple_of(prime) {
            remaining /= prime;
            exponent += 1;
        }
        if exponent & 1 != 0 {
            kernel.push(prime);
        }
    }
    if remaining > 1 {
        if !is_prime_u64_deterministic(remaining) {
            return None;
        }
        kernel.push(remaining);
    }
    if kernel.len() > 3 {
        return None;
    }
    let mut large_primes = Vec::with_capacity(kernel.len());
    for prime in kernel {
        if prime <= base.algebraic_bound {
            return None;
        }
        let root = match_lp_root(prime, a, b, n, polynomial, base, lp_cache)?;
        large_primes.push((prime, root));
    }
    Some(ClassifiedNorm { mask, large_primes })
}

const MAX_SQUARE_ROOT_BITS: u64 = 131_072;
const ATTEMPT_BASE_MS: u64 = 1_500;
const ATTEMPT_SCALE_SHIFT: u32 = 13;
const ATTEMPT_CEILING_MS: u64 = 12_000;
const SUBSET_GENERATOR_LIMIT: usize = 10;
const SUBSET_TRIAL_LIMIT: usize = 4;

fn attempt_budget_ms(alpha_bits: u64) -> u64 {
    let scaled = (alpha_bits.saturating_mul(alpha_bits)) >> ATTEMPT_SCALE_SHIFT;
    ATTEMPT_BASE_MS
        .saturating_add(scaled)
        .min(ATTEMPT_CEILING_MS)
}

fn attempt_deadline(
    deadline: Option<Instant>,
    started: Instant,
    alpha_bits: u64,
) -> Option<Instant> {
    let budget = started + std::time::Duration::from_millis(attempt_budget_ms(alpha_bits));
    Some(match deadline {
        Some(limit) if limit < budget => limit,
        _ => budget,
    })
}

fn direct_congruence_of_squares(
    n: &Integer,
    evaluation_point: &Integer,
    relations: &[GeneralRelation],
    dependency: &[usize],
    deadline: Option<Instant>,
) -> Option<Integer> {
    if dependency.is_empty() {
        return None;
    }
    if deadline.is_some_and(|limit| Instant::now() >= limit) {
        return None;
    }
    let mut rational = Integer::from(1);
    let mut evaluated = Integer::from(1);
    for index in dependency {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return None;
        }
        let relation = relations.get(*index)?;
        rational *= &relation.rational_value;
        if relation.algebraic_value.coefficients.len() < 2 {
            return None;
        }
        let constant = &relation.algebraic_value.coefficients[0];
        let slope = &relation.algebraic_value.coefficients[1];
        evaluated *= &(constant - slope * evaluation_point);
    }
    let rational = crate::integer::integer_abs(&rational);
    let rational_root = integer_sqrt(&rational);
    if &rational_root * &rational_root != rational {
        return None;
    }
    let evaluated = crate::integer::integer_abs(&evaluated);
    let evaluated_root = integer_sqrt(&evaluated);
    if &evaluated_root * &evaluated_root != evaluated {
        return None;
    }
    for x_sign in [Integer::from(1), Integer::from(-1)] {
        let x = (&evaluated_root * &x_sign) % n;
        for y_sign in [Integer::from(1), Integer::from(-1)] {
            let y = (&rational_root * &y_sign) % n;
            if (&x * &x - &y * &y) % n != Integer::from(0) {
                continue;
            }
            let divisor = integer_gcd(&(&x - &y), n);
            if divisor > Integer::from(1) && divisor < *n && n % &divisor == Integer::from(0) {
                return Some(divisor);
            }
            let divisor = integer_gcd(&(&x + &y), n);
            if divisor > Integer::from(1) && divisor < *n && n % &divisor == Integer::from(0) {
                return Some(divisor);
            }
        }
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn extract_factor(
    n: &Integer,
    field: &NumberField,
    evaluation_point: &Integer,
    relations: &[GeneralRelation],
    dependency: &[usize],
    units: &[Element],
    base: &GeneralBase,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
    units_enabled: bool,
) -> Option<Integer> {
    if let Some(divisor) =
        direct_congruence_of_squares(n, evaluation_point, relations, dependency, deadline)
    {
        return Some(divisor);
    }
    let mut rational = Integer::from(1);
    let mut algebraic = Element::one(field.degree);
    for index in dependency {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return None;
        }
        let relation = relations.get(*index)?;
        rational *= &relation.rational_value;
        algebraic = field.multiply(&algebraic, &relation.algebraic_value);
    }
    let alpha_bits = algebraic
        .coefficients
        .iter()
        .map(|coefficient| integer_abs(coefficient).bits())
        .max()
        .unwrap_or(0);
    if alpha_bits > MAX_SQUARE_ROOT_BITS {
        return None;
    }
    let started = Instant::now();
    let attempt_deadline = attempt_deadline(deadline, started, alpha_bits);
    let deadline = attempt_deadline;
    let rational = crate::integer::integer_abs(&rational);
    let rational_root = integer_sqrt(&rational);
    if &rational_root * &rational_root != rational {
        return None;
    }
    let corrected =
        match unit_corrected_square_root(field, &algebraic, units, deadline, stop, units_enabled) {
            Some(corrected) => corrected,
            None => {
                if let Some(result) = ideal_corrected_square_root(
                    field, &algebraic, relations, dependency, base, deadline, stop,
                ) {
                    result
                } else {
                    CorrectedSquareRoot {
                        root: algebraic_square_root(field, &algebraic, deadline, stop)?,
                        unit: Element::one(field.degree),
                        power: 0,
                    }
                }
            }
        };
    let multiplier = field.evaluate(&corrected.unit, evaluation_point);
    let targets = multiplier_square_roots(&rational_root, &multiplier, corrected.power);
    for sign in [Integer::from(1), Integer::from(-1)] {
        let mut root = corrected.root.clone();
        for coefficient in &mut root.coefficients {
            *coefficient *= &sign;
        }
        let algebraic_value = field.evaluate(&root, evaluation_point);
        let x = &algebraic_value % n;
        for target in &targets {
            let y = target % n;
            if (&x * &x - &y * &y) % n != Integer::from(0) {
                continue;
            }
            let divisor = integer_gcd(&(&x - &y), n);
            if divisor > Integer::from(1) && divisor < *n && n % &divisor == Integer::from(0) {
                return Some(divisor);
            }
            let divisor = integer_gcd(&(&x + &y), n);
            if divisor > Integer::from(1) && divisor < *n && n % &divisor == Integer::from(0) {
                return Some(divisor);
            }
        }
    }
    None
}

fn ideal_corrected_square_root(
    field: &NumberField,
    value: &Element,
    relations: &[GeneralRelation],
    dependency: &[usize],
    base: &GeneralBase,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> Option<CorrectedSquareRoot> {
    let mut columns = HashMap::new();
    for (prime, root, column) in &base.algebraic_roots {
        if *column >= base.rational_primes.len() && *column < base.sign_column {
            columns.insert(*column, (*prime, *root));
        }
    }
    let mut key_set = HashSet::new();
    for index in dependency {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return None;
        }
        let relation = relations.get(*index)?;
        for (word_index, word) in relation.mask.iter().enumerate() {
            let mut rest = *word;
            while rest != 0 {
                let column = word_index * 64 + rest.trailing_zeros() as usize;
                if let Some(key) = columns.get(&column) {
                    key_set.insert(*key);
                }
                rest &= rest - 1;
            }
        }
    }
    let mut keys: Vec<(u64, u64)> = key_set.into_iter().collect();
    keys.sort_unstable();
    if keys.is_empty() {
        return None;
    }
    let mut generators: Vec<(u64, u64, Element, usize)> = Vec::with_capacity(keys.len());
    for (prime, root) in keys {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return None;
        }
        if let Some(vector) = shortest_vector(
            ideal_basis(field.degree, prime, root, field, deadline, stop),
            deadline,
            stop,
        ) {
            let size = usize::try_from(
                vector
                    .iter()
                    .map(|coefficient| integer_abs(coefficient).bits())
                    .max()
                    .unwrap_or(0),
            )
            .unwrap_or(usize::MAX);
            generators.push((prime, root, Element::from_coefficients(vector), size));
        }
    }
    generators.sort_by(|left, right| {
        left.3
            .cmp(&right.3)
            .then(left.0.cmp(&right.0))
            .then(left.1.cmp(&right.1))
    });
    let mut generator = Element::one(field.degree);
    for (_, _, element, _) in &generators {
        generator = field.multiply(&generator, element);
    }
    let generator_square = field.multiply(&generator, &generator);
    if let Some(inverse) = element_inverse(field, &generator_square) {
        let target = field.multiply(value, &inverse);
        if let Some(root) = algebraic_square_root(field, &target, deadline, stop)
            && field.multiply(&root, &root) == target
        {
            return Some(CorrectedSquareRoot {
                root,
                unit: generator,
                power: 2,
            });
        }
    }
    let mut target = value.clone();
    for power in 1..=2 {
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            return None;
        }
        target = field.multiply(&target, &generator);
        if let Some(root) = algebraic_square_root(field, &target, deadline, stop)
            && field.multiply(&root, &root) == target
        {
            return Some(CorrectedSquareRoot {
                root,
                unit: Element::one(field.degree),
                power,
            });
        }
    }
    let value_norm = element_norm(field, value);
    let capped = generators.len().min(SUBSET_GENERATOR_LIMIT);
    if capped >= 2 {
        let full = (1_usize << capped) - 1;
        let mut masks: Vec<usize> = (1..full).collect();
        masks.sort_by_key(|mask| (mask.count_ones(), *mask));
        let budgeted = if value_norm.is_some() && capped <= 6 {
            masks.len()
        } else {
            SUBSET_TRIAL_LIMIT
        };
        for (trials, mask) in masks.into_iter().enumerate() {
            if trials >= budgeted {
                break;
            }
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            let mut partial = Element::one(field.degree);
            for (index, (_, _, element, _)) in generators.iter().enumerate().take(capped) {
                if mask & (1_usize << index) != 0 {
                    partial = field.multiply(&partial, element);
                }
            }
            if let Some(norm) = value_norm.as_ref()
                && let Some(partial_norm) = element_norm(field, &partial)
                && !partial_norm.is_zero()
            {
                let square_norm = &partial_norm * &partial_norm;
                let divisible = (norm % &square_norm).is_zero();
                let magnitude = integer_abs(&(norm / &square_norm));
                let root = integer_sqrt(&magnitude);
                if !divisible || &root * &root != magnitude {
                    continue;
                }
            }
            let square = field.multiply(&partial, &partial);
            if let Some(inverse) = element_inverse(field, &square) {
                let target = field.multiply(value, &inverse);
                if let Some(root) = algebraic_square_root(field, &target, deadline, stop)
                    && field.multiply(&root, &root) == target
                {
                    return Some(CorrectedSquareRoot {
                        root,
                        unit: partial,
                        power: 2,
                    });
                }
            }
        }
    }
    None
}

struct CorrectedSquareRoot {
    root: Element,
    unit: Element,
    power: usize,
}

fn unit_search_bound(degree: usize) -> i64 {
    if degree <= 3 {
        16_i64
    } else if degree == 4 {
        6_i64
    } else if degree == 5 {
        2_i64
    } else {
        0_i64
    }
}

fn find_units(field: &NumberField, bound: i64, limit: usize) -> Vec<Element> {
    if !matches!(field.degree, 3..=5) || bound < 1 || limit == 0 {
        return Vec::new();
    }
    let range: Vec<Integer> = (-bound..=bound).map(Integer::from).collect();
    let mut current = vec![Integer::from(0); field.degree];
    let mut units = Vec::new();
    collect_units(field, &range, 0, &mut current, limit, &mut units);
    units
}

fn collect_units(
    field: &NumberField,
    range: &[Integer],
    index: usize,
    current: &mut Vec<Integer>,
    limit: usize,
    units: &mut Vec<Element>,
) {
    if units.len() >= limit {
        return;
    }
    if index == current.len() {
        if current.iter().all(|value| value.is_zero()) {
            return;
        }
        let candidate = Element::from_coefficients(current.clone());
        if element_norm(field, &candidate)
            .is_some_and(|norm| norm.is_one() || norm == Integer::from(-1))
        {
            units.push(candidate);
        }
        return;
    }
    for value in range {
        current[index] = value.clone();
        collect_units(field, range, index + 1, current, limit, units);
        if units.len() >= limit {
            return;
        }
    }
}

fn element_norm(field: &NumberField, element: &Element) -> Option<Integer> {
    if field.degree > 6 || element.coefficients.len() != field.degree {
        return None;
    }
    exact_determinant(&multiplication_matrix(field, element))
}

fn multiplication_matrix(field: &NumberField, element: &Element) -> Vec<Vec<Integer>> {
    let degree = field.degree;
    let mut columns = Vec::with_capacity(degree);
    for column in 0..degree {
        let mut basis = vec![Integer::from(0); degree];
        basis[column] = Integer::from(1);
        columns.push(
            field
                .multiply(element, &Element::from_coefficients(basis))
                .coefficients,
        );
    }
    (0..degree)
        .map(|row| columns.iter().map(|column| column[row].clone()).collect())
        .collect()
}

fn exact_determinant(matrix: &[Vec<Integer>]) -> Option<Integer> {
    let size = matrix.len();
    if size == 0 || matrix.iter().any(|row| row.len() != size) {
        return None;
    }
    let mut work = matrix.to_vec();
    let mut sign = Integer::from(1);
    let mut divisor = Integer::from(1);
    for index in 0..size {
        let mut pivot = index;
        while pivot < size && work[pivot][index].is_zero() {
            pivot += 1;
        }
        if pivot == size {
            return Some(Integer::from(0));
        }
        if pivot != index {
            work.swap(pivot, index);
            sign = -sign;
        }
        let diagonal = work[index][index].clone();
        for row in (index + 1)..size {
            for column in (index + 1)..size {
                let numerator =
                    &work[row][column] * &diagonal - &work[row][index] * &work[index][column];
                work[row][column] = &numerator / &divisor;
            }
            work[row][index] = Integer::from(0);
        }
        divisor = diagonal;
    }
    Some(sign * work[size - 1][size - 1].clone())
}

fn element_inverse(field: &NumberField, element: &Element) -> Option<Element> {
    if field.degree > 6 || element.coefficients.len() != field.degree {
        return None;
    }
    let matrix = multiplication_matrix(field, element);
    let determinant = exact_determinant(&matrix)?;
    if determinant.is_zero() {
        return None;
    }
    let mut coefficients = vec![Integer::from(0); field.degree];
    for column in 0..field.degree {
        let mut replaced = matrix.clone();
        for (row, values) in replaced.iter_mut().enumerate() {
            values[column] = if row == column {
                Integer::from(1)
            } else {
                Integer::from(0)
            };
        }
        let numerator = exact_determinant(&replaced)?;
        if (&numerator % &determinant) != Integer::from(0) {
            return None;
        }
        coefficients[column] = numerator / &determinant;
    }
    Some(Element::from_coefficients(coefficients))
}

fn unit_corrected_targets(
    field: &NumberField,
    value: &Element,
    units: &[Element],
) -> Vec<(Element, Element, usize)> {
    let mut multipliers: Vec<Element> = vec![Element::one(field.degree)];
    for unit in units {
        if !multipliers.contains(unit) {
            multipliers.push(unit.clone());
        }
    }
    let mut rows: Vec<Vec<Element>> = Vec::with_capacity(multipliers.len());
    for unit in &multipliers {
        let mut row = Vec::with_capacity(UNIT_POWER_LIMIT);
        row.push(value.clone());
        for power in 1..UNIT_POWER_LIMIT {
            let previous = row[power - 1].clone();
            row.push(field.multiply(&previous, unit));
        }
        rows.push(row);
    }
    let mut attempted: Vec<Element> = Vec::new();
    let mut targets = Vec::new();
    for power in 0..UNIT_POWER_LIMIT {
        for (index, row) in rows.iter().enumerate() {
            if attempted.len() >= UNIT_TARGET_LIMIT {
                break;
            }
            let scaled = row[power].clone();
            if attempted.contains(&scaled) {
                continue;
            }
            attempted.push(scaled.clone());
            targets.push((scaled, multipliers[index].clone(), power));
        }
    }
    targets
}

fn algebraic_square_root(
    field: &NumberField,
    value: &Element,
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
) -> Option<Element> {
    if value.coefficients.len() == field.degree {
        let norm = field.norm(value);
        if norm < Integer::from(0) {
            return None;
        }
        let root = integer_sqrt(&norm);
        if &root * &root != norm {
            return None;
        }
    }
    if let Some(root) = crate::cubic::find_square_root(field, value, deadline, stop) {
        return Some(root);
    }
    if field.degree == 3 {
        return None;
    }
    crate::hensel::find_square_root(field, value, deadline, stop)
}

fn unit_corrected_square_root(
    field: &NumberField,
    value: &Element,
    units: &[Element],
    deadline: Option<Instant>,
    stop: Option<&AtomicBool>,
    units_enabled: bool,
) -> Option<CorrectedSquareRoot> {
    if let Some(root) = algebraic_square_root(field, value, deadline, stop) {
        return Some(CorrectedSquareRoot {
            root,
            unit: Element::one(field.degree),
            power: 0,
        });
    }
    if !units_enabled {
        return None;
    }
    for (scaled, unit, power) in unit_corrected_targets(field, value, units) {
        if let Some(root) = algebraic_square_root(field, &scaled, deadline, stop)
            && field.multiply(&root, &root) == scaled
        {
            return Some(CorrectedSquareRoot { root, unit, power });
        }
    }
    None
}

fn multiplier_square_roots(
    rational_root: &Integer,
    multiplier: &Integer,
    power: usize,
) -> Vec<Integer> {
    let mut value = rational_root * rational_root;
    for _ in 0..power {
        value *= multiplier;
    }
    let root = integer_sqrt(&value);
    if &root * &root != value {
        return Vec::new();
    }
    vec![root.clone(), -root]
}

pub(crate) fn integer_root(value: &Integer, degree: usize) -> Integer {
    if degree == 0 {
        return Integer::from(1);
    }
    if degree == 1 {
        return value.clone();
    }
    if value <= &Integer::from(0) {
        return Integer::from(0);
    }
    if value <= &Integer::from(1) {
        return value.clone();
    }
    let shift = value.bits().div_ceil(u64::try_from(degree).unwrap_or(1)) + 1;
    let mut lo = Integer::from(0);
    let mut hi = Integer::from(1) << shift;
    let exponent = u32::try_from(degree).unwrap_or(u32::MAX);
    while &lo + Integer::from(1) < hi {
        let mid = (&lo + &hi) / Integer::from(2);
        if mid.pow(exponent) <= *value {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

fn integer_pow(value: &Integer, degree: usize) -> Integer {
    value.pow(degree as u32)
}

fn is_perfect_power(value: &Integer, degree: usize) -> bool {
    if value <= &Integer::from(1) {
        return false;
    }
    let root = integer_root(value, degree);
    root.pow(degree as u32) == *value
}

fn discriminant(polynomial: &PolynomialPair) -> Integer {
    polynomial_discriminant(&polynomial.coefficients)
}

fn mark_general_residue(scores: &mut [u32], max_a: usize, prime: u64, residue: u64, weight: u32) {
    let start = ((u128::from(residue) + max_a as u128) % u128::from(prime)) as usize;
    let step = usize::try_from(prime).unwrap_or(scores.len()).max(1);
    let mut index = start;
    while index < scores.len() {
        scores[index] = scores[index].saturating_add(weight);
        if step >= scores.len() - index {
            break;
        }
        index += step;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn negative_residue_maps_to_prime_minus_one() {
        assert_eq!(residue_mod(&Integer::from(-1), 7), 6);
        assert_eq!(residue_mod(&Integer::from(-8), 5), 2);
        assert_eq!(residue_mod(&Integer::from(0), 7), 0);
        assert_eq!(residue_mod(&Integer::from(8), 7), 1);
    }

    #[test]
    fn signed_coefficient_residues_match_euclidean_remainder() {
        let coefficients = [Integer::from(-1), Integer::from(0), Integer::from(1)];
        let residues: Vec<u64> = coefficients
            .iter()
            .map(|value| residue_mod(value, 7))
            .collect();
        assert_eq!(residues, vec![6, 0, 1]);
    }

    #[test]
    fn polynomial_discriminant_matches_known_values() {
        let cubic_pure = vec![
            Integer::from(-2),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ];
        assert_eq!(polynomial_discriminant(&cubic_pure), Integer::from(-108));
        let cubic_general = vec![
            Integer::from(1),
            Integer::from(1),
            Integer::from(0),
            Integer::from(1),
        ];
        assert_eq!(polynomial_discriminant(&cubic_general), Integer::from(-31));
        let quadratic = vec![Integer::from(1), Integer::from(0), Integer::from(1)];
        assert_eq!(polynomial_discriminant(&quadratic), Integer::from(-4));
        let quartic_pure = vec![
            Integer::from(2),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ];
        assert_eq!(polynomial_discriminant(&quartic_pure), Integer::from(2048));
    }

    #[test]
    fn selected_polynomial_has_no_small_index_divisor() {
        let n = Integer::from(100160063);
        for degree in 3..=5 {
            let polynomial = select_polynomial(&n, degree, 64).expect("polynomial");
            assert!(has_unit_index(&polynomial));
        }
    }

    #[test]
    fn irreducibility_gate_rescues_quartic_hidden_from_small_primes() {
        let hidden = PolynomialPair {
            degree: 4,
            coefficients: vec![
                Integer::from(105),
                Integer::from(15),
                Integer::from(6),
                Integer::from(9),
                Integer::from(1),
            ],
            m: Integer::from(10),
        };
        assert!(is_irreducible_candidate(&hidden));
    }

    #[test]
    fn irreducibility_gate_still_rejects_reducible_quartic() {
        let reducible = PolynomialPair {
            degree: 4,
            coefficients: vec![
                Integer::from(2),
                Integer::from(0),
                Integer::from(3),
                Integer::from(0),
                Integer::from(1),
            ],
            m: Integer::from(10),
        };
        assert!(!is_irreducible_candidate(&reducible));
    }

    #[test]
    fn irreducibility_gate_rejects_reducible_cubic() {
        let reducible = PolynomialPair {
            degree: 3,
            coefficients: vec![
                Integer::from(1),
                Integer::from(1),
                Integer::from(1),
                Integer::from(1),
            ],
            m: Integer::from(10),
        };
        assert!(!is_irreducible_candidate(&reducible));
    }

    #[test]
    fn irreducibility_gate_accepts_irreducible_cubic() {
        let irreducible = PolynomialPair {
            degree: 3,
            coefficients: vec![
                Integer::from(1),
                Integer::from(1),
                Integer::from(0),
                Integer::from(1),
            ],
            m: Integer::from(10),
        };
        assert!(is_irreducible_candidate(&irreducible));
    }

    #[test]
    fn number_field_construction_does_not_prove_irreducibility() {
        let coefficients = vec![
            Integer::from(2),
            Integer::from(0),
            Integer::from(3),
            Integer::from(0),
            Integer::from(1),
        ];
        assert!(NumberField::new(4, coefficients.clone()).is_some());
        let candidate = PolynomialPair {
            degree: 4,
            coefficients,
            m: Integer::from(10),
        };
        assert!(!is_irreducible_candidate(&candidate));
    }

    #[test]
    fn selected_polynomials_are_irreducible_over_rationals() {
        let n = Integer::from(100160063);
        for degree in 3..=5 {
            let polynomial = select_polynomial(&n, degree, 64).expect("polynomial");
            assert_eq!(polynomial.degree, degree);
            assert_eq!(polynomial.evaluate(&polynomial.m), n);
            assert!(is_irreducible_candidate(&polynomial));
            assert!(has_unit_index(&polynomial));
            NumberField::new(degree, polynomial.coefficients.clone()).expect("number field");
        }
    }

    fn pure_cubic(constant: i64) -> NumberField {
        NumberField::new(
            3,
            vec![
                Integer::from(constant),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ],
        )
        .expect("cubic field")
    }

    fn cubic(values: &[i64]) -> Element {
        Element::from_coefficients(values.iter().copied().map(Integer::from).collect())
    }

    fn pure_cubic_norm(values: &[i64], constant: i64) -> Integer {
        let a = Integer::from(constant);
        let cube = |value: i64| Integer::from(value).pow(3);
        let linear = |value: i64| Integer::from(value);
        cube(values[0]) - &a * cube(values[1])
            + (&a * &a) * cube(values[2])
            + Integer::from(3 * constant)
                * linear(values[0])
                * linear(values[1])
                * linear(values[2])
    }

    fn minor_determinant(matrix: &[Vec<Integer>], skip_row: usize, skip_column: usize) -> Integer {
        let reduced: Vec<Vec<Integer>> = matrix
            .iter()
            .enumerate()
            .filter(|(row, _)| *row != skip_row)
            .map(|(_, row)| {
                row.iter()
                    .enumerate()
                    .filter(|(column, _)| *column != skip_column)
                    .map(|(_, value)| value.clone())
                    .collect()
            })
            .collect();
        exact_determinant(&reduced).expect("minor determinant")
    }

    #[test]
    fn determinant_norm_matches_pure_cubic_form() {
        let field = pure_cubic(2);
        for values in [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 1],
            [1, -1, 1],
            [-3, 5, 2],
            [7, -4, -1],
        ] {
            assert_eq!(
                element_norm(&field, &cubic(&values)),
                Some(pure_cubic_norm(&values, 2))
            );
        }
        let other = pure_cubic(7);
        for values in [[1, 2, 0], [0, -1, 3], [4, 4, -1]] {
            assert_eq!(
                element_norm(&other, &cubic(&values)),
                Some(pure_cubic_norm(&values, 7))
            );
        }
    }

    #[test]
    fn norm_is_multiplicative() {
        let field = pure_cubic(2);
        let left = cubic(&[3, 1, 0]);
        let right = cubic(&[-2, 4, 1]);
        let product = field.multiply(&left, &right);
        let left_norm = element_norm(&field, &left).expect("left norm");
        let right_norm = element_norm(&field, &right).expect("right norm");
        assert_eq!(element_norm(&field, &product), Some(left_norm * right_norm));
    }

    #[test]
    fn quartic_norm_comes_from_determinant() {
        let field = NumberField::new(
            4,
            vec![
                Integer::from(2),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ],
        )
        .expect("quartic field");
        let theta = Element::from_coefficients(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
        ]);
        let one_plus_theta = Element::from_coefficients(vec![
            Integer::from(1),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
        ]);
        assert_eq!(element_norm(&field, &theta), Some(Integer::from(2)));
        assert_eq!(
            element_norm(&field, &one_plus_theta),
            Some(Integer::from(3))
        );
        let mixed = Element::from_coefficients(vec![
            Integer::from(2),
            Integer::from(-1),
            Integer::from(3),
            Integer::from(1),
        ]);
        assert_eq!(element_norm(&field, &mixed), Some(Integer::from(662)));
    }

    #[test]
    fn cubic_unit_has_integral_inverse() {
        let field = pure_cubic(2);
        let unit = cubic(&[1, -1, 1]);
        assert_eq!(element_norm(&field, &unit), Some(Integer::from(1)));
        let matrix = multiplication_matrix(&field, &unit);
        let size = matrix.len();
        let inverse = Element::from_coefficients(
            (0..size)
                .map(|row| {
                    let minor = minor_determinant(&matrix, 0, row);
                    if row % 2 == 0 { minor } else { -minor }
                })
                .collect(),
        );
        assert_eq!(field.multiply(&unit, &inverse), Element::one(3));
    }

    #[test]
    fn unit_search_finds_cubic_unit_of_norm_one() {
        let field = pure_cubic(2);
        let units = find_units(&field, 2, 8);
        assert!(units.contains(&cubic(&[1, -1, 1])));
        assert!(!units.is_empty());
        for unit in &units {
            assert!(unit.coefficients.iter().any(|value| !value.is_zero()));
            assert!(
                element_norm(&field, unit)
                    .is_some_and(|norm| norm.is_one() || norm == Integer::from(-1))
            );
        }
    }

    #[test]
    fn unit_correction_returns_verified_root() {
        let field = pure_cubic(2);
        let root = cubic(&[3, 1, 0]);
        let value = field.multiply(&root, &root);
        let units = find_units(&field, 2, 8);
        let corrected = unit_corrected_square_root(&field, &value, &units, None, None, true)
            .expect("algebraic square root");
        assert_eq!(corrected.power, 0);
        assert_eq!(corrected.unit, Element::one(3));
        assert_eq!(field.multiply(&corrected.root, &corrected.root), value);
    }

    #[test]
    fn unit_targets_cover_the_small_power_range() {
        let field = pure_cubic(2);
        let value = cubic(&[9, 6, 1]);
        let unit = cubic(&[1, -1, 1]);
        let targets = unit_corrected_targets(&field, &value, std::slice::from_ref(&unit));
        assert_eq!(targets.len(), UNIT_POWER_LIMIT);
        let mut expected = value.clone();
        for (index, (scaled, multiplier, power)) in targets.iter().enumerate() {
            assert_eq!(*scaled, expected);
            assert_eq!(*power, index);
            if index == 0 {
                assert_eq!(*multiplier, Element::one(3));
            } else {
                assert_eq!(*multiplier, unit);
            }
            expected = field.multiply(&expected, &unit);
        }
    }

    #[test]
    fn unit_targets_deduplicate_the_identity_multiplier() {
        let field = pure_cubic(2);
        let value = cubic(&[9, 6, 1]);
        let units = vec![Element::one(3), cubic(&[1, 0, 0]), cubic(&[1, -1, 1])];
        let targets = unit_corrected_targets(&field, &value, &units);
        let mut scaled_values: Vec<Vec<Integer>> = targets
            .iter()
            .map(|(scaled, _, _)| scaled.coefficients.clone())
            .collect();
        let total = scaled_values.len();
        scaled_values.sort();
        scaled_values.dedup();
        assert_eq!(scaled_values.len(), total);
        assert_eq!(targets[0].0, value);
    }

    #[test]
    fn multiplier_square_roots_respect_the_unit_power() {
        let field = pure_cubic(2);
        let unit = cubic(&[1, -1, 1]);
        let point = Integer::from(7);
        let multiplier = field.evaluate(&unit, &point);
        let rational_root = Integer::from(5);
        let identity = multiplier_square_roots(&rational_root, &multiplier, 0);
        let mirrored = -rational_root.clone();
        assert_eq!(identity, vec![rational_root.clone(), mirrored]);
        assert!(multiplier_square_roots(&rational_root, &multiplier, 1).is_empty());
        let squared = &rational_root * &rational_root * &multiplier * &multiplier;
        let expected = integer_sqrt(&squared);
        let negative = -expected.clone();
        assert_eq!(
            multiplier_square_roots(&rational_root, &multiplier, 2),
            vec![expected.clone(), negative]
        );
        assert_eq!(&expected * &expected, squared);
    }

    fn pure_quintic(constant: i64) -> NumberField {
        NumberField::new(
            5,
            vec![
                Integer::from(constant),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ],
        )
        .expect("quintic field")
    }

    fn quintic(values: &[i64]) -> Element {
        Element::from_coefficients(values.iter().copied().map(Integer::from).collect())
    }

    #[test]
    fn quintic_unit_search_stays_bounded() {
        let field = pure_quintic(7);
        let bound = unit_search_bound(5);
        assert!((1..=3).contains(&bound));
        let width = (2 * bound + 1) as u64;
        assert!(width.pow(5) <= 20_000);
        let units = find_units(&field, bound, UNIT_SEARCH_LIMIT);
        assert!(!units.is_empty());
        assert!(units.len() <= UNIT_SEARCH_LIMIT);
        for unit in &units {
            assert!(
                element_norm(&field, unit)
                    .is_some_and(|norm| norm.is_one() || norm == Integer::from(-1))
            );
        }
    }

    #[test]
    fn quintic_unit_correction_succeeds() {
        let field = pure_quintic(7);
        let beta = quintic(&[1, 2, 3, 4, 5]);
        let value = field.multiply(&beta, &beta);
        let units = find_units(&field, unit_search_bound(5), UNIT_SEARCH_LIMIT);
        assert!(!units.is_empty());
        let corrected = unit_corrected_square_root(&field, &value, &units, None, None, true)
            .expect("quintic square root");
        assert_eq!(field.multiply(&corrected.root, &corrected.root), value);
    }

    #[test]
    fn quintic_unit_correction_recovers_unit_multiple() {
        let field = pure_quintic(7);
        let beta = quintic(&[1, 2, 3, 4, 5]);
        let square = field.multiply(&beta, &beta);
        let negated: Vec<Integer> = square
            .coefficients
            .iter()
            .map(|value| -value.clone())
            .collect();
        let value = Element::from_coefficients(negated);
        assert!(algebraic_square_root(&field, &value, None, None).is_none());
        let units = find_units(&field, unit_search_bound(5), UNIT_SEARCH_LIMIT);
        let minus_one = Element::from_coefficients(vec![
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ]);
        assert!(units.contains(&minus_one));
        let corrected = unit_corrected_square_root(&field, &value, &units, None, None, true)
            .expect("unit corrected quintic root");
        assert!(corrected.power > 0);
        let mut expected = value.clone();
        for _ in 0..corrected.power {
            expected = field.multiply(&expected, &corrected.unit);
        }
        assert_eq!(field.multiply(&corrected.root, &corrected.root), expected);
        assert_eq!(expected, square);
    }

    #[test]
    fn quintic_direct_sqrt_fails_but_unit_scan_succeeds() {
        let field = pure_quintic(7);
        assert_eq!(unit_search_bound(5), 2);
        assert_eq!(UNIT_SEARCH_LIMIT, 4);
        let beta = quintic(&[1, 2, 3, 4, 5]);
        let square = field.multiply(&beta, &beta);
        let negated: Vec<Integer> = square
            .coefficients
            .iter()
            .map(|value| -value.clone())
            .collect();
        let value = Element::from_coefficients(negated);
        assert!(algebraic_square_root(&field, &value, None, None).is_none());
        assert!(unit_corrected_square_root(&field, &value, &[], None, None, false).is_none());
        let units = find_units(&field, unit_search_bound(5), UNIT_SEARCH_LIMIT);
        assert!(!units.is_empty());
        assert!(units.len() <= UNIT_SEARCH_LIMIT);
        let corrected = unit_corrected_square_root(&field, &value, &units, None, None, true)
            .expect("unit scan root");
        assert!(corrected.power > 0);
        let mut expected = value.clone();
        for _ in 0..corrected.power {
            expected = field.multiply(&expected, &corrected.unit);
        }
        assert_eq!(field.multiply(&corrected.root, &corrected.root), expected);
        assert_eq!(expected, square);
    }

    #[test]
    fn gnfs_failure_taxonomy_labels_are_distinct() {
        assert_eq!(scan_failure_taxonomy(false, false), "sqrt_failed");
        assert_eq!(scan_failure_taxonomy(true, false), "unit_correction_failed");
        assert_eq!(scan_failure_taxonomy(false, true), "found");
        assert_eq!(scan_failure_taxonomy(true, true), "found");
        let labels = [
            "insufficient_relations",
            "rank_deficient",
            "sqrt_failed",
            "unit_correction_failed",
        ];
        let mut sorted = labels.to_vec();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), labels.len());
        assert_eq!(unit_search_bound(5), 2);
        let field = pure_quintic(7);
        let units = find_units(&field, unit_search_bound(5), UNIT_SEARCH_LIMIT);
        assert!(!units.is_empty());
        assert!(units.len() <= UNIT_SEARCH_LIMIT);
        let detail = format!(
            "units_tried={} bound={} found=false outcome={}",
            units.len(),
            unit_search_bound(5),
            scan_failure_taxonomy(true, false)
        );
        assert!(detail.contains("units_tried="));
        assert!(detail.contains("unit_correction_failed"));
    }

    #[test]
    fn zero_memory_budget_rejects_live_estimate() {
        let options = FactorOptions {
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        let need = crate::options::estimate_live_bytes(600, 900, 10, 5);
        assert!(need > 0);
        let outcome = live_budget_exceeded(&options, 600, 900, 10, 5);
        assert!(
            matches!(outcome, Some(FactorError::ResourceLimit(_))),
            "zero budget must reject a positive live estimate"
        );
        assert!(live_budget_exceeded(&FactorOptions::default(), 10, 15, 0, 0).is_none());
    }

    #[test]
    fn zero_memory_budget_rejects_relation_matrix() {
        let n = Integer::from(100160063);
        let options = FactorOptions {
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        let cancelled = AtomicBool::new(false);
        let outcome = gnfs_divisor(&n, &options, &cancelled, None);
        assert!(
            matches!(outcome, Err(FactorError::ResourceLimit(_))),
            "zero budget must fail before sieving"
        );
    }

    #[test]
    fn default_memory_budget_covers_small_matrix() {
        let options = FactorOptions::default();
        assert!(
            crate::options::estimate_relation_matrix_megabytes(600) as u64
                <= options.memory_budget_bytes()
        );
    }

    #[test]
    fn total_bytes_estimate_exceeds_tiny_budget() {
        let options = FactorOptions {
            max_memory_mb: Some(1),
            ..FactorOptions::default()
        };
        let budget = options.memory_budget_bytes();
        assert!(estimate_total_bytes(200_000, 300_000) > budget);
        assert!(estimate_total_bytes(10, 15) <= FactorOptions::default().memory_budget_bytes());
    }

    #[test]
    fn integer_root_matches_floor_on_hang_adversaries() {
        assert_eq!(integer_root(&Integer::from(1000), 3), Integer::from(10));
        assert_eq!(integer_root(&Integer::from(1330), 3), Integer::from(10));
        for value in [1000, 1330] {
            let n = Integer::from(value);
            let root = integer_root(&n, 3);
            assert!(root.pow(3) <= n);
            assert!((&root + Integer::from(1)).pow(3) > n);
        }
    }

    #[test]
    fn integer_root_recovers_large_exact_cubes() {
        let base = Integer::from(123456789);
        let cube = base.pow(3);
        assert_eq!(integer_root(&cube, 3), base);
        assert!(is_perfect_power(&cube, 3));
        let wide_base = Integer::from(123456789) << 60_u32;
        let wide_cube = wide_base.pow(3);
        assert_eq!(integer_root(&wide_cube, 3), wide_base);
        assert!(is_perfect_power(&wide_cube, 3));
    }

    fn factor_style_corrected_root(value: &Integer, exponent: u32) -> Integer {
        let mut root = integer_root(value, exponent as usize);
        while root.pow(exponent) > *value {
            root -= Integer::from(1);
        }
        while (&root + Integer::from(1)).pow(exponent) <= *value {
            root += Integer::from(1);
        }
        root
    }

    #[test]
    fn integer_root_brackets_two_pow_128_minus_one() {
        let n = (Integer::from(1) << 128_u32) - Integer::from(1);
        let root = factor_style_corrected_root(&n, 3);
        assert!(root.pow(3) <= n);
        assert!((&root + Integer::from(1)).pow(3) > n);
        assert!(!is_perfect_power(&n, 3));
    }

    #[test]
    fn integer_root_hang_regression_pins_complete() {
        let wide_base = Integer::from(123456789) << 60_u32;
        let wide_cube = wide_base.pow(3);
        assert_eq!(factor_style_corrected_root(&wide_cube, 3), wide_base);
        let semi = (Integer::from(1) << 64_u32) + Integer::from(13);
        let other = (Integer::from(1) << 64_u32) + Integer::from(29);
        let n = &semi * &other;
        for exponent in [3_u32, 5, 7, 31, 127] {
            let root = factor_style_corrected_root(&n, exponent);
            assert!(root.pow(exponent) <= n);
            assert!((&root + Integer::from(1)).pow(exponent) > n);
        }
        assert!(!is_perfect_power(&n, 3));
        assert!(!is_perfect_power(&n, 5));
    }

    #[test]
    fn total_bytes_estimate_covers_drain_growth() {
        let tiny = FactorOptions {
            max_memory_mb: Some(1),
            ..FactorOptions::default()
        }
        .memory_budget_bytes();
        let ample = FactorOptions::default().memory_budget_bytes();
        assert!(estimate_total_bytes(200_000, 300_016) > tiny);
        assert!(estimate_total_bytes(600, 916) <= ample);
        assert!(estimate_total_bytes(200_000, 300_016) >= estimate_total_bytes(200_000, 300_000));
    }
}
