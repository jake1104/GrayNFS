use std::collections::{HashMap, HashSet};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::FactorError;
use crate::general::{GeneralBase, GeneralRelation, PartialBank};
use crate::integer::Integer;
use crate::number_field::Element;
use crate::options::FactorOptions;
use crate::polynomial::PolynomialPair;

pub const CHECKPOINT_VERSION: u32 = 1;

type MixedPendingEntry = (((u64, u64), u64), usize);
type TriadPendingEntry = ((Vec<(u64, u64)>, u64), usize);

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StoredPolynomial {
    degree: usize,
    coefficients: Vec<String>,
    m: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StoredBase {
    rational_primes: Vec<u64>,
    algebraic_roots: Vec<(u64, u64, usize)>,
    modulus_primes: Vec<u64>,
    discriminant: String,
    algebraic_bound: u64,
    large_prime_bound: u64,
    sign_column: usize,
    lp_columns: usize,
    algebraic_index: Vec<(u64, Vec<(u64, usize)>)>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StoredRelation {
    rational_value: String,
    algebraic_value: Vec<String>,
    mask: Vec<u64>,
    large_prime: Option<(u64, u64)>,
    double_edge: Option<[(u64, u64); 2]>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StoredBank {
    singles: Vec<StoredRelation>,
    single_sources: Vec<(i64, i64)>,
    single_aux: Vec<u64>,
    single_lps: Vec<Vec<(u64, u64)>>,
    collector_pending: Vec<((u64, u64), usize)>,
    rat_pending: Vec<(u64, usize)>,
    mixed_pending: Vec<MixedPendingEntry>,
    triad_pending: Vec<TriadPendingEntry>,
    doubles: Vec<StoredRelation>,
    double_sources: Vec<(i64, i64)>,
    double_aux: Vec<u64>,
    double_parent: Vec<((u64, u64), (u64, u64))>,
    double_members: Vec<((u64, u64), Vec<usize>)>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StoredPosition {
    attempt_a: i64,
    attempt_b: i64,
    passes_done: usize,
    usable_total: usize,
    finished: bool,
    consumed: Vec<(i64, i64)>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StoredOptions {
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    large_prime_bound: usize,
    seed: u64,
    enable_double_large_primes: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CheckpointPayload {
    n: String,
    degree: usize,
    options: StoredOptions,
    polynomial: StoredPolynomial,
    base: StoredBase,
    complete: Vec<StoredRelation>,
    bank: StoredBank,
    position: StoredPosition,
}

#[derive(Serialize, Deserialize)]
struct CheckpointFile {
    version: u32,
    payload: CheckpointPayload,
    checksum: String,
}

pub struct CollectResume {
    pub(crate) complete: Vec<GeneralRelation>,
    pub(crate) bank: PartialBank,
    pub(crate) consumed: HashSet<(i64, i64)>,
    pub(crate) attempt_a: i64,
    pub(crate) attempt_b: i64,
    pub(crate) passes_done: usize,
    pub(crate) usable_total: usize,
    pub(crate) finished: bool,
}

impl CheckpointPayload {
    pub fn is_finished(&self) -> bool {
        self.position.finished
    }
}

pub fn checkpoint_checksum(payload_bytes: &[u8]) -> String {
    let mut hash: u64 = 0xcbf29ce484222325;
    for byte in payload_bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    format!("{hash:016x}")
}

pub fn save_gnfs_checkpoint(path: &Path, payload: &CheckpointPayload) -> Result<(), FactorError> {
    let payload_bytes = serde_json::to_vec(payload)
        .map_err(|error| FactorError::Checkpoint(format!("checkpoint encode failed: {error}")))?;
    let checksum = checkpoint_checksum(&payload_bytes);
    let file = CheckpointFile {
        version: CHECKPOINT_VERSION,
        payload: payload.clone(),
        checksum,
    };
    let bytes = serde_json::to_vec(&file)
        .map_err(|error| FactorError::Checkpoint(format!("checkpoint encode failed: {error}")))?;
    std::fs::write(path, bytes)
        .map_err(|error| FactorError::Checkpoint(format!("checkpoint write failed: {error}")))?;
    Ok(())
}

pub fn load_gnfs_checkpoint(path: &Path) -> Result<CheckpointPayload, FactorError> {
    let bytes = std::fs::read(path)
        .map_err(|error| FactorError::Checkpoint(format!("checkpoint unreadable: {error}")))?;
    let file: CheckpointFile = serde_json::from_slice(&bytes)
        .map_err(|error| FactorError::Checkpoint(format!("checkpoint corrupt: {error}")))?;
    if file.version != CHECKPOINT_VERSION {
        return Err(FactorError::Checkpoint(format!(
            "checkpoint version mismatch: found {}, supported {CHECKPOINT_VERSION}",
            file.version
        )));
    }
    let payload_bytes = serde_json::to_vec(&file.payload)
        .map_err(|error| FactorError::Checkpoint(format!("checkpoint corrupt: {error}")))?;
    let checksum = checkpoint_checksum(&payload_bytes);
    if checksum != file.checksum {
        return Err(FactorError::Checkpoint(
            "checkpoint checksum mismatch: file is corrupt or tampered".to_string(),
        ));
    }
    Ok(file.payload)
}

fn encode_integer(value: &Integer) -> String {
    value.to_string()
}

fn decode_integer(text: &str) -> Result<Integer, FactorError> {
    text.parse::<Integer>().map_err(|_| {
        FactorError::Checkpoint("checkpoint corrupt: invalid integer encoding".to_string())
    })
}

fn store_relation(relation: &GeneralRelation) -> StoredRelation {
    StoredRelation {
        rational_value: encode_integer(&relation.rational_value),
        algebraic_value: relation
            .algebraic_value
            .coefficients
            .iter()
            .map(encode_integer)
            .collect(),
        mask: relation.mask.clone(),
        large_prime: relation.large_prime,
        double_edge: relation.double_edge,
    }
}

fn restore_relation(stored: &StoredRelation) -> Result<GeneralRelation, FactorError> {
    let mut coefficients = Vec::with_capacity(stored.algebraic_value.len());
    for text in &stored.algebraic_value {
        coefficients.push(decode_integer(text)?);
    }
    Ok(GeneralRelation {
        rational_value: decode_integer(&stored.rational_value)?,
        algebraic_value: Element::from_coefficients(coefficients),
        mask: stored.mask.clone(),
        large_prime: stored.large_prime,
        double_edge: stored.double_edge,
    })
}

fn sorted_pairs<K: Ord + Clone, V: Clone>(map: &HashMap<K, V>) -> Vec<(K, V)> {
    let mut entries: Vec<(K, V)> = map
        .iter()
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect();
    entries.sort_by(|left, right| left.0.cmp(&right.0));
    entries
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn payload_from_state(
    n: &Integer,
    degree: usize,
    options: &FactorOptions,
    polynomial: &PolynomialPair,
    base: &GeneralBase,
    complete: &[GeneralRelation],
    bank: &PartialBank,
    consumed: &HashSet<(i64, i64)>,
    attempt_a: i64,
    attempt_b: i64,
    passes_done: usize,
    usable_total: usize,
    finished: bool,
) -> CheckpointPayload {
    let mut consumed_sorted: Vec<(i64, i64)> = consumed.iter().copied().collect();
    consumed_sorted.sort();
    let mut algebraic_index = sorted_pairs(&base.algebraic_index);
    for (_, roots) in algebraic_index.iter_mut() {
        roots.sort();
    }
    CheckpointPayload {
        n: encode_integer(n),
        degree,
        options: StoredOptions {
            rational_base_bound: options.rational_base_bound,
            algebraic_base_bound: options.algebraic_base_bound,
            sieve_a_bound: options.sieve_a_bound,
            sieve_b_bound: options.sieve_b_bound,
            large_prime_bound: options.large_prime_bound,
            seed: options.seed,
            enable_double_large_primes: options.enable_double_large_primes,
        },
        polynomial: StoredPolynomial {
            degree: polynomial.degree,
            coefficients: polynomial.coefficients.iter().map(encode_integer).collect(),
            m: encode_integer(&polynomial.m),
        },
        base: StoredBase {
            rational_primes: base.rational_primes.clone(),
            algebraic_roots: base.algebraic_roots.clone(),
            modulus_primes: base.modulus_primes.clone(),
            discriminant: encode_integer(&base.discriminant),
            algebraic_bound: base.algebraic_bound,
            large_prime_bound: base.large_prime_bound,
            sign_column: base.sign_column,
            lp_columns: base.lp_columns,
            algebraic_index,
        },
        complete: complete.iter().map(store_relation).collect(),
        bank: StoredBank {
            singles: bank.singles.iter().map(store_relation).collect(),
            single_sources: bank.single_sources.clone(),
            single_aux: bank.single_aux.clone(),
            single_lps: bank.single_lps.clone(),
            collector_pending: sorted_pairs(&bank.collector.pending),
            rat_pending: sorted_pairs(&bank.rat_pending),
            mixed_pending: sorted_pairs(&bank.mixed_pending),
            triad_pending: sorted_pairs(&bank.triad_pending),
            doubles: bank.doubles.iter().map(store_relation).collect(),
            double_sources: bank.double_sources.clone(),
            double_aux: bank.double_aux.clone(),
            double_parent: sorted_pairs(&bank.double_parent),
            double_members: sorted_pairs(&bank.double_members),
        },
        position: StoredPosition {
            attempt_a,
            attempt_b,
            passes_done,
            usable_total,
            finished,
            consumed: consumed_sorted,
        },
    }
}

fn restore_map<K: Eq + std::hash::Hash, V>(entries: Vec<(K, V)>) -> HashMap<K, V> {
    entries.into_iter().collect()
}

pub(crate) fn validate_and_restore(
    payload: &CheckpointPayload,
    n: &Integer,
    degree: usize,
    options: &FactorOptions,
    polynomial: &PolynomialPair,
    base: &GeneralBase,
) -> Result<CollectResume, FactorError> {
    if payload.n != encode_integer(n) {
        return Err(FactorError::Checkpoint(
            "checkpoint n mismatch: resume requires the same input".to_string(),
        ));
    }
    if payload.degree != degree {
        return Err(FactorError::Checkpoint(
            "checkpoint degree mismatch: resume requires the same degree".to_string(),
        ));
    }
    let stored_poly = &payload.polynomial;
    let live_coeffs: Vec<String> = polynomial.coefficients.iter().map(encode_integer).collect();
    if stored_poly.degree != polynomial.degree
        || stored_poly.coefficients != live_coeffs
        || stored_poly.m != encode_integer(&polynomial.m)
    {
        return Err(FactorError::Checkpoint(
            "checkpoint polynomial mismatch: resume requires the same polynomial".to_string(),
        ));
    }
    let stored_options = &payload.options;
    if stored_options.rational_base_bound != options.rational_base_bound
        || stored_options.algebraic_base_bound != options.algebraic_base_bound
        || stored_options.sieve_a_bound != options.sieve_a_bound
        || stored_options.sieve_b_bound != options.sieve_b_bound
        || stored_options.large_prime_bound != options.large_prime_bound
        || stored_options.seed != options.seed
        || stored_options.enable_double_large_primes != options.enable_double_large_primes
    {
        return Err(FactorError::Checkpoint(
            "checkpoint options mismatch: resume requires the same sieve bounds, large-prime bound, double-prime flag and seed".to_string(),
        ));
    }
    let stored_base = &payload.base;
    let mut live_index = sorted_pairs(&base.algebraic_index);
    for (_, roots) in live_index.iter_mut() {
        roots.sort();
    }
    if stored_base.rational_primes != base.rational_primes
        || stored_base.algebraic_roots != base.algebraic_roots
        || stored_base.modulus_primes != base.modulus_primes
        || stored_base.discriminant != encode_integer(&base.discriminant)
        || stored_base.algebraic_bound != base.algebraic_bound
        || stored_base.large_prime_bound != base.large_prime_bound
        || stored_base.sign_column != base.sign_column
        || stored_base.lp_columns != base.lp_columns
        || stored_base.algebraic_index != live_index
    {
        return Err(FactorError::Checkpoint(
            "checkpoint factor-base mismatch: resume requires the same factor base".to_string(),
        ));
    }
    let mut complete = Vec::with_capacity(payload.complete.len());
    for stored in &payload.complete {
        complete.push(restore_relation(stored)?);
    }
    let stored_bank = &payload.bank;
    let mut singles = Vec::with_capacity(stored_bank.singles.len());
    for stored in &stored_bank.singles {
        singles.push(restore_relation(stored)?);
    }
    let mut doubles = Vec::with_capacity(stored_bank.doubles.len());
    for stored in &stored_bank.doubles {
        doubles.push(restore_relation(stored)?);
    }
    let bank = PartialBank {
        singles,
        single_sources: stored_bank.single_sources.clone(),
        single_aux: stored_bank.single_aux.clone(),
        single_lps: stored_bank.single_lps.clone(),
        collector: crate::general::LargePrimeCollector {
            pending: restore_map(stored_bank.collector_pending.clone()),
        },
        rat_pending: restore_map(stored_bank.rat_pending.clone()),
        mixed_pending: restore_map(stored_bank.mixed_pending.clone()),
        triad_pending: restore_map(stored_bank.triad_pending.clone()),
        doubles,
        double_sources: stored_bank.double_sources.clone(),
        double_aux: stored_bank.double_aux.clone(),
        double_parent: restore_map(stored_bank.double_parent.clone()),
        double_members: restore_map(stored_bank.double_members.clone()),
    };
    if bank.singles.len() != bank.single_sources.len()
        || bank.singles.len() != bank.single_aux.len()
        || bank.singles.len() != bank.single_lps.len()
        || bank.doubles.len() != bank.double_sources.len()
        || bank.doubles.len() != bank.double_aux.len()
    {
        return Err(FactorError::Checkpoint(
            "checkpoint corrupt: partial-bank lengths disagree".to_string(),
        ));
    }
    let consumed: HashSet<(i64, i64)> = payload.position.consumed.iter().copied().collect();
    Ok(CollectResume {
        complete,
        bank,
        consumed,
        attempt_a: payload.position.attempt_a,
        attempt_b: payload.position.attempt_b,
        passes_done: payload.position.passes_done,
        usable_total: payload.position.usable_total,
        finished: payload.position.finished,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::general::PartialBank;
    use crate::number_field::Element;

    fn stored_relation_fixture(rational: i64, mask: Vec<u64>) -> StoredRelation {
        StoredRelation {
            rational_value: rational.to_string(),
            algebraic_value: vec![rational.to_string(), "1".to_string()],
            mask,
            large_prime: None,
            double_edge: None,
        }
    }

    fn payload_fixture() -> CheckpointPayload {
        CheckpointPayload {
            n: "100160063".to_string(),
            degree: 3,
            options: StoredOptions {
                rational_base_bound: 2000,
                algebraic_base_bound: 2000,
                sieve_a_bound: 700,
                sieve_b_bound: 120,
                large_prime_bound: 50000,
                seed: 7,
                enable_double_large_primes: false,
            },
            polynomial: StoredPolynomial {
                degree: 3,
                coefficients: vec![
                    "1".to_string(),
                    "0".to_string(),
                    "0".to_string(),
                    "2".to_string(),
                ],
                m: "464".to_string(),
            },
            base: StoredBase {
                rational_primes: vec![2, 3, 5],
                algebraic_roots: vec![(7, 1, 3)],
                modulus_primes: vec![7],
                discriminant: "-27".to_string(),
                algebraic_bound: 2000,
                large_prime_bound: 50000,
                sign_column: 4,
                lp_columns: 0,
                algebraic_index: vec![(7, vec![(1, 3)])],
            },
            complete: vec![stored_relation_fixture(11, vec![0b101])],
            bank: StoredBank {
                singles: vec![stored_relation_fixture(13, vec![0b011])],
                single_sources: vec![(1, 5)],
                single_aux: vec![1],
                single_lps: vec![vec![(101, 3)]],
                collector_pending: vec![((101, 3), 0)],
                rat_pending: Vec::new(),
                mixed_pending: Vec::new(),
                triad_pending: Vec::new(),
                doubles: Vec::new(),
                double_sources: Vec::new(),
                double_aux: Vec::new(),
                double_parent: Vec::new(),
                double_members: Vec::new(),
            },
            position: StoredPosition {
                attempt_a: 700,
                attempt_b: 120,
                passes_done: 1,
                usable_total: 1,
                finished: false,
                consumed: vec![(1, 5)],
            },
        }
    }

    fn temp_path(name: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "graynfs-checkpoint-{name}-{}-{}.json",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|span| span.as_nanos())
                .unwrap_or(0)
        ));
        path
    }

    #[test]
    fn checkpoint_roundtrip_preserves_collector_state() {
        let payload = payload_fixture();
        let path = temp_path("roundtrip");
        save_gnfs_checkpoint(&path, &payload).expect("save succeeds");
        let loaded = load_gnfs_checkpoint(&path).expect("load succeeds");
        let _ = std::fs::remove_file(&path);
        assert_eq!(loaded.n, payload.n);
        assert_eq!(loaded.degree, payload.degree);
        assert_eq!(loaded.complete.len(), 1);
        assert_eq!(loaded.complete[0].rational_value, "11");
        assert_eq!(loaded.complete[0].algebraic_value, vec!["11", "1"]);
        assert_eq!(loaded.bank.singles.len(), 1);
        assert_eq!(loaded.bank.single_sources, vec![(1, 5)]);
        assert_eq!(loaded.bank.single_aux, vec![1]);
        assert_eq!(loaded.bank.single_lps, vec![vec![(101, 3)]]);
        assert_eq!(loaded.bank.collector_pending, vec![((101, 3), 0)]);
        assert_eq!(loaded.position.attempt_a, 700);
        assert_eq!(loaded.position.attempt_b, 120);
        assert_eq!(loaded.position.passes_done, 1);
        assert_eq!(loaded.position.usable_total, 1);
        assert!(!loaded.position.finished);
        assert_eq!(loaded.position.consumed, vec![(1, 5)]);
        assert_eq!(loaded.base.rational_primes, vec![2, 3, 5]);
        assert_eq!(loaded.options.seed, 7);
    }

    #[test]
    fn checkpoint_checksum_tamper_is_refused() {
        let payload = payload_fixture();
        let path = temp_path("tamper");
        save_gnfs_checkpoint(&path, &payload).expect("save succeeds");
        let mut text = std::fs::read_to_string(&path).expect("read succeeds");
        let target = "\"rational_value\":\"11\"";
        let replacement = "\"rational_value\":\"12\"";
        assert!(text.contains(target));
        text = text.replace(target, replacement);
        std::fs::write(&path, text).expect("rewrite succeeds");
        let error = load_gnfs_checkpoint(&path).expect_err("tampered file is refused");
        let detail = error.to_string();
        let _ = std::fs::remove_file(&path);
        assert!(
            detail.contains("checksum"),
            "tamper refusal must name the checksum: {detail}"
        );
    }

    #[test]
    fn checkpoint_version_mismatch_is_refused() {
        let payload = payload_fixture();
        let path = temp_path("version");
        save_gnfs_checkpoint(&path, &payload).expect("save succeeds");
        let text = std::fs::read_to_string(&path).expect("read succeeds");
        assert!(text.contains("\"version\":1"));
        let text = text.replacen("\"version\":1", "\"version\":2", 1);
        std::fs::write(&path, text).expect("rewrite succeeds");
        let error = load_gnfs_checkpoint(&path).expect_err("version-mismatched file is refused");
        let detail = error.to_string();
        let _ = std::fs::remove_file(&path);
        assert!(
            detail.contains("version"),
            "version refusal must name the version: {detail}"
        );
    }

    #[test]
    fn checkpoint_truncated_file_is_refused() {
        let payload = payload_fixture();
        let path = temp_path("truncated");
        save_gnfs_checkpoint(&path, &payload).expect("save succeeds");
        let mut bytes = std::fs::read(&path).expect("read succeeds");
        bytes.truncate(bytes.len() / 2);
        std::fs::write(&path, bytes).expect("rewrite succeeds");
        let error = load_gnfs_checkpoint(&path).expect_err("truncated file is refused");
        let detail = error.to_string();
        let _ = std::fs::remove_file(&path);
        assert!(
            detail.contains("corrupt") || detail.contains("checksum"),
            "truncation refusal must name corruption: {detail}"
        );
    }

    #[test]
    fn checkpoint_state_export_covers_bank_and_coordinates() {
        let consumed: HashSet<(i64, i64)> = [(2_i64, 9_i64)].into_iter().collect();
        let bank = PartialBank::new();
        let polynomial = PolynomialPair {
            degree: 3,
            coefficients: vec![
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(2),
            ],
            m: Integer::from(464),
        };
        let base = GeneralBase {
            rational_primes: vec![2, 3],
            rational_indices: [(2_u64, 0_usize), (3, 1)].into_iter().collect(),
            algebraic_roots: Vec::new(),
            algebraic_index: HashMap::new(),
            modulus_primes: Vec::new(),
            discriminant: Integer::from(-27),
            algebraic_bound: 2000,
            large_prime_bound: 50000,
            large_prime_roots: HashMap::new(),
            sign_column: 2,
            lp_columns: 0,
        };
        let relation = GeneralRelation {
            rational_value: Integer::from(11),
            algebraic_value: Element::from_coefficients(vec![Integer::from(11), Integer::from(1)]),
            mask: vec![0b101],
            large_prime: Some((101, 3)),
            double_edge: None,
        };
        let options = FactorOptions {
            seed: 7,
            ..FactorOptions::default()
        };
        let n: Integer = "100160063".parse().expect("valid integer");
        let payload = payload_from_state(
            &n,
            3,
            &options,
            &polynomial,
            &base,
            std::slice::from_ref(&relation),
            &bank,
            &consumed,
            700,
            120,
            1,
            1,
            false,
        );
        assert_eq!(payload.position.consumed, vec![(2, 9)]);
        assert_eq!(payload.complete.len(), 1);
        assert_eq!(payload.complete[0].large_prime, Some((101, 3)));
        assert_eq!(payload.bank.singles.len(), 0);
        assert_eq!(payload.base.sign_column, 2);
    }
}
