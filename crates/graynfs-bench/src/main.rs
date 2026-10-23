use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use clap::Parser;
use graynfs_core::prime::{
    AKS_PROOF_MAX_BITS, is_baillie_psw_prime, is_frobenius_probable_prime, is_prime_aks_proof,
    is_probable_prime, is_probable_prime_with_options,
};
use graynfs_core::{FactorOptions, Integer, factor_with_progress, probable_prime};
use num_traits::Num;
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug)]
#[command(
    name = "graynfs-bench",
    version,
    about = "Benchmark GrayNFS across factorization difficulty tiers; repro267 archives belong under logs/<task-id>/"
)]
struct Cli {
    #[arg(long, default_value = "full")]
    tier: String,
    #[arg(long, default_value_t = 30)]
    timeout_seconds: u64,
    #[arg(long)]
    threads: Option<usize>,
    #[arg(long, default_value = "adhoc")]
    task_id: String,
    #[arg(long, default_value = "logs/graynfs-benchmark.json")]
    output: PathBuf,
    #[arg(long)]
    case_filter: Option<String>,
    #[arg(long, default_value_t = 40)]
    random_count: usize,
    #[arg(long, default_value_t = 0x9e37_79b9_7f4a_7c15)]
    random_seed: u64,
    #[arg(long, default_value_t = 16)]
    random_min_factor_bits: usize,
    #[arg(long, default_value_t = 24)]
    random_max_factor_bits: usize,
    #[arg(long)]
    compare_baseline: Option<PathBuf>,
    #[arg(long)]
    compare_candidate: Option<PathBuf>,
}

#[derive(Clone, Debug)]
struct Case {
    name: String,
    value: Integer,
    expected: Vec<Integer>,
    bits: usize,
    repeats: usize,
    wide: bool,
    deadline_secs: Option<u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct StageSample {
    stage: String,
    current: u64,
    total: u64,
    message: String,
    start_ms: f64,
    end_ms: f64,
    duration_ms: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Observation {
    name: String,
    input: String,
    bits: usize,
    method: String,
    factors: Vec<String>,
    expected: Vec<String>,
    verified: bool,
    total_ms: f64,
    status: String,
    error: Option<String>,
    stages: Vec<StageSample>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct BenchArchive {
    task_id: String,
    run_id: String,
    commit_or_id: String,
    commit_sha: Option<String>,
    tier: String,
    timeout_seconds: u64,
    threads_requested: Option<usize>,
    threads_effective: usize,
    repeats_total: usize,
    case_filter: Option<String>,
    random_count: usize,
    random_seed: u64,
    random_min_factor_bits: usize,
    random_max_factor_bits: usize,
    timestamp_unix_secs: u64,
    observations: Vec<Observation>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    if cli.compare_baseline.is_some() || cli.compare_candidate.is_some() {
        return run_compare(&cli);
    }
    let threads_effective = cli.threads.unwrap_or(
        std::thread::available_parallelism()
            .map(|count| count.get())
            .unwrap_or(1),
    );
    let timestamp_unix_secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_secs())
        .unwrap_or(0);
    let mut cases = if cli.tier == "random" {
        random_cases(&cli)?
    } else {
        build_cases(&cli.tier)?
    };
    if let Some(filter) = cli.case_filter.as_deref() {
        cases.retain(|case| case.name.contains(filter));
        if cases.is_empty() {
            return Err(format!("unknown case: {filter}").into());
        }
    }
    let mut observations = Vec::new();
    for case in &cases {
        for _ in 0..case.repeats {
            observations.push(run_case(case, cli.timeout_seconds, threads_effective));
        }
    }
    print_summary(&observations);
    if let Some(parent) = cli.output.parent() {
        fs::create_dir_all(parent)?;
    }
    let archive = BenchArchive {
        task_id: cli.task_id.clone(),
        run_id: build_run_id(&cli.task_id),
        commit_or_id: cli.task_id.clone(),
        commit_sha: resolve_commit_sha(),
        tier: cli.tier.clone(),
        timeout_seconds: cli.timeout_seconds,
        threads_requested: cli.threads,
        threads_effective,
        repeats_total: observations.len(),
        case_filter: cli.case_filter.clone(),
        random_count: cli.random_count,
        random_seed: cli.random_seed,
        random_min_factor_bits: cli.random_min_factor_bits,
        random_max_factor_bits: cli.random_max_factor_bits,
        timestamp_unix_secs,
        observations,
    };
    fs::write(&cli.output, serde_json::to_vec_pretty(&archive)?)?;
    println!("json: {}", cli.output.display());
    Ok(())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ParityProvenance {
    run_id: String,
    commit_or_id: String,
    commit_sha: Option<String>,
    observation_count: usize,
    verified_count: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ParityMember {
    input: String,
    expected: Vec<String>,
    baseline_factors: Vec<String>,
    candidate_factors: Vec<String>,
    baseline_total_ms: f64,
    candidate_total_ms: f64,
    baseline_runs: usize,
    candidate_runs: usize,
    baseline_verified: bool,
    candidate_verified: bool,
    status: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ParityEnvelope {
    task_id: String,
    run_id: String,
    commit_or_id: String,
    commit_sha: Option<String>,
    baseline: ParityProvenance,
    candidate: ParityProvenance,
    members: Vec<ParityMember>,
    #[serde(default)]
    primality_members: Vec<PrimalityMember>,
    #[serde(default)]
    primality_status: String,
    status: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct PrimalityMember {
    input: String,
    expected_is_prime: bool,
    probable: bool,
    bpsw: bool,
    frobenius_gated: bool,
    aks: Option<bool>,
    status: String,
    run_id: String,
    commit_or_id: String,
    commit_sha: Option<String>,
}

fn primality_corpus() -> Vec<(String, bool)> {
    let mut corpus = Vec::new();
    for digits in [
        "561",
        "1105",
        "1729",
        "2465",
        "2821",
        "6601",
        "8911",
        "10585",
        "15841",
        "29341",
        "41041",
        "46657",
        "52633",
        "62745",
        "63973",
        "75361",
        "2047",
        "3215031751",
        "2152302898747",
        "3474749660383",
        "341550071728321",
        "3825123056546413051",
    ] {
        corpus.push((digits.to_string(), false));
    }
    for digits in [
        "2",
        "3",
        "5",
        "7",
        "11",
        "13",
        "17",
        "19",
        "23",
        "29",
        "31",
        "37",
        "41",
        "43",
        "47",
        "101",
        "997",
        "1000003",
        "2147483647",
        "170141183460469231731687303715884105727",
    ] {
        corpus.push((digits.to_string(), true));
    }
    corpus
}

fn aks_composite_feasible(value: &Integer) -> bool {
    value.bits() <= 43
}

fn aks_prime_feasible(value: &Integer) -> bool {
    value.bits() <= 17
}

fn primality_member_status(
    expected_is_prime: bool,
    probable: bool,
    bpsw: bool,
    frobenius_gated: bool,
    aks: Option<bool>,
) -> &'static str {
    let mut ok = probable == expected_is_prime
        && bpsw == expected_is_prime
        && frobenius_gated == expected_is_prime;
    if let Some(proof) = aks {
        ok = ok && proof == expected_is_prime;
    }
    status_for(ok)
}

fn build_primality_members(
    run_id: &str,
    commit_or_id: &str,
    commit_sha: Option<String>,
) -> Vec<PrimalityMember> {
    let gated = FactorOptions {
        use_frobenius: true,
        ..FactorOptions::default()
    };
    let mut members = Vec::new();
    for (input, expected_is_prime) in primality_corpus() {
        let value = Integer::from_str_radix(&input, 10).expect("primality input parses");
        let probable = is_probable_prime(&value);
        let bpsw = is_baillie_psw_prime(&value);
        let frobenius_gated = is_probable_prime_with_options(&value, &gated);
        let frobenius_direct = is_frobenius_probable_prime(&value);
        let combined_gated = frobenius_gated && frobenius_direct;
        let aks = if value.bits() > AKS_PROOF_MAX_BITS {
            None
        } else if expected_is_prime {
            if aks_prime_feasible(&value) {
                Some(is_prime_aks_proof(&value))
            } else {
                None
            }
        } else if aks_composite_feasible(&value) {
            Some(is_prime_aks_proof(&value))
        } else {
            None
        };
        let status =
            primality_member_status(expected_is_prime, probable, bpsw, combined_gated, aks)
                .to_string();
        members.push(PrimalityMember {
            input,
            expected_is_prime,
            probable,
            bpsw,
            frobenius_gated: combined_gated,
            aks,
            status,
            run_id: run_id.to_string(),
            commit_or_id: commit_or_id.to_string(),
            commit_sha: commit_sha.clone(),
        });
    }
    members
}

fn parity_corpus() -> Vec<(String, Vec<String>)> {
    vec![
        (
            "770817968652704265841785071059".to_string(),
            vec!["825459858613693".to_string(), "933804303878863".to_string()],
        ),
        (
            "3426790303412814276981061".to_string(),
            vec!["1771204694653".to_string(), "1934722911337".to_string()],
        ),
        (
            "1155285037654380637".to_string(),
            vec!["1074741869".to_string(), "1074941873".to_string()],
        ),
        (
            "166871515843890537614148105919789085907512212182703922830333764223351416812060867"
                .to_string(),
            vec![
                "121680281986369904489".to_string(),
                "1371393237423486234387340775961402466614689310566958802486603".to_string(),
            ],
        ),
    ]
}

fn parity_cases() -> Vec<Case> {
    let bits = [100_usize, 82, 61, 267];
    let names = [
        "parity_100bit",
        "parity_82bit",
        "parity_deg3skip",
        "parity_demo267",
    ];
    parity_corpus()
        .into_iter()
        .zip(bits)
        .zip(names)
        .map(|(((input, expected), bit_count), name)| {
            let value = Integer::from_str_radix(&input, 10).expect("parity input parses");
            let factors = expected
                .iter()
                .map(|factor| Integer::from_str_radix(factor, 10).expect("parity factor parses"))
                .collect();
            Case {
                name: name.to_string(),
                value,
                expected: factors,
                bits: bit_count,
                repeats: 5,
                wide: false,
                deadline_secs: None,
            }
        })
        .collect()
}

fn provenance_of(archive: &BenchArchive) -> ParityProvenance {
    ParityProvenance {
        run_id: archive.run_id.clone(),
        commit_or_id: archive.commit_or_id.clone(),
        commit_sha: archive.commit_sha.clone(),
        observation_count: archive.observations.len(),
        verified_count: archive
            .observations
            .iter()
            .filter(|observation| observation.verified)
            .count(),
    }
}

fn observations_for<'a>(archive: &'a BenchArchive, input: &str) -> Vec<&'a Observation> {
    archive
        .observations
        .iter()
        .filter(|observation| observation.input == input)
        .collect()
}

fn median_of(mut values: Vec<f64>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|first, second| first.total_cmp(second));
    Some(values[values.len() / 2])
}

fn compare_parity_envelopes(baseline: &BenchArchive, candidate: &BenchArchive) -> ParityEnvelope {
    let mut members = Vec::new();
    for (input, expected) in parity_corpus() {
        let befores = observations_for(baseline, &input);
        let afters = observations_for(candidate, &input);
        let mut member = ParityMember {
            input: input.clone(),
            expected: sorted_strings(&expected),
            baseline_factors: Vec::new(),
            candidate_factors: Vec::new(),
            baseline_total_ms: f64::NAN,
            candidate_total_ms: f64::NAN,
            baseline_runs: befores.len(),
            candidate_runs: afters.len(),
            baseline_verified: false,
            candidate_verified: false,
            status: status_for(false).to_string(),
        };
        if !befores.is_empty() && !afters.is_empty() {
            let before_totals: Vec<f64> = befores
                .iter()
                .map(|observation| observation.total_ms)
                .collect();
            let after_totals: Vec<f64> = afters
                .iter()
                .map(|observation| observation.total_ms)
                .collect();
            member.baseline_total_ms = median_of(before_totals).unwrap_or(f64::NAN);
            member.candidate_total_ms = median_of(after_totals).unwrap_or(f64::NAN);
            member.baseline_factors = sorted_strings(&befores[0].factors);
            member.candidate_factors = sorted_strings(&afters[0].factors);
            member.baseline_verified = befores.iter().all(|observation| observation.verified);
            member.candidate_verified = afters.iter().all(|observation| observation.verified);
            let seen: Vec<Vec<String>> = befores
                .iter()
                .chain(afters.iter())
                .map(|observation| sorted_strings(&observation.factors))
                .collect();
            let factors_identical = seen.iter().all(|factors| *factors == seen[0]);
            let product_holds = member.baseline_verified
                && member.candidate_verified
                && seen[0] == sorted_strings(&expected);
            let total_within_old = member.candidate_total_ms <= member.baseline_total_ms;
            member.status =
                status_for(factors_identical && product_holds && total_within_old).to_string();
        }
        members.push(member);
    }
    let parity_status = members.iter().all(|member| member.status == "ok");
    let run_id = build_run_id("v0.1.5.2.7");
    let commit_or_id = "v0.1.5.2.7".to_string();
    let commit_sha = resolve_commit_sha();
    let primality_members = build_primality_members(&run_id, &commit_or_id, commit_sha.clone());
    let primality_status = primality_members.iter().all(|member| member.status == "ok");
    let primality_report = status_for(primality_status).to_string();
    let status = status_for(parity_status && primality_status).to_string();
    ParityEnvelope {
        task_id: "v0.1.5.2.7".to_string(),
        run_id,
        commit_or_id,
        commit_sha,
        baseline: provenance_of(baseline),
        candidate: provenance_of(candidate),
        members,
        primality_members,
        primality_status: primality_report,
        status,
    }
}

fn run_compare(cli: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let baseline_path = cli
        .compare_baseline
        .as_ref()
        .ok_or("compare_baseline is required with compare_candidate")?;
    let candidate_path = cli
        .compare_candidate
        .as_ref()
        .ok_or("compare_candidate is required with compare_baseline")?;
    let baseline: BenchArchive = serde_json::from_str(&fs::read_to_string(baseline_path)?)?;
    let candidate: BenchArchive = serde_json::from_str(&fs::read_to_string(candidate_path)?)?;
    let envelope = compare_parity_envelopes(&baseline, &candidate);
    println!(
        "input,expected,baseline_factors,candidate_factors,baseline_ms,candidate_ms,baseline_verified,candidate_verified,status"
    );
    for member in &envelope.members {
        println!(
            "{},{},{},{},{:.3},{:.3},{},{},{}",
            member.input,
            member.expected.join("+"),
            member.baseline_factors.join("+"),
            member.candidate_factors.join("+"),
            member.baseline_total_ms,
            member.candidate_total_ms,
            member.baseline_verified,
            member.candidate_verified,
            member.status,
        );
    }
    println!(
        "baseline_run={},candidate_run={},status={},primality_status={}",
        envelope.baseline.run_id,
        envelope.candidate.run_id,
        envelope.status,
        envelope.primality_status,
    );
    println!(
        "input,expected_is_prime,probable,bpsw,frobenius_gated,aks,status,run_id,commit_or_id"
    );
    for member in &envelope.primality_members {
        println!(
            "{},{},{},{},{},{},{},{},{}",
            member.input,
            member.expected_is_prime,
            member.probable,
            member.bpsw,
            member.frobenius_gated,
            member
                .aks
                .map(|proof| proof.to_string())
                .unwrap_or("deferred".to_string()),
            member.status,
            member.run_id,
            member.commit_or_id,
        );
    }
    if let Some(parent) = cli.output.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&cli.output, serde_json::to_vec_pretty(&envelope)?)?;
    println!("json: {}", cli.output.display());
    if envelope.status != "ok" {
        return Err("parity envelope verification_failed".into());
    }
    Ok(())
}

static RUN_ID_SEQ: AtomicU64 = AtomicU64::new(0);

fn build_run_id(task_id: &str) -> String {
    let nanos_since_epoch = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_nanos())
        .unwrap_or(0);
    let seq = RUN_ID_SEQ.fetch_add(1, Ordering::Relaxed);
    format!(
        "{}-{}-{}-{}",
        task_id,
        nanos_since_epoch,
        std::process::id(),
        seq
    )
}

fn status_for(verified: bool) -> &'static str {
    if verified {
        "ok"
    } else {
        "verification_failed"
    }
}

fn resolve_commit_sha() -> Option<String> {
    for key in ["GITHUB_SHA", "GRAYNFS_COMMIT"] {
        match std::env::var(key) {
            Ok(value) if !value.trim().is_empty() => return Some(value.trim().to_string()),
            _ => continue,
        }
    }
    None
}

fn run_case(case: &Case, timeout_seconds: u64, threads: usize) -> Observation {
    let options = if case.wide {
        let scale = wide_scale(case.bits);
        FactorOptions {
            rational_base_bound: 2_000 * scale,
            algebraic_base_bound: 2_000 * scale,
            sieve_a_bound: 700 * scale,
            sieve_b_bound: 120 * scale,
            threads: Some(threads),
            max_duration: Some(Duration::from_secs(
                case.deadline_secs.unwrap_or(timeout_seconds),
            )),
            allow_fallback: true,
            ..FactorOptions::default()
        }
    } else {
        FactorOptions {
            threads: Some(threads),
            max_duration: Some(Duration::from_secs(
                case.deadline_secs.unwrap_or(timeout_seconds),
            )),
            allow_fallback: true,
            ..FactorOptions::default()
        }
    };
    let started = Instant::now();
    let mut raw_events = Vec::new();
    let result = factor_with_progress(&case.value, &options, |event| {
        raw_events.push((
            started.elapsed().as_secs_f64() * 1_000.0,
            event.stage.clone(),
            event.current,
            event.total,
            event.message.clone(),
        ));
        Ok(())
    });
    let total_ms = started.elapsed().as_secs_f64() * 1_000.0;
    let stages = finish_stages(raw_events, total_ms);
    let mut observation = Observation {
        name: case.name.clone(),
        input: case.value.to_string(),
        bits: case.bits,
        method: String::new(),
        factors: Vec::new(),
        expected: case.expected.iter().map(ToString::to_string).collect(),
        verified: false,
        total_ms,
        status: "ok".to_string(),
        error: None,
        stages,
    };
    match result {
        Ok(report) => {
            observation.method = report.method;
            observation.factors = report.factors;
            let product =
                parse_decimal_product(&observation.factors).unwrap_or_else(|_| Integer::from(0));
            observation.verified = product == case.value
                && sorted_strings(&observation.factors) == sorted_strings(&observation.expected);
            observation.status = status_for(observation.verified).to_string();
        }
        Err(error) => {
            observation.status = "error".to_string();
            observation.error = Some(error.to_string());
        }
    }
    observation
}

fn finish_stages(
    raw_events: Vec<(f64, String, u64, u64, String)>,
    total_ms: f64,
) -> Vec<StageSample> {
    let mut stages = Vec::with_capacity(raw_events.len());
    for (index, (start_ms, stage, current, total, message)) in raw_events.iter().enumerate() {
        let end_ms = raw_events
            .get(index + 1)
            .map(|event| event.0)
            .unwrap_or(total_ms);
        stages.push(StageSample {
            stage: stage.clone(),
            current: *current,
            total: *total,
            message: message.clone(),
            start_ms: *start_ms,
            end_ms,
            duration_ms: (end_ms - start_ms).max(0.0),
        });
    }
    stages
}

fn print_summary(observations: &[Observation]) {
    println!(
        "case,bits,status,verified,method,total_ms,preprocessing_ms,polynomial_ms,factor_base_ms,collect_ms,relations_ms,rank_ms,dependencies_ms,square_root_ms"
    );
    for observation in observations {
        println!(
            "{},{},{},{},{},{:.3},{},{},{},{},{},{},{},{}",
            observation.name,
            observation.bits,
            observation.status,
            observation.verified,
            observation.method,
            observation.total_ms,
            duration_for(&observation.stages, "preprocessing"),
            duration_for(&observation.stages, "polynomial"),
            duration_for(&observation.stages, "base"),
            duration_for(&observation.stages, "collect"),
            duration_for(&observation.stages, "relations"),
            duration_for(&observation.stages, "rank"),
            duration_for(&observation.stages, "dependencies"),
            duration_for(&observation.stages, "square_root"),
        );
    }
}

fn duration_for(stages: &[StageSample], kind: &str) -> String {
    let total = stages
        .iter()
        .filter(|stage| match kind {
            "preprocessing" => stage.stage == "preprocessing",
            "polynomial" => stage.message == "polynomial selected",
            "base" => stage.message == "factor base built",
            "collect" => stage.message == "sieving the quadratic base",
            "relations" => stage.message == "relations collected",
            "rank" => stage.message == "matrix rank computed",
            "dependencies" => stage.message == "relation dependencies found",
            "square_root" => stage.message == "verified a congruence of squares",
            _ => false,
        })
        .map(|stage| stage.duration_ms)
        .sum::<f64>();
    if total == 0.0 {
        "-".to_string()
    } else {
        format!("{total:.3}")
    }
}

fn wide_scale(bits: usize) -> usize {
    bits.div_ceil(16)
}

fn build_cases(tier: &str) -> Result<Vec<Case>, Box<dyn std::error::Error>> {
    if tier == "repro267" {
        return Ok(vec![repro_267_case()]);
    }
    if tier == "parity" {
        return Ok(parity_cases());
    }
    if !matches!(tier, "quick" | "full" | "stress") {
        return Err(format!("unknown tier: {tier}").into());
    }
    let mut cases = vec![
        power_case("pow2_20", 2, 20, 3),
        power_case("pow2_40", 2, 40, 3),
        mixed_case("pow3_15_5_10", &[(3, 15), (5, 10)], 3),
        power_case("square_997_4", 997, 4, 3),
        mixed_case(
            "smooth_transition",
            &[(2, 8), (3, 6), (5, 3), (10007, 3)],
            3,
        ),
    ];
    if tier == "quick" {
        cases.truncate(4);
        return Ok(cases);
    }
    cases.push(semiprime_case("balanced_32", 32, 100, 1, true));
    cases.push(semiprime_case("balanced_40", 40, 100, 1, true));
    cases.push(semiprime_case("balanced_48", 48, 100, 1, true));
    cases.push(semiprime_case("balanced_64", 64, 100, 1, true));
    if tier == "stress" {
        cases.push(semiprime_case("stress_80", 80, 100, 1, true));
        cases.push(semiprime_case("stress_96", 96, 100, 1, true));
        cases.push(semiprime_case("stress_128", 128, 100, 1, true));
        cases.push(repro_267_case());
    }
    Ok(cases)
}

fn power_case(name: &str, base: u64, exponent: u32, repeats: usize) -> Case {
    let value = Integer::from(base).pow(exponent);
    let expected = (0..exponent)
        .map(|_| Integer::from(base))
        .collect::<Vec<_>>();
    Case {
        name: name.to_string(),
        bits: (u64::from(exponent) * u64::from(base.ilog2() + 1)).max(1) as usize,
        value,
        expected,
        repeats,
        wide: false,
        deadline_secs: None,
    }
}

fn mixed_case(name: &str, factors: &[(u64, u32)], repeats: usize) -> Case {
    let mut value = Integer::from(1);
    let mut expected = Vec::new();
    for (base, exponent) in factors {
        value *= Integer::from(*base).pow(*exponent);
        for _ in 0..*exponent {
            expected.push(Integer::from(*base));
        }
    }
    Case {
        name: name.to_string(),
        bits: value.to_string().len() * 4,
        value,
        expected,
        repeats,
        wide: false,
        deadline_secs: None,
    }
}

fn semiprime_case(name: &str, bits: usize, offset: u64, repeats: usize, wide: bool) -> Case {
    let first = next_prime(bits, offset);
    let second = next_prime_above(&(&first + Integer::from(2)));
    Case {
        name: name.to_string(),
        bits: bits * 2,
        value: &first * &second,
        expected: vec![first, second],
        repeats,
        wide,
        deadline_secs: None,
    }
}

fn repro_267_case() -> Case {
    let value = Integer::from_str_radix(
        "166871515843890537614148105919789085907512212182703922830333764223351416812060867",
        10,
    )
    .expect("repro267 value parses");
    let first =
        Integer::from_str_radix("121680281986369904489", 10).expect("repro267 first factor parses");
    let second = Integer::from_str_radix(
        "1371393237423486234387340775961402466614689310566958802486603",
        10,
    )
    .expect("repro267 second factor parses");
    Case {
        name: "repro_267".to_string(),
        bits: 267,
        value,
        expected: vec![first, second],
        repeats: 1,
        wide: false,
        deadline_secs: Some(60),
    }
}

fn next_prime(bits: usize, offset: u64) -> Integer {
    next_prime_above(&((Integer::from(1) << (bits - 1)) | Integer::from(offset)))
}

fn next_prime_above(lower: &Integer) -> Integer {
    let mut candidate = lower | Integer::from(1);
    loop {
        if probable_prime(&candidate) {
            return candidate;
        }
        candidate += Integer::from(2);
    }
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut mixed = *state;
    mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    mixed ^ (mixed >> 31)
}

fn random_semiprime(name: &str, first_bits: usize, second_bits: usize, state: &mut u64) -> Case {
    let first_span = first_bits.saturating_sub(1).min(40);
    let first_mask = (1_u64 << first_span) - 1;
    let first = next_prime(first_bits, splitmix64(state) & first_mask);
    let floor = Integer::from(1) << second_bits.saturating_sub(1);
    let above = &first + Integer::from(2);
    let start = if above >= floor { above } else { floor };
    let second = next_prime_above(&start);
    Case {
        name: name.to_string(),
        bits: first_bits + second_bits,
        value: &first * &second,
        expected: vec![first, second],
        repeats: 1,
        wide: false,
        deadline_secs: None,
    }
}

fn random_cases(cli: &Cli) -> Result<Vec<Case>, Box<dyn std::error::Error>> {
    let low = cli.random_min_factor_bits.max(8);
    let high = cli.random_max_factor_bits.max(low);
    if cli.random_count == 0 {
        return Ok(Vec::new());
    }
    let mut state = cli.random_seed;
    let mut cases = Vec::with_capacity(cli.random_count);
    for index in 0..cli.random_count {
        let span = high - low + 1;
        let first_bits = low + (splitmix64(&mut state) as usize) % span;
        let second_bits = low + (splitmix64(&mut state) as usize) % span;
        cases.push(random_semiprime(
            &format!("random_{index:03}_{first_bits}x{second_bits}"),
            first_bits,
            second_bits,
            &mut state,
        ));
    }
    Ok(cases)
}

fn parse_decimal_product(values: &[String]) -> Result<Integer, String> {
    values.iter().try_fold(Integer::from(1), |acc, value| {
        Ok(acc * Integer::from_str_radix(value, 10).map_err(|error| error.to_string())?)
    })
}

fn sorted_strings(values: &[String]) -> Vec<String> {
    let mut result = values.to_vec();
    result.sort();
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    static ENV_LOCK: Mutex<()> = Mutex::new(());

    fn saved_env(key: &str) -> Option<String> {
        std::env::var(key).ok()
    }

    fn restore_env(key: &str, saved: Option<String>) {
        unsafe {
            match saved {
                Some(value) => std::env::set_var(key, value),
                None => std::env::remove_var(key),
            }
        }
    }

    #[test]
    fn run_id_uniqueness_format() {
        let first = build_run_id("taskx");
        let second = build_run_id("taskx");
        assert_ne!(first, second);
        let mut seqs = Vec::new();
        for id in [&first, &second] {
            assert!(id.starts_with("taskx-"));
            let rest = &id["taskx-".len()..];
            let parts: Vec<&str> = rest.split('-').collect();
            assert_eq!(parts.len(), 3);
            assert!(parts[0].parse::<u128>().is_ok());
            assert_eq!(parts[1], std::process::id().to_string());
            let seq = parts[2].parse::<u64>().expect("seq parses");
            seqs.push(seq);
        }
        assert_eq!(seqs[1], seqs[0] + 1);
    }

    #[test]
    fn verification_failed_status() {
        assert_eq!(status_for(true), "ok");
        assert_eq!(status_for(false), "verification_failed");
        let factors = vec!["6".to_string()];
        let expected = vec!["7".to_string()];
        let product = parse_decimal_product(&factors).expect("product parses");
        let verified =
            product == Integer::from(42) && sorted_strings(&factors) == sorted_strings(&expected);
        assert!(!verified);
        assert_eq!(status_for(verified), "verification_failed");
    }

    #[test]
    fn commit_sha_injected_from_env() {
        let _guard = ENV_LOCK.lock().expect("env lock poisons");
        let saved_github = saved_env("GITHUB_SHA");
        let saved_gray = saved_env("GRAYNFS_COMMIT");
        unsafe {
            std::env::remove_var("GITHUB_SHA");
        }
        unsafe {
            std::env::set_var("GRAYNFS_COMMIT", "abc123def456");
        }
        let resolved = resolve_commit_sha();
        assert_eq!(resolved, Some("abc123def456".to_string()));
        restore_env("GITHUB_SHA", saved_github);
        restore_env("GRAYNFS_COMMIT", saved_gray);
    }

    #[test]
    fn commit_sha_absent_without_env() {
        let _guard = ENV_LOCK.lock().expect("env lock poisons");
        let saved_github = saved_env("GITHUB_SHA");
        let saved_gray = saved_env("GRAYNFS_COMMIT");
        unsafe {
            std::env::remove_var("GITHUB_SHA");
            std::env::remove_var("GRAYNFS_COMMIT");
        }
        let resolved = resolve_commit_sha();
        assert_eq!(resolved, None);
        restore_env("GITHUB_SHA", saved_github);
        restore_env("GRAYNFS_COMMIT", saved_gray);
    }

    #[test]
    fn commit_sha_trims_whitespace_and_prefers_github() {
        let _guard = ENV_LOCK.lock().expect("env lock poisons");
        let saved_github = saved_env("GITHUB_SHA");
        let saved_gray = saved_env("GRAYNFS_COMMIT");
        unsafe {
            std::env::set_var("GITHUB_SHA", "  abc123  ");
            std::env::set_var("GRAYNFS_COMMIT", "  xyz  ");
        }
        let resolved = resolve_commit_sha();
        assert_eq!(resolved, Some("abc123".to_string()));
        restore_env("GITHUB_SHA", saved_github);
        restore_env("GRAYNFS_COMMIT", saved_gray);
    }

    #[test]
    fn commit_sha_skips_empty_github_sha() {
        let _guard = ENV_LOCK.lock().expect("env lock poisons");
        let saved_github = saved_env("GITHUB_SHA");
        let saved_gray = saved_env("GRAYNFS_COMMIT");
        unsafe {
            std::env::set_var("GITHUB_SHA", "");
            std::env::set_var("GRAYNFS_COMMIT", "  xyz  ");
        }
        let resolved = resolve_commit_sha();
        assert_eq!(resolved, Some("xyz".to_string()));
        restore_env("GITHUB_SHA", saved_github);
        restore_env("GRAYNFS_COMMIT", saved_gray);
    }

    fn stub_observation(
        input: &str,
        factors: &[&str],
        verified: bool,
        total_ms: f64,
    ) -> Observation {
        Observation {
            name: format!("stub_{input}"),
            input: input.to_string(),
            bits: 64,
            method: "stub".to_string(),
            factors: factors.iter().map(|value| value.to_string()).collect(),
            expected: factors.iter().map(|value| value.to_string()).collect(),
            verified,
            total_ms,
            status: status_for(verified).to_string(),
            error: None,
            stages: Vec::new(),
        }
    }

    fn stub_archive(run_id: &str, entries: Vec<(&str, Vec<&str>, bool, f64)>) -> BenchArchive {
        BenchArchive {
            task_id: "v0.1.5.1.7".to_string(),
            run_id: run_id.to_string(),
            commit_or_id: "v0.1.5.1.7".to_string(),
            commit_sha: None,
            tier: "parity".to_string(),
            timeout_seconds: 120,
            threads_requested: None,
            threads_effective: 1,
            repeats_total: entries.len(),
            case_filter: None,
            random_count: 0,
            random_seed: 0,
            random_min_factor_bits: 0,
            random_max_factor_bits: 0,
            timestamp_unix_secs: 0,
            observations: entries
                .into_iter()
                .map(|(input, factors, verified, total_ms)| {
                    stub_observation(input, &factors, verified, total_ms)
                })
                .collect(),
        }
    }

    fn stub_pair(
        candidate_input: &str,
        candidate_factors: Vec<&str>,
        candidate_verified: bool,
        candidate_total_ms: f64,
    ) -> (BenchArchive, BenchArchive) {
        let corpus = parity_corpus();
        let baseline_entries = corpus
            .iter()
            .map(|(input, expected)| {
                (
                    input.as_str(),
                    expected.iter().map(String::as_str).collect(),
                    true,
                    100.0,
                )
            })
            .collect();
        let candidate_entries = corpus
            .iter()
            .map(|(input, expected)| {
                if input == candidate_input {
                    (
                        input.as_str(),
                        candidate_factors.clone(),
                        candidate_verified,
                        candidate_total_ms,
                    )
                } else {
                    (
                        input.as_str(),
                        expected.iter().map(String::as_str).collect(),
                        true,
                        50.0,
                    )
                }
            })
            .collect();
        (
            stub_archive("baseline-run", baseline_entries),
            stub_archive("candidate-run", candidate_entries),
        )
    }

    #[test]
    fn parity_envelope_rejects_mismatched_factors() {
        let corpus = parity_corpus();
        let target = corpus[0].0.clone();
        let (baseline, candidate) = stub_pair(&target, vec!["11", "17"], true, 50.0);
        let envelope = compare_parity_envelopes(&baseline, &candidate);
        let member = envelope
            .members
            .iter()
            .find(|member| member.input == target)
            .expect("target member present");
        assert_eq!(member.status, "verification_failed");
        assert_eq!(envelope.status, "verification_failed");
    }

    #[test]
    fn parity_envelope_rejects_slower_candidate_total() {
        let corpus = parity_corpus();
        let target = corpus[1].0.clone();
        let expected: Vec<String> = corpus[1].1.clone();
        let borrowed: Vec<&str> = expected.iter().map(String::as_str).collect();
        let (baseline, candidate) = stub_pair(&target, borrowed, true, 200.0);
        let envelope = compare_parity_envelopes(&baseline, &candidate);
        let member = envelope
            .members
            .iter()
            .find(|member| member.input == target)
            .expect("target member present");
        assert_eq!(member.status, "verification_failed");
        assert_eq!(envelope.status, "verification_failed");
    }

    #[test]
    fn parity_envelope_accepts_identical_faster_factors() {
        let corpus = parity_corpus();
        let baseline_entries = corpus
            .iter()
            .map(|(input, expected)| {
                (
                    input.as_str(),
                    expected.iter().map(String::as_str).collect(),
                    true,
                    100.0,
                )
            })
            .collect();
        let candidate_entries = corpus
            .iter()
            .map(|(input, expected)| {
                (
                    input.as_str(),
                    expected.iter().map(String::as_str).collect(),
                    true,
                    50.0,
                )
            })
            .collect();
        let baseline = stub_archive("baseline-run", baseline_entries);
        let candidate = stub_archive("candidate-run", candidate_entries);
        let envelope = compare_parity_envelopes(&baseline, &candidate);
        assert_eq!(envelope.status, "ok");
        assert!(envelope.members.iter().all(|member| member.status == "ok"));
        assert_eq!(envelope.baseline.run_id, "baseline-run");
        assert_eq!(envelope.candidate.run_id, "candidate-run");
    }

    #[test]
    fn parity_envelope_compares_median_totals() {
        let corpus = parity_corpus();
        let target = corpus[0].0.clone();
        let expected: Vec<String> = corpus[0].1.clone();
        let borrowed: Vec<&str> = expected.iter().map(String::as_str).collect();
        let mut baseline_entries = Vec::new();
        let mut candidate_entries = Vec::new();
        for (input, factors) in &corpus {
            if *input == target {
                continue;
            }
            baseline_entries.push((
                input.as_str(),
                factors.iter().map(String::as_str).collect(),
                true,
                100.0,
            ));
            candidate_entries.push((
                input.as_str(),
                factors.iter().map(String::as_str).collect(),
                true,
                50.0,
            ));
        }
        baseline_entries.push((target.as_str(), borrowed.clone(), true, 100.0));
        candidate_entries.push((target.as_str(), borrowed.clone(), true, 500.0));
        candidate_entries.push((target.as_str(), borrowed.clone(), true, 50.0));
        candidate_entries.push((target.as_str(), borrowed, true, 50.0));
        let baseline = stub_archive("baseline-run", baseline_entries);
        let candidate = stub_archive("candidate-run", candidate_entries);
        let envelope = compare_parity_envelopes(&baseline, &candidate);
        let member = envelope
            .members
            .iter()
            .find(|member| member.input == target)
            .expect("target member present");
        assert_eq!(member.candidate_total_ms, 50.0);
        assert_eq!(member.status, "ok");
        assert_eq!(envelope.status, "ok");
    }

    #[test]
    fn median_of_matches_partial_cmp_order_on_finite_inputs() {
        assert_eq!(median_of(Vec::new()), None);
        assert_eq!(median_of(vec![2.0]), Some(2.0));
        assert_eq!(median_of(vec![3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(median_of(vec![1.0, 2.0, 3.0, 4.0]), Some(3.0));
        let boundary = vec![
            f64::MIN,
            -1e300,
            -1.0,
            -f64::MIN_POSITIVE,
            -0.0,
            1.0,
            f64::MIN_POSITIVE,
            1e300,
            f64::MAX,
        ];
        let mut by_partial = boundary.clone();
        by_partial.sort_by(|first, second| {
            first
                .partial_cmp(second)
                .expect("boundary total_ms is finite")
        });
        let mut by_total = boundary.clone();
        by_total.sort_by(|first, second| first.total_cmp(second));
        assert_eq!(by_total, by_partial);
        assert_eq!(
            median_of(boundary.clone()),
            Some(by_partial[by_partial.len() / 2])
        );
        assert_eq!(median_of(vec![50.0, 500.0, 50.0]), Some(50.0));
    }

    #[test]
    fn parity_envelope_rejects_missing_member() {
        let corpus = parity_corpus();
        let baseline_entries = corpus
            .iter()
            .map(|(input, expected)| {
                (
                    input.as_str(),
                    expected.iter().map(String::as_str).collect(),
                    true,
                    100.0,
                )
            })
            .collect();
        let mut candidate_entries = Vec::new();
        for (index, (input, expected)) in corpus.iter().enumerate() {
            if index == 0 {
                continue;
            }
            candidate_entries.push((
                input.as_str(),
                expected.iter().map(String::as_str).collect(),
                true,
                50.0,
            ));
        }
        let baseline = stub_archive("baseline-run", baseline_entries);
        let candidate = stub_archive("candidate-run", candidate_entries);
        let envelope = compare_parity_envelopes(&baseline, &candidate);
        assert_eq!(envelope.status, "verification_failed");
    }

    #[test]
    fn primality_corpus_has_required_members() {
        let corpus = primality_corpus();
        assert_eq!(corpus.len(), 42);
        for digits in [
            "561", "1105", "1729", "2465", "2821", "6601", "8911", "10585", "15841", "29341",
            "41041", "46657", "52633", "62745", "63973", "75361",
        ] {
            assert!(
                corpus
                    .iter()
                    .any(|(input, expected)| input == digits && !expected)
            );
        }
        for digits in [
            "2047",
            "3215031751",
            "2152302898747",
            "3474749660383",
            "341550071728321",
            "3825123056546413051",
        ] {
            assert!(
                corpus
                    .iter()
                    .any(|(input, expected)| input == digits && !expected)
            );
        }
        for digits in [
            "2",
            "3",
            "5",
            "7",
            "11",
            "13",
            "17",
            "19",
            "23",
            "29",
            "31",
            "37",
            "41",
            "43",
            "47",
            "101",
            "997",
            "1000003",
            "2147483647",
            "170141183460469231731687303715884105727",
        ] {
            assert!(
                corpus
                    .iter()
                    .any(|(input, expected)| input == digits && *expected)
            );
        }
    }

    #[test]
    fn primality_member_status_rejects_mismatch() {
        assert_eq!(
            primality_member_status(false, false, false, false, Some(false)),
            "ok"
        );
        assert_eq!(
            primality_member_status(true, true, true, true, Some(true)),
            "ok"
        );
        assert_eq!(
            primality_member_status(false, false, false, false, None),
            "ok"
        );
        assert_eq!(
            primality_member_status(false, true, false, false, Some(false)),
            "verification_failed"
        );
        assert_eq!(
            primality_member_status(false, false, false, false, Some(true)),
            "verification_failed"
        );
        assert_eq!(
            primality_member_status(true, false, true, true, None),
            "verification_failed"
        );
    }

    #[test]
    fn primality_envelope_records_real_gates() {
        let members = build_primality_members("run-x", "v0.1.5.2.7", None);
        assert_eq!(members.len(), primality_corpus().len());
        assert!(members.iter().all(|member| member.status == "ok"));
        for member in &members {
            assert_eq!(member.run_id, "run-x");
            assert_eq!(member.commit_or_id, "v0.1.5.2.7");
            assert_eq!(member.expected_is_prime, member.probable);
            assert_eq!(member.expected_is_prime, member.bpsw);
            assert_eq!(member.expected_is_prime, member.frobenius_gated);
            if let Some(proof) = member.aks {
                assert_eq!(proof, member.expected_is_prime);
            }
        }
        let by_input = |digits: &str| {
            members
                .iter()
                .find(|member| member.input == digits)
                .expect("member present")
        };
        assert!(!by_input("561").expected_is_prime);
        assert!(!by_input("3215031751").expected_is_prime);
        assert!(by_input("1000003").expected_is_prime);
        assert!(by_input("170141183460469231731687303715884105727").expected_is_prime);
        assert_eq!(
            by_input("170141183460469231731687303715884105727").aks,
            None
        );
    }

    #[test]
    fn primality_compare_folds_into_envelope_status() {
        let corpus = parity_corpus();
        let baseline_entries = corpus
            .iter()
            .map(|(input, expected)| {
                (
                    input.as_str(),
                    expected.iter().map(String::as_str).collect(),
                    true,
                    100.0,
                )
            })
            .collect();
        let candidate_entries = corpus
            .iter()
            .map(|(input, expected)| {
                (
                    input.as_str(),
                    expected.iter().map(String::as_str).collect(),
                    true,
                    50.0,
                )
            })
            .collect();
        let baseline = stub_archive("baseline-run", baseline_entries);
        let candidate = stub_archive("candidate-run", candidate_entries);
        let envelope = compare_parity_envelopes(&baseline, &candidate);
        assert_eq!(envelope.primality_status, "ok");
        assert!(
            envelope
                .primality_members
                .iter()
                .all(|member| member.status == "ok")
        );
        assert_eq!(envelope.status, "ok");
    }
}
