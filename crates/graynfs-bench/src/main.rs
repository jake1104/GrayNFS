use std::fs;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use clap::Parser;
use graynfs_core::{FactorOptions, Integer, factor_with_progress, probable_prime};
use num_traits::Num;
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(
    name = "graynfs-bench",
    version,
    about = "Benchmark GrayNFS across factorization difficulty tiers"
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

#[derive(Clone, Debug, Serialize)]
struct StageSample {
    stage: String,
    current: u64,
    total: u64,
    message: String,
    start_ms: f64,
    end_ms: f64,
    duration_ms: f64,
}

#[derive(Clone, Debug, Serialize)]
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

#[derive(Clone, Debug, Serialize)]
struct BenchArchive {
    task_id: String,
    tier: String,
    timeout_seconds: u64,
    threads_requested: Option<usize>,
    threads_effective: usize,
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
        tier: cli.tier.clone(),
        timeout_seconds: cli.timeout_seconds,
        threads_requested: cli.threads,
        threads_effective,
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
