use clap::{CommandFactory, Parser, Subcommand};
use graynfs_core::batch::{batch_job_options, batch_slots};
use graynfs_core::{
    CancellationToken, FactorError, FactorOptions, FactorReport, ProgressEvent, StageTiming,
    factor_with_timings_cancelled, parse_decimal,
};
use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, mpsc};

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser, Debug)]
#[command(name = "graynfs", version, about = "Factor integers with GNFS")]
struct Cli {
    number: Option<String>,
    #[arg(
        long,
        help = "Force a fixed polynomial degree instead of the automatic race"
    )]
    degree: Option<usize>,
    #[arg(
        long,
        default_value_t = 2_000,
        help = "Bound of the rational factor base"
    )]
    rational_base_bound: usize,
    #[arg(
        long,
        default_value_t = 2_000,
        help = "Bound of the algebraic factor base"
    )]
    algebraic_base_bound: usize,
    #[arg(long, default_value_t = 700, help = "Sieve bound on a")]
    sieve_a_bound: usize,
    #[arg(long, default_value_t = 120, help = "Sieve bound on b")]
    sieve_b_bound: usize,
    #[arg(
        long,
        default_value_t = 50_000,
        help = "Upper bound of a residual large prime"
    )]
    large_prime_bound: usize,
    #[arg(
        long,
        help = "Seed for the deterministic RNG (absent = 0x9e3779b97f4a7c15)"
    )]
    seed: Option<u64>,
    #[arg(
        long,
        default_value_t = 60,
        help = "Overall budget in seconds; on expiry the run ends with a clean ResourceLimit"
    )]
    timeout_seconds: u64,
    #[arg(
        long,
        help = "Process-wide CPU-time budget in seconds (allocated wall-clock multiplied by threads); on expiry the run ends with a clean ResourceLimit"
    )]
    max_cpu_seconds: Option<u64>,
    #[arg(
        long,
        help = "Cap estimated sieve and matrix memory in MiB before allocating (unset; 1 GiB budget)"
    )]
    max_memory_mb: Option<usize>,
    #[arg(long, help = "Disable the Pollard-Rho fallback")]
    no_fallback: bool,
    #[arg(
        long,
        help = "Skip trial division and the Fermat/p-1/p+1/ECM probes (square/power checks still run)"
    )]
    skip_preprocessing: bool,
    #[arg(long, help = "Print a per-stage timing breakdown")]
    timings: bool,
    #[arg(
        long,
        help = "Allow the general siever to keep relations with two or three large primes and pair them cyclically"
    )]
    double_large_primes: bool,
    #[arg(long, help = "Skip the Fermat close-factor probe")]
    no_fermat: bool,
    #[arg(long, help = "Skip Pollard p-1")]
    no_pm1: bool,
    #[arg(long, help = "Skip Williams p+1")]
    no_pp1: bool,
    #[arg(long, help = "Skip Lenstra ECM")]
    no_ecm: bool,
    #[arg(
        long,
        help = "Disable the SNFS special-form path inside GNFS (generic polynomials only)"
    )]
    no_snfs: bool,
    #[arg(
        long,
        help = "Route degree 2 to the classic single-polynomial QS path instead of SIQS"
    )]
    no_mpqs: bool,
    #[arg(
        long,
        help = "Worker thread count for the parallel stages (defaults to all logical cores)"
    )]
    threads: Option<usize>,
    #[arg(
        long,
        help = "Mark timing records as performance-mode runs without changing factorization behavior"
    )]
    perf_mode: bool,
    #[arg(
        long,
        help = "Write GNFS collector checkpoints to PATH during long runs"
    )]
    checkpoint: Option<PathBuf>,
    #[arg(long, help = "Resume a GNFS run from the checkpoint file at PATH")]
    resume: Option<PathBuf>,
    #[arg(long, help = "Emit the result as a single JSON object on stdout")]
    json: bool,
    #[arg(long, help = "Stream per-event factorization progress to stderr")]
    progress: bool,
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Batch {
        #[arg(long, help = "Path to a text file holding one integer per line")]
        input: PathBuf,
        #[arg(long, help = "Emit JSON Lines records (the default batch encoding)")]
        jsonl: bool,
        #[arg(
            long,
            help = "Record per-line failures and keep going instead of aborting"
        )]
        continue_on_error: bool,
        #[arg(
            long,
            default_value_t = 1,
            help = "Batch worker count; records still emit in input order"
        )]
        jobs: usize,
    },
}

fn options_from_cli(cli: &Cli) -> FactorOptions {
    FactorOptions {
        degree: cli.degree,
        rational_base_bound: cli.rational_base_bound,
        algebraic_base_bound: cli.algebraic_base_bound,
        sieve_a_bound: cli.sieve_a_bound,
        sieve_b_bound: cli.sieve_b_bound,
        large_prime_bound: cli.large_prime_bound,
        seed: cli.seed.unwrap_or(FactorOptions::default().seed),
        max_duration: Some(std::time::Duration::from_secs(cli.timeout_seconds)),
        max_cpu_seconds: cli.max_cpu_seconds,
        max_memory_mb: cli.max_memory_mb,
        allow_fallback: !cli.no_fallback,
        skip_preprocessing: cli.skip_preprocessing,
        enable_double_large_primes: cli.double_large_primes,
        use_fermat: !cli.no_fermat,
        use_pm1: !cli.no_pm1,
        use_pp1: !cli.no_pp1,
        use_ecm: !cli.no_ecm,
        use_snfs: !cli.no_snfs,
        use_mpqs: !cli.no_mpqs,
        threads: cli.threads,
        perf_mode: cli.perf_mode,
        checkpoint_path: cli.checkpoint.clone(),
        resume_path: cli.resume.clone(),
        ..FactorOptions::default()
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match &cli.command {
        Some(Commands::Batch {
            input,
            continue_on_error,
            jobs,
            ..
        }) => run_batch(&cli, input, *continue_on_error, *jobs),
        None => run_single(&cli),
    }
}

struct BatchRecord {
    input: String,
    factors: Vec<String>,
    status: String,
    method: Option<String>,
    elapsed_ms: u128,
    error: Option<String>,
}

fn batch_record_value(record: &BatchRecord) -> serde_json::Value {
    serde_json::json!({
        "input": record.input,
        "factors": record.factors,
        "status": record.status,
        "method": record.method,
        "elapsed_ms": record.elapsed_ms,
        "error": record.error,
    })
}

fn factor_batch_line(
    raw: &str,
    options: &FactorOptions,
    token: &CancellationToken,
    show_progress: bool,
) -> BatchRecord {
    factor_batch_line_with_sink(raw, options, token, &|event| {
        if show_progress {
            eprintln!("{}", format_progress(event));
        }
    })
}

fn factor_batch_line_with_sink(
    raw: &str,
    options: &FactorOptions,
    token: &CancellationToken,
    on_progress: &(dyn Fn(&ProgressEvent) + Sync),
) -> BatchRecord {
    let started = std::time::Instant::now();
    let number = match parse_decimal(raw) {
        Ok(value) => value,
        Err(error) => {
            let (_, detail) = error_parts(&error);
            return BatchRecord {
                input: raw.to_string(),
                factors: Vec::new(),
                status: "error".to_string(),
                method: None,
                elapsed_ms: started.elapsed().as_millis(),
                error: Some(detail),
            };
        }
    };
    let (result, _) = factor_with_timings_cancelled(&number, options, token, |event| {
        if token.is_cancelled() {
            return Err(FactorError::Cancelled);
        }
        on_progress(event);
        Ok(())
    });
    match result {
        Ok(report) => BatchRecord {
            input: raw.to_string(),
            factors: report.factors,
            status: "ok".to_string(),
            method: Some(report.method),
            elapsed_ms: started.elapsed().as_millis(),
            error: None,
        },
        Err(error) => {
            let (category, detail) = error_parts(&error);
            BatchRecord {
                input: raw.to_string(),
                factors: Vec::new(),
                status: "error".to_string(),
                method: None,
                elapsed_ms: started.elapsed().as_millis(),
                error: Some(format!("{category}: {detail}")),
            }
        }
    }
}

fn run_batch_lines(
    reader: impl BufRead,
    options: &FactorOptions,
    token: &CancellationToken,
    show_progress: bool,
    continue_on_error: bool,
    emit: &mut dyn FnMut(&BatchRecord),
) -> std::io::Result<bool> {
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let record = factor_batch_line(trimmed, options, token, show_progress);
        let failed = record.status == "error";
        let cancelled = batch_record_is_cancelled(&record);
        emit(&record);
        if cancelled || (failed && !continue_on_error) {
            return Ok(false);
        }
        if token.is_cancelled() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn batch_record_is_cancelled(record: &BatchRecord) -> bool {
    record
        .error
        .as_deref()
        .is_some_and(|text| text.starts_with("cancelled"))
}

fn run_batch_parallel(
    reader: impl BufRead + Send,
    options: &FactorOptions,
    token: &CancellationToken,
    show_progress: bool,
    continue_on_error: bool,
    jobs: usize,
    emit: &mut dyn FnMut(&BatchRecord),
) -> std::io::Result<bool> {
    let workers = batch_slots(options, jobs).max(1);
    let scoped = batch_job_options(options, jobs);
    let (task_tx, task_rx) =
        mpsc::sync_channel::<(usize, String)>(workers.saturating_mul(2).max(1));
    let (result_tx, result_rx) = mpsc::channel::<(usize, BatchRecord)>();
    let task_rx = Arc::new(Mutex::new(task_rx));
    let stop = Arc::new(AtomicBool::new(false));
    let dispatched = Arc::new(AtomicUsize::new(0));
    let dispatcher_done = Arc::new(AtomicBool::new(false));
    let io_error = Arc::new(Mutex::new(None::<String>));
    let progress_lock = Arc::new(Mutex::new(()));
    let mut completed = true;
    std::thread::scope(|scope| {
        {
            let task_tx = task_tx;
            let stop = Arc::clone(&stop);
            let dispatched = Arc::clone(&dispatched);
            let dispatcher_done = Arc::clone(&dispatcher_done);
            let io_error = Arc::clone(&io_error);
            scope.spawn(move || {
                let mut lines = reader.lines();
                let mut seq = 0usize;
                'dispatch: loop {
                    if stop.load(Ordering::SeqCst) {
                        break;
                    }
                    let line = match lines.next() {
                        None => break,
                        Some(Ok(text)) => text,
                        Some(Err(error)) => {
                            if let Ok(mut slot) = io_error.lock() {
                                *slot = Some(error.to_string());
                            }
                            break;
                        }
                    };
                    let trimmed = line.trim().to_string();
                    if trimmed.is_empty() {
                        continue;
                    }
                    let mut pending = Some((seq, trimmed));
                    seq = seq.saturating_add(1);
                    while let Some(job) = pending.take() {
                        if stop.load(Ordering::SeqCst) {
                            break 'dispatch;
                        }
                        match task_tx.try_send(job) {
                            Ok(()) => {
                                dispatched.fetch_add(1, Ordering::SeqCst);
                            }
                            Err(mpsc::TrySendError::Full(job)) => {
                                pending = Some(job);
                                std::thread::sleep(std::time::Duration::from_millis(1));
                            }
                            Err(mpsc::TrySendError::Disconnected(_)) => {
                                break 'dispatch;
                            }
                        }
                    }
                }
                dispatcher_done.store(true, Ordering::SeqCst);
            });
        }
        for _ in 0..workers {
            let task_rx = Arc::clone(&task_rx);
            let result_tx = result_tx.clone();
            let scoped = scoped.clone();
            let token = token.clone();
            let progress_lock = Arc::clone(&progress_lock);
            let stop = Arc::clone(&stop);
            scope.spawn(move || {
                while let Ok((seq, line)) = task_rx.lock().expect("batch task lock").recv() {
                    if stop.load(Ordering::SeqCst) {
                        continue;
                    }
                    let record = if show_progress {
                        factor_batch_line_with_sink(&line, &scoped, &token, &|event| {
                            if let Ok(_guard) = progress_lock.lock() {
                                eprintln!("{}", format_progress(event));
                            }
                        })
                    } else {
                        factor_batch_line_with_sink(&line, &scoped, &token, &|_| {})
                    };
                    if result_tx.send((seq, record)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(result_tx);
        let mut pending: BTreeMap<usize, BatchRecord> = BTreeMap::new();
        let mut next_seq = 0usize;
        loop {
            while let Some(record) = pending.remove(&next_seq) {
                let failed = record.status == "error";
                let cancelled = batch_record_is_cancelled(&record);
                emit(&record);
                next_seq = next_seq.saturating_add(1);
                if cancelled || (failed && !continue_on_error) || token.is_cancelled() {
                    completed = false;
                    stop.store(true, Ordering::SeqCst);
                    for _ in result_rx.iter() {}
                    return;
                }
            }
            if dispatcher_done.load(Ordering::SeqCst)
                && next_seq >= dispatched.load(Ordering::SeqCst)
                && pending.is_empty()
            {
                break;
            }
            match result_rx.recv() {
                Ok((seq, record)) => {
                    pending.insert(seq, record);
                }
                Err(_) => break,
            }
        }
        while let Some(record) = pending.remove(&next_seq) {
            let failed = record.status == "error";
            let cancelled = batch_record_is_cancelled(&record);
            emit(&record);
            next_seq = next_seq.saturating_add(1);
            if cancelled || (failed && !continue_on_error) || token.is_cancelled() {
                completed = false;
                stop.store(true, Ordering::SeqCst);
                break;
            }
        }
        if token.is_cancelled() {
            completed = false;
        }
    });
    if let Ok(slot) = io_error.lock()
        && let Some(detail) = slot.as_ref()
    {
        return Err(std::io::Error::other(detail.clone()));
    }
    Ok(completed)
}

fn batch_exit_code(failed: bool, cancelled: bool) -> ExitCode {
    if failed || cancelled {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}

fn run_batch(cli: &Cli, input: &PathBuf, continue_on_error: bool, jobs: usize) -> ExitCode {
    if jobs == 0 {
        eprintln!("error: input: --jobs must be at least one");
        return ExitCode::FAILURE;
    }
    if cli.checkpoint.is_some() || cli.resume.is_some() {
        eprintln!("error: input: --checkpoint and --resume are single-number options");
        return ExitCode::FAILURE;
    }
    let file = match std::fs::File::open(input) {
        Ok(handle) => handle,
        Err(error) => {
            eprintln!(
                "error: input: cannot open batch file {}: {error}",
                input.display()
            );
            return ExitCode::FAILURE;
        }
    };
    let options = options_from_cli(cli);
    let show_progress = cli.progress;
    let token = CancellationToken::new();
    let interrupt = token.clone();
    let _ = ctrlc::set_handler(move || {
        interrupt.cancel();
    });
    let reader = std::io::BufReader::new(file);
    let mut failed = false;
    let mut emit = |record: &BatchRecord| {
        if record.status == "error" {
            failed = true;
        }
        print_json(&batch_record_value(record));
    };
    let completed_ok = if jobs <= 1 {
        run_batch_lines(
            reader,
            &options,
            &token,
            show_progress,
            continue_on_error,
            &mut emit,
        )
    } else {
        run_batch_parallel(
            reader,
            &options,
            &token,
            show_progress,
            continue_on_error,
            jobs,
            &mut emit,
        )
    };
    if let Err(error) = completed_ok {
        eprintln!("error: io: cannot read batch line: {error}");
        return ExitCode::FAILURE;
    }
    batch_exit_code(failed, token.is_cancelled())
}

fn run_single(cli: &Cli) -> ExitCode {
    let number_text = match &cli.number {
        Some(value) => value.clone(),
        None => Cli::command()
            .error(
                clap::error::ErrorKind::MissingRequiredArgument,
                "the following required arguments were not provided:\n  <NUMBER>",
            )
            .exit(),
    };
    let number = match parse_decimal(&number_text) {
        Ok(value) => value,
        Err(error) => {
            if cli.json {
                print_json(&json_error_value(&error, None));
            } else {
                report_error(&error);
            }
            return ExitCode::FAILURE;
        }
    };
    let options = options_from_cli(cli);
    let show_progress = cli.progress;
    let token = CancellationToken::new();
    let interrupt = token.clone();
    let _ = ctrlc::set_handler(move || {
        interrupt.cancel();
    });
    let (result, timings) = factor_with_timings_cancelled(&number, &options, &token, |event| {
        if token.is_cancelled() {
            return Err(FactorError::Cancelled);
        }
        if show_progress {
            eprintln!("{}", format_progress(event));
        }
        Ok(())
    });
    match result {
        Ok(report) => {
            if cli.json {
                let selected = cli.timings.then_some(report.timings.as_slice());
                print_json(&json_report_value(&report, selected));
            } else {
                println!("N: {}", report.input);
                println!("method: {}", report.method);
                println!(
                    "degree: {}",
                    report
                        .degree
                        .map_or_else(|| "n/a".to_string(), |value| value.to_string())
                );
                println!("factors: {}", report.factors.join(" * "));
                println!("relations: {}", report.relations);
                println!("rank: {}", report.rank);
                println!("elapsed_ms: {}", report.elapsed_ms);
                if cli.timings {
                    print_timings(&report.timings);
                }
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            if cli.json {
                let selected = cli.timings.then_some(timings.as_slice());
                print_json(&json_error_value(&error, selected));
            } else {
                if cli.timings {
                    print_timings(&timings);
                }
                report_error(&error);
            }
            ExitCode::FAILURE
        }
    }
}

fn print_json(value: &serde_json::Value) {
    println!("{value}");
}

fn format_progress(event: &ProgressEvent) -> String {
    format!(
        "progress: stage={} current={} total={} message={}",
        event.stage, event.current, event.total, event.message
    )
}

fn timings_json(entries: &[StageTiming]) -> serde_json::Value {
    serde_json::Value::Array(
        entries
            .iter()
            .map(|entry| {
                serde_json::json!({
                    "stage": entry.stage,
                    "degree": entry.degree,
                    "elapsed_ms": entry.elapsed_ms,
                    "detail": entry.detail,
                })
            })
            .collect::<Vec<_>>(),
    )
}

fn json_report_value(report: &FactorReport, timings: Option<&[StageTiming]>) -> serde_json::Value {
    let mut value = serde_json::json!({
        "input": report.input,
        "method": report.method,
        "degree": report.degree,
        "factors": report.factors,
        "factor_confidence": report.factor_confidence,
        "polynomial": report.polynomial,
        "root": report.root,
        "relations": report.relations,
        "rank": report.rank,
        "elapsed_ms": report.elapsed_ms,
    });
    if let Some(entries) = timings {
        value["timings"] = timings_json(entries);
    }
    value
}

fn json_error_value(error: &FactorError, timings: Option<&[StageTiming]>) -> serde_json::Value {
    let (category, detail) = error_parts(error);
    let mut value = serde_json::json!({
        "error": {
            "category": category,
            "detail": detail,
        }
    });
    if let Some(entries) = timings {
        value["timings"] = timings_json(entries);
    }
    value
}

fn print_timings(entries: &[StageTiming]) {
    for entry in entries {
        let degree = entry
            .degree
            .map_or_else(|| "n/a".to_string(), |value| value.to_string());
        println!(
            "stage: {:<20} degree: {:<4} elapsed_ms: {:<8} {}",
            entry.stage, degree, entry.elapsed_ms, entry.detail
        );
    }
}

fn error_parts(error: &FactorError) -> (&'static str, String) {
    match error {
        FactorError::InvalidInput(detail) => ("input", detail.clone()),
        FactorError::ResourceLimit(detail) => ("resource-limit", detail.clone()),
        FactorError::Cancelled => ("cancelled", "the operation was cancelled".to_string()),
        FactorError::Checkpoint(detail) => ("checkpoint", detail.clone()),
        other => ("internal", other.to_string()),
    }
}

fn report_error(error: &FactorError) {
    let (category, detail) = error_parts(error);
    eprintln!("error: {category}: {detail}");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_report() -> FactorReport {
        FactorReport {
            input: "100160063".to_string(),
            factors: vec!["10007".to_string(), "10009".to_string()],
            factor_confidence: vec!["proven_prime".to_string(), "proven_prime".to_string()],
            method: "gnfs".to_string(),
            degree: Some(3),
            polynomial: None,
            root: None,
            relations: 604,
            rank: 293,
            elapsed_ms: 120,
            timings: vec![StageTiming {
                stage: "preprocessing".to_string(),
                degree: Some(2),
                elapsed_ms: 1,
                detail: "checked".to_string(),
            }],
        }
    }

    #[test]
    fn json_report_has_identity_fields_and_round_trips() {
        let report = sample_report();
        let value = json_report_value(&report, None);
        assert_eq!(value["input"], "100160063");
        assert_eq!(value["method"], "gnfs");
        assert_eq!(value["degree"], 3);
        assert_eq!(value["factors"], serde_json::json!(["10007", "10009"]));
        assert_eq!(
            value["factor_confidence"],
            serde_json::json!(["proven_prime", "proven_prime"])
        );
        assert!(value["polynomial"].is_null());
        assert!(value["root"].is_null());
        assert_eq!(value["relations"], 604);
        assert_eq!(value["rank"], 293);
        assert_eq!(value["elapsed_ms"], 120);
        assert!(value.get("timings").is_none());
        let text = value.to_string();
        let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed, value);
    }

    #[test]
    fn json_report_emits_confidence_polynomial_root_additive_only() {
        let report = sample_report();
        let value = json_report_value(&report, None);
        let mut keys: Vec<&str> = value
            .as_object()
            .expect("report object")
            .keys()
            .map(String::as_str)
            .collect();
        keys.sort_unstable();
        assert_eq!(
            keys,
            vec![
                "degree",
                "elapsed_ms",
                "factor_confidence",
                "factors",
                "input",
                "method",
                "polynomial",
                "rank",
                "relations",
                "root"
            ]
        );
        let mut with_poly = sample_report();
        with_poly.polynomial = Some("x^2 + 1".to_string());
        with_poly.root = Some("12345".to_string());
        let poly_value = json_report_value(&with_poly, None);
        assert_eq!(poly_value["polynomial"], "x^2 + 1");
        assert_eq!(poly_value["root"], "12345");
        assert_eq!(
            poly_value["factor_confidence"],
            serde_json::json!(["proven_prime", "proven_prime"])
        );
    }

    #[test]
    fn json_report_includes_timings_only_when_requested() {
        let report = sample_report();
        let value = json_report_value(&report, Some(report.timings.as_slice()));
        let timings = value["timings"].as_array().expect("timings array");
        assert_eq!(timings.len(), 1);
        assert_eq!(timings[0]["stage"], "preprocessing");
        assert_eq!(timings[0]["degree"], 2);
        assert_eq!(timings[0]["elapsed_ms"], 1);
        assert_eq!(timings[0]["detail"], "checked");
    }

    #[test]
    fn json_report_passes_scheduler_stages_with_stable_keys() {
        let mut report = sample_report();
        report.timings = vec![
            StageTiming {
                stage: "reclaim".to_string(),
                degree: None,
                elapsed_ms: 0,
                detail: "reclaim=none siqs_width=700 ecm_digits=15 workers=16".to_string(),
            },
            StageTiming {
                stage: "joint".to_string(),
                degree: None,
                elapsed_ms: 0,
                detail: "shares=fallback winner=ecm banked=0".to_string(),
            },
            StageTiming {
                stage: "cascade_projection".to_string(),
                degree: Some(2),
                elapsed_ms: 0,
                detail: "cost_ecm_ms=10 cost_siqs_ms=20 lower_cost=siqs".to_string(),
            },
            StageTiming {
                stage: "ecm".to_string(),
                degree: None,
                elapsed_ms: 5,
                detail: "divisor=7 curves=2 ecm_ms=5".to_string(),
            },
        ];
        let value = json_report_value(&report, Some(report.timings.as_slice()));
        let timings = value["timings"].as_array().expect("timings array");
        assert_eq!(timings.len(), 4);
        let stages: Vec<&str> = timings
            .iter()
            .map(|entry| entry["stage"].as_str().expect("stage string"))
            .collect();
        assert_eq!(
            stages,
            vec!["reclaim", "joint", "cascade_projection", "ecm"]
        );
        for entry in timings {
            let object = entry.as_object().expect("timing object");
            let mut keys: Vec<&str> = object.keys().map(String::as_str).collect();
            keys.sort_unstable();
            assert_eq!(keys, vec!["degree", "detail", "elapsed_ms", "stage"]);
        }
        assert_eq!(timings[1]["detail"], "shares=fallback winner=ecm banked=0");
        assert_eq!(value["method"], "gnfs");
    }

    #[test]
    fn json_report_serializes_nullable_degree() {
        let mut report = sample_report();
        report.degree = None;
        let value = json_report_value(&report, None);
        assert!(value["degree"].is_null());
    }

    #[test]
    fn json_error_carries_category_and_detail() {
        let value = json_error_value(&FactorError::ResourceLimit("deadline".to_string()), None);
        assert_eq!(value["error"]["category"], "resource-limit");
        assert_eq!(value["error"]["detail"], "deadline");
        assert!(value.get("timings").is_none());
    }

    #[test]
    fn json_error_includes_timings_only_when_requested() {
        let error = FactorError::ResourceLimit("deadline".to_string());
        let entries = vec![StageTiming {
            stage: "memory_shadow".to_string(),
            degree: Some(3),
            elapsed_ms: 0,
            detail: "site=general-collect-precheck verdict=reject".to_string(),
        }];
        let with = json_error_value(&error, Some(entries.as_slice()));
        let timings = with["timings"].as_array().expect("timings array");
        assert_eq!(timings.len(), 1);
        assert_eq!(timings[0]["stage"], "memory_shadow");
        let without = json_error_value(&error, None);
        assert!(without.get("timings").is_none());
    }

    #[test]
    fn error_parts_maps_input() {
        let (category, detail) = error_parts(&FactorError::InvalidInput("bad digit".to_string()));
        assert_eq!(category, "input");
        assert_eq!(detail, "bad digit");
    }

    #[test]
    fn error_parts_maps_resource_limit() {
        let (category, detail) = error_parts(&FactorError::ResourceLimit("deadline".to_string()));
        assert_eq!(category, "resource-limit");
        assert_eq!(detail, "deadline");
    }

    #[test]
    fn error_parts_maps_cancelled() {
        let (category, detail) = error_parts(&FactorError::Cancelled);
        assert_eq!(category, "cancelled");
        assert!(!detail.is_empty());
    }

    #[test]
    fn error_parts_maps_internal_variants() {
        let variants = [
            FactorError::NoPolynomial,
            FactorError::SearchExhausted,
            FactorError::InsufficientRelations {
                found: 1,
                required: 2,
            },
            FactorError::AlgebraicSquareRootFailed,
            FactorError::LinearAlgebra("rank".to_string()),
            FactorError::Factorization("boom".to_string()),
            FactorError::UnsupportedDegree(3),
        ];
        for error in &variants {
            let (category, detail) = error_parts(error);
            assert_eq!(category, "internal");
            assert!(!detail.is_empty());
        }
    }

    #[test]
    fn error_parts_maps_checkpoint() {
        let (category, detail) = error_parts(&FactorError::Checkpoint(
            "checkpoint checksum mismatch".to_string(),
        ));
        assert_eq!(category, "checkpoint");
        assert!(detail.contains("checksum"));
    }

    #[test]
    fn json_error_shape_has_category_and_detail_for_each_category() {
        let cases = [
            (FactorError::InvalidInput("bad digit".to_string()), "input"),
            (
                FactorError::ResourceLimit("deadline".to_string()),
                "resource-limit",
            ),
            (FactorError::Cancelled, "cancelled"),
            (FactorError::NoPolynomial, "internal"),
        ];
        for (error, category) in &cases {
            let value = json_error_value(error, None);
            assert_eq!(value["error"]["category"], *category);
            assert!(value["error"]["detail"].is_string());
            assert!(value.get("input").is_none());
        }
    }

    #[test]
    fn exit_codes_pin_success_zero_and_failure_nonzero() {
        assert_eq!(ExitCode::SUCCESS, ExitCode::from(0));
        assert_ne!(ExitCode::FAILURE, ExitCode::from(0));
    }

    #[test]
    fn batch_continue_golden_three_lines_one_failure() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        let reader = std::io::BufReader::new(std::io::Cursor::new("15\nabc\n21\n"));
        let mut records = Vec::new();
        let completed = run_batch_lines(reader, &options, &token, false, true, &mut |record| {
            records.push(batch_record_value(record));
        })
        .expect("streaming batch reads");
        assert!(completed);
        assert_eq!(records.len(), 3);
        assert_eq!(records[0]["input"], "15");
        assert_eq!(records[0]["status"], "ok");
        assert_eq!(records[0]["factors"], serde_json::json!(["3", "5"]));
        assert!(records[0]["method"].is_string());
        assert!(records[0]["error"].is_null());
        assert_eq!(records[1]["input"], "abc");
        assert_eq!(records[1]["status"], "error");
        assert_eq!(records[1]["factors"], serde_json::json!([]));
        assert!(records[1]["error"].is_string());
        assert_eq!(records[2]["input"], "21");
        assert_eq!(records[2]["status"], "ok");
        assert_eq!(records[2]["factors"], serde_json::json!(["3", "7"]));
        assert!(records[2]["error"].is_null());
        for value in &records {
            let mut keys: Vec<&str> = value
                .as_object()
                .expect("batch record object")
                .keys()
                .map(String::as_str)
                .collect();
            keys.sort_unstable();
            assert_eq!(
                keys,
                vec![
                    "elapsed_ms",
                    "error",
                    "factors",
                    "input",
                    "method",
                    "status"
                ]
            );
            let text = value.to_string();
            let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
            assert_eq!(&parsed, value);
        }
    }

    #[test]
    fn batch_aborts_on_first_failure_without_continue_flag() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        let reader = std::io::BufReader::new(std::io::Cursor::new("15\nabc\n21\n"));
        let mut records = Vec::new();
        let completed = run_batch_lines(reader, &options, &token, false, false, &mut |record| {
            records.push(batch_record_value(record));
        })
        .expect("streaming batch reads");
        assert!(!completed);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0]["status"], "ok");
        assert_eq!(records[0]["factors"], serde_json::json!(["3", "5"]));
        assert_eq!(records[1]["status"], "error");
        assert!(records[1]["error"].is_string());
    }

    #[test]
    fn batch_streaming_skips_blank_lines() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        let reader = std::io::BufReader::new(std::io::Cursor::new("15\n\n   \n21\n"));
        let mut records = Vec::new();
        let completed = run_batch_lines(reader, &options, &token, false, true, &mut |record| {
            records.push(batch_record_value(record));
        })
        .expect("streaming batch reads");
        assert!(completed);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0]["input"], "15");
        assert_eq!(records[0]["status"], "ok");
        assert_eq!(records[1]["input"], "21");
        assert_eq!(records[1]["status"], "ok");
    }

    #[test]
    fn cli_batch_subcommand_parses_input_jsonl_continue_flags() {
        let cli = Cli::try_parse_from([
            "graynfs",
            "batch",
            "--input",
            "numbers.txt",
            "--jsonl",
            "--continue-on-error",
        ])
        .expect("parse batch subcommand");
        match cli.command.expect("batch command present") {
            Commands::Batch {
                input,
                jsonl,
                continue_on_error,
                jobs,
            } => {
                assert_eq!(input, PathBuf::from("numbers.txt"));
                assert!(jsonl);
                assert!(continue_on_error);
                assert_eq!(jobs, 1);
            }
        }
        assert!(cli.number.is_none());
        let bare = Cli::try_parse_from(["graynfs", "batch", "--input", "numbers.txt"])
            .expect("parse batch without option flags");
        match bare.command.expect("batch command present") {
            Commands::Batch {
                jsonl,
                continue_on_error,
                ..
            } => {
                assert!(!jsonl);
                assert!(!continue_on_error);
            }
        }
    }

    #[test]
    fn cli_batch_jobs_flag_defaults_to_one_and_parses_explicit() {
        let flagged =
            Cli::try_parse_from(["graynfs", "batch", "--input", "numbers.txt", "--jobs", "4"])
                .expect("parse batch with jobs");
        match flagged.command.expect("batch command present") {
            Commands::Batch { jobs, .. } => {
                assert_eq!(jobs, 4);
            }
        }
    }

    #[test]
    fn batch_jobs_zero_is_rejected_before_file_open() {
        let cli = Cli::try_parse_from([
            "graynfs",
            "batch",
            "--input",
            "definitely-missing-numbers.txt",
            "--jobs",
            "0",
        ])
        .expect("parse batch with zero jobs");
        match &cli.command {
            Some(Commands::Batch {
                input,
                continue_on_error,
                jobs,
                ..
            }) => {
                assert_eq!(
                    run_batch(&cli, input, *continue_on_error, *jobs),
                    ExitCode::FAILURE
                );
            }
            None => panic!("batch command present"),
        }
    }

    #[test]
    fn batch_parallel_jobs_one_matches_sequential() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        let input = "15\nabc\n21\n";
        let mut sequential = Vec::new();
        let reader = std::io::BufReader::new(std::io::Cursor::new(input));
        let done_seq = run_batch_lines(reader, &options, &token, false, true, &mut |record| {
            sequential.push(batch_record_value(record));
        })
        .expect("sequential batch reads");
        let mut parallel = Vec::new();
        let reader = std::io::BufReader::new(std::io::Cursor::new(input));
        let done_par =
            run_batch_parallel(reader, &options, &token, false, true, 1, &mut |record| {
                parallel.push(batch_record_value(record));
            })
            .expect("parallel batch reads");
        assert_eq!(done_seq, done_par);
        assert_eq!(sequential.len(), parallel.len());
        assert_eq!(sequential.len(), 3);
        for (left, right) in sequential.iter().zip(parallel.iter()) {
            assert_eq!(left["input"], right["input"]);
            assert_eq!(left["status"], right["status"]);
            assert_eq!(left["factors"], right["factors"]);
            assert_eq!(left["method"], right["method"]);
            assert_eq!(left["error"], right["error"]);
        }
    }

    #[test]
    fn batch_parallel_preserves_input_order() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        let lines = ["91", "abc", "15", "21", "33", "35", "39", "49", "77"];
        let input = lines.join("\n") + "\n";
        let mut records = Vec::new();
        let done = run_batch_parallel(
            std::io::BufReader::new(std::io::Cursor::new(input)),
            &options,
            &token,
            false,
            true,
            4,
            &mut |record| {
                records.push(batch_record_value(record));
            },
        )
        .expect("parallel batch reads");
        assert!(done);
        assert_eq!(records.len(), lines.len());
        for (record, expected) in records.iter().zip(lines.iter()) {
            assert_eq!(record["input"], *expected);
        }
        assert_eq!(records[2]["factors"], serde_json::json!(["3", "5"]));
        assert_eq!(records[1]["status"], "error");
        assert!(records[1]["error"].is_string());
    }

    #[test]
    fn batch_parallel_aborts_on_first_failure_without_continue_flag() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        let mut records = Vec::new();
        let done = run_batch_parallel(
            std::io::BufReader::new(std::io::Cursor::new("15\nabc\n21\n")),
            &options,
            &token,
            false,
            false,
            4,
            &mut |record| {
                records.push(batch_record_value(record));
            },
        )
        .expect("parallel batch reads");
        assert!(!done);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0]["status"], "ok");
        assert_eq!(records[0]["factors"], serde_json::json!(["3", "5"]));
        assert_eq!(records[1]["status"], "error");
    }

    #[test]
    fn batch_parallel_over_budget_queues_instead_of_aborting() {
        let options = FactorOptions {
            threads: Some(2),
            max_memory_mb: Some(0),
            ..FactorOptions::default()
        };
        assert_eq!(batch_slots(&options, 4), 1);
        let token = CancellationToken::new();
        let mut records = Vec::new();
        let done = run_batch_parallel(
            std::io::BufReader::new(std::io::Cursor::new("15\n21\n33\n")),
            &options,
            &token,
            false,
            true,
            4,
            &mut |record| {
                records.push(batch_record_value(record));
            },
        )
        .expect("queued batch reads");
        assert!(done);
        assert_eq!(records.len(), 3);
        assert_eq!(records[0]["input"], "15");
        assert_eq!(records[0]["status"], "ok");
        assert_eq!(records[0]["factors"], serde_json::json!(["3", "5"]));
        assert_eq!(records[1]["factors"], serde_json::json!(["3", "7"]));
        assert_eq!(records[2]["factors"], serde_json::json!(["3", "11"]));
    }

    #[test]
    fn batch_parallel_precancelled_token_aborts_after_first_record() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        token.cancel();
        let mut records = Vec::new();
        let done = run_batch_parallel(
            std::io::BufReader::new(std::io::Cursor::new("15\n21\n33\n")),
            &options,
            &token,
            false,
            true,
            4,
            &mut |record| {
                records.push(batch_record_value(record));
            },
        )
        .expect("parallel batch reads");
        assert!(!done);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0]["status"], "error");
        assert!(
            records[0]["error"]
                .as_str()
                .is_some_and(|text| text.starts_with("cancelled"))
        );
    }

    #[test]
    fn cli_checkpoint_resume_flags_default_absent_and_parse_explicit() {
        let bare = Cli::try_parse_from(["graynfs", "100160063"]).expect("parse without flags");
        assert_eq!(bare.checkpoint, None);
        assert_eq!(bare.resume, None);
        let flagged = Cli::try_parse_from([
            "graynfs",
            "100160063",
            "--checkpoint",
            "collect.json",
            "--resume",
            "resume.json",
        ])
        .expect("parse with checkpoint flags");
        assert_eq!(flagged.checkpoint, Some(PathBuf::from("collect.json")));
        assert_eq!(flagged.resume, Some(PathBuf::from("resume.json")));
        let options = options_from_cli(&flagged);
        assert_eq!(options.checkpoint_path, Some(PathBuf::from("collect.json")));
        assert_eq!(options.resume_path, Some(PathBuf::from("resume.json")));
        let plain = options_from_cli(&bare);
        assert_eq!(plain.checkpoint_path, None);
        assert_eq!(plain.resume_path, None);
    }

    #[test]
    fn cli_batch_rejects_checkpoint_and_resume_flags() {
        let cli = Cli::try_parse_from([
            "graynfs",
            "--checkpoint",
            "collect.json",
            "batch",
            "--input",
            "numbers.txt",
        ])
        .expect("parse batch with checkpoint flag");
        match &cli.command {
            Some(Commands::Batch {
                input,
                continue_on_error,
                jobs,
                ..
            }) => {
                assert_eq!(
                    run_batch(&cli, input, *continue_on_error, *jobs),
                    ExitCode::FAILURE
                );
            }
            None => panic!("batch command present"),
        }
        let cli = Cli::try_parse_from([
            "graynfs",
            "--resume",
            "resume.json",
            "batch",
            "--input",
            "numbers.txt",
        ])
        .expect("parse batch with resume flag");
        match &cli.command {
            Some(Commands::Batch {
                input,
                continue_on_error,
                jobs,
                ..
            }) => {
                assert_eq!(
                    run_batch(&cli, input, *continue_on_error, *jobs),
                    ExitCode::FAILURE
                );
            }
            None => panic!("batch command present"),
        }
    }

    #[test]
    fn cli_single_number_still_parses_as_positional() {
        let cli = Cli::try_parse_from(["graynfs", "100160063"]).expect("parse single number");
        assert_eq!(cli.number.as_deref(), Some("100160063"));
        assert!(cli.command.is_none());
    }

    #[test]
    fn cli_no_squfof_flag_is_rejected_as_unknown() {
        let result = Cli::try_parse_from(["graynfs", "100160063", "--no-squfof"]);
        assert!(result.is_err());
    }

    #[test]
    fn cli_progress_defaults_off_and_parses_when_set() {
        let bare = Cli::try_parse_from(["graynfs", "100160063"]).expect("parse without progress");
        assert!(!bare.progress);
        let flagged = Cli::try_parse_from(["graynfs", "100160063", "--progress"])
            .expect("parse with progress");
        assert!(flagged.progress);
    }

    #[test]
    fn format_progress_keeps_stage_current_total_message_shape() {
        let event = ProgressEvent {
            stage: "collect".to_string(),
            current: 12,
            total: 200,
            message: "sieving the quadratic base".to_string(),
        };
        assert_eq!(
            format_progress(&event),
            "progress: stage=collect current=12 total=200 message=sieving the quadratic base"
        );
    }

    #[test]
    fn cli_parses_seed_default_absent_and_explicit_round_trip() {
        let bare = Cli::try_parse_from(["graynfs", "100160063"]).expect("parse without seed");
        assert_eq!(bare.seed, None);
        assert_eq!(options_from_cli(&bare).seed, FactorOptions::default().seed);
        let flagged = Cli::try_parse_from(["graynfs", "100160063", "--seed", "111"])
            .expect("parse with seed");
        assert_eq!(flagged.seed, Some(111));
        assert_eq!(options_from_cli(&flagged).seed, 111);
    }

    #[test]
    fn cli_parses_max_memory_mb() {
        let cli = Cli::try_parse_from(["graynfs", "100160063", "--max-memory-mb", "512"])
            .expect("parse with memory cap");
        assert_eq!(cli.max_memory_mb, Some(512));
        let options = options_from_cli(&cli);
        assert_eq!(options.max_memory_mb, Some(512));
        let bare = Cli::try_parse_from(["graynfs", "100160063"]).expect("parse without memory cap");
        assert_eq!(bare.max_memory_mb, None);
        let bare_options = options_from_cli(&bare);
        assert_eq!(bare_options.max_memory_mb, None);
    }

    #[test]
    fn cli_parses_max_cpu_seconds() {
        let cli = Cli::try_parse_from(["graynfs", "100160063", "--max-cpu-seconds", "30"])
            .expect("parse with cpu cap");
        assert_eq!(cli.max_cpu_seconds, Some(30));
        let options = options_from_cli(&cli);
        assert_eq!(options.max_cpu_seconds, Some(30));
        assert_eq!(options.cpu_budget_ms(), Some(30_000));
        let bare = Cli::try_parse_from(["graynfs", "100160063"]).expect("parse without cpu cap");
        assert_eq!(bare.max_cpu_seconds, None);
        let bare_options = options_from_cli(&bare);
        assert_eq!(bare_options.max_cpu_seconds, None);
        assert_eq!(bare_options.cpu_budget_ms(), None);
    }

    #[test]
    fn batch_exit_code_truth_table() {
        assert_eq!(batch_exit_code(false, false), ExitCode::SUCCESS);
        assert_eq!(batch_exit_code(true, false), ExitCode::FAILURE);
        assert_eq!(batch_exit_code(false, true), ExitCode::FAILURE);
        assert_eq!(batch_exit_code(true, true), ExitCode::FAILURE);
    }

    #[test]
    fn batch_precancelled_token_aborts_after_first_record() {
        let options = FactorOptions::default();
        let token = CancellationToken::new();
        token.cancel();
        let reader = std::io::BufReader::new(std::io::Cursor::new("15\n21\n33\n"));
        let mut records = Vec::new();
        let completed = run_batch_lines(reader, &options, &token, false, true, &mut |record| {
            records.push(batch_record_value(record));
        })
        .expect("streaming batch reads");
        assert!(!completed);
        assert_eq!(records.len(), 1);
    }
}
