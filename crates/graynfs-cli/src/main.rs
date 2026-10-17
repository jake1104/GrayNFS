use clap::Parser;
use graynfs_core::{
    FactorError, FactorOptions, FactorReport, StageTiming, factor_with_timings, parse_decimal,
};
use std::process::ExitCode;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser, Debug)]
#[command(name = "graynfs", version, about = "Factor integers with GNFS")]
struct Cli {
    number: String,
    #[arg(
        long,
        help = "Force a fixed polynomial degree instead of the automatic race"
    )]
    degree: Option<usize>,
    #[arg(long, default_value_t = 2_000)]
    rational_base_bound: usize,
    #[arg(long, default_value_t = 2_000)]
    algebraic_base_bound: usize,
    #[arg(long, default_value_t = 700)]
    sieve_a_bound: usize,
    #[arg(long, default_value_t = 120)]
    sieve_b_bound: usize,
    #[arg(long, default_value_t = 50_000)]
    large_prime_bound: usize,
    #[arg(long, default_value_t = 60)]
    timeout_seconds: u64,
    #[arg(
        long,
        help = "Cap estimated sieve and matrix memory in MiB before allocating (default 1024)"
    )]
    max_memory_mb: Option<usize>,
    #[arg(long)]
    no_fallback: bool,
    #[arg(long)]
    skip_preprocessing: bool,
    #[arg(long, help = "Print a per-stage timing breakdown")]
    timings: bool,
    #[arg(
        long,
        help = "Allow the general siever to keep relations with two or three large primes and pair them cyclically"
    )]
    double_large_primes: bool,
    #[arg(long)]
    no_fermat: bool,
    #[arg(long)]
    no_squfof: bool,
    #[arg(long)]
    no_pm1: bool,
    #[arg(long)]
    no_pp1: bool,
    #[arg(long)]
    no_ecm: bool,
    #[arg(long)]
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
    #[arg(long, help = "Emit the result as a single JSON object on stdout")]
    json: bool,
}

fn options_from_cli(cli: &Cli) -> FactorOptions {
    FactorOptions {
        degree: cli.degree,
        rational_base_bound: cli.rational_base_bound,
        algebraic_base_bound: cli.algebraic_base_bound,
        sieve_a_bound: cli.sieve_a_bound,
        sieve_b_bound: cli.sieve_b_bound,
        large_prime_bound: cli.large_prime_bound,
        max_duration: Some(std::time::Duration::from_secs(cli.timeout_seconds)),
        max_memory_mb: cli.max_memory_mb,
        allow_fallback: !cli.no_fallback,
        skip_preprocessing: cli.skip_preprocessing,
        enable_double_large_primes: cli.double_large_primes,
        use_fermat: !cli.no_fermat,
        use_squfof: !cli.no_squfof,
        use_pm1: !cli.no_pm1,
        use_pp1: !cli.no_pp1,
        use_ecm: !cli.no_ecm,
        use_mpqs: !cli.no_mpqs,
        threads: cli.threads,
        perf_mode: cli.perf_mode,
        ..FactorOptions::default()
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let number = match parse_decimal(&cli.number) {
        Ok(value) => value,
        Err(error) => {
            if cli.json {
                print_json(&json_error_value(&error));
            } else {
                report_error(&error);
            }
            return ExitCode::FAILURE;
        }
    };
    let options = options_from_cli(&cli);
    let (result, timings) = factor_with_timings(&number, &options, |_| Ok(()));
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
                print_json(&json_error_value(&error));
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

fn json_report_value(report: &FactorReport, timings: Option<&[StageTiming]>) -> serde_json::Value {
    let mut value = serde_json::json!({
        "input": report.input,
        "method": report.method,
        "degree": report.degree,
        "factors": report.factors,
        "relations": report.relations,
        "rank": report.rank,
        "elapsed_ms": report.elapsed_ms,
    });
    if let Some(entries) = timings {
        let array = entries
            .iter()
            .map(|entry| {
                serde_json::json!({
                    "stage": entry.stage,
                    "degree": entry.degree,
                    "elapsed_ms": entry.elapsed_ms,
                    "detail": entry.detail,
                })
            })
            .collect::<Vec<_>>();
        value["timings"] = serde_json::Value::Array(array);
    }
    value
}

fn json_error_value(error: &FactorError) -> serde_json::Value {
    let (category, detail) = error_parts(error);
    serde_json::json!({
        "error": {
            "category": category,
            "detail": detail,
        }
    })
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
        assert_eq!(value["relations"], 604);
        assert_eq!(value["rank"], 293);
        assert_eq!(value["elapsed_ms"], 120);
        assert!(value.get("timings").is_none());
        let text = value.to_string();
        let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed, value);
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
    fn json_report_serializes_nullable_degree() {
        let mut report = sample_report();
        report.degree = None;
        let value = json_report_value(&report, None);
        assert!(value["degree"].is_null());
    }

    #[test]
    fn json_error_carries_category_and_detail() {
        let value = json_error_value(&FactorError::ResourceLimit("deadline".to_string()));
        assert_eq!(value["error"]["category"], "resource-limit");
        assert_eq!(value["error"]["detail"], "deadline");
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
            let value = json_error_value(error);
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
}
