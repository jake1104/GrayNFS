use std::time::Duration;

use graynfs_core::{FactorOptions, FactorReport, StageTiming, factor_with_timings, parse_decimal};
use serde_json::{Value, json};
use wasm_bindgen::prelude::*;

const MIN_BUDGET_MS: u32 = 1_000;
const MAX_BUDGET_MS: u32 = 10_000;
const MEMORY_LIMIT_MB: usize = 512;
const SEED: u64 = 0x9e3779b97f4a7c15;

#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[wasm_bindgen]
pub fn factorize(decimal: &str, budget_ms: u32) -> String {
    let budget = budget_ms.clamp(MIN_BUDGET_MS, MAX_BUDGET_MS);
    let number = match parse_decimal(decimal) {
        Ok(number) => number,
        Err(error) => return failure(decimal, &error.to_string()),
    };
    let options = FactorOptions {
        threads: Some(1),
        max_duration: Some(Duration::from_millis(u64::from(budget))),
        max_memory_mb: Some(MEMORY_LIMIT_MB),
        seed: SEED,
        ..FactorOptions::default()
    };
    let (result, timings) = factor_with_timings(&number, &options, |_| Ok(()));
    match result {
        Ok(report) => success(decimal, &report, &timings),
        Err(error) => failure(decimal, &error.to_string()),
    }
}

fn success(input: &str, report: &FactorReport, timings: &[StageTiming]) -> String {
    let stages: Vec<Value> = timings
        .iter()
        .map(|timing| {
            json!({
                "stage": &timing.stage,
                "degree": timing.degree,
                "ms": millis(timing.elapsed_ms),
                "detail": &timing.detail,
            })
        })
        .collect();
    json!({
        "ok": true,
        "input": input,
        "factors": &report.factors,
        "factorConfidence": &report.factor_confidence,
        "method": &report.method,
        "degree": report.degree,
        "polynomial": &report.polynomial,
        "root": &report.root,
        "relations": report.relations,
        "rank": report.rank,
        "elapsedMs": millis(report.elapsed_ms),
        "stages": stages,
    })
    .to_string()
}

fn failure(input: &str, error: &str) -> String {
    json!({
        "ok": false,
        "input": input,
        "error": error,
    })
    .to_string()
}

fn millis(value: u128) -> u64 {
    value.min(u128::from(u64::MAX)) as u64
}
