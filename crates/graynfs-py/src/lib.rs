use graynfs_core::{
    FactorError, FactorOptions, StageTiming, factor, factor_with_progress, factor_with_timings,
    parse_decimal, probable_prime, probable_prime_fast,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn error_category(error: &FactorError) -> &'static str {
    match error {
        FactorError::InvalidInput(_) => "input",
        FactorError::ResourceLimit(_) => "resource-limit",
        FactorError::Cancelled => "cancelled",
        FactorError::Checkpoint(_) => "checkpoint",
        _ => "internal",
    }
}

fn error_detail(error: &FactorError) -> String {
    match error {
        FactorError::InvalidInput(detail)
        | FactorError::ResourceLimit(detail)
        | FactorError::Checkpoint(detail) => detail.clone(),
        FactorError::Cancelled => "the operation was cancelled".to_string(),
        other => other.to_string(),
    }
}

fn py_error_message(error: &FactorError) -> String {
    format!("{}: {}", error_category(error), error_detail(error))
}

fn to_py_error(error: graynfs_core::FactorError) -> PyErr {
    PyValueError::new_err(py_error_message(&error))
}

fn stage_timing_to_py(py: Python<'_>, entry: StageTiming) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("stage", entry.stage)?;
    dict.set_item("degree", entry.degree)?;
    dict.set_item("elapsed_ms", entry.elapsed_ms)?;
    dict.set_item("detail", entry.detail)?;
    Ok(dict.unbind())
}

#[allow(clippy::too_many_arguments)]
fn build_options(
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    max_memory_mb: Option<usize>,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: usize,
    enable_double_large_primes: bool,
    seed: Option<u64>,
    max_cpu_seconds: Option<u64>,
    use_snfs: bool,
    checkpoint_path: Option<String>,
    resume_path: Option<String>,
) -> FactorOptions {
    FactorOptions {
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        large_prime_bound,
        seed: seed.unwrap_or(FactorOptions::default().seed),
        threads,
        max_duration: timeout_seconds.map(std::time::Duration::from_secs),
        max_memory_mb,
        allow_fallback,
        skip_preprocessing,
        enable_double_large_primes,
        use_fermat,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        perf_mode,
        max_cpu_seconds,
        use_snfs,
        checkpoint_path: checkpoint_path.map(std::path::PathBuf::from),
        resume_path: resume_path.map(std::path::PathBuf::from),
        ..FactorOptions::default()
    }
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, max_memory_mb=None, allow_fallback=true, use_fermat=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false, skip_preprocessing=false, large_prime_bound=50000, enable_double_large_primes=false, seed=None, max_cpu_seconds=None, use_snfs=true, checkpoint_path=None, resume_path=None))]
fn factor_decimal(
    py: Python<'_>,
    n: &str,
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    max_memory_mb: Option<usize>,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: usize,
    enable_double_large_primes: bool,
    seed: Option<u64>,
    max_cpu_seconds: Option<u64>,
    use_snfs: bool,
    checkpoint_path: Option<String>,
    resume_path: Option<String>,
) -> PyResult<Vec<String>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        max_memory_mb,
        allow_fallback,
        use_fermat,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
        skip_preprocessing,
        large_prime_bound,
        enable_double_large_primes,
        seed,
        max_cpu_seconds,
        use_snfs,
        checkpoint_path,
        resume_path,
    );
    let report = py
        .detach(|| factor(&value, &options))
        .map_err(to_py_error)?;
    Ok(report.factors)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, callback=None, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, max_memory_mb=None, allow_fallback=true, use_fermat=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false, skip_preprocessing=false, large_prime_bound=50000, enable_double_large_primes=false, seed=None, max_cpu_seconds=None, use_snfs=true, checkpoint_path=None, resume_path=None))]
fn factor_with_progress_decimal(
    py: Python<'_>,
    n: &str,
    callback: Option<Py<PyAny>>,
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    max_memory_mb: Option<usize>,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: usize,
    enable_double_large_primes: bool,
    seed: Option<u64>,
    max_cpu_seconds: Option<u64>,
    use_snfs: bool,
    checkpoint_path: Option<String>,
    resume_path: Option<String>,
) -> PyResult<Vec<String>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        max_memory_mb,
        allow_fallback,
        use_fermat,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
        skip_preprocessing,
        large_prime_bound,
        enable_double_large_primes,
        seed,
        max_cpu_seconds,
        use_snfs,
        checkpoint_path.clone(),
        resume_path.clone(),
    );
    if let Some(callback) = callback {
        let report = py
            .detach(|| {
                factor_with_progress(&value, &options, |event| {
                    Python::attach(|py| {
                        let event_dict = PyDict::new(py);
                        event_dict
                            .set_item("stage", &event.stage)
                            .map_err(|error| FactorError::Factorization(error.to_string()))?;
                        event_dict
                            .set_item("current", event.current)
                            .map_err(|error| FactorError::Factorization(error.to_string()))?;
                        event_dict
                            .set_item("total", event.total)
                            .map_err(|error| FactorError::Factorization(error.to_string()))?;
                        event_dict
                            .set_item("message", &event.message)
                            .map_err(|error| FactorError::Factorization(error.to_string()))?;
                        callback
                            .call1(py, (event_dict,))
                            .map_err(|error| FactorError::Factorization(error.to_string()))
                            .and_then(|outcome| {
                                if outcome.is_none(py)
                                    || outcome.is_truthy(py).map_err(|error| {
                                        FactorError::Factorization(error.to_string())
                                    })?
                                {
                                    Ok(())
                                } else {
                                    Err(FactorError::Cancelled)
                                }
                            })
                    })
                })
            })
            .map_err(to_py_error)?;
        Ok(report.factors)
    } else {
        factor_decimal(
            py,
            n,
            degree,
            rational_base_bound,
            algebraic_base_bound,
            sieve_a_bound,
            sieve_b_bound,
            timeout_seconds,
            max_memory_mb,
            allow_fallback,
            use_fermat,
            use_pm1,
            use_pp1,
            use_ecm,
            use_mpqs,
            threads,
            perf_mode,
            skip_preprocessing,
            large_prime_bound,
            enable_double_large_primes,
            seed,
            max_cpu_seconds,
            use_snfs,
            checkpoint_path,
            resume_path,
        )
    }
}

#[pyfunction]
fn is_probable_prime_decimal(py: Python<'_>, n: &str) -> PyResult<bool> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    Ok(py.detach(|| probable_prime(&value)))
}

#[pyfunction]
#[pyo3(signature = (n, rounds=64))]
fn is_probable_prime_fast_decimal(py: Python<'_>, n: &str, rounds: u64) -> PyResult<bool> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    Ok(py.detach(|| probable_prime_fast(&value, rounds)))
}

#[pyfunction]
#[pyo3(signature = (n, method="aks"))]
fn is_prime_proof_decimal(py: Python<'_>, n: &str, method: &str) -> PyResult<bool> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    py.detach(|| match method {
        "bpsw" => Ok(graynfs_core::prime::is_baillie_psw_prime(&value)),
        "frobenius" => Ok(graynfs_core::prime::is_frobenius_probable_prime(&value)),
        "aks" => Ok(graynfs_core::prime::is_prime_aks_proof(&value)),
        "bpsw_frobenius" => {
            let options_frobenius = FactorOptions {
                use_frobenius: true,
                ..FactorOptions::default()
            };
            Ok(graynfs_core::prime::is_probable_prime_with_options(
                &value,
                &options_frobenius,
            ))
        }
        _ => Err(PyValueError::new_err(format!(
            "input: unknown proof method: {method}"
        ))),
    })
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, max_memory_mb=None, allow_fallback=true, use_fermat=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false, skip_preprocessing=false, large_prime_bound=50000, enable_double_large_primes=false, seed=None, max_cpu_seconds=None, use_snfs=true, checkpoint_path=None, resume_path=None))]
fn stage_timings_decimal(
    py: Python<'_>,
    n: &str,
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    max_memory_mb: Option<usize>,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: usize,
    enable_double_large_primes: bool,
    seed: Option<u64>,
    max_cpu_seconds: Option<u64>,
    use_snfs: bool,
    checkpoint_path: Option<String>,
    resume_path: Option<String>,
) -> PyResult<Vec<Py<PyAny>>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        max_memory_mb,
        allow_fallback,
        use_fermat,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
        skip_preprocessing,
        large_prime_bound,
        enable_double_large_primes,
        seed,
        max_cpu_seconds,
        use_snfs,
        checkpoint_path,
        resume_path,
    );
    let report = py
        .detach(|| factor(&value, &options))
        .map_err(to_py_error)?;
    report
        .timings
        .into_iter()
        .map(|entry| Ok(stage_timing_to_py(py, entry)?.into_any()))
        .collect()
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, max_memory_mb=None, allow_fallback=true, use_fermat=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false, skip_preprocessing=false, large_prime_bound=50000, enable_double_large_primes=false, seed=None, max_cpu_seconds=None, use_snfs=true, checkpoint_path=None, resume_path=None))]
fn factor_report_decimal(
    py: Python<'_>,
    n: &str,
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    max_memory_mb: Option<usize>,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: usize,
    enable_double_large_primes: bool,
    seed: Option<u64>,
    max_cpu_seconds: Option<u64>,
    use_snfs: bool,
    checkpoint_path: Option<String>,
    resume_path: Option<String>,
) -> PyResult<Py<PyDict>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        max_memory_mb,
        allow_fallback,
        use_fermat,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
        skip_preprocessing,
        large_prime_bound,
        enable_double_large_primes,
        seed,
        max_cpu_seconds,
        use_snfs,
        checkpoint_path,
        resume_path,
    );
    let (result, timings) = py.detach(|| factor_with_timings(&value, &options, |_| Ok(())));
    let report = result.map_err(to_py_error)?;
    let dict = PyDict::new(py);
    dict.set_item("factors", report.factors)?;
    dict.set_item("factor_confidence", report.factor_confidence)?;
    dict.set_item("polynomial", report.polynomial)?;
    dict.set_item("root", report.root)?;
    dict.set_item("method", report.method)?;
    dict.set_item("degree", report.degree)?;
    dict.set_item("relations", report.relations)?;
    dict.set_item("rank", report.rank)?;
    dict.set_item("elapsed_ms", report.elapsed_ms)?;
    let entries: PyResult<Vec<Py<PyDict>>> = timings
        .into_iter()
        .map(|entry| stage_timing_to_py(py, entry))
        .collect();
    dict.set_item("timings", entries?)?;
    Ok(dict.unbind())
}

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_options(max_memory_mb: Option<usize>) -> FactorOptions {
        build_options(
            None,
            2000,
            2000,
            700,
            120,
            None,
            max_memory_mb,
            true,
            true,
            true,
            true,
            true,
            true,
            None,
            false,
            false,
            50_000,
            false,
            None,
            None,
            true,
            None,
            None,
        )
    }

    fn checkpoint_options() -> FactorOptions {
        build_options(
            None,
            2000,
            2000,
            700,
            120,
            None,
            None,
            true,
            true,
            true,
            true,
            true,
            true,
            None,
            false,
            false,
            50_000,
            false,
            None,
            None,
            true,
            Some("collect.json".to_string()),
            Some("resume.json".to_string()),
        )
    }

    #[test]
    fn build_options_carries_memory_budget() {
        let options = base_options(Some(512));
        assert_eq!(options.max_memory_mb, Some(512));
        assert_eq!(options.memory_budget_bytes(), 512 * 1024 * 1024);
    }

    #[test]
    fn build_options_defaults_to_unbounded_budget() {
        let options = base_options(None);
        assert_eq!(options.max_memory_mb, None);
    }

    #[test]
    fn build_options_enables_snfs_path() {
        assert!(base_options(None).use_snfs);
    }

    #[test]
    fn build_options_carries_cpu_budget() {
        let capped = build_options(
            None,
            2000,
            2000,
            700,
            120,
            None,
            None,
            true,
            true,
            true,
            true,
            true,
            true,
            None,
            false,
            false,
            50_000,
            false,
            None,
            Some(30),
            true,
            None,
            None,
        );
        assert_eq!(capped.max_cpu_seconds, Some(30));
        assert_eq!(capped.cpu_budget_ms(), Some(30_000));
        let open = base_options(None);
        assert_eq!(open.max_cpu_seconds, None);
        assert_eq!(open.cpu_budget_ms(), None);
    }

    #[test]
    fn build_options_carries_checkpoint_paths() {
        let options = checkpoint_options();
        assert_eq!(
            options.checkpoint_path,
            Some(std::path::PathBuf::from("collect.json"))
        );
        assert_eq!(
            options.resume_path,
            Some(std::path::PathBuf::from("resume.json"))
        );
        let plain = base_options(None);
        assert_eq!(plain.checkpoint_path, None);
        assert_eq!(plain.resume_path, None);
    }

    #[test]
    fn error_category_maps_checkpoint() {
        let error = FactorError::Checkpoint("checkpoint checksum mismatch".to_string());
        assert_eq!(error_category(&error), "checkpoint");
        assert!(py_error_message(&error).starts_with("checkpoint: "));
    }

    #[test]
    fn error_category_maps_input() {
        let error = FactorError::InvalidInput("bad digit".to_string());
        assert_eq!(error_category(&error), "input");
        assert!(py_error_message(&error).starts_with("input: "));
    }

    #[test]
    fn error_category_maps_resource_limit() {
        let error = FactorError::ResourceLimit("deadline".to_string());
        assert_eq!(error_category(&error), "resource-limit");
        assert!(py_error_message(&error).starts_with("resource-limit: "));
    }

    #[test]
    fn error_category_maps_cancelled() {
        let error = FactorError::Cancelled;
        assert_eq!(error_category(&error), "cancelled");
        assert!(py_error_message(&error).starts_with("cancelled: "));
    }

    #[test]
    fn error_category_maps_internal_variants() {
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
            assert_eq!(error_category(error), "internal");
            assert!(py_error_message(error).starts_with("internal: "));
        }
    }

    #[test]
    fn py_error_message_carries_detail() {
        let error = FactorError::InvalidInput("invalid decimal integer: abc".to_string());
        assert_eq!(
            py_error_message(&error),
            "input: invalid decimal integer: abc"
        );
    }
}

#[pymodule]
fn _graynfs(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(factor_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(factor_with_progress_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(stage_timings_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(factor_report_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(is_probable_prime_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(is_probable_prime_fast_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(is_prime_proof_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    Ok(())
}
