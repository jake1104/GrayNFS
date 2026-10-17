use graynfs_core::{
    FactorError, FactorOptions, factor, factor_with_progress, parse_decimal, probable_prime,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn to_py_error(error: graynfs_core::FactorError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[allow(clippy::too_many_arguments)]
fn build_options(
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    allow_fallback: bool,
    use_fermat: bool,
    use_squfof: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
) -> FactorOptions {
    FactorOptions {
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        max_duration: timeout_seconds.map(std::time::Duration::from_secs),
        allow_fallback,
        use_fermat,
        use_squfof,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
        ..FactorOptions::default()
    }
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, allow_fallback=true, use_fermat=true, use_squfof=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false))]
fn factor_decimal(
    py: Python<'_>,
    n: &str,
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    allow_fallback: bool,
    use_fermat: bool,
    use_squfof: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
) -> PyResult<Vec<String>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        allow_fallback,
        use_fermat,
        use_squfof,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
    );
    let report = py
        .detach(|| factor(&value, &options))
        .map_err(to_py_error)?;
    Ok(report.factors)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, callback=None, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, allow_fallback=true, use_fermat=true, use_squfof=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false))]
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
    allow_fallback: bool,
    use_fermat: bool,
    use_squfof: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
) -> PyResult<Vec<String>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        allow_fallback,
        use_fermat,
        use_squfof,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
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
                            .map_err(|error| FactorError::Factorization(error.to_string()))?;
                        Ok(())
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
            allow_fallback,
            use_fermat,
            use_squfof,
            use_pm1,
            use_pp1,
            use_ecm,
            use_mpqs,
            threads,
            perf_mode,
        )
    }
}

#[pyfunction]
fn is_probable_prime_decimal(py: Python<'_>, n: &str) -> PyResult<bool> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    Ok(py.detach(|| probable_prime(&value)))
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (n, degree=None, rational_base_bound=2000, algebraic_base_bound=2000, sieve_a_bound=700, sieve_b_bound=120, timeout_seconds=None, allow_fallback=true, use_fermat=true, use_squfof=true, use_pm1=true, use_pp1=true, use_ecm=true, use_mpqs=true, threads=None, perf_mode=false))]
fn stage_timings_decimal(
    py: Python<'_>,
    n: &str,
    degree: Option<usize>,
    rational_base_bound: usize,
    algebraic_base_bound: usize,
    sieve_a_bound: usize,
    sieve_b_bound: usize,
    timeout_seconds: Option<u64>,
    allow_fallback: bool,
    use_fermat: bool,
    use_squfof: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: Option<usize>,
    perf_mode: bool,
) -> PyResult<Vec<Py<PyAny>>> {
    let value = parse_decimal(n).map_err(to_py_error)?;
    let options = build_options(
        degree,
        rational_base_bound,
        algebraic_base_bound,
        sieve_a_bound,
        sieve_b_bound,
        timeout_seconds,
        allow_fallback,
        use_fermat,
        use_squfof,
        use_pm1,
        use_pp1,
        use_ecm,
        use_mpqs,
        threads,
        perf_mode,
    );
    let report = py
        .detach(|| factor(&value, &options))
        .map_err(to_py_error)?;
    let entries = report.timings.into_iter().map(|entry| {
        let dict = PyDict::new(py);
        dict.set_item("stage", entry.stage).unwrap();
        dict.set_item("degree", entry.degree).unwrap();
        dict.set_item("elapsed_ms", entry.elapsed_ms).unwrap();
        dict.set_item("detail", entry.detail).unwrap();
        dict.into_any().unbind()
    });
    Ok(entries.collect())
}

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn _graynfs(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(factor_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(factor_with_progress_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(stage_timings_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(is_probable_prime_decimal, module)?)?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    Ok(())
}
