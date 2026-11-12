from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Optional

from . import _graynfs

__all__ = [
    "GrayNFSError",
    "factor",
    "factor_batch",
    "factor_report",
    "factor_with_progress",
    "is_prime_proof",
    "is_probable_prime",
    "is_probable_prime_fast",
    "stage_timings",
    "version",
]

ProgressCallback = Callable[[dict[str, Any]], Optional[bool]]


class GrayNFSError(ValueError):
    def __init__(self, category: str, detail: str) -> None:
        self.category = category
        self.detail = detail
        super().__init__(category + ": " + detail)


_GRAY_ERROR_PREFIXES = ("input: ", "resource-limit: ", "cancelled: ", "internal: ")


def _to_gray_error(exc: ValueError) -> ValueError:
    message = str(exc)
    for prefix in _GRAY_ERROR_PREFIXES:
        if message.startswith(prefix):
            category, _, detail = message.partition(": ")
            return GrayNFSError(category, detail)
    return exc


def _validate_number(value: int) -> str:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("n must be an integer")
    if value < 1:
        raise ValueError("input: n must be a positive integer")
    return str(value)


def _validate_seed(seed: int | None) -> None:
    if seed is None:
        return
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("seed must be an integer or None")
    if seed < 0:
        raise ValueError("input: seed must be a non-negative integer")


def factor(
    n: int,
    *,
    degree: int | None = None,
    rational_base_bound: int = 2000,
    algebraic_base_bound: int = 2000,
    sieve_a_bound: int = 700,
    sieve_b_bound: int = 120,
    timeout_seconds: int | None = None,
    max_memory_mb: int | None = None,
    allow_fallback: bool = True,
    use_fermat: bool = True,
    use_pm1: bool = True,
    use_pp1: bool = True,
    use_ecm: bool = True,
    use_mpqs: bool = True,
    threads: int | None = None,
    perf_mode: bool = False,
    skip_preprocessing: bool = False,
    large_prime_bound: int = 50000,
    enable_double_large_primes: bool = False,
    seed: int | None = None,
    max_cpu_seconds: int | None = None,
    use_snfs: bool = True,
) -> list[int]:
    """Factor n into ascending prime factors; keywords mirror the documented surface (None means auto or unset, True/False toggles keep or skip a stage)."""
    try:
        _validate_seed(seed)
        return [
            int(value)
            for value in _graynfs.factor_decimal(
                _validate_number(n),
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
            )
        ]
    except ValueError as exc:
        raise _to_gray_error(exc)


def factor_with_progress(
    n: int,
    callback: ProgressCallback,
    *,
    degree: int | None = None,
    rational_base_bound: int = 2000,
    algebraic_base_bound: int = 2000,
    sieve_a_bound: int = 700,
    sieve_b_bound: int = 120,
    timeout_seconds: int | None = None,
    max_memory_mb: int | None = None,
    allow_fallback: bool = True,
    use_fermat: bool = True,
    use_pm1: bool = True,
    use_pp1: bool = True,
    use_ecm: bool = True,
    use_mpqs: bool = True,
    threads: int | None = None,
    perf_mode: bool = False,
    skip_preprocessing: bool = False,
    large_prime_bound: int = 50000,
    enable_double_large_primes: bool = False,
    seed: int | None = None,
    max_cpu_seconds: int | None = None,
    use_snfs: bool = True,
) -> list[int]:
    """Factor n like factor, but call callback({stage, current, total, message}) at each stage boundary; the callback is required here (the raw binding accepts None and falls back to a plain factor). Return a falsy value (e.g. False) from the callback to cancel the operation (surfaced as GrayNFSError with category=cancelled); truthy values and None continue."""
    if not callable(callback):
        raise TypeError("callback must be callable")
    try:
        _validate_seed(seed)
        return [
            int(value)
            for value in _graynfs.factor_with_progress_decimal(
                _validate_number(n),
                callback,
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
            )
        ]
    except ValueError as exc:
        raise _to_gray_error(exc)


def is_probable_prime(n: int) -> bool:
    """BPSW probable-prime gate (deterministic for u64, probable otherwise; not a proof)."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 0:
        return False
    return bool(_graynfs.is_probable_prime_decimal(str(n)))


def is_probable_prime_fast(n: int, rounds: int = 64) -> bool:
    """Miller-Rabin fast filter with the given round count; probable only, never a proof."""
    try:
        if isinstance(n, bool) or not isinstance(n, int):
            raise TypeError("n must be an integer")
        if isinstance(rounds, bool) or not isinstance(rounds, int):
            raise TypeError("rounds must be an integer")
        if rounds < 0:
            raise ValueError("input: rounds must be a non-negative integer")
        return bool(_graynfs.is_probable_prime_fast_decimal(str(n), rounds))
    except ValueError as exc:
        raise _to_gray_error(exc)


def is_prime_proof(n: int, method: str = "aks") -> bool:
    """Deterministic proof by default (method="aks", 64-bit cap, off the hot factoring path); "bpsw", "frobenius", and "bpsw_frobenius" are probable-only (probable, not a proof); over-64-bit False means unproven, not composite."""
    try:
        if isinstance(n, bool) or not isinstance(n, int):
            raise TypeError("n must be an integer")
        if not isinstance(method, str):
            raise TypeError("method must be a string")
        return bool(_graynfs.is_prime_proof_decimal(str(n), method))
    except ValueError as exc:
        raise _to_gray_error(exc)


def stage_timings(
    n: int,
    *,
    degree: int | None = None,
    rational_base_bound: int = 2000,
    algebraic_base_bound: int = 2000,
    sieve_a_bound: int = 700,
    sieve_b_bound: int = 120,
    timeout_seconds: int | None = None,
    max_memory_mb: int | None = None,
    allow_fallback: bool = True,
    use_fermat: bool = True,
    use_pm1: bool = True,
    use_pp1: bool = True,
    use_ecm: bool = True,
    use_mpqs: bool = True,
    threads: int | None = None,
    perf_mode: bool = False,
    skip_preprocessing: bool = False,
    large_prime_bound: int = 50000,
    enable_double_large_primes: bool = False,
    seed: int | None = None,
    max_cpu_seconds: int | None = None,
    use_snfs: bool = True,
) -> list[dict[str, Any]]:
    """Return per-stage timing records [{stage, degree, elapsed_ms, detail}] for factoring n without printing."""
    try:
        _validate_seed(seed)
        return [
            {
                "stage": str(entry["stage"]),
                "degree": None if entry["degree"] is None else int(entry["degree"]),
                "elapsed_ms": int(entry["elapsed_ms"]),
                "detail": str(entry["detail"]),
            }
            for entry in _graynfs.stage_timings_decimal(
                _validate_number(n),
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
            )
        ]
    except ValueError as exc:
        raise _to_gray_error(exc)


def version() -> str:
    return str(_graynfs.version())


def _validate_max_workers(max_workers: int | None) -> int | None:
    if max_workers is None:
        return None
    if isinstance(max_workers, bool) or not isinstance(max_workers, int):
        raise TypeError("max_workers must be an integer or None")
    if max_workers < 1:
        raise ValueError("input: max_workers must be a positive integer")
    return max_workers


def _batch_slots(
    *,
    max_memory_mb: int | None,
    rational_base_bound: int,
    algebraic_base_bound: int,
    sieve_a_bound: int,
    per_job_threads: int,
    max_workers: int,
) -> int:
    budget_bytes = (1024 if max_memory_mb is None else max_memory_mb) * 1024 * 1024
    per_job_bytes = (rational_base_bound + algebraic_base_bound) * 64 + (
        2 * sieve_a_bound + 1
    ) * 68 * max(1, per_job_threads)
    if per_job_bytes <= 0:
        return max_workers
    if per_job_bytes > budget_bytes:
        return 1
    return max(1, min(max_workers, budget_bytes // per_job_bytes))


def factor_batch(
    ns: list[int],
    *,
    fail_fast: bool = True,
    degree: int | None = None,
    rational_base_bound: int = 2000,
    algebraic_base_bound: int = 2000,
    sieve_a_bound: int = 700,
    sieve_b_bound: int = 120,
    timeout_seconds: int | None = None,
    max_memory_mb: int | None = None,
    allow_fallback: bool = True,
    use_fermat: bool = True,
    use_pm1: bool = True,
    use_pp1: bool = True,
    use_ecm: bool = True,
    use_mpqs: bool = True,
    threads: int | None = None,
    perf_mode: bool = False,
    skip_preprocessing: bool = False,
    large_prime_bound: int = 50000,
    enable_double_large_primes: bool = False,
    seed: int | None = None,
    max_cpu_seconds: int | None = None,
    use_snfs: bool = True,
    max_workers: int | None = None,
) -> list[list[int] | GrayNFSError]:
    """Factor each input in ns with the shared factor keywords; fail_fast=True raises on the first failure while False records the GrayNFSError in place. max_workers=None (default) factors sequentially; an integer selects that many bounded parallel workers while preserving input order, dividing threads across jobs and queuing jobs so total memory stays within the memory budget."""
    try:
        resolved_workers = _validate_max_workers(max_workers)
    except ValueError as exc:
        raise _to_gray_error(exc)
    if resolved_workers is None or resolved_workers == 1:
        return _factor_batch_sequential(
            ns,
            fail_fast=fail_fast,
            degree=degree,
            rational_base_bound=rational_base_bound,
            algebraic_base_bound=algebraic_base_bound,
            sieve_a_bound=sieve_a_bound,
            sieve_b_bound=sieve_b_bound,
            timeout_seconds=timeout_seconds,
            max_memory_mb=max_memory_mb,
            allow_fallback=allow_fallback,
            use_fermat=use_fermat,
            use_pm1=use_pm1,
            use_pp1=use_pp1,
            use_ecm=use_ecm,
            use_mpqs=use_mpqs,
            threads=threads,
            perf_mode=perf_mode,
            skip_preprocessing=skip_preprocessing,
            large_prime_bound=large_prime_bound,
            enable_double_large_primes=enable_double_large_primes,
            seed=seed,
            max_cpu_seconds=max_cpu_seconds,
            use_snfs=use_snfs,
        )
    return _factor_batch_parallel(
        ns,
        fail_fast=fail_fast,
        degree=degree,
        rational_base_bound=rational_base_bound,
        algebraic_base_bound=algebraic_base_bound,
        sieve_a_bound=sieve_a_bound,
        sieve_b_bound=sieve_b_bound,
        timeout_seconds=timeout_seconds,
        max_memory_mb=max_memory_mb,
        allow_fallback=allow_fallback,
        use_fermat=use_fermat,
        use_pm1=use_pm1,
        use_pp1=use_pp1,
        use_ecm=use_ecm,
        use_mpqs=use_mpqs,
        threads=threads,
        perf_mode=perf_mode,
        skip_preprocessing=skip_preprocessing,
        large_prime_bound=large_prime_bound,
        enable_double_large_primes=enable_double_large_primes,
        seed=seed,
        max_cpu_seconds=max_cpu_seconds,
        use_snfs=use_snfs,
        max_workers=resolved_workers,
    )


def _factor_batch_sequential(
    ns: list[int],
    *,
    fail_fast: bool,
    degree: int | None,
    rational_base_bound: int,
    algebraic_base_bound: int,
    sieve_a_bound: int,
    sieve_b_bound: int,
    timeout_seconds: int | None,
    max_memory_mb: int | None,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: int | None,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: int,
    enable_double_large_primes: bool,
    seed: int | None,
    max_cpu_seconds: int | None,
    use_snfs: bool,
) -> list[list[int] | GrayNFSError]:
    results: list[list[int] | GrayNFSError] = []
    for value in ns:
        try:
            results.append(
                factor(
                    value,
                    degree=degree,
                    rational_base_bound=rational_base_bound,
                    algebraic_base_bound=algebraic_base_bound,
                    sieve_a_bound=sieve_a_bound,
                    sieve_b_bound=sieve_b_bound,
                    timeout_seconds=timeout_seconds,
                    max_memory_mb=max_memory_mb,
                    allow_fallback=allow_fallback,
                    use_fermat=use_fermat,
                    use_pm1=use_pm1,
                    use_pp1=use_pp1,
                    use_ecm=use_ecm,
                    use_mpqs=use_mpqs,
                    threads=threads,
                    perf_mode=perf_mode,
                    skip_preprocessing=skip_preprocessing,
                    large_prime_bound=large_prime_bound,
                    enable_double_large_primes=enable_double_large_primes,
                    seed=seed,
                    max_cpu_seconds=max_cpu_seconds,
                    use_snfs=use_snfs,
                )
            )
        except GrayNFSError as exc:
            if fail_fast:
                raise
            results.append(exc)
    return results


def _factor_batch_parallel(
    ns: list[int],
    *,
    fail_fast: bool,
    degree: int | None,
    rational_base_bound: int,
    algebraic_base_bound: int,
    sieve_a_bound: int,
    sieve_b_bound: int,
    timeout_seconds: int | None,
    max_memory_mb: int | None,
    allow_fallback: bool,
    use_fermat: bool,
    use_pm1: bool,
    use_pp1: bool,
    use_ecm: bool,
    use_mpqs: bool,
    threads: int | None,
    perf_mode: bool,
    skip_preprocessing: bool,
    large_prime_bound: int,
    enable_double_large_primes: bool,
    seed: int | None,
    max_cpu_seconds: int | None,
    use_snfs: bool,
    max_workers: int,
) -> list[list[int] | GrayNFSError]:
    values = list(ns)
    total_workers = threads if threads is not None else (os.cpu_count() or 1)
    per_job_threads = max(1, total_workers // max_workers)
    slots = _batch_slots(
        max_memory_mb=max_memory_mb,
        rational_base_bound=rational_base_bound,
        algebraic_base_bound=algebraic_base_bound,
        sieve_a_bound=sieve_a_bound,
        per_job_threads=per_job_threads,
        max_workers=max_workers,
    )

    def _run(value: int) -> list[int] | GrayNFSError:
        try:
            return factor(
                value,
                degree=degree,
                rational_base_bound=rational_base_bound,
                algebraic_base_bound=algebraic_base_bound,
                sieve_a_bound=sieve_a_bound,
                sieve_b_bound=sieve_b_bound,
                timeout_seconds=timeout_seconds,
                max_memory_mb=max_memory_mb,
                allow_fallback=allow_fallback,
                use_fermat=use_fermat,
                use_pm1=use_pm1,
                use_pp1=use_pp1,
                use_ecm=use_ecm,
                use_mpqs=use_mpqs,
                threads=per_job_threads,
                perf_mode=perf_mode,
                skip_preprocessing=skip_preprocessing,
                large_prime_bound=large_prime_bound,
                enable_double_large_primes=enable_double_large_primes,
                seed=seed,
                max_cpu_seconds=max_cpu_seconds,
                use_snfs=use_snfs,
            )
        except GrayNFSError as exc:
            return exc

    with ThreadPoolExecutor(max_workers=slots) as executor:
        futures = [executor.submit(_run, value) for value in values]
        results: list[list[int] | GrayNFSError] = []
        for future in futures:
            outcome = future.result()
            if isinstance(outcome, GrayNFSError):
                if fail_fast:
                    raise outcome
                results.append(outcome)
            else:
                results.append(outcome)
        return results


def factor_report(
    n: int,
    *,
    degree: int | None = None,
    rational_base_bound: int = 2000,
    algebraic_base_bound: int = 2000,
    sieve_a_bound: int = 700,
    sieve_b_bound: int = 120,
    timeout_seconds: int | None = None,
    max_memory_mb: int | None = None,
    allow_fallback: bool = True,
    use_fermat: bool = True,
    use_pm1: bool = True,
    use_pp1: bool = True,
    use_ecm: bool = True,
    use_mpqs: bool = True,
    threads: int | None = None,
    perf_mode: bool = False,
    skip_preprocessing: bool = False,
    large_prime_bound: int = 50000,
    enable_double_large_primes: bool = False,
    seed: int | None = None,
    max_cpu_seconds: int | None = None,
    use_snfs: bool = True,
) -> dict[str, Any]:
    """Factor n once and return {factors, factor_confidence, polynomial, root, method, degree, relations, rank, elapsed_ms, timings}."""
    try:
        _validate_seed(seed)
        raw = _graynfs.factor_report_decimal(
            _validate_number(n),
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
        )
    except ValueError as exc:
        raise _to_gray_error(exc)
    return {
        "factors": [int(value) for value in raw["factors"]],
        "factor_confidence": [str(value) for value in raw["factor_confidence"]],
        "polynomial": None if raw["polynomial"] is None else str(raw["polynomial"]),
        "root": None if raw["root"] is None else str(raw["root"]),
        "method": str(raw["method"]),
        "degree": None if raw["degree"] is None else int(raw["degree"]),
        "relations": int(raw["relations"]),
        "rank": int(raw["rank"]),
        "elapsed_ms": int(raw["elapsed_ms"]),
        "timings": [
            {
                "stage": str(entry["stage"]),
                "degree": None if entry["degree"] is None else int(entry["degree"]),
                "elapsed_ms": int(entry["elapsed_ms"]),
                "detail": str(entry["detail"]),
            }
            for entry in raw["timings"]
        ],
    }
