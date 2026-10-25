from __future__ import annotations

from typing import Any, Callable

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

ProgressCallback = Callable[[dict[str, Any]], None]


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
    use_squfof: bool = True,
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
) -> list[int]:
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
                use_squfof,
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
    use_squfof: bool = True,
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
) -> list[int]:
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
                use_squfof,
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
    use_squfof: bool = True,
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
) -> list[dict[str, Any]]:
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
                use_squfof,
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
            )
        ]
    except ValueError as exc:
        raise _to_gray_error(exc)


def version() -> str:
    return str(_graynfs.version())


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
    use_squfof: bool = True,
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
                    use_squfof=use_squfof,
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
                )
            )
        except GrayNFSError as exc:
            if fail_fast:
                raise
            results.append(exc)
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
    use_squfof: bool = True,
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
) -> dict[str, Any]:
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
            use_squfof,
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
        )
    except ValueError as exc:
        raise _to_gray_error(exc)
    return {
        "factors": [int(value) for value in raw["factors"]],
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
