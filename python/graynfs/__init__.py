from __future__ import annotations

from typing import Any, Callable

from . import _graynfs

__all__ = [
    "GrayNFSError",
    "factor",
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
) -> list[int]:
    try:
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
) -> list[int]:
    if not callable(callback):
        raise TypeError("callback must be callable")
    try:
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
            )
        ]
    except ValueError as exc:
        raise _to_gray_error(exc)


def is_probable_prime(n: int) -> bool:
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 0:
        return False
    return bool(_graynfs.is_probable_prime_decimal(str(n)))


def is_probable_prime_fast(n: int, rounds: int = 64) -> bool:
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


def is_prime_proof(n: int, method: str = "bpsw") -> bool:
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
) -> list[dict[str, Any]]:
    try:
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
            )
        ]
    except ValueError as exc:
        raise _to_gray_error(exc)


def version() -> str:
    return str(_graynfs.version())


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
) -> dict[str, Any]:
    try:
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
