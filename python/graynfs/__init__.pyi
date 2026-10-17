from typing import Any, Callable

ProgressCallback = Callable[[dict[str, Any]], None]

def factor(
    n: int,
    *,
    degree: int | None = ...,
    rational_base_bound: int = ...,
    algebraic_base_bound: int = ...,
    sieve_a_bound: int = ...,
    sieve_b_bound: int = ...,
    timeout_seconds: int | None = ...,
    allow_fallback: bool = ...,
    use_fermat: bool = ...,
    use_squfof: bool = ...,
    use_pm1: bool = ...,
    use_pp1: bool = ...,
    use_ecm: bool = ...,
    use_mpqs: bool = ...,
    threads: int | None = ...,
    perf_mode: bool = ...,
) -> list[int]: ...
def factor_with_progress(
    n: int,
    callback: ProgressCallback,
    *,
    degree: int | None = ...,
    rational_base_bound: int = ...,
    algebraic_base_bound: int = ...,
    sieve_a_bound: int = ...,
    sieve_b_bound: int = ...,
    timeout_seconds: int | None = ...,
    allow_fallback: bool = ...,
    use_fermat: bool = ...,
    use_squfof: bool = ...,
    use_pm1: bool = ...,
    use_pp1: bool = ...,
    use_ecm: bool = ...,
    use_mpqs: bool = ...,
    threads: int | None = ...,
    perf_mode: bool = ...,
) -> list[int]: ...
def stage_timings(
    n: int,
    *,
    degree: int | None = ...,
    rational_base_bound: int = ...,
    algebraic_base_bound: int = ...,
    sieve_a_bound: int = ...,
    sieve_b_bound: int = ...,
    timeout_seconds: int | None = ...,
    allow_fallback: bool = ...,
    use_fermat: bool = ...,
    use_squfof: bool = ...,
    use_pm1: bool = ...,
    use_pp1: bool = ...,
    use_ecm: bool = ...,
    use_mpqs: bool = ...,
    threads: int | None = ...,
    perf_mode: bool = ...,
) -> list[dict[str, Any]]: ...
def is_probable_prime(n: int) -> bool: ...
def version() -> str: ...
