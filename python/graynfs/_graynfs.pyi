from typing import Any, Callable

def factor_decimal(
    n: str,
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
) -> list[str]: ...
def factor_with_progress_decimal(
    n: str,
    callback: Callable[[dict[str, Any]], None] | None = ...,
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
) -> list[str]: ...
def stage_timings_decimal(
    n: str,
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
def is_probable_prime_decimal(n: str) -> bool: ...
def is_probable_prime_fast_decimal(n: str, rounds: int = ...) -> bool: ...
def is_prime_proof_decimal(n: str, method: str = ...) -> bool: ...
def version() -> str: ...
