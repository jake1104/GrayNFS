from typing import Any

import graynfs


def test_factor_small_composite() -> None:
    assert graynfs.factor(15) == [3, 5]


def test_probable_prime() -> None:
    assert graynfs.is_probable_prime(97)
    assert not graynfs.is_probable_prime(91)


def test_invalid_input() -> None:
    try:
        graynfs.factor(0)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_factor_auto_degree() -> None:
    assert graynfs.factor(100160063) == [10007, 10009]


def test_factor_explicit_degree() -> None:
    factors = graynfs.factor(
        100160063,
        degree=3,
        algebraic_base_bound=500,
        sieve_a_bound=1600,
        sieve_b_bound=700,
        timeout_seconds=120,
    )
    assert factors == [10007, 10009]


def test_factor_with_progress() -> None:
    events: list[dict] = []
    factors = graynfs.factor_with_progress(
        100160063,
        events.append,
        degree=3,
        algebraic_base_bound=500,
        sieve_a_bound=1600,
        sieve_b_bound=700,
        timeout_seconds=120,
        use_fermat=False,
        use_pm1=False,
        use_pp1=False,
        use_ecm=False,
    )
    assert factors == [10007, 10009]
    assert events
    assert events[-1]["stage"] == "gnfs"
    assert events[-1]["current"] == events[-1]["total"]


def test_factor_with_progress_rejects_non_callable() -> None:
    callback: Any = "not callable"
    try:
        graynfs.factor_with_progress(15, callback)
    except TypeError:
        pass
    else:
        raise AssertionError("expected TypeError")


def test_stage_timings() -> None:
    rows = graynfs.stage_timings(
        100160063,
        degree=3,
        algebraic_base_bound=500,
        sieve_a_bound=1600,
        sieve_b_bound=700,
        timeout_seconds=120,
        use_fermat=False,
        use_pm1=False,
        use_pp1=False,
        use_ecm=False,
    )
    assert rows
    stages = {row["stage"] for row in rows}
    assert "relation_collection" in stages
    assert "total" in stages
    for row in rows:
        assert set(row) == {"stage", "degree", "elapsed_ms", "detail"}
        assert isinstance(row["stage"], str)
        assert row["degree"] is None or isinstance(row["degree"], int)
        assert isinstance(row["elapsed_ms"], int)
        assert isinstance(row["detail"], str)


def test_stage_timings_empty_for_trial_division() -> None:
    assert graynfs.stage_timings(15) == []


def test_fermat_shortcircuit_and_opt_out() -> None:
    rows = graynfs.stage_timings(100160063, degree=3, timeout_seconds=120)
    stages = {row["stage"] for row in rows}
    assert "fermat" in stages
    assert "relation_collection" not in stages
    rows_nofermat = graynfs.stage_timings(
        100160063, degree=3, timeout_seconds=120, use_fermat=False, use_pm1=False, use_pp1=False, use_ecm=False
    )
    stages_nofermat = {row["stage"] for row in rows_nofermat}
    assert "fermat" not in stages_nofermat
    assert "relation_collection" in stages_nofermat


def test_threads_explicit():
    assert graynfs.factor(100160063, threads=1) == [10007, 10009]
    assert graynfs.factor(100160063, threads=4) == [10007, 10009]


def test_version() -> None:
    assert graynfs.version() == "0.1.3"


def test_factor_rejects_non_integer() -> None:
    for bad in ["15", 15.0, True, None]:
        try:
            graynfs.factor(bad)  # type: ignore[arg-type]
        except TypeError:
            pass
        else:
            raise AssertionError("expected TypeError for %r" % (bad,))


def test_factor_rejects_non_positive() -> None:
    for bad in [0, -1, -100160063]:
        try:
            graynfs.factor(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError for %r" % (bad,))


def test_is_probable_prime_rejects_non_integer() -> None:
    for bad in ["97", 97.0, True, None]:
        try:
            graynfs.is_probable_prime(bad)  # type: ignore[arg-type]
        except TypeError:
            pass
        else:
            raise AssertionError("expected TypeError for %r" % (bad,))


def test_factor_timeout_surfaces_value_error() -> None:
    try:
        graynfs.factor(
            100160063,
            degree=3,
            timeout_seconds=0,
            use_fermat=False,
            use_pm1=False,
            use_pp1=False,
            use_ecm=False,
        )
    except ValueError as exc:
        assert "deadline" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_factor_max_memory_passthrough() -> None:
    assert graynfs.factor(15, max_memory_mb=None) == [3, 5]
    try:
        graynfs.factor(
            100160063,
            degree=3,
            algebraic_base_bound=500,
            sieve_a_bound=1600,
            sieve_b_bound=700,
            timeout_seconds=120,
            max_memory_mb=0,
            use_fermat=False,
            use_pm1=False,
            use_pp1=False,
            use_ecm=False,
        )
    except ValueError as exc:
        assert "resource-limit" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_factor_with_progress_event_shape() -> None:
    events: list[dict] = []
    assert graynfs.factor_with_progress(15, events.append) == [3, 5]
    assert events
    for event in events:
        assert set(event) == {"stage", "current", "total", "message"}
        assert isinstance(event["stage"], str)
        assert isinstance(event["message"], str)


def test_stage_timings_rejects_bad_input() -> None:
    try:
        graynfs.stage_timings(0)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")
    try:
        graynfs.stage_timings("15")  # type: ignore[arg-type]
    except TypeError:
        pass
    else:
        raise AssertionError("expected TypeError")
