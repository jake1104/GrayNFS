from typing import Any

import pytest

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
    assert graynfs.version() == "0.2.0"


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


def test_factor_report_fields() -> None:
    report = graynfs.factor_report(100160063)
    assert report["factors"] == [10007, 10009]
    assert isinstance(report["method"], str)
    assert report["degree"] is None or isinstance(report["degree"], int)
    assert isinstance(report["relations"], int)
    assert isinstance(report["rank"], int)
    assert isinstance(report["elapsed_ms"], int)
    assert report["timings"]
    for row in report["timings"]:
        assert set(row) == {"stage", "degree", "elapsed_ms", "detail"}
        assert isinstance(row["stage"], str)
        assert row["degree"] is None or isinstance(row["degree"], int)
        assert isinstance(row["elapsed_ms"], int)
        assert isinstance(row["detail"], str)


def test_factor_report_matches_factor_single_call() -> None:
    options: Any = {
        "degree": 3,
        "algebraic_base_bound": 500,
        "sieve_a_bound": 1600,
        "sieve_b_bound": 700,
        "timeout_seconds": 120,
    }
    assert graynfs.factor_report(100160063, **options)["factors"] == graynfs.factor(
        100160063, **options
    )


def test_factor_report_timeout_surfaces_value_error() -> None:
    try:
        graynfs.factor_report(
            100160063,
            degree=3,
            timeout_seconds=0,
            use_fermat=False,
            use_pm1=False,
            use_pp1=False,
            use_ecm=False,
        )
    except ValueError as exc:
        assert str(exc).startswith("resource-limit: ")
        assert "deadline" in str(exc)
    else:
        raise AssertionError("expected ValueError")


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


def test_gray_error_value_error_compatibility() -> None:
    assert issubclass(graynfs.GrayNFSError, ValueError)
    with pytest.raises(ValueError):
        graynfs.factor(0)
    with pytest.raises(graynfs.GrayNFSError):
        graynfs.factor(0)


def test_gray_error_categories() -> None:
    try:
        graynfs.factor(0)
    except graynfs.GrayNFSError as exc:
        assert exc.category == "input"
        assert exc.detail == "n must be a positive integer"
    else:
        raise AssertionError("expected GrayNFSError")
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
    except graynfs.GrayNFSError as exc:
        assert exc.category == "resource-limit"
        assert "deadline" in exc.detail
    else:
        raise AssertionError("expected GrayNFSError")


def test_gray_error_message_identity() -> None:
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
    except graynfs.GrayNFSError as wrapped:
        wrapped_message = str(wrapped)
        assert wrapped.category == "resource-limit"
        assert wrapped_message == wrapped.category + ": " + wrapped.detail
    else:
        raise AssertionError("expected GrayNFSError")
    try:
        graynfs._graynfs.factor_decimal(
            "100160063",
            3,
            2000,
            2000,
            700,
            120,
            0,
            None,
            True,
            False,
            True,
            False,
            False,
            False,
            True,
            None,
            False,
        )
    except ValueError as raw:
        assert type(raw) is ValueError
        assert str(raw) == wrapped_message
    else:
        raise AssertionError("expected ValueError")


def test_gray_error_callback_failure_maps_internal() -> None:
    marker_message = "boom-from-callback"

    def failing_callback(event: dict) -> None:
        raise ValueError(marker_message)

    try:
        graynfs.factor_with_progress(15, failing_callback)
    except graynfs.GrayNFSError as exc:
        assert exc.category == "internal"
        assert marker_message in exc.detail
        assert marker_message in str(exc)
        assert str(exc) == exc.category + ": " + exc.detail
    else:
        raise AssertionError("expected GrayNFSError")
    try:
        graynfs._graynfs.factor_with_progress_decimal(
            "15",
            failing_callback,
            None,
            2000,
            2000,
            700,
            120,
            None,
            None,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            None,
            False,
        )
    except ValueError as raw:
        assert type(raw) is ValueError
        assert marker_message in str(raw)
    else:
        raise AssertionError("expected ValueError")


def test_to_gray_error_maps_known_prefixes() -> None:
    cases = (
        ("input: bad digit", "input", "bad digit"),
        ("resource-limit: deadline exceeded", "resource-limit", "deadline exceeded"),
        ("cancelled: stopped", "cancelled", "stopped"),
        ("internal: boom", "internal", "boom"),
        ("input: invalid decimal integer: abc", "input", "invalid decimal integer: abc"),
    )
    for message, category, detail in cases:
        mapped = graynfs._to_gray_error(ValueError(message))
        assert type(mapped) is graynfs.GrayNFSError
        assert mapped.category == category
        assert mapped.detail == detail
        assert str(mapped) == message


def test_to_gray_error_passes_unknown_through() -> None:
    original = ValueError("boom-from-callback")
    assert graynfs._to_gray_error(original) is original
    empty = ValueError("")
    assert graynfs._to_gray_error(empty) is empty
