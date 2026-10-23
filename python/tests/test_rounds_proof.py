import inspect

import pytest

import graynfs


def test_existing_signature_unchanged() -> None:
    assert list(inspect.signature(graynfs.is_probable_prime).parameters) == ["n"]
    with pytest.raises(TypeError):
        graynfs.is_probable_prime(97, rounds=64)  # type: ignore[call-arg]


def test_fast_rounds_tunable() -> None:
    assert graynfs.is_probable_prime_fast(101, rounds=1)
    assert graynfs.is_probable_prime_fast(101, rounds=1000)
    assert not graynfs.is_probable_prime_fast(91, rounds=1000)
    assert not graynfs.is_probable_prime_fast(5459, rounds=1000)


def test_fast_default_rounds() -> None:
    assert graynfs.is_probable_prime_fast(101)
    assert not graynfs.is_probable_prime_fast(91)
    assert not graynfs.is_probable_prime_fast(-5)


def test_fast_rejects_bad_rounds() -> None:
    with pytest.raises(TypeError):
        graynfs.is_probable_prime_fast(101, rounds="64")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        graynfs.is_probable_prime_fast(101, rounds=-1)
    with pytest.raises(TypeError):
        graynfs.is_probable_prime_fast("101", rounds=64)  # type: ignore[arg-type]


def test_proof_accepts_prime_all_methods() -> None:
    assert graynfs.is_prime_proof(101) == True
    for method in ("bpsw", "frobenius", "aks", "bpsw_frobenius"):
        assert graynfs.is_prime_proof(101, method=method) is True
    assert graynfs.is_prime_proof(-5) is False


def test_proof_rejects_composites_all_methods() -> None:
    for method in ("bpsw", "frobenius", "aks", "bpsw_frobenius"):
        assert graynfs.is_prime_proof(91, method=method) is False
        assert graynfs.is_prime_proof(561, method=method) is False
        assert graynfs.is_prime_proof(5459, method=method) is False


def test_proof_rejects_unknown_method() -> None:
    with pytest.raises(ValueError):
        graynfs.is_prime_proof(101, method="miller-rabin")
    with pytest.raises(TypeError):
        graynfs.is_prime_proof(101, method=64)  # type: ignore[arg-type]


def test_proof_has_no_rounds() -> None:
    assert list(inspect.signature(graynfs.is_prime_proof).parameters) == ["n", "method"]
    with pytest.raises(TypeError):
        graynfs.is_prime_proof(101, rounds=64)  # type: ignore[call-arg]


def test_proof_maps_input_error_category() -> None:
    try:
        graynfs.is_prime_proof(101, method="miller-rabin")
    except graynfs.GrayNFSError as exc:
        assert exc.category == "input"
    else:
        raise AssertionError("expected GrayNFSError")
