import inspect
import json

import pytest

import graynfs


def test_proof_method_routing_matrix_small_inputs() -> None:
    """aks proves while bpsw, frobenius, and bpsw_frobenius stay probable-only."""
    for method in ("aks", "bpsw", "frobenius", "bpsw_frobenius"):
        assert graynfs.is_prime_proof(101, method=method) is True
        assert graynfs.is_prime_proof(91, method=method) is False
        assert graynfs.is_prime_proof(561, method=method) is False
    assert graynfs.is_probable_prime(101) is True
    assert graynfs.is_probable_prime(91) is False


def test_proof_in_cap_boundary_composite_is_composite() -> None:
    """2**64 - 1 sits inside the cap, so aks False agrees with the probable gate."""
    cap_edge = (1 << 64) - 1
    assert cap_edge.bit_length() == 64
    assert graynfs.is_prime_proof(cap_edge, method="aks") is False
    assert graynfs.is_probable_prime(cap_edge) is False


def test_proof_in_cap_boundary_prime_stays_probable() -> None:
    """Largest 64-bit prime is in-cap and probable; aks itself is not run here."""
    largest_64bit_prime = 18446744073709551557
    assert largest_64bit_prime.bit_length() == 64
    assert graynfs.is_probable_prime(largest_64bit_prime) is True


def test_proof_over_cap_composite_is_composite() -> None:
    """2**64 is above the cap and composite under every method."""
    over_cap_composite = 1 << 64
    assert over_cap_composite.bit_length() == 65
    assert graynfs.is_prime_proof(over_cap_composite, method="aks") is False
    assert graynfs.is_probable_prime(over_cap_composite) is False
    for method in ("bpsw", "frobenius", "bpsw_frobenius"):
        assert graynfs.is_prime_proof(over_cap_composite, method=method) is False


def test_proof_over_cap_prime_is_unproven_not_composite() -> None:
    """M89 is prime but above the cap, so aks False pairs with probable True."""
    mersenne89 = (1 << 89) - 1
    assert mersenne89.bit_length() == 89
    assert graynfs.is_prime_proof(mersenne89, method="aks") is False
    assert graynfs.is_probable_prime(mersenne89) is True
    for method in ("bpsw", "frobenius", "bpsw_frobenius"):
        assert graynfs.is_prime_proof(mersenne89, method=method) is True


def test_proof_unknown_method_names_valid_set_and_unproven_rule() -> None:
    """Unknown methods fail as input errors that restate the routing contract."""
    assert list(inspect.signature(graynfs.is_prime_proof).parameters) == ["n", "method"]
    assert inspect.signature(graynfs.is_prime_proof).parameters["method"].default == "aks"
    with pytest.raises(graynfs.GrayNFSError) as excinfo:
        graynfs.is_prime_proof(101, method="miller-rabin")
    assert excinfo.value.category == "input"
    assert "aks|bpsw|frobenius|bpsw_frobenius" in excinfo.value.detail
    assert "unproven, not composite" in excinfo.value.detail
    with pytest.raises(TypeError):
        graynfs.is_prime_proof(101, method=64)  # type: ignore[arg-type]


def test_proof_report_keys_stay_additive_only() -> None:
    """factor_report keys may grow but must keep every documented entry."""
    report = graynfs.factor_report(15)
    assert set(
        [
            "factors",
            "factor_confidence",
            "polynomial",
            "root",
            "method",
            "degree",
            "relations",
            "rank",
            "elapsed_ms",
            "timings",
        ]
    ) <= set(report.keys())
    assert report["factors"] == [3, 5]
    assert len(report["factor_confidence"]) == len(report["factors"])
    assert isinstance(report["timings"], list)
    assert report["timings"]
    for entry in report["timings"]:
        assert set(["stage", "degree", "elapsed_ms", "detail"]) <= set(entry.keys())
    parsed = json.loads(json.dumps(report))
    assert set(
        [
            "factors",
            "factor_confidence",
            "polynomial",
            "root",
            "method",
            "degree",
            "relations",
            "rank",
            "elapsed_ms",
            "timings",
        ]
    ) <= set(parsed.keys())
    for entry in parsed["timings"]:
        assert set(["stage", "degree", "elapsed_ms", "detail"]) <= set(entry.keys())
