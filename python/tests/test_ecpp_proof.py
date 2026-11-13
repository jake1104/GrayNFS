import graynfs


def test_ecpp_proves_small_inputs_like_u64_gate() -> None:
    assert graynfs.is_prime_proof(101, method="ecpp") is True
    assert graynfs.is_prime_proof(2, method="ecpp") is True
    assert graynfs.is_prime_proof(91, method="ecpp") is False
    assert graynfs.is_prime_proof(561, method="ecpp") is False


def test_ecpp_proves_over_cap_prime_where_aks_stays_unproven() -> None:
    mersenne89 = (1 << 89) - 1
    assert mersenne89.bit_length() == 89
    assert graynfs.is_prime_proof(mersenne89, method="aks") is False
    assert graynfs.is_prime_proof(mersenne89, method="ecpp") is True
    assert graynfs.is_probable_prime(mersenne89) is True


def test_ecpp_declines_over_cap_composites() -> None:
    assert graynfs.is_prime_proof(1 << 64, method="ecpp") is False
    assert graynfs.is_prime_proof((1 << 89) - 1 + 2, method="ecpp") is False
