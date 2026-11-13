use graynfs_core::prime::{AKS_PROOF_MAX_BITS, aks_calls, frobenius_calls, is_prime_u64};
use graynfs_core::{Integer, is_prime_ecpp_proof, prove_ecpp, verify_ecpp_certificate};
use num_bigint::BigInt;
use num_traits::Num;

static ECPP_PROOF_SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn proof_guard() -> std::sync::MutexGuard<'static, ()> {
    ECPP_PROOF_SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

fn decimal(text: &str) -> Integer {
    BigInt::from_str_radix(text, 10).unwrap()
}

#[test]
fn aks_cap_stays_at_64_bits() {
    assert_eq!(AKS_PROOF_MAX_BITS, 64);
    assert!(is_prime_u64(2));
    let wide = Integer::from(1) << 70;
    assert!(!graynfs_core::prime::is_prime_aks_proof(&wide));
}

#[test]
fn ecpp_proof_delegates_u64_range() {
    let _guard = proof_guard();
    assert!(is_prime_ecpp_proof(&decimal("101")));
    assert!(is_prime_ecpp_proof(&decimal("9223372036854775837")));
    assert!(!is_prime_ecpp_proof(&decimal("561")));
    assert!(!is_prime_ecpp_proof(&decimal("2047")));
    assert!(!is_prime_ecpp_proof(&decimal("4")));
}

#[test]
fn ecpp_proof_certifies_mid_range_prime() {
    let _guard = proof_guard();
    let prime = decimal("618970019642690137449562111");
    assert_eq!(prime.bits(), 89);
    let snap_aks = aks_calls();
    let snap_frobenius = frobenius_calls();
    assert!(is_prime_ecpp_proof(&prime));
    let cert = prove_ecpp(&prime).expect("mid-range prime must certify");
    assert!(verify_ecpp_certificate(&cert));
    assert_eq!(aks_calls(), snap_aks);
    assert_eq!(frobenius_calls(), snap_frobenius);
}

#[test]
fn ecpp_proof_declines_composites_and_out_of_scope() {
    let _guard = proof_guard();
    assert!(!is_prime_ecpp_proof(&decimal(
        "79228162422030617973593122087"
    )));
    assert!(!is_prime_ecpp_proof(&decimal(
        "170141183460469231731687303715884105728"
    )));
    assert!(prove_ecpp(&decimal("101")).is_none());
    assert!(prove_ecpp(&(Integer::from(1) << 129)).is_none());
    assert!(!is_prime_ecpp_proof(&(Integer::from(1) << 129)));
}
