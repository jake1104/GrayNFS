use graynfs_core::prime::{
    aks_calls, frobenius_calls, is_baillie_psw_prime, is_frobenius_probable_prime,
    is_prime_aks_proof, is_prime_u64, is_probable_prime,
};
use graynfs_core::{
    Integer, certificate_to_json, prove_ecpp, verify_ecpp_certificate, verify_ecpp_json,
};
use num_bigint::BigInt;
use num_traits::Num;

const CORPUS_PRIMES: [(&str, u64); 12] = [
    ("9223372036854775837", 64),
    ("2361183241434822606859", 72),
    ("604462909807314587353111", 80),
    ("154742504910672534362390567", 88),
    ("39614081257132168796771975177", 96),
    ("10141204801825835211973625643089", 104),
    ("2596148429267413814265248164610099", 112),
    ("664613997892457936451903530140172297", 120),
    ("170141183460469231731687303715884105757", 128),
    ("618970019642690137449562111", 89),
    ("162259276829213363391578010288127", 107),
    ("170141183460469231731687303715884105727", 127),
];

const CORPUS_COMPOSITES: [&str; 9] = [
    "340282366920938463463374607431768211455",
    "340282366920938461286658806734041124249",
    "79228162422030616971593122087",
    "79228162237563176810023223171",
    "1000112004278059472142857",
    "170141183460469231731687303715884105728",
    "1000000000000000000000000000058",
    "2047",
    "3215031751",
];

fn decimal(text: &str) -> Integer {
    BigInt::from_str_radix(text, 10).unwrap()
}

static ECPP_SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn ecpp_guard() -> std::sync::MutexGuard<'static, ()> {
    ECPP_SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

#[test]
fn aks_authority_on_small_inputs() {
    let _guard = ecpp_guard();
    for small in [101_u64, 127, 1009, 7919] {
        assert!(is_prime_aks_proof(&decimal(&small.to_string())));
        assert!(is_prime_u64(small));
    }
    for small in [91_u64, 561, 2047, 3277] {
        assert!(!is_prime_aks_proof(&decimal(&small.to_string())));
        assert!(!is_probable_prime(&decimal(&small.to_string())));
    }
}

#[test]
fn cross_verification_matches_bpsw_gate() {
    let _guard = ecpp_guard();
    eprintln!("bits | value-prefix | bpsw | aks | ecpp | ecpp_ms");
    for (text, bits) in CORPUS_PRIMES {
        let prime = decimal(text);
        assert_eq!(prime.bits(), bits, "corpus bit length drift: {text}");
        let bpsw = is_probable_prime(&prime)
            && is_baillie_psw_prime(&prime)
            && is_frobenius_probable_prime(&prime);
        assert!(bpsw, "BPSW must accept corpus prime {text}");
        let start = std::time::Instant::now();
        let snap_aks = aks_calls();
        let snap_frobenius = frobenius_calls();
        let cert = prove_ecpp(&prime).expect("ECPP must certify corpus prime");
        let elapsed = start.elapsed().as_millis();
        assert!(verify_ecpp_certificate(&cert));
        let rendered = certificate_to_json(&cert);
        assert!(verify_ecpp_json(&prime, &rendered).expect("schema must hold"));
        assert_eq!(aks_calls(), snap_aks);
        assert_eq!(frobenius_calls(), snap_frobenius);
        let aks_text = if bits > 64 {
            assert!(
                !is_prime_aks_proof(&prime),
                "over-cap aks must stay False at {text}"
            );
            "false"
        } else {
            "skip"
        };
        let prefix: String = text.chars().take(12).collect();
        eprintln!("{bits} | {prefix} | true | {aks_text} | true | {elapsed}");
    }
    for text in CORPUS_COMPOSITES {
        let composite = decimal(text);
        assert!(
            !is_probable_prime(&composite),
            "BPSW must reject composite {text}"
        );
        assert!(
            !is_baillie_psw_prime(&composite),
            "BPSW stack must reject {text}"
        );
        assert!(
            !is_frobenius_probable_prime(&composite),
            "Frobenius must reject {text}"
        );
        if composite.bits() > 64 {
            assert!(!is_prime_aks_proof(&composite));
        }
        let snap_aks = aks_calls();
        let snap_frobenius = frobenius_calls();
        assert!(
            prove_ecpp(&composite).is_none(),
            "ECPP must decline composite {text}"
        );
        assert_eq!(aks_calls(), snap_aks);
        assert_eq!(frobenius_calls(), snap_frobenius);
        let prefix: String = text.chars().take(12).collect();
        eprintln!(
            "{} | {prefix} | false | false | false | -",
            composite.bits()
        );
    }
}
