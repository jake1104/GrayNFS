use graynfs_core::{
    CancellationToken, FactorError, FactorOptions, factor, factor_with_cancellation, parse_decimal,
};

#[test]
fn token_starts_uncancelled_and_cancel_sets_flag() {
    let token = CancellationToken::new();
    assert!(!token.is_cancelled());
    token.cancel();
    assert!(token.is_cancelled());
}

#[test]
fn token_clone_shares_cancellation_state() {
    let token = CancellationToken::new();
    let shared = token.clone();
    assert!(!shared.is_cancelled());
    token.cancel();
    assert!(shared.is_cancelled());
}

#[test]
fn pre_cancelled_token_factor_ends_cancelled() {
    let token = CancellationToken::new();
    token.cancel();
    let number = parse_decimal("100160063").expect("valid input");
    let outcome = factor_with_cancellation(&number, &FactorOptions::default(), &token, |_| Ok(()));
    assert!(
        matches!(outcome, Err(FactorError::Cancelled)),
        "pre-cancelled token must end as Cancelled, got {outcome:?}"
    );
}

#[test]
fn progress_driven_cancel_ends_cancelled() {
    let token = CancellationToken::new();
    let number = parse_decimal("100160063").expect("valid input");
    let outcome = factor_with_cancellation(&number, &FactorOptions::default(), &token, |_| {
        token.cancel();
        Ok(())
    });
    assert!(
        matches!(outcome, Err(FactorError::Cancelled)),
        "mid-run cancel must end as Cancelled, got {outcome:?}"
    );
}

#[test]
fn live_token_factor_matches_plain_factor() {
    let token = CancellationToken::new();
    let number = parse_decimal("100160063").expect("valid input");
    let options = FactorOptions::default();
    let expected = factor(&number, &options).expect("plain factor succeeds");
    let actual = factor_with_cancellation(&number, &options, &token, |_| Ok(()))
        .expect("token factor succeeds");
    assert_eq!(actual.factors, expected.factors);
    assert_eq!(actual.method, expected.method);
}
