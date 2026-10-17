use thiserror::Error;

#[derive(Debug, Error)]
pub enum FactorError {
    #[error("invalid input: {0}")]
    InvalidInput(String),
    #[error("no suitable polynomial was found")]
    NoPolynomial,
    #[error("the configured search space was exhausted")]
    SearchExhausted,
    #[error("insufficient relations: found {found}, required at least {required}")]
    InsufficientRelations { found: usize, required: usize },
    #[error("the algebraic square root could not be constructed")]
    AlgebraicSquareRootFailed,
    #[error("linear algebra failed: {0}")]
    LinearAlgebra(String),
    #[error("the operation was cancelled")]
    Cancelled,
    #[error("resource limit exceeded: {0}")]
    ResourceLimit(String),
    #[error("integer factorization failed: {0}")]
    Factorization(String),
    #[error("unsupported polynomial degree {0}; the current implementation requires degree 2")]
    UnsupportedDegree(usize),
}
