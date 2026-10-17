use crate::error::FactorError;
use crate::integer::{Integer, integer_sqrt};
use crate::options::FactorOptions;

#[derive(Clone, Debug)]
pub struct PolynomialPair {
    pub degree: usize,
    pub coefficients: Vec<Integer>,
    pub m: Integer,
}

impl PolynomialPair {
    pub fn constant(&self) -> Integer {
        self.coefficients
            .first()
            .cloned()
            .unwrap_or_else(|| Integer::from(0))
    }

    pub fn homogeneous(&self, a: &Integer, b: &Integer) -> Integer {
        let degree = self.coefficients.len() - 1;
        let mut a_powers = vec![Integer::from(1); degree + 1];
        for index in 1..=degree {
            a_powers[index] = &a_powers[index - 1] * a;
        }
        let mut b_powers = vec![Integer::from(1); degree + 1];
        for index in 1..=degree {
            b_powers[index] = &b_powers[index - 1] * b;
        }
        let mut result = Integer::from(0);
        for (index, coefficient) in self.coefficients.iter().enumerate() {
            result += coefficient * &a_powers[index] * &b_powers[degree - index];
        }
        result
    }

    pub fn evaluate(&self, value: &Integer) -> Integer {
        let mut result = Integer::from(0);
        let mut power = Integer::from(1);
        for coefficient in &self.coefficients {
            result += coefficient * &power;
            power *= value;
        }
        result
    }

    pub fn f_at_m(&self) -> Integer {
        self.evaluate(&self.m)
    }

    pub fn describe(&self) -> String {
        format!(
            "x^{} + {} with root {} modulo N",
            self.degree,
            self.constant(),
            self.m
        )
    }
}

pub fn select_polynomial(
    n: &Integer,
    options: &FactorOptions,
) -> Result<PolynomialPair, FactorError> {
    if options.degree != Some(2) {
        return Err(FactorError::UnsupportedDegree(options.degree.unwrap_or(2)));
    }
    if n <= &Integer::from(1) {
        return Err(FactorError::InvalidInput(
            "N must be greater than one".to_string(),
        ));
    }
    let mut m = integer_sqrt(n);
    if &m * &m == *n {
        m -= Integer::from(1);
    }
    let mut attempts = 0_usize;
    let constant = loop {
        if m <= Integer::from(0) {
            return Err(FactorError::NoPolynomial);
        }
        let candidate = n - &m * &m;
        if candidate > Integer::from(0) && !is_square(&candidate) {
            break candidate;
        }
        m -= Integer::from(1);
        attempts += 1;
        if attempts > 32 {
            return Err(FactorError::NoPolynomial);
        }
    };
    Ok(PolynomialPair {
        degree: 2,
        coefficients: vec![constant, Integer::from(0), Integer::from(1)],
        m,
    })
}

fn is_square(value: &Integer) -> bool {
    if value <= &Integer::from(0) {
        return false;
    }
    let root = integer_sqrt(value);
    &root * &root == *value
}
