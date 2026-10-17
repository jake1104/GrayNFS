use crate::poly::ModPolynomial;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExtensionField {
    prime: u64,
    modulus: ModPolynomial,
    size: u64,
}

impl ExtensionField {
    pub fn new(prime: u64, modulus: ModPolynomial) -> Option<Self> {
        let degree = modulus.degree();
        if prime < 2 || degree < 1 || modulus.is_constant() {
            return None;
        }
        if modulus.coefficients.iter().any(|value| *value >= prime) {
            return None;
        }
        if !crate::poly::is_irreducible(&modulus, prime) {
            return None;
        }
        let size = prime.checked_pow(degree as u32)?;
        if size < 2 {
            return None;
        }
        Some(Self {
            prime,
            modulus,
            size,
        })
    }

    pub fn prime(&self) -> u64 {
        self.prime
    }

    pub fn size(&self) -> u64 {
        self.size
    }

    pub fn degree(&self) -> usize {
        self.modulus.degree()
    }

    pub fn modulus(&self) -> &ModPolynomial {
        &self.modulus
    }

    pub fn zero(&self) -> Vec<u64> {
        vec![0_u64; self.degree()]
    }

    pub fn one(&self) -> Vec<u64> {
        let mut value = self.zero();
        if let Some(first) = value.first_mut() {
            *first = 1;
        }
        value
    }

    pub fn is_zero(&self, value: &[u64]) -> bool {
        value
            .iter()
            .all(|coefficient| *coefficient % self.prime == 0)
    }

    pub fn equal(&self, left: &[u64], right: &[u64]) -> bool {
        self.normalize(left) == self.normalize(right)
    }

    pub fn from_scalar(&self, value: u64) -> Vec<u64> {
        let mut result = self.zero();
        if let Some(first) = result.first_mut() {
            *first = value % self.prime;
        }
        result
    }

    pub fn normalize(&self, value: &[u64]) -> Vec<u64> {
        let degree = self.degree();
        let mut result = vec![0_u64; degree];
        for (index, coefficient) in value.iter().take(degree).enumerate() {
            result[index] = coefficient % self.prime;
        }
        result
    }

    pub fn mul(&self, left: &[u64], right: &[u64]) -> Vec<u64> {
        let first = ModPolynomial::new(self.normalize(left), self.prime);
        let second = ModPolynomial::new(self.normalize(right), self.prime);
        let product = first.mul(&second);
        self.normalize(&product.rem(&self.modulus).coefficients)
    }

    pub fn pow(&self, value: &[u64], mut exponent: u64) -> Vec<u64> {
        let mut result = self.one();
        let mut base = self.normalize(value);
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = self.mul(&result, &base);
            }
            base = self.mul(&base, &base);
            exponent >>= 1;
        }
        result
    }

    pub fn inverse(&self, value: &[u64]) -> Option<Vec<u64>> {
        if self.is_zero(value) {
            return None;
        }
        Some(self.pow(value, self.size.checked_sub(2)?))
    }

    pub fn frobenius(&self, value: &[u64]) -> Vec<u64> {
        self.pow(value, self.prime)
    }
}

pub fn is_square(field: &ExtensionField, value: &[u64]) -> bool {
    if field.is_zero(value) {
        return true;
    }
    let exponent = (field.size() - 1) / 2;
    field.equal(&field.pow(value, exponent), &field.one())
}

pub fn non_residue(field: &ExtensionField) -> Option<Vec<u64>> {
    let mut candidate = 2_u64;
    while candidate < field.prime() {
        let value = field.from_scalar(candidate);
        if !is_square(field, &value) {
            return Some(field.normalize(&value));
        }
        candidate += 1;
    }
    if field.degree() >= 2 {
        let mut lead = 0_u64;
        while lead < field.prime() {
            let mut value = field.zero();
            value[0] = lead % field.prime();
            value[1] = 1;
            if !is_square(field, &value) {
                return Some(field.normalize(&value));
            }
            lead += 1;
        }
    }
    None
}

pub fn square_root(field: &ExtensionField, value: &[u64]) -> Option<Vec<u64>> {
    let target = field.normalize(value);
    if field.is_zero(&target) {
        return Some(field.zero());
    }
    let mut twos = 0_u32;
    let mut odd_part = field.size() - 1;
    while odd_part.is_multiple_of(2) {
        odd_part /= 2;
        twos += 1;
    }
    if !is_square(field, &target) {
        return None;
    }
    if twos == 1 {
        return Some(field.pow(&target, (field.size() + 1) / 4));
    }
    let residue = non_residue(field)?;
    let mut current = field.pow(&residue, odd_part);
    let mut reduced = field.pow(&target, odd_part);
    let mut root = field.pow(&target, odd_part.div_ceil(2));
    let mut width = twos;
    loop {
        if field.equal(&reduced, &field.one()) {
            return Some(root);
        }
        let mut order = 1_u32;
        let mut probe = field.mul(&reduced, &reduced);
        while order < width && !field.equal(&probe, &field.one()) {
            probe = field.mul(&probe, &probe);
            order += 1;
        }
        if order >= width {
            return None;
        }
        let mut factor = current.clone();
        for _ in 0..(width - order - 1) {
            factor = field.mul(&factor, &factor);
        }
        reduced = field.mul(&reduced, &field.mul(&factor, &factor));
        root = field.mul(&root, &factor);
        current = field.mul(&factor, &factor);
        width = order;
    }
}

pub fn multiplicative_order(field: &ExtensionField, value: &[u64]) -> u64 {
    let mut order = 1_u64;
    let mut current = field.normalize(value);
    while !field.equal(&current, &field.one()) {
        current = field.mul(&current, value);
        order += 1;
        if order > field.size() {
            break;
        }
    }
    order
}

pub fn verify_square_root(field: &ExtensionField, value: &[u64], root: &[u64]) -> bool {
    let squared = field.mul(root, root);
    field.equal(&squared, &field.normalize(value))
}

#[cfg(test)]
mod tests {
    use super::{ExtensionField, is_square, multiplicative_order, non_residue, square_root};
    use crate::poly::ModPolynomial;

    fn field(prime: u64, degree: usize) -> ExtensionField {
        assert!(degree > 0);
        let total = (prime as usize).pow(degree as u32);
        for index in 0..total {
            let mut coefficients = vec![0_u64; degree + 1];
            let mut remaining = index;
            for coefficient in coefficients.iter_mut().take(degree) {
                *coefficient = (remaining % prime as usize) as u64;
                remaining /= prime as usize;
            }
            coefficients[degree] = 1;
            let candidate = ModPolynomial::new(coefficients, prime);
            if crate::poly::is_irreducible(&candidate, prime) {
                return ExtensionField::new(prime, candidate).expect("field");
            }
        }
        panic!("no irreducible of degree {degree} over {prime}");
    }

    fn assert_field(prime: u64, degree: usize) {
        let extension = field(prime, degree);
        for value in 0..prime {
            let target = extension.from_scalar(value);
            match square_root(&extension, &target) {
                Some(root) => {
                    let squared = extension.mul(&root, &root);
                    assert!(
                        extension.equal(&squared, &target),
                        "prime {prime} value {value} root mismatch"
                    );
                }
                None => assert!(
                    !is_square(&extension, &target),
                    "prime {prime} value {value} should be a square"
                ),
            }
        }
        let residue = non_residue(&extension).expect("non residue");
        assert!(!is_square(&extension, &residue), "residue must be square");
        let degree = extension.degree();
        let mut found = false;
        'outer: for lead in 1..prime {
            for second in 0..prime {
                if degree < 2 && second != 0 {
                    continue;
                }
                let mut value = extension.zero();
                value[0] = lead;
                if degree > 1 {
                    value[1] = second;
                }
                if is_square(&extension, &value) && multiplicative_order(&extension, &value) > 1 {
                    found = true;
                    break 'outer;
                }
            }
        }
        assert!(found, "no generator of order above one");
    }

    #[test]
    fn square_roots_over_prime_fields() {
        assert_field(5, 1);
        assert_field(7, 1);
        assert_field(13, 1);
    }

    #[test]
    fn square_roots_over_extension_fields() {
        assert_field(5, 2);
        assert_field(3, 3);
        assert_field(7, 4);
        assert_field(3, 6);
    }

    #[test]
    fn field_rejects_invalid_shape() {
        assert!(ExtensionField::new(5, ModPolynomial::new(vec![1], 5)).is_none());
        assert!(ExtensionField::new(1, ModPolynomial::new(vec![0, 1], 1)).is_none());
    }

    #[test]
    fn scalar_field_matches_modular_arithmetic() {
        let extension = field(11, 1);
        for value in 0..11_u64 {
            let target = extension.from_scalar(value);
            if let Some(root) = square_root(&extension, &target) {
                assert_eq!(extension.mul(&root, &root), extension.normalize(&target));
            }
        }
    }

    #[test]
    fn inverse_is_consistent() {
        let extension = field(7, 2);
        for value in 1..7_u64 {
            let target = extension.from_scalar(value);
            let inverse = extension.inverse(&target).expect("inverse");
            assert_eq!(extension.mul(&target, &inverse), extension.one());
        }
    }

    #[test]
    fn frobenius_is_identity_on_prime_field() {
        let extension = field(11, 1);
        for value in 0..11_u64 {
            let target = extension.from_scalar(value);
            assert_eq!(extension.frobenius(&target), extension.normalize(&target));
        }
    }
}
