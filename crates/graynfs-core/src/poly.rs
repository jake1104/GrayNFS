#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModPolynomial {
    pub coefficients: Vec<u64>,
    pub modulus: u64,
}

impl ModPolynomial {
    pub fn new(coefficients: Vec<u64>, modulus: u64) -> Self {
        let mut result = Self {
            coefficients,
            modulus,
        };
        result.normalize();
        result
    }

    pub fn zero(modulus: u64) -> Self {
        Self {
            coefficients: vec![0],
            modulus,
        }
    }

    pub fn one(modulus: u64) -> Self {
        Self {
            coefficients: vec![1],
            modulus,
        }
    }

    pub fn degree(&self) -> usize {
        if self.coefficients.len() <= 1 {
            0
        } else {
            self.coefficients.len() - 1
        }
    }

    pub fn normalize(&mut self) {
        if self.coefficients.is_empty() {
            self.coefficients.push(0);
        }
        if self.modulus == 0 {
            return;
        }
        while self.coefficients.len() > 1 && self.coefficients.last() == Some(&0) {
            self.coefficients.pop();
        }
        for value in &mut self.coefficients {
            *value %= self.modulus;
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        let length = self.coefficients.len().max(other.coefficients.len());
        let mut coefficients: Vec<u64> = vec![0; length];
        for (index, coefficient) in coefficients.iter_mut().enumerate() {
            let left = self.coefficients.get(index).copied().unwrap_or(0);
            let right = other.coefficients.get(index).copied().unwrap_or(0);
            *coefficient =
                ((u128::from(left) + u128::from(right)) % (u128::from(self.modulus))) as u64;
        }
        Self::new(coefficients, self.modulus)
    }

    pub fn sub(&self, other: &Self) -> Self {
        let length = self.coefficients.len().max(other.coefficients.len());
        let mut coefficients: Vec<u64> = vec![0; length];
        for (index, coefficient) in coefficients.iter_mut().enumerate() {
            let left = u128::from(self.coefficients.get(index).copied().unwrap_or(0));
            let right = u128::from(other.coefficients.get(index).copied().unwrap_or(0));
            *coefficient =
                ((left + u128::from(self.modulus) - right) % u128::from(self.modulus)) as u64;
        }
        Self::new(coefficients, self.modulus)
    }

    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero(self.modulus);
        }
        let mut coefficients = vec![0_u64; self.coefficients.len() + other.coefficients.len() - 1];
        for (left_index, left) in self.coefficients.iter().enumerate() {
            for (right_index, right) in other.coefficients.iter().enumerate() {
                let product = u128::from(*left) * u128::from(*right)
                    + u128::from(coefficients[left_index + right_index]);
                coefficients[left_index + right_index] =
                    (product % u128::from(self.modulus)) as u64;
            }
        }
        Self::new(coefficients, self.modulus)
    }

    pub fn pow(&self, mut exponent: u64) -> Self {
        let mut result = Self::one(self.modulus);
        let mut base = self.clone();
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exponent >>= 1;
        }
        result
    }

    pub fn divrem(&self, divisor: &Self) -> (Self, Self) {
        assert!(!divisor.is_zero(), "polynomial division by zero");
        let mut remainder = self.clone();
        let divisor_degree = divisor.degree();
        let mut quotient_length = self
            .degree()
            .checked_sub(divisor_degree)
            .map(|value| value + 1)
            .unwrap_or(1);
        let mut quotient = vec![0_u64; quotient_length];
        let leading = divisor.coefficients[divisor_degree] as u128;
        let modulus = self.modulus as u128;
        while !remainder.is_zero() && remainder.degree() >= divisor_degree {
            let shift = remainder.degree() - divisor_degree;
            let inverse = mod_inverse(leading as u64, self.modulus);
            let scale = (u128::from(remainder.coefficients[remainder.degree()])
                * u128::from(inverse))
                % modulus;
            quotient[shift] = scale as u64;
            for index in 0..=divisor_degree {
                let product = (u128::from(divisor.coefficients[index]) * scale) % modulus;
                let target = shift + index;
                remainder.coefficients[target] =
                    ((u128::from(remainder.coefficients[target]) + modulus - product) % modulus)
                        as u64;
            }
            remainder.normalize();
        }
        quotient_length = quotient_length.max(1);
        (
            Self::new(quotient[..quotient_length].to_vec(), self.modulus),
            remainder,
        )
    }

    pub fn gcd(&self, other: &Self) -> Self {
        let mut left = self.clone();
        let mut right = other.clone();
        while !right.is_zero() {
            let remainder = left.divrem(&right).1;
            left = right;
            right = remainder;
        }
        if left.coefficients.is_empty() {
            return Self::one(self.modulus);
        }
        let leading = left.coefficients[left.degree()] as u128;
        let scale = mod_inverse(leading as u64, self.modulus);
        for coefficient in &mut left.coefficients {
            *coefficient =
                ((u128::from(*coefficient) * u128::from(scale)) % u128::from(self.modulus)) as u64;
        }
        left
    }

    pub fn eval(&self, point: u64) -> u64 {
        let mut result = 0_u64;
        for coefficient in self.coefficients.iter().rev() {
            result = ((u128::from(result) * u128::from(point) + u128::from(*coefficient))
                % u128::from(self.modulus)) as u64;
        }
        result
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|value| *value == 0)
    }

    pub fn roots(&self) -> Vec<u64> {
        let monic = self.monic();
        let degree = monic.degree();
        if degree == 0 {
            return Vec::new();
        }
        if degree == 1 {
            return monic.linear_root().into_iter().collect();
        }
        let common = monic.gcd(&monic.derivative());
        let radical = if common.is_constant() {
            monic
        } else {
            monic.divrem(&common).0
        };
        let x = ModPolynomial::new(vec![0, 1], self.modulus);
        let all_roots = x.pow_mod(self.modulus, &radical).sub(&x);
        let linear_part = radical.gcd(&all_roots).monic();
        let mut roots = Vec::new();
        split_linear_factors(&linear_part, &mut roots);
        roots.sort_unstable();
        roots.dedup();
        roots
    }

    fn linear_root(&self) -> Option<u64> {
        if self.degree() != 1 {
            return None;
        }
        let leading = *self.coefficients.get(1)?;
        if leading == 0 {
            return None;
        }
        let inverse = mod_inverse(leading, self.modulus);
        let width = u128::from(self.modulus);
        let constant = u128::from(*self.coefficients.first()?);
        Some(((width * width - constant * u128::from(inverse)) % width) as u64)
    }

    pub fn is_constant(&self) -> bool {
        self.degree() == 0
    }

    pub fn constant(&self) -> u64 {
        self.coefficients.first().copied().unwrap_or(0)
    }

    pub fn rem(&self, divisor: &Self) -> Self {
        self.divrem(divisor).1
    }

    pub fn monic(&self) -> Self {
        if self.is_zero() {
            return Self::zero(self.modulus);
        }
        let leading = self.coefficients.last().copied().unwrap_or(1);
        if leading == 1 {
            return self.clone();
        }
        let scale = mod_inverse(leading, self.modulus);
        let coefficients = self
            .coefficients
            .iter()
            .map(|coefficient| {
                ((u128::from(*coefficient) * u128::from(scale)) % u128::from(self.modulus)) as u64
            })
            .collect();
        Self::new(coefficients, self.modulus)
    }

    pub fn pow_mod(&self, mut exponent: u64, modulus: &Self) -> Self {
        if modulus.is_zero() {
            return Self::zero(self.modulus);
        }
        let mut result = Self::one(self.modulus);
        let mut base = self.rem(modulus);
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result.mul(&base).rem(modulus);
            }
            base = base.mul(&base).rem(modulus);
            exponent >>= 1;
        }
        result
    }

    pub fn derivative(&self) -> Self {
        let degree = self.degree();
        if degree == 0 {
            return Self::zero(self.modulus);
        }
        let mut coefficients = vec![0_u64; degree];
        for (index, coefficient) in coefficients.iter_mut().enumerate() {
            let raw = self.coefficients.get(index + 1).copied().unwrap_or(0);
            *coefficient = ((u128::from(raw) * u128::from((index + 1) as u64))
                % u128::from(self.modulus)) as u64;
        }
        Self::new(coefficients, self.modulus)
    }
}

pub fn null_space(matrix: &[Vec<u64>], columns: usize, prime: u64) -> Vec<Vec<u64>> {
    let width = u128::from(prime);
    let mut rows: Vec<Vec<u64>> = matrix.to_vec();
    for row in rows.iter_mut() {
        row.resize(columns, 0);
    }
    let mut pivots: Vec<Option<usize>> = vec![None; columns];
    let mut pivot_row = 0;
    for column in 0..columns {
        if pivot_row >= rows.len() {
            break;
        }
        let mut selected = None;
        for (offset, row) in rows.iter().enumerate().skip(pivot_row) {
            if row[column] != 0 {
                selected = Some(offset);
                break;
            }
        }
        let Some(row) = selected else {
            continue;
        };
        rows.swap(pivot_row, row);
        let inverse = mod_inverse(rows[pivot_row][column], prime);
        for value in rows[pivot_row].iter_mut() {
            *value = ((u128::from(*value) * u128::from(inverse)) % width) as u64;
        }
        let normalized = rows[pivot_row].clone();
        for (index, row) in rows.iter_mut().enumerate() {
            if index == pivot_row || row[column] == 0 {
                continue;
            }
            let factor = row[column];
            for (position, value) in row.iter_mut().enumerate().take(columns) {
                let left = u128::from(*value);
                let right = u128::from(factor) * u128::from(normalized[position]);
                *value = ((left + width * width - right) % width) as u64;
            }
        }
        pivots[column] = Some(pivot_row);
        pivot_row += 1;
    }
    let mut basis: Vec<Vec<u64>> = Vec::new();
    for free_column in 0..columns {
        if pivots[free_column].is_some() {
            continue;
        }
        let mut vector = vec![0_u64; columns];
        vector[free_column] = 1;
        for (column, pivot) in pivots.iter().enumerate() {
            if let Some(row) = pivot {
                let value = rows[*row][free_column];
                vector[column] = ((width - u128::from(value)) % width) as u64;
            }
        }
        basis.push(vector);
    }
    basis
}

fn split_linear_factors(polynomial: &ModPolynomial, roots: &mut Vec<u64>) {
    let degree = polynomial.degree();
    if degree == 0 {
        return;
    }
    if degree == 1 {
        if let Some(root) = polynomial.linear_root() {
            roots.push(root);
        }
        return;
    }
    if let Some((left, right)) = cauchy_split(polynomial) {
        split_linear_factors(&left, roots);
        split_linear_factors(&right, roots);
        return;
    }
    for point in 0..polynomial.modulus {
        if polynomial.eval(point) != 0 {
            continue;
        }
        let factor = ModPolynomial::new(
            vec![
                (polynomial.modulus - point % polynomial.modulus) % polynomial.modulus,
                1,
            ],
            polynomial.modulus,
        );
        roots.push(point);
        let mut remainder = polynomial.monic();
        for _ in 0..degree {
            if remainder.degree() == 0 {
                break;
            }
            remainder = remainder.divrem(&factor).1;
        }
        split_linear_factors(&remainder.monic(), roots);
        return;
    }
}

fn cauchy_split(polynomial: &ModPolynomial) -> Option<(ModPolynomial, ModPolynomial)> {
    let exponent = polynomial.modulus.div_ceil(2);
    for offset in 0..16_u64 {
        let probe = ModPolynomial::new(vec![offset, 1], polynomial.modulus);
        let reduced = probe.rem(polynomial);
        let image = probe.pow_mod(exponent, polynomial).sub(&reduced);
        let divisor = polynomial.gcd(&image);
        if !divisor.is_constant() && divisor.degree() < polynomial.degree() {
            return Some((divisor.monic(), polynomial.divrem(&divisor).0.monic()));
        }
    }
    None
}

pub fn factor_squarefree(polynomial: &ModPolynomial) -> Vec<ModPolynomial> {
    let prime = polynomial.modulus;
    let monic = polynomial.monic();
    let degree = monic.degree();
    if degree == 0 {
        return Vec::new();
    }
    if !monic.derivative().gcd(&monic).is_constant() {
        return Vec::new();
    }
    if degree == 1 {
        return vec![monic];
    }
    let x = ModPolynomial::new(vec![0, 1], prime);
    let x_to_p = x.pow_mod(prime, &monic);
    let mut matrix: Vec<Vec<u64>> = Vec::with_capacity(degree);
    let mut power = ModPolynomial::one(prime);
    let mut generator = ModPolynomial::one(prime);
    for _ in 0..degree {
        let image = power.sub(&generator);
        let mut row: Vec<u64> = image.coefficients.clone();
        row.resize(degree, 0);
        row.truncate(degree);
        for (index, coefficient) in row.iter_mut().enumerate() {
            let value = power.coefficients.get(index).copied().unwrap_or(0) % prime;
            let other = generator.coefficients.get(index).copied().unwrap_or(0) % prime;
            *coefficient = ((u128::from(*coefficient) + u128::from(prime) + u128::from(prime)
                - u128::from(value)
                - u128::from(other))
                % u128::from(prime)) as u64;
        }
        matrix.push(row);
        power = x_to_p.mul(&power).rem(&monic);
        generator = x.mul(&generator).rem(&monic);
    }
    let mut factors: Vec<ModPolynomial> = vec![monic.clone()];
    for vector in null_space(&matrix, degree, prime) {
        if vector.iter().all(|value| *value == 0) {
            continue;
        }
        let probe = ModPolynomial::new(vector, prime);
        let mut next: Vec<ModPolynomial> = Vec::new();
        for current in factors {
            let divisor = current.gcd(&probe);
            if divisor.is_constant() || divisor.degree() == current.degree() {
                next.push(current);
                continue;
            }
            let quotient = current.divrem(&divisor).0.monic();
            next.push(divisor.monic());
            next.push(quotient);
        }
        factors = next;
    }
    factors.retain(|factor| factor.degree() > 0);
    factors
}

pub fn is_irreducible(modulus: &ModPolynomial, prime: u64) -> bool {
    let degree = modulus.degree();
    if degree < 1 || modulus.is_zero() {
        return false;
    }
    let x = ModPolynomial::new(vec![0, 1], prime);
    let total = match prime.checked_pow(degree as u32) {
        Some(value) => value,
        None => return false,
    };
    if x.pow_mod(total, modulus) != x.rem(modulus) {
        return false;
    }
    let mut remaining = degree;
    let mut step = 2_usize;
    while step * step <= remaining {
        if remaining.is_multiple_of(step) {
            if !separable(&x, modulus, prime, degree, step) {
                return false;
            }
            while remaining.is_multiple_of(step) {
                remaining /= step;
            }
        }
        step += 1;
    }
    remaining <= 1 || separable(&x, modulus, prime, degree, remaining)
}

fn separable(
    x: &ModPolynomial,
    modulus: &ModPolynomial,
    prime: u64,
    degree: usize,
    factor: usize,
) -> bool {
    let reduced = degree / factor;
    let Some(exponent) = prime.checked_pow(reduced as u32) else {
        return false;
    };
    let shifted = x.pow_mod(exponent, modulus).sub(&x.rem(modulus));
    shifted.gcd(modulus).degree() == 0
}

pub fn crt_reconstruct(
    factors: &[ModPolynomial],
    residues: &[ModPolynomial],
) -> Option<ModPolynomial> {
    let mut modulus = factors.first()?.clone();
    let mut result = residues.first()?.rem(&modulus);
    for index in 1..factors.len() {
        let next = factors[index].clone();
        let target = residues[index].rem(&next);
        let reduced = modulus.rem(&next);
        let size = reduced.modulus.checked_pow(next.degree() as u32)?;
        let inverse = reduced.pow_mod(size.checked_sub(2)?, &next);
        let difference = target.sub(&result.rem(&next));
        let adjustment = difference.mul(&inverse).rem(&next);
        let combined = modulus.mul(&next);
        result = result.add(&modulus.mul(&adjustment)).rem(&combined);
        modulus = combined;
    }
    Some(result)
}

pub fn mod_inverse(value: u64, modulus: u64) -> u64 {
    let mut t = 0_i128;
    let mut new_t = 1_i128;
    let mut r = modulus as i128;
    let mut new_r = value as i128;
    while new_r != 0 {
        let quotient = r / new_r;
        let temp_t = t - quotient * new_t;
        t = new_t;
        new_t = temp_t;
        let temp_r = r - quotient * new_r;
        r = new_r;
        new_r = temp_r;
    }
    if t < 0 {
        t += modulus as i128;
    }
    (t as u64) % modulus
}

pub fn homogeneous(coefficients: &[u128], a: &u128, b: &u128) -> u128 {
    let degree = coefficients.len() - 1;
    let mut result = 0_u128;
    for coefficient in coefficients {
        result = result * a + coefficient * b;
    }
    let _ = degree;
    result
}
