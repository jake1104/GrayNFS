use num_traits::Zero;

use crate::integer::Integer;

const NORM_EXPANSION_LIMIT: usize = 8;

fn determinant(matrix: &[Vec<Integer>]) -> Integer {
    let size = matrix.len();
    if size == 0 {
        return Integer::from(1);
    }
    if size == 1 {
        return matrix[0][0].clone();
    }
    if size == 2 {
        return &matrix[0][0] * &matrix[1][1] - &matrix[0][1] * &matrix[1][0];
    }
    let mut total = Integer::from(0);
    for column in 0..size {
        let minor: Vec<Vec<Integer>> = matrix
            .iter()
            .skip(1)
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(index, _)| *index != column)
                    .map(|(_, entry)| entry.clone())
                    .collect()
            })
            .collect();
        let cofactor = determinant(&minor);
        if column % 2 == 0 {
            total += &matrix[0][column] * &cofactor;
        } else {
            total -= &matrix[0][column] * &cofactor;
        }
    }
    total
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Element {
    pub coefficients: Vec<Integer>,
}

impl Element {
    pub fn zero(degree: usize) -> Self {
        Self {
            coefficients: vec![Integer::from(0); degree],
        }
    }

    pub fn one(degree: usize) -> Self {
        let mut coefficients = vec![Integer::from(0); degree];
        if !coefficients.is_empty() {
            coefficients[0] = Integer::from(1);
        }
        Self { coefficients }
    }

    pub fn from_coefficients(coefficients: Vec<Integer>) -> Self {
        Self { coefficients }
    }
}

#[derive(Clone, Debug)]
pub struct NumberField {
    pub degree: usize,
    pub polynomial: Vec<Integer>,
}

impl NumberField {
    pub fn new(degree: usize, polynomial: Vec<Integer>) -> Option<Self> {
        if degree < 2 || polynomial.len() != degree + 1 {
            return None;
        }
        if polynomial[degree] != Integer::from(1) {
            return None;
        }
        Some(Self { degree, polynomial })
    }

    pub fn one(&self) -> Element {
        Element::one(self.degree)
    }

    pub fn multiply(&self, left: &Element, right: &Element) -> Element {
        let degree = self.degree;
        let mut convolution = vec![Integer::from(0); degree * 2 - 1];
        for i in 0..degree {
            for j in 0..degree {
                convolution[i + j] += &left.coefficients[i] * &right.coefficients[j];
            }
        }
        for index in (degree..convolution.len()).rev() {
            let value = convolution[index].clone();
            if value.is_zero() {
                continue;
            }
            for j in 0..degree {
                convolution[index - degree + j] -= &value * &self.polynomial[j];
            }
        }
        convolution.truncate(degree);
        Element::from_coefficients(convolution)
    }

    pub fn evaluate(&self, element: &Element, point: &Integer) -> Integer {
        let mut result = Integer::from(0);
        for coefficient in element.coefficients.iter().rev() {
            result = result * point + coefficient;
        }
        result
    }

    pub fn homogeneous(&self, a: &Integer, b: &Integer) -> Integer {
        let degree = self.polynomial.len() - 1;
        let mut result = self.polynomial[degree].clone();
        let mut b_power = Integer::from(1);
        for index in (0..degree).rev() {
            b_power *= b;
            result = result * a + &self.polynomial[index] * &b_power;
        }
        result
    }

    pub fn norm(&self, element: &Element) -> Integer {
        let degree = self.degree;
        if element.coefficients.len() != degree || degree > NORM_EXPANSION_LIMIT {
            return Integer::from(0);
        }
        let mut matrix: Vec<Vec<Integer>> = Vec::with_capacity(degree);
        for row in 0..degree {
            let mut entries = Vec::with_capacity(degree);
            for column in 0..degree {
                let mut basis = vec![Integer::from(0); degree];
                basis[column] = Integer::from(1);
                entries.push(
                    self.multiply(&Element::from_coefficients(basis), element)
                        .coefficients[row]
                        .clone(),
                );
            }
            matrix.push(entries);
        }
        determinant(&matrix)
    }

    pub fn discriminant(&self) -> Integer {
        let degree = Integer::from(self.degree);
        let constant = self.polynomial[0].clone();
        let mut discriminant = degree.pow(self.degree as u32);
        discriminant *= constant.pow((self.degree - 1) as u32);
        if self.degree % 2 == 1 {
            discriminant = -discriminant;
        }
        discriminant
    }
}
