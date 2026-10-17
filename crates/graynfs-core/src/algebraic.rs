use crate::integer::Integer;

use num_traits::Zero;

use crate::prime::integer_sqrt;

#[derive(Clone, Debug)]
pub struct QuadraticElement {
    pub x: Integer,
    pub y: Integer,
}

#[derive(Clone, Debug)]
pub struct QuadraticRing {
    pub constant: Integer,
    pub linear: Integer,
}

impl QuadraticRing {
    pub fn new(constant: Integer) -> Self {
        Self {
            constant,
            linear: Integer::from(0),
        }
    }

    pub fn with_linear(constant: Integer, linear: Integer) -> Self {
        Self { constant, linear }
    }

    pub fn discriminant(&self) -> Integer {
        &self.linear * &self.linear - Integer::from(4) * &self.constant
    }

    pub fn element(&self, x: Integer, y: Integer) -> QuadraticElement {
        QuadraticElement { x, y }
    }

    pub fn multiply(&self, left: &QuadraticElement, right: &QuadraticElement) -> QuadraticElement {
        let cross = &left.x * &right.y + &left.y * &right.x;
        let product = &left.y * &right.y;
        QuadraticElement {
            x: &left.x * &right.x - &self.constant * &product,
            y: cross - &self.linear * &product,
        }
    }

    pub fn norm(&self, value: &QuadraticElement) -> Integer {
        &value.x * &value.x - &self.linear * (&value.x * &value.y)
            + &self.constant * (&value.y * &value.y)
    }

    pub fn evaluate(&self, value: &QuadraticElement, point: &Integer) -> Integer {
        &value.x + &(value.y.clone() * point)
    }

    pub fn square_root(&self, value: &QuadraticElement) -> Option<QuadraticElement> {
        if self.linear.is_zero() {
            return self.square_root_diagonal(value);
        }
        self.square_root_general(value)
    }

    fn square_root_diagonal(&self, value: &QuadraticElement) -> Option<QuadraticElement> {
        let norm = self.norm(value);
        if norm < Integer::from(0) {
            return None;
        }
        let norm_root = integer_sqrt(&norm);
        if &norm_root * &norm_root != norm {
            return None;
        }
        let x = &value.x;
        let y = &value.y;
        if y.is_zero() {
            if x < &Integer::from(0) {
                return None;
            }
            let root = integer_sqrt(x);
            if &root * &root != *x {
                return None;
            }
            return Some(QuadraticElement {
                x: root,
                y: Integer::from(0),
            });
        }
        let two = Integer::from(2);
        let first = (&norm_root + x) / &two;
        let second_numerator = &norm_root - x;
        let second_denominator = &two * &self.constant;
        if first < Integer::from(0)
            || second_numerator < Integer::from(0)
            || second_denominator <= Integer::from(0)
        {
            return None;
        }
        if (&norm_root + x) % &two != Integer::from(0)
            || &second_numerator % &second_denominator != Integer::from(0)
        {
            return None;
        }
        let second = second_numerator / &second_denominator;
        let mut u = integer_sqrt(&first);
        let mut v = integer_sqrt(&second);
        if &u * &u != first || &v * &v != second {
            return None;
        }
        if &(&u * &v) * &two != *y {
            v = -v;
        }
        if &(&u * &v) * &two != *y {
            return None;
        }
        if u < Integer::from(0) {
            u = -u;
        }
        Some(QuadraticElement { x: u, y: v })
    }

    fn square_root_general(&self, value: &QuadraticElement) -> Option<QuadraticElement> {
        let norm = self.norm(value);
        if norm < Integer::from(0) {
            return None;
        }
        let norm_root = integer_sqrt(&norm);
        if &norm_root * &norm_root != norm {
            return None;
        }
        let discriminant = self.discriminant();
        if discriminant.is_zero() {
            return None;
        }
        let x = &value.x;
        let y = &value.y;
        if y.is_zero() {
            let two = Integer::from(2);
            if x % &two != Integer::from(0) {
                return None;
            }
            let half = x / &two;
            if half < Integer::from(0) {
                return None;
            }
            let root = integer_sqrt(&half);
            if &root * &root != half {
                return None;
            }
            return Some(QuadraticElement {
                x: &root * &two,
                y: Integer::from(0),
            });
        }
        for sign in [Integer::from(1), Integer::from(-1)] {
            let numerator =
                Integer::from(2) * x - &self.linear * y - Integer::from(2) * (&norm_root * &sign);
            if &numerator % &discriminant != Integer::from(0) {
                continue;
            }
            let squared_b = numerator / &discriminant;
            if squared_b < Integer::from(0) {
                continue;
            }
            let b_root = integer_sqrt(&squared_b);
            if &b_root * &b_root != squared_b || b_root.is_zero() {
                continue;
            }
            for b_sign in [Integer::from(1), Integer::from(-1)] {
                let b = &b_root * &b_sign;
                let two_b = Integer::from(2) * &b;
                let a_numerator = y + &self.linear * (&b * &b);
                if &a_numerator % &two_b != Integer::from(0) {
                    continue;
                }
                let a = a_numerator / &two_b;
                if &a * &a - &self.constant * (&b * &b) != *x {
                    continue;
                }
                if Integer::from(2) * (&a * &b) - &self.linear * (&b * &b) != *y {
                    continue;
                }
                let a = if a < Integer::from(0) { -a } else { a };
                return Some(QuadraticElement { x: a, y: b });
            }
        }
        None
    }
}
