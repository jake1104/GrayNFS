use crate::integer::Integer;

pub struct Montgomery {
    modulus: Vec<u64>,
    modulus_integer: Integer,
    limbs: usize,
    inv_prime: u64,
    r_squared: Vec<u64>,
}

fn limbs_of(value: &Integer, limbs: usize) -> Vec<u64> {
    let (_, bytes) = value.to_bytes_le();
    let mut out = vec![0_u64; limbs];
    for (index, chunk) in bytes.chunks(8).enumerate() {
        if index >= limbs {
            break;
        }
        let mut word = 0_u64;
        for (position, byte) in chunk.iter().enumerate() {
            word |= (*byte as u64) << (8 * position);
        }
        out[index] = word;
    }
    out
}

fn integer_of(limbs: &[u64]) -> Integer {
    let mut bytes = Vec::with_capacity(limbs.len() * 8);
    for word in limbs {
        bytes.extend_from_slice(&word.to_le_bytes());
    }
    Integer::from_bytes_le(num_bigint::Sign::Plus, &bytes)
}

impl Montgomery {
    pub fn new(modulus: &Integer) -> Option<Self> {
        if modulus <= &Integer::from(1) {
            return None;
        }
        if modulus % Integer::from(2) != Integer::from(1) {
            return None;
        }
        let (_, bytes) = modulus.to_bytes_le();
        let limbs = bytes.len().div_ceil(8);
        if limbs == 0 {
            return None;
        }
        let first = limbs_of(modulus, limbs)[0];
        let mut inverse = 1_u64;
        for _ in 0..6 {
            inverse = inverse.wrapping_mul(2u64.wrapping_sub(first.wrapping_mul(inverse)));
        }
        let shift = Integer::from(1) << (64 * limbs as u32);
        let r_squared = limbs_of(&(shift.clone() * shift % modulus), limbs);
        Some(Self {
            modulus: limbs_of(modulus, limbs),
            modulus_integer: modulus.clone(),
            limbs,
            inv_prime: inverse.wrapping_neg(),
            r_squared,
        })
    }

    pub fn limbs(&self) -> usize {
        self.limbs
    }

    pub fn convert_in(&self, value: &Integer) -> Vec<u64> {
        let reduced = limbs_of(&(value % &self.modulus_integer), self.limbs);
        self.multiply(&reduced, &self.r_squared)
    }

    pub fn convert_out(&self, limbs: &[u64]) -> Integer {
        let one = limbs_of(&Integer::from(1), self.limbs);
        integer_of(&self.multiply(limbs, &one))
    }

    pub fn multiply_scratch_len(&self) -> usize {
        2 * self.limbs + 2
    }

    pub fn multiply_into(&self, out: &mut [u64], left: &[u64], right: &[u64], middle: &mut [u64]) {
        debug_assert_eq!(out.len(), self.limbs);
        debug_assert_eq!(middle.len(), 2 * self.limbs + 2);
        match self.limbs {
            2 => self.cios_into::<2>(out, left, right, middle),
            3 => self.cios_into::<3>(out, left, right, middle),
            4 => self.cios_into::<4>(out, left, right, middle),
            5 => self.cios_into::<5>(out, left, right, middle),
            _ => self.cios_into_dyn(out, left, right, middle),
        }
    }

    fn cios_into<const LIMBS: usize>(
        &self,
        out: &mut [u64],
        left: &[u64],
        right: &[u64],
        middle: &mut [u64],
    ) {
        debug_assert_eq!(out.len(), LIMBS);
        debug_assert_eq!(self.limbs, LIMBS);
        debug_assert_eq!(middle.len(), 2 * LIMBS + 2);
        for slot in middle.iter_mut() {
            *slot = 0;
        }
        for i in 0..LIMBS {
            let mut carry = 0_u128;
            for j in 0..LIMBS {
                let total = middle[i + j] as u128 + (left[i] as u128) * (right[j] as u128) + carry;
                middle[i + j] = total as u64;
                carry = total >> 64;
            }
            let mut position = i + LIMBS;
            let mut extra = carry;
            while extra != 0 {
                let total = middle[position] as u128 + extra;
                middle[position] = total as u64;
                extra = total >> 64;
                position += 1;
            }
            let factor = middle[i].wrapping_mul(self.inv_prime);
            let mut second = 0_u128;
            for j in 0..LIMBS {
                let total =
                    middle[i + j] as u128 + (factor as u128) * (self.modulus[j] as u128) + second;
                middle[i + j] = total as u64;
                second = total >> 64;
            }
            let mut position = i + LIMBS;
            let mut extra = second;
            while extra != 0 {
                let total = middle[position] as u128 + extra;
                middle[position] = total as u64;
                extra = total >> 64;
                position += 1;
            }
        }
        debug_assert_eq!(middle[2 * LIMBS + 1], 0);
        debug_assert!(middle[2 * LIMBS] <= 1);
        for (index, slot) in out.iter_mut().enumerate() {
            *slot = middle[index + LIMBS];
        }
        if middle[2 * LIMBS] != 0 || Self::at_least(out, &self.modulus) {
            Self::subtract_in_place(out, &self.modulus);
        }
    }

    fn cios_into_dyn(&self, out: &mut [u64], left: &[u64], right: &[u64], middle: &mut [u64]) {
        debug_assert_eq!(out.len(), self.limbs);
        debug_assert_eq!(middle.len(), 2 * self.limbs + 2);
        for slot in middle.iter_mut() {
            *slot = 0;
        }
        let limbs = self.limbs;
        for i in 0..limbs {
            let mut carry = 0_u128;
            for j in 0..limbs {
                let total = middle[i + j] as u128 + (left[i] as u128) * (right[j] as u128) + carry;
                middle[i + j] = total as u64;
                carry = total >> 64;
            }
            let mut position = i + limbs;
            let mut extra = carry;
            while extra != 0 {
                let total = middle[position] as u128 + extra;
                middle[position] = total as u64;
                extra = total >> 64;
                position += 1;
            }
            let factor = middle[i].wrapping_mul(self.inv_prime);
            let mut second = 0_u128;
            for j in 0..limbs {
                let total =
                    middle[i + j] as u128 + (factor as u128) * (self.modulus[j] as u128) + second;
                middle[i + j] = total as u64;
                second = total >> 64;
            }
            let mut position = i + limbs;
            let mut extra = second;
            while extra != 0 {
                let total = middle[position] as u128 + extra;
                middle[position] = total as u64;
                extra = total >> 64;
                position += 1;
            }
        }
        debug_assert_eq!(middle[2 * limbs + 1], 0);
        debug_assert!(middle[2 * limbs] <= 1);
        for (index, slot) in out.iter_mut().enumerate() {
            *slot = middle[index + limbs];
        }
        if middle[2 * limbs] != 0 || Self::at_least(out, &self.modulus) {
            Self::subtract_in_place(out, &self.modulus);
        }
    }

    pub fn multiply(&self, left: &[u64], right: &[u64]) -> Vec<u64> {
        let mut middle = vec![0_u64; self.multiply_scratch_len()];
        let mut out = vec![0_u64; self.limbs];
        self.multiply_into(&mut out, left, right, &mut middle);
        out
    }

    fn at_least(left: &[u64], right: &[u64]) -> bool {
        for (a, b) in left.iter().rev().zip(right.iter().rev()) {
            if a != b {
                return a > b;
            }
        }
        true
    }

    fn subtract_in_place(left: &mut [u64], right: &[u64]) {
        let mut borrow = 0_i128;
        for (a, b) in left.iter_mut().zip(right.iter()) {
            let total = *a as i128 - *b as i128 - borrow;
            if total < 0 {
                *a = (total + (1_i128 << 64)) as u64;
                borrow = 1;
            } else {
                *a = total as u64;
                borrow = 0;
            }
        }
    }

    pub fn add_into(&self, out: &mut [u64], left: &[u64], right: &[u64]) {
        debug_assert_eq!(out.len(), self.limbs);
        let mut carry = 0_u64;
        for index in 0..self.limbs {
            let total = left[index] as u128 + right[index] as u128 + carry as u128;
            out[index] = total as u64;
            carry = (total >> 64) as u64;
        }
        if carry != 0 || Self::at_least(out, &self.modulus) {
            Self::subtract_in_place(out, &self.modulus);
        }
    }

    pub fn add(&self, left: &[u64], right: &[u64]) -> Vec<u64> {
        let mut out = vec![0_u64; self.limbs];
        self.add_into(&mut out, left, right);
        out
    }

    pub fn subtract_into(&self, out: &mut [u64], left: &[u64], right: &[u64]) {
        debug_assert_eq!(out.len(), self.limbs);
        let mut borrow = 0_i128;
        for index in 0..self.limbs {
            let total = left[index] as i128 - right[index] as i128 - borrow;
            if total < 0 {
                out[index] = (total + (1_i128 << 64)) as u64;
                borrow = 1;
            } else {
                out[index] = total as u64;
                borrow = 0;
            }
        }
        if borrow != 0 {
            let mut carry = 0_u64;
            for (slot, m) in out.iter_mut().zip(self.modulus.iter()) {
                let total = *slot as u128 + *m as u128 + carry as u128;
                *slot = total as u64;
                carry = (total >> 64) as u64;
            }
        }
    }

    pub fn subtract(&self, left: &[u64], right: &[u64]) -> Vec<u64> {
        let mut out = vec![0_u64; self.limbs];
        self.subtract_into(&mut out, left, right);
        out
    }
}

#[cfg(test)]
mod montgomery_tests {
    use super::*;
    use num_traits::Euclid;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn random_integer(state: &mut u64, bits: u32) -> Integer {
        let mut acc = Integer::from(0);
        let mut remaining = bits;
        while remaining > 0 {
            let take = remaining.min(64);
            let word = xorshift(state) & (u64::MAX >> (64 - take));
            acc = (acc << take) + Integer::from(word);
            remaining -= take;
        }
        acc
    }

    #[test]
    fn operations_match_reference() {
        let mut state = 0x243F6A8885A308D3_u64;
        for modulus_bits in [3_u32, 17, 65, 127, 129, 267, 300] {
            let modulus = loop {
                let candidate = random_integer(&mut state, modulus_bits) | Integer::from(1);
                if candidate > Integer::from(1) {
                    break candidate;
                }
            };
            let field = Montgomery::new(&modulus).expect("odd modulus");
            for _ in 0..200 {
                let left = random_integer(&mut state, modulus_bits + 2);
                let right = random_integer(&mut state, modulus_bits + 2);
                let left_mont = field.convert_in(&left);
                let right_mont = field.convert_in(&right);
                let product = field.convert_out(&field.multiply(&left_mont, &right_mont));
                assert_eq!(product, (left.clone() * &right) % &modulus);
                let sum = field.convert_out(&field.add(&left_mont, &right_mont));
                assert_eq!(sum, (left.clone() + &right) % &modulus);
                let difference = field.convert_out(&field.subtract(&left_mont, &right_mont));
                assert_eq!(difference, (left.clone() - &right).rem_euclid(&modulus));
            }
        }
    }

    #[test]
    fn unrolled_n34_matches_dyn() {
        let mut state = 0x9E3779B97F4A7C15_u64;
        for modulus_bits in [129_u32, 140, 160, 192, 193, 200, 250, 256] {
            let modulus = random_integer(&mut state, modulus_bits)
                | (Integer::from(1) << (modulus_bits - 1))
                | Integer::from(1);
            let field = Montgomery::new(&modulus).expect("odd modulus");
            let limbs = field.limbs();
            assert!(limbs == 3 || limbs == 4);
            for _ in 0..100 {
                let left_raw = random_integer(&mut state, modulus_bits + 2) % &modulus;
                let right_raw = random_integer(&mut state, modulus_bits + 2) % &modulus;
                let left_mont = field.convert_in(&left_raw);
                let right_mont = field.convert_in(&right_raw);
                let mut out_kernel = vec![0_u64; limbs];
                let mut out_dyn = vec![0_u64; limbs];
                let mut mid_kernel = vec![0_u64; 2 * limbs + 2];
                let mut mid_dyn = vec![0_u64; 2 * limbs + 2];
                match limbs {
                    3 => field.cios_into::<3>(
                        &mut out_kernel,
                        &left_mont,
                        &right_mont,
                        &mut mid_kernel,
                    ),
                    4 => field.cios_into::<4>(
                        &mut out_kernel,
                        &left_mont,
                        &right_mont,
                        &mut mid_kernel,
                    ),
                    _ => unreachable!(),
                }
                field.cios_into_dyn(&mut out_dyn, &left_mont, &right_mont, &mut mid_dyn);
                assert_eq!(out_kernel, out_dyn);
                let dispatch = field.multiply(&left_mont, &right_mont);
                assert_eq!(dispatch, out_dyn);
            }
        }
    }

    #[test]
    fn dispatch_covers_n34() {
        let src: String = include_str!("montgomery.rs")
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect();
        for limbs in [3, 4] {
            let arm: String = format!("{limbs}=>self.cios_into::<{limbs}>(out,left,right,middle)");
            assert!(src.contains(&arm), "missing dispatch arm {limbs}");
        }
    }

    #[test]
    fn conversions_round_trip() {
        let mut state = 0x452821E638D01377_u64;
        let modulus = random_integer(&mut state, 267) | Integer::from(1);
        let field = Montgomery::new(&modulus).expect("odd modulus");
        for _ in 0..100 {
            let value = random_integer(&mut state, 300);
            let back = field.convert_out(&field.convert_in(&value));
            assert_eq!(back, value % &modulus);
        }
    }

    #[test]
    fn rejects_even_modulus() {
        assert!(Montgomery::new(&Integer::from(100)).is_none());
        assert!(Montgomery::new(&Integer::from(1)).is_none());
    }
}
