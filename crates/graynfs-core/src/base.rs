use std::collections::HashMap;

use crate::integer::{Integer, IntegerExt};

use crate::error::FactorError;
use crate::polynomial::PolynomialPair;
use crate::prime::{mod_sqrt, sieve_primes};

#[derive(Clone, Debug)]
pub struct AlgebraicPrime {
    pub p: u64,
    pub root: u64,
    pub id: usize,
}

#[derive(Clone, Debug)]
pub struct FactorBase {
    pub rational_primes: Vec<u64>,
    pub rational_indices: HashMap<u64, usize>,
    pub algebraic_primes: Vec<AlgebraicPrime>,
    pub algebraic_by_prime: HashMap<u64, Vec<usize>>,
    pub algebraic_modulus_primes: Vec<u64>,
    pub sign_column: usize,
}

impl FactorBase {
    pub fn build(
        polynomial: &PolynomialPair,
        n: &Integer,
        rational_bound: usize,
        algebraic_bound: usize,
    ) -> Result<Self, FactorError> {
        let rational_primes: Vec<u64> = sieve_primes(rational_bound)
            .into_iter()
            .filter(|prime| n % *prime != Integer::from(0))
            .collect();
        let mut rational_indices = HashMap::with_capacity(rational_primes.len());
        for (index, prime) in rational_primes.iter().enumerate() {
            rational_indices.insert(*prime, index);
        }
        let mut algebraic_primes = Vec::new();
        let mut algebraic_by_prime: HashMap<u64, Vec<usize>> = HashMap::new();
        let mut algebraic_modulus_primes = Vec::new();
        let constant = polynomial.constant();
        for prime in sieve_primes(algebraic_bound) {
            if n % prime == Integer::from(0) {
                continue;
            }
            algebraic_modulus_primes.push(prime);
            let constant_mod = (&constant % Integer::from(prime)).to_u64_saturating();
            let residue = (prime - constant_mod) % prime;
            let Some(root) = mod_sqrt(residue, prime) else {
                continue;
            };
            let roots = if root == 0 {
                vec![0]
            } else {
                let other = prime - root;
                if other == root {
                    vec![root]
                } else {
                    vec![root, other]
                }
            };
            for root in roots {
                let id = rational_primes.len() + algebraic_primes.len();
                algebraic_by_prime
                    .entry(prime)
                    .or_default()
                    .push(algebraic_primes.len());
                algebraic_primes.push(AlgebraicPrime { p: prime, root, id });
            }
        }
        let sign_column = rational_primes.len() + algebraic_primes.len();
        Ok(Self {
            rational_primes,
            rational_indices,
            algebraic_primes,
            algebraic_by_prime,
            algebraic_modulus_primes,
            sign_column,
        })
    }

    pub fn column_count(&self) -> usize {
        self.sign_column + 1
    }

    pub fn toggle_rational(&self, mask: &mut [u64], prime: u64) -> Result<bool, FactorError> {
        let Some(&index) = self.rational_indices.get(&prime) else {
            return Ok(false);
        };
        mask[index / 64] ^= 1_u64 << (index % 64);
        Ok(true)
    }

    pub fn toggle_algebraic(&self, mask: &mut [u64], prime: u64, a: i64, b: i64) {
        self.toggle_algebraic_u64(mask, prime, a as u64, b as u64);
    }

    pub fn toggle_algebraic_u64(&self, mask: &mut [u64], prime: u64, a: u64, b: u64) {
        let Some(indices) = self.algebraic_by_prime.get(&prime) else {
            return;
        };
        let a_mod = a % prime;
        let b_mod = b % prime;
        for index in indices {
            let root = self.algebraic_primes[*index].root;
            if (a_mod + prime - b_mod * root % prime).is_multiple_of(prime) {
                let column = self.algebraic_primes[*index].id;
                mask[column / 64] ^= 1_u64 << (column % 64);
            }
        }
    }

    pub fn toggle_algebraic_bigint(&self, mask: &mut [u64], prime: u64, a: &Integer, b: &Integer) {
        let Some(indices) = self.algebraic_by_prime.get(&prime) else {
            return;
        };
        let a_mod = (a % Integer::from(prime)).to_u64_saturating();
        let b_mod = (b % Integer::from(prime)).to_u64_saturating();
        for index in indices {
            let root = self.algebraic_primes[*index].root;
            if (a_mod + prime - b_mod * root % prime).is_multiple_of(prime) {
                let column = self.algebraic_primes[*index].id;
                mask[column / 64] ^= 1_u64 << (column % 64);
            }
        }
    }

    pub fn contains_prime(&self, prime: u64) -> bool {
        self.rational_indices.contains_key(&prime)
    }

    pub fn has_algebraic_root(&self, prime: u64) -> bool {
        self.algebraic_by_prime.contains_key(&prime)
    }

    pub fn for_siqs(n: &Integer, rational_bound: usize, algebraic_bound: usize) -> Self {
        let rational_primes: Vec<u64> = sieve_primes(rational_bound)
            .into_iter()
            .filter(|prime| n % *prime != Integer::from(0))
            .collect();
        let mut rational_indices = HashMap::with_capacity(rational_primes.len());
        for (index, prime) in rational_primes.iter().enumerate() {
            rational_indices.insert(*prime, index);
        }
        let mut algebraic_primes = Vec::new();
        let mut algebraic_by_prime: HashMap<u64, Vec<usize>> = HashMap::new();
        let mut algebraic_modulus_primes = Vec::new();
        for prime in sieve_primes(algebraic_bound) {
            if n % Integer::from(prime) == Integer::from(0) {
                continue;
            }
            let residue = (n % Integer::from(prime)).to_u64_saturating();
            let Some(root) = mod_sqrt(residue, prime) else {
                continue;
            };
            algebraic_modulus_primes.push(prime);
            let roots = if root == 0 {
                vec![0]
            } else {
                let other = prime - root;
                if other == root {
                    vec![root]
                } else {
                    vec![root, other]
                }
            };
            for root in roots {
                let id = rational_primes.len() + algebraic_primes.len();
                algebraic_by_prime
                    .entry(prime)
                    .or_default()
                    .push(algebraic_primes.len());
                algebraic_primes.push(AlgebraicPrime { p: prime, root, id });
            }
        }
        let sign_column = rational_primes.len() + algebraic_primes.len();
        Self {
            rational_primes,
            rational_indices,
            algebraic_primes,
            algebraic_by_prime,
            algebraic_modulus_primes,
            sign_column,
        }
    }
}
