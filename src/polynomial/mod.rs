use ff::PrimeField;

use crate::{
    domains::Domain,
};

pub mod arithmetic;

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial<F: PrimeField>(Vec<F>);

// Construc impl block
impl<F: PrimeField> Polynomial<F> {

    // Assumes that v.len() is a power of 2. 
    pub(crate) fn from_vec(v: Vec<F>) -> Self {
        Self(v)
    }
 
    pub fn eval_single(&self, point: &F) -> F {
        let mut result: F = self.leading_coefficient();
        for i in (0..self.len() - 1).rev() {
            result = (result * *point) + self.0[i];
        }
        result
    }

    pub(crate) fn pop_zeros(&mut self) {
        for i in (0..self.len()).rev() {
            if self.0[i] == F::ZERO {
                self.0.pop();
            }
            else {
                break;
            }
        }
        if !self.len().is_power_of_two() { 
            let mut x = 0usize;
            let mut y = self.len();
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            self.pad_to_base(x + 1); 
        }
    }

    pub fn vanishing_polynomial(roots: u128) -> Self {
        let mut v = Vec::new();
        v.push(F::ZERO - F::ONE);
        let mut len = vec![F::ZERO; (roots - 1) as usize];
        v.append(&mut len);
        v.push(F::ONE);

        // Pad the vector to the next largest power of 2
        let mut additional_zeros = vec![F::ZERO; v.len().next_power_of_two() - v.len()];
        v.append(&mut additional_zeros);

        Self::from_vec(v)
    }

    pub fn pad_to_base(&mut self, new_base: usize) {
        let size = 1<<new_base;
        assert!(size > self.len());
        let mut extension = vec![F::ZERO; size - self.len()];
        self.0.append(&mut extension);
    }
}

// Utility/organizational functions
impl<F: PrimeField> Polynomial<F> {

    // Dealing with coefficients of the polynomial:
    pub fn coefficients(&self) -> Vec<F> {
        self.0.clone()
    }

    pub(crate) fn leading_coefficient(&self) -> F {
        let index = self.len();
        self.0[index-1]
    }

    pub(crate) fn coefficient_at(&self, index: usize) -> F {
        self.0[index]
    }

    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn log_n(&self) -> usize {
        let mut x = 0usize;
        let mut y = self.len().next_power_of_two();
        while y != 1 {
            y >>= 1;
            x += 1;
        }
        x
    }
}

// Testing impl block
impl<F: PrimeField> Polynomial<F> {

    // O(N^2). For testing only
    pub fn eval_at_naive(&self, x: &F) -> F {
        let mut e = *x;
        let mut ret = self.0[0];
        for i in 1..self.len() {
            ret += e * self.0[i];
            e *= *x;
        }
        ret
    }

    pub fn divisor_polynomial(roots: u64) -> Self {
        let gen_base: F = Domain::new_for_size(roots).unwrap().generator;
        let mut polynomial = Polynomial::from_vec(vec![F::ZERO - F::ONE, F::ONE]); // (x - 1)
        let mut gen = gen_base;

        for _ in 1..roots {
            let new_poly = Polynomial::from_vec(vec![F::ZERO - gen, F::ONE]); // (x - g^i)
            polynomial = polynomial * new_poly;
            gen *= gen_base;
        }
        polynomial.pop_zeros();
        polynomial
    }

}
