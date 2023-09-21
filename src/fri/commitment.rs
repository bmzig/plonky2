use ff::PrimeField;
use blake3::{Hasher};
use primitive_types::U256;

use crate::{
    polynomial::Polynomial,
    fri::{FriCommitment},
    fft::serial,
    constants::*,
    domains::Domain,
};

impl<F: PrimeField> Polynomial<F> {

    pub(crate) fn commitment(&self) -> FriCommitment<F> {

        // FRI commitment is evaluation of a polynomial across the dp-th roots of unity, where d is
        // the degree of the polynomial and p is the FRI_BLOWUP_FACTOR constant.

        let mut evaluations = self.coefficients();
        let log_n = {
            let mut x = 0usize;
            let mut y = evaluations.len();
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            x
        };
        let extended_log_n = log_n + FRI_BLOWUP_LOG;
        let omega = Domain::root_with_order_unchecked((FRI_BLOWUP_FACTOR * evaluations.len()) as u64);
        
        evaluations.append(&mut vec![F::ZERO; (1<<extended_log_n) - (1<<log_n)]);
        serial::serial_fft(evaluations.as_mut_slice(), &omega, extended_log_n as u32);

        // Create the merkle tree out of the evaluations of the dp-th roots of unity:
        // f(w^0) --
        //          | H(f1|f2) --
        // f(w^1) --            |
        //                      | H(H(f0|f1)|H(f2|f3)) --
        // f(w^2) --            |                       |
        //          | H(f2|f3) --                       |
        // f(w^3) --                                    |
        //                                              | Com(f)
        // f(w^4) --                                    |-------
        //          | H(f4|f5) --                       |
        // f(w^5) --            |                       |
        //                      | H(H(f4|f5)|H(f6|f7)) --
        // f(w^6) --            |
        //          | H(f6|f7) -- 
        // f(w^7) -- 

        let mut hash_vector = Vec::new();
        for i in (0..evaluations.len()).step_by(2) {
            let mut hash = Hasher::new();
            hash.update(evaluations[i].to_repr().as_ref());
            hash.update(evaluations[i+1].to_repr().as_ref());
            hash_vector.push(hash.finalize());
        }

        for _ in 0..(extended_log_n-1) {
            let mut new_hash_vector = Vec::new();
            for i in (0..hash_vector.len()).step_by(2) {
                let mut hasher = Hasher::new();
                hasher.update(hash_vector[i].as_bytes().as_slice());
                hasher.update(hash_vector[i+1].as_bytes().as_slice());
                new_hash_vector.push(hasher.finalize());
            }
            hash_vector = new_hash_vector;
        }
        assert!(hash_vector.len() == 1);
        FriCommitment::new(hash_vector[0])
    }
}

impl<F: PrimeField> FriCommitment<F> {
    pub fn interpret_as_element(&self) -> F {
        crate::utils::field_element_from_bytes(self.value().as_bytes().as_slice())
    }

    pub fn interpret_as_root_of_unity(&self, domain_size: u64) -> F {
        let random_exponent = U256::from_big_endian(self.value().as_bytes().as_slice()).low_u64();
        let base: F = crate::domains::Domain::root_with_order_unchecked(domain_size);
        base.pow([random_exponent])
    }
}
