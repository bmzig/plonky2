use ff::PrimeField;
use blake3::{Hasher, Hash};

use crate::{
    polynomial::Polynomial,
    fri::{FriCommitment, FriChallenge},
    constants::*,
};

impl<F: PrimeField> Polynomial<F> {

    pub fn fold_full(&self) -> (Vec<FriCommitment<F>>, Vec<Self>) {
        let mut target_length = self.len()/2;
        let log_n = {
            let mut x = 0usize;
            let mut y = self.len();
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            x
        };

        let mut commitment_vector = Vec::with_capacity(log_n);
        let mut polynomial_vector = Vec::with_capacity(log_n);

        let mut com = self.commitment();
        let mut r: F = com.interpret_as_element();

        let mut folded = vec![F::ZERO; target_length];
        let mut c = 0usize;
        for element in folded.iter_mut().take(target_length) { *element = self.coefficient_at(c) + (self.coefficient_at(c+1) * r); c += 2; }

        let mut intermediate = Polynomial::from_vec(folded);
        for _i in 0..(log_n-1) {
            polynomial_vector.push(intermediate.clone());

            com = intermediate.commitment();
            r = com.interpret_as_element();
            target_length = intermediate.len()/2;
            let mut folded = vec![F::ZERO; target_length];
            let mut c = 0usize;
            for element in folded.iter_mut().take(target_length) { *element = intermediate.coefficient_at(c) + (intermediate.coefficient_at(c+1) * r); c += 2; }
            intermediate = Polynomial::from_vec(folded);
            
            commitment_vector.push(com.clone());
        }
        commitment_vector.push(intermediate.commitment());
        polynomial_vector.push(intermediate);

        (commitment_vector, polynomial_vector)
    }
}

impl<F: PrimeField> FriChallenge<F> {

    pub(crate) fn query_check(&self, top_commitment: &FriCommitment<F>, random_root_of_unity: &F) -> Hash {

        let mut target = random_root_of_unity.square();
        let alpha: F = top_commitment.interpret_as_element();
        let even = (self.positive_evaluation() + self.negative_evaluation()) * F::from(2).invert().unwrap();
        let odd = (self.positive_evaluation() - self.negative_evaluation()) * (F::from(2) * random_root_of_unity).invert().unwrap();
        let mut assembled = even + (alpha * odd);

        for i in 0..self.fold_queries().len() {

            let alpha: F = self.commitment_vector()[i].interpret_as_element();
            let even = (assembled + self.fold_queries()[i]) * F::from(2).invert().unwrap();
            let odd = (assembled - self.fold_queries()[i]) * (F::from(2) * target).invert().unwrap();
            assembled = even + (alpha * odd);
            target = target.square()
        }

        let mut evals: [Hash; FRI_BLOWUP_FACTOR/2] = [Hash::from(ZERO_BYTES); FRI_BLOWUP_FACTOR/2];
        for eval in evals.iter_mut().take(FRI_BLOWUP_FACTOR/2) {
            let mut hasher = Hasher::new();
            hasher.update(assembled.to_repr().as_ref());
            hasher.update(assembled.to_repr().as_ref());
            *eval = hasher.finalize();
        }

        for _ in 0..FRI_BLOWUP_LOG-1 {
            for (c, i) in (0..evals.len()).step_by(2).enumerate() {
                let mut hasher = Hasher::new();
                hasher.update(evals[i].as_bytes().as_slice());
                hasher.update(evals[i+1].as_bytes().as_slice());
                evals[c] = hasher.finalize();
            }
        }

        evals[0]
    }
}
