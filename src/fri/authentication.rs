use crate::{
    fri::{AuthenticationPath, AuthenticationHash},
    fft::serial,
    domains::Domain,
    constants::*,
    polynomial::Polynomial,
};

use ff::PrimeField;
use blake3::{Hasher, Hash};

impl<F: PrimeField> AuthenticationPath<F> {

    pub fn new(first_evaluation: F, second_evaluation: F, authentication_path: Vec<AuthenticationHash>) -> Self {
        Self {
            first_evaluation,
            second_evaluation,
            authentication_path,
        }
    }

    pub fn derive_root(&self) -> Hash {
        let mut hasher = Hasher::new();
        hasher.update(self.first_evaluation.to_repr().as_ref());
        hasher.update(self.second_evaluation.to_repr().as_ref());
        let mut target = hasher.finalize();

        /*
         * Reconstructs the root from log_d hashes:
         *           H2 --
         *               |
         *               | -- H3 --
         * F1 --         |        |
         *      | -- H1 --        |
         * F2 --                  |-- Com(f(x))
         *                        |
         *                        |
         *                    H4 --
        */

        for i in 0..self.authentication_path.len() {
            let mut hasher = Hasher::new();
            if self.authentication_path[i].is_first {
                hasher.update(self.authentication_path[i].hash_ref());
                hasher.update(target.as_bytes().as_slice());
            }
            else {
                hasher.update(target.as_bytes().as_slice());
                hasher.update(self.authentication_path[i].hash_ref());
            }
            target = hasher.finalize();
        }
        target
    }

    pub(crate) fn contains_evaluation(&self, evaluation: &F) -> bool {
        (self.first_evaluation == *evaluation) || (self.second_evaluation == *evaluation)
    }
}

impl AuthenticationHash {

    pub fn new(hash: Hash, is_first: bool) -> Self {
        Self {
            hash,
            is_first,
        }
    }

    pub fn hash_ref(&self) -> &[u8] {
        self.hash.as_bytes().as_ref()
    }
}

impl<F: PrimeField> Polynomial<F> {
    pub(crate) fn authentication_path_for(&self, root: &F) -> AuthenticationPath<F> {

        let target = self.eval_single(root);

        let log_n = self.log_n() + FRI_BLOWUP_LOG;
        let honest_base_generator = Domain::root_with_order_unchecked((self.len() * FRI_BLOWUP_FACTOR) as u64);
        let mut evaluations = self.coefficients();
        evaluations.append(&mut vec![F::ZERO; (1<<log_n) - (1<<self.log_n())]);

        serial::serial_fft(evaluations.as_mut_slice(), &honest_base_generator, log_n as u32);

        let mut hash_vec = Vec::with_capacity((self.len()/2)-1);
        let mut authentication_vec = Vec::with_capacity(log_n);
        let mut first_evaluation: Option<F> = None;
        let mut second_evaluation: Option<F> = None;
        let mut index = 0;
        for i in (0..evaluations.len()).step_by(2) {
            if ((evaluations[i] == target) || (evaluations[i+1] == target)) && first_evaluation.is_none() {
                first_evaluation = Some(evaluations[i]);
                second_evaluation = Some(evaluations[i+1]);
                index = i/2;
            }
            else {
                let mut hasher = Hasher::new();
                hasher.update(evaluations[i].to_repr().as_ref());
                hasher.update(evaluations[i+1].to_repr().as_ref());
                hash_vec.push(hasher.finalize());
            }
        }
        assert!(first_evaluation.is_some(), "Polynomial does not have root.");

        for _i in 0..(log_n-2) {

            // I should find a better way to do this because vec.remove is O(n)
            if index & 1 == 0 { authentication_vec.push(AuthenticationHash::new(hash_vec.remove(index), false)); }
            else { authentication_vec.push(AuthenticationHash::new(hash_vec.remove(index-1), true)); }

            let mut new_hash_vec = Vec::with_capacity(hash_vec.len()/2);
            for j in (0..hash_vec.len()).step_by(2) {
                let mut hasher = Hasher::new();
                hasher.update(hash_vec[j].as_bytes().as_slice());
                hasher.update(hash_vec[j+1].as_bytes().as_slice());
                new_hash_vec.push(hasher.finalize());
            }
            hash_vec = new_hash_vec;
            index /= 2;
        }

        let is_first = index != 0;

        assert!(hash_vec.len() == 1);
        authentication_vec.push(AuthenticationHash::new(hash_vec[0], is_first));
        AuthenticationPath::new(first_evaluation.unwrap(), second_evaluation.unwrap(), authentication_vec)
    }
}
