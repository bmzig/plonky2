use ff::PrimeField;

use crate::{
    fri::FriCommitment,
    stark::{FriChallenge},
};

mod fft;
mod fri;
mod polynomial;
mod utils;
mod field;
mod domains;
mod plonk;
mod constants;
mod stark;

/*
#[derive(Debug)]
pub struct PlonkyProof {
}
*/

#[derive(Debug)]
pub struct FriProof<F: PrimeField> {
    w_com: FriCommitment<F>,
    fri_challenge: FriChallenge<F>
}

impl<F: PrimeField> FriProof<F> {

    pub(crate) fn w_com(&self) -> &FriCommitment<F> {
        &self.w_com
    }

    pub(crate) fn fri_challenge(&self) -> &FriChallenge<F> {
        &self.fri_challenge
    }

}

#[cfg(test)]
mod plonky2 {

    use super::*;

    use crate::{
        FriProof,
        polynomial::Polynomial,
        field::goldilocks::Goldilocks,
        stark::FriChallenge,
        domains::Domain,
        constants::*,
    };

    use rand::{Rng, RngCore};

    use ff::{Field};

    fn random_root_of_unity<R: RngCore, F: PrimeField>(rng: &mut R, size: u64) -> F {
        let random_exponent = rng.gen::<u64>();
        let root: F = Domain::root_with_order_unchecked(size);
        root.pow([random_exponent])
    }

    #[test]
    fn fri_semi_interactive_proof() {
        // Setup: prover has secret "knowledge" polynomial f(x). Verifier knows that if prover is
        // honest, then f(x) must have degree <= 4 and must evaluate to y for some input w. 
        let f_x = Polynomial::from_vec(vec![Goldilocks::ONE, Goldilocks::from(5), Goldilocks::from(5), Goldilocks::ONE, Goldilocks::from(10), Goldilocks::from(9), Goldilocks::ZERO, Goldilocks::from(88)]);
        let mut rng = rand::thread_rng();
        let target_degree: u64 = f_x.len().next_power_of_two() as u64;
        let _log_n = target_degree.next_power_of_two().ilog2() as u32;
        
        // (P) Prover sends commitment to f(x) to verifier
        let _f_commitment = f_x.commitment();

        // (V) Verifier samples a random r in the field and gives this value to the prover.
        let random_r = Goldilocks::random(rng.clone());

        // (P) Prover makes polynomial w(x). Gives commitment to verifier.
        let w_x = f_x.shift_polynomial(random_r);
        let w_commitment = w_x.commitment();
        
        // (P) Prover folds w_x. Prover claims to make log_n folds before ending at a constant
        // function. Verifier must now check that prover 1) made log_n folds, and 2) that all folds
        // were done honestly, which is done by checking each commitment on every layer of the
        // fold. This part is noninteractive. The folding is done using randomness derived from
        // commitments of every intermediate polynomial. Commitments are given to the verifier.
        let (commitment_vector, polynomial_vector) = w_x.fold_full();

        // (V) In the first portion of the proof, the verifier queries a random root of unity...
        let random_root_of_unity: Goldilocks = random_root_of_unity(&mut rng, target_degree * FRI_BLOWUP_FACTOR as u64);

        // ... (P) asks the prover for an authentication path for a random root of unity and its
        // negative counterpart for w_x.
        let positive_authentication_path = w_x.authentication_path_for(&random_root_of_unity);
        let negative_authentication_path = w_x.authentication_path_for(&-random_root_of_unity);
        let positive_evaluation = w_x.eval_single(&random_root_of_unity);
        let negative_evaluation = w_x.eval_single(&-random_root_of_unity);
        
        // ... (V) checks whether this authentication path is consistent with the commitment...
        assert_eq!(w_commitment.value(), positive_authentication_path.derive_root());
        assert_eq!(w_commitment.value(), negative_authentication_path.derive_root());

        // ... (P) requests the rest of the required queries for the fold.
        let mut auth_vec = Vec::with_capacity(polynomial_vector.len());
        let mut query_vec = Vec::with_capacity(polynomial_vector.len());
        let mut target = random_root_of_unity.square();
        for i in 0..polynomial_vector.len()-1 {
            auth_vec.push(polynomial_vector[i].authentication_path_for(&-target));
            query_vec.push(polynomial_vector[i].eval_single(&-target));
            target = target.square();
        }
        
        // auth_vec.len() == query_vec.len() and commitment_vec.len() == polynomial_vec.len()
        assert!(auth_vec.len() == query_vec.len());
        assert!(commitment_vector.len() == polynomial_vector.len());
        
        let fri_challenge = FriChallenge::new(
            positive_evaluation, 
            negative_evaluation, 
            positive_authentication_path.clone(), 
            negative_authentication_path.clone(), 
            auth_vec.clone(), 
            query_vec, 
            commitment_vector.clone()
        );
       
        // (V) Now, the verifier has everything needed to check the proof. The verifier first checks
        // the folds with the values received ...
        let should_be_constant_function = fri_challenge.query_check(&w_commitment, &random_root_of_unity);

        assert!(should_be_constant_function == commitment_vector.last().unwrap().value());

        // ... (V) and then verifies that the commitments are consistent with the values obtained
        // from the query phase above.
        let _target = random_root_of_unity.square();

        assert!(positive_authentication_path.contains_evaluation(&positive_evaluation));
        assert!(negative_authentication_path.contains_evaluation(&negative_evaluation));

        for i in 0..auth_vec.len() {
            assert!(auth_vec[i].contains_evaluation(&fri_challenge.fold_queries()[i]));
            assert_eq!(commitment_vector[i].value(), auth_vec[i].derive_root());
        }

        // If everything passes, then the prover either found a collision in the hash function,
        // successfully grinded Fiat Shamir, or overwhelmingly most likely, I made a mistake in the code.
    }

    #[test]
    fn fri_noninteractive_proof() {

        // This time, the prover has a public polynomial f(x) and wants to convince somebody that
        // f(w) = y for some w. The verifier wants to check this, but does not want to run f(w)
        // themselves. Therefore, the prover just makes a single STARK proof for this polynomial.
        let f_x = Polynomial::from_vec(vec![Goldilocks::ONE, Goldilocks::from(5000), Goldilocks::from(50), Goldilocks::ONE, -Goldilocks::from(10), -Goldilocks::from(9), Goldilocks::ZERO, Goldilocks::from(88)]);

        let stark_proof = FriProof::evaluation_proof(&f_x, None);

        let result = stark_proof.verify();

        assert!(result.is_valid())
    }
}
