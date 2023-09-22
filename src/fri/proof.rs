use ff::PrimeField;

use crate::{
    FriProof,
    fri::{FriCommitment, FriChallenge, VerificationResult},
    constants::*,
    polynomial::Polynomial,
};

impl<F: PrimeField> FriProof<F> {

    pub fn new(w_com: FriCommitment<F>, fri_challenge: FriChallenge<F>) -> Self {
        Self {
            w_com,
            fri_challenge,
        }
    }

    pub fn evaluation_proof(f_x: &Polynomial<F>, r: Option<F>) -> Self {
        
        // Prover makes w_x out of f_x and randomness.
        let w_x = {
            if let Some(x) = r { f_x.shift_polynomial(x) }
            else { f_x.shift_polynomial(f_x.commitment().interpret_as_element()) }
        };
        let w_commitment = w_x.commitment();

        // Prover folds w_x. Has vector with intermediate polynomials to use as a utility and a
        // commitment vector to send to the prover.
        let (commitment_vector, polynomial_vector) = w_x.fold_full();

        // Prover commits to the entire fold by interpreting the constant function commitment as
        // a root of unity.
        let random_root_of_unity = commitment_vector
            .last()
            .expect("Commitment vector empty.")
            .interpret_as_root_of_unity((f_x.len().next_power_of_two() * FRI_BLOWUP_FACTOR) as u64);

        // Prover makes authentication paths for w_x and evaluates values accordingly.
        let positive_authentication_path = w_x.authentication_path_for(&random_root_of_unity);
        let negative_authentication_path = w_x.authentication_path_for(&-random_root_of_unity);
        let positive_evaluation = w_x.eval_single(&random_root_of_unity);
        let negative_evaluation = w_x.eval_single(&-random_root_of_unity);

        // Prover makes queries and sources authentication paths for the rest of the fold
        let mut auth_vec = Vec::with_capacity(polynomial_vector.len());
        let mut query_vec = Vec::with_capacity(polynomial_vector.len());

        let mut target = random_root_of_unity.square();
        for polynomial in polynomial_vector.iter().take(polynomial_vector.len()-1) {
            auth_vec.push(polynomial.authentication_path_for(&-target));
            query_vec.push(polynomial.eval_single(&-target));
            target = target.square();
        }

        let fri_challenge = FriChallenge::new(
            positive_evaluation, 
            negative_evaluation, 
            positive_authentication_path, 
            negative_authentication_path, 
            auth_vec, 
            query_vec, 
            commitment_vector
        );

        Self::new(w_commitment, fri_challenge)
    }
    
    pub fn verify(&self) -> VerificationResult {

        // Check that the queries are consistent with the authentication paths
        if !self.fri_challenge().positive_authentication_path().contains_evaluation(&self.fri_challenge().positive_evaluation()) { return VerificationResult::InvalidProof; }
        if !self.fri_challenge().negative_authentication_path().contains_evaluation(&self.fri_challenge().negative_evaluation()) { return VerificationResult::InvalidProof; }
        for i in 0..self.fri_challenge().authentication_paths().len() {
            if !self.fri_challenge().authentication_paths()[i].contains_evaluation(&self.fri_challenge().fold_queries()[i]) { return VerificationResult::InvalidProof; }
        }

        // Check that the authentication paths are consistent with the commitments
        
        if self.w_com().value() != self.fri_challenge().positive_authentication_path().derive_root() { return VerificationResult::InvalidProof; }
        if self.w_com().value() != self.fri_challenge().negative_authentication_path().derive_root() { return VerificationResult::InvalidProof; }

        for i in 0..self.fri_challenge().authentication_paths().len() {
            if self.fri_challenge().authentication_paths()[i].derive_root() != self.fri_challenge().commitment_vector()[i].value() { return VerificationResult::InvalidProof; }
        }
        // Check that the fold is proper
        
        let should_be_root = self
            .fri_challenge()
            .commitment_vector()
            .last()
            .unwrap()
            .interpret_as_root_of_unity(1<<(self.fri_challenge().commitment_vector().len() + FRI_BLOWUP_LOG));

        let should_be_constant_function = self.fri_challenge().query_check(self.w_com(), &should_be_root);
        if should_be_constant_function != self.fri_challenge().commitment_vector().last().unwrap().value() { return VerificationResult::InvalidProof; }
        
        VerificationResult::ValidProof

    }
}
