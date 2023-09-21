use crate::{
    fri::{FriCommitment, AuthenticationPath},
    polynomial::Polynomial,
};

use ff::PrimeField;


mod fold;
mod proof;

#[derive(Debug)]
pub struct FriChallenge<F: PrimeField> {
    positive_evaluation: F,
    negative_evaluation: F,
    positive_authentication_path: AuthenticationPath<F>,
    negative_authentication_path: AuthenticationPath<F>,
    authentication_vector: Vec<AuthenticationPath<F>>,
    fold_queries: Vec<F>,
    commitment_vector: Vec<FriCommitment<F>>,
}

#[derive(Debug, Eq, PartialEq)]
pub enum VerificationResult {
    ValidProof,
    InvalidProof,
}

impl<F: PrimeField> FriChallenge<F> {
    pub fn new(
            positive_evaluation: F, 
            negative_evaluation: F, 
            positive_authentication_path: AuthenticationPath<F>,
            negative_authentication_path: AuthenticationPath<F>,
            authentication_vector: Vec<AuthenticationPath<F>>,
            fold_queries: Vec<F>, 
            commitment_vector: Vec<FriCommitment<F>>
        ) -> Self {
        
            Self {
                positive_evaluation,
                negative_evaluation,
                positive_authentication_path,
                negative_authentication_path,
                authentication_vector,
                fold_queries,
                commitment_vector,
            }
    }

    pub fn positive_evaluation(&self) -> F {
        self.positive_evaluation
    }

    pub fn negative_evaluation(&self) -> F {
        self.negative_evaluation
    }

    pub fn fold_queries(&self) -> &Vec<F> {
        &self.fold_queries
    }

    pub fn commitment_vector(&self) -> &Vec<FriCommitment<F>> {
        &self.commitment_vector
    }

    pub(crate) fn positive_authentication_path(&self) -> &AuthenticationPath<F> {
        &self.positive_authentication_path
    }

    pub(crate) fn negative_authentication_path(&self) -> &AuthenticationPath<F> {
        &self.negative_authentication_path
    }

    pub(crate) fn authentication_paths(&self) -> &Vec<AuthenticationPath<F>> {
        &self.authentication_vector
    }

}

impl<F: PrimeField> Polynomial<F> {

    // The polynomial w(x) s.t. w(x)=(f(x)-v)(x-r)^-1
    pub fn shift_polynomial(&self, r: F) -> Self {
        let eval = self.eval_single(&r);
        let numerator = self.sub_constant(eval);
        let denominator = Polynomial::from_vec(vec![F::ZERO-r, F::ONE]);
        let (w_x, _) = numerator.long_division(&denominator);
        w_x
    }
}

impl VerificationResult {
    pub fn is_valid(&self) -> bool {
        *self == VerificationResult::ValidProof
    }
}

