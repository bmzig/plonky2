use blake3::Hash;
use ff::PrimeField;

use std::marker::PhantomData;

mod commitment;
mod authentication;

#[derive(Debug, Clone)]
pub struct FriCommitment<F: PrimeField>(Hash, PhantomData<F>);

#[derive(Debug, Clone)]
pub struct AuthenticationHash {
    hash: Hash,
    is_first: bool,
}

#[derive(Debug, Clone)]
pub struct AuthenticationPath<F: PrimeField> {
    first_evaluation: F,
    second_evaluation: F,
    authentication_path: Vec<AuthenticationHash>
}

impl<F: PrimeField> FriCommitment<F> {
    pub fn new(h: Hash) -> Self {
        Self(h, PhantomData)
    }

    pub fn value(&self) -> Hash {
        self.0
    }

    pub fn next_value(&self) -> Hash {
        blake3::hash(self.0.as_bytes().as_slice())
    }
}
