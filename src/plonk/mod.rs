use ff::PrimeField;

use crate::{
    FriProof,
    polynomial::Polynomial,
    fri::FriCommitment,
};

mod protocols;
mod circuit;
mod proofs;

pub struct Evaluation<F: PrimeField> {
    eval: F,
    eval_proof: FriProof<F>,
}

pub struct Circuit<F: PrimeField> {
    selector: Polynomial<F>,
    gates: Vec<Gate<F>>,
}

pub struct Gate<F: PrimeField> {
    wires: (F, F),
}

pub struct ZeroTestProof<F: PrimeField> {
    f_r: Evaluation<F>,
    q_r: Evaluation<F>,
}

pub struct ProductCheckProof<F: PrimeField> {
    end_eval: Evaluation<F>, // Com(t)
    t_r: Evaluation<F>, // Com(t)
    t_wr: Evaluation<F>, // Com(t)
    q_r: Evaluation<F>, // Com(q)
    f_wr: Evaluation<F>, // Com(f)
}

pub struct RationalProductCheckProof<F: PrimeField> {
    end_eval: Evaluation<F>, // Com(t)
    t_r: Evaluation<F>,
    t_wr: Evaluation<F>,
    q_r: Evaluation<F>, // Com(q)
    g_wr: Evaluation<F>, // Com(g)
    f_wr: Evaluation<F>, // Com(f)
}

pub struct PermutationCheckProof<F: PrimeField> {
    f_com: FriCommitment<F>,
    end_eval: Evaluation<F>, // Com(t)
    t_r: Evaluation<F>,
    t_wr: Evaluation<F>,
    q_r: Evaluation<F>, // Com(q)
    g_wr: Evaluation<F>, // Com(g)
    f_wr: Evaluation<F>, // Com(f)
}

pub struct PrescribedPermutationCheckProof<F: PrimeField> {
    f_com: FriCommitment<F>,
    g_com: FriCommitment<F>,
    end_eval: Evaluation<F>, // Com(t)
    t_r: Evaluation<F>,
    t_wr: Evaluation<F>,
    q_r: Evaluation<F>, // Com(q)
    g_wr: Evaluation<F>, // Com(g)
    f_wr: Evaluation<F>, // Com(f)
    w_wr: Evaluation<F>, // Com(w)
}
