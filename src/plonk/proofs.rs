use ff::PrimeField;

use crate::{
    FriProof,
    plonk::{
        ZeroTestProof, 
        ProductCheckProof, 
        Evaluation, 
        RationalProductCheckProof, 
        PermutationCheckProof, 
        PrescribedPermutationCheckProof
    },
    fri::{FriCommitment, VerificationResult}
};

/*
macro_rules! validate {
    ($($eval:ident),+ $(,)?) => {{
        $(
            assert!($eval.check().is_valid());
        )
    }}
}
*/

impl<F: PrimeField> Evaluation<F> {
    
    pub fn new(eval: F, eval_proof: FriProof<F>) -> Self {
        Self {
            eval,
            eval_proof,
        }
    }

    pub fn check(&self) -> VerificationResult {
        if self.eval_proof.verify().is_valid() { return VerificationResult::ValidProof; }
        VerificationResult::InvalidProof
    }

    pub fn evaluation(&self) -> F {
        self.eval
    }
}

impl<F: PrimeField> ZeroTestProof<F> {

    pub fn new(f_r: Evaluation<F>, q_r: Evaluation<F>) -> Self {
        Self {
            f_r,
            q_r,
        }
    }

    pub fn verify(&self) -> VerificationResult {
        if !self.f_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.q_r.check().is_valid() { return VerificationResult::InvalidProof; }

        let vp = F::ONE;

        if self.f_r.evaluation() != self.q_r.evaluation() * vp { return VerificationResult::InvalidProof; }
        VerificationResult::ValidProof
    }
}

impl<F: PrimeField> ProductCheckProof<F> {

    pub fn new(
        end_eval: Evaluation<F>,
        t_r: Evaluation<F>,
        t_wr: Evaluation<F>,
        q_r: Evaluation<F>,
        f_wr: Evaluation<F>,
    ) -> Self {
        Self {
            end_eval,
            t_r,
            t_wr,
            q_r,
            f_wr,
        }
    }

    pub fn verify(&self) -> VerificationResult {
        if !self.end_eval.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.q_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.f_wr.check().is_valid() { return VerificationResult::InvalidProof; }

        let vp = F::ONE;

        if self.end_eval.evaluation() != F::ONE { return VerificationResult::InvalidProof; }
        let lhs = self.t_wr.evaluation() - (self.t_r.evaluation() * self.f_wr.evaluation());
        let rhs = self.q_r.evaluation() * vp;
        if lhs != rhs { return VerificationResult::InvalidProof; }

        VerificationResult::ValidProof
    }
}

impl<F: PrimeField> RationalProductCheckProof<F> {

    pub fn new(
        end_eval: Evaluation<F>,
        t_r: Evaluation<F>,
        t_wr: Evaluation<F>,
        q_r: Evaluation<F>,
        g_wr: Evaluation<F>,
        f_wr: Evaluation<F>,
    ) -> Self {
        Self {
            end_eval,
            t_r,
            t_wr,
            q_r,
            g_wr,
            f_wr,
        }
    }

    pub fn verify(&self) -> VerificationResult {
        if !self.end_eval.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.q_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.g_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.f_wr.check().is_valid() { return VerificationResult::InvalidProof; }

        let vp = F::ONE;

        if self.end_eval.evaluation() != F::ONE { return VerificationResult::InvalidProof; }
        let lhs = (self.t_wr.evaluation() * self.g_wr.evaluation()) - (self.t_r.evaluation() * self.f_wr.evaluation());
        let rhs = self.q_r.evaluation() * vp;
        if lhs != rhs { return VerificationResult::InvalidProof; }

        VerificationResult::ValidProof
    }
}

impl<F: PrimeField> PermutationCheckProof<F> {

    pub fn new(
        f_com: FriCommitment<F>,
        end_eval: Evaluation<F>,
        t_r: Evaluation<F>,
        t_wr: Evaluation<F>,
        q_r: Evaluation<F>,
        g_wr: Evaluation<F>,
        f_wr: Evaluation<F>,
    ) -> Self {
        Self {
            f_com,
            end_eval,
            t_r,
            t_wr,
            q_r,
            g_wr,
            f_wr,
        }
    }

    pub fn verify(&self) -> VerificationResult {

        if !self.end_eval.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.q_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.g_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.f_wr.check().is_valid() { return VerificationResult::InvalidProof; }

        let vp = F::ONE;
        let r = self.f_com.interpret_as_element();

        // I might need to change the evaluation in the "protocols" file to w^r instead of wr.
        let g = r - self.g_wr.evaluation();
        let f = r - self.f_wr.evaluation();

        if self.end_eval.evaluation() != F::ONE { return VerificationResult::InvalidProof; }
        let lhs = (self.t_wr.evaluation() * g) - (self.t_r.evaluation() * f);
        let rhs = self.q_r.evaluation() * vp;
        if lhs != rhs { return VerificationResult::InvalidProof; }

        VerificationResult::ValidProof
    }
}

impl<F: PrimeField> PrescribedPermutationCheckProof<F> {

    pub fn new(
        f_com: FriCommitment<F>, 
        g_com: FriCommitment<F>,
        end_eval: Evaluation<F>,
        t_r: Evaluation<F>,
        t_wr: Evaluation<F>,
        q_r: Evaluation<F>,
        f_wr: Evaluation<F>,
        g_wr: Evaluation<F>,
        w_wr: Evaluation<F>
    ) -> Self {
        Self {
            f_com,
            g_com,
            end_eval,
            t_r,
            t_wr,
            q_r,
            f_wr,
            g_wr,
            w_wr,
        }
    }

    pub fn verify(&self) -> VerificationResult {

        if !self.end_eval.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.t_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.q_r.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.g_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.f_wr.check().is_valid() { return VerificationResult::InvalidProof; }
        if !self.w_wr.check().is_valid() { return VerificationResult::InvalidProof; }

        let vp = F::ONE;
        let r = self.f_com.interpret_as_element();
        let s = self.g_com.interpret_as_element();
        let a = F::ONE;

        // I might need to change the evaluation in the "protocols" file to w^r instead of wr.
        let f = r - (s * self.w_wr.evaluation()) - self.f_wr.evaluation();
        let g = r - (s * a) - self.g_wr.evaluation();

        if self.end_eval.evaluation() != F::ONE { return VerificationResult::InvalidProof; }
        let lhs = (self.t_wr.evaluation() * g) - (self.t_r.evaluation() * f);
        let rhs = self.q_r.evaluation() * vp;
        if lhs != rhs { return VerificationResult::InvalidProof; }

        VerificationResult::ValidProof
    }
}
