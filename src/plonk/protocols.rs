use ff::PrimeField;

use crate::{
    FriProof,
    plonk::{
        ZeroTestProof, 
        ProductCheckProof, 
        Evaluation, 
        RationalProductCheckProof,
        PermutationCheckProof,
        PrescribedPermutationCheckProof,
    },
    polynomial::Polynomial,
    domains::Domain,
    fft::serial,
};

impl<F: PrimeField> Polynomial<F> {

    // <------------------------------------------------------------------------------------------->
    // For the zero test, a prover wants to convince the verifier that they have knowledge of some
    // polynomial p(x) which is zero for some subset omega in F_p (in practice, the subset is some 
    // set of primitive roots of unity. If this is true, then the prover can divide out the polynomial 
    // constructed from those roots. The prover then commits to q(x) = p(x)/Z(x) (using KZG, FRI,
    // Bulletproofs, etc.). Verifier samples a random r (or use of Fiat-Shamir for noninteractive), and
    // prover queries q(r) and p(r) and sends field elements to prover. KZG is just checking q(r)Z(r)
    // = p(r), FRI checks Merkle authentication paths, etc.
    // <------------------------------------------------------------------------------------------->
    pub fn zero_test(&self, z: &Polynomial<F>) -> ZeroTestProof<F> {
        let (q_x, _) = self.long_division(z);
        
        let r: F = q_x
            .commitment()
            .interpret_as_element();

        let f_eval = self.eval_single(&r);
        let q_eval = q_x.eval_single(&r);

        let f_eval_proof = FriProof::evaluation_proof(self, Some(r));
        let q_eval_proof = FriProof::evaluation_proof(&q_x, Some(r));

        ZeroTestProof::new(
            Evaluation::new(f_eval, f_eval_proof),
            Evaluation::new(q_eval, q_eval_proof)
        )

    }

    // <------------------------------------------------------------------------------------------->
    // For the product check, the prover wants to convince the verifier that they have knowledge of
    // some polynomial k-degree polynomial f(x) such that the series product of f(x) from 0 to k is
    // 1. It is assumed that the verifier has a commitment to f. To begin, the prover interpolates 
    // t(x) from points of the k-degree polynomial of knowledge f(x) such that t(1) = f(1), t(w) = 
    // f(1)f(w), t(w^2) = f(1)f(w)f(w^2), ... , t(w^k-1) = f(1)f(w)...f(w^k-1). There is an important 
    // recurrence relation here: For all points on the set of omega, t(wx) = t(x) * f(wx) (even true 
    // for t(w^k-1) since w is a primitive root of unity).
    //
    // The important things to notice here is that 1. t(w^k-1) should be 1, and 2. t(wx) - t(x) * f(wx)
    // should be 0. The round begins with the prover constructing another polynomial h(x) = t(w*x)
    // - t(x) * f(w*x). (From here it is similar to the zero test). The prover creates a q(x)
    // = h(x)/Z(x), and commits to q(x) and t(x). The verifier samples a random r in F_p and queries
    // t(x) at w^k-1, r, and w, and q(x) at r and f(x) at wr. The verifier will accept if and only if
    // t(w^k-1) == 1, t(wr) - t(r)f(wr) == q(r)(r^k - 1), and all the commitment checks are valid.
    // <------------------------------------------------------------------------------------------->
    pub fn product_check(&self) -> ProductCheckProof<F> {

        let mut evaluations = self.coefficients();

        let size = self.len().next_power_of_two();
        let log_n = {
            let mut x = 0u32;
            let mut y = size;
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            x
        };

        let omega = Domain::root_with_order_unchecked(size as u64);
        serial::serial_fft(evaluations.as_mut_slice(), &omega, log_n);

        let mut t_x = vec![F::ZERO; size];
        let mut target = F::ONE;
        for i in 0..size {
            target *= evaluations[i];
            t_x[i] = target;
        }

        let t_end = t_x.last().unwrap().clone();

        serial::serial_ifft(t_x.as_mut_slice(), &omega, log_n);
        let t_x = Polynomial::from_vec(t_x);
        let t_end_proof = FriProof::evaluation_proof(&t_x, Some(omega.pow([size as u64 - 1])));

        let r = t_x.commitment().interpret_as_element();

        let t_r = t_x.eval_single(&r);
        let t_r_proof = FriProof::evaluation_proof(&t_x, Some(r));
        let t_wr = t_x.eval_single(&(omega * r));
        let t_wr_proof = FriProof::evaluation_proof(&t_x, Some(omega * r));

        let vanishing = Polynomial::vanishing_polynomial(size as u128);
        let (q_x, _) = t_x.long_division(&vanishing);
        let q_r = q_x.eval_single(&r);
        let q_r_proof = FriProof::evaluation_proof(&q_x, Some(r));

        let f_wr = self.eval_single(&(omega * r));
        let f_wr_proof = FriProof::evaluation_proof(self, Some(omega * r));

        ProductCheckProof::new(
            Evaluation::new(t_end, t_end_proof),
            Evaluation::new(t_r, t_r_proof),
            Evaluation::new(t_wr, t_wr_proof),
            Evaluation::new(q_r, q_r_proof),
            Evaluation::new(f_wr, f_wr_proof)
        )
    }

    // <------------------------------------------------------------------------------------------->
    // The product check for rational functions is largely the same. The prover wants to show that they
    // know polynomials f(x) and q(x) such that series product from 0 to k of f(w^i)/g(w^i) = 1. The
    // verifier first takes commitments to f(x) and g(x) through KZG, FRI, etc. Like before, the prover
    // interpolates a polynomial t(x) such that t(1) = f(1)/g(1), t(w) = f(w)/g(w), t(w^2)
    // = f(w^2)/g(w^2), ... , t(w^k-1) = f(w^k-1)/g(w^k-1). Like last time, if constructed honestly,
    // then t(w^k-1) = 1 and t(wx) * g(wx) = t(x) * f(wx) for all x in the subset omega.
    // <------------------------------------------------------------------------------------------->
    pub fn product_check_rational(&self, denominator: &Polynomial<F>) -> RationalProductCheckProof<F> {

        let mut numerator_evaluations = self.coefficients();
        let mut denominator_evaluations = denominator.coefficients();

        let size = self.len().next_power_of_two();
        let log_n = {
            let mut x = 0u32;
            let mut y = size;
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            x
        };

        let omega = Domain::root_with_order_unchecked(size as u64);
        serial::serial_fft(numerator_evaluations.as_mut_slice(), &omega, log_n);
        serial::serial_fft(denominator_evaluations.as_mut_slice(), &omega, log_n);

        let mut t_x = vec![F::ZERO; size];
        let mut target = F::ONE;
        for i in 0..size {
            target *= numerator_evaluations[i] * denominator_evaluations[i].invert().unwrap();
            t_x[i] = target;
        }

        let t_end = t_x.last().unwrap().clone();

        serial::serial_ifft(t_x.as_mut_slice(), &omega, log_n);
        let t_x = Polynomial::from_vec(t_x);
        let t_end_proof = FriProof::evaluation_proof(&t_x, Some(omega.pow([size as u64 - 1])));

        let r = t_x.commitment().interpret_as_element();

        let t_r = t_x.eval_single(&r);
        let t_r_proof = FriProof::evaluation_proof(&t_x, Some(r));
        let t_wr = t_x.eval_single(&(omega * r));
        let t_wr_proof = FriProof::evaluation_proof(&t_x, Some(omega * r));

        let vanishing = Polynomial::vanishing_polynomial(size as u128);
        let (q_x, _) = t_x.long_division(&vanishing);
        let q_r = q_x.eval_single(&r);
        let q_r_proof = FriProof::evaluation_proof(&q_x, Some(r));

        let g_wr = denominator.eval_single(&(omega * r));
        let g_wr_proof = FriProof::evaluation_proof(denominator, Some(omega * r));
        let f_wr = self.eval_single(&(omega * r));
        let f_wr_proof = FriProof::evaluation_proof(self, Some(omega * r));

        RationalProductCheckProof::new(
            Evaluation::new(t_end, t_end_proof),
            Evaluation::new(t_r, t_r_proof),
            Evaluation::new(t_wr, t_wr_proof),
            Evaluation::new(q_r, q_r_proof),
            Evaluation::new(g_wr, g_wr_proof),
            Evaluation::new(f_wr, f_wr_proof)
        )

    }

    // <------------------------------------------------------------------------------------------->
    // In this protocol, the prover wants to show that for a subset omega in F_p (for example, the set
    // of the nth primitive roots of unity), the evaluations of f(omega) is a permutation of the
    // evaluations of g(omega). We assume that the prover commited to f and g through KZG, FRI, etc.
    //
    // For the protocol, we contruct f_hat(x) and g_hat(x). Both are constructed by interpolating off
    // of the roots of (x-f(a)) and (x-g(a)) for a being an element in omega (so that we can use FFT).
    // As you can see, both f_hat(x) and g_hat(x) should be equal to each other if they are indeed
    // permutations of each other. Now, the prover and the verifier can engage in the product check
    // protocol and prove that f_hat(x)/g_hat(x) = 1 for all x in omega.
    // <------------------------------------------------------------------------------------------->
    pub fn permutation_check(&self, permutation: &Polynomial<F>) -> PermutationCheckProof<F> {

        let f_commitment = self.commitment();
        let r = f_commitment.interpret_as_element();

        let mut f_evals = self.coefficients();
        let mut g_evals = permutation.coefficients();

        let size = self.len().next_power_of_two();
        let log_n = {
            let mut x = 0u32;
            let mut y = size;
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            x
        };

        let omega = Domain::root_with_order_unchecked(size as u64);
        serial::serial_fft(f_evals.as_mut_slice(), &omega, log_n);
        serial::serial_fft(g_evals.as_mut_slice(), &omega, log_n);

        let mut t_x = vec![F::ZERO; size];
        let mut target = F::ONE;
        for i in 0..size {
            target *= (r - f_evals[i]) * (r - g_evals[i]).invert().unwrap();
            t_x[i] = target;
        }

        let t_end = t_x.last().unwrap().clone();

        serial::serial_ifft(t_x.as_mut_slice(), &omega, log_n);
        let t_x = Polynomial::from_vec(t_x);
        
        let t_end_proof = FriProof::evaluation_proof(&t_x, Some(omega.pow([size as u64 - 1])));

        let t_r = t_x.eval_single(&r);
        let t_r_proof = FriProof::evaluation_proof(&t_x, Some(r));
        let t_wr = t_x.eval_single(&(omega * r));
        let t_wr_proof = FriProof::evaluation_proof(&t_x, Some(omega * r));

        let vanishing = Polynomial::vanishing_polynomial(size as u128);
        let (q_x, _) = t_x.long_division(&vanishing);
        let q_r = q_x.eval_single(&r);
        let q_r_proof = FriProof::evaluation_proof(&q_x, Some(r));

        let g_wr = permutation.eval_single(&(r * omega));
        let g_wr_proof = FriProof::evaluation_proof(permutation, Some(omega * r));

        let f_wr = self.eval_single(&(r * omega));
        let f_wr_proof = FriProof::evaluation_proof(self, Some(omega * r));

        PermutationCheckProof::new(
            f_commitment,
            Evaluation::new(t_end, t_end_proof),
            Evaluation::new(t_r, t_r_proof),
            Evaluation::new(t_wr, t_wr_proof),
            Evaluation::new(q_r, q_r_proof),
            Evaluation::new(g_wr, g_wr_proof),
            Evaluation::new(f_wr, f_wr_proof)
        )
        
    }

    // <------------------------------------------------------------------------------------------->
    // In this protocol, given a known permutation W: w_1 -> w_2, the prover wants to convince the
    // verifier that they have polynomials f(x) and g(x) which satisfy f(y) = g(W(y)) for all y in
    // omega. Essentially, this is saying that not only are f(x) and g(x) permutations of each other,
    // but the two polynomials are linked together through a specific permutation W. We assume that the
    // verifier has commitments to f, g, and W through methods stated before.
    //
    // We need to do this in linear time, which is achieved by using a clever product check.
    // Essentially, if f(x) is a prescribed permutation of g(x), then for an element a in omega, 
    // (W(a), f(a)) will be a permutation of (a, g(a)).
    //
    // For example, let W(w^0) = w^2, W(w^1) = w^0, and W(w^2) = w^1, then...
    // (w^0, g(w^0)), (w^1, g(w^1)), (w^2, g(w^2)) -> i.e. the values in the right tuple
    // (w^2, f(w^0)), (w^0, f(w^1)), (w^1, f(w^2)) -> i.e. the values in the left tuple
    // will be permutations of each other. 
    //
    // We can think of constructing two polynomials, f_hat(x, y) and g_hat(x, y), where f_hat(x,y) is
    // equal to (x - y * W(a) - f(a)) and g_hat(x, y) is (x - y * a - g(a)), but this would complicate
    // the algorithm, so we do not actually do this. Instead, we perform a product check on the two
    // polynomials. The verifier samples two random values in the field, s, and r, and the verifier and
    // prover go through the product check, except the prover calculates the top polynomial (i.e.
    // f_hat(x, y)) as (r - s * W(x) - f(x)) and the bottom polynomial as (r - s * x - g(x)). Both
    // these polynomials are univariate, and if the product check passes, then f(x) is a prescribed
    // permutation of g(x) with high probability, since both polynomials are equal at this random
    // point.
    // <------------------------------------------------------------------------------------------->
    pub fn prescribed_permutation_check(&self, permutation: &Polynomial<F>, rules: &Polynomial<F>) -> PrescribedPermutationCheckProof<F> {

        let f_commitment = self.commitment();
        let g_commitment = permutation.commitment();
        let r = f_commitment.interpret_as_element();
        let s = g_commitment.interpret_as_element();

        let mut f_evals = self.coefficients();
        let mut g_evals = permutation.coefficients();
        let mut rules_evals = rules.coefficients();

        let size = self.len().next_power_of_two();
        let log_n = {
            let mut x = 0u32;
            let mut y = size;
            while y != 1 {
                y >>= 1;
                x += 1;
            }
            x
        };

        let omega = Domain::root_with_order_unchecked(size as u64);
        serial::serial_fft(f_evals.as_mut_slice(), &omega, log_n);
        serial::serial_fft(g_evals.as_mut_slice(), &omega, log_n);
        serial::serial_fft(rules_evals.as_mut_slice(), &omega, log_n);

        let mut t_x = vec![F::ZERO; size];
        let mut target = F::ONE;
        let mut g = F::ONE;
        for i in 0..size {
            target *= (r - (s * rules_evals[i]) - f_evals[i]) * (r - (s * g) - g_evals[i]).invert().unwrap();
            g *= omega;
            t_x[i] = target;
        }

        let t_end = t_x.last().unwrap().clone();

        serial::serial_ifft(t_x.as_mut_slice(), &omega, log_n);
        let t_x = Polynomial::from_vec(t_x);
        
        let t_end_proof = FriProof::evaluation_proof(&t_x, Some(omega.pow([size as u64 - 1])));

        let t_r = t_x.eval_single(&r);
        let t_r_proof = FriProof::evaluation_proof(&t_x, Some(r));
        let t_wr = t_x.eval_single(&(omega * r));
        let t_wr_proof = FriProof::evaluation_proof(&t_x, Some(omega * r));

        let vanishing = Polynomial::vanishing_polynomial(size as u128);
        let (q_x, _) = t_x.long_division(&vanishing);
        let q_r = q_x.eval_single(&r);
        let q_r_proof = FriProof::evaluation_proof(&q_x, Some(r));

        let g_wr = permutation.eval_single(&(r * omega));
        let g_wr_proof = FriProof::evaluation_proof(permutation, Some(omega * r));

        let f_wr = self.eval_single(&(r * omega));
        let f_wr_proof = FriProof::evaluation_proof(self, Some(omega * r));

        let w_wr = rules.eval_single(&(r * omega));
        let w_wr_proof = FriProof::evaluation_proof(rules, Some(omega * r));

        PrescribedPermutationCheckProof::new(
            f_commitment,
            g_commitment,
            Evaluation::new(t_end, t_end_proof),
            Evaluation::new(t_r, t_r_proof),
            Evaluation::new(t_wr, t_wr_proof),
            Evaluation::new(q_r, q_r_proof),
            Evaluation::new(g_wr, g_wr_proof),
            Evaluation::new(f_wr, f_wr_proof),
            Evaluation::new(w_wr, w_wr_proof)
        )
    }
}
