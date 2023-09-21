use ff::PrimeField;

use crate::{
    plonk::{Circuit, Gate},
};

impl<F: PrimeField> Circuit<F> {

    pub fn new(gates: Vec<Gate<F>>) -> Self {
        unimplemented!();
    }

}
