use ff::PrimeField;

/*
 * I copied this from Matter Labs
 * https://github.com/matter-labs/hodor/blob/master/src/domains/mod.rs
*/

#[derive(Debug)]
pub enum SynthesisError {
    Error,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct Domain<F: PrimeField> {
    pub size: u64,
    pub power_of_two: u64,
    pub generator: F,
}

impl<F: PrimeField> Domain<F> {

    pub(crate) fn root_with_order_unchecked(order: u64) -> F {
        let mut size = order.next_power_of_two();
        let mut log_n = 0;
        while size != 1 {
            size >>= 1;
            log_n += 1;
        }
        let max_power_of_two = F::S;
        let mut g = F::ROOT_OF_UNITY;
        for _ in log_n..(max_power_of_two-1) {
            g = g.square();
        }
        g
    }

    pub fn new_for_size(size: u64) -> Result<Self, SynthesisError> {
        let size = size.next_power_of_two();
        let mut power_of_two = 0;
        let mut k = size;
        while k != 1 {
            k >>= 1;
            power_of_two += 1;
        }
        let max_power_of_two = F::S as u64;
        if power_of_two > max_power_of_two {
            return Err(SynthesisError::Error);
        }

        let mut generator = F::ROOT_OF_UNITY;
        for _ in power_of_two..max_power_of_two {
            generator = generator.square();
        }

        Ok(Self {
            size,
            power_of_two,
            generator,
        })
    }
}
