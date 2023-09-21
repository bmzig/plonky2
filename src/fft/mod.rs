pub(crate) mod serial;

#[cfg(test)]
mod fft_tests {

    use super::*;
    use ff::Field;
    use crate::{
        field::goldilocks::Goldilocks,
        domains::Domain,
        constants::*,
    };

    #[test]
    #[ignore]
    fn benchmark() {
        let log_n = 8;
        let base = 1<<log_n;
        let size = FRI_BLOWUP_FACTOR * base;

        let omega = Domain::new_for_size(size as u64).unwrap().generator;

        let rng = rand::thread_rng();

        let mut a = (0..base).map(|_| Goldilocks::random(rng.clone())).collect::<Vec<_>>();
        let _b = a.clone();
        let now = std::time::Instant::now();
        serial::serial_fft(a.as_mut_slice(), &omega, log_n);
        let after = std::time::Instant::now();
        println!("Serial FFT took {:?}", after - now);
        let now = std::time::Instant::now();
        // recursive::recursive_fft(_b, &omega);
        let after = std::time::Instant::now();
        println!("Recursive FFT took {:?}", after - now);
    }
}
