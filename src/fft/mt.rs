use ff::PrimeField;

// TODO: Make multi-threaded FFT.
pub(crate) fn mt_fft<F: PrimeField>(_s: Vec<F>, _g: &F, _threads: u16) -> Vec<F> {
    unimplemented!();
}

pub(crate) fn mt_ifft<F: PrimeField>(_s: Vec<F>, _g: &F, _threads: u16) -> Vec<F> {
    unimplemented!();
    // let invlen = F::from_u128(s.len() as u128).invert().unwrap();
    // mt_fft(s, g, threads)
}

