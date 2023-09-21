use ff::PrimeField;

pub(crate) fn recursive_fft<F: PrimeField>(s: Vec<F>, g: &F) -> Vec<F> {
    let split_vec = | x: Vec<F> | -> (Vec<F>, Vec<F>) {
        let mut ev = Vec::new();
        let mut ov = Vec::new();
        for i in 0..x.len() {
            if i % 2 == 0 { ev.push(x[i]); }
            else { ov.push(x[i]); }
        }
        (ev, ov)
    };

    if s.len() == 1 { return s; }
    let (ov, ev) = split_vec(s);
    let l = recursive_fft(ev, &g.pow([2u64]));
    let r = recursive_fft(ov, &g.pow([2u64]));
    let mut o = Vec::new();
    let mut a = Vec::new();
    for (i, (x,y)) in r.into_iter().zip(l).enumerate() {
        let yr = y*g.pow([i as u64]);
        o.push(x+yr);
        a.push(x-yr);
    }
    o.append(&mut a);
    o
}

pub(crate) fn recursive_ifft<F: PrimeField>(s: Vec<F>, g: &F) -> Vec<F> {
    let invlen = F::from_u128(s.len() as u128).invert().unwrap();
    recursive_fft(s, &g.invert().unwrap()).into_iter().map(|x| x * invlen).collect::<Vec<F>>()
}
