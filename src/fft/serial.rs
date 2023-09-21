use ff::PrimeField;

pub(crate) fn serial_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32) {
    
    #[inline(always)]
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow([(n / (2*m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = F::ONE;
            for j in 0..m {

                let mut t = a[(k+j+m) as usize];
                t.mul_assign(&w);
                let mut tmp = a[(k+j) as usize];
                tmp.sub_assign(&t);
                a[(k+j+m) as usize] = tmp;
                a[(k+j) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2*m;
        }

        m *= 2;
    }
}

pub(crate) fn serial_ifft<F: PrimeField>(s: &mut [F], omega: &F, log_n: u32) {
    let invlen = F::from_u128(s.len() as u128).invert().unwrap();
    serial_fft(s, &omega.invert().unwrap(), log_n);
    for item in s.iter_mut() {
        *item *= invlen;
    }
}

