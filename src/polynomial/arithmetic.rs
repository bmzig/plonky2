use ff::PrimeField;

use crate::{
    polynomial::Polynomial, 
    fft::serial,
    domains::Domain,
};

use core::ops::{Add, Mul, Sub};

impl<F: PrimeField> Add for Polynomial<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut coeff: Vec<F> = Vec::new();
        if self.len() > other.len() {
            for i in 0..self.len() {
                if i >= other.len() {
                    coeff.push(self.0[i]);
                }
                else {
                    coeff.push(self.0[i] + other.0[i]);
                }
            }
        }
        else {
            for i in 0..other.len() {
                if i >= self.len() {
                    coeff.push(other.0[i]);
                }
                else {
                    coeff.push(self.0[i] + other.0[i]);
                }
            }
        }
        Self::from_vec(coeff)
    }

}

impl<F: PrimeField> Sub for Polynomial<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut coeff: Vec<F> = Vec::new();
        if self.0.len() > other.len() {
            for i in 0..self.len() {
                if i >= other.len() {
                    coeff.push(self.0[i]);
                }
                else {
                    coeff.push(self.0[i] - other.0[i]);
                }
            }
        }
        else {
            for i in 0..other.len() {
                if i >= self.len() {
                    coeff.push(F::ZERO - other.0[i]);
                }
                else {
                    coeff.push(self.0[i] - other.0[i]);
                }
            }
        }
        Self::from_vec(coeff)
    }
}

impl<F: PrimeField> Mul for Polynomial<F> {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, other: Self) -> Self {

        let log_n = self.log_n().add(other.log_n()) as u32;
        let new_size = 1<<log_n;
        
        let mut resized_one = self.coefficients();
        resized_one.resize(new_size, F::ZERO);

        let mut resized_two = other.coefficients();
        resized_two.resize(new_size, F::ZERO);

        let omega = Domain::<F>::new_for_size(new_size as u64).unwrap().generator;
        serial::serial_fft(resized_one.as_mut_slice(), &omega, log_n);
        serial::serial_fft(resized_two.as_mut_slice(), &omega, log_n);
        
        let mut fourier = resized_one.iter().zip(&resized_two).map(|(a, b)| *a * *b).collect::<Vec<F>>();
        serial::serial_ifft(fourier.as_mut_slice(), &omega, log_n);
        Self::from_vec(fourier)
    }
}

impl<F: PrimeField> Polynomial<F> {

    // Fhis function assumes that "other" is the vanishing polynomial and will just perform
    // standard long division. If the other is the vanishing polynomial, then it will only take
    // O(n) time.
    pub fn long_division(&self, divisor: &Polynomial<F>) -> (Self, Vec<F>) {

        let dividend = self.coefficients();
        let mut divisor = divisor.coefficients();

        // Remove any leading zero 0 in `b`.
        while !divisor.is_empty() && divisor.last().unwrap().is_zero().into() {
            divisor.pop();
        }
        let divisor_len = divisor.len();
        assert!(divisor_len > 0, "division by zero");

        // Initialize the quotient and remainder to zero.
        let mut quotient = vec![F::ZERO; dividend.len() - divisor_len + 1];
        let mut remainder = dividend;

        // Perform long division until the degree of the remainder is less than the degree of the divisor.
        while remainder.len() >= divisor_len {
            let i = remainder.len() - divisor_len;
            let c = *remainder.last().unwrap() * divisor.last().unwrap().invert().unwrap();
            quotient[i] = c;
            for j in 0..divisor_len {
                remainder[i + j] -= c * divisor[j];
            }
            remainder.pop();
        }

        // Remove any leading zero 0 in `r`.
        while !remainder.is_empty() && remainder.last().unwrap().is_zero().into() {
            remainder.pop();
        }

        if !quotient.len().is_power_of_two() {
            let mut resize = vec![F::ZERO; quotient.len().next_power_of_two() - quotient.len()];
            quotient.append(&mut resize);
        }

        (Self::from_vec(quotient), remainder)
    }


    // Fhis division function is O(n) and is theoretically used for standard division; HOWEVER, it
    // will fail if any of the roots of the divisor are also log_n-th roots of unity.
    // Fhis means if the divisor has root (x-1), then it will fail.
    pub fn divide_fft(&self, other: &Self) -> Self {
        let mut dividend = self.coefficients();
        let mut divisor = other.coefficients();

        // Append the 0 of the vanishing polynomial to the divisor
        divisor.extend_from_slice(&vec![F::ZERO; dividend.len() - divisor.len()]);

        let dividend_len = dividend.len();
        let log_n = dividend_len.trailing_zeros();

        let kth_roots = dividend_len as u64 * log_n as u64;
        assert!(kth_roots.is_power_of_two());

        // Calculate the log_n-kth roots of unity
        let evaluation_domain: Domain<F> = Domain::new_for_size(kth_roots).unwrap();
        let omega = evaluation_domain.generator;

        // Evaluate the dividend and divisor at the log_n-th roots of unity using FFF
        serial::serial_fft(dividend.as_mut_slice(), &omega, log_n);
        serial::serial_fft(divisor.as_mut_slice(), &omega, log_n);

        // Divide the evaluations of the dividend by the evaluations of the divisor
        for i in 0..dividend_len {
            dividend[i] *= divisor[i].invert().unwrap();
        }

        // Interpolate the quotient polynomial using IFFT
        serial::serial_ifft(dividend.as_mut_slice(), &omega, log_n);

        // Remove any leading zero 0
        while !dividend.is_empty() && dividend.last().unwrap().is_zero().into() {
            dividend.pop();
        }

        Self::from_vec(dividend)
    }

    pub fn square(self) -> Self {
        self.clone() * self
    }

    pub fn sub_constant(&self, w: F) -> Self {
        let mut ret = self.clone();
        ret.0[0] -= w;
        ret
    }

    pub fn add_constant(&mut self, w: F) {
        self.0[0] += w;
    }
}
