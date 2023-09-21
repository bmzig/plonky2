use ff::PrimeField;
use primitive_types::U256;

pub(crate) fn field_element_from_bytes<F: PrimeField>(bytes: &[u8]) -> F {
    let mut repr = F::Repr::default();

    let repr_size = repr.as_ref().len();
    let mut parsed_bytes: [u8; 32] = [0u8; 32];
    let modulus = U256::from_str_radix(F::MODULUS, 16).unwrap();
    (U256::from_big_endian(bytes) % modulus).to_little_endian(parsed_bytes.as_mut_slice());
    let copy = parsed_bytes.chunks_exact(repr_size).next().expect("Repr is larger than 256 bits.");

    repr.as_mut()[..repr_size].copy_from_slice(&copy[..repr_size]);
    
    F::from_repr(repr).unwrap()
}

#[cfg(test)]
mod utils_tests {
    use super::*;
    use crate::field::Fp;

    use ff::Field;

    #[test]
    fn from_uniform_bytes() {
        let mut one: [u8; 32] = [0; 32];
        one[31] = 0x01;
        let element_one: Fp = field_element_from_bytes(one.as_slice());
        assert_eq!(element_one, Fp::ONE);
        
        let max: [u8; 32] = [
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x11, 
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
        ];
        let element_max: Fp = field_element_from_bytes(max.as_slice());
        assert_eq!(element_max, Fp::ZERO-Fp::ONE);

        let wraparound: [u8; 32] = [
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x11, 
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 
        ];
        let element_wraparound: Fp = field_element_from_bytes(wraparound.as_slice());
        assert_eq!(element_wraparound, Fp::from(14)); // f - 1 = e
    }
}
