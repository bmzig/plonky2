pub(crate) mod goldilocks;

use ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "3618502788666131213697322783095070105623107215331596699973092056135872020481"] // 2^251 + 2^196 + 2^192 + 1
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fp([u64; 4]);
