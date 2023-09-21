use ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "18446744069414584321"] // 2^64 - 2^32 + 1
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct Goldilocks([u64; 2]);
