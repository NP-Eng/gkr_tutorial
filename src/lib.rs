pub mod data_structures;
#[allow(non_upper_case_globals)]
pub mod gkr;
pub mod oracles;
pub mod parties;
mod tests;
mod utils;

use ark_ff::fields::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "101"] // MODIFY
#[generator = "2"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;
