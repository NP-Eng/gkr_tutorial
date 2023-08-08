pub mod data_structures;
#[allow(non_upper_case_globals)]
pub mod gkr;
pub mod oracles;
pub mod parties;
mod tests;
mod utils;
pub use ark_bn254::Fr;
use ark_ff::fields::{Fp64, MontBackend, MontConfig};

