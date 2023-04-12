
mod polynomial;
mod parties;
mod tests;

use ark_ff::{fields::{Fp64, MontBackend, MontConfig}};

use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term}, DenseMVPolynomial,
};

use parties::run_sumcheck_protocol;


#[derive(MontConfig)]
#[modulus = "101"] // MODIFY
#[generator = "2"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

fn main() {

    // MODIFY
    let pol = SparsePolynomial::from_coefficients_vec(
        3,
        vec![
            (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
            (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
        ],
    );

    let verbose = true; // MODIFY

    if run_sumcheck_protocol(pol, verbose) {
        println!("Verifier trusts!");
    } else {
        println!("Verifier aborted");
    }

}
