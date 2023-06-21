// mod gkr;
mod parties;
mod polynomial;
// mod tests;
mod utils;

use ark_ff::fields::{Fp64, MontBackend, MontConfig};

use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, DenseMultilinearExtension, MultilinearExtension,
};

use parties::run_sumcheck_protocol;

#[derive(MontConfig)]
#[modulus = "101"] // MODIFY
#[generator = "2"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

fn main() {

    // eval vector in LE indexing
    let pol = DenseMultilinearExtension::from_evaluations_vec(
        2,
        vec![1, 6, 5, 0]
            .into_iter()
            .map(|x| Fq::from(x as u64))
            .collect(),
    );

    // write some assertions for pol evaluations
    assert_eq!(pol.evaluate(&[Fq::from(0), Fq::from(0)]).unwrap(), Fq::from(1));
    assert_eq!(pol.evaluate(&[Fq::from(1), Fq::from(0)]).unwrap(), Fq::from(6));
    assert_eq!(pol.evaluate(&[Fq::from(2), Fq::from(0)]).unwrap(), Fq::from(11));
    assert_eq!(pol.evaluate(&[Fq::from(0), Fq::from(1)]).unwrap(), Fq::from(5));

    run_sumcheck_protocol(pol);
}
