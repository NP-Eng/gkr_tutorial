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
    // MODIFY
    // pol = 1 - 10 * x1 * x2 + 5 * x1 + 4 * x2 ; hypersum: 12
    // r = 1
    // for each b = {0, 1}: f(x1, b)
    // f(x1, 0) = 1 + 5 x1 ; f(0, 0) = 1; f(1, 0) = 6; f(2, 0) = 11
    // f(x1, 1) = 5 - 5 x1 ; f(0, 1) = 5; f(1, 1) = 0; f(2, 1) = -5
    // interpolate a poly which takes the value 6 at 0 and 6 at 1

    // f(48, x2)
    // for t = 0,1,2
    // f(48, 0) = 39 (mod 101); f(48, 1) = 68 (mod 101)

    // let pol2: DenseMultilinearExtension::from_evaluations_vec(
    //     2,
    //     vec![0, 2, 1, 3]
    //         .into_iter()
    //         .map(|x| Fq::from(x as u64))
    //         .collect(),
    // );

    let pol = DenseMultilinearExtension::from_evaluations_vec(
        2,
        vec![1, 6, 5, 0]
            .into_iter()
            .map(|x| Fq::from(x as u64))
            .collect(),
    );

    // write some assertions for pol evaluations
    //
    assert_eq!(pol.evaluate(&[Fq::from(0), Fq::from(0)]).unwrap(), Fq::from(1));
    assert_eq!(pol.evaluate(&[Fq::from(1), Fq::from(0)]).unwrap(), Fq::from(6));
    assert_eq!(pol.evaluate(&[Fq::from(2), Fq::from(0)]).unwrap(), Fq::from(11));
    assert_eq!(pol.evaluate(&[Fq::from(0), Fq::from(1)]).unwrap(), Fq::from(5));

    run_sumcheck_protocol(pol);
}
