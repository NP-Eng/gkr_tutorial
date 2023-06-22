// mod gkr;
mod parties;
// mod tests;
mod utils;

use ark_ff::fields::{Fp64, MontBackend, MontConfig};

use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, DenseMultilinearExtension, MultilinearExtension,
};

use parties::run_sumcheck_protocol;

use crate::utils::interpolate_uni_poly;

#[derive(MontConfig)]
#[modulus = "101"] // MODIFY
#[generator = "2"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

fn main() {

    
    let val = crate::utils::interpolate_uni_poly(
        &vec![23, 24, 25].into_iter().map(|e| Fq::from(e as u64)).collect::<Vec<Fq>>(), Fq::from(20)
    );
    println!("VAL {val}");

    // eval vector in LE indexing
    // f = 1 - 10 * x1 * x2 + 5 * x1 + 4 * x2
    let f = DenseMultilinearExtension::from_evaluations_vec(
        2,
        vec![1, 6, 5, 0]
            .into_iter()
            .map(|x| Fq::from(x as u64))
            .collect(),
    );

    // g = 3 + x3 + 4 * x3 * x4 + x4
    let g = DenseMultilinearExtension::from_evaluations_vec(
        2,
            vec![3, 4, 4, 9]
            .into_iter()
            .map(|x| Fq::from(x as u64))
            .collect(),
    );

    // g = 3 + 2 * x1 + 4 * x1 * x2
    let g = DenseMultilinearExtension::from_evaluations_vec(
        2,
            vec![3, 5, 3, 9]
            .into_iter()
            .map(|x| Fq::from(x as u64))
            .collect(),
    );
    
    // write some assertions for pol evaluations
    assert_eq!(f.evaluate(&[Fq::from(0), Fq::from(0)]).unwrap(), Fq::from(1));
    assert_eq!(f.evaluate(&[Fq::from(1), Fq::from(0)]).unwrap(), Fq::from(6));
    assert_eq!(f.evaluate(&[Fq::from(2), Fq::from(0)]).unwrap(), Fq::from(11));
    assert_eq!(f.evaluate(&[Fq::from(0), Fq::from(1)]).unwrap(), Fq::from(5));

    run_sumcheck_protocol(f, g); // expected: 48
}

// if we implement this with 2 generic polys and GKR uses it with the same twice, we are wasting e.g. table computation
// gkr: two different points; sumcheck: same point?
// it is not algorithm three! that is for the product of two polys on the *same* variables
// needed to instantiate MLE
// failure of interpolation