// mod gkr;
mod parties;
// mod tests;
mod utils;

use ark_ff::fields::{Fp64, MontBackend, MontConfig};

use ark_poly::{
    DenseMultilinearExtension, MultilinearExtension,
};

use parties::run_sumcheck_protocol;

#[derive(MontConfig)]
#[modulus = "32771"] // MODIFY
#[generator = "2"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

fn main() {
    
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

    /*     // g = 3 + 2 * x1 + 4 * x1 * x2
    let g = DenseMultilinearExtension::from_evaluations_vec(
        2,
            vec![3, 5, 3, 9]
            .into_iter()
            .map(|x| Fq::from(x as u64))
            .collect(),
    ); */
    
    // write some assertions for pol evaluations
    assert_eq!(f.evaluate(&[Fq::from(0), Fq::from(0)]).unwrap(), Fq::from(1));
    assert_eq!(f.evaluate(&[Fq::from(1), Fq::from(0)]).unwrap(), Fq::from(6));
    assert_eq!(f.evaluate(&[Fq::from(2), Fq::from(0)]).unwrap(), Fq::from(11));
    assert_eq!(f.evaluate(&[Fq::from(0), Fq::from(1)]).unwrap(), Fq::from(5));

    run_sumcheck_protocol(f, g); // expected: 48 */
}
