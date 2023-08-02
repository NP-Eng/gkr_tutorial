use ark_ff::PrimeField;

use ark_poly::{MultilinearExtension, SparseMultilinearExtension};

/// Like `Product`, but f2, f3 aren't known to the verifier, only their evals at a certain point
#[derive(Clone, Debug)]
pub struct GKRVerifierProduct<F: PrimeField>(pub SparseMultilinearExtension<F>, pub F, pub F);

impl<F: PrimeField> From<(SparseMultilinearExtension<F>, F, F)> for GKRVerifierProduct<F> {
    fn from(tuple: (SparseMultilinearExtension<F>, F, F)) -> Self {
        Self(tuple.0, tuple.1, tuple.2)
    }
}

#[derive(Clone, Debug)]
pub struct Product<F: PrimeField, MLE: MultilinearExtension<F>>(
    pub SparseMultilinearExtension<F>,
    pub MLE,
    pub MLE,
);

impl<F: PrimeField, MLE: MultilinearExtension<F>> From<(SparseMultilinearExtension<F>, MLE, MLE)>
    for Product<F, MLE>
{
    fn from(tuple: (SparseMultilinearExtension<F>, MLE, MLE)) -> Self {
        assert_eq!(
            tuple.1.num_vars(),
            tuple.2.num_vars(),
            "f2 and f3 must have the same number of variables"
        );
        Self(tuple.0, tuple.1, tuple.2)
    }
}

#[derive(Clone, Debug)]
pub struct GKRVerifierSumOfProducts<F: PrimeField> {
    // for add:
    // 1st: [alpha * add(g1, x, y) + beta * add(g2, x, y)] * f2(x) * 1 +
    // 2nd: [alpha * add(g1, x, y) + beta * add(g2, x, y)] * 1 * f3(y) +
    // 3rd [alpha * mul(g1, x, y) + beta * mul(g2, x, y)] * f2(x) * f3(y)
    pub terms: Vec<GKRVerifierProduct<F>>,
}

impl<F: PrimeField> From<Vec<GKRVerifierProduct<F>>> for GKRVerifierSumOfProducts<F> {
    fn from(terms: Vec<GKRVerifierProduct<F>>) -> Self {
        Self { terms }
    }
}

#[derive(Clone, Debug)]
pub struct SumOfProducts<F: PrimeField, MLE: MultilinearExtension<F>> {
    // for add:
    // 1st: [alpha * add(g1, x, y) + beta * add(g2, x, y)] * f2(x) * 1 +
    // 2nd: [alpha * add(g1, x, y) + beta * add(g2, x, y)] * 1 * f3(y) +
    // 3rd [alpha * mul(g1, x, y) + beta * mul(g2, x, y)] * f2(x) * f3(y)
    pub terms: Vec<Product<F, MLE>>,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> From<Vec<Product<F, MLE>>>
    for SumOfProducts<F, MLE>
{
    fn from(terms: Vec<Product<F, MLE>>) -> Self {
        Self { terms }
    }
}
