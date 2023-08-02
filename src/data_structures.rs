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

/// Product of a Sparse MLE and two other MLEs.
/// E.g. for GKR, Sparse MLE is either add/mul MLE, and the other two are W_{i+1}
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

/// Sum of products as seen by the verifier.
/// For example, for GKR, we have a sum of products of 3 elements:
/// 1st: add_i(g, x, y) * W_{i+1}(x) * 1
/// 2nd: add_i(g, x, y) * 1 * W_{i+1}(y)
/// 3rd  mul_i(g, x, y) * W_{i+1}(x) * W_{i+1}(y)
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

/// Sum of products as seen by the prover.
/// For example, for GKR, we have a sum of products of 3 elements:
/// 1st: add(g, x, y) * f2(x) * 1
/// 2nd: add(g, x, y) * 1 * f3(y)
/// 3rd  mul(g, x, y) * f2(x) * f3(y)
#[derive(Clone, Debug)]
pub struct SumOfProducts<F: PrimeField, MLE: MultilinearExtension<F>> {
    pub terms: Vec<Product<F, MLE>>,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> From<Vec<Product<F, MLE>>>
    for SumOfProducts<F, MLE>
{
    fn from(terms: Vec<Product<F, MLE>>) -> Self {
        Self { terms }
    }
}
