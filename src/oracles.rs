use ark_ff::PrimeField;

use ark_poly::MultilinearExtension;

use crate::parties::{Product, SumOfProducts};

pub trait Oracle<F> {
    fn divinate(&self, x: &[F], y: &[F]) -> F;
    fn num_rounds(&self) -> usize;
    fn trust_message(&self) -> String;
}

pub(crate) struct PolyOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: SumOfProducts<F, MLE>,
    g: Vec<F>,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> PolyOracle<F, MLE> {
    pub(crate) fn new(sum_of_products: SumOfProducts<F, MLE>, g: Vec<F>) -> Self {
        Self { sum_of_products, g }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for PolyOracle<F, MLE> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {
        let mut zxy: Vec<F> = self.g.clone();
        zxy.extend(x);
        zxy.extend(y);

        let mut sum = F::zero();

        for Product(f1, f2, f3) in &self.sum_of_products.terms {
            sum += f1.evaluate(&zxy).unwrap() * f2.evaluate(x).unwrap() * f3.evaluate(y).unwrap()
        }

        sum
    }

    // returns the number of variables of EACH of the two polynomials
    // the total number of variables of the product is twice that
    fn num_rounds(&self) -> usize {
        let Product(_, f2, _) = &self.sum_of_products.terms[0];
        2 * f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}

pub(crate) struct CombinedPolyOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: SumOfProducts<F, MLE>,
    g1: Vec<F>,
    g2: Vec<F>,
    alpha: F,
    beta: F,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> CombinedPolyOracle<F, MLE> {
    pub(crate) fn new(
        sum_of_products: SumOfProducts<F, MLE>,
        g1: Vec<F>,
        g2: Vec<F>,
        alpha: F,
        beta: F,
    ) -> Self {
        Self {
            sum_of_products,
            g1,
            g2,
            alpha,
            beta,
        }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for CombinedPolyOracle<F, MLE> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {
        let mut sums = Vec::new();

        for (coeff, g) in vec![(self.alpha, self.g1.clone()), (self.beta, self.g2.clone())] {
            let mut zxy: Vec<F> = g.clone();
            zxy.extend(x);
            zxy.extend(y);

            let mut sum = F::zero();

            for Product(f1, f2, f3) in &self.sum_of_products.terms {
                sum +=
                    f1.evaluate(&zxy).unwrap() * f2.evaluate(x).unwrap() * f3.evaluate(y).unwrap()
            }

            sums.push(coeff * sum);
        }

        sums[0] + sums[1]
    }

    // returns the number of variables of EACH of the two polynomials
    // the total number of variables of the product is twice that
    fn num_rounds(&self) -> usize {
        let Product(_, f2, _) = &self.sum_of_products.terms[0];
        2 * f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}

pub(crate) struct GKROracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: SumOfProducts<F, MLE>,
    g1: Vec<F>,
    g2: Vec<F>,
    alpha: F,
    beta: F,
    w_u_i: F,
    w_v_i: F,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> GKROracle<F, MLE> {
    pub(crate) fn new(
        sum_of_products: SumOfProducts<F, MLE>,
        g1: Vec<F>,
        g2: Vec<F>,
        alpha: F,
        beta: F,
        w_u_i: F,
        w_v_i: F,
    ) -> Self {
        Self {
            sum_of_products,
            g1,
            g2,
            alpha,
            beta,
            w_u_i,
            w_v_i,
        }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for GKROracle<F, MLE> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {
        let mut sums = Vec::new();

        for (coeff, g) in vec![(self.alpha, self.g1.clone()), (self.beta, self.g2.clone())] {
            let mut zxy: Vec<F> = g.clone();
            zxy.extend(x);
            zxy.extend(y);

            let mut sum = F::zero();

            // almost the same as in `CombinedPolyOracle`, but instead of using f2, f3, we use the provided values w_u_i, w_v_i.
            for Product(f1, _, _) in &self.sum_of_products.terms {
                sum += f1.evaluate(&zxy).unwrap() * self.w_u_i * self.w_v_i
            }

            sums.push(coeff * sum);
        }

        sums[0] + sums[1]
    }

    // returns the number of variables of the SparseMLE, minus the number of bits needed to represent "current" layer.
    fn num_rounds(&self) -> usize {
        let Product(f1, _, _) = &self.sum_of_products.terms[0];
        f1.num_vars() - self.g1.len()
    }

    fn trust_message(&self) -> String {
        format!("assuming that f(b*) and f(c*) are from the prover")
    }
}
