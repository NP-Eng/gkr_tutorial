
#[cfg(test)]
mod tests {

    use ark_ff::fields::{Fp64, MontBackend, MontConfig};

    use ark_poly::{
        polynomial::multivariate::{SparsePolynomial, SparseTerm, Term}, DenseMVPolynomial,
    };
    
    use crate::{polynomial::*, parties::run_sumcheck_protocol};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_sq_mul() {
        assert_eq!(pow_sq_m(&Fq::from(4), 19), Fq::from(13));
    }

    #[test]
    fn test_specialise() {
        let g = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(4), SparseTerm::new(vec![(0, 3)])),
                (Fq::from(9), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(5), SparseTerm::new(vec![(2, 1)])),
                (Fq::from(7), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Fq::from(2), SparseTerm::new(vec![(1, 2), (2, 3)])),
                (Fq::from(1), SparseTerm::new(vec![])),
            ],
        );

        let f = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(4), SparseTerm::new(vec![(0, 3)])),
                (Fq::from(12), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(4), SparseTerm::new(vec![(0, 1)])),
                (Fq::from(16), SparseTerm::new(vec![])),
            ],
        );

        assert_eq!(specialise(&g, 2, Fq::from(3)), f);

    }

    #[test]
    fn test_univariate() {
        
        let g = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(4), SparseTerm::new(vec![(0, 3)])),
                (Fq::from(12), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(4), SparseTerm::new(vec![(0, 1)])),
                (Fq::from(16), SparseTerm::new(vec![])),
            ],
        );

        let f = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(4), SparseTerm::new(vec![(0, 0)])),
                (Fq::from(12), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(4), SparseTerm::new(vec![(1, 1), (2, 0)])),
                (Fq::from(16), SparseTerm::new(vec![])),
            ],
        );

        assert!(univariate_index(&g).is_err());
        assert_eq!(univariate_index(&f), Ok(Some(1)));
    }

    #[test]
    fn test_degree() {
        
        let g = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(4), SparseTerm::new(vec![(0, 3)])),
                (Fq::from(9), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(5), SparseTerm::new(vec![(2, 1)])),
                (Fq::from(7), SparseTerm::new(vec![(0, 1), (2, 1), (0, 3)])),
                (Fq::from(2), SparseTerm::new(vec![(1, 2), (2, 3)])),
                (Fq::from(1), SparseTerm::new(vec![])),
            ],
        );

        assert_eq!(partial_degree(&g, 0), 4);
        assert_eq!(partial_degree(&g, 1), 2);
        assert_eq!(partial_degree(&g, 2), 3);
    }

    #[test]
    fn test_univariate_eval() {
        let pol = UnivariatePolynomial::new(
            &vec![Fq::from(1), Fq::from(8), Fq::from(0), Fq::from(3), Fq::from(0)]
        );

        assert_eq!(pol.eval(Fq::from(2)), Fq::from(7));
    }

    #[test]
    fn test_univariate_from_multivariate() {
        let g = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(4), SparseTerm::new(vec![(0, 0)])),
                (Fq::from(12), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(1), SparseTerm::new(vec![(1, 2)])),
                (Fq::from(7), SparseTerm::new(vec![(1, 4)])),
                (Fq::from(4), SparseTerm::new(vec![(1, 1), (2, 0)])),
                (Fq::from(16), SparseTerm::new(vec![])),
            ],
        );

        let c = UnivariatePolynomial::from_sparse_multivariate(&g).unwrap();
        println!("{:?}", c.coefficients());

        //assert_eq!(UnivariatePolynomial::from_sparse_multivariate(&g).unwrap(), UnivariatePolynomial::new(vec!([Fq::from(3), Fq::from(4), Fq::from(13), Fq::zero(), Fq::from(7)])));
    }

    #[test]
    fn test_univariate_display() {
        let p1 = UnivariatePolynomial::new(
            &vec![Fq::from(1), Fq::from(8), Fq::from(0), Fq::from(3)]
        );

        let p2 = UnivariatePolynomial::new(
            &vec![Fq::from(0), Fq::from(1), Fq::from(8), Fq::from(0), Fq::from(3)]
        );

        let p3 = UnivariatePolynomial::new(
            &vec![Fq::from(0), Fq::from(0), Fq::from(1), Fq::from(8), Fq::from(0), Fq::from(3)]
        );

        assert_eq!(p1.to_string(), "3 * x^3 + 8 * x + 1");
        assert_eq!(p2.to_string(), "3 * x^4 + 8 * x^2 + 1 * x");
        assert_eq!(p3.to_string(), "3 * x^5 + 8 * x^3 + 1 * x^2");
    }

    #[test]
    fn test_protocol() {
        let pol = SparsePolynomial::from_coefficients_vec(
            3,
            vec![
                (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
                (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
            ],
        );
        assert!(run_sumcheck_protocol(pol, true))
    }

}
