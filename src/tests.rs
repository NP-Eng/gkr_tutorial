#[cfg(test)]
mod tests {
    use crate::{
        data_structures::{Product, SumOfProducts},
        parties::{
            initialise_phase_1 as init_phase_1_original,
            initialise_phase_2 as init_phase_2_original, run_sumcheck_protocol,
            run_sumcheck_protocol_combined, run_sumcheck_protocol_combined_multiprod,
        },
        utils::usize_to_zxy,
        Fr as Fq,
    };
    use ark_ff::PrimeField;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
    use ark_std::UniformRand;
    use halo2_base::halo2_proofs::halo2curves::bn256::Fr as Fr2;

    fn initialise_phase_1<F: PrimeField, MLE: MultilinearExtension<F>>(
        f_1: &SparseMultilinearExtension<F>,
        f_3: &MLE,
        g: &[F],
    ) -> Vec<F> {
        init_phase_1_original(f_1, f_3, g, g, F::one(), F::zero())
    }

    fn initialise_phase_2<F: PrimeField>(
        f_1: &SparseMultilinearExtension<F>,
        g: &[F],
        u: &[F],
    ) -> Vec<F> {
        init_phase_2_original(f_1, g, g, u, F::one(), F::zero())
    }

    #[test]
    fn test_int_to_zxy() {
        let v = 2;
        let (z0, x0, y0) = (3, 0, 1);

        let vp = 1 << v;
        let idx = z0 * (vp * vp) + x0 * vp + y0;

        assert_eq!((z0, x0, y0), usize_to_zxy(idx, v));
    }

    #[test]
    fn tests_phase_1_should_be_0() {
        let mut test_rng = ark_std::test_rng();
        // f(x1, x2, x3, x4, x5, x6) = f(g1, g2, x1, x2, y1, y2)
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        // when first 2 variables are 1, 1, then f1 should be 0 no matter what f3 is
        let g = vec![Fq::from(1u64), Fq::from(1u64)];

        // f3 is a random 2-variate MLE
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            (0..4).map(|_| Fq::rand(&mut test_rng)).collect(),
        );
        let actual = initialise_phase_1(&f1, &f3, &g);

        // computed by hand
        let expected = vec![
            Fq::from(0u64),
            Fq::from(0u64),
            Fq::from(0u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn tests_phase_1() {
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        // f3 would be a constant 2-variate MLE
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![3, 3, 3, 3]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let g = vec![Fq::from(0u64), Fq::from(0u64)];

        let actual = initialise_phase_1(&f1, &f3, &g);

        // computed by hand
        let expected = vec![
            Fq::from(6u64),
            Fq::from(0u64),
            Fq::from(3u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn tests_phase_1_full() {
        // f(x1, x2, x3, x4, x5, x6) = f(g1, g2, x1, x2, y1, y2)
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        // when first 2 variables are 1, 1, then f1 should be 0 no matter what f3 is
        let g = vec![Fq::from(0u64), Fq::from(0u64)];

        // eval vector in LE indexing
        // f3(y1, y2) = 1 - 10 * y1 * y2 + 5 * y1 + 4 * y2
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );
        let actual = initialise_phase_1(&f1, &f3, &g);

        // computed by hand
        let expected = vec![
            Fq::from(10u64),
            Fq::from(0u64),
            Fq::from(1u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn tests_phase_1_full_non_zero_g() {
        // f(x1, x2, x3, x4, x5, x6) = f(g1, g2, x1, x2, y1, y2)
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(2u64), Fq::from(3u64)];

        // eval vector in LE indexing
        // f3(y1, y2) = 1 - 10 * y1 * y2 + 5 * y1 + 4 * y2
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );
        let actual = initialise_phase_1(&f1, &f3, &g);

        // computed by hand
        let expected = vec![
            Fq::from(20u64),
            Fq::from(0u64),
            Fq::from(2u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn tests_phase_1_max() {
        // f1 = x1 * (1 - 10 * g1 * g2 + 5 + 4 * x2 + y1 + g2 * y2 - 7 * y2 * x2)
        let points = vec![
            (4, Fq::from(6_u64)),
            (5, Fq::from(6_u64)),
            (6, Fq::from(6_u64)),
            (7, Fq::from(97_u64)),
            (12, Fq::from(10_u64)),
            (13, Fq::from(10_u64)),
            (14, Fq::from(10_u64)),
            (20, Fq::from(7_u64)),
            (21, Fq::from(7_u64)),
            (22, Fq::from(7_u64)),
            (23, Fq::from(98_u64)),
            (28, Fq::from(11_u64)),
            (29, Fq::from(11_u64)),
            (30, Fq::from(11_u64)),
            (31, Fq::from(1_u64)),
            (36, Fq::from(6_u64)),
            (37, Fq::from(6_u64)),
            (38, Fq::from(7_u64)),
            (39, Fq::from(98_u64)),
            (44, Fq::from(3_u64)),
            (45, Fq::from(3_u64)),
            (46, Fq::from(4_u64)),
            (47, Fq::from(95_u64)),
            (52, Fq::from(7_u64)),
            (53, Fq::from(7_u64)),
            (54, Fq::from(8_u64)),
            (55, Fq::from(99_u64)),
            (60, Fq::from(4_u64)),
            (61, Fq::from(4_u64)),
            (62, Fq::from(5_u64)),
            (63, Fq::from(96_u64)),
        ];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        // when first 2 variables are 1, 1, then f1 should be 0 no matter what f3 is
        let g = vec![Fq::from(0u64), Fq::from(0u64)];

        // eval vector in LE indexing
        // f3(y1, y2) = 1 - 10 * y1 * y2 + 5 * y1 + 4 * y2
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );
        let actual = initialise_phase_1(&f1, &f3, &g);

        // computed by hand
        let expected = vec![
            Fq::from(0u64),
            Fq::from(78u64),
            Fq::from(0u64),
            Fq::from(91u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn tests_phase_2_max() {
        // f1 = x1 * (1 - 10 * g1 * g2 + 5 + 4 * x2 + y1 + g2 * y2 - 7 * y2 * x2)
        let points = vec![
            (4, Fq::from(6_u64)),
            (5, Fq::from(6_u64)),
            (6, Fq::from(6_u64)),
            (7, Fq::from(97_u64)),
            (12, Fq::from(10_u64)),
            (13, Fq::from(10_u64)),
            (14, Fq::from(10_u64)),
            (20, Fq::from(7_u64)),
            (21, Fq::from(7_u64)),
            (22, Fq::from(7_u64)),
            (23, Fq::from(98_u64)),
            (28, Fq::from(11_u64)),
            (29, Fq::from(11_u64)),
            (30, Fq::from(11_u64)),
            (31, Fq::from(1_u64)),
            (36, Fq::from(6_u64)),
            (37, Fq::from(6_u64)),
            (38, Fq::from(7_u64)),
            (39, Fq::from(98_u64)),
            (44, Fq::from(3_u64)),
            (45, Fq::from(3_u64)),
            (46, Fq::from(4_u64)),
            (47, Fq::from(95_u64)),
            (52, Fq::from(7_u64)),
            (53, Fq::from(7_u64)),
            (54, Fq::from(8_u64)),
            (55, Fq::from(99_u64)),
            (60, Fq::from(4_u64)),
            (61, Fq::from(4_u64)),
            (62, Fq::from(5_u64)),
            (63, Fq::from(96_u64)),
        ];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(80u64), Fq::from(6u64)];
        let u = vec![Fq::from(27u64), Fq::from(12u64)];

        let actual = initialise_phase_2::<Fq>(&f1, &g, &u);

        // computed by hand
        let expected = vec![
            Fq::from(27_u64),
            Fq::from(54_u64),
            Fq::from(42_u64),
            Fq::from(69_u64),
        ];

        assert_eq!(actual, expected);
    }

    #[test]
    fn tests_sumcheck_basic() {
        // f(x1, x2, x3, x4, x5, x6) = f(g1, g2, x1, x2, y1, y2)
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(2u64), Fq::from(3u64)];

        let f2 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );
        // eval vector in LE indexing
        // f3(y1, y2) = 1 - 10 * y1 * y2 + 5 * y1 + 4 * y2
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        run_sumcheck_protocol(f1, f2, f3, &g);
    }

    #[test]
    fn tests_sumcheck_randomized() {
        let mut test_rng = ark_std::test_rng();
        // degree_x is the degree of both f2 (in x) and f3 (in y)
        // only go up to degree 7, otherwise tests start taking long
        for degree_x in 3..8 {
            // total degree of f1
            let d = 3 * degree_x;
            // create a random sparse MLE with 2^(d-2) evaluation points (25% density)
            let num_evals: usize = 1 << (d - 2);
            let points = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f1 = SparseMultilinearExtension::from_evaluations(d, &points);
            // generate random g
            let g = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate a random f2 from Vec<Fq>
            let f2 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            // same for f3
            let f3 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );

            run_sumcheck_protocol(f1, f2, f3, &g);
        }
    }

    #[test]
    fn tests_sumcheck_trivial_combination() {
        // same parameters as the test tests_sumcheck, but using the linear-combination interface:
        // first with coefficients (1, 0), then (0, 1)

        // f(x1, x2, x3, x4, x5, x6) = f(g1, g2, x1, x2, y1, y2)
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(2u64), Fq::from(3u64)];

        let f2 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );
        // eval vector in LE indexing
        // f3(y1, y2) = 1 - 10 * y1 * y2 + 5 * y1 + 4 * y2
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let zero_vector = vec![Fq::from(0_u64); g.len()];
        let ones_vector = vec![Fq::from(0_u64); g.len()];

        run_sumcheck_protocol_combined(
            f1.clone(),
            f2.clone(),
            f3.clone(),
            &g,
            &zero_vector,
            Fq::from(1u64),
            Fq::from(0_u64),
        );
        run_sumcheck_protocol_combined(
            f1,
            f2,
            f3,
            &ones_vector,
            &g,
            Fq::from(0_u64),
            Fq::from(1u64),
        );
    }

    #[test]
    fn tests_sumcheck_randomized_combination() {
        let mut test_rng = ark_std::test_rng();
        // degree_x is the degree of both f2 (in x) and f3 (in y)
        // only go up to degree 7, otherwise tests start taking long
        for degree_x in 3..8 {
            // total degree of f1
            let d = 3 * degree_x;
            // create a random sparse MLE with 2^(d-2) evaluation points (25% density)
            let num_evals: usize = 1 << (d - 2);
            let points = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f1 = SparseMultilinearExtension::from_evaluations(d, &points);
            // generate random g1
            let g1 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate random g2
            let g2 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate a random f2 from Vec<Fq>
            let f2 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            // same for f3
            let f3 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );

            run_sumcheck_protocol_combined(
                f1.clone(),
                f2.clone(),
                f3.clone(),
                &g1,
                &g2,
                Fq::rand(&mut test_rng),
                Fq::rand(&mut test_rng),
            );
        }
    }

    #[test]
    fn tests_sumcheck_randomized_combination_multiprod() {
        let mut test_rng = ark_std::test_rng();
        // degree_x is the degree of both f2 (in x) and f3 (in y)
        // only go up to degree 7, otherwise tests start taking long
        for degree_x in 3..8 {
            // total degree of f1
            let d = 3 * degree_x;
            // create a random sparse MLE with 2^(d-2) evaluation points (25% density)
            let num_evals: usize = 1 << (d - 2);
            let points = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f11 = SparseMultilinearExtension::from_evaluations(d, &points);
            // generating second sparse f1
            let _points2 = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f12 = SparseMultilinearExtension::from_evaluations(d, &points);
            // generate random g1
            let g1 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate random g2
            let g2 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();

            // generating first f2 from random coeffs
            let f21 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            // same for first f3
            let f31 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            // same for second f2
            let f22 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            // same for second f3
            let f32 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );

            let sum_of_products = SumOfProducts {
                terms: vec![
                    Product(f11.clone(), f21.clone(), f31.clone()),
                    Product(f12.clone(), f22.clone(), f32.clone()),
                ],
            };

            run_sumcheck_protocol_combined_multiprod(
                sum_of_products,
                &g1,
                &g2,
                Fq::rand(&mut test_rng),
                Fq::rand(&mut test_rng),
            );
        }
    }
}
