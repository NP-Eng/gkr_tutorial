#[cfg(test)]
mod tests {
    use crate::{
        parties::{initialise_phase_1, initialise_phase_2, precompute, run_sumcheck_protocol},
        utils::usize_to_zxy,
        Fq,
    };
    use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
    use ark_std::UniformRand;

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
            (0..4)
                .into_iter()
                .map(|_| Fq::rand(&mut test_rng))
                .collect(),
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
        let mut test_rng = ark_std::test_rng();
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
        let mut test_rng = ark_std::test_rng();
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
        let mut test_rng = ark_std::test_rng();
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
        let mut test_rng = ark_std::test_rng();
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

        let actual = initialise_phase_2(&f1, &g, &u);

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
    fn tests_sumcheck() {
        let mut test_rng = ark_std::test_rng();
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
}
