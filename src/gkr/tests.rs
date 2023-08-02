#[cfg(test)]
mod test {

    use crate::{
        gkr::{
            parties::{Prover, Verifier},
            *,
        },
        utils::test_sponge,
    };
    use ark_bls12_381::Fq;
    use ark_crypto_primitives::sponge::Absorb;

    fn make_test_circuit<F: PrimeField + Absorb>() -> UniformCircuit<F, 2> {
        // example from Thaler's book p. 60, bottom - one gate changed for addition

        // mul, layer 0 (output)
        let mul0_0 = Wiring::new(0, 0, 1);
        let mul0_1 = Wiring::new(1, 2, 3);

        let layer_0 = Layer::new(Vec::new(), vec![mul0_0, mul0_1]);

        // add, layer 0
        // empty

        // mul, layer 1
        let mul1_0 = Wiring::new(0, 0, 0);
        let mul1_1 = Wiring::new(1, 1, 1);
        let mul1_2 = Wiring::new(2, 1, 2);

        // add, layer 1
        let add1_3 = Wiring::new(3, 3, 3);

        let layer_1 = Layer::new(vec![add1_3], vec![mul1_0, mul1_1, mul1_2]);

        UniformCircuit::<F, 2>::new(vec![layer_0, layer_1])
    }

    #[test]
    fn simple_circuit() {
        let circuit = make_test_circuit::<Fq>();
        let computed_out = circuit.evaluate(
            vec![3, 2, 3, 1]
                .iter()
                .map(|x| Fq::from(*x as u64))
                .collect(),
        );

        assert_eq!(
            *computed_out.last().unwrap(),
            vec![36, 12, 18, 18]
                .iter()
                .map(|x| Fq::from(*x as u64))
                .collect::<Vec<Fq>>()
        );
    }

    #[test]
    fn layer_to_mles() {
        let circuit = make_test_circuit::<Fq>();
        for layer in circuit.layers.iter() {
            let _: [SparseMultilinearExtension<Fq>; 2] = layer.into();
        }
    }

    #[test]
    fn test_gkr() {
        // Need to use the same sponge, since it's initialized with random values
        let sponge = test_sponge();

        let mut gkr_prover = Prover::<Fq, 2>::new(make_test_circuit(), sponge.clone());
        let circuit_input = vec![3, 2, 3, 1]
            .iter()
            .map(|x| Fq::from(*x as u64))
            .collect();
        let transcript = gkr_prover.run(circuit_input);

        let mut gkr_verifier =
            Verifier::<Fq, 2>::new(make_test_circuit(), sponge.clone(), transcript.proof);
        assert!(gkr_verifier.run());
    }
}
