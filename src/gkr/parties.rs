use std::{marker::PhantomData, vec};

use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};
use halo2_base::utils::ScalarField;
use poseidon_native::Poseidon;

use super::UniformCircuit;
use crate::utils::{test_sponge, SumcheckProof};
use crate::{
    oracles::GKROracle,
    parties::{Prover as SumcheckProver, Verifier as SumcheckVerifier},
};

#[derive(Clone, Debug, Default)]
pub struct GKRProof<F: PrimeField> {
    pub values: Vec<Vec<F>>,
    pub sumcheck_proofs: Vec<SumcheckProof<F>>,
}

/// Each party will instantiate a `Transcript`.
/// The prover will feed into the transcript the GKR layer claims, and separately the sumcheck messages, thus creating the `GKRProof`.
/// The verifier will feed into the transcript the values received in `GKRProof`.
#[derive(Clone, Debug)]
pub struct Transcript<F: PrimeField + Absorb, F2: ScalarField> {
    pub proof: GKRProof<F>,
    pub sponge: Poseidon<F2, 3, 2>,
}

impl<F: PrimeField + Absorb, F2: ScalarField> Transcript<F, F2> {
    fn new() -> Self {
        Self {
            proof: GKRProof::default(),
            sponge: test_sponge(),
        }
    }
    fn update(&mut self, elems: &[F], n: usize) -> Vec<F> {
        for elem in elems.iter() {
            let bytes = BigInteger::to_bytes_le(&elem.into_bigint());
            let halo2_base_elem = F2::from_bytes_le(&bytes);
            self.sponge.update(&[halo2_base_elem]);
        }
        self.proof.values.push(elems.to_vec());

        let mut squeezed = vec![];
        for i in 0..n {
            let halo2_out = self.sponge.squeeze();
            let out = F::from_le_bytes_mod_order(&halo2_out.to_bytes_le());
            squeezed.push(out);
        }
        squeezed
    }

    pub fn new_from_proof(proof: GKRProof<F>) -> Transcript<F, F2> {
        Self {
            proof,
            sponge: test_sponge(),
        }
    }
}
pub struct Prover<F: PrimeField + Absorb, F2: ScalarField, const s: usize> {
    circuit: UniformCircuit<F, s>,
    transcript: Transcript<F, F2>,
    phantom: PhantomData<F2>,
}

impl<F: PrimeField + Absorb, F2: ScalarField, const s: usize> Prover<F, F2, s> {
    pub fn new(circuit: UniformCircuit<F, s>) -> Self {
        Self {
            circuit,
            transcript: Transcript::new(),
            phantom: PhantomData,
        }
    }
}

impl<F: PrimeField + Absorb, F2: ScalarField, const s: usize> Prover<F, F2, s> {
    pub fn run(&mut self, x: Vec<F>) -> GKRProof<F> {
        // number of layers
        let d = self.circuit.layers.len();

        let identically_one =
            DenseMultilinearExtension::from_evaluations_vec(s, vec![F::one(); 1 << s]);

        // compute all the layers
        let w = self.circuit.evaluate(x);
        // go over the circuit layers in reverse index order for w, since the input layer is at index 0
        let r_0 = self.transcript.update(&w[d - 1], s);

        let mut alpha = F::one();
        let mut beta = F::zero();

        let mut u_i = r_0.clone();
        let mut v_i = r_0;

        for i in 0..d {
            let w_iplus1_mle =
                DenseMultilinearExtension::from_evaluations_vec(s, w[d - i - 1].clone());

            let [add_i_mle, mul_i_mle]: [SparseMultilinearExtension<F>; 2] =
                (&self.circuit.layers[i]).into();

            // TODO lots of cloning here, can we do better?
            let f_i_1 = (
                add_i_mle.clone(),
                w_iplus1_mle.clone(),
                identically_one.clone(),
            );
            let f_i_2 = (add_i_mle, identically_one.clone(), w_iplus1_mle.clone());
            let f_i_3 = (mul_i_mle, w_iplus1_mle.clone(), w_iplus1_mle.clone());

            let mut sumcheck_prover: SumcheckProver<F, F2, DenseMultilinearExtension<F>> =
                SumcheckProver::new(vec![f_i_1.into(), f_i_2.into(), f_i_3.into()].into());
            let (sumcheck_proof, random_sumcheck_challenges) =
                sumcheck_prover.run(&u_i, &v_i, alpha, beta);

            // the first half of the transcript is b*, the second is c*
            let (b_star, c_star) = random_sumcheck_challenges.split_at((1 << s) / 2);

            let w_b_star = w_iplus1_mle.evaluate(b_star).unwrap();
            let w_c_star = w_iplus1_mle.evaluate(c_star).unwrap();

            u_i = b_star.to_vec();
            v_i = c_star.to_vec();

            // get the random scalars from the verifier to compute the linear combination
            let temp_scalars = self.transcript.update(&[w_b_star, w_c_star], 2);

            alpha = temp_scalars[0];
            beta = temp_scalars[1];

            // add the sumcheck transcript to the GKR transcript
            self.transcript.proof.sumcheck_proofs.push(sumcheck_proof);
        }

        self.transcript.proof.clone()
    }
}

pub struct Verifier<F: PrimeField + Absorb, F2: ScalarField, const s: usize> {
    circuit: UniformCircuit<F, s>,
    transcript: Transcript<F, F2>,
    phantom: PhantomData<F2>,
}

// impl run for the verifier
impl<F: PrimeField + Absorb, F2: ScalarField, const s: usize> Verifier<F, F2, s> {
    pub fn new(circuit: UniformCircuit<F, s>, proof: GKRProof<F>) -> Self {
        Self {
            circuit,
            transcript: Transcript::new_from_proof(proof),
            phantom: PhantomData,
        }
    }

    pub fn run(&mut self) -> bool {
        let d = self.circuit.layers.len();

        let claim = self.transcript.proof.values[0].clone();

        let r_0 = self.transcript.update(&claim, s);

        // verify the transcript
        let mut u_i = r_0.clone();
        let mut v_i = r_0;

        let mut alpha = F::one();
        let mut beta = F::zero();

        for i in 0..d {
            let [add_i_mle, mul_i_mle]: [SparseMultilinearExtension<F>; 2] =
                (&self.circuit.layers[i]).into();

            let w_u_i = self.transcript.proof.values[i + 1][0];
            let w_v_i = self.transcript.proof.values[i + 1][1];

            let f_i_1 = (add_i_mle.clone(), w_u_i, F::one());
            let f_i_2 = (add_i_mle, F::one(), w_v_i);
            let f_i_3 = (mul_i_mle, w_u_i, w_v_i);
            let sum_of_prods = vec![f_i_1.into(), f_i_2.into(), f_i_3.into()].into();

            let mut sumcheck_verifier: SumcheckVerifier<F, F2, GKROracle<F>> =
                SumcheckVerifier::new(
                    GKROracle::new(sum_of_prods, u_i, v_i, alpha, beta),
                    self.transcript.proof.sumcheck_proofs[i].clone(),
                );

            let (result, random_challenges) = sumcheck_verifier.run();
            assert!(result, "sumcheck failed at round {}", i);

            let (tv1, tv2) = random_challenges.split_at((1 << s) / 2);
            u_i = tv1.to_vec();
            v_i = tv2.to_vec();

            let temp_scalars = self.transcript.update(&[w_u_i, w_v_i], 2);

            alpha = temp_scalars[0];
            beta = temp_scalars[1];
        }
        true
    }
}
