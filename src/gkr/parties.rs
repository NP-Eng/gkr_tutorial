use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};

use super::{parties::Transcript as SC_Transcript, UniformCircuit};
use crate::parties::{Product, Prover as SumcheckProver};

#[derive(Clone, Debug)]
pub struct Transcript<F: PrimeField + Absorb> {
    pub values: Vec<Vec<F>>,
    pub challenges: Vec<F>,
    pub sumcheck_transcripts: Vec<SC_Transcript<F>>,
}

impl<F: PrimeField + Absorb> Transcript<F> {
    fn new() -> Self {
        Self {
            values: Vec::new(),
            challenges: Vec::new(),
            sumcheck_transcripts: Vec::new(),
        }
    }
    fn update(&mut self, elems: Vec<F>, sponge: &mut PoseidonSponge<F>, n: usize) -> Vec<F> {
        for elem in elems.iter() {
            sponge.absorb(elem);
        }
        self.values.push(elems);

        // challenges simlated via the FS transform
        let rs = sponge.squeeze_field_elements::<F>(n);
        self.challenges.extend(rs.clone());
        rs
    }
    fn verify(&self, sponge: &mut PoseidonSponge<F>) -> bool {
        // not useful
        for (msg, h) in self.values.iter().zip(self.challenges.iter()) {
            for elem in msg.iter() {
                sponge.absorb(elem);
            }

            if sponge.squeeze_field_elements::<F>(1)[0] != *h {
                return false;
            }
        }

        true
    }
}
pub struct Prover<F: PrimeField + Absorb, const d: usize> {
    circuit: UniformCircuit<F, d>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

// new prover
impl<F: PrimeField + Absorb, const d: usize> Prover<F, d> {
    pub fn new(circuit: UniformCircuit<F, d>, sponge: PoseidonSponge<F>) -> Self {
        Self {
            circuit,
            sponge,
            transcript: Transcript::new(),
        }
    }
}

impl<F: PrimeField + Absorb, const d: usize> Prover<F, d> {
    fn run(&mut self, x: Vec<F>) -> Transcript<F> {
        let k = d;
        let identically_one =
            DenseMultilinearExtension::from_evaluations_vec(k, vec![F::one(); 1 << k]);

        // compute all the layers
        let w = self.circuit.evaluate(x);
        // TODO provably need to go in reverse index order for w, since the input layer is at index 0
        let r_0 = self.transcript.update(w[0].clone(), &mut self.sponge, 1);

        let mut alpha = F::one();
        let mut beta = F::zero();

        let mut u_i = r_0.clone();
        let mut v_i = r_0;

        for i in 0..k {
            let w_iplus1_mle = DenseMultilinearExtension::from_evaluations_vec(k, w[i + 1].clone());

            let [add_i_mle, mul_i_mle]: [SparseMultilinearExtension<F>; 2] =
                (&self.circuit.layers[i + 1]).into();

            // TODO lots of cloning here, can we do better?
            let f_i_1 = (
                add_i_mle.clone(),
                w_iplus1_mle.clone(),
                identically_one.clone(),
            );
            let f_i_2 = (add_i_mle, identically_one.clone(), w_iplus1_mle.clone());
            let f_i_3 = (mul_i_mle, w_iplus1_mle.clone(), w_iplus1_mle.clone());

            let mut sumcheck_prover = SumcheckProver::new(
                vec![f_i_1.into(), f_i_2.into(), f_i_3.into()].into(),
                self.sponge.clone(),
            );
            let sumcheck_transcript = sumcheck_prover.run(&u_i, &v_i, alpha, beta);
            // the first half of the transcript is b*, the second is c*
            let (b_star, c_star) = sumcheck_transcript.challenges.split_at(k / 2);

            let w_b_star = w_iplus1_mle.evaluate(&b_star).unwrap();
            let w_c_star = w_iplus1_mle.evaluate(&c_star).unwrap();

            u_i = b_star.to_vec();
            v_i = c_star.to_vec();

            // get the random scalars from the verifier to compute the linear combination
            let temp_scalars =
                self.transcript
                    .update([w_b_star, w_c_star].to_vec(), &mut self.sponge, 2);

            alpha = temp_scalars[0];
            beta = temp_scalars[1];
        }
        self.transcript.clone()
    }
}

struct Verifier<F: PrimeField + Absorb, const d: usize> {
    circuit: UniformCircuit<F, d>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

// impl run for the verifier
impl<F: PrimeField + Absorb, const d: usize> Verifier<F, d> {
    fn run(&mut self) -> bool {
        true
    }
}
