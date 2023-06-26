use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};

use super::{parties::Transcript as SC_Transcript, UniformCircuit};
use crate::parties::Prover as SumcheckProver;

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
struct Prover<F: PrimeField + Absorb> {
    circuit: UniformCircuit<F>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

impl<F: PrimeField + Absorb> Prover<F> {
    fn run(&mut self, x: Vec<F>) -> Transcript<F> {
        let k = self.circuit.layers[0].num_vars;
        let identically_one =
            DenseMultilinearExtension::from_evaluations_vec(k, vec![F::one(); 1 << k]);

        // compute all the layers
        let w = self.circuit.evaluate(x);
        let r_0 = self.transcript.update(w[0].clone(), &mut self.sponge, 1);

        // init random scalars to 1, 1
        let mut random_scalars = vec![F::one(), F::one()];

        for i in 0..k {
            let alpha = random_scalars[0];
            let beta = random_scalars[1];

            let w_iplus1_mle = DenseMultilinearExtension::from_evaluations_vec(k, w[i + 1]);

            let [add_i_mle, mul_i_mle]: [SparseMultilinearExtension<F>; 2] =
                self.circuit.layers[i + 1].into();

            // TODO multiply by alpha and beta
            let f_i_1 = (add_i_mle, w_iplus1_mle, identically_one);
            let f_i_2 = (add_i_mle, identically_one, w_iplus1_mle);
            let f_i_3 = (mul_i_mle, w_iplus1_mle, w_iplus1_mle);

            let sumcheck_prover = SumcheckProver::new(vec![f_i_1, f_i_2, f_i_3]);
            let sumcheck_transcript = sumcheck_prover.run(&r_0);
            // the first half of the transcript is b*, the second is c*
            let (b_star, c_star) = sumcheck_transcript.challenges.split_at(k / 2);

            let w_b_star = w_iplus1_mle.evaluate(&b_star).unwrap();
            let w_c_star = w_iplus1_mle.evaluate(&c_star).unwrap();

            // get the random scalars from the verifier to compute the linear combination
            random_scalars =
                self.transcript
                    .update([w_b_star, w_c_star].to_vec(), &mut self.sponge, 2);
        }
        self.transcript
    }
}

struct Verifier<F: PrimeField + Absorb> {
    circuit: UniformCircuit<F>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

// impl run for the verifier
impl<F: PrimeField + Absorb> Verifier<F> {
    fn run(&mut self) -> bool {
        true
    }
}

mod test {

    #[test]
    fn run_prover() {}
}
