use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension};

use super::{parties::Transcript as SC_Transcript, UniformCircuit};
use crate::polynomial::UnivariatePolynomial;

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
        let w_0 = self.circuit.evaluate(x);

        let r_0 = self.transcript.update(w_0.clone(), &mut self.sponge, 1);
        let w0_mle = DenseMultilinearExtension::from_evaluations_vec(
            self.circuit.layers[0].num_vars,
            w_0
        );

        let m = w0_mle.evaluate(&r_0); // TODO possibly unnecessary

        let d = self.circuit.layers.len();


        f(x, y) = f1(x) + f2(y);
        mle(f);

        eval_list = Vec::new();

        for c in 0..2**k0 {
            for b 0..2**k0 {
                
            }
        }
    
        w_i_b_plus_w_i_c = DenseMultilinearExtension::from_evaluations_vec(
            add_i_ri 
        )

        for i in 0..d {

            // compute wi+1
            // w_iplus1_mle = DenseMultilinearExtension::from_evaluations_vec(w_iplus1) + 

            let f_i_ri = 
                DenseMultilinearExtension::from_evaluations_vec(add_i_ri) *
                    w_iplus1_mle + w_iplus1_mle


        }

        unimplemented!()
    }
}

struct Verifier<F: PrimeField + Absorb> {
    circuit: UniformCircuit<F>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

mod test {

    #[test]
    fn run_prover() {}
}
