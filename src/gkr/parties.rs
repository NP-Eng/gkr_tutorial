use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::PrimeField;

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

        self.transcript.update(w_0, &mut self.sponge, 1);

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
