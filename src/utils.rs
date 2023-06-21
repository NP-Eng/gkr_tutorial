use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    Absorb, CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_std::test_rng;

use crate::polynomial::UnivariatePolynomial;

pub(crate) fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;
    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];
    let mut v = Vec::new();
    let mut ark_rng = test_rng();
    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();
        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}

#[derive(Clone, Debug)]
pub struct Transcript<F: PrimeField + Absorb> {
    pub values: Vec<Vec<F>>,
    pub challenges: Vec<F>,
}

impl<F: PrimeField + Absorb> Transcript<F> {
    pub fn new() -> Self {
        Self {
            values: Vec::new(),
            challenges: Vec::new(),
        }
    }
    
    pub fn update(&mut self, elements: &[F], sponge: &mut PoseidonSponge<F>) -> F {

        for elem in elements.iter() {
            sponge.absorb(elem);
        }
        self.values.push(elements.to_vec());

        let r = sponge.squeeze_field_elements::<F>(1)[0];
        self.challenges.push(r);
        r
    }

    pub fn verify(&self, sponge: &mut PoseidonSponge<F>) -> bool {
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
