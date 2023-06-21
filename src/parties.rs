use std::marker::PhantomData;

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::PrimeField;

use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm},
    Polynomial,
};
use ark_poly::{
    DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, SparseMultilinearExtension,
};

use crate::utils::{interpolate_uni_poly, test_sponge, Transcript};

pub struct Prover<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> {
    g: MLE,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

pub trait Oracle<F> {
    fn divinate(&self, x: &[F]) -> F;
    fn num_vars(&self) -> usize;
}

struct PolyOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    poly: MLE,
    phantom: PhantomData<F>,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> PolyOracle<F, MLE> {
    fn new(poly: MLE) -> Self {
        Self {
            poly,
            phantom: PhantomData,
        }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for PolyOracle<F, MLE> {
    fn divinate(&self, x: &[F]) -> F {
        self.poly.evaluate(&x.to_vec()).unwrap()
    }

    fn num_vars(&self) -> usize {
        self.poly.num_vars()
    }
}

pub struct Verifier<F: PrimeField + Absorb, O: Oracle<F>> {
    oracle: O,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

impl<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> Prover<F, MLE> {
    fn run(&mut self) -> Transcript<F> {
        let v = self.g.num_vars();

        let mut a = self.g.to_evaluations();

        let mut f = vec![vec![F::zero(); 1 << (v - 1)]; 3];

        // first round; claimed_value contains a constant univariate polynomial corresponding to the value the Prover wants to prove
        let claimed_value = a.iter().sum();

        self.transcript
            .update(&[claimed_value], &mut self.sponge);

        // preparing index conversion big-endian -> little-endian
        // we do this because `DenseMultilinearExtension::from_evaluations_slice` expects little-endian notation for the evaluation slice
        let le_indices: Vec<usize> = (0usize..(1 << v))
            .into_iter()
            .map(|i| i.reverse_bits() >> (usize::BITS as usize - v))
            .collect();

        for i in 1..=v {
            // r-th round

            // Algorithm 1

            for b in 0..(1 << (v - i)) {
                f[0][b] = a[le_indices[b]];
                f[1][b] = a[le_indices[b + (1 << (v - i))]];
                f[2][b] = -a[le_indices[b]] + a[le_indices[b + (1 << (v - i))]] * F::from(2u64);
            }

            let values: Vec<F> = (0..=2)
                .map(|t| f[t].iter().take(1 << (v - i)).sum())
                .collect();

            let r_i = self.transcript.update(&values, &mut self.sponge);

            for b in 0..(1 << (v - i)) {
                a[le_indices[b]] =
                    a[le_indices[b]] * (F::one() - r_i) + a[le_indices[b + (1 << (v - i))]] * r_i;
            }
        }

        println!("Prover finished successfully");

        self.transcript.clone()
    }
}

impl<F: PrimeField + Absorb, O: Oracle<F>> Verifier<F, O> {
    fn run(&mut self) -> bool {
        self.transcript.verify(&mut self.sponge);

        let v = self.oracle.num_vars();

        // this is actually a univariate polynomial
        let mut last_pol = Vec::from(&self.transcript.values[0][..]);

        let mut last_sent_scalar = self.transcript.challenges[0];

        // round numbering in line with the book's (necessarily inconsistent variable indexing!)
        for r in 1..=v {
            let new_pol_evals = self.transcript.values[r].clone();

            let claimed_sum = new_pol_evals[0] + new_pol_evals[1];
            // let new_pol = DenseMultilinearExtension::from_evaluations_slice(1, &new_pol_evals[..2]);
            let new_pol = Vec::from(&new_pol_evals[..]);

            if claimed_sum != interpolate_uni_poly(&last_pol, last_sent_scalar) {
                println!("Verifier found inconsistent evaluation in round {r}: received univariate polynomial is f_n = {new_pol:?},\n  \
                          previous one is f_o = {last_pol:?}. The equation f_n(0) + f_n(1) = f_o({last_sent_scalar}) fails to hold. Aborting."
                        );
                return false;
            }

            // cannot fail to unwrap due to round count
            last_sent_scalar = self.transcript.challenges[r];

            last_pol = new_pol;
        }

        if interpolate_uni_poly(&last_pol, last_sent_scalar)
            != self.oracle.divinate(&self.transcript.challenges[1..])
        {
            println!(
                "Verifier found inconsistent evaluation in the oracle call: \
                        received univariate polynomial {last_pol:?} evaluates to {},\n  \
                        whereas original multivariate one evaluates to {} on {:?}. Aborting.",
                interpolate_uni_poly(&last_pol, last_sent_scalar),
                self.oracle.divinate(&self.transcript.challenges[1..]),
                self.transcript.challenges[1..].to_vec(),
            );
            return false;
        }

        println!(
            "Verifier finished successfully and is confident in the value {}",
            self.transcript.values[0][0]
        );

        true
    }
}

// run the protocol and return true iff the verifier does not abort
pub fn run_sumcheck_protocol<F: PrimeField + Absorb, MLE: MultilinearExtension<F>>(pol: MLE) {
    println!("Initiated sumcheck protocol on the polynomial {pol:?}");

    let mut prover = Prover {
        g: pol.clone(),
        sponge: test_sponge(),
        transcript: Transcript::new(),
    };

    let transcript = prover.run();

    let mut verifier = Verifier {
        oracle: PolyOracle::new(pol),
        sponge: test_sponge(),
        transcript,
    };

    assert!(verifier.run());
}
