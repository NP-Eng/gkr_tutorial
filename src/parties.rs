use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::PrimeField;

use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm},
    Polynomial,
};

use crate::polynomial::{partial_degree, specialise, UnivariatePolynomial};
use crate::utils::test_sponge;

pub struct Prover<F: PrimeField + Absorb> {
    g: SparsePolynomial<F, SparseTerm>,
    verbose: bool,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

pub struct Verifier<F: PrimeField + Absorb> {
    g: SparsePolynomial<F, SparseTerm>,
    verbose: bool,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

#[derive(Clone, Debug)]
pub struct Transcript<F: PrimeField + Absorb + Absorb> {
    pub values: Vec<Vec<F>>,
    pub challenges: Vec<F>,
}

impl<F: PrimeField + Absorb> Transcript<F> {
    fn new() -> Self {
        Self {
            values: Vec::new(),
            challenges: Vec::new(),
        }
    }
    fn update(&mut self, poly: UnivariatePolynomial<F>, sponge: &mut PoseidonSponge<F>) {
        let coeffs = poly.coefficients().to_vec();

        for elem in coeffs.iter() {
            sponge.absorb(elem);
        }
        self.values.push(coeffs);

        let r = sponge.squeeze_field_elements::<F>(1)[0];
        self.challenges.push(r);
    }
    fn verify(&self, sponge: &mut PoseidonSponge<F>) -> bool {
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

impl<F: PrimeField + Absorb> Prover<F> {
    fn run(&mut self) -> Transcript<F> {
        let v = self.g.num_vars;
        let mut hyperpolys = vec![self.g.clone()];

        // sum g over each bit hypercube of dimension i = 1, 2, ... corresponding to the last i variables
        for i in (0..v).rev() {
            let last = hyperpolys.last().unwrap();
            hyperpolys.push(specialise(&last, i, F::zero()) + specialise(&last, i, F::one()));
        }

        // first round; h contains a constant univariate polynomial corresponding to the value the Prover wants to prove
        let h = UnivariatePolynomial::from_sparse_multivariate(&hyperpolys.pop().unwrap());

        self.transcript.update(h.unwrap(), &mut self.sponge);

        for _ in 1..=v {
            // r-th round
            let mut round_pol = hyperpolys.pop().unwrap();

            for (i, val) in self.transcript.challenges.iter().skip(1).enumerate() {
                round_pol = specialise(&round_pol, i, *val);
            }

            let uni_round_pol: UnivariatePolynomial<F> = UnivariatePolynomial::from_sparse_multivariate(&round_pol).unwrap();

            self.transcript.update(uni_round_pol, &mut self.sponge);
            
        }

        println!("Prover finished successfully");

        self.transcript.clone()
    }
}

impl<F: PrimeField + Absorb> Verifier<F> {
    fn run(&mut self) -> bool {
        self.transcript.verify(&mut self.sponge);

        let v = self.g.num_vars;

        let mut last_pol = UnivariatePolynomial::new(&self.transcript.values[0]);

        let mut last_sent_scalar = self.transcript.challenges[1];

        // round numbering in line with the book's (necessarily inconsistent variable indexing!)
        for r in 1..=v {
            let new_pol = UnivariatePolynomial::new(&self.transcript.values[r]);

            if new_pol.degree() > partial_degree(&self.g, r - 1) {
                println!("Verifier found inconsistent degrees in round {r}: the received unviariate polynomial {new_pol} has degree {}\n  \
                          whereas the original one has degree {} in variable x_{}. Aborting.",
                          new_pol.degree(), partial_degree(&self.g, r - 1), r - 1
                        );
                return false;
            }

            if new_pol.eval(F::zero()) + new_pol.eval(F::one()) != last_pol.eval(last_sent_scalar) {
                println!("Verifier found inconsistent evaluation in round {r}: received univariate polynomial is f_n = {new_pol},\n  \
                          previous one is f_o = {last_pol}. The equation f_n(0) + f_n(1) = f_o({last_sent_scalar}) fails to hold. Aborting."
                        );
                return false;
            }

            // cannot fail to unwrap due to round count
            last_sent_scalar = self.transcript.challenges[r];

            last_pol = new_pol;
        }

        if last_pol.eval(last_sent_scalar) != self.g.evaluate(&self.transcript.challenges[1..].to_vec()) {
            println!("Verifier found inconsistent evaluation in the oracle call: received univariate polynomial {last_pol} evaluates to {},\n  \
                      whereas original multivariate one evaluates to {} on {:?}. Aborting.",
                      last_pol.eval(last_sent_scalar), self.g.evaluate(&self.transcript.challenges), self.transcript.challenges
                    );
            return false;
        }

        println!("Verifier finished successfully and is confident in the value {}", self.transcript.values[0][0]);

        true
    }
}

// run the protocol and return true iff the verifier does not abort
pub fn run_sumcheck_protocol<F: PrimeField + Absorb>(
    pol: SparsePolynomial<F, SparseTerm>,
    verbose: bool,
) {
    println!("Initiated sumcheck protocol on the polynomial {pol:?}");

    let mut prover = Prover {
        g: pol.clone(),
        verbose: verbose,
        sponge: test_sponge(),
        transcript: Transcript::new(),
    };

    let transcript = prover.run();

    let mut verifier = Verifier {
        g: pol,
        verbose: verbose,
        sponge: test_sponge(),
        transcript,
    };

    assert!(verifier.run());

}
