
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

use ark_ff::Field;

use ark_poly::{polynomial::multivariate::{SparsePolynomial, SparseTerm}, Polynomial};

use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use crate::polynomial::{UnivariatePolynomial, specialise, partial_degree};

pub struct Prover<F: Field> {
    tx: Sender<UnivariatePolynomial<F>>,
    rx: Receiver<F>,
    g: SparsePolynomial<F, SparseTerm>,
    verbose: bool,
}

pub struct Verifier<F: Field> {
    tx: Sender<F>,
    rx: Receiver<UnivariatePolynomial<F>>,
    g: SparsePolynomial<F, SparseTerm>,
    verbose: bool,
}

impl<F: Field> Prover<F> {
    fn run(&self) {

        let v = self.g.num_vars;
        let mut hyperpolys = vec![self.g.clone()];

        // sum g over each bit hypercube of dimension i = 1, 2, ... corresponding to the last i variables
        for i in (0..v).rev() {
            let last = hyperpolys.last().unwrap();
            hyperpolys.push(specialise(&last, i, F::zero()) + specialise(&last, i, F::one()));
        }

        let mut ver_values = Vec::new();

        // first round
        let h = UnivariatePolynomial::from_sparse_multivariate(&hyperpolys.pop().unwrap());

        // h contains a constant univariate polynomial corresponding to the value the Prover wants to prove
        if self.tx.send(h.unwrap()).is_err() {
            println!("Prover failed to send claimed value. Aborting.");
            return
        }

        for r in 1..v + 1 {
            // r-th round
            let mut round_pol = hyperpolys.pop().unwrap();
    
            for (i, v) in ver_values.iter().enumerate() {
                round_pol = specialise(&round_pol, i, *v);
            }

            let uni_round_pol = UnivariatePolynomial::from_sparse_multivariate(&round_pol).unwrap();

            if self.tx.send(uni_round_pol).is_err() {
                println!("Prover failed to send polynomial (the channel might have been closed). Aborting.");
                return
            }

            if r < v {
                let ri = match self.rx.recv() {
                    Err(e) => {
                        println!("Prover failed to receive value: {}. Aborting.", e);
                        return
                    },
                    Ok(val) => {
                        if self.verbose {
                            println!(" - Prover round {}: received value {}", r, val);
                        }
                        val
                    },
                };

                ver_values.push(ri);
            }

        }

        println!("Prover finished successfully");

    }
}

impl<F: Field> Verifier<F> {
    fn run(&self) -> bool {

        let v = self.g.num_vars;
        let rng = &mut ChaCha20Rng::from_entropy();

        let mut scalars = Vec::new();

        for _ in 0..v {
            scalars.push(F::rand(rng));
        } 

        // MODIFY
        // in order to simulate the example on page 36 of Thaler's book, uncomment
        /* let two = F::one() + F::one(); // the Field trait does not have a F::from(2) method, unlike some of its implementors
        let three = two + F::one();
        let six = three + three;
        let scalars = vec![two, three, six]; */

        let mut scalars_iter = scalars.iter();
    
        let mut last_pol = match self.rx.recv() {
            Err(e) => {
                println!("Verifier failed to receive value: {}. Aborting.", e);
                return false;
            },
            Ok(pol) => pol,
        };

        // pol is constant here (even if the Prover cheats and sends a non-constant polynomial, this step is still sound)
        let claimed_val = last_pol.eval(F::zero());

        let mut last_sent_scalar = F::zero(); // dummy value for the first round (the evaluated polynomial is constant)

        // round numbering in line with the book's (necessarily inconsistent variable indexing!)
        for r in 1..v + 1 {

            let new_pol = match self.rx.recv() {
                Err(_) => {
                    println!("Verifier failed to receive polynomial in round {}. Aborting.", r);
                    return false;
                },
                Ok(pol) => {
                    if self.verbose {
                        println!(" - Verifier round {}: received polynomial {}", r, pol);
                    }
                    pol
                },
            };

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
            last_sent_scalar = *scalars_iter.next().unwrap();
            
            if r < v {
                if self.tx.send(last_sent_scalar).is_err() {
                    println!("Verifier failed to send value (the channel might have been closed). Aborting.");
                    return false;
                }
            }

            last_pol = new_pol;

        }

        if last_pol.eval(last_sent_scalar) != self.g.evaluate(&scalars) {
            println!("Verifier found inconsistent evaluation in the oracle call: received univariate polynomial {last_pol} evaluates to {},\n  \
                      whereas original multivariate one evaluates to {} on {:?}. Aborting.",
                      last_pol.eval(last_sent_scalar), self.g.evaluate(&scalars), scalars
                    );
            return false;
        }

        println!("Verifier finished successfully and is confident in the value {claimed_val}");

        true

    }
}

// run the protocol and return true iff the verifier does not abort
pub fn run_sumcheck_protocol<F: Field>(pol: SparsePolynomial<F, SparseTerm>, verbose: bool) -> bool {

    println!("Initiated sumcheck protocol on the polynomial {pol:?}");

    let (tx_pol, rx_pol) = mpsc::channel::<UnivariatePolynomial<F>>();
    let (tx_val, rx_val) = mpsc::channel::<F>();

    let prover = Prover {
        tx: tx_pol,
        rx: rx_val,
        g: pol.clone(),
        verbose: verbose,
    };

    let verifier = Verifier {
        tx: tx_val,
        rx: rx_pol,
        g: pol,
        verbose: verbose,
    };

    let thread_prover = thread::spawn(move || prover.run());

    let thread_verifier = thread::spawn(move || verifier.run());

    thread_prover.join();

    thread_verifier.join().unwrap()

}
