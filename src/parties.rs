use std::{marker::PhantomData, thread::sleep_ms};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;

use ark_poly::{MultilinearExtension, SparseMultilinearExtension};
use rand_chacha::rand_core::le;

use crate::utils::{interpolate_uni_poly, test_sponge, usize_to_zxy, Transcript};

pub struct Prover<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> {
    f1: SparseMultilinearExtension<F>,
    f2: MLE,
    f3: MLE,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

pub trait Oracle<F> {
    fn divinate(&self, x: &[F], y: &[F]) -> F;
    fn num_vars(&self) -> usize;
}

struct PolyOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    f1: SparseMultilinearExtension<F>,
    f2: MLE,
    f3: MLE,
    g: Vec<F>,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> PolyOracle<F, MLE> {
    fn new(f1: SparseMultilinearExtension<F>, f2: MLE, f3: MLE, g: Vec<F>) -> Self {
        Self { f1, f2, f3, g }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for PolyOracle<F, MLE> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {
        let mut zxy: Vec<F> = Vec::from(self.g.clone());
        zxy.extend(x);
        zxy.extend(y);

        self.f1.evaluate(&zxy).unwrap()
            * self.f2.evaluate(&x).unwrap()
            * self.f3.evaluate(&y).unwrap()
    }

    // returns the number of variables of EACH of the two polynomials
    // the total number of variables of the product is twice that
    fn num_vars(&self) -> usize {
        2 * self.f2.num_vars()
    }
}

pub struct Verifier<F: PrimeField + Absorb, O: Oracle<F>> {
    oracle: O,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

pub(crate) fn to_le_indices(v: usize) -> Vec<usize> {
    // preparing index conversion big-endian -> little-endian
    // we do this because `DenseMultilinearExtension::from_evaluations_slice` expects little-endian notation for the evaluation slice
    (0usize..(1 << v))
        .into_iter()
        .map(|i| i.reverse_bits() >> (usize::BITS as usize - v))
        .collect()
}

impl<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> Prover<F, MLE> {
    fn new(f1: SparseMultilinearExtension<F>, f2: MLE, f3: MLE, sponge: PoseidonSponge<F>) -> Self {
        // let v = f2.num_vars();
        // assert_eq!(
        //     v,
        //     g.num_vars(),
        //     "number of variables mismatched: {v} vs {}",
        //     g.num_vars()
        // );

        Self {
            f1,
            f2,
            f3,
            sponge,
            transcript: Transcript::new(),
        }
    }
    fn sumcheck_prod(&mut self, A_f: &mut Vec<F>, A_g: &mut Vec<F>, v: usize) -> F {
        // first round; claimed_value contains the value the Prover wants to prove
        let claimed_value = (0..1 << v).map(|i| A_f[i] * A_g[i]).sum();

        // self.transcript.update(&[claimed_value], &mut self.sponge);

        let mut f_values = vec![vec![F::zero(); 1 << (v - 1)]; 3];
        let mut g_values = vec![vec![F::zero(); 1 << (v - 1)]; 3];

        let le_indices = to_le_indices(v);
        for i in 1..=v {
            // i-th round

            // Algorithm 1
            for b in 0..(1 << (v - i)) {
                f_values[0][b] = A_f[le_indices[b]];
                f_values[1][b] = A_f[le_indices[b + (1 << (v - i))]];
                f_values[2][b] =
                    -A_f[le_indices[b]] + A_f[le_indices[b + (1 << (v - i))]] * F::from(2u64);

                g_values[0][b] = A_g[le_indices[b]];
                g_values[1][b] = A_g[le_indices[b + (1 << (v - i))]];
                g_values[2][b] =
                    -A_g[le_indices[b]] + A_g[le_indices[b + (1 << (v - i))]] * F::from(2u64);
            }

            // Algorithm 3
            let values: Vec<F> = (0..=2)
                .map(|t| ((0..1 << (v - i)).map(|b| f_values[t][b] * g_values[t][b])).sum())
                .collect();

            let r_i = self.transcript.update(&values, &mut self.sponge);

            // Algorithm 1, part 2
            for b in 0..(1 << (v - i)) {
                A_f[le_indices[b]] = A_f[le_indices[b]] * (F::one() - r_i)
                    + A_f[le_indices[b + (1 << (v - i))]] * r_i;
                A_g[le_indices[b]] = A_g[le_indices[b]] * (F::one() - r_i)
                    + A_g[le_indices[b + (1 << (v - i))]] * r_i;
            }
        }

        claimed_value
    }

    fn run(&mut self, g: &[F]) -> Transcript<F> {
        // Algorithm 6. Sumcheck GKR
        let mut A_h = initialise_phase_1(&self.f1, &self.f3, g);

        // phase 1
        let mut A_f2 = self.f2.to_evaluations();
        let A_f3 = self.f3.to_evaluations();

        let out = self.sumcheck_prod(&mut A_h, &mut A_f2, self.f2.num_vars());

        // the final claim comes from the first run, since the second is summing over a "different" `f1` (with `x` fixed to some random `u`)
        self.transcript.set_claim(out);

        // phase 2
        let u = &self.transcript.challenges[..];

        let f2_u = self.f2.evaluate(&u).unwrap();

        let mut A_f1 = initialise_phase_2::<F>(&self.f1, g, &u);
        let mut A_f3_f2_u = A_f3.iter().map(|x| *x * f2_u).collect::<Vec<F>>();

        self.sumcheck_prod(&mut A_f1, &mut A_f3_f2_u, self.f2.num_vars());

        println!("Prover finished successfully");

        self.transcript.clone()
    }
}

impl<F: PrimeField + Absorb, O: Oracle<F>> Verifier<F, O> {
    fn run(&mut self) -> bool {
        self.transcript.verify(&mut self.sponge);

        let v = self.oracle.num_vars();

        // last_pol is initialised to a singleton list containing the claimed value
        // in subsequent rounds, it contains the values of the would-be-sent polynomial at 0, 1 and 2
        let mut last_pol = vec![self.transcript.claim.unwrap()];

        // dummy scalar to be used when evaluating the constant poly in the first iteration
        let mut last_sent_scalar = F::one();

        // round numbering in line with Thaler's book (necessarily inconsistent variable indexing!)
        for r in 0..v {
            let new_pol_evals = self.transcript.values[r].clone();
            let claimed_sum = new_pol_evals[0] + new_pol_evals[1];
            let new_pol = Vec::from(&new_pol_evals[..]);

            if claimed_sum != interpolate_uni_poly(&last_pol, last_sent_scalar) {
                println!("Verifier found inconsistent evaluation in round {r}: received univariate polynomial is f_n = {new_pol:?},\n  \
                          previous one is f_o = {last_pol:?}. The equation f_n(0) + f_n(1) = f_o({last_sent_scalar}) fails to hold. Aborting."
                        );
                return false;
            }
            last_sent_scalar = self.transcript.challenges[r];

            last_pol = new_pol;
        }

        if interpolate_uni_poly(&last_pol, last_sent_scalar)
            != self.oracle.divinate(
                &self.transcript.challenges[..(v / 2)],
                &self.transcript.challenges[(v / 2)..],
            )
        {
            println!(
                "Verifier found inconsistent evaluation in the oracle call: \
                        received univariate polynomial {last_pol:?} evaluates to {},\n  \
                        whereas original multivariate one evaluates to {} on {:?}. Aborting.",
                interpolate_uni_poly(&last_pol, last_sent_scalar),
                self.oracle.divinate(
                    &self.transcript.challenges[..(v / 2)],
                    &self.transcript.challenges[(v / 2)..],
                ),
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

pub(crate) fn precompute<F: PrimeField>(vals: &[F]) -> Vec<F> {
    // TODO since we're cloning the new table, maybe the table sizes can be optimized
    // or keep 2 fixed tables and swap pointers in each iteration
    let n = vals.len();

    let mut table = vec![F::zero(); 1 << n];
    let mut new_table: Vec<F> = vec![F::zero(); 1 << n];

    table[0] = F::one();

    for i in 0..n {
        for b in 0..(1 << i) {
            new_table[2 * b] = table[b] * (F::one() - vals[i]);
            new_table[2 * b + 1] = table[b] * vals[i];
        }
        table = new_table.clone();
    }

    new_table
}

pub(crate) fn initialise_phase_1<F: PrimeField, MLE: MultilinearExtension<F>>(
    f_1: &SparseMultilinearExtension<F>,
    f_3: &MLE,
    g: &[F],
) -> Vec<F> {
    let v = f_3.num_vars();
    let le_indices_f1 = to_le_indices(f_1.num_vars);
    let le_indices_f3 = to_le_indices(v);

    let table_g = precompute(g);
    let mut ahg = vec![F::zero(); 1 << v];

    for (idx_le, val) in f_1.evaluations.iter() {
        let idx = le_indices_f1[*idx_le];
        let (z, x, y) = usize_to_zxy(idx, v);
        ahg[le_indices_f3[x]] += table_g[z] * val * f_3.to_evaluations()[le_indices_f3[y]];
    }

    ahg
}

pub(crate) fn initialise_phase_2<F: PrimeField>(
    f_1: &SparseMultilinearExtension<F>,
    g: &[F],
    u: &[F],
) -> Vec<F> {
    let v = g.len();
    let le_indices_f_1 = to_le_indices(f_1.num_vars);
    let le_indices_g = to_le_indices(g.len());

    let table_g = precompute(g);
    let table_u = precompute(u);

    let mut af1 = vec![F::zero(); 1 << v];

    for (idx_le, val) in f_1.evaluations.iter() {
        let idx = le_indices_f_1[*idx_le];
        let (z, x, y) = usize_to_zxy(idx, v);
        af1[le_indices_g[y]] += table_g[z] * table_u[x] * val;
    }

    af1
}

// run the protocol and return true iff the verifier does not abort
pub fn run_sumcheck_protocol<F: PrimeField + Absorb, MLE: MultilinearExtension<F>>(
    f1: SparseMultilinearExtension<F>,
    f2: MLE,
    f3: MLE,
    g: &[F],
) {
    // println!(
    //     "Initiated sumcheck protocol on the product of the polynomials \n\t{f:?}\nand\n\t{g:?}"
    // );

    let mut prover = Prover::new(f1.clone(), f2.clone(), f3.clone(), test_sponge());
    // let g = vec![F::one(); 1 << f2.num_vars()];

    let transcript = prover.run(&g);

    println!("CLAIMED VALUE {:?}", transcript.claim);

    let mut verifier = Verifier {
        oracle: PolyOracle::new(f1, f2, f3, g.to_vec()), // TODO: decide if passing & is enough
        sponge: test_sponge(),
        transcript,
    };

    assert!(verifier.run());
}
