use std::iter::zip;

use ark_bls12_381::g2;
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;

use ark_poly::{MultilinearExtension, SparseMultilinearExtension};

use crate::utils::{interpolate_uni_poly, test_sponge, usize_to_zxy, Transcript};

#[derive(Clone, Debug)]
pub struct Product<F: PrimeField, MLE: MultilinearExtension<F>>(
    pub SparseMultilinearExtension<F>,
    pub MLE,
    pub MLE,
);

#[derive(Clone, Debug)]
pub struct SumOfProducts<F: PrimeField, MLE: MultilinearExtension<F>> {
    // for add:
    // 1st: [alpha * add(g1, x, y) + beta * add(g2, x, y)] * f2(x) * 1 +
    // 2nd: [alpha * add(g1, x, y) + beta * add(g2, x, y)] * 1 * f3(y) +
    // 3rd [alpha * mul(g1, x, y) + beta * mul(g2, x, y)] * f2(x) * f3(y)
    pub terms: Vec<Product<F, MLE>>,
}

pub struct Prover<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> {
    sum_of_products: SumOfProducts<F, MLE>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

pub trait Oracle<F> {
    fn divinate(&self, x: &[F], y: &[F]) -> F;
    fn num_rounds(&self) -> usize;
    fn trust_message(&self) -> String;
}

struct PolyOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: SumOfProducts<F, MLE>,
    g: Vec<F>,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> PolyOracle<F, MLE> {
    fn new(sum_of_products: SumOfProducts<F, MLE>, g: Vec<F>) -> Self {
        Self { sum_of_products, g }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for PolyOracle<F, MLE> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {
        let mut zxy: Vec<F> = Vec::from(self.g.clone());
        zxy.extend(x);
        zxy.extend(y);

        let mut sum = F::zero();

        for Product(f1, f2, f3) in &self.sum_of_products.terms {
            sum += f1.evaluate(&zxy).unwrap() * f2.evaluate(&x).unwrap() * f3.evaluate(&y).unwrap()
        }

        sum
    }

    // returns the number of variables of EACH of the two polynomials
    // the total number of variables of the product is twice that
    fn num_rounds(&self) -> usize {
        let Product(_, f2, _) = &self.sum_of_products.terms[0];
        2 * f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}

struct CombinedPolyOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: SumOfProducts<F, MLE>,
    g1: Vec<F>,
    g2: Vec<F>,
    alpha: F,
    beta: F
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> CombinedPolyOracle<F, MLE> {
    fn new(sum_of_products: SumOfProducts<F, MLE>, g1: Vec<F>, g2: Vec<F>, alpha: F, beta: F) -> Self {
        Self { sum_of_products, g1, g2, alpha, beta }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for CombinedPolyOracle<F, MLE> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {

        let mut sums = Vec::new();

        for (coeff, g) in vec![(self.alpha, self.g1.clone()), (self.beta, self.g2.clone())] {
            let mut zxy: Vec<F> = Vec::from(g.clone());
            zxy.extend(x);
            zxy.extend(y);

            let mut sum = F::zero();

            for Product(f1, f2, f3) in &self.sum_of_products.terms {
                sum += f1.evaluate(&zxy).unwrap() * f2.evaluate(&x).unwrap() * f3.evaluate(&y).unwrap()
            }

            sums.push(coeff * sum);
        }

        sums[0] + sums[1]
    }

    // returns the number of variables of EACH of the two polynomials
    // the total number of variables of the product is twice that
    fn num_rounds(&self) -> usize {
        let Product(_, f2, _) = &self.sum_of_products.terms[0];
        2 * f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}

struct GKROracle<F: PrimeField> {
    f1: SparseMultilinearExtension<F>,
    g: Vec<F>,
    yb: F,
    yc: F,
    num_rounds: usize,
}

impl<F: PrimeField> GKROracle<F> {
    fn new(f1: SparseMultilinearExtension<F>, g: Vec<F>, yb: F, yc: F, num_rounds: usize) -> Self {
        Self {
            f1,
            g,
            yb,
            yc,
            num_rounds,
        }
    }
}

impl<F: PrimeField> Oracle<F> for GKROracle<F> {
    fn divinate(&self, x: &[F], y: &[F]) -> F {
        let mut zxy: Vec<F> = Vec::from(self.g.clone());
        zxy.extend(x);
        zxy.extend(y);

        self.f1.evaluate(&zxy).unwrap() * self.yb * self.yc
    }

    // returns the number of variables of EACH of the two polynomials
    // the total number of variables of the product is twice that
    fn num_rounds(&self) -> usize {
        self.num_rounds
    }

    fn trust_message(&self) -> String {
        format!(
            "assuming that f(b*) = {:?} and f(c*) = {:?}",
            self.yb, self.yc
        )
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
    pub fn new(
        sum_of_products: SumOfProducts<F, MLE>,
        sponge: PoseidonSponge<F>,
    ) -> Self {
        Self {
            sum_of_products,
            sponge,
            transcript: Transcript::new(),
        }
    }

    #[allow(non_snake_case)]
    fn sumcheck_prod(&mut self, A_fs: &mut Vec<Vec<F>>, A_gs: &mut Vec<Vec<F>>, v: usize) -> F {
        // first round; claimed_value contains the value the Prover wants to prove
        let num_summands = A_fs.len();
        assert_eq!(num_summands, A_gs.len());

        let claimed_value = (0..num_summands)
            .map(|summand| {
                (0..1 << v)
                    .map(|i| A_fs[summand][i] * A_gs[summand][i])
                    .sum::<F>()
            })
            .sum();

        let le_indices = to_le_indices(v);

        
        for i in 1..=v {
            let mut values = vec![F::zero(); 3];
            
            for summand in 0..num_summands {
                let mut f_values = vec![vec![F::zero(); 1 << (v - 1)]; 3];
                let mut g_values = vec![vec![F::zero(); 1 << (v - 1)]; 3];
                let A_f = &A_fs[summand];
                let A_g = &A_gs[summand];
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
                // evaluations at t = 0,1,2
                let temp_values: Vec<F> = (0..=2)
                    .map(|t| ((0..1 << (v - i)).map(|b| f_values[t][b] * g_values[t][b])).sum())
                    .collect();

                values[0] += temp_values[0];
                values[1] += temp_values[1];
                values[2] += temp_values[2];
            }

            let r_i = self.transcript.update(&values, &mut self.sponge);

            // Algorithm 1, part 2
            for summand in 0..num_summands {
                let A_f = &mut A_fs[summand];
                let A_g = &mut A_gs[summand];

                for b in 0..(1 << (v - i)) {
                    A_f[le_indices[b]] =
                        A_f[le_indices[b]] * (F::one() - r_i) + A_f[le_indices[b + (1 << (v - i))]] * r_i;
                    A_g[le_indices[b]] =
                        A_g[le_indices[b]] * (F::one() - r_i) + A_g[le_indices[b + (1 << (v - i))]] * r_i;
                }
            }
        }

        claimed_value
    }

    // Algorithm 6 (Sumcheck GKR) adapted to a sum of products
    #[allow(non_snake_case)]
    pub fn run(&mut self, g1: &[F], g2: &[F], alpha: F, beta: F) -> Transcript<F> {
        let Product(_, pol, _) = &self.sum_of_products.terms[0];
        let k: usize = pol.num_vars();

        let mut A_hs: Vec<Vec<F>> = Vec::new();
        let mut A_f2s: Vec<Vec<F>> = Vec::new();
        
        // phase 1
        for Product(f1, f2, f3) in self.sum_of_products.terms.iter() {
            // TODO can we share precomputation for the first two products?
            A_hs.push(initialise_phase_1(&f1, f3, g1, g2, alpha, beta));
            A_f2s.push(f2.to_evaluations());
        }
        
        let out = self.sumcheck_prod(&mut A_hs, &mut A_f2s, k);
        
        // the final claim comes from the first run, since the second is summing over a "different" `f1` (with `x` fixed to some random `u`)
        self.transcript.set_claim(out);
        
        // phase 2
        let u = &self.transcript.challenges[..];
        
        let mut A_f3s: Vec<Vec<F>> = Vec::new();
        let mut A_f1s: Vec<Vec<F>> = Vec::new();
        // let f2s_u: Vec<F> = Vec::new();
        let mut A_f3_f2_us: Vec<Vec<F>> = Vec::new();
        for Product(f1, f2, f3) in &self.sum_of_products.terms {
            let f2_u = f2.evaluate(&u).unwrap();
            let A_f3 = f3.to_evaluations();
            
            A_f1s.push(initialise_phase_2::<F>(&f1, g1, g2, &u, alpha, beta));
            A_f3_f2_us.push(A_f3.iter().map(|x| *x * f2_u).collect::<Vec<F>>());
        }

        self.sumcheck_prod(&mut A_f1s, &mut A_f3_f2_us, k);

        println!("Prover finished successfully");

        self.transcript.clone()
    }
}

impl<F: PrimeField + Absorb, O: Oracle<F>> Verifier<F, O> {
    fn run(&mut self) -> bool {
        self.transcript.verify(&mut self.sponge);

        let v = self.oracle.num_rounds();

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
            "Verifier finished successfully and is confident in the value {:?} {}",
            self.transcript.claim,
            self.oracle.trust_message()
        );

        true
    }
}

pub(crate) fn old_precompute<F: PrimeField>(vals: &[F]) -> Vec<F> {
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

pub(crate) fn precompute<F: PrimeField>(vals: &[F]) -> Vec<F> {
    let n = vals.len();

    let mut table_1 = vec![F::zero(); 1 << n];
    let mut table_2: Vec<F> = vec![F::zero(); 1 << n];

    let mut old_table = &mut table_1;
    let mut new_table = &mut table_2;

    old_table[0] = F::one();

    for i in 0..n {
        for b in 0..(1 << i) {
            new_table[2 * b] = old_table[b] * (F::one() - vals[i]);
            new_table[2 * b + 1] = old_table[b] * vals[i];
        }
        std::mem::swap(&mut old_table, &mut new_table);
    }

    if n % 2 == 0 {
        table_1
    } else {
        table_2
    }
}

pub(crate) fn initialise_phase_1<F: PrimeField, MLE: MultilinearExtension<F>>(
    f_1: &SparseMultilinearExtension<F>,
    f_3: &MLE,
    g1: &[F],
    g2: &[F],
    alpha: F,
    beta: F,
) -> Vec<F> {
    let v = f_3.num_vars();
    let le_indices_f1 = to_le_indices(f_1.num_vars);
    let le_indices_f3 = to_le_indices(v);

    let table_g1 = precompute(g1);
    let table_g2 = precompute(g2);
    let table_g = table_g1
        .iter()
        .zip(table_g2.iter())
        .map(|(x, y)| alpha * *x + beta * *y)
        .collect::<Vec<F>>();

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
    g1: &[F],
    g2: &[F],
    u: &[F],
    alpha: F,
    beta: F,
) -> Vec<F> {
    let v = g1.len();
    let le_indices_f_1 = to_le_indices(f_1.num_vars);
    let le_indices_g: Vec<usize> = to_le_indices(v);

    let table_g1 = precompute(g1);
    let table_g2 = precompute(g2);
    let table_u = precompute(u);

    let mut af1 = vec![F::zero(); 1 << v];

    for (idx_le, val) in f_1.evaluations.iter() {
        let idx = le_indices_f_1[*idx_le];
        let (z, x, y) = usize_to_zxy(idx, v);
        af1[le_indices_g[y]] += (alpha * table_g1[z] + beta * table_g2[z]) * table_u[x] * val;
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

    let simple_sum = SumOfProducts {
        terms: vec![Product(f1, f2, f3)],
    };
    
    let mut prover = Prover::new(simple_sum.clone(), test_sponge());

    let transcript = prover.run(&g, &g, F::one(), F::zero());

    let mut verifier = Verifier {
        oracle: PolyOracle::new(simple_sum, g.to_vec()), // TODO: decide if passing & is enough
        sponge: test_sponge(),
        transcript,
    };

    assert!(verifier.run());
}


// run the protocol and return true iff the verifier does not abort
pub fn run_sumcheck_protocol_combined<F: PrimeField + Absorb, MLE: MultilinearExtension<F>>(
    f1: SparseMultilinearExtension<F>,
    f2: MLE,
    f3: MLE,
    g1: &[F],
    g2: &[F],
    alpha: F,
    beta: F
) {

    let simple_sum = SumOfProducts {
        terms: vec![Product(f1, f2, f3)],
    };
    
    let mut prover = Prover::new(simple_sum.clone(), test_sponge());

    let transcript = prover.run(&g1, &g2, alpha, beta);

    let mut verifier = Verifier {
        oracle: CombinedPolyOracle::new(
            simple_sum,
            g1.to_vec(),
            g2.to_vec(),
            alpha,
            beta
        ),
        sponge: test_sponge(),
        transcript,
    };

    assert!(verifier.run());
}

// run the protocol and return true iff the verifier does not abort
pub fn run_sumcheck_protocol_combined_multiprod<F: PrimeField + Absorb, MLE: MultilinearExtension<F>>(
    sum_of_products: SumOfProducts<F, MLE>,
    g1: &[F],
    g2: &[F],
    alpha: F,
    beta: F
) {
    let mut prover = Prover::new(sum_of_products.clone(), test_sponge());

    let transcript = prover.run(&g1, &g2, alpha, beta);

    let mut verifier = Verifier {
        oracle: CombinedPolyOracle::new(
            sum_of_products,
            g1.to_vec(),
            g2.to_vec(),
            alpha,
            beta
        ),
        sponge: test_sponge(),
        transcript,
    };

    assert!(verifier.run());
}