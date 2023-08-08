use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    Absorb, CryptographicSponge,
};
use ark_ff::{BigInteger, PrimeField};
use halo2_base::utils::ScalarField;
use poseidon_native::Poseidon;

pub(crate) fn test_sponge<F: ScalarField>() -> Poseidon<F, 3, 2> {
    Poseidon::<F, 3, 2>::new(8, 57)
}

/// `SumcheckProof` can be safely shared with the verifier in full.
#[derive(Clone, Debug, Default)]
pub struct SumcheckProof<F: PrimeField> {
    pub values: Vec<Vec<F>>,
    pub claim: Option<F>,
}

/// Each party will instantiate a `Transcript`.
/// The prover will feed into the transcript the values of the polynomial, thus creating the `SumcheckProof`.
/// The verifier will feed into the transcript the values received in `SumcheckProof`.
#[derive(Clone)]
pub struct Transcript<F: PrimeField + Absorb, F2: ScalarField> {
    pub proof: SumcheckProof<F>,
    pub sponge: Poseidon<F2, 3, 2>,
}

impl<F: PrimeField + Absorb, F2: ScalarField> Transcript<F, F2> {
    pub fn new() -> Self {
        Self {
            proof: SumcheckProof::default(),
            sponge: test_sponge(),
        }
    }

    pub fn new_from_proof(proof: SumcheckProof<F>) -> Self {
        Self {
            proof,
            sponge: test_sponge(),
        }
    }

    pub fn update(&mut self, elements: &[F]) -> F {
        for elem in elements.iter() {
            let bytes = BigInteger::to_bytes_le(&elem.into_bigint());
            let halo2_base_elem = F2::from_bytes_le(&bytes);
            self.sponge.update(&[halo2_base_elem]);
        }
        self.proof.values.push(elements.to_vec());

        let halo2_out = self.sponge.squeeze();
        let out = F::from_le_bytes_mod_order(&halo2_out.to_bytes_le());
        out
    }

    pub fn set_claim(&mut self, claimed_value: F) {
        self.proof.claim = Some(claimed_value);
    }
}

// separate a bit string represented as a usize into three values,
// the last of which correspond to
#[inline]
pub fn usize_to_zxy(idx: usize, v: usize) -> (usize, usize, usize) {
    let vp = 1 << v;

    let y = idx % vp;
    let idx_zx = idx / vp;
    let x = idx_zx % vp;
    let z = idx_zx / vp;

    (z, x, y)
}

/// interpolate the *unique* univariate polynomial of degree *at most*
/// p_i.len()-1 passing through the y-values in p_i at x = 0,..., p_i.len()-1
/// and evaluate this  polynomial at `eval_at`. In other words, efficiently compute
///  \sum_{i=0}^{len p_i - 1} p_i[i] * (\prod_{j!=i} (eval_at - j)/(i-j))
pub(crate) fn interpolate_uni_poly<F: PrimeField>(p_i: &[F], eval_at: F) -> F {
    let len = p_i.len();

    let mut evals = vec![];

    let mut prod = eval_at;
    evals.push(eval_at);

    //`prod = \prod_{j} (eval_at - j)`
    // we return early if 0 <= eval_at <  len, i.e. if the desired value has been passed
    let mut check = F::zero();
    for i in 1..len {
        if eval_at == check {
            return p_i[i - 1];
        }
        check += F::one();

        let tmp = eval_at - check;
        evals.push(tmp);
        prod *= tmp;
    }

    if eval_at == check {
        return p_i[len - 1];
    }

    let mut res = F::zero();
    // we want to compute \prod (j!=i) (i-j) for a given i
    //
    // we start from the last step, which is
    //  denom[len-1] = (len-1) * (len-2) *... * 2 * 1
    // the step before that is
    //  denom[len-2] = (len-2) * (len-3) * ... * 2 * 1 * -1
    // and the step before that is
    //  denom[len-3] = (len-3) * (len-4) * ... * 2 * 1 * -1 * -2
    //
    // i.e., for any i, the one before this will be derived from
    //  denom[i-1] = - denom[i] * (len-i) / i
    //
    // that is, we only need to store
    // - the last denom for i = len-1, and
    // - the ratio between the current step and the last step, which is the
    //   product of -(len-i) / i from all previous steps and we store
    //   this product as a fraction number to reduce field divisions.

    // We know
    //  - 2^61 < factorial(20) < 2^62
    //  - 2^122 < factorial(33) < 2^123
    // so we will be able to compute the ratio
    //  - for len <= 20 with i64
    //  - for len <= 33 with i128
    //  - for len >  33 with BigInt
    if p_i.len() <= 20 {
        let last_denom = F::from(u64_factorial(len - 1));
        let mut ratio_numerator = 1i64;
        let mut ratio_enumerator = 1u64;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                -F::from((-ratio_numerator) as u64)
            } else {
                F::from(ratio_numerator as u64)
            };

            res += p_i[i] * prod * F::from(ratio_enumerator)
                / (last_denom * ratio_numerator_f * evals[i]);

            // compute ratio for the next step which is current_ratio * -(len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i64 - i as i64);
                ratio_enumerator *= i as u64;
            }
        }
    } else if p_i.len() <= 33 {
        let last_denom = F::from(u128_factorial(len - 1));
        let mut ratio_numerator = 1i128;
        let mut ratio_enumerator = 1u128;

        for i in (0..len).rev() {
            let ratio_numerator_f = if ratio_numerator < 0 {
                -F::from((-ratio_numerator) as u128)
            } else {
                F::from(ratio_numerator as u128)
            };

            res += p_i[i] * prod * F::from(ratio_enumerator)
                / (last_denom * ratio_numerator_f * evals[i]);

            // compute ratio for the next step which is current_ratio * -(len-i)/i
            if i != 0 {
                ratio_numerator *= -(len as i128 - i as i128);
                ratio_enumerator *= i as u128;
            }
        }
    } else {
        // since we are using field operations, we can merge
        // `last_denom` and `ratio_numerator` into a single field element.
        let mut denom_up = field_factorial::<F>(len - 1);
        let mut denom_down = F::one();

        for i in (0..len).rev() {
            res += p_i[i] * prod * denom_down / (denom_up * evals[i]);

            // compute denom for the next step is -current_denom * (len-i)/i
            if i != 0 {
                denom_up *= -F::from((len - i) as u64);
                denom_down *= F::from(i as u64);
            }
        }
    }

    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u128_factorial(a: usize) -> u128 {
    let mut res = 1u128;
    for i in 1..=a {
        res *= i as u128;
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn u64_factorial(a: usize) -> u64 {
    let mut res = 1u64;
    for i in 1..=a {
        res *= i as u64;
    }
    res
}

/// compute the factorial(a) = 1 * 2 * ... * a
#[inline]
fn field_factorial<F: PrimeField>(a: usize) -> F {
    let mut res = F::one();
    for i in 1..=a {
        res *= F::from(i as u64);
    }
    res
}

// test for conversions between ScalarField and PrimeField
#[cfg(test)]
mod tests {
    use ark_ff::BigInteger;

    #[test]
    fn test_conversion_ark_to_halo() {
        use ark_bn254::Fr;
        use ark_ff::PrimeField;

        use halo2_base::halo2_proofs::halo2curves::bn256::Fr as Fr2;
        use halo2_base::utils::ScalarField;

        use std::ops::Neg;

        for i in 0u64..10 {
            let elem = Fr::from(i);
            let elem2 = Fr2::from(i);
            let bytes_elem = BigInteger::to_bytes_le(&elem.into_bigint());
            let bytes_elem2 = elem2.to_bytes_le();
            assert_eq!(bytes_elem, bytes_elem2);
        }

        for i in 0u64..10 {
            let elem = Fr::from(i).neg();
            let elem2 = -Fr2::from(i);
            let bytes_elem = BigInteger::to_bytes_le(&elem.into_bigint());
            let bytes_elem2 = elem2.to_bytes_le();
            assert_eq!(bytes_elem, bytes_elem2);
        }
    }
}
