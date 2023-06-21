use std::{fmt, iter::zip};

use ark_ff::Field;

use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};

const SIZE_USIZE: usize = 8 * std::mem::size_of::<usize>();

// raise the element to a given power using the square-and-multiply method
pub fn pow_sq_m<F: Field>(base: &F, pow: usize) -> F {
    let mut pow_two = base.to_owned();

    if base.is_zero() || base.is_one() || pow == 1 {
        return pow_two;
    }

    let mut acc = F::one();

    for b in (0..(SIZE_USIZE as u32) - pow.leading_zeros()).map(|n| (pow >> n) & 1) {
        if b == 1 {
            acc *= pow_two;
        }

        pow_two.square_in_place(); // does one trailing squaring more than necessary
    }

    acc
}

// specialise the variable with index i to the value val in pol
// returns a (new, owned) polynomial with one fewer variable than the original one
// if the index is not smaller than the number of variables, do nothing
pub fn specialise<F: Field>(
    pol: &SparsePolynomial<F, SparseTerm>,
    i: usize,
    val: F,
) -> SparsePolynomial<F, SparseTerm> {
    if i >= pol.num_vars {
        return pol.clone();
    }

    let mut monomials = Vec::new();

    for (coeff, term) in &pol.terms {
        let mut spec_term = Vec::new();
        let mut new_c = F::from(*coeff);

        for (id, pow) in zip(term.vars(), term.powers()) {
            if id == i {
                new_c *= pow_sq_m(&val, pow);
            } else {
                spec_term.push((id, pow));
            }
        }

        monomials.push((new_c, SparseTerm::new(spec_term)));
    }

    SparsePolynomial::from_coefficients_vec(pol.num_vars, monomials)
}

// if exactly one variable appears with non-zero exponents, return Ok(Some(its index));
// if constant, return Ok(None)
// otherwise, return Err(())
pub fn univariate_index<F: Field>(
    pol: &SparsePolynomial<F, SparseTerm>,
) -> Result<Option<usize>, ()> {
    let mut found = None;

    for (_, term) in &pol.terms {
        for (id, pow) in zip(term.vars(), term.powers()) {
            // checking non-zero might be unncessary by construction of SparsePolynomial
            if pow != 0 {
                match found {
                    None => found = Some(id),
                    Some(i) => {
                        if id != i {
                            return Err(());
                        }
                    }
                }
            }
        }
    }

    if found.is_none() {
        Ok(None)
    } else {
        Ok(found)
    }
}

// compute the degree of the given polynomial in its i-th variable
pub fn partial_degree<F: Field>(pol: &SparsePolynomial<F, SparseTerm>, i: usize) -> usize {
    if i >= pol.num_vars {
        return 0;
    }

    pol.terms
        .iter()
        .map(
            |(_, term)| {
                zip(term.vars(), term.powers())
                    .map(|(id, exp)| if id == i { exp } else { 0 })
                    .sum()
            }, // the constructor of SparsePolynomial already groups tuples corresponding to the same variable, but we add just to be safe
        )
        .max()
        .unwrap()
}

pub struct UnivariatePolynomial<F: Field>(Vec<F>);

impl<F: Field> UnivariatePolynomial<F> {
    // coefficients in *increasing* degree order. zeros are dropped.
    pub fn new(coeffs: &[F]) -> Self {
        if coeffs.len() == 0 {
            return Self(vec![F::zero()]);
        }

        let mut i = coeffs.len() - 1;

        while i > 0 && coeffs[i] == F::zero() {
            i -= 1;
        }

        Self(coeffs[..i + 1].to_owned())
        // TODO just change to Self(coeffs.to_owned())
    }

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn coefficients(&self) -> &[F] {
        &self.0
    }

    pub fn eval(&self, val: F) -> F {
        let mut res = F::zero();
        let mut pows = F::one();

        for c in self.0.iter() {
            res += pows * c;
            pows *= val;
        }

        res
    }

    pub fn from_sparse_multivariate(pol: &SparsePolynomial<F, SparseTerm>) -> Option<Self> {
        match univariate_index(&pol) {
            Err(()) => None,
            Ok(None) => Some(Self::new(&[pol.terms.first().unwrap().0])),
            Ok(Some(i)) => {
                let mut coeffs = vec![F::zero(); partial_degree(&pol, i) + 1];

                for (c, term) in pol.terms.iter() {
                    if term.is_constant() {
                        coeffs[0] += c;
                    } else {
                        for (id, pow) in zip(term.vars(), term.powers()) {
                            if id == i {
                                coeffs[pow] += c;
                            }
                        }
                    }
                }

                Some(Self::new(&coeffs))
            }
        }
    }
    fn interpolate(coords: &[(F, F)]) -> Self {
        
        //let linear = map(x)
        //for i in 
        unimplemented!();
    }
}

impl<F: Field> fmt::Display for UnivariatePolynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut o = 0;

        while o < self.0.len() && self.0[o] == F::zero() {
            o += 1;
        }

        for (d, c) in self.0[o + 1..].iter().enumerate().rev() {
            if *c != F::zero() {
                let rd = d + o + 1;
                if rd > 1 {
                    write!(f, "{c} * x^{} + ", rd)?
                } else {
                    // d = 1, as it can never be 0 in this loop
                    write!(f, "{} * x + ", c)?
                }
            }
        }

        let c = self.0[o];

        if o > 1 {
            write!(f, "{} * x^{}", c, o)
        } else if o == 1 {
            write!(f, "{} * x", c)
        } else {
            write!(f, "{}", c)
        }
    }
}
