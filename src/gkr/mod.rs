use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};
use std::marker::PhantomData;

use ark_ff::PrimeField;

use crate::parties::to_le_indices;

pub mod parties;
mod tests;

#[derive(Debug)]
pub struct Wiring<const D: usize> {
    curr_index: usize,
    left: usize,
    right: usize,
}

impl<const D: usize> Wiring<D> {
    fn new(curr_index: usize, left: usize, right: usize) -> Self {
        assert!(curr_index < 1 << D);
        assert!(left < 1 << D);
        assert!(right < 1 << D);
        Self {
            curr_index,
            left,
            right,
        }
    }
}

#[derive(Debug)]
pub struct Layer<const D: usize> {
    // pub num_vars: usize,
    pub add: Vec<Wiring<D>>,
    pub mul: Vec<Wiring<D>>,
}

impl<const D: usize> Layer<D> {
    fn new(add: Vec<Wiring<D>>, mul: Vec<Wiring<D>>) -> Self {
        // let num_vars = &add.len() + &mul.len();
        Self { add, mul }
    }
}

// TODO - implement the circuit conversion into MLEs
impl<F: PrimeField, const D: usize> Into<[SparseMultilinearExtension<F>; 2]> for &Layer<D> {
    fn into(self) -> [SparseMultilinearExtension<F>; 2] {
        // Assume uniform circuit sizes for now
        let le_indices = to_le_indices(3 * D);
        // let D = self.num_vars / 3;
        println!("D: {}", D);
        println!("le_indices: {:?}", le_indices);
        println!("add: {:?}", self.add);
        let add_indices: Vec<usize> = self
            .add
            .iter()
            .map(|w| {
                // construct the index from curr_index, left, right
                let index_be: usize = w.curr_index << (2 * D) + w.left << (D) + w.right;
                le_indices[index_be]
                // index
            })
            .collect();
        let add_evals: Vec<(usize, F)> = add_indices.iter().map(|i| (*i, F::one())).collect();
        let add_mle = SparseMultilinearExtension::from_evaluations(3 * D, &add_evals);
        // same for mul
        let mul_indices: Vec<usize> = self
            .mul
            .iter()
            .map(|w| {
                // construct the index from curr_index, left, right
                let index_be: usize = (w.curr_index << (2 * D)) + (w.left << (D)) + w.right;
                le_indices[index_be]
                // index
            })
            .collect();
        let mul_evals: Vec<(usize, F)> = mul_indices.iter().map(|i| (*i, F::one())).collect();
        let mul_mle = SparseMultilinearExtension::from_evaluations(3 * D, &mul_evals);
        // return both
        [add_mle, mul_mle]
    }
}

pub struct UniformCircuit<F, const D: usize> {
    layers: Vec<Layer<D>>,
    phantom: PhantomData<F>,
}

impl<F: PrimeField, const D: usize> UniformCircuit<F, D> {
    fn new(layers: Vec<Layer<D>>) -> Self {
        Self {
            layers,
            phantom: PhantomData::<F>,
        }

        // test indices
    }

    // not as in the GKR protocol - plain circuit evaluation
    fn evaluate(&self, x: Vec<F>) -> Vec<Vec<F>> {
        let mut evals = Vec::new();

        let mut last_layer = x;

        for layer in self.layers.iter().rev() {
            let mut new_layer: Vec<F> = vec![F::zero(); layer.add.len() + layer.mul.len()];

            // handle addition
            for Wiring {
                curr_index,
                left,
                right,
            } in layer.add.iter()
            {
                new_layer[*curr_index] = last_layer[*left] + last_layer[*right];
            }

            // handle mul
            for Wiring {
                curr_index,
                left,
                right,
            } in layer.mul.iter()
            {
                new_layer[*curr_index] = last_layer[*left] * last_layer[*right];
            }

            evals.push(new_layer.clone());

            last_layer = new_layer;
        }

        evals
    }
}
