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
pub struct Wiring<const s: usize> {
    curr_index: usize,
    left: usize,
    right: usize,
}

impl<const s: usize> Wiring<s> {
    fn new(curr_index: usize, left: usize, right: usize) -> Self {
        assert!(curr_index < 1 << s);
        assert!(left < 1 << s);
        assert!(right < 1 << s);
        Self {
            curr_index,
            left,
            right,
        }
    }
}

#[derive(Debug)]
pub struct Layer<const s: usize> {
    // pub num_vars: usize,
    pub add: Vec<Wiring<s>>,
    pub mul: Vec<Wiring<s>>,
}

impl<const s: usize> Layer<s> {
    fn new(add: Vec<Wiring<s>>, mul: Vec<Wiring<s>>) -> Self {
        // let num_vars = &add.len() + &mul.len();
        Self { add, mul }
    }
}

// TODO - implement the circuit conversion into MLEs
impl<F: PrimeField, const s: usize> Into<[SparseMultilinearExtension<F>; 2]> for &Layer<s> {
    fn into(self) -> [SparseMultilinearExtension<F>; 2] {
        // Assume uniform circuit sizes for now
        let le_indices = to_le_indices(3 * s);
        let add_indices: Vec<usize> = self
            .add
            .iter()
            .map(|w| {
                // construct the index from curr_index, left, right
                let index_be: usize = (w.curr_index << (2 * s)) + (w.left << (s)) + w.right;
                le_indices[index_be]
                // index
            })
            .collect();
        let add_evals: Vec<(usize, F)> = add_indices.iter().map(|i| (*i, F::one())).collect();
        let add_mle = SparseMultilinearExtension::from_evaluations(3 * s, &add_evals);
        // same for mul
        let mul_indices: Vec<usize> = self
            .mul
            .iter()
            .map(|w| {
                // construct the index from curr_index, left, right
                let index_be: usize = (w.curr_index << (2 * s)) + (w.left << (s)) + w.right;
                le_indices[index_be]
                // index
            })
            .collect();
        let mul_evals: Vec<(usize, F)> = mul_indices.iter().map(|i| (*i, F::one())).collect();
        let mul_mle = SparseMultilinearExtension::from_evaluations(3 * s, &mul_evals);
        // return both
        [add_mle, mul_mle]
    }
}

pub struct UniformCircuit<F, const s: usize> {
    layers: Vec<Layer<s>>,
    phantom: PhantomData<F>,
}

impl<F: PrimeField, const s: usize> UniformCircuit<F, s> {
    /// Pads each layer to have 1<<s elements, so that we have a uniform circuit representation. We do this by adding padding addition gates
    pub fn new(layers: Vec<Layer<s>>) -> Self {
        let mut layers = layers;
        for layer in layers.iter_mut() {
            let num_gates = layer.add.len() + layer.mul.len();
            let diff = (1 << s) - (num_gates);
            for i in 0..diff {
                layer.add.push(Wiring::new(i + num_gates, 0, 0));
            }
        }

        Self {
            layers,
            phantom: PhantomData::<F>,
        }
    }

    /// Plain circuit evaluation given the input x
    /// Asserts that each layer is uniform with 1<<s evaluations
    pub fn evaluate(&self, x: Vec<F>) -> Vec<Vec<F>> {
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

            assert_eq!(new_layer.len(), 1 << s, "non-uniform circuit");
            evals.push(new_layer.clone());

            last_layer = new_layer;
        }

        evals
    }
}
