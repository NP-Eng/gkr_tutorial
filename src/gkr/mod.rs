use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};
use std::marker::PhantomData;

use ark_ff::PrimeField;

use crate::parties::to_le_indices;

// pub mod parties;
// mod tests;

#[derive(Debug)]
pub struct Wiring<const d: usize> {
    curr_index: usize,
    left: usize,
    right: usize,
}

impl<const d: usize> Wiring<d> {
    fn new(curr_index: usize, left: usize, right: usize) -> Self {
        assert!(curr_index < 1 << d);
        assert!(left < 1 << d);
        assert!(right < 1 << d);
        Self {
            curr_index,
            left,
            right,
        }
    }
}

#[derive(Debug)]
pub struct Layer<const d: usize> {
    // pub num_vars: usize,
    pub add: Vec<Wiring<d>>,
    pub mul: Vec<Wiring<d>>,
}

impl<const d: usize> Layer<d> {
    fn new(add: Vec<Wiring<d>>, mul: Vec<Wiring<d>>) -> Self {
        // let num_vars = &add.len() + &mul.len();
        Self { add, mul }
    }
}

// TODO - implement the circuit conversion into MLEs
impl<F: PrimeField, const d: usize> Into<[SparseMultilinearExtension<F>; 2]> for &Layer<d> {
    fn into(self) -> [SparseMultilinearExtension<F>; 2] {
        // Assume uniform circuit sizes for now
        let le_indices = to_le_indices(3 * d);
        // let d = self.num_vars / 3;
        println!("d: {}", d);
        println!("le_indices: {:?}", le_indices);
        println!("add: {:?}", self.add);
        let add_indices: Vec<usize> = self
            .add
            .iter()
            .map(|w| {
                // construct the index from curr_index, left, right
                let index_be: usize = w.curr_index << (2 * d) + w.left << (d) + w.right;
                le_indices[index_be]
                // index
            })
            .collect();
        let add_evals: Vec<(usize, F)> = add_indices.iter().map(|i| (*i, F::one())).collect();
        let add_mle = SparseMultilinearExtension::from_evaluations(3 * d, &add_evals);
        // same for mul
        let mul_indices: Vec<usize> = self
            .mul
            .iter()
            .map(|w| {
                // construct the index from curr_index, left, right
                let index_be: usize = (w.curr_index << (2 * d)) + (w.left << (d)) + w.right;
                le_indices[index_be]
                // index
            })
            .collect();
        let mul_evals: Vec<(usize, F)> = mul_indices.iter().map(|i| (*i, F::one())).collect();
        let mul_mle = SparseMultilinearExtension::from_evaluations(3 * d, &mul_evals);
        // return both
        [add_mle, mul_mle]
    }
}

pub struct UniformCircuit<F, const d: usize> {
    layers: Vec<Layer<d>>,
    phantom: PhantomData<F>,
}

impl<F: PrimeField, const d: usize> UniformCircuit<F, d> {
    fn new(layers: Vec<Layer<d>>) -> Self {
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

#[cfg(test)]
mod test {

    use super::*;
    use ark_bls12_381::Fq;

    #[test]
    fn simple_circuit() {
        // example from Thaler's book p. 60, bottom - one gate changed for addition

        // mul, layer 0 (output)
        let mul0_0 = Wiring::new(0, 0, 1);
        let mul0_1 = Wiring::new(1, 2, 3);

        let layer_0 = Layer::new(Vec::new(), vec![mul0_0, mul0_1]);

        // add, layer 0
        // empty

        // mul, layer 1
        let mul1_0 = Wiring::new(0, 0, 0);
        let mul1_1 = Wiring::new(1, 1, 1);
        let mul1_2 = Wiring::new(2, 1, 2);

        // add, layer 1
        let add1_3 = Wiring::new(3, 3, 3);

        let layer_1 = Layer::new(vec![add1_3], vec![mul1_0, mul1_1, mul1_2]);

        let circuit = UniformCircuit::<Fq, 2>::new(vec![layer_0, layer_1]);

        let computed_out = circuit.evaluate(
            vec![3, 2, 3, 1]
                .iter()
                .map(|x| Fq::from(*x as u64))
                .collect(),
        );

        assert_eq!(
            *computed_out.last().unwrap(),
            vec![36, 12]
                .iter()
                .map(|x| Fq::from(*x as u64))
                .collect::<Vec<Fq>>()
        );
    }

    #[test]
    fn layer_to_mles() {
        let mul0_0 = Wiring::new(0, 0, 1);
        let mul0_1 = Wiring::new(1, 2, 3);

        let layer_0 = Layer::<2>::new(Vec::new(), vec![mul0_0, mul0_1]);

        let mles: [SparseMultilinearExtension<Fq>; 2] = (&layer_0).into();
    }
}
