
# Sum-check :repeat:

A toy implementation of the sum-check protocol for zero-knowledge proofs of polynomial sums over binary hypercubes. We follow the description in Thaler's [Proofs, Arguments, and Zero-Knowledge](https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf) (section 4.1).

> **Important disclaimer: this project is *completely unsuitable for use in real-world cryptographic applications*. Its purposes are purely illustrative.**

## Usage

The interface is very crude at the moment: in order to execute the protocol, you may modify the hard-coded polynomial, modulus and primitive root  (the $p$ in $\mathbb{F}_p$ and a generator of $\mathbb{F}_p^\ast$) in `main.rs`, compile and execute. If the flag `verbose` is set to `true` (its default value), the prover and verifier will print some information messages during execution. These three parameters are marked with the comment `// MODIFY`.

The polynomial is given as in [this example](https://docs.rs/ark-poly/latest/ark_poly/polynomial/multivariate/struct.SparsePolynomial.html#examples-2) from the `ark_poly` crate.

In order to simulate the exact execution in the example from page 36 in the aforementioned book, leave the original polynomial from `main.rs` unchanged and uncomment the only block comment in the file `parties.rs`, which is also marked with the line comment `// MODIFY`. This prescribes which values the verifier should send to the prover (these are chosen randomly during a regular execution).

## Tests

A few tests are included for various functions. The last of them, `test_protocol` tests the entire protocol for the polynomial in the example from page 36 cited above (in order to simulate that exact execution, follow the instructions at the end of section [Usage]).

The tests can be run with `cargo test` (or `cargo test -- --nocapture` if you don't want the output to be suppresed).

## Dependencies

- `ark-ff` for finite field arithmetic
- `ark-poly` for multivariate polynomials (an ad-hoc implementation of univariate polynomials is included in `polynomials.rs`)
- `rand` and `rand_chacha` for (cryptografically secure) choice of random field elements by the verifier.
