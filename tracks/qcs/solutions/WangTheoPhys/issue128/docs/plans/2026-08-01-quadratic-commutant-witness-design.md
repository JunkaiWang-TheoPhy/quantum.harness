# Quadratic Commutant Witness Design

**Date:** 2026-08-01

**Status:** Approved as the next gate in the gauge-aware spectral-certificate
design and authorized by the user's instruction to continue.

## Objective

Turn the candidate polynomial commutant

```text
W = H^2 - tau(H^2) I - tau(H^3)/tau(H^2) H
```

into an exact, reproducible certificate for the periodic 12-by-12 spin-1/2
isotropic Heisenberg model and the fixed five-copy fourth-order Suzuki formula.
Here `tau` is normalized Hilbert--Schmidt trace.  The certificate must prove
that `W` commutes with `H`, is orthogonal to both calibration directions `I`
and `H`, and has nonzero exact pairing with the leading logarithm defect `E5`.

## Claim boundary

A nonzero `tau(W E5)` proves a retiming-robust **leading-order spectral
obstruction**:

```text
E5 not in image(i ad_H) + span{I,H}.
```

It does not by itself prove a finite-step eigenphase lower bound at total time
one.  Promotion to that claim still requires a local logarithm branch theorem
and a certified all-order remainder small enough to preserve the leading
pairing.

## Exact construction

1. Lift every canonical two-by-two-cell symplectic Pauli density through all
   even translations on the periodic 12-by-12 torus.
2. Validate the lift by reconstructing all 864 Pauli terms of the full
   Heisenberg Hamiltonian from the four matching densities.
3. Form the exact rational Pauli map of `H^2` and obtain trace moments by Pauli
   orthogonality, without constructing dense matrices or `H^3`.
4. Generate the exact cubic-field `E5` density and compute `tau(H E5)` and
   `tau(H^2 E5)` by translated sparse lookup.
5. Emit a deterministic JSON certificate containing the formula identifier,
   lattice size, exact moments, exact cubic pairings, witness coefficients,
   algebraic proof obligations, and source hashes.
6. Recompute the certificate in an explicit slow verifier.  Keep ordinary
   unit tests fast by checking the lifting/moment kernel on small exact inputs
   and validating the stored certificate fail-closed.

## Expected exact identities

For the 12-by-12 periodic lattice,

```text
tau(H^2) = 54
tau(H^3) = -27
W = H^2 - 54 I + H/2
```

The exact cubic pairing is nonzero if any coordinate in the rational basis
`(1, alpha, alpha^2)`, `alpha^3 = 4`, is nonzero.

## Files and tests

- Create `src/trottercert/commutant_witness.py` for torus lifting, exact Pauli
  products, trace moments, certificate construction, and lightweight payload
  verification.
- Create `scripts/certify_quadratic_commutant_witness.py` for deterministic
  generation and full recomputation verification.
- Create `tests/test_commutant_witness.py` using test-first development.
- Create `docs/experiments/processor-obstruction/quadratic-commutant-witness.json`.
- Update the gauge-aware audit and its tests so the calibrated status becomes
  `leading_order_no_go`, while the finite-step status remains `inconclusive`.

No HPC work is authorized by this design.
