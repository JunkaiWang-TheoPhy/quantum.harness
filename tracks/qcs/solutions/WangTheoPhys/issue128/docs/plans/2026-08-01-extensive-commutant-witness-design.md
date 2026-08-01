# Extensive Commutant-Witness Theorem Design

**Date:** 2026-08-01

**Status:** Approved by the user's instruction to proceed after the proposed
system-size, normalization, and finite-step theorem sequence.

## Objective

Determine whether the exact 12-by-12 witness

```text
W_L = H_L^2 - tau(H_L^2) I + H_L/2
```

is one member of an extensive family on periodic even square lattices.  If it
is, certify the system-size dependence of its leading-defect pairing and its
Hilbert--Schmidt normalization, then formulate the exact remainder inequality
needed to promote the leading-order spectral no-go to a finite-step theorem.

## Claim boundary

The first phase may prove only an exact family of leading-order quotient
obstructions.  It must not claim a finite-step eigenphase lower bound until a
local-log branch and an all-order remainder have been bounded in the same
dual pairing.  Small tori whose periodic identifications alias a local E5
cluster must be identified and excluded rather than silently folded into the
stable-volume formula.

## Approaches considered

### A. Exact size scan followed by linked-cluster proof (selected)

Generate the exact E5 density once, lift it onto several even tori, and compute
`tau(H^2)`, `tau(H^3)`, `tau(H E5)`, `tau(H^2 E5)`, `tau(W E5)`, and
`tau(W^2)`.  Detect the first alias-free size from support geometry and test
whether each extensive quantity divided by the number of two-by-two cells is
constant.  Convert stable observations into exact local counting identities.

This route produces both discovery data and an independently checkable
theorem object.  It is computationally modest because no dense matrix or HPC
job is required.

### B. Direct 12-by-12 finite-step remainder

Normalize the existing witness immediately and combine it with the current
operator-norm ledger.  This is faster to state but risks a volume-loose bound
and would leave the result tied to one lattice size.

### C. Cross-model witness survey

Repeat the quotient construction for other Hamiltonians.  This could broaden
the paper, but it is premature before the Heisenberg witness has a clean
system-size theorem and a finite-step promotion criterion.

## Exact data model

For each admissible even `L`, emit exact rational or cubic-field coordinates
for:

```text
N = L^2
M = N/4                         number of 2x2 cells
tau(H^2)
tau(H^3)
tau(H^4) = ||H^2||_HS^2
tau(H E5)
tau(H^2 E5)
tau(W E5)
tau(W^2)
rho_L = tau(W E5)^2 / tau(W^2)
```

The squared normalized obstruction `rho_L` stays in the cubic field and
avoids introducing uncertified square roots.  It is the exact squared dual
lower bound obtained from the normalized witness `W/sqrt(tau(W^2))`.

## Aliasing and normalization checks

Before evaluating a size, determine every E5 Pauli term's coordinate diameter.
Reject a torus if two distinct coordinates in any term coincide modulo `L`.
For accepted sizes:

1. lift all four matching densities and reproduce the full Hamiltonian;
2. construct `H^2` using cancellation of anticommuting ordered pairs;
3. obtain `tau(H^4)` from the Pauli coefficient norm of `H^2`;
4. obtain all E5 pairings by sparse translated lookup; and
5. verify the polynomial commutant and calibration orthogonality identities.

## Expected analytic structure

Existing exact data at `L=4` and `L=12` suggest

```text
tau(H^2) = 3 N / 8,
tau(H^3) = -3 N / 16,
```

so the Hamiltonian coefficient in `W_L` is universally `+1/2`.  The size scan
will test, not assume, the corresponding stable formula for `tau(W E5)` and
the polynomial in `N` governing `tau(W^2)`.

## Finite-step promotion criterion

For total time one and `r` product-formula steps, write the calibrated
effective logarithm defect schematically as

```text
Delta_r = E5/r^4 + R_r.
```

The exact leading obstruction survives whenever

```text
|tau(W_L R_r)| < |tau(W_L E5)| / r^4.
```

A sufficient norm-based condition is

```text
||R_r||_HS < |tau(W_L E5)| /
             (r^4 sqrt(tau(W_L^2))).
```

The implementation will emit this threshold exactly as a squared algebraic
quantity.  It will not substitute the current operator-norm ledger unless the
normalizations and logarithm branch are proved compatible.

## Deliverables

- an exact multi-size scan and verifier;
- a deterministic JSON artifact with alias decisions and all moments;
- tests for moment identities, alias rejection, cubic normalization, and
  mutation resistance;
- a theorem note separating verified finite-size identities, conjectured
  stable formulas, and analytically proved counting relations;
- an updated go/no-go decision for the finite-step remainder phase.

No shared HPC job is authorized in this design.
