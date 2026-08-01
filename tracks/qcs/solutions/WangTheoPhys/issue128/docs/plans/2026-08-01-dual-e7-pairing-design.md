# Dual-Only E7 Pairing Design

**Date:** 2026-08-01

**Status:** Approved under the user's continuing instruction to pursue the
finite-step theorem without submitting new HPC work.

## Objective

Compute the exact cubic-field scalar

```text
tau(W_L E7_L),  W_L = H_L^2 - (3L^2/8) I + H_L/2,
```

without constructing or storing the complete E7 Pauli map.  Use the result to
replace the loose generic E7 majorant in the dual-pairing remainder ledger.

## Approaches

1. **Full physical E7 map on shared HPC:** maximally reusable but evaluates
   millions of coefficients that are orthogonal to `W`; deferred.
2. **Deterministic local dual contraction (selected):** stream every exact
   seventh-degree word, retain only support at most four, contract directly
   with `H` and `H^2`, and reduce exact scalar shards.
3. **Improve the generic E7 l1 majorant:** cheap but would need roughly three
   orders of magnitude of improvement and still would not exploit exact
   commutant orthogonality.

## Algorithm

The seventh-degree free logarithm has 16,380 words grouped into 4,096 common
six-letter suffixes.  For each suffix group:

1. evaluate its shared nested-commutator tail once;
2. evaluate the possible outer colors;
3. discard Pauli strings of support greater than four, because neither `H`
   nor `H^2` contains them;
4. map one representative onto an alias-free `L=6` torus;
5. use even-translation invariance to multiply the target coefficient by the
   nine two-by-two cells rather than repeating nine lookups;
6. multiply by the exact free-word cubic coefficient and common rational
   Dynkin/dyadic denominator; and
7. release the group cache before continuing.

Assign sorted suffix-group ordinal `g` to shard `g mod shard_count`.  Four
local shards keep each process under the ten-minute local-compute target.

## Verification

- A degree-five run must reproduce the already certified `tau(H E5)`,
  `tau(H^2 E5)`, and `tau(W E5)` values.
- Every shard records its exact group ordinals, word count, nonzero-word count,
  retained-term count, source hashes, and three cubic partial sums.
- The reducer requires identical configuration/source hashes, disjoint group
  sets, and exact coverage of all suffix groups.
- Forward and reverse shard reduction orders must agree exactly.
- The final E7 artifact remains a scalar dual certificate; it must not claim
  the full E7 operator map or an all-order finite-step theorem.

## Claim boundary

An exact `tau(W E7)` closes only the first omitted logarithm coefficient.  The
finite-step status remains inconclusive until an E9-and-higher dual remainder
and a compatible logarithm branch are certified.  No shared HPC job is
authorized.
