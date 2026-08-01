# Gauge-Aware Spectral Certificate Design

**Date:** 2026-08-01

**Status:** Approved for implementation by the user's "做吧" instruction after
the venue-readiness review.

## Objective

Build the first exact, fail-closed prototype for distinguishing three kinds of
leading product-formula defect:

1. a commutator term that changes only the eigenbasis;
2. an identity or Hamiltonian-parallel term that is removable by global-phase
   or time/energy calibration; and
3. a residual term that is a genuine spectral obstruction.

The prototype must either construct an exact correction (a primal
certificate) or return an exact orthogonal witness (a dual certificate).  It
must not turn a restricted-processor obstruction into a spectral no-go claim.

## Current evidence and claim boundary

The frozen fourth-order endpoint-processor audit proves

```text
<H, E5> != 0, while <H, i[Q,H]> = 0 for every Hermitian Q.
```

This proves that the complete leading defect cannot be removed by endpoint
conjugation when the target time is fixed.  It does not by itself survive an
allowed time recalibration because `H` belongs to the calibration direction.

The same audit proves that support-at-most-four processors cannot cancel the
support-six part of the operator defect.  That is a restricted-processor
operator no-go, not automatically a spectral lower bound.  A spectral no-go
requires a witness in the commutant of `H`, orthogonal to both `I` and `H`, or
an equivalent exact block-diagonal obstruction.

## Approaches considered

### A. Documentation-only correction

Amend the obstruction report to distinguish fixed-time and gauge-aware
claims.  This is cheap and necessary, but it produces no reusable theorem
object and no executable evidence for a paper.

### B. Exact primal/dual gauge prototype (selected)

Implement an exact rational linear-algebra layer for the quotient

```text
Hermitian operators /
(image(i ad_H) + span{I,H}).
```

For small exact matrices, generate the full Hermitian commutator image
automatically.  For coefficient-space problems, accept an explicit generator
span and label it as either complete or restricted.  Return either an exact
gauge reconstruction or a primitive integer witness.  This directly tests the
paper's proposed normal-form logic and exposes the distinction between full
spectral and restricted local claims.

### C. Immediate physical E7/HPC enumeration

Run the shared 16,380-word s10/s11 evaluator now.  This could provide a useful
coefficient but would not resolve the local-log remainder, the calibration
gauge, or the spectral-versus-support claim boundary.  It is deferred until
the theorem prototype leaves a plausible finite-step budget.

## Architecture

### Exact coefficient-space solver

Create a focused module that receives a target vector and named gauge
generators over exact SymPy rationals.  It will:

1. select an independent generator basis exactly;
2. compute the Hilbert--Schmidt orthogonal projection;
3. return the exact residual;
4. convert a nonzero rational residual to a primitive integer witness; and
5. verify orthogonality and the nonzero target pairing.

The result carries one of two explicit scopes:

- `full_commutator_image`: orthogonality supports a spectral obstruction;
- `restricted_processor_span`: orthogonality supports only a restricted
  processor no-go.

No caller may promote the second scope to the first without a separate span
completeness identifier.

### Small-matrix reference construction

For an exact Hermitian matrix `H`, construct a real Hermitian basis consisting
of diagonal, symmetric off-diagonal, and imaginary antisymmetric matrix units.
The vectors of `i[K,H]` over this complete basis generate the full commutator
image.  Add `I` and `H` as optional gauge directions, decompose an exact defect
`E`, and reconstruct the witness matrix from its coordinates.

The reference audit will cover three decisive cases:

1. `E` off diagonal: removable by a commutator;
2. `E` parallel to `H`: obstructed at fixed time but removable after time
   calibration; and
3. a diagonal `E` orthogonal to `I` and `H`: a calibration-robust spectral
   obstruction.

### Issue-128 claim audit

Create a deterministic JSON report tied by SHA-256 to the existing exact
processor-obstruction artifact.  It will preserve the valid fixed-time no-go,
mark the current gauge-aware spectral status `inconclusive`, and explain that
no independent commutant witness orthogonal to `I,H` has yet been attached.
The report verifier will reject a forged `spectral_no_go` status.

## Data flow

```text
exact H,E or exact coefficient vectors
          |
          v
full/restricted gauge generator span + {I,H}
          |
          v
exact projection --------------------> primal reconstruction
          |
          +--------------------------> primitive dual witness
                                          |
                                          v
                       scope-aware verifier and JSON audit
```

## Error handling and fail-closed rules

- Reject nonexact floating-point entries.
- Reject dimension mismatches and empty targets.
- Reject non-Hermitian matrix inputs.
- Reject a `full_commutator_image` claim without automatic full-basis
  construction or an explicit completeness identifier.
- Never label a zero residual as an obstruction.
- Never label a restricted-span residual as a spectral obstruction.
- Preserve fixed-time and calibrated-time decisions as separate fields.
- Bind every Issue-128 audit to the source artifact digest.

## Testing

Unit tests will check exact projection, primitive-witness normalization,
orthogonality, primal reconstruction, non-Hermitian rejection, float
rejection, and scope enforcement.  Matrix tests will verify the three
reference cases above.  Artifact tests will mutate the source hash, status,
and claim scope and require verification failure.

The focused suite must pass before the existing non-slow Issue-128 suite is
run.  No HPC job is authorized by this design.

## Success and promotion gate

This phase succeeds when the repository contains a reproducible, exact report
that demonstrates all three gauge cases and accurately downgrades the current
retiming-robust no-go status to `inconclusive` without weakening the valid
fixed-time result.

The next scientific gate is then concrete: find a nonzero commutant witness
orthogonal to `I,H`, or close the processed local-log remainder tightly enough
to justify physical E7.  Only the latter can authorize a new shared HPC run.
