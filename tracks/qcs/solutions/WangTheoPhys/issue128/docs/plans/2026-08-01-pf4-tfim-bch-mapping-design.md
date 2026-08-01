# Exact PF4--TFIM BCH Mapping Design

## Purpose

Replace the unverified, plan-stated PF4 quadratic form with an independently
derived exact free-trace identity for the actual five-copy fourth-order Suzuki
formula.  Use that identity to decide, without numerical fitting, whether the
periodic transverse-field Ising family has a leading-order obstruction to a
pure endpoint-conjugation processor at fixed time and fixed Hamiltonian
normalization.

This result is a Paper B family theorem.  It is not the finite-step Heisenberg
certificate, and it does not establish a gauge-aware affine obstruction under
arbitrary additions of `a I + b H`.

## Approaches considered

### A. Exact cyclic free-word derivation and certificate promotion

Expand the real-generator logarithm of the exact cubic Suzuki stages through
degree five, multiply by `H=A+B`, and quotient degree-six words by cyclic trace
equivalence.  Independently construct the cyclic classes of the nested
commutators

\[
C=[A,[A,B]],\qquad D=[B,[B,A]],
\]

and prove equality coefficient by coefficient.  Cross-check the same identity
on exact rational Hermitian matrices before changing any claim status.

This is the selected approach.  It uses the formula actually certified by the
repository, exposes the incorrect old coefficients, and produces a compact
proof object that can be checked without physical Pauli enumeration.

### B. Preserve the conditional TFIM artifact

Record the newly observed coefficients only as a diagnostic and retain the
current `blocked_pending_pf4_bch_mapping` claim.  This is safe but wastes an
exact identity already available from the committed cubic-series engine and
does not satisfy the Paper B family-theorem gate.

### C. Build a full gauge-aware TFIM commutant witness

Project the PF4 defect away from `span{I,H}` inside every degenerate energy
block and construct a trace-norm-normalized witness for the resulting quotient.
This would address the larger affine orbit, but it is a separate theorem with
additional degeneracy and normalization work.  It is intentionally excluded
from this implementation.

## Exact identity

Let `alpha` be the positive real root of `alpha^3=4`, and let the five-copy
Suzuki scale be

\[
u=(4-\alpha)^{-1}=\frac4{15}+\frac1{15}\alpha+
\frac1{60}\alpha^2.
\]

For the stage convention implemented by
`fourth_order_suzuki_cubic_stages(2)`, write

\[
\log S_4(h)=h(A+B)+h^5E_5+O(h^7).
\]

The candidate identity to prove exactly is

\[
\operatorname{Tr}((A+B)E_5)
=\gamma\operatorname{Tr}\!\left(
C^2-4CD+\frac83D^2
\right),
\]

where

\[
\gamma=
\frac{37}{900000}
+\frac{313}{14400000}\alpha
+\frac{29}{1800000}\alpha^2>0.
\]

The equality is a cyclic free-word identity: after canonical cyclic rotation,
both sides have the same ten degree-six word classes and the same exact cubic
coefficient in every class.  Positivity of `gamma` follows directly from the
positive coordinates and `alpha>0`; it is not inferred from floating-point
evaluation.

For physical unitary stages `exp(-i h A)` and `exp(-i h B)`, homogeneity gives
the same Hermitian `E5` coefficient in the principal effective Hamiltonian,
because `(-i)^5=-i` and multiplication of the logarithm by `i/h` removes that
factor.

## TFIM specialization and theorem boundary

For the periodic even-length family

\[
H=A+B,\qquad
A=h\sum_i X_i,\qquad
B=j\sum_i Z_iZ_{i+1},
\]

exact Pauli counting already proves

\[
\frac{\operatorname{Tr}(C^2)}{d}=128Lh^4j^2,\quad
\frac{\operatorname{Tr}(CD)}{d}=0,\quad
\frac{\operatorname{Tr}(D^2)}{d}=128Lh^2j^4.
\]

Therefore

\[
\frac{\operatorname{Tr}(HE_5)}{d}
=128L\gamma h^2j^2\left(h^2+\frac83j^2\right).
\]

It is strictly positive exactly when `h*j != 0`.  If either coupling is zero,
the two fragments commute and the obstruction vanishes.

For every Hermitian endpoint generator `Q`, trace cyclicity gives

\[
\operatorname{Tr}(H\,i[Q,H])=0.
\]

Consequently the nonzero pairing excludes `E5=-i[Q,H]` and gives a rigorous
operator-norm lower bound through normalized trace duality.  The committed
certificate may therefore claim:

- exact PF4 BCH/free-trace mapping;
- an all-even-length TFIM leading-order obstruction to pure endpoint
  conjugation at fixed time and normalization;
- an exact algebraic lower-bound coefficient using the existing Pauli-l1
  upper bound on the normalized trace norm of `H`.

It must continue to state:

- no finite-step TFIM no-go is proved;
- no obstruction modulo `a I + b H` is proved;
- no total-time eigenphase lower bound modulo `2 pi` is proved;
- the Heisenberg E9 certificate is a separate finite-step result.

## Software architecture

Create `src/trottercert/pf4_bch_mapping.py` as the small trusted algebraic
kernel.  It owns generic exact cubic word addition, multiplication,
commutators, cyclic canonicalization, the two sides of the identity, and a
strict verifier.  It depends only on the committed `cubic_field.py` and
`cubic_local.py`; it must not depend on the currently untracked
`rational_words.py` or `processed_kernels.py` files.

Modify `trace_obstruction.py` so its canonical PF4 identity record contains the
proved coefficients and exact gamma coordinates.  Keep the existing PF2
identity unchanged.  Replace the old conditional quadratic evaluator with the
correct core `C2 - 4 CD + (8/3) D2`, and provide a separate function returning
the exact cubic Suzuki pairing.

Modify `tfim_obstruction.py` and `certify_tfim_obstruction.py` to serialize the
exact rational moments, cubic pairing, algebraic lower-bound coefficient,
commuting exception, and narrowly scoped promoted status.  Regenerate the
canonical JSON artifact only from the corrected implementation.

## Verification strategy

1. Assert the free-word logarithm has no degrees two through four and obtain
   the degree-five map from the exact cubic stages.
2. Compare all ten cyclic degree-six classes of both sides exactly.
3. Mutate each rational coefficient and each gamma coordinate and require the
   verifier to fail.
4. Evaluate both sides on seeded exact rational symmetric matrices in
   dimensions two and three.
5. Recheck the TFIM Pauli formulas at lengths 4, 6, and 8.
6. Verify the closed family formula for positive, negative, zero, and rational
   couplings.
7. Require artifact verification to reject the old coefficients, conditional
   status, forged promoted status, and source-hash drift.
8. Run the focused PF4/TFIM suites in a clean checkout.  The repository-wide
   suite remains a separate final gate because current HEAD contains committed
   tests whose implementation dependencies are still untracked concurrent
   work.

## Failure policy

Any disagreement between free-word classes, matrix evaluation, Pauli counting,
or frozen artifact recomputation leaves the claim conditional.  No verifier
may accept a status change based only on a recomputed payload digest.  The old
unverified coefficient set is evidence of why the algebraic identity, not the
research plan text, is authoritative.
