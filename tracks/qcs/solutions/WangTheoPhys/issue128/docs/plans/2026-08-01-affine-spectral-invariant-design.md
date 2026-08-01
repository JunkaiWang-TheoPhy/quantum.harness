# Affine spectral-invariant closure design

**Status:** Approved as a rigor amendment to the finite-step Paper B gate.

## Problem

The exact E5/E7/E9 contraction and compatible E11-and-higher tail certify a
nonzero pairing

```text
p = tau(W D),
W = H^2 - 54 I + H/2,
D = log(S(1/97))/(1/97) - H.
```

This excludes the tangent gauge
`image(i ad_H) + span(I,H)`.  By itself it does not exclude the entire
nonlinear affine unitary orbit `a I + b U H U^dagger`: a finite conjugation can
have nonzero pairing with a witness that annihilates the tangent space.  A
publication-level finite-step eigenphase claim therefore needs one more
nonlinear remainder gate.

## Considered routes

1. **Stop at a finite-step local-log obstruction.**  This is sound and needs no
   new mathematics, but it cannot be advertised as a finite-step spectral
   no-go.
2. **Use an affine spectral moment invariant.**  This is the selected route.
   It reuses the same witness pairing and the already certified principal-log
   norm envelope, adds only exact rational arithmetic, and has a positive
   feasibility margin on the frozen instance.
3. **Enclose the full spectrum.**  A direct `2^144`-dimensional spectral
   computation is neither feasible nor necessary.

## Exact invariant

For Hermitian `A`, define its centered part and normalized moments by

```text
C(A) = A - tau(A) I,
m_k(A) = tau(C(A)^k).
```

For the frozen Hamiltonian, `m_2(H)=54` and `m_3(H)=-27`.  Define

```text
Phi_H(A) = 54^3 m_3(A)^2 - 27^2 m_2(A)^3.
```

If `A=a I+b U H U^dagger`, then

```text
m_2(A)=b^2 m_2(H),
m_3(A)=b^3 m_3(H),
```

and hence `Phi_H(A)=0` exactly, for either sign of `b`.  Thus a certified
nonzero value of `Phi_H(log(S(h))/h)` excludes endpoint conjugation, global
phase, and global time/energy calibration at the actual finite step.

## Linear term and nonlinear remainder

Write `A=H+D` and `E=D-tau(D)I`.  Set

```text
d2 = 2 tau(H E) + tau(E^2),
d3 = 3 tau(H^2 E) + 3 tau(H E^2) + tau(E^3).
```

Expanding the invariant gives the exact linear identity

```text
Phi_H(H+D)_linear = -162 * 54^3 * tau(W D).
```

Let `delta` certify `||D|| <= delta`.  Centering gives
`||E|| <= epsilon = 2 delta`.  With the already certified rational cap
`sqrt(tau(H^2)) <= 15/2`, define

```text
L2 = 15 epsilon,
N2 = epsilon^2,
A2 = L2 + N2,

L3 = 162 epsilon,
N3 = (45/2) epsilon^2 + epsilon^3,
A3 = L3 + N3.
```

Triangle inequality and normalized Schatten inequalities give the rational
remainder cap

```text
R_Phi = 54^3 [A3^2 + 2*27*N3]
      + 27^2 [3*54^2*N2 + 3*54*A2^2 + A2^3].
```

The dual certificate provides a per-cell lower bound `m_dual` on
`|tau(WD)|/(N/4)`.  For 36 cells the finite-step affine-spectral gate is

```text
162 * 54^3 * 36 * m_dual - R_Phi > 0.
```

The tail artifact bounds the one-step logarithm defect.  The combiner must use

```text
delta = 97 * envelope.log_defect
```

for the effective Hamiltonian `log(S(1/97))/(1/97)`; using the one-step value
without this factor is a verifier error.

## Certificate boundary

The existing signed dual margin remains a separately reported field.  It may
establish `certified_local_log_obstruction`.  The public finite-step spectral
status changes to `certified_affine_spectral_obstruction` only when the exact
lower endpoint of the affine-invariant margin is strictly positive.  A zero or
negative endpoint leaves the spectral claim `inconclusive`, even if the dual
pairing is nonzero.

The production JSON binds E5, E7, E9, the dual-log tail, the manifest index,
the effective-log norm cap, all implementation sources, and both arithmetic
margins.  Verifiers reject a missing factor of 97, a lowered norm cap, a
hand-edited claim, or any source/digest drift.

## Feasibility result

Before E9, the exact E5+E7 dual margin is approximately
`6.2741778358e-11` per cell.  Using the committed log-defect cap gives

```text
epsilon ~= 2.8446418973e-6,
R_Phi ~= 3.5254820251e-2,
linear invariant lower bound ~= 5.7617660331e-2.
```

The conservative nonlinear cap consumes about `61.2%` of the pre-E9 linear
bound, leaving a plausible positive gate.  These decimals are diagnostics;
the certificate uses only exact rational endpoints and the verified E9 value.
