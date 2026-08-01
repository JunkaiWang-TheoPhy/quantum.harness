# Dual local-log tail via the even right generator

**Date:** 2026-08-01

**Status:** Approved as the concrete E11-and-higher route under the user's
standing `做吧` and `继续` instructions.  This refines, rather than weakens,
the fail-closed claim boundary in the earlier dual-E9 design.

## Objective

Certify, for the frozen periodic `12 x 12` Heisenberg instance and the
principal logarithm at `t = 1/97`, a compatible bound

```text
|sum_(odd d >= 11) q_d / 97^(d-1)| <= B_log_tail,

q_d = tau(W E_d)/(N/4),
W = H^2 - tau(H^2) I + H/2.
```

The target is the same normalized trace pairing used by the exact E5, E7,
and E9 contractions.  A right-generator operator-norm error may be used only
as a proved intermediate inequality; it may not be substituted for the local
logarithm tail.

The finite-step claim remains inconclusive until the exact E9 artifact is
available and its signed contribution leaves strictly more margin than the
verified E11-and-higher bound.

## Fixed-instance scope

The finite-step target is the manuscript's periodic `12 x 12` benchmark:

```text
N = 144 sites,
C = N/4 = 36 two-by-two cells,
tau(W^2) = 23085/4,
max nonzero Pauli coefficient of W = 1/4.
```

The exact E5 and E7 pairings are stable local densities, but the branch and
logarithm-stability part of this certificate is deliberately a fixed-instance
statement.  It does not claim a finite-step theorem uniform in lattice size.

## Approaches considered

### A. Direct polymer or Collatz majorant for `log S`

Expand the logarithm as a connected polymer series and track support return
with a nonnegative transition matrix.  This is the most natural route to a
uniform-in-volume result.  It also requires a new linked-cluster theorem,
support aggregation proof, and degree-uniform domination certificate.  A
support-only absolute recurrence loses the factorial cancellation between
stage exponentials and commutator growth and is numerically useless at the
target.

This route remains appropriate for a later scaling paper, but it is not the
smallest proof surface for the fixed `12 x 12` submission theorem.

### B. Exact E11/E13 coefficients followed by an observed ratio

The manifest machinery can be extended to higher odd degrees.  Exact higher
coefficients would be valuable cross-checks, but no finite list proves an
all-order tail.  An empirical ratio is therefore prohibited as the trusted
closure rule.

### C. Even right generator plus a certified `dexp` correction (selected)

Use time symmetry to isolate the even part of the right generator.  Its
degree-ten-and-higher dual tail has a small, already structured local
majorant.  Then prove, rather than assume, the conversion from that object to
the derivative of the principal logarithm.  Commutation of `W` with `H`
forces the conversion error to contain two product-formula defects, so it
starts at order ten in the derivative and is quantitatively much smaller than
the direct right-generator tail.

This route reuses checked local-growth lemmas, keeps the final quantity exactly
aligned with E5/E7/E9, and produces a short independent rational verifier.

## Algebraic reduction

Let

```text
S(t) = exp(A(t)),
K_R(t) = S'(t) S(t)^(-1).
```

Time symmetry gives `A(-t) = -A(t)`.  Differentiating
`S(-t) = S(t)^(-1)` shows that the even part of the right generator is

```text
K_even(t) = (K_R(t) + K_R(-t))/2
          = sinh(ad_A(t))/ad_A(t) applied to A'(t).
```

Consequently

```text
K_even - A'
  = sum_(k >= 1) ad_A^(2k)(A') / (2k+1)!.
```

The normalized trace pairing is invariant under moving commutators between
factors.  Since `[W,H] = 0`, write

```text
A(t) = t H + R(t).
```

For every `k >= 1`,

```text
|tau(W ad_A^(2k)(A'))|
 <= ||ad_A^(2k-1)(W)||_2 ||[A,A']||_2
 <= 2 ||R|| ||W||_2 (2||A||)^(2k-2) ||[A,A']||_2.
```

The first commutator on `W` is `[R,W]`, not `[tH,W]`.  Moreover

```text
[A,A'] = [tH,R'] + [R,H] + [R,R'],
R' = A' - H.
```

Thus both factors vanish with the formula defect.  The scalar series obeys,
for `||A|| <= 3/2`,

```text
sum_(k >= 1) (2||A||)^(2k-2)/(2k+1)! <= 10/33.
```

The bound follows by taking the first term `1/3!` and bounding successive
term ratios by `9/20`.

## Centered branch and logarithm stability

Each matching Hamiltonian may be shifted by a scalar without changing any
commutator or pairing with `W`.  The centered matching radius is `N/4 = 36`,
and the centered sum has the conservative operator bound

```text
||H_centered|| <= N = 144.
```

At `t = 1/97`, the exact-evolution phase radius is therefore
`144/97 < 3/2`.

Let `B_avg(t)` be the verified average right-generator defect and define

```text
epsilon_step(t) = t B_avg(t).
```

Duhamel gives

```text
||S(t) - exp(t H)|| <= epsilon_step(t).
```

The relative principal logarithm is bounded without floating point by

```text
||log(exp(-tH) S(t))||
 <= epsilon_step/(1-epsilon_step).
```

Along the corresponding unitary geodesic, all centered logarithms remain in
the radius `3/2` ball.  On normalized Hilbert--Schmidt space the inverse
`dexp` multiplier satisfies

```text
||dexp_A^(-1)||_(2 -> 2)
 <= ||A||/sin(||A||)
 <= 1/(1-||A||^2/6)
 <= 8/5.
```

The last two inequalities use the alternating sine lower bound on
`[0,3/2]`.  Hence

```text
||R(t)||
 <= (8/5) epsilon_step/(1-epsilon_step).
```

This also proves the required centered principal-log branch independently of
the older path-length-only bound.

## Bound for `R'`

On Hilbert--Schmidt space,

```text
A' = f(ad_A) K_R,
f(z) = z/(exp(z)-1).
```

For imaginary `z` with `|z| <= 3`, the certificate uses the rational scalar
bounds

```text
|f(z)| <= 8/5,
|(f(z)-1)/z| <= 1.
```

The first is the same sine bound as above.  The second will be checked by a
small interval-polynomial lemma, not accepted as an artifact input.  Since
`ad_A(H) = ad_R(H)`, this yields

```text
||R'||_2
 <= (8/5) ||K_R-H||_2 + 2 ||R|| ||H||_2.
```

The exact moments give `||H||_2 = sqrt(54) < 15/2`.  No extensive
Hilbert--Schmidt estimate is inserted by hand.

## Direct even-generator tail

For one stage with absolute coefficient `a` and later-stage absolute prefix
`p`, a degree-`n` conjugation term has:

- Taylor factor `a p^n/n!`;
- local nested-commutator cell bound `(3/2)(n+1)!`;
- dual coefficient bound `max |W_P| = 1/4`; and
- average time factor `1/(n+1)`.

The factorials and the cell normalization cancel exactly.  Summing all
degrees `n >= 10` gives the rational compatible bound

```text
B_direct(t)
 = (3/8) sum_stages a (p t)^10/(1-p t).
```

This is algebraically the same scalar expression as
`defect_tail_site_bound(..., first_omitted_degree=10)`, but the certificate
must document and verify the new dual-pairing normalization.  It is not
labeled an operator-norm substitution.

At `t=1/97`, the current exact-rational implementation evaluates to

```text
B_direct = 6.448912075284801e-12  (outward decimal summary).
```

## `dexp` correction tail

Use the existing exact D4--D7 local majorants and the degree-eight-and-higher
right-generator majorant to construct two monotone functions:

```text
B_point(s) >= ||K_R(s)-H||,
B_avg(s)   >= (1/s) integral_0^s ||K_R(u)-H|| du.
```

For the geometric stage tail, the pointwise scalar sum is recomputed as

```text
sum_(n >= q) (n+1) x^n
 = x^q ((q+1)-q x)/(1-x)^2.
```

The verified fourth-order cancellations imply the monotone power laws on
`0 <= s <= t`:

```text
B_point(s) <= (s/t)^4 B_point(t),
||R(s)||   <= (s/t)^5 ||R(t)||,
||R'(s)||_2 <= (s/t)^4 ||R'(t)||_2,
||[A,A'](s)||_2 <= (s/t)^5 C_comm(t).
```

With

```text
C_comm(t) = 2 (t ||H_centered|| ||R'||_2
                 + ||R|| ||H||_2
                 + ||R|| ||R'||_2),
```

integration of the order-ten product gives

```text
B_dexp
 <= (||W||_2/C) * 2 ||R|| C_comm * (10/33) / 11.
```

Use the exact moment bound

```text
||W||_2/C = sqrt(23085/4)/36 < 17/8.
```

The rational pre-audit at `t=1/97` gives the following outward decimal
summaries:

| quantity | bound |
|---|---:|
| average global `K-H` defect | `8.889505847599637e-7` |
| one-step unitary defect | `9.164439018143955e-9` |
| pointwise global `K-H` defect | `5.478749822410207e-6` |
| `||R||` | `1.466310256340944e-8` |
| `||R'||_2` | `8.985946254307472e-6` |
| `C_comm` | `2.689986970136198e-5` |
| `B_dexp` | `4.618047053958782e-14` |
| `B_direct + B_dexp` | `6.495092545992439e-12` |

The combined pre-audit is about `9.4%` of the current E5+E7 signed magnitude
at `r=97`.  This is a feasibility result, not the final gate: exact E9 can
change the signed margin and must be incorporated before promotion.

### Implementation provenance refinement

The initial feasibility audit regenerated generic 30-digit/24-digit
right-generator constants from the live source tree.  Final artifact work
found that two of those older source files carried unrelated uncommitted cache
and scan additions.  The production certificate therefore does not bind or
consume that dirty worktree state.  It binds the committed schema-v3
D5-integrated main certificate and exactly recovers its accepted D4--D7
density caps from the recorded contributions.  The new module independently
rebuilds the formula's 31 interval stages at the main certificate's frozen
12-digit precision.  This both restores clean-checkout reproducibility and
tightens the conversion term to the values above.

## Certificate architecture

Add one scientific module with pure exact-rational functions for:

1. the pointwise and average stage-tail sums;
2. centered branch/logarithm-stability gates;
3. `R`, `R'`, and commutator bounds;
4. direct dual right-generator tail;
5. `dexp` correction and total E11+ bound; and
6. state-by-state verification of every scalar precondition.

Add one CLI that emits canonical JSON containing:

- formula and fixed-instance identifiers;
- source SHA-256 values;
- exact moment identities and rational square-root caps;
- exact D4--D7 constants recovered from the hash-bound schema-v3 main
  certificate and independently rebuilt interval stages;
- every intermediate rational bound;
- the direct, conversion, and total tail values;
- the frozen E5/E7 artifact hashes;
- an optional exact E9 artifact hash and signed-margin decision; and
- a fail-closed claim block.

Without E9 the artifact must say

```text
tail_bound = certified,
finite_step_status = inconclusive,
missing = exact dual E9 pairing.
```

With E9, the verifier recomputes the exact cubic E5+E7+E9 interval and promotes
only if the total E11+ rational bound is strictly below its distance from zero.

## Independent verification and mutation policy

The verifier must regenerate all mathematical values.  It rejects:

- changed lattice size, step count, formula, or source hashes;
- a forged `max |W_P|`, moment, square-root cap, or branch radius;
- a tail ratio outside its convergence region;
- changed D4--D7 constants or pointwise/average sums;
- a lowered `R`, `R'`, commutator, direct-tail, or conversion bound;
- an E9 digest mismatch or malformed exact cubic value; and
- any promoted claim whose strict signed-margin inequality fails.

Tests include exact rational regressions, independent small scalar-series
sums, endpoint power-law checks, analytic scalar-multiplier checks, mutation
attacks on every claim-bearing input, and the full non-slow Issue-128 suite.

## Success criterion

This phase succeeds when a source-regenerating verifier proves

```text
B_log_tail = B_direct + B_dexp
```

for the principal local logarithm and, after exact E9 reduction, proves

```text
B_log_tail
  < distance_to_zero(q5/97^4 + q7/97^6 + q9/97^8).
```

Only that strict inequality authorizes promotion of the fixed-instance
finite-step spectral obstruction.  E9 computation, scheduler completion, or a
small unverified decimal by itself does not.
