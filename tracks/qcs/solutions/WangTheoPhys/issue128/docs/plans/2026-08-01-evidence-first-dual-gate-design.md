# Issue 128 evidence-first dual-gate design

**Date:** 2026-08-01
**Route:** Paper A Quantum completion plus Paper B PRX Quantum sprint
**Status:** Route approved; exact compute setup awaiting final ratification

## 1. Decision

Run three scientific lanes in parallel, but promote their claims independently:

1. certify or reject the strengthened 353-step control used in the informal
   `10591/2851` comparison;
2. extend the frozen Heisenberg compiler to a preregistered XXZ family and seek
   a positive certified transfer away from the isotropic point;
3. replace the unverified PF4 quadratic form with a free-trace identity derived
   from the actual five-copy Suzuki formula.

No lane may borrow another lane's status.  In particular, a plausible 353-step
number does not become a baseline certificate, an unsupported XXZ row does not
become transfer evidence, and exact TFIM Pauli moments do not become a PF4
operator obstruction without the BCH-to-trace mapping.

## 2. Frozen scientific setup

### 2.1 Paper A flagship and strengthened control

Use the existing flagship convention

```text
H = sum_<i,j> (X_i X_j + Y_i Y_j + Z_i Z_j) / 4
lattice = 12 x 12 periodic square lattice
formula = five-copy fourth-order Suzuki formula
target time T = 1
operator-norm tolerance epsilon = 10^-6
resource model G(r) = 30 r + 1 merged group exponentials
```

The current certified candidate is `r=95`, `G=2851`.  The strengthened control
claim is exactly `r=353`, `G=10591`; both numbers must be derived from a complete
formula evaluation and independently verified.  Until then it remains
`missing_evidence`, and the only formal comparison is `11791/2851` against the
published `r=393` baseline.

### 2.2 Paper A XXZ transfer

Use a staged finite-lattice program.  The compiler-validation pilot is

```text
H(Delta) = sum_<i,j> (X_i X_j + Y_i Y_j + Delta Z_i Z_j) / 4
lattice = 4 x 4 periodic square lattice (direct finite-torus evaluation)
formula, T, epsilon, and resource model = the flagship values above
Delta grid = 0, 1/4, 1/2, 1, 3/2, 2, 4
```

The first two nonisotropic pilot points are `Delta=1/2` and `Delta=2`, one
easy-plane and one easy-axis instance.  The 4 by 4 implementation must evaluate
the finite torus directly or reject detected wraparound aliasing; it may not
multiply an infinite-lattice or unit-cell density by 16.  After the direct
finite implementation and verifier pass, the publication transfer run scales
the same frozen rules to the 12 by 12 flagship geometry.

The compiler and all acceptance rules are frozen before evaluating any new
anisotropy.  `Delta=1` is an identity reduction to the existing Heisenberg
certificate, not independent transfer evidence.  A positive transfer result
requires at least one `Delta != 1` row with a complete certificate accepted by
the independent verifier.  Unsupported and inconclusive rows remain part of
the frozen release.

### 2.3 Paper B PF4 mapping and TFIM family

Use the algebraic convention

```text
S2(t) = exp(t A/2) exp(t B) exp(t A/2)
S4(t) = S2(u t)^2 S2((1-4u)t) S2(u t)^2
u = 1 / (4 - alpha),  alpha^3 = 4
log S4(t) = t(A+B) + t^5 E5 + O(t^7)
C = [A,[A,B]]
D = [B,[B,A]]
```

The historical free-word discovery suggests

```text
Tr((A+B) E5) = gamma Tr(C^2 - 4 C D + (8/3) D^2),
```

not the currently stored plan form with coefficients `1/2, 14/3, 4/3`.
The new coefficients and the cubic-field scalar `gamma` are hypotheses until
derived from the exact Suzuki word and independently checked.  The derivation
must solve for all cyclic coefficients without presupposing either the new or
legacy vector; the legacy vector is retained only as a negative regression.
Because the TFIM family has `Tr(CD)=0` and the two diagonal vectors differ only
by an overall factor, TFIM instances alone cannot distinguish the two formulas.

The target TFIM family is

```text
A = h sum_i X_i
B = j sum_i Z_i Z_(i+1)
boundary = periodic
domain = even length L >= 4 and exact rational h,j
```

The already certified identities

```text
Tr(C^2)/d = 128 L h^4 j^2,
Tr(CD)/d = 0,
Tr(D^2)/d = 128 L h^2 j^4
```

remain valid inputs.  They become an operator obstruction only after the exact
PF4 mapping, the sign/nonvanishing of `gamma`, and the gauge implication are
proved.

## 3. Evidence and verification architecture

### Lane A1: strengthened control

- Reconstruct the complete published/control theorem evaluation from frozen
  Hamiltonian, Suzuki stages, theorem center, norm rule, and cost model.
- Produce canonical JSON with exact rationals, source hashes, adjacent-step
  evidence, and the derived `r=353`, `G=10591` only if they follow.
- Verify it with a second implementation that does not trust submitted totals.
- If the exact result differs, publish the derived result and retire 10591.

### Lane A2: grouped XXZ

- Add a parameterized exact Hamiltonian and fragment identity.
- On the 4 by 4 torus, prove that the four checkerboard matchings contain eight
  disjoint bonds each and cover all 32 periodic nearest-neighbor bonds.
- Separate discovery grouping from the trusted certificate checker.
- Certify Pauli coverage, coefficients, translation multiplicities,
  anticommutation, finite-step ledger, and resources for each supported Delta.
- Promote merged group exponentials first; retain CNOT cost as missing evidence
  until a general XXZ bond-synthesis certificate exists.
- Freeze all failures as well as successes.

#### Production-ledger execution design

The direct finite-torus theorem is unchanged, but its production witness uses
an exact D4 orbit compression so that verification remains practical.  The
raw theorem stream is enumerated in full and hashed canonically.  Every
actual five-letter word retains its own exact positive theorem weight; weights
are never inferred from an orbit representative and are never combined with a
Pauli coefficient map.  Only the unweighted commutator polynomial is evaluated
for the lexicographically canonical representative of each frozen 4 by 4 D4
orbit.  The frozen site permutation transports Pauli masks and proves equality
of coefficient magnitudes and commutation graphs for the other orbit members.

Representative blocks use `deterministic_pair_only_bitset_v1`: terms are
ordered by decreasing exact magnitude and mask, and each unmatched term is
paired with the earliest unmatched anticommuting partner.  Groups are therefore
singletons or pairs and can be independently replayed without trusting a graph
heuristic.  The theorem constants are still evaluated as the literal sums

```text
K_grouped = sum_actual_words w_k U_rep(k) / 5!
K_triangle = sum_actual_words w_k ||C_rep(k)||_1 / 5!.
```

This is the selected production path.  Evaluating pair-only groups for every
actual block is a correct but slower fallback.  Multi-member greedy grouping is
reserved for a second pass over the canonical representatives only when the
pair-only adjacent-step result does not improve the same-Delta baseline enough
to pass the preregistered pilot gate.

A bounded prefix is profiling evidence only.  It cannot produce a certified
finite-step result.  Production status requires full raw-stream coverage,
complete actual-word-to-representative coverage, exact representative maps,
deterministic group replay, strict accepted/previous-step inequalities, and an
independent standard-library verifier.

### Lane B: PF4 BCH mapping

- Implementation 1 expands the exact cubic-field Suzuki word through degree
  five and reduces `Tr((A+B)E5)` modulo cyclic trace identities.
- Implementation 2 independently expands seeded Hermitian matrices and solves
  for the trace polynomial coefficients without importing implementation 1.
- Exact small rational matrices, commuting cases, `A=0`, `B=0`, swap/order
  conventions, and sign changes are mandatory adversarial fixtures.
- Only a canonical identity record with both implementations agreeing can
  replace `algebraic_form_only_unverified_bch_mapping`.

## 4. Promotion gates

### Paper A strengthened-baseline gate

Pass only if the full authoritative calculation and independent verifier both
derive `353` and `10591`.  Otherwise retain `complete_formula=missing_evidence`.

### Paper A positive-transfer gate

Pass only if at least one preregistered `Delta != 1` row is certified without
post-hoc rule changes.  The `Delta=1` identity row alone does not pass this gate.

### Paper B PF4-family gate

Pass only if all of the following hold:

- the actual five-copy Suzuki BCH/free-trace mapping is exact;
- the cubic-field prefactor is proved nonzero with the needed sign;
- the TFIM substitution gives a strict family statement on its declared
  domain;
- the physical `A,B -> -iA,-iB` sign convention is bridged explicitly;
- the gauge/correctability implication is stated and proved separately with a
  witness orthogonal to both `I` and `H`; pairing with `H` alone establishes
  only a fixed-target obstruction and cannot survive the allowed retiming
  gauge;
- mutation tests cannot promote the old conditional artifact by rehashing it.

Finite-step PRX Quantum promotion remains a separate gate and cannot follow
from a leading-order PF4 identity alone.

## 5. Execution and ownership

The three lanes use disjoint file allowlists and separate commits.  Paper A
artifacts enter the Paper A release only after their gates pass.  PF4/TFIM
artifacts remain Paper B-owned.  Existing dirty manuscript, processed-kernel,
and HPC files are outside these edits and must not be staged incidentally.
