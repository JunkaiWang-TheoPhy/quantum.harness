# Compatible dual local-log tail ledger

## Scope

This ledger certifies the E11-and-higher dual local-log tail for the frozen
periodic `12 x 12` isotropic Heisenberg benchmark at `r=97`.  The normalized
quantity is

```text
tau(W R_log)/(N/4),
W = H^2 - tau(H^2) I + H/2.
```

The certificate is fixed-instance and principal-branch.  It does not claim a
finite-step theorem uniform in lattice size.  It also does not promote the
finite-step obstruction yet: exact E9 is still missing.

The canonical mathematical record is
`docs/experiments/processor-obstruction/dual-log-tail.json`, whose SHA-256 is

```text
11b2b1595b39d3a4421124adc1e67136c7cb777376406ca368cedfb8abc8ad39
```

## Identity connecting the two proof objects

Write `S(t)=exp(A(t))` and `K_R(t)=S'(t)S(t)^(-1)`.  Time symmetry gives
`A(-t)=-A(t)`, so

```text
K_even(t) = (K_R(t)+K_R(-t))/2
          = sinh(ad_A(t))/ad_A(t) applied to A'(t).
```

Thus

```text
K_even-A' = sum_(k>=1) ad_A^(2k)(A')/(2k+1)!.
```

Because `[W,H]=0`, setting `A=tH+R` makes the first commutator on `W`
equal to `[R,W]`.  Moving one more commutator through the normalized trace
gives the certified Hilbert--Schmidt inequality

```text
|tau(W ad_A^(2k)(A'))|
 <= 2 ||R|| ||W||_2 (2||A||)^(2k-2) ||[A,A']||_2.
```

This is the essential opening: the conversion is quadratic in formula
defects rather than a generic first-order change of norm.

## Exact fixed inputs

The verifier regenerates and checks:

| input | exact value |
|---|---:|
| sites `N` | `144` |
| two-by-two cells `N/4` | `36` |
| steps | `97` |
| centered `H` operator cap | `144` |
| `tau(H^2)` | `54` |
| `H` Hilbert--Schmidt cap | `15/2` |
| `tau(W^2)` | `23085/4` |
| `||W||_2/(N/4)` cap | `17/8` |
| maximum nonzero Pauli coefficient of `W` | `1/4` |
| centered log-radius cap | `3/2` |
| inverse-`dexp` cap | `8/5` |
| multiplier-difference cap | `1` |
| even-series cap | `10/33` |

The exact E5 extensive-witness input has SHA-256

```text
5e5ba831312109df42897c147a85dc4bb77e6670fec4fc6d5befae24bdd8416d
```

and the exact E7 dual-pairing input has SHA-256

```text
9d68771af90096d0ad96bdf6a0861c10d8ec39a864b417d772e51e0c6a9d5d7b
```

## Right-generator and logarithm envelope

Every exact value below is stored as a canonical integer pair at the matching
JSON pointer under `/envelope`.  The decimal column is a non-normative
floating summary; the final column records the exact denominator length as a
quick audit of the outward rational computation.

| JSON field | outward decimal summary | denominator digits |
|---|---:|---:|
| `average_generator_defect` | `2.1469354409384142e-6` | 1,196 |
| `one_step_unitary_defect` | `2.2133355061220766e-8` | 1,198 |
| `pointwise_generator_defect` | `1.1875617201006646e-5` | 2,136 |
| `relative_log_defect` | `2.2133355551106182e-8` | 1,198 |
| `centered_exact_phase_radius` | `1.4845360824742269` | 2 |
| `centered_log_radius_bound` | `1.4845361046075825` | 1,200 |
| `log_defect` | `3.5413368881769893e-8` | 1,199 |
| `log_derivative_defect_hs` | `1.9532188054837184e-5` | 3,333 |
| `log_commutator_defect_hs` | `5.8523677790783468e-5` | 4,531 |

The centered-log value is strictly below `3/2`; hence the artifact proves the
principal-log stability gate used by the `dexp` conversion.

## Compatible E11-and-higher bound

For a stage coefficient `a`, later-stage absolute prefix `p`, and omitted
generator degree `n>=10`, the conjugation Taylor factor, connected local
commutator growth, time average, and dual normalization combine as

```text
(a p^n/n!) * ((3/2)(n+1)!) * (1/(n+1)) * (1/4)
= (3/8) a p^n.
```

Therefore the direct even-generator contribution is

```text
B_direct = (3/8) sum_stages a (p/97)^10/(1-p/97).
```

The old D8+ operator tail was not substituted.  Its scalar stage sum was
rederived with the exact `1/4` Pauli coefficient and per-cell dual
normalization, and the difference between `K_even` and `A'` was certified as a
separate `dexp` correction.

| JSON field | outward decimal summary | denominator digits |
|---|---:|---:|
| `direct_even_generator_tail` | `6.4489120752848012e-12` | 1,263 |
| `dexp_correction_tail` | `2.4265048229124384e-13` | 5,730 |
| `total_log_tail` | `6.6915625575760451e-12` | 5,796 |

The exact total is the sum of the two exact rational pairs.  The current
E5+E7 signed magnitude at `r=97` is approximately
`6.923687090362501e-11`, so the tail pre-audit consumes about `9.7%` of that
pre-E9 margin.  This comparison is only a feasibility diagnostic.  The exact
E9 contribution can change the final signed margin and must be verified before
promotion.

## Source binding

The artifact binds these implementation SHA-256 values:

| source | SHA-256 |
|---|---|
| `scripts/certify_dual_log_tail.py` | `9a5122ecab7e69cc675943f08f38b73ec364cd30f4dfd9c5844a996db5aa1eb0` |
| `src/trottercert/dual_log_tail.py` | `413f8d71a3f7b1b03721e43b3ea3e8da85002138395ed14ca0acc02dc2afa0cf` |
| `src/trottercert/refined_error.py` | `6fd62e1d045de4adfcc407463882fedf5fef30c34b6e793c052ddad62e8f4ebf` |
| `src/trottercert/rigorous_fourth.py` | `e1b626c4464fd9d8a3d4aa1ef88844974abc80f04f14bc01feaf860b1582310b` |

Python's default 4,300-digit integer-to-text safety limit is smaller than the
largest exact denominators in this artifact.  The dedicated generator and
verifier explicitly enable unbounded integer conversion only inside their
bounded, source-regenerated certificate process.  The JSON is 59 KiB; this is
not an unbounded-input parsing path.

## Claim boundary

The canonical claim is exactly

```text
dual_e11_plus_tail = certified
finite_step_status = inconclusive
missing = exact dual E9 pairing
```

No scheduler completion, sampled coefficient, or decimal comparison may
change that status.  Promotion requires a hash-bound exact E9 reduction and
the strict exact-cubic interval gate

```text
B_log_tail
  < distance_to_zero(q5/97^4 + q7/97^6 + q9/97^8).
```

## Verification

From the Issue-128 directory:

```bash
export PYTHONPATH=src:.
python scripts/certify_dual_log_tail.py \
  --verify docs/experiments/processor-obstruction/dual-log-tail.json

pytest -q \
  tests/test_dual_log_tail.py \
  tests/test_dual_log_tail_certificate.py \
  tests/test_refined_error.py \
  tests/test_dual_e9_pairing.py

pytest -q -m 'not slow'
```

The verifier rebuilds the stages, D4--D7 constants, all intermediate
rationals, source hashes, E5/E7 hashes, bounds, and claim.  Lowering a bound or
changing any claim-bearing input causes strict regeneration mismatch.

The final focused run completed with `35 passed, 4 deselected` in 165.41
seconds.  The complete non-slow Issue-128 suite completed with `270 passed, 16
deselected` in 393.56 seconds.  Ruff import/typing checks, bytecode compilation,
and `git diff --check` also passed before publication.
