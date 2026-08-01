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
5d43135f3b1c3bf66499535868a6a4c20e02ad35d68e31b36792a8d3c5accec0
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

The D4--D7 density caps are recovered exactly from the already verified
schema-v3 D5-integrated main certificate, whose SHA-256 is

```text
ec60a7458f90fa2206cf1f58e05da3ea5c9b99062087c2e5bdd4ee06cff9811f
```

For each degree `j=4,...,7`, the verifier inverts the recorded contribution
at its certified 95-step point by

```text
d_j_site = contribution_j * (j+1) * 95^j / 144.
```

It independently rebuilds the 31 interval stages at the certificate's frozen
12-digit coefficient precision.  This avoids depending on unrelated dirty
working-tree cache or scan additions in the older generator modules.

## Right-generator and logarithm envelope

Every exact value below is stored as a canonical integer pair at the matching
JSON pointer under `/envelope`.  The decimal column is a non-normative
floating summary; the final column records the exact denominator length as a
quick audit of the outward rational computation.

| JSON field | outward decimal summary | denominator digits |
|---|---:|---:|
| `average_generator_defect` | `8.8895058475996365e-7` | 504 |
| `one_step_unitary_defect` | `9.1644390181439548e-9` | 506 |
| `pointwise_generator_defect` | `5.4787498224102065e-6` | 899 |
| `relative_log_defect` | `9.1644391021308979e-9` | 506 |
| `centered_exact_phase_radius` | `1.4845360824742269` | 2 |
| `centered_log_radius_bound` | `1.4845360916386658` | 508 |
| `log_defect` | `1.4663102563409437e-8` | 507 |
| `log_derivative_defect_hs` | `8.9859462543074722e-6` | 1,404 |
| `log_commutator_defect_hs` | `2.6899869701361978e-5` | 1,911 |

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
| `direct_even_generator_tail` | `6.4489120754528513e-12` | 535 |
| `dexp_correction_tail` | `4.6180470539587816e-14` | 2,417 |
| `total_log_tail` | `6.4950925459924386e-12` | 2,447 |

The exact total is the sum of the two exact rational pairs.  The current
E5+E7 signed magnitude at `r=97` is approximately
`6.923687090362501e-11`, so the tail pre-audit consumes about `9.4%` of that
pre-E9 margin.  This comparison is only a feasibility diagnostic.  The exact
E9 contribution can change the final signed margin and must be verified before
promotion.

## Source binding

The artifact binds these implementation SHA-256 values:

| source | SHA-256 |
|---|---|
| `scripts/certify_dual_log_tail.py` | `cc6dd59789589107b5f77b4d61e4257260094a4c6dd10726a6aff286a4fbec2b` |
| `src/trottercert/dual_log_tail.py` | `ac404b6676b20f74b859681ab301beae7b17e3e4b8f54529b8ee6db693d3a61e` |
| `src/trottercert/intervals.py` | `7b81e9123f9643860962047fec9ffa0b38ec5773fb0b9c56a9a436bc370eeee8` |

The largest exact denominator has 2,447 digits, below Python's default
integer-to-text safety limit.  The canonical JSON is 26,793 bytes.

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

The final focused dual-tail run completed with `27 passed` in 1.14 seconds.
A detached, genuine Git worktree at the then-current branch head ran the 27
dual-tail tests plus all seven committed delivery-package tests with `34
passed` in 16.73 seconds.  The shared dirty worktree's non-slow run reached
`298 passed, 16 deselected` and four failures, all in an uncommitted concurrent
`test_publication_scope.py` suite whose ownership manifest did not yet classify
an untracked `CITATION.cff`; those files are outside this certificate change.
Ruff import/typing checks, bytecode compilation, and `git diff --check` passed.
