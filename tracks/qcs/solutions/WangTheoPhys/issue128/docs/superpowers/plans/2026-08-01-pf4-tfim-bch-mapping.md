# Exact PF4--TFIM BCH Mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prove the exact degree-five cyclic free-trace identity for the actual five-copy Suzuki PF4 formula and promote the periodic TFIM artifact only to a rigorously scoped leading-order endpoint-conjugation obstruction.

**Architecture:** A dependency-light cubic free-word module constructs both sides of the trace identity and verifies their ten cyclic degree-six classes exactly.  The trace-obstruction layer exposes the corrected quadratic core and exact cubic Suzuki prefactor; the TFIM layer combines those with existing exact Pauli counts and freezes a schema-v2 certificate.  Claim fields distinguish fixed-normalization endpoint conjugation from affine gauge, finite-step, and total-time eigenphase statements.

**Tech Stack:** Python 3.12, `fractions.Fraction`, exact `Cubic` arithmetic in `Q(alpha)/(alpha^3-4)`, exact free-word series, exact Pauli sums, canonical JSON, SHA-256, pytest, Ruff.

## Global Constraints

- Preserve every unrelated modified or untracked file in the primary worktree.
- Implement in `/tmp/issue128-pf4-bch-20260801` on branch `codex/issue128-pf4-bch`.
- Do not import the currently untracked `rational_words.py` or `processed_kernels.py` from the trusted PF4 mapping module.
- Use the exact stages returned by `fourth_order_suzuki_cubic_stages(2)`; never reconstruct Suzuki coefficients from decimal floats.
- The authoritative identity is equality of all canonical cyclic word classes, not agreement on sampled matrices.
- The matrix tests are independent cross-checks and cannot replace the free-word identity.
- Promote only the leading-order, fixed-time, fixed-normalization obstruction to pure endpoint conjugation.
- Keep affine `a I + b H`, finite-step TFIM, and total-time eigenphase claims explicitly unclaimed.
- The Heisenberg E9 run and its finite-step certificate remain a separate proof chain.
- Every generated JSON artifact must be canonical, source-hash-bound, atomically published, and mutation-tested.
- Every commit message must include `Co-authored-by: OmX <omx@oh-my-codex.dev>`.

---

### Task 1: Exact cubic cyclic free-word identity

**Files:**
- Create: `src/trottercert/pf4_bch_mapping.py`
- Create: `tests/test_pf4_bch_mapping.py`

**Interfaces:**
- Consumes: `Cubic`, `fourth_order_suzuki_cubic_stages`, and `cubic_formula_log_series`.
- Produces: `Word`, `CubicWordPolynomial`, `pf4_suzuki_gamma() -> Cubic`, `pf4_suzuki_trace_polynomials() -> tuple[CubicWordPolynomial, CubicWordPolynomial]`, `cyclic_trace_classes(polynomial: Mapping[Word, Cubic]) -> CubicWordPolynomial`, and `verify_pf4_suzuki_trace_identity() -> None`.

- [ ] **Step 1: Write the exact class-equality test**

```python
from fractions import Fraction

from trottercert.cubic_field import Cubic
from trottercert.pf4_bch_mapping import (
    cyclic_trace_classes,
    pf4_suzuki_gamma,
    pf4_suzuki_trace_polynomials,
    verify_pf4_suzuki_trace_identity,
)


def test_pf4_suzuki_trace_identity_matches_all_cyclic_classes() -> None:
    left, right = pf4_suzuki_trace_polynomials()
    left_classes = cyclic_trace_classes(left)
    right_classes = cyclic_trace_classes(right)

    assert len(left_classes) == 10
    assert left_classes == right_classes
    assert pf4_suzuki_gamma() == Cubic(
        Fraction(37, 900000),
        Fraction(313, 14400000),
        Fraction(29, 1800000),
    )
    verify_pf4_suzuki_trace_identity()
```

- [ ] **Step 2: Run the test and verify the missing-module failure**

Run:

```bash
cd tracks/qcs/solutions/WangTheoPhys/issue128
PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py::test_pf4_suzuki_trace_identity_matches_all_cyclic_classes
```

Expected: collection fails with `ModuleNotFoundError: trottercert.pf4_bch_mapping`.

- [ ] **Step 3: Implement exact polynomial primitives and the two sides**

Create `pf4_bch_mapping.py` with the following public structure and exact constants:

```python
from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction

from .cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from .cubic_local import cubic_formula_log_series

Word = tuple[int, ...]
CubicWordPolynomial = dict[Word, Cubic]


def pf4_suzuki_gamma() -> Cubic:
    return Cubic(
        Fraction(37, 900000),
        Fraction(313, 14400000),
        Fraction(29, 1800000),
    )


def _add(left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]) -> CubicWordPolynomial:
    result = dict(left)
    for word, coefficient in right.items():
        updated = result.get(word, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            result.pop(word, None)
        else:
            result[word] = updated
    return result


def _scale(polynomial: Mapping[Word, Cubic], scalar: Cubic | Fraction | int) -> CubicWordPolynomial:
    factor = Cubic.coerce(scalar)
    return {
        word: coefficient * factor
        for word, coefficient in polynomial.items()
        if coefficient * factor != Cubic.zero()
    }


def _multiply(left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]) -> CubicWordPolynomial:
    result: CubicWordPolynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            result = _add(
                result,
                {word: left_coefficient * right_coefficient},
            )
    return result


def _commutator(left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]) -> CubicWordPolynomial:
    return _add(_multiply(left, right), _scale(_multiply(right, left), -1))


def cyclic_trace_classes(polynomial: Mapping[Word, Cubic]) -> CubicWordPolynomial:
    classes: CubicWordPolynomial = {}
    for word, coefficient in polynomial.items():
        if not word:
            representative = ()
        else:
            representative = min(
                word[offset:] + word[:offset] for offset in range(len(word))
            )
        classes = _add(classes, {representative: coefficient})
    return classes


def pf4_suzuki_trace_polynomials() -> tuple[CubicWordPolynomial, CubicWordPolynomial]:
    a = {(0,): Cubic.one()}
    b = {(1,): Cubic.one()}
    hamiltonian = _add(a, b)
    logarithm = cubic_formula_log_series(
        fourth_order_suzuki_cubic_stages(2), 5
    )
    if logarithm[1] != hamiltonian or any(logarithm[d] for d in (2, 3, 4)):
        raise ArithmeticError("exact Suzuki stages are not fourth order")
    e5 = logarithm[5]
    c = _commutator(a, _commutator(a, b))
    d = _commutator(b, _commutator(b, a))
    quadratic = _add(
        _multiply(c, c),
        _add(
            _scale(_multiply(c, d), -4),
            _scale(_multiply(d, d), Fraction(8, 3)),
        ),
    )
    return _multiply(hamiltonian, e5), _scale(
        quadratic, pf4_suzuki_gamma()
    )


def verify_pf4_suzuki_trace_identity() -> None:
    left, right = pf4_suzuki_trace_polynomials()
    left_classes = cyclic_trace_classes(left)
    right_classes = cyclic_trace_classes(right)
    if len(left_classes) != 10 or left_classes != right_classes:
        raise ArithmeticError("exact PF4 cyclic free-trace identity failed")
```

- [ ] **Step 4: Add exact mutation and fourth-order tests**

Add tests that independently rebuild the right side with one coefficient changed and require class inequality:

```python
@pytest.mark.parametrize(
    ("c2", "cd", "d2"),
    (
        (Fraction(2), Fraction(-4), Fraction(8, 3)),
        (Fraction(1), Fraction(-3), Fraction(8, 3)),
        (Fraction(1), Fraction(-4), Fraction(7, 3)),
    ),
)
def test_pf4_trace_identity_rejects_mutated_quadratic_coefficients(
    c2: Fraction, cd: Fraction, d2: Fraction
) -> None:
    assert (c2, cd, d2) != (Fraction(1), Fraction(-4), Fraction(8, 3))
```

Implement the test helper with its own short word-polynomial construction so it does not call the right-side builder under test.  Also assert the exact logarithm degree maps at degrees two, three, and four are empty.

- [ ] **Step 5: Add exact rational matrix cross-checks**

For seeded symmetric integer matrices in dimensions two and three, evaluate each raw word polynomial with `dtype=object`, take the ordinary trace, and compare the three `Cubic` coordinates exactly:

```python
@pytest.mark.parametrize("dimension", (2, 3))
def test_pf4_trace_identity_on_exact_symmetric_matrices(dimension: int) -> None:
    rng = np.random.default_rng(12840 + dimension)
    raw_a = rng.integers(-3, 4, size=(dimension, dimension))
    raw_b = rng.integers(-3, 4, size=(dimension, dimension))
    a = np.asarray(raw_a + raw_a.T, dtype=object)
    b = np.asarray(raw_b + raw_b.T, dtype=object)
    left, right = pf4_suzuki_trace_polynomials()

    assert _evaluate_trace(left, a, b) == _evaluate_trace(right, a, b)
```

The local `_evaluate_trace` helper must start with an exact identity matrix, multiply matrices according to each word, and accumulate every matrix entry as a `Cubic` value without converting to float.

- [ ] **Step 6: Run and lint the new algebraic kernel**

Run:

```bash
PYTHONPATH=src:. python -m pytest -q tests/test_pf4_bch_mapping.py
python -m ruff check src/trottercert/pf4_bch_mapping.py \
  tests/test_pf4_bch_mapping.py
```

Expected: all tests pass and Ruff reports `All checks passed!`.

- [ ] **Step 7: Commit the exact identity kernel**

```bash
git add src/trottercert/pf4_bch_mapping.py tests/test_pf4_bch_mapping.py
git commit -m "feat(issue128): prove exact Suzuki PF4 trace identity" \
  -m "Co-authored-by: OmX <omx@oh-my-codex.dev>"
```

### Task 2: Correct and bind the canonical PF4 trace record

**Files:**
- Modify: `src/trottercert/trace_obstruction.py`
- Modify: `tests/test_trace_obstruction.py`
- Test: `tests/test_pf4_bch_mapping.py`

**Interfaces:**
- Consumes: `pf4_suzuki_gamma()` and `verify_pf4_suzuki_trace_identity()` from Task 1.
- Produces: corrected `pf4_trace_quadratic_form(...) -> Fraction`, new `pf4_suzuki_trace_pairing(...) -> Cubic`, and canonical `pf4_trace_identity_record()` with status `proved_exact_suzuki_pf4_free_trace_identity`.

- [ ] **Step 1: Change the tests to the corrected core and proved status**

Replace the old coefficient expectation with:

```python
def test_pf4_quadratic_form_preserves_exact_derived_coefficients() -> None:
    value = pf4_trace_quadratic_form(6, 5, 3, Fraction(7, 11))
    expected = Fraction(7, 11) * (
        Fraction(6) - 4 * 5 + Fraction(8, 3) * 3
    )
    assert value == expected
```

Add:

```python
def test_pf4_exact_suzuki_pairing_includes_cubic_gamma() -> None:
    core = pf4_trace_quadratic_form(6, 0, 3, 1)
    assert pf4_suzuki_trace_pairing(6, 0, 3) == pf4_suzuki_gamma() * core
```

Update the canonical-record assertions to require coefficients `1`, `-4`,
`8/3`, the three gamma coordinates, and status
`proved_exact_suzuki_pf4_free_trace_identity`.

- [ ] **Step 2: Run the new standalone mapping tests and record the known collection boundary**

Run:

```bash
PYTHONPATH=src:. python -m pytest -q tests/test_pf4_bch_mapping.py
PYTHONPATH=src:. python -m pytest -q tests/test_trace_obstruction.py
```

Expected: the mapping suite passes.  At current clean HEAD the second command may fail during collection because the already committed test imports the still-untracked concurrent modules `processed_kernels.py` and `rational_words.py`; record that as an external repository-boundary failure, not as evidence for or against this task.  Re-run it after those files are committed by their owner.

- [ ] **Step 3: Implement the corrected evaluator and exact pairing**

In `trace_obstruction.py`, use:

```python
def pf4_trace_quadratic_form(
    trace_c2: Fraction | int,
    trace_cd: Fraction | int,
    trace_d2: Fraction | int,
    gamma: Fraction | int,
) -> Fraction:
    c2 = _exact_fraction(trace_c2, "trace_c2")
    cd = _exact_fraction(trace_cd, "trace_cd")
    d2 = _exact_fraction(trace_d2, "trace_d2")
    scale = _exact_fraction(gamma, "gamma")
    return scale * (c2 - 4 * cd + Fraction(8, 3) * d2)


def pf4_suzuki_trace_pairing(
    trace_c2: Fraction | int,
    trace_cd: Fraction | int,
    trace_d2: Fraction | int,
) -> Cubic:
    verify_pf4_suzuki_trace_identity()
    return pf4_suzuki_gamma() * pf4_trace_quadratic_form(
        trace_c2, trace_cd, trace_d2, 1
    )
```

Import `Cubic`, `pf4_suzuki_gamma`, and
`verify_pf4_suzuki_trace_identity` from committed modules only.

- [ ] **Step 4: Replace the PF4 canonical record**

Build the record with:

```python
formula=(
    "Tr((A+B)E5) = gamma*(Tr(C^2) - 4*Tr(CD) "
    "+ (8/3)*Tr(D^2))"
),
assumptions=(
    "A and B are finite-dimensional Hermitian matrices",
    "E5 is the degree-five logarithm term of the exact five-copy Suzuki PF4 formula",
    "C=[A,[A,B]] and D=[B,[B,A]]",
    "trace is cyclic",
    "the identity is verified in the exact cubic coefficient field alpha^3=4",
),
exact_coefficients=(
    ("trace_c2", Fraction(1)),
    ("trace_cd", Fraction(-4)),
    ("trace_d2", Fraction(8, 3)),
    ("gamma_a0", Fraction(37, 900000)),
    ("gamma_a1", Fraction(313, 14400000)),
    ("gamma_a2", Fraction(29, 1800000)),
),
identity_status="proved_exact_suzuki_pf4_free_trace_identity",
```

Call `verify_pf4_suzuki_trace_identity()` before returning the record so a
future drift in the exact stages fails closed.

- [ ] **Step 5: Run focused tests and mutation checks**

Run:

```bash
PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py \
  tests/test_tfim_obstruction.py
python -m ruff check src/trottercert/trace_obstruction.py \
  src/trottercert/pf4_bch_mapping.py \
  tests/test_pf4_bch_mapping.py
```

Expected: all runnable focused tests pass; no code path retains the old
`1/2, 14/3, 4/3` coefficient set.

- [ ] **Step 6: Commit the corrected trace record**

```bash
git add src/trottercert/trace_obstruction.py tests/test_trace_obstruction.py
git commit -m "fix(issue128): bind PF4 record to exact BCH identity" \
  -m "Co-authored-by: OmX <omx@oh-my-codex.dev>"
```

### Task 3: Promote the narrowly scoped periodic TFIM family certificate

**Files:**
- Modify: `src/trottercert/tfim_obstruction.py`
- Modify: `scripts/certify_tfim_obstruction.py`
- Modify: `tests/test_tfim_obstruction.py`
- Modify: `docs/experiments/processor-obstruction/tfim-family-obstruction.json`

**Interfaces:**
- Consumes: corrected rational core and cubic Suzuki pairing from Task 2.
- Produces: rational keys `trace_c2_over_d`, `trace_d2_over_d`, `trace_cd_over_d`, `trace_h2_over_d`, `quadratic_core_over_d`, `hamiltonian_norm_upper`; cubic keys `trace_h_e5_over_d`, `endpoint_conjugation_lower_bound`; canonical schema-v2 artifact and strict verifier.

- [ ] **Step 1: Write the corrected TFIM formula tests**

For lengths 4, 6, and 8 require:

```python
core = (
    moments["trace_c2_over_d"]
    - 4 * moments["trace_cd_over_d"]
    + Fraction(8, 3) * moments["trace_d2_over_d"]
)
assert moments["quadratic_core_over_d"] == core
assert moments["trace_h_e5_over_d"] == pf4_suzuki_gamma() * core
assert moments["hamiltonian_norm_upper"] == length * (abs(h) + abs(j))
assert moments["endpoint_conjugation_lower_bound"] == (
    pf4_suzuki_gamma() * core / moments["hamiltonian_norm_upper"]
)
```

Add parameterized rational couplings `(2,3)`, `(-2,3)`, `(2,-3)`, `(0,3)`,
`(2,0)`, and `(0,0)`.  Require a strictly positive cubic lower-bound
coordinate triple exactly when `h*j != 0`; require zero in every commuting
case.

- [ ] **Step 2: Run the tests against the old schema**

Run:

```bash
PYTHONPATH=src:. python -m pytest -q tests/test_tfim_obstruction.py
```

Expected: failures show the old quadratic coefficients, missing cubic fields,
and conditional claim status.

- [ ] **Step 3: Extend the TFIM exact result schema**

In `tfim_obstruction.py`, define:

```python
RATIONAL_RESULT_KEYS = (
    "trace_c2_over_d",
    "trace_d2_over_d",
    "trace_cd_over_d",
    "trace_h2_over_d",
    "quadratic_core_over_d",
    "hamiltonian_norm_upper",
)
CUBIC_RESULT_KEYS = (
    "trace_h_e5_over_d",
    "endpoint_conjugation_lower_bound",
)
```

Return a `dict[str, Fraction | Cubic]`.  Compute the exact pairing with
`pf4_suzuki_trace_pairing`.  For nonzero `hamiltonian_norm_upper`, set the
lower bound to `pf4_suzuki_gamma() * abs(core) / norm_upper`; for the zero
Hamiltonian require `core == 0` and return `Cubic.zero()`.

- [ ] **Step 4: Upgrade the artifact contract to schema version two**

In `certify_tfim_obstruction.py`:

- add `src/trottercert/pf4_bch_mapping.py` and
  `src/trottercert/cubic_local.py` to `SOURCE_PATHS`;
- serialize rational results as `[numerator, denominator]`;
- serialize cubic results as three canonical rational pairs;
- require `verify_pf4_suzuki_trace_identity()` during build and verify;
- set the family domain to
  `even length >= 4; exact rational h,j; strict endpoint-conjugation obstruction iff h*j != 0; commuting exception otherwise`;
- freeze formulas
  `quadratic_core_over_d = trace_c2_over_d - 4*trace_cd_over_d + 8*trace_d2_over_d/3` and
  `trace_h_e5_over_d = gamma*quadratic_core_over_d`;
- set `schema_version` to `2`;
- use the following exact claim object:

```python
EXPECTED_CLAIM = {
    "trace_moment_status": "certified_exact_pauli_counting",
    "pf4_bch_mapping_status": "proved_exact_suzuki_pf4_free_trace_identity",
    "endpoint_conjugation_status": (
        "certified_leading_order_fixed_time_fixed_normalization"
    ),
    "operator_lower_bound_status": "certified_exact_algebraic_leading_coefficient",
    "affine_gauge_status": "not_claimed",
    "finite_step_no_go": "not_claimed",
    "total_time_eigenphase_status": "not_claimed",
    "promotion_status": "promoted_scoped_family_theorem",
}
```

- [ ] **Step 5: Add fail-closed artifact mutation tests**

Mutate and reject each of the following independently:

- `schema_version` from `2` to `1`;
- quadratic coefficient `trace_cd` from `[-4,1]` to `[14,3]`;
- `gamma_coordinates[2]`;
- one numerator in `trace_h_e5_over_d`;
- `affine_gauge_status` to `proved`;
- `finite_step_no_go` to `proved`;
- `total_time_eigenphase_status` to `proved`;
- one source SHA;
- canonical JSON whitespace and a duplicate key.

Every mutation must fail for its semantic gate before payload digest mismatch
where the verifier has a dedicated check.

- [ ] **Step 6: Regenerate and verify the canonical artifact**

Run:

```bash
artifact=docs/experiments/processor-obstruction/tfim-family-obstruction.json
temporary=$(mktemp -d)/tfim-family-obstruction.json
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py --build "$temporary"
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py --verify "$temporary"
cp "$temporary" "$artifact"
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py --verify "$artifact"
```

Expected: both verifications print `verified=...`; the tracked artifact is a
single canonical JSON line with schema version two.

- [ ] **Step 7: Run the complete runnable PF4/TFIM gate**

Run:

```bash
PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py \
  tests/test_tfim_obstruction.py
python -m ruff check \
  src/trottercert/pf4_bch_mapping.py \
  src/trottercert/trace_obstruction.py \
  src/trottercert/tfim_obstruction.py \
  scripts/certify_tfim_obstruction.py \
  tests/test_pf4_bch_mapping.py \
  tests/test_tfim_obstruction.py
```

Expected: all tests pass and Ruff reports no findings.

- [ ] **Step 8: Commit the scoped family certificate**

```bash
git add src/trottercert/tfim_obstruction.py \
  scripts/certify_tfim_obstruction.py \
  tests/test_tfim_obstruction.py \
  docs/experiments/processor-obstruction/tfim-family-obstruction.json
git commit -m "feat(issue128): certify scoped PF4 TFIM obstruction" \
  -m "Co-authored-by: OmX <omx@oh-my-codex.dev>"
```

### Task 4: Correct publication text and verify the release boundary

**Files:**
- Modify: `docs/plans/2026-08-01-balanced-dual-track-publication-design.md`
- Modify: `docs/experiments/processor-obstruction/README.md`
- Create: `docs/report/pf4-tfim-bch-ledger.md`
- Modify: `artifacts/publication/paper-a-file-ownership.json`
- Test: `tests/test_publication_scope.py`

**Interfaces:**
- Consumes: the verified schema-v2 TFIM artifact and canonical PF4 identity record.
- Produces: a human-readable theorem ledger and an ownership classification that keeps all PF4/TFIM files in Paper B and out of the Paper A release.

- [ ] **Step 1: Replace the old PF4 coefficients in publication plans**

Search:

```bash
rg -n "14/3|14\\s*\\*.*CD|algebraic_form_only_unverified_bch_mapping|blocked_pending_pf4_bch_mapping" \
  docs src scripts tests
```

For every theorem or plan statement, replace the old unverified form with
`C^2 - 4 C D + (8/3) D^2`, cite the exact gamma coordinates, and label the
result `leading-order endpoint-conjugation obstruction at fixed time and fixed normalization`.
Historical plan text may say the old form was superseded, but it must not
continue presenting those coefficients as valid.

- [ ] **Step 2: Write the theorem ledger with exact scope**

Create `docs/report/pf4-tfim-bch-ledger.md` containing:

1. the exact stage convention and coefficient field;
2. the ten-class cyclic free-word equality;
3. the gamma coordinates and positivity argument;
4. the TFIM Pauli-counting formulas and closed family expression;
5. the duality inequality excluding `E5=-i[Q,H]`;
6. the exact artifact path, payload SHA, verifier command, and source hashes;
7. explicit nonclaims for affine gauge, finite step, and total-time eigenphases;
8. a note that the former `1/2,14/3,4/3` plan coefficients were never certified and are superseded.

Do not insert a decimal-only theorem coefficient; any decimal may appear only
as a clearly marked readability diagnostic beside the exact cubic value.

- [ ] **Step 3: Classify every new or changed Paper B file**

In `paper-a-file-ownership.json`, add or update exact entries with
`{"ownership":"paper-b","release_included":false}` for:

```text
docs/plans/2026-08-01-pf4-tfim-bch-mapping-design.md
docs/superpowers/plans/2026-08-01-pf4-tfim-bch-mapping.md
docs/report/pf4-tfim-bch-ledger.md
src/trottercert/pf4_bch_mapping.py
tests/test_pf4_bch_mapping.py
src/trottercert/trace_obstruction.py
tests/test_trace_obstruction.py
src/trottercert/tfim_obstruction.py
scripts/certify_tfim_obstruction.py
tests/test_tfim_obstruction.py
docs/experiments/processor-obstruction/tfim-family-obstruction.json
```

Preserve every existing ownership entry and canonical key order.

- [ ] **Step 4: Run focused verification and ownership audit**

Run:

```bash
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py \
  --verify docs/experiments/processor-obstruction/tfim-family-obstruction.json
PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py \
  tests/test_tfim_obstruction.py
PYTHONPATH=src:. python scripts/audit_publication_scope.py \
  artifacts/publication/paper-a-file-ownership.json
```

Expected: the certificate and focused tests pass.  If the ownership audit
reports paths belonging to still-untracked concurrent Paper A work, preserve
those files and record the external boundary; do not delete or adopt them.

- [ ] **Step 5: Perform a clean-checkout gate**

Create a detached worktree from the final PF4 commit and run:

```bash
PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py \
  tests/test_tfim_obstruction.py
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py \
  --verify docs/experiments/processor-obstruction/tfim-family-obstruction.json
```

Expected: both commands pass without access to the dirty primary worktree.  A
repository-wide test run is required at the final Issue-128 completion audit
after the owners of `reference_verify.py`, `processed_kernels.py`, and
`rational_words.py` commit their corresponding source boundary.

- [ ] **Step 6: Commit publication documentation separately**

```bash
git add docs/plans/2026-08-01-balanced-dual-track-publication-design.md \
  docs/experiments/processor-obstruction/README.md \
  docs/report/pf4-tfim-bch-ledger.md
git add -f docs/superpowers/plans/2026-08-01-pf4-tfim-bch-mapping.md
git commit -m "docs(issue128): state scoped PF4 TFIM theorem" \
  -m "Co-authored-by: OmX <omx@oh-my-codex.dev>"
```

Commit the ownership manifest only after confirming it does not conflict with
newer primary-worktree edits:

```bash
git add artifacts/publication/paper-a-file-ownership.json
git commit -m "chore(issue128): classify PF4 TFIM proof files" \
  -m "Co-authored-by: OmX <omx@oh-my-codex.dev>"
```

## Self-review

- **Spec coverage:** Tasks 1 and 2 prove and bind the exact BCH/free-trace identity; Task 3 specializes it to the periodic TFIM family and freezes the scoped certificate; Task 4 corrects publication text, records nonclaims, classifies ownership, and supplies clean-checkout evidence.
- **Placeholder scan:** The plan contains no `TBD`, `TODO`, deferred implementation, or unnamed test instruction.  The only deferred repository-wide gate names the exact untracked concurrent dependencies and is part of the final Issue-128 audit, not this theorem's focused acceptance test.
- **Type consistency:** `pf4_suzuki_gamma()` and `pf4_suzuki_trace_pairing()` return `Cubic`; rational moment/core fields return `Fraction`; schema-v2 serialization distinguishes the two; every later task uses the same names.
- **Claim consistency:** No task promotes affine gauge, finite-step TFIM, or total-time eigenphase claims.  The E9 finite-step proof remains separate.

## Execution choice

The user selected inline execution with “做吧”.  Proactive subagents are
disabled by the active collaboration policy, and the optional
`superpowers:executing-plans` skill is not installed in this session.  Execute
the checklist directly in the isolated worktree with TDD and commit after each
reviewable task.
