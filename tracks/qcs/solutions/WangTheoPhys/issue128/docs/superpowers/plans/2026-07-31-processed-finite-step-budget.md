# Processed Finite-Step Budget Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the exact processed-kernel free-Lie audit into rigorous local Pauli bounds for rationalization residuals and the proof-only processor, then expose a fail-closed finite-step budget contract for future physical E7 and all-order remainder certificates.

**Architecture:** Apply the Dynkin--Specht--Wever projection to exact rational word polynomials using the existing four-matching symplectic local evaluator. Keep the already-computable degree-three/five residual and processor bounds separate from the not-yet-computed physical degree-seven coefficient map and all-order remainder. A budget is accepted only when every required component is supplied and the sum is below tolerance.

**Tech Stack:** Python 3.12, `Fraction`, existing `SymplecticDyadicLocalDensityEvaluator`, exact rational exponential enclosure, pytest 8.

## Global Constraints

- Preserve the dirty worktree and add focused files only.
- Frozen benchmark: periodic `12 x 12` isotropic Heisenberg, `N = 144`, `T = 1`, tolerance `10^-6`, four matching fragments.
- Primary integer gates: `s10` at `r = 47` and `39`; `s11` at `r = 43` and `35`.
- Every physical coefficient and bound is an exact `Fraction`; floating point may appear only in printed diagnostics.
- The degree-seven physical site density and the all-order remainder are required external proof inputs. Neither may default to zero.
- Do not infer a global logarithm branch from a formal BCH series. The budget object records the local-log theorem assumption explicitly and cannot by itself certify it.
- Do not prepare or submit HPC jobs in this plan.

---

### Task 1: Exact Free-Lie to Local-Pauli Projection

**Files:**
- Create: `src/trottercert/rational_lie_local.py`
- Test: `tests/test_rational_lie_local.py`

**Interfaces:**
- Consumes: `WordPolynomial` and the existing symplectic local evaluator.
- Produces: `RationalLocalLieTerms`, `project_lie_polynomial_to_cell(polynomial)`, and `cell_l1(polynomial)`.

- [x] **Step 1: Write failing projection tests**

```python
from fractions import Fraction

from trottercert.rational_lie_local import project_lie_polynomial_to_cell


def test_zero_and_commutator_projection() -> None:
    assert not project_lie_polynomial_to_cell({})
    projected = project_lie_polynomial_to_cell(
        {(0, 1): Fraction(1), (1, 0): Fraction(-1)}
    )
    assert projected
    assert sum(abs(value) for value in projected.values()) > 0


def test_projection_is_exactly_linear() -> None:
    left = {(0, 1): Fraction(1), (1, 0): Fraction(-1)}
    right = {(2, 3): Fraction(2), (3, 2): Fraction(-2)}
    combined = {**left, **right}
    assert project_lie_polynomial_to_cell(combined) == {
        **project_lie_polynomial_to_cell(left),
        **project_lie_polynomial_to_cell(right),
    }
```

Use disjoint matching-color supports in the linearity test; if canonical Pauli keys overlap, compare by exact coefficient addition instead of dictionary unpacking.

- [x] **Step 2: Run the test and confirm missing-module failure**

Run: `python -m pytest -q tests/test_rational_lie_local.py`

Expected: collection fails with `ModuleNotFoundError`.

- [x] **Step 3: Implement the exact DSW projection**

Require every nonzero word to have one common positive degree.  For each word,
use the existing evaluator's right-nested commutator and multiply by

```text
word_coefficient / (degree * 2^denominator_exponent).
```

Canonicalize every symplectic Pauli term to the shared `2 x 2` unit cell with
`canonicalize_symplectic_unit_cell`, combine equal terms, and remove exact
zeros.

- [x] **Step 4: Run the projection tests**

Run: `python -m pytest -q tests/test_rational_lie_local.py`

Expected: all tests pass in under one second.

### Task 2: Residual and Processor Local Bounds

**Files:**
- Create: `src/trottercert/processed_local_bounds.py`
- Test: `tests/test_processed_local_bounds.py`

**Interfaces:**
- Consumes: `EffectiveOrderAudit` and Task 1 projection.
- Produces: `ProcessedLocalBounds`, `build_processed_local_bounds(kernel)`, `rational_exp_upper(x, terms=32)`, and `processor_distance_upper(bounds, n_sites, steps)`.

- [x] **Step 1: Write failing exact-bound tests**

```python
from fractions import Fraction

from trottercert.processed_kernels import published_effective_order_six_kernel
from trottercert.processed_local_bounds import (
    build_processed_local_bounds,
    processor_distance_upper,
    rational_exp_upper,
)


def test_rational_exp_upper_encloses_simple_values() -> None:
    assert rational_exp_upper(Fraction()) == 1
    assert rational_exp_upper(Fraction(1, 10)) > Fraction(11, 10)


def test_s10_and_s11_local_bounds_are_nonzero_and_ordered() -> None:
    s10 = build_processed_local_bounds(
        published_effective_order_six_kernel("s10")
    )
    s11 = build_processed_local_bounds(
        published_effective_order_six_kernel("s11")
    )
    assert 0 < s10.degree3_residual_site_l1 < Fraction(1, 10**25)
    assert 0 < s11.degree3_residual_site_l1 < s10.degree3_residual_site_l1
    assert s11.r2_site_l1 < s10.r2_site_l1
    assert processor_distance_upper(s10, 144, 39) < Fraction(1, 10)
    assert processor_distance_upper(s11, 144, 35) < Fraction(1, 1000)
```

- [x] **Step 2: Run the tests and confirm missing-module failure**

Run: `python -m pytest -q tests/test_processed_local_bounds.py`

Expected: collection fails with `ModuleNotFoundError`.

- [x] **Step 3: Implement exact local fields**

Store cell and site l1 values for `R2`, `R4`, processed degree-three residual,
and processed degree-five residual.  Site density is exact cell density divided
by four.

- [x] **Step 4: Implement a rational exponential upper enclosure**

For nonnegative rational `x`, sum `sum(x^k/k!, k=0..m)` and bound the rest by

```text
next_term / (1 - x/(m+2))
```

when `x/(m+2) < 1`.  Reject values outside that convergence condition.  Bound
the processor rotation by

```text
rho <= exp(N * (R2_site/r^2 + R4_site/r^4)) - 1.
```

- [x] **Step 5: Run the local-bound tests**

Run: `python -m pytest -q tests/test_processed_local_bounds.py`

Expected: all tests pass; both candidates remain below their processor-distance thresholds.

### Task 3: Fail-Closed Finite-Step Budget Contract

**Files:**
- Create: `src/trottercert/processed_budget.py`
- Test: `tests/test_processed_budget.py`

**Interfaces:**
- Consumes: Task 2 local bounds.
- Produces: `ProcessedFiniteStepInputs`, `ProcessedFiniteStepBudget`, and `evaluate_processed_budget(inputs)`.

- [x] **Step 1: Write failing budget tests**

```python
from fractions import Fraction

import pytest

from trottercert.processed_budget import (
    ProcessedFiniteStepInputs,
    evaluate_processed_budget,
)


def test_budget_requires_degree7_and_remainder_proof_inputs() -> None:
    with pytest.raises(ValueError, match="degree-seven"):
        ProcessedFiniteStepInputs(
            kernel_name="s10", steps=39, n_sites=144,
            tolerance=Fraction(1, 10**6), degree7_site_l1=None,
            higher_order_remainder=None, local_log_theorem_id=None,
        )


def test_budget_sums_exact_log_contributions() -> None:
    inputs = ProcessedFiniteStepInputs(
        kernel_name="s10", steps=39, n_sites=144,
        tolerance=Fraction(1, 10**6),
        degree7_site_l1=Fraction(1, 10**6),
        higher_order_remainder=Fraction(1, 10**8),
        local_log_theorem_id="test-only",
    )
    result = evaluate_processed_budget(inputs)
    assert result.total == (
        result.degree3_residual + result.degree5_residual
        + result.degree7 + result.higher_order_remainder
    )
```

- [x] **Step 2: Run the tests and confirm missing-module failure**

Run: `python -m pytest -q tests/test_processed_budget.py`

Expected: collection fails with `ModuleNotFoundError`.

- [x] **Step 3: Implement exact repeated-log budget arithmetic**

For total time one, use

```text
degree3 residual = N * C3_site / r^2
degree5 residual = N * C5_site / r^4
degree7 contribution = N * C7_site / r^6.
```

Add the externally certified all-higher-degree remainder.  Return `accepted =
total <= tolerance`, remaining budget, processor distance, and the theorem ID.
Document that acceptance is conditional on the cited local-log theorem and
does not establish that theorem.

- [x] **Step 4: Run the budget tests**

Run: `python -m pytest -q tests/test_processed_budget.py`

Expected: all tests pass.

### Task 4: Emit the Pre-E7 Local Budget Report

**Files:**
- Create: `scripts/report_processed_local_budget.py`
- Create: `docs/experiments/processed-finite-step-budget/pre-e7-local-budget.json`
- Create: `docs/experiments/processed-finite-step-budget/README.md`
- Test: `tests/test_processed_local_budget_report.py`

**Interfaces:**
- Consumes: Tasks 1–3.
- Produces: canonical JSON for `s10@39`, `s10@47`, `s11@35`, and `s11@43`, with exact residual and processor bounds but explicit missing `degree7` and `higher_remainder` fields.

- [x] **Step 1: Write report-schema and missing-proof tests**

Test that every point contains exact rational pairs for residuals and processor distance, and that its status is exactly `awaiting_degree7_and_local_log_remainder`.  Mutating that status to `accepted` must make the verifier fail.

- [x] **Step 2: Implement report generation and verification**

The CLI must print one progress line per candidate/step point and write sorted compact JSON.  The verifier recomputes Task 2 fields and rejects any non-null degree-seven or remainder value until a theorem-bound artifact is explicitly linked by SHA-256.

- [x] **Step 3: Emit and verify the report**

Run:

```bash
PYTHONPATH=src python -u -m scripts.report_processed_local_budget \
  --output docs/experiments/processed-finite-step-budget/pre-e7-local-budget.json
python -m pytest -q tests/test_processed_local_budget_report.py
```

Expected: four points emitted; residual and processor values verified; no point labeled accepted.

- [x] **Step 4: Run focused and full regression suites**

Run the four new test files, then `python -m pytest -q`.  Record counts and
runtime in the README.

## Self-Review

- Exact residuals, processor rotation, future E7 input, and future remainder input are separate fields.
- The budget cannot silently replace missing physical or theorem evidence by zero.
- The global-log branch problem is stated as an external local-log theorem requirement.
- No HPC, E9, manuscript, or submission action is included.
