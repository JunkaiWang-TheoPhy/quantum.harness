# Dual Local-Log Tail via Even Generator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a source-regenerating exact-rational certificate for the E11-and-higher dual local-log tail of the fixed periodic `12 x 12` Issue-128 benchmark at `r=97`.

**Architecture:** A small pure module first certifies the scalar stage sums and analytic multiplier caps.  A fixed-instance layer then combines the existing D4--D7 right-generator constants, the centered principal-log stability estimate, and the even-generator `dexp` correction into one immutable rational record.  A separate CLI serializes and verifies a canonical JSON artifact while keeping the finite-step claim inconclusive until exact E9 is supplied by a later promotion phase.

**Tech Stack:** Python 3.11+, standard-library `dataclasses`, `fractions.Fraction`, `hashlib`, and `json`; existing `trottercert` cubic/interval/right-generator modules; pytest; Ruff.

## Global Constraints

- Work only in `tracks/qcs/solutions/WangTheoPhys/issue128` except for repository-level commands.
- Preserve every unrelated dirty or untracked file in the shared worktree.
- Use exact `Fraction` arithmetic for every claim-bearing value; decimal output is commentary only.
- Recompute the fourth-order constants from the frozen interval Suzuki stages; do not trust JSON-supplied constants.
- Fix `n_sites=144`, `cells=36`, and `steps=97`; reject any artifact with different values.
- The final mathematical object is `tau(W R_log)/(N/4)`, not the global right-generator operator norm.
- Keep `finite_step_status=inconclusive` and `missing=exact dual E9 pairing` throughout this plan.
- Do not edit the frozen exact E7 implementation or artifact.
- Do not infer an all-order ratio from E5, E7, or any sampled coefficient.
- Every commit message uses Lore format and ends with `Co-authored-by: OmX <omx@oh-my-codex.dev>`.

---

## File map

- Create `src/trottercert/dual_log_tail.py`: exact scalar sums, analytic gates, fixed-instance records, and the complete tail computation.
- Create `tests/test_dual_log_tail.py`: unit tests for scalar identities, rational gates, monotonicity, normalization, and the final bound.
- Create `scripts/certify_dual_log_tail.py`: canonical artifact builder, source binding, verifier, and CLI.
- Create `tests/test_dual_log_tail_certificate.py`: JSON round trips and fail-closed mutation attacks.
- Create `docs/experiments/processor-obstruction/dual-log-tail.json`: source-bound tail artifact with no E9 claim promotion.
- Create `docs/report/dual-log-tail-ledger.md`: human-readable derivation, exact-input ledger, outward decimals, and claim boundary.

The later E9 integration will receive its own plan because it consumes an
external HPC artifact not produced here and changes the theorem/audit claim.

---

### Task 1: Exact scalar stage sums and analytic caps

**Files:**
- Create: `src/trottercert/dual_log_tail.py`
- Create: `tests/test_dual_log_tail.py`

**Interfaces:**
- Consumes: `Sequence[IntervalStage]` from `trottercert.rigorous_fourth`.
- Produces: `average_stage_tail(stages, steps, first_omitted_degree) -> Fraction`.
- Produces: `pointwise_stage_tail(stages, steps, first_omitted_degree) -> Fraction`.
- Produces: `verify_analytic_caps(log_radius_cap: Fraction, dexp_cap: Fraction, multiplier_difference_cap: Fraction, even_series_cap: Fraction) -> None`.

- [x] **Step 1: Write failing scalar-tail tests**

Add imports and direct finite-sum reference functions to
`tests/test_dual_log_tail.py`:

```python
from fractions import Fraction

import pytest

from trottercert.dual_log_tail import (
    average_stage_tail,
    pointwise_stage_tail,
    verify_analytic_caps,
)
from trottercert.rigorous_fourth import (
    fourth_order_suzuki_interval_stages,
)


def test_geometric_stage_tail_matches_long_finite_sum() -> None:
    stages, _ = fourth_order_suzuki_interval_stages(
        4, decimal_digits=30
    )
    steps = 97
    first = 10
    prefix = Fraction()
    average_reference = Fraction()
    pointwise_reference = Fraction()
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        average_reference += coefficient * sum(
            ratio**degree for degree in range(first, 160)
        )
        pointwise_reference += coefficient * sum(
            (degree + 1) * ratio**degree
            for degree in range(first, 160)
        )
        prefix += coefficient
    assert average_stage_tail(stages, steps, first) > average_reference
    assert pointwise_stage_tail(stages, steps, first) > pointwise_reference
    assert average_stage_tail(stages, steps, first) - average_reference < Fraction(1, 10**100)
    assert pointwise_stage_tail(stages, steps, first) - pointwise_reference < Fraction(1, 10**98)


def test_stage_tail_rejects_invalid_convergence_region() -> None:
    stages, _ = fourth_order_suzuki_interval_stages(4, decimal_digits=30)
    with pytest.raises(ValueError, match="convergence"):
        average_stage_tail(stages, 1, 10)
    with pytest.raises(ValueError, match="omitted degree"):
        pointwise_stage_tail(stages, 97, -1)
```

- [x] **Step 2: Run the focused tests and confirm the import failure**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail.py
```

Expected: collection fails because `trottercert.dual_log_tail` does not exist.

- [x] **Step 3: Implement the two exact stage sums**

Create `src/trottercert/dual_log_tail.py` with:

```python
from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

from .rigorous_fourth import IntervalStage


def _validate_tail_inputs(
    steps: int,
    first_omitted_degree: int,
) -> None:
    if not isinstance(steps, int) or isinstance(steps, bool) or steps < 1:
        raise ValueError("steps must be a positive integer")
    if (
        not isinstance(first_omitted_degree, int)
        or isinstance(first_omitted_degree, bool)
        or first_omitted_degree < 0
    ):
        raise ValueError("first omitted degree must be nonnegative")


def average_stage_tail(
    stages: Sequence[IntervalStage],
    steps: int,
    first_omitted_degree: int,
) -> Fraction:
    _validate_tail_inputs(steps, first_omitted_degree)
    prefix = Fraction()
    total = Fraction()
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("stage tail is outside its convergence region")
        total += coefficient * ratio**first_omitted_degree / (1 - ratio)
        prefix += coefficient
    return total


def pointwise_stage_tail(
    stages: Sequence[IntervalStage],
    steps: int,
    first_omitted_degree: int,
) -> Fraction:
    _validate_tail_inputs(steps, first_omitted_degree)
    prefix = Fraction()
    total = Fraction()
    degree = first_omitted_degree
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("stage tail is outside its convergence region")
        total += (
            coefficient
            * ratio**degree
            * ((degree + 1) - degree * ratio)
            / (1 - ratio) ** 2
        )
        prefix += coefficient
    return total
```

- [x] **Step 4: Add failing analytic-cap tests**

Append:

```python
def test_selected_analytic_caps_are_proved_by_rational_gates() -> None:
    verify_analytic_caps(
        log_radius_cap=Fraction(3, 2),
        dexp_cap=Fraction(8, 5),
        multiplier_difference_cap=Fraction(1),
        even_series_cap=Fraction(10, 33),
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("log_radius_cap", Fraction(151, 100), "log radius"),
        ("dexp_cap", Fraction(159, 100), "dexp cap"),
        ("multiplier_difference_cap", Fraction(99, 100), "multiplier"),
        ("even_series_cap", Fraction(3, 10), "even-series"),
    ],
)
def test_analytic_caps_reject_unsafe_claims(
    field: str,
    value: Fraction,
    message: str,
) -> None:
    arguments = {
        "log_radius_cap": Fraction(3, 2),
        "dexp_cap": Fraction(8, 5),
        "multiplier_difference_cap": Fraction(1),
        "even_series_cap": Fraction(10, 33),
    }
    arguments[field] = value
    with pytest.raises(ValueError, match=message):
        verify_analytic_caps(**arguments)
```

- [x] **Step 5: Implement the exact analytic gates**

Add:

```python
LOG_RADIUS_CAP = Fraction(3, 2)
DEXP_CAP = Fraction(8, 5)
MULTIPLIER_DIFFERENCE_CAP = Fraction(1)
EVEN_SERIES_CAP = Fraction(10, 33)


def verify_analytic_caps(
    *,
    log_radius_cap: Fraction,
    dexp_cap: Fraction,
    multiplier_difference_cap: Fraction,
    even_series_cap: Fraction,
) -> None:
    radius = Fraction(log_radius_cap)
    if radius > LOG_RADIUS_CAP:
        raise ValueError("log radius exceeds the proved interval")
    # sin(x)/x >= 1-x^2/6 on [0, 3/2].
    required_dexp = 1 / (1 - radius**2 / 6)
    if Fraction(dexp_cap) < required_dexp:
        raise ValueError("dexp cap is below the rational sine bound")
    # For y=2x and |x|<=3/2, the cotangent remainder satisfies
    # |x*cot(x)-1|/|x| <= x/(3*(1-x^2/6)) <= 4/5.
    # Hence |(f(iy)-1)/(iy)|^2 <= 1/4+(2/5)^2 < 1.
    if Fraction(multiplier_difference_cap) < 1:
        raise ValueError("multiplier-difference cap is below the proved cap")
    # The k=1 term is 1/3! and every later ratio is <=9/20.
    required_even_series = Fraction(1, 6) / (1 - Fraction(9, 20))
    if Fraction(even_series_cap) < required_even_series:
        raise ValueError("even-series cap is below the ratio majorant")
```

- [x] **Step 6: Run Task-1 tests and static checks**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail.py
ruff check --select I,UP035 src/trottercert/dual_log_tail.py tests/test_dual_log_tail.py
python -m compileall -q src/trottercert/dual_log_tail.py tests/test_dual_log_tail.py
git diff --check
```

Expected: all commands exit zero.

- [x] **Step 7: Commit Task 1**

Stage only the two Task-1 files and commit with subject:

```text
feat(issue128): add rational dual-tail scalar kernel
```

The Lore body records the two closed-form sums, the rational analytic lemmas,
the focused test command, constraints, and confidence, followed by the OmX
co-author trailer.

---

### Task 2: Fixed-instance right-generator and logarithm envelope

**Files:**
- Modify: `src/trottercert/dual_log_tail.py`
- Modify: `tests/test_dual_log_tail.py`

**Interfaces:**
- Consumes: `RefinedFourthOrderConstants` from `trottercert.refined_error`.
- Produces: frozen `DualTailGeometry` and `DualTailEnvelope` dataclasses.
- Produces: `issue128_dual_tail_geometry() -> DualTailGeometry`.
- Produces: `build_dual_tail_envelope(constants: RefinedFourthOrderConstants, geometry: DualTailGeometry) -> DualTailEnvelope`.

- [x] **Step 1: Write failing fixed-geometry and moment tests**

Append:

```python
from trottercert.dual_log_tail import (
    build_dual_tail_envelope,
    issue128_dual_tail_geometry,
)
from trottercert.refined_error import build_refined_fourth_order_constants


def test_issue128_geometry_has_exact_moment_caps() -> None:
    geometry = issue128_dual_tail_geometry()
    assert geometry.n_sites == 144
    assert geometry.cells == 36
    assert geometry.steps == 97
    assert geometry.h_hs_squared == 54
    assert geometry.w_hs_squared == Fraction(23085, 4)
    assert geometry.h_hs_cap**2 >= geometry.h_hs_squared
    assert (
        geometry.w_hs_per_cell_cap**2
        >= geometry.w_hs_squared / geometry.cells**2
    )
    assert geometry.max_w_pauli_coefficient == Fraction(1, 4)


def test_fixed_envelope_regression_and_branch_gate() -> None:
    constants = build_refined_fourth_order_constants(
        decimal_digits=30,
        quantization_digits=24,
    )
    result = build_dual_tail_envelope(
        constants, issue128_dual_tail_geometry()
    )
    assert result.centered_exact_phase_radius < Fraction(3, 2)
    assert result.centered_log_radius_bound < Fraction(3, 2)
    assert Fraction(2146, 10**9) < result.average_generator_defect
    assert result.average_generator_defect < Fraction(2148, 10**9)
    assert Fraction(1187, 10**8) < result.pointwise_generator_defect
    assert result.pointwise_generator_defect < Fraction(1188, 10**8)
    assert Fraction(354, 10**10) < result.log_defect
    assert result.log_defect < Fraction(355, 10**10)
```

- [x] **Step 2: Run the new tests and confirm missing interfaces**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail.py -k 'geometry or envelope'
```

Expected: import or attribute failure for the new records/functions.

- [x] **Step 3: Add immutable geometry and envelope records**

Add:

```python
from dataclasses import dataclass

from .refined_error import RefinedFourthOrderConstants


@dataclass(frozen=True, slots=True)
class DualTailGeometry:
    n_sites: int
    cells: int
    steps: int
    centered_h_op_cap: Fraction
    h_hs_squared: Fraction
    h_hs_cap: Fraction
    w_hs_squared: Fraction
    w_hs_per_cell_cap: Fraction
    max_w_pauli_coefficient: Fraction
    log_radius_cap: Fraction
    dexp_cap: Fraction
    multiplier_difference_cap: Fraction
    even_series_cap: Fraction


@dataclass(frozen=True, slots=True)
class DualTailEnvelope:
    average_generator_defect: Fraction
    one_step_unitary_defect: Fraction
    pointwise_generator_defect: Fraction
    relative_log_defect: Fraction
    centered_exact_phase_radius: Fraction
    centered_log_radius_bound: Fraction
    log_defect: Fraction
    log_derivative_defect_hs: Fraction
    log_commutator_defect_hs: Fraction


def issue128_dual_tail_geometry() -> DualTailGeometry:
    return DualTailGeometry(
        n_sites=144,
        cells=36,
        steps=97,
        centered_h_op_cap=Fraction(144),
        h_hs_squared=Fraction(54),
        h_hs_cap=Fraction(15, 2),
        w_hs_squared=Fraction(23085, 4),
        w_hs_per_cell_cap=Fraction(17, 8),
        max_w_pauli_coefficient=Fraction(1, 4),
        log_radius_cap=LOG_RADIUS_CAP,
        dexp_cap=DEXP_CAP,
        multiplier_difference_cap=MULTIPLIER_DIFFERENCE_CAP,
        even_series_cap=EVEN_SERIES_CAP,
    )
```

- [x] **Step 4: Implement geometry validation and the full envelope**

Add private `_verify_geometry` checks for every exact fixed value and cap.
Then implement:

```python
def build_dual_tail_envelope(
    constants: RefinedFourthOrderConstants,
    geometry: DualTailGeometry,
) -> DualTailEnvelope:
    _verify_geometry(geometry)
    verify_analytic_caps(
        log_radius_cap=geometry.log_radius_cap,
        dexp_cap=geometry.dexp_cap,
        multiplier_difference_cap=geometry.multiplier_difference_cap,
        even_series_cap=geometry.even_series_cap,
    )
    t = Fraction(1, geometry.steps)
    n_sites = Fraction(geometry.n_sites)
    direct_degrees = (
        constants.d4_site * t**4 / 5
        + constants.d5_site * t**5 / 6
        + constants.d6_site * t**6 / 7
        + constants.d7_site * t**7 / 8
    )
    average_defect = n_sites * (
        direct_degrees
        + Fraction(3, 8)
        * average_stage_tail(constants.stages, geometry.steps, 8)
    )
    pointwise_defect = n_sites * (
        constants.d4_site * t**4
        + constants.d5_site * t**5
        + constants.d6_site * t**6
        + constants.d7_site * t**7
        + Fraction(3, 8)
        * pointwise_stage_tail(constants.stages, geometry.steps, 8)
    )
    one_step = t * average_defect
    if one_step >= 1:
        raise ArithmeticError("unitary defect is outside the logarithm gate")
    relative_log = one_step / (1 - one_step)
    exact_radius = t * geometry.centered_h_op_cap
    centered_radius = exact_radius + relative_log
    if centered_radius > geometry.log_radius_cap:
        raise ArithmeticError("centered principal-log radius gate failed")
    log_defect = geometry.dexp_cap * relative_log
    log_derivative_defect = (
        geometry.dexp_cap * pointwise_defect
        + 2
        * geometry.multiplier_difference_cap
        * log_defect
        * geometry.h_hs_cap
    )
    commutator_defect = 2 * (
        exact_radius * log_derivative_defect
        + log_defect * geometry.h_hs_cap
        + log_defect * log_derivative_defect
    )
    return DualTailEnvelope(
        average_generator_defect=average_defect,
        one_step_unitary_defect=one_step,
        pointwise_generator_defect=pointwise_defect,
        relative_log_defect=relative_log,
        centered_exact_phase_radius=exact_radius,
        centered_log_radius_bound=centered_radius,
        log_defect=log_defect,
        log_derivative_defect_hs=log_derivative_defect,
        log_commutator_defect_hs=commutator_defect,
    )
```

- [x] **Step 5: Add mutation-style geometry tests**

Use `dataclasses.replace` to lower each cap or change each fixed integer and
assert `_verify_geometry` is reached through `build_dual_tail_envelope`:

```python
from dataclasses import replace


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("n_sites", 143, "site count"),
        ("cells", 35, "cell count"),
        ("steps", 96, "step count"),
        ("h_hs_squared", Fraction(53), "H moment"),
        ("w_hs_squared", Fraction(23084, 4), "W moment"),
        ("max_w_pauli_coefficient", Fraction(1, 5), "W Pauli"),
    ],
)
def test_geometry_mutations_fail_closed(
    field: str,
    value: object,
    message: str,
) -> None:
    constants = build_refined_fourth_order_constants(
        decimal_digits=30, quantization_digits=24
    )
    forged = replace(issue128_dual_tail_geometry(), **{field: value})
    with pytest.raises(ValueError, match=message):
        build_dual_tail_envelope(constants, forged)
```

- [x] **Step 6: Run Task-2 tests and static checks**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail.py
ruff check --select I,UP035 src/trottercert/dual_log_tail.py tests/test_dual_log_tail.py
python -m compileall -q src/trottercert/dual_log_tail.py tests/test_dual_log_tail.py
git diff --check
```

Expected: all commands exit zero.

- [x] **Step 7: Commit Task 2**

Commit only the module and unit test with subject:

```text
feat(issue128): certify centered dual-log envelope
```

The Lore body must name the fixed instance, exact moments, branch gate,
pointwise/average distinction, test commands, and claim constraint.

---

### Task 3: Complete compatible E11-and-higher bound

**Files:**
- Modify: `src/trottercert/dual_log_tail.py`
- Modify: `tests/test_dual_log_tail.py`

**Interfaces:**
- Consumes: `DualTailEnvelope` and `DualTailGeometry` from Task 2.
- Produces: frozen `DualLogTailBound`.
- Produces: `certify_issue128_dual_log_tail(constants: RefinedFourthOrderConstants | None = None) -> DualLogTailBound`.

- [x] **Step 1: Write the failing final-bound regression**

Append:

```python
from trottercert.dual_log_tail import certify_issue128_dual_log_tail


def test_complete_e11_tail_is_compatible_and_below_preaudit_gate() -> None:
    result = certify_issue128_dual_log_tail()
    assert Fraction(6448, 10**15) < result.direct_even_generator_tail
    assert result.direct_even_generator_tail < Fraction(6450, 10**15)
    assert Fraction(242, 10**15) < result.dexp_correction_tail
    assert result.dexp_correction_tail < Fraction(244, 10**15)
    assert result.total_log_tail == (
        result.direct_even_generator_tail + result.dexp_correction_tail
    )
    assert Fraction(669, 10**14) < result.total_log_tail
    assert result.total_log_tail < Fraction(670, 10**14)
    assert result.first_omitted_log_degree == 11
    assert result.claim_scope == "fixed_12x12_r97_dual_local_log"
```

- [x] **Step 2: Run the regression and confirm the missing interface**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail.py::test_complete_e11_tail_is_compatible_and_below_preaudit_gate
```

Expected: import or attribute failure.

- [x] **Step 3: Implement the final immutable record and composition**

Add:

```python
@dataclass(frozen=True, slots=True)
class DualLogTailBound:
    geometry: DualTailGeometry
    envelope: DualTailEnvelope
    first_omitted_generator_degree: int
    first_omitted_log_degree: int
    direct_even_generator_tail: Fraction
    dexp_correction_tail: Fraction
    total_log_tail: Fraction
    claim_scope: str


def certify_issue128_dual_log_tail(
    constants: RefinedFourthOrderConstants | None = None,
) -> DualLogTailBound:
    if constants is None:
        from .refined_error import build_refined_fourth_order_constants

        constants = build_refined_fourth_order_constants(
            decimal_digits=30,
            quantization_digits=24,
        )
    geometry = issue128_dual_tail_geometry()
    envelope = build_dual_tail_envelope(constants, geometry)
    direct = (
        Fraction(3, 8)
        * average_stage_tail(
            constants.stages,
            geometry.steps,
            first_omitted_degree=10,
        )
    )
    correction = (
        geometry.w_hs_per_cell_cap
        * 2
        * envelope.log_defect
        * envelope.log_commutator_defect_hs
        * geometry.even_series_cap
        / 11
    )
    total = direct + correction
    if total <= 0:
        raise ArithmeticError("dual local-log tail must be positive")
    return DualLogTailBound(
        geometry=geometry,
        envelope=envelope,
        first_omitted_generator_degree=10,
        first_omitted_log_degree=11,
        direct_even_generator_tail=direct,
        dexp_correction_tail=correction,
        total_log_tail=total,
        claim_scope="fixed_12x12_r97_dual_local_log",
    )
```

- [x] **Step 4: Add normalization and monotonicity tests**

Append:

```python
def test_direct_tail_has_the_exact_dual_normalization() -> None:
    result = certify_issue128_dual_log_tail()
    constants = build_refined_fourth_order_constants(
        decimal_digits=30, quantization_digits=24
    )
    raw_cell = (
        Fraction(3, 2)
        * average_stage_tail(constants.stages, 97, 10)
    )
    assert result.direct_even_generator_tail == raw_cell * Fraction(1, 4)


def test_total_tail_decreases_at_larger_step_counts() -> None:
    constants = build_refined_fourth_order_constants(
        decimal_digits=30, quantization_digits=24
    )
    base = issue128_dual_tail_geometry()
    at_97 = build_dual_tail_envelope(constants, base)
    at_98 = build_dual_tail_envelope(constants, replace(base, steps=98))
    assert at_98.average_generator_defect < at_97.average_generator_defect
    assert at_98.log_defect < at_97.log_defect
```

For this test only, relax `_verify_geometry` into two functions: one validates
the fixed certificate geometry and one validates the mathematical geometry.
`build_dual_tail_envelope` uses the mathematical validator (`steps >= 97`),
while `certify_issue128_dual_log_tail` invokes the fixed validator and still
requires exactly `97`.  The JSON verifier in Task 4 also uses the fixed
validator.

- [x] **Step 5: Run Task-3 tests and the existing refined-tail regressions**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail.py tests/test_refined_error.py
ruff check --select I,UP035 src/trottercert/dual_log_tail.py tests/test_dual_log_tail.py
python -m compileall -q src/trottercert/dual_log_tail.py tests/test_dual_log_tail.py
git diff --check
```

Expected: all commands exit zero.

- [x] **Step 6: Commit Task 3**

Commit only the two Task-3 files with subject:

```text
feat(issue128): close compatible E11+ dual tail
```

The Lore body records the direct and conversion bounds separately, the exact
normalization identity, fixed-scope constraint, tests, and confidence.

---

### Task 4: Canonical artifact, verifier, and evidence ledger

**Files:**
- Create: `scripts/certify_dual_log_tail.py`
- Create: `tests/test_dual_log_tail_certificate.py`
- Create: `docs/experiments/processor-obstruction/dual-log-tail.json`
- Create: `docs/report/dual-log-tail-ledger.md`

**Interfaces:**
- Consumes: `certify_issue128_dual_log_tail()` from Task 3.
- Consumes: frozen E5 artifact `docs/experiments/processor-obstruction/extensive-commutant-witness.json` only for provenance and the current pre-E9 comparison.
- Consumes: frozen E7 artifact `docs/experiments/processor-obstruction/dual-e7-pairing.json` only for provenance and the current pre-E9 comparison.
- Produces: `build_payload() -> dict[str, object]` and `verify_payload(payload: Mapping[str, object]) -> None`.
- Produces: canonical schema-v1 `issue128_dual_log_tail` JSON.

- [x] **Step 1: Write failing payload and mutation tests**

Create `tests/test_dual_log_tail_certificate.py`:

```python
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.certify_dual_log_tail import build_payload, verify_payload

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "docs/experiments/processor-obstruction/dual-log-tail.json"
)


def test_payload_round_trip_is_source_regenerated() -> None:
    payload = build_payload()
    verify_payload(payload)
    assert payload["kind"] == "issue128_dual_log_tail"
    assert payload["claim"] == {
        "dual_e11_plus_tail": "certified",
        "finite_step_status": "inconclusive",
        "missing": "exact dual E9 pairing",
    }
    assert payload["bounds"]["total_log_tail"] == [
        build_payload()["bounds"]["total_log_tail"][0],
        build_payload()["bounds"]["total_log_tail"][1],
    ]


def test_checked_in_artifact_is_canonical_and_current() -> None:
    encoded = ARTIFACT.read_bytes()
    payload = json.loads(encoded)
    verify_payload(payload)
    assert encoded == (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    assert payload == build_payload()


@pytest.mark.parametrize(
    ("path", "message"),
    [
        (("instance", "steps"), "payload regeneration mismatch"),
        (("moments", "w_hs_squared"), "payload regeneration mismatch"),
        (("analytic_caps", "dexp"), "payload regeneration mismatch"),
        (("envelope", "log_defect"), "payload regeneration mismatch"),
        (("bounds", "direct_even_generator_tail"), "payload regeneration mismatch"),
        (("bounds", "dexp_correction_tail"), "payload regeneration mismatch"),
        (("claim", "finite_step_status"), "payload regeneration mismatch"),
    ],
)
def test_payload_mutations_fail_closed(
    path: tuple[str, str],
    message: str,
) -> None:
    forged = copy.deepcopy(build_payload())
    parent, child = path
    value = forged[parent][child]
    forged[parent][child] = "forged" if isinstance(value, str) else [0, 1]
    with pytest.raises(ValueError, match=message):
        verify_payload(forged)
```

- [x] **Step 2: Run tests and confirm the missing script**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail_certificate.py
```

Expected: collection fails because `scripts.certify_dual_log_tail` does not
exist.

- [x] **Step 3: Implement canonical serialization and source binding**

Create `scripts/certify_dual_log_tail.py` with:

```python
from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

from trottercert.dual_log_tail import certify_issue128_dual_log_tail

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/dual-log-tail.json"
)
E5_PATH = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
E7_PATH = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
)
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/dual_log_tail.py",
    ISSUE_ROOT / "src/trottercert/refined_error.py",
    ISSUE_ROOT / "src/trottercert/rigorous_fourth.py",
    Path(__file__).resolve(),
)


def _pair(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
```

Serialize every field of `DualTailGeometry`, `DualTailEnvelope`, and
`DualLogTailBound` explicitly as canonical rational pairs.  Do not serialize
dataclasses with `asdict`, because field additions must require a conscious
schema update.  Include E5/E7 file SHA-256 values and source paths relative to
`ISSUE_ROOT`.

Set the exact claim block to:

```python
{
    "dual_e11_plus_tail": "certified",
    "finite_step_status": "inconclusive",
    "missing": "exact dual E9 pairing",
}
```

Implement `verify_payload` as strict equality with a regenerated
`build_payload()` after first checking that the input is a mapping.  The CLI
supports mutually exclusive `--output PATH` and `--verify PATH`.  `--output`
writes to a sibling `.pending.<pid>` path, rereads and verifies it, then uses
`Path.replace` to publish.  It refuses to overwrite an existing final output
unless `--replace` is explicitly present; `--replace` is allowed only for the
checked-in default artifact.

- [x] **Step 4: Run payload tests and generate the artifact**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail_certificate.py -k 'not checked_in'
PYTHONPATH=src:. python scripts/certify_dual_log_tail.py --output docs/experiments/processor-obstruction/dual-log-tail.json
PYTHONPATH=src:. python scripts/certify_dual_log_tail.py --verify docs/experiments/processor-obstruction/dual-log-tail.json
PYTHONPATH=src:. pytest -q tests/test_dual_log_tail_certificate.py
```

Expected: the first focused subset passes, the generator and verifier print
success, and the complete certificate test passes.

- [x] **Step 5: Write the evidence ledger from regenerated values**

Create `docs/report/dual-log-tail-ledger.md` with these sections and exact
content sources:

1. `Scope`: fixed periodic `12 x 12`, `r=97`, principal local log, E11+ only.
2. `Identity`: the displayed `K_even=sinh(ad_A)/ad_A A'` relation and the
   two-defect commutator inequality from the design.
3. `Exact inputs`: N, cells, H/W moments, rational caps, formula stage/source
   hashes, and E5/E7 artifact hashes copied from the generated JSON.
4. `Rational ledger`: every pair in `envelope` and `bounds`, with outward
   decimals generated by `float(Fraction(*pair))` and labeled non-normative.
5. `Claim boundary`: E11+ certified; exact E9 missing; no finite-step
   promotion; no all-size theorem.
6. `Verification`: the exact CLI commands above and the focused/full test
   commands below.

The report must state explicitly that the old D8+ operator tail was not
substituted: the direct scalar expression was rederived with the `1/4` Pauli
coefficient and per-cell normalization, and the `dexp` difference was bounded
separately.

- [x] **Step 6: Run focused and full verification**

Run:

```bash
ruff check --select I,UP035 \
  src/trottercert/dual_log_tail.py \
  scripts/certify_dual_log_tail.py \
  tests/test_dual_log_tail.py \
  tests/test_dual_log_tail_certificate.py
python -m compileall -q \
  src/trottercert/dual_log_tail.py \
  scripts/certify_dual_log_tail.py \
  tests/test_dual_log_tail.py \
  tests/test_dual_log_tail_certificate.py
PYTHONPATH=src:. pytest -q \
  tests/test_dual_log_tail.py \
  tests/test_dual_log_tail_certificate.py \
  tests/test_refined_error.py \
  tests/test_dual_e9_pairing.py
PYTHONPATH=src:. pytest -q -m 'not slow'
git diff --check
```

Expected: all commands exit zero.  Record exact test counts and elapsed time
in the commit body and evidence ledger.

- [x] **Step 7: Commit Task 4**

Stage only the four Task-4 files and any Task-3 files changed solely to repair
a verifier-discovered bug.  Commit with subject:

```text
cert(issue128): emit compatible E11+ dual tail
```

The Lore body records artifact/source hashes, exact test counts, claim
boundary, no-HPC constraint for this local phase, and confidence.

---

## Self-review

- Spec coverage: Tasks 1--3 cover both scalar tails, analytic caps, fixed
  moments, branch stability, `R/R'`, commutator conversion, normalization, and
  the total bound.  Task 4 covers source regeneration, provenance, canonical
  serialization, mutation rejection, evidence, and the fail-closed claim.
- Deferred scope: exact E9 ingestion and theorem promotion are intentionally
  excluded because the required HPC artifact does not yet exist.  The current
  plan still produces a complete, independently useful E11+ certificate.
- Placeholder scan: every task names concrete interfaces, test code, failure
  modes, implementation equations, commands, and commit boundaries; no
  unfinished marker remains.
- Type consistency: Task 2 defines `DualTailGeometry` and `DualTailEnvelope`;
  Task 3 consumes them and defines `DualLogTailBound`; Task 4 consumes only the
  Task-3 public function and serializes those named fields.
- Claim safety: no task can set `finite_step_status` to anything other than
  `inconclusive`; the future E9 promotion requires a separate reviewed change.

## Executed-plan provenance refinement

During final artifact review, whole-file hashes exposed unrelated uncommitted
cache and center-scan additions in two legacy generator modules.  The executed
implementation was tightened after the four planned tasks: it now consumes
the committed, verified schema-v3 D5-integrated main certificate, recovers its
D4--D7 density caps exactly from the recorded contributions, and independently
rebuilds the 31 interval stages inside `dual_log_tail.py`.  The final bound is
`6.495092545992439e-12`, tighter than the plan's generic-source pre-audit, and
the canonical verifier no longer depends on either dirty legacy file.
