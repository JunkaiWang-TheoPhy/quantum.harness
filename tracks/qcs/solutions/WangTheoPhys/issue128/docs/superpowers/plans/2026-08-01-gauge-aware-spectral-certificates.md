# Gauge-Aware Spectral Certificates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an exact primal/dual solver for the calibrated spectral gauge, demonstrate it on complete small-matrix commutator images, and emit a fail-closed Issue-128 claim audit.

**Architecture:** A metric-aware exact linear-algebra core computes the quotient by a supplied gauge span and labels full versus restricted claims explicitly. A small-matrix adapter constructs the complete Hermitian commutator image automatically, while a deterministic CLI binds the current processor-obstruction artifact to a calibrated-time claim audit and three exact reference examples.

**Tech Stack:** Python 3.11+, `fractions.Fraction`, SymPy 1.13+, dataclasses, enums, canonical JSON, SHA-256, pytest 8+

## Global Constraints

- Work only under `tracks/qcs/solutions/WangTheoPhys/issue128` on the existing `codex/issue128-hpc-six-hour` branch.
- Preserve all unrelated dirty-worktree changes and stage only files named by the active task.
- Accept exact integers, `Fraction`, and exact SymPy expressions; reject Python floats and `sympy.Float` atoms.
- Keep fixed-time, calibrated-time, full-commutator-image, and restricted-processor claims distinct.
- Do not authorize or submit HPC work in this plan.
- Every JSON status must be recomputed by a verifier; a mutated status or source digest must fail.

---

### Task 1: Exact Metric-Aware Gauge Quotient

**Files:**
- Create: `src/trottercert/spectral_gauge.py`
- Create: `tests/test_spectral_gauge.py`

**Interfaces:**
- Consumes: exact coordinate vectors and optional diagonal Hilbert--Schmidt metric weights.
- Produces: `GaugeScope`, `GaugeDecomposition`, `decompose_gauge`, `verify_gauge_decomposition`, and `primitive_integer_vector`.

- [ ] **Step 1: Write failing tests for primal, dual, metric, and scope behavior**

```python
from fractions import Fraction

import pytest

from trottercert.spectral_gauge import (
    GaugeScope,
    decompose_gauge,
    verify_gauge_decomposition,
)


def test_exact_primal_reconstruction() -> None:
    result = decompose_gauge(
        (Fraction(2), Fraction(3)),
        ((1, 0), (0, 1)),
        scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
        completeness_id="test-full-span",
    )
    assert result.status == "removable"
    assert result.residual == (0, 0)
    verify_gauge_decomposition(result)


def test_full_and_restricted_dual_claims_differ() -> None:
    full = decompose_gauge(
        (1, 2, 3),
        ((1, 0, 0),),
        scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
        completeness_id="test-full-span",
    )
    restricted = decompose_gauge(
        (1, 2, 3),
        ((1, 0, 0),),
        scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
    )
    assert full.status == "spectral_obstruction"
    assert restricted.status == "restricted_obstruction"
    assert full.primitive_witness == restricted.primitive_witness == (0, 2, 3)
    verify_gauge_decomposition(full)
    verify_gauge_decomposition(restricted)


def test_weighted_metric_changes_projection() -> None:
    result = decompose_gauge(
        (1, 0),
        ((1, 1),),
        metric=(1, 2),
        scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
    )
    assert result.projection == (Fraction(1, 3), Fraction(1, 3))
    assert result.residual == (Fraction(2, 3), Fraction(-1, 3))
    assert result.primitive_witness == (2, -1)
    verify_gauge_decomposition(result)


def test_exact_and_scope_inputs_fail_closed() -> None:
    with pytest.raises(ValueError, match="exact"):
        decompose_gauge(
            (0.1, 1),
            (),
            scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
        )
    with pytest.raises(ValueError, match="completeness"):
        decompose_gauge(
            (1, 0),
            (),
            scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
        )
```

- [ ] **Step 2: Run the tests and confirm the missing-module failure**

Run:

```bash
python -m pytest -q tests/test_spectral_gauge.py
```

Expected: collection fails with `ModuleNotFoundError: trottercert.spectral_gauge`.

- [ ] **Step 3: Implement exact vector validation and primitive normalization**

Implement these public definitions:

```python
class GaugeScope(str, Enum):
    FULL_COMMUTATOR_IMAGE = "full_commutator_image"
    RESTRICTED_PROCESSOR_SPAN = "restricted_processor_span"


@dataclass(frozen=True)
class GaugeDecomposition:
    target: tuple[sp.Expr, ...]
    independent_generators: tuple[tuple[sp.Expr, ...], ...]
    pivot_indices: tuple[int, ...]
    metric: tuple[sp.Expr, ...]
    coefficients: tuple[sp.Expr, ...]
    projection: tuple[sp.Expr, ...]
    residual: tuple[sp.Expr, ...]
    primitive_witness: tuple[int, ...] | None
    witness_target_pairing: sp.Expr
    scope: GaugeScope
    completeness_id: str | None
    status: str


def primitive_integer_vector(values: Sequence[sp.Expr]) -> tuple[int, ...] | None:
    """Clear rational denominators, divide by the integer gcd, and fix sign."""


def decompose_gauge(
    target: Sequence[ExactScalar],
    generators: Sequence[Sequence[ExactScalar]],
    *,
    metric: Sequence[ExactScalar] | None = None,
    scope: GaugeScope,
    completeness_id: str | None = None,
) -> GaugeDecomposition:
    """Project `target` onto an exact gauge span and emit a primal or dual result."""
```

Validation must reject booleans, Python floats, non-real expressions,
`sympy.Float` atoms, nonpositive metric weights, empty targets, and dimension
mismatches.  `FULL_COMMUTATOR_IMAGE` must require a nonempty
`completeness_id`.

- [ ] **Step 4: Implement exact weighted projection and verification**

Use an exact generator matrix `G`, obtain independent columns from RREF pivot
indices, and compute

```python
gram = basis.T * sp.diag(*metric) * basis
rhs = basis.T * sp.diag(*metric) * target_column
coefficients = gram.inv() * rhs
projection = basis * coefficients
residual = target_column - projection
```

For an empty gauge span, the projection is zero.  Set statuses exactly as:

```python
if residual.is_zero_matrix:
    status = "removable"
elif scope is GaugeScope.FULL_COMMUTATOR_IMAGE:
    status = "spectral_obstruction"
else:
    status = "restricted_obstruction"
```

Implement:

```python
def verify_gauge_decomposition(result: GaugeDecomposition) -> None:
    """Recompute projection, residual, weighted orthogonality, and claim scope."""
```

The verifier must require `witness_target_pairing != 0` for either obstruction
status and exact weighted orthogonality to every independent generator.

- [ ] **Step 5: Run focused tests**

Run:

```bash
python -m pytest -q tests/test_spectral_gauge.py
```

Expected: all tests pass.

- [ ] **Step 6: Commit Task 1**

```bash
git add \
  tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/spectral_gauge.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_spectral_gauge.py
git commit -m "feat(issue128): add exact spectral gauge quotient"
```

---

### Task 2: Complete Small-Matrix Commutator Image

**Files:**
- Modify: `src/trottercert/spectral_gauge.py`
- Modify: `tests/test_spectral_gauge.py`

**Interfaces:**
- Consumes: exact square Hermitian SymPy matrices `H` and `E`.
- Produces: `MatrixGaugeDecomposition`, `hermitian_basis`, `hermitian_coordinates`, `matrix_from_hermitian_coordinates`, and `decompose_matrix_spectral_gauge`.

- [ ] **Step 1: Add failing matrix reference tests**

```python
import sympy as sp

from trottercert.spectral_gauge import decompose_matrix_spectral_gauge


def test_off_diagonal_defect_is_commutator_removable() -> None:
    hamiltonian = sp.diag(-1, 1)
    defect = sp.Matrix([[0, 1], [1, 0]])
    result = decompose_matrix_spectral_gauge(hamiltonian, defect)
    assert result.gauge.status == "removable"
    assert result.witness_matrix is None


def test_time_calibration_absorbs_hamiltonian_parallel_defect() -> None:
    hamiltonian = sp.diag(-1, 1)
    fixed_time = decompose_matrix_spectral_gauge(
        hamiltonian,
        hamiltonian,
        include_time_calibration=False,
    )
    calibrated = decompose_matrix_spectral_gauge(
        hamiltonian,
        hamiltonian,
        include_time_calibration=True,
    )
    assert fixed_time.gauge.status == "spectral_obstruction"
    assert calibrated.gauge.status == "removable"


def test_commutant_witness_survives_phase_and_time_gauge() -> None:
    hamiltonian = sp.diag(-1, 0, 2)
    defect = sp.diag(2, -3, 1)
    result = decompose_matrix_spectral_gauge(hamiltonian, defect)
    assert result.gauge.status == "spectral_obstruction"
    assert result.witness_matrix is not None
    assert sp.simplify(result.witness_matrix * hamiltonian - hamiltonian * result.witness_matrix) == sp.zeros(3)
    assert sp.trace(result.witness_matrix) == 0
    assert sp.trace(result.witness_matrix * hamiltonian) == 0


def test_matrix_inputs_reject_floats_and_nonhermitian_values() -> None:
    with pytest.raises(ValueError, match="exact"):
        decompose_matrix_spectral_gauge(sp.diag(0.1, 1), sp.eye(2))
    with pytest.raises(ValueError, match="Hermitian"):
        decompose_matrix_spectral_gauge(
            sp.Matrix([[0, 1], [0, 0]]),
            sp.eye(2),
        )
```

- [ ] **Step 2: Run the matrix tests and confirm missing-interface failures**

Run:

```bash
python -m pytest -q tests/test_spectral_gauge.py -k 'matrix or diagonal or calibration or commutant'
```

Expected: import or attribute failure for `decompose_matrix_spectral_gauge`.

- [ ] **Step 3: Implement the exact Hermitian basis and coordinate map**

Use this ordered real basis for dimension `d`:

1. `|j><j|` for every diagonal index;
2. `|j><k| + |k><j|` for `j < k`;
3. `i|j><k| - i|k><j|` for `j < k`.

Its Hilbert--Schmidt metric is diagonal with weight one on diagonal elements
and two on both off-diagonal families.  Implement exact round-trip checks in:

Implement functions with these exact signatures:

```python
def hermitian_basis(dimension: int) -> tuple[sp.Matrix, ...]:
    """Return the ordered real Hermitian basis described above."""


def hermitian_coordinates(matrix: sp.Matrix) -> tuple[sp.Expr, ...]:
    """Return diagonal, real-off-diagonal, then imaginary coordinates."""


def matrix_from_hermitian_coordinates(
    coordinates: Sequence[ExactScalar],
    dimension: int,
) -> sp.Matrix:
    """Reconstruct and validate an exact Hermitian matrix."""
```

The implementation body of `hermitian_basis` must append `sp.eye(d)[j,j]`
matrix units explicitly, followed by the two off-diagonal families.  The
coordinate function reads `matrix[j,j]`, `re(matrix[j,k])`, and
`im(matrix[j,k])` in the same order.  Reconstruction is the exact linear
combination of those coordinates and basis matrices and must round-trip to the
input under `sp.simplify`.

- [ ] **Step 4: Implement automatic full-image decomposition**

Add:

```python
@dataclass(frozen=True)
class MatrixGaugeDecomposition:
    hamiltonian: sp.ImmutableMatrix
    defect: sp.ImmutableMatrix
    include_global_phase: bool
    include_time_calibration: bool
    gauge: GaugeDecomposition
    projection_matrix: sp.ImmutableMatrix
    residual_matrix: sp.ImmutableMatrix
    witness_matrix: sp.ImmutableMatrix | None


def decompose_matrix_spectral_gauge(
    hamiltonian: sp.MatrixBase,
    defect: sp.MatrixBase,
    *,
    include_global_phase: bool = True,
    include_time_calibration: bool = True,
) -> MatrixGaugeDecomposition:
```

Generate every `i * (K*H - H*K)` from the complete Hermitian basis, append
`I` and `H` according to the two flags, and call `decompose_gauge` with:

```python
scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
completeness_id="exact-full-hermitian-basis-v1",
metric=(1,) * d + (2,) * (d * (d - 1)),
```

Reconstruct the projection, residual, and primitive witness matrices.  Verify
the witness commutes with `H` whenever the result is a spectral obstruction.

- [ ] **Step 5: Run all Task 1--2 tests**

Run:

```bash
python -m pytest -q tests/test_spectral_gauge.py
```

Expected: all tests pass.

- [ ] **Step 6: Commit Task 2**

```bash
git add \
  tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/spectral_gauge.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_spectral_gauge.py
git commit -m "feat(issue128): construct complete matrix gauge image"
```

---

### Task 3: Fail-Closed Gauge-Aware Issue-128 Audit

**Files:**
- Create: `scripts/audit_gauge_aware_obstruction.py`
- Create: `tests/test_gauge_aware_obstruction_audit.py`
- Create: `docs/experiments/processor-obstruction/gauge-aware-audit.json`
- Modify: `docs/experiments/processor-obstruction/README.md`

**Interfaces:**
- Consumes: `docs/experiments/processor-obstruction/exact-obstruction.json` and Task 2's exact matrix examples.
- Produces: `build_payload(source: Path) -> dict[str, object]`, `verify_payload(payload, source) -> None`, and CLI options `--source`, `--output`, `--verify`.

- [ ] **Step 1: Write failing artifact and mutation tests**

```python
import copy
import json
from pathlib import Path

import pytest

from scripts.audit_gauge_aware_obstruction import build_payload, verify_payload


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/experiments/processor-obstruction/exact-obstruction.json"


def test_current_obstruction_is_fixed_time_only() -> None:
    payload = build_payload(SOURCE)
    assert payload["fixed_time_endpoint_processor"]["status"] == "no_go"
    assert payload["calibrated_spectral_obstruction"]["status"] == "inconclusive"
    assert payload["reference_examples"]["hamiltonian_parallel"]["calibrated_status"] == "removable"
    assert payload["reference_examples"]["independent_commutant"]["status"] == "spectral_obstruction"
    verify_payload(payload, SOURCE)


def test_audit_rejects_forged_status_and_digest() -> None:
    payload = build_payload(SOURCE)
    forged = copy.deepcopy(payload)
    forged["calibrated_spectral_obstruction"]["status"] = "spectral_no_go"
    with pytest.raises(ValueError, match="calibrated"):
        verify_payload(forged, SOURCE)
    forged = copy.deepcopy(payload)
    forged["source"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="digest"):
        verify_payload(forged, SOURCE)
```

- [ ] **Step 2: Run the test and confirm missing-module failure**

Run:

```bash
python -m pytest -q tests/test_gauge_aware_obstruction_audit.py
```

Expected: collection fails for the missing audit module.

- [ ] **Step 3: Implement canonical payload construction and verification**

The payload must contain these top-level fields:

```python
{
    "schema_version": 1,
    "kind": "issue128_gauge_aware_processor_obstruction_audit",
    "source": {"path": str, "sha256": str},
    "allowed_calibration_gauge": ["image(i ad_H)", "span(I)", "span(H)"],
    "fixed_time_endpoint_processor": {
        "status": "no_go",
        "witness": "H",
        "source_nonzero": True,
    },
    "calibrated_spectral_obstruction": {
        "status": "inconclusive",
        "reason": "the attached witness H lies in the allowed time-calibration direction",
    },
    "restricted_support_processor": {
        "status": "operator_no_go_only",
        "spectral_implication": "not_established",
    },
    "reference_examples": {
        "off_diagonal": {"status": "removable"},
        "hamiltonian_parallel": {
            "fixed_time_status": "spectral_obstruction",
            "calibrated_status": "removable",
        },
        "independent_commutant": {
            "status": "spectral_obstruction",
            "orthogonal_to_identity": True,
            "orthogonal_to_hamiltonian": True,
            "commutes_with_hamiltonian": True,
        },
    },
}
```

`build_payload` must read the source's exact cubic overlap and require its
`nonzero` flag.  `verify_payload` must rebuild the entire payload and compare
canonical JSON bytes; it must never trust submitted status strings.

- [ ] **Step 4: Generate and verify the canonical artifact**

Run:

```bash
PYTHONPATH=src python -u -m scripts.audit_gauge_aware_obstruction \
  --source docs/experiments/processor-obstruction/exact-obstruction.json \
  --output docs/experiments/processor-obstruction/gauge-aware-audit.json

PYTHONPATH=src python -u -m scripts.audit_gauge_aware_obstruction \
  --source docs/experiments/processor-obstruction/exact-obstruction.json \
  --output docs/experiments/processor-obstruction/gauge-aware-audit.json \
  --verify
```

Expected: the first command writes sorted canonical JSON; the second prints a
single success line containing the source digest and exits zero.

- [ ] **Step 5: Document the calibrated claim boundary**

Append a `## Calibration-aware refinement` section to the existing README.
State explicitly:

```text
The exact <H,E5> witness proves the fixed-target-time endpoint no-go.  Because
H is itself an allowed time-calibration direction, this witness alone does not
prove a retiming-robust spectral obstruction.  The current calibrated status
is INCONCLUSIVE.  The support-six result remains a restricted local-processor
operator no-go and is not promoted to an eigenphase lower bound.
```

Link the new JSON digest and give both reproduction commands.

- [ ] **Step 6: Run focused tests**

Run:

```bash
python -m pytest -q \
  tests/test_spectral_gauge.py \
  tests/test_gauge_aware_obstruction_audit.py
```

Expected: all tests pass.

- [ ] **Step 7: Commit Task 3**

```bash
git add \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/audit_gauge_aware_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_gauge_aware_obstruction_audit.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/gauge-aware-audit.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/README.md
git commit -m "feat(issue128): audit calibrated spectral obstruction"
```

---

### Task 4: Regression and Claim Audit

**Files:**
- Modify only if a test exposes a defect: files introduced in Tasks 1--3.

**Interfaces:**
- Consumes: all Task 1--3 code and the existing Issue-128 non-slow suite.
- Produces: a clean focused result, a clean full non-slow result, and a final no-HPC decision.

- [ ] **Step 1: Run source-level static checks**

Run:

```bash
python -m compileall -q src/trottercert/spectral_gauge.py \
  scripts/audit_gauge_aware_obstruction.py
```

Expected: exit zero with no output.

- [ ] **Step 2: Run the focused scientific suite**

Run:

```bash
python -m pytest -q \
  tests/test_spectral_gauge.py \
  tests/test_gauge_aware_obstruction_audit.py \
  tests/test_effective_processor.py \
  tests/test_spectral_processing.py
```

Expected: all tests pass and no test labels a restricted span as a spectral
obstruction.

- [ ] **Step 3: Run the complete non-slow Issue-128 suite**

Run:

```bash
python -m pytest -q
```

Expected: all selected tests pass; slow tests remain deselected by project
configuration.

- [ ] **Step 4: Audit changed-file scope and artifact reproducibility**

Run:

```bash
git diff --check
git status --short
PYTHONPATH=src python -u -m scripts.audit_gauge_aware_obstruction \
  --source docs/experiments/processor-obstruction/exact-obstruction.json \
  --output docs/experiments/processor-obstruction/gauge-aware-audit.json \
  --verify
```

Expected: no whitespace errors; the audit verifies; unrelated pre-existing
dirty files remain unstaged and untouched.

- [ ] **Step 5: Record the promotion decision**

Do not request HPC.  Record in the handoff that the next gate is either:

1. an exact commutant witness orthogonal to `I,H`; or
2. a rigorous processed local-log remainder leaving a plausible physical-E7
   budget at `s11@35` or `s10@39`.
