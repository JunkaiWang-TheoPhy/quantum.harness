# Processed Spectral Local Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a local, exact-rational promotion gate for the published `s8`, `s10`, and `s11` effective-order-six kernels, including processor order audits, spectral/QPE bookkeeping, and a hash-bound payload suitable for deciding whether a shared E7 HPC run is justified.

**Architecture:** Keep published coefficient provenance, exact free-word algebra, processor construction, and spectral certification in separate modules. Treat each published decimal string as the exact rational number it denotes, record the resulting nonzero order-condition residuals, and set the unneeded sixth-degree processor coordinate to zero. The output is a canonical JSON audit; it is a local proof input, not yet a physical Pauli E7 certificate.

**Tech Stack:** Python 3.12, standard-library `fractions`, `decimal`, `hashlib`, and `json`; SymPy 1.13 for exact sparse linear algebra; NumPy/SciPy only for explicitly labeled dense diagnostics; pytest 8.

## Global Constraints

- Preserve the current dirty worktree and do not modify unrelated existing files.
- The frozen physical benchmark remains periodic `12 x 12` isotropic Heisenberg, `T = 1`, operator-norm tolerance `10^-6`, split into four matching fragments.
- `s10` and `s11` are primary candidates; `s8` is the backup and processed-control candidate.
- Parse every published coefficient from its decimal source string into `Fraction`; never pass through binary floating point in a proof path.
- Bind provenance to arXiv `2404.04340`, Zenodo DOI `10.5281/zenodo.10814897`, source file `metemprar.py`, and MD5 `35a37af90096f371e7c82366cc7fad8b`.
- Carry rationalization residuals explicitly; do not label the rationalized kernels exact algebraic sixth-order formulas.
- The proof-only processor is `R(h) = h^2 R2 + h^4 R4`; set `R6 = 0` for the frozen local gate.
- Do not submit or prepare a new HPC job in this plan. Promotion to HPC requires both primary kernels to pass the exact and dense gates.
- Use canonical compact sorted JSON and rational `[numerator, denominator]` pairs for persistent proof data.
- Run focused tests first; run the full non-slow issue-128 suite only after all tasks pass.

---

### Task 1: Exact Published Kernel Registry

**Files:**
- Create: `src/trottercert/processed_kernels.py`
- Test: `tests/test_processed_kernels.py`

**Interfaces:**
- Consumes: no new project interface.
- Produces: `RationalStage`, `ProcessedKernel`, `published_effective_order_six_kernel(name)`, `alternating_kernel_stages(kernel, n_fragments=4)`, `merged_groups_per_step(kernel, n_fragments=4)`, and `maximum_steps_for_group_budget(kernel, group_budget, n_fragments=4)`.

- [x] **Step 1: Write the failing registry tests**

```python
from fractions import Fraction

import pytest

from trottercert.processed_kernels import (
    alternating_kernel_stages,
    maximum_steps_for_group_budget,
    merged_groups_per_step,
    published_effective_order_six_kernel,
)


@pytest.mark.parametrize(
    ("name", "stage_count", "beat_steps", "fivefold_steps"),
    [("s8", 49, 59, 49), ("s10", 61, 47, 39), ("s11", 67, 43, 35)],
)
def test_published_kernel_resource_gates(
    name: str, stage_count: int, beat_steps: int, fivefold_steps: int
) -> None:
    kernel = published_effective_order_six_kernel(name)
    stages = alternating_kernel_stages(kernel)
    assert len(stages) == stage_count
    assert sum((stage.coefficient for stage in stages), Fraction()) == 4
    assert merged_groups_per_step(kernel) == stage_count
    assert maximum_steps_for_group_budget(kernel, 2850) == beat_steps
    assert maximum_steps_for_group_budget(kernel, 2358) == fivefold_steps


def test_registry_rejects_unknown_kernel() -> None:
    with pytest.raises(ValueError, match="unknown processed kernel"):
        published_effective_order_six_kernel("s9")
```

- [x] **Step 2: Run the focused tests and confirm the missing-module failure**

Run: `python -m pytest -q tests/test_processed_kernels.py`

Expected: collection fails with `ModuleNotFoundError: trottercert.processed_kernels`.

- [x] **Step 3: Implement the exact registry and alternating four-fragment expansion**

```python
@dataclass(frozen=True)
class RationalStage:
    fragment_index: int
    coefficient: Fraction


@dataclass(frozen=True)
class ProcessedKernel:
    name: str
    seed_coefficients: tuple[Fraction, ...]
    source_strings: tuple[str, ...]
    source_doi: str
    source_file: str
    source_md5: str

    @property
    def composition_coefficients(self) -> tuple[Fraction, ...]:
        center = Fraction(1, 2) - sum(self.seed_coefficients, Fraction())
        half = self.seed_coefficients + (center,)
        return half + tuple(reversed(half))
```

Implement `alternating_kernel_stages` in the left-to-right product convention of Eq. (1.8): iterate the full coefficient list in reverse, use fragment order `(3,2,1,0)` for `chi` and `(0,1,2,3)` for `chi*`, then merge adjacent equal fragment indices exactly.

- [x] **Step 4: Run the registry tests**

Run: `python -m pytest -q tests/test_processed_kernels.py`

Expected: `4 passed`.

- [ ] **Step 5: Commit the registry**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/processed_kernels.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_processed_kernels.py
git commit -m "feat(issue128): register exact processed kernels"
```

### Task 2: Exact Free-Word and Lyndon Algebra

**Files:**
- Create: `src/trottercert/rational_words.py`
- Test: `tests/test_rational_words.py`

**Interfaces:**
- Consumes: `RationalStage` from Task 1.
- Produces: `Word`, `WordPolynomial`, `add_polynomials`, `scale_polynomial`, `multiply_polynomials`, `commutator_polynomial`, `rational_formula_log_series`, `lyndon_words`, `standard_lyndon_bracket`, and `word_l1`.

- [x] **Step 1: Write exact algebra tests**

```python
from fractions import Fraction

from trottercert.processed_kernels import RationalStage
from trottercert.rational_words import (
    commutator_polynomial,
    lyndon_words,
    rational_formula_log_series,
    standard_lyndon_bracket,
    word_l1,
)


def test_commutator_and_degree_four_lyndon_dimension() -> None:
    left = {(0,): Fraction(1)}
    right = {(1,): Fraction(1)}
    assert commutator_polynomial(left, right) == {
        (0, 1): Fraction(1),
        (1, 0): Fraction(-1),
    }
    basis = lyndon_words(4, 4)
    assert len(basis) == 60
    assert all(standard_lyndon_bracket(word) for word in basis)


def test_exact_strang_log_has_no_degree_two() -> None:
    stages = (
        RationalStage(0, Fraction(1, 2)),
        RationalStage(1, Fraction(1)),
        RationalStage(0, Fraction(1, 2)),
    )
    logarithm = rational_formula_log_series(stages, 3)
    assert word_l1(logarithm[2]) == 0
    assert word_l1(logarithm[3]) > 0
```

- [x] **Step 2: Run the tests and confirm the missing-module failure**

Run: `python -m pytest -q tests/test_rational_words.py`

Expected: collection fails with `ModuleNotFoundError: trottercert.rational_words`.

- [x] **Step 3: Implement sparse exact series multiplication**

Represent a homogeneous polynomial as `dict[tuple[int, ...], Fraction]`. Build the product formula degree by degree by multiplying each current series by one single-generator exponential:

```python
for target_degree in range(order, -1, -1):
    for power in range(1, target_degree + 1):
        scalar = stage.coefficient**power / factorial(power)
        suffix = (stage.fragment_index,) * power
        for word, coefficient in product[target_degree - power].items():
            updated[target_degree][word + suffix] += coefficient * scalar
```

Construct `log(I + X)` as `sum((-1)^(k+1) X^k/k, k=1..order)` using exact truncated multiplication.

- [x] **Step 4: Implement deterministic Lyndon enumeration and standard bracketing**

A word is Lyndon when it is strictly lexicographically smaller than every nontrivial rotation. For a Lyndon word of length greater than one, choose its longest proper Lyndon suffix `v`, write `w = uv`, and return `[bracket(u), bracket(v)]` via `commutator_polynomial`.

- [x] **Step 5: Run the exact algebra tests**

Run: `python -m pytest -q tests/test_rational_words.py`

Expected: `2 passed`.

- [ ] **Step 6: Commit the exact word algebra**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/rational_words.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_rational_words.py
git commit -m "feat(issue128): add exact free-word algebra"
```

### Task 3: Rational Effective-Order Processor Audit

**Files:**
- Create: `src/trottercert/effective_processor.py`
- Test: `tests/test_effective_processor.py`

**Interfaces:**
- Consumes: `ProcessedKernel`, `alternating_kernel_stages`, and exact polynomial operations from Tasks 1–2.
- Produces: `ProcessorCoordinates`, `EffectiveOrderAudit`, `least_squares_lie_image`, and `audit_effective_order_six(kernel)`.

- [x] **Step 1: Write the processor audit tests**

```python
from fractions import Fraction

import pytest

from trottercert.effective_processor import audit_effective_order_six
from trottercert.processed_kernels import published_effective_order_six_kernel


@pytest.mark.parametrize(
    ("name", "maximum_degree_three_l1", "maximum_degree_five_l1"),
    [
        ("s8", Fraction(1, 10**20), Fraction(1, 10**20)),
        ("s10", Fraction(1, 10**24), Fraction(1, 10**24)),
        ("s11", Fraction(1, 10**26), Fraction(1, 10**26)),
    ],
)
def test_rationalized_kernels_have_recorded_small_residuals(
    name: str,
    maximum_degree_three_l1: Fraction,
    maximum_degree_five_l1: Fraction,
) -> None:
    audit = audit_effective_order_six(
        published_effective_order_six_kernel(name)
    )
    assert audit.processed_degree_three_l1 > 0
    assert audit.processed_degree_five_l1 > 0
    assert audit.processed_degree_three_l1 < maximum_degree_three_l1
    assert audit.processed_degree_five_l1 < maximum_degree_five_l1
    assert len(audit.r2.coordinates) == 6
    assert len(audit.r4.coordinates) == 60
    assert audit.processed_degree_seven_l1 > 0
```

- [x] **Step 2: Run the tests and confirm the missing-module failure**

Run: `python -m pytest -q tests/test_effective_processor.py`

Expected: collection fails with `ModuleNotFoundError: trottercert.effective_processor`.

- [x] **Step 3: Implement exact least-squares projection onto a Lie-image space**

Build a sparse row map over the union of word supports. For basis polynomials `L_j`, solve the exact normal equations

```text
(M^T M) x = M^T target
```

with `sympy.linsolve`. Set any free symbols to zero and convert every result back to `Fraction`. Return both the coordinate tuple and the exact residual polynomial; never discard a rationalization residual.

- [x] **Step 4: Implement the degree-three and degree-five processor equations**

Use `A = log_kernel[1]`, `B3 = log_kernel[3]`, and `B5 = log_kernel[5]`. Solve

```text
B3 + [R2,A] = residual3
B5 + [R4,A] + [R2,B3] + 1/2 [R2,[R2,A]] = residual5
```

in the standard degree-two and degree-four Lyndon bases. Set `R6 = 0` and construct the processed degree-seven polynomial as

```text
B7 + [R2,B5] + [R4,B3]
 + 1/2 [R2,[R2,B3]]
 + 1/2 [R2,[R4,A]]
 + 1/2 [R4,[R2,A]]
 + 1/6 [R2,[R2,[R2,A]]].
```

- [x] **Step 5: Run the processor audit tests with timing**

Run: `/usr/bin/time -p python -m pytest -q tests/test_effective_processor.py`

Expected: `3 passed`; wall time below ten minutes and peak memory below 16 GiB. If it exceeds either local threshold, stop and optimize sparse normal-equation construction before continuing.

- [ ] **Step 6: Commit the processor audit**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/effective_processor.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_effective_processor.py
git commit -m "feat(issue128): audit effective-order processors exactly"
```

### Task 4: Spectral and QPE Bound Primitives

**Files:**
- Create: `src/trottercert/spectral_processing.py`
- Test: `tests/test_spectral_processing.py`

**Interfaces:**
- Consumes: no exact-word implementation details.
- Produces: `SpectralPhaseCertificate`, `spectral_phase_certificate(one_step_error, repetitions)`, and `target_projector_loss(processor_distance, effective_hamiltonian_error, gap)`.

- [x] **Step 1: Write boundary and monotonicity tests**

```python
from fractions import Fraction
from math import asin

import pytest

from trottercert.spectral_processing import (
    spectral_phase_certificate,
    target_projector_loss,
)


def test_repeated_spectral_phase_certificate() -> None:
    result = spectral_phase_certificate(Fraction(1, 1000), 35)
    assert result.repeated_error == Fraction(35, 1000)
    assert result.phase_radius == pytest.approx(2 * asin(0.035 / 2))


def test_phase_certificate_and_projector_loss_fail_closed() -> None:
    with pytest.raises(ValueError, match="below two"):
        spectral_phase_certificate(Fraction(1, 10), 20)
    assert target_projector_loss(Fraction(1, 100), Fraction(1, 1000), Fraction(1, 5)) == Fraction(3, 100)
    with pytest.raises(ValueError, match="positive gap"):
        target_projector_loss(Fraction(), Fraction(), Fraction())
```

- [x] **Step 2: Run the tests and confirm the missing-module failure**

Run: `python -m pytest -q tests/test_spectral_processing.py`

Expected: collection fails with `ModuleNotFoundError: trottercert.spectral_processing`.

- [x] **Step 3: Implement fail-closed theorem bookkeeping**

Store the exact repeated norm bound `min(2, repetitions * one_step_error)` only when the uncapped value is strictly below two; otherwise raise. Return the floating phase radius `2*asin(float(repeated_error)/2)`. Define the conservative projector loss as

```text
2 * processor_distance + 2 * effective_hamiltonian_error / gap.
```

Reject negative inputs and nonpositive gaps.

- [x] **Step 4: Run the spectral primitive tests**

Run: `python -m pytest -q tests/test_spectral_processing.py`

Expected: `2 passed`.

- [ ] **Step 5: Commit the spectral primitives**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/spectral_processing.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_spectral_processing.py
git commit -m "feat(issue128): add processed spectral bound primitives"
```

### Task 5: Canonical Local-Gate Payload and Independent Verification

**Files:**
- Create: `scripts/build_processed_kernel_audit.py`
- Create: `scripts/verify_processed_kernel_audit.py`
- Create: `tests/test_processed_kernel_audit.py`
- Create: `artifacts/processed-kernel-audit/.gitkeep`

**Interfaces:**
- Consumes: Tasks 1–4.
- Produces: schema-v1 `issue128_processed_kernel_local_audit` canonical JSON and a verifier that recomputes coefficients, stages, resource gates, processor coordinates, residuals, and hashes from source code.

- [x] **Step 1: Write payload round-trip and corruption tests**

```python
from pathlib import Path

import pytest

from scripts.build_processed_kernel_audit import build_payload, write_payload
from scripts.verify_processed_kernel_audit import verify_payload


def test_processed_kernel_payload_round_trip_and_corruption(tmp_path: Path) -> None:
    path = tmp_path / "audit.json"
    payload = build_payload(("s10", "s11"))
    write_payload(path, payload)
    verified = verify_payload(path)
    assert verified == ("s10", "s11")
    payload["kernels"][0]["stage_count"] = 60
    write_payload(path, payload)
    with pytest.raises(ValueError, match="stage count"):
        verify_payload(path)
```

- [x] **Step 2: Run the test and confirm the missing-module failure**

Run: `python -m pytest -q tests/test_processed_kernel_audit.py`

Expected: collection fails because the builder/verifier modules do not exist.

- [x] **Step 3: Implement the canonical builder**

For each kernel include:

```json
{
  "name": "s10",
  "source_strings": [],
  "composition_coefficients": [],
  "stage_count": 61,
  "resource_gates": {"beat_current_steps": 47, "fivefold_steps": 39},
  "processor": {"r2_lyndon": [], "r4_lyndon": [], "r6_policy": "zero"},
  "residuals": {"degree3": [], "degree5": []},
  "processed_degree7": [],
  "word_l1": {}
}
```

Encode every word polynomial as a lexicographically sorted list `[word_as_string, numerator, denominator]`. Include the source DOI/file/MD5, Python source SHA-256 digests, schema version, and explicit statement `rationalized_from_published_decimal_strings`.

- [x] **Step 4: Implement the recomputing verifier**

The verifier must reject schema, provenance, coefficient, order residual, stage count, resource arithmetic, word ordering, or source-hash drift. It must call the registry and processor audit rather than trusting submitted l1 totals.

- [x] **Step 5: Run the payload tests**

Run: `python -m pytest -q tests/test_processed_kernel_audit.py`

Expected: `1 passed`.

- [x] **Step 6: Emit and verify the frozen local payload**

Run:

```bash
PYTHONPATH=src python -u -m scripts.build_processed_kernel_audit \
  --kernels s10 s11 \
  --output artifacts/processed-kernel-audit/s10-s11-local-audit.json
PYTHONPATH=src python -u -m scripts.verify_processed_kernel_audit \
  artifacts/processed-kernel-audit/s10-s11-local-audit.json
```

Expected: both commands exit zero; verifier prints both candidate names, exact degree-three/five residual l1 values, processed degree-seven word l1 values, and resource gates.

- [ ] **Step 7: Commit the payload tooling and artifact**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_processed_kernel_audit.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/verify_processed_kernel_audit.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_processed_kernel_audit.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/processed-kernel-audit/s10-s11-local-audit.json
git commit -m "feat(issue128): freeze processed kernel local audit"
```

### Task 6: High-Precision Dense Spectral Discovery Gate

**Files:**
- Create: `scripts/screen_processed_spectra.py`
- Create: `tests/test_processed_spectral_screen.py`
- Create: `docs/experiments/processed-spectral-screen/README.md`

**Interfaces:**
- Consumes: Task 1 registry and stage convention.
- Produces: deterministic diagnostic JSON containing operator-error and eigenphase-error slopes for `s8`, `s10`, and `s11`. This output is explicitly non-rigorous and cannot enter a certificate bound.

- [x] **Step 1: Write a deterministic reduced-screen test**

```python
from scripts.screen_processed_spectra import run_random_fragment_screen


def test_processed_spectral_screen_separates_operator_and_phase_order() -> None:
    result = run_random_fragment_screen(
        kernels=("s10", "s11"), dimension=6, repetitions=(3, 4, 5, 6, 8), seed=128
    )
    assert 1.7 < result["s10"]["operator_slope"] < 2.3
    assert 5.3 < result["s10"]["phase_slope"] < 6.7
    assert 3.5 < result["s11"]["operator_slope"] < 4.5
    assert 5.3 < result["s11"]["phase_slope"] < 6.7
```

- [x] **Step 2: Run the test and confirm the missing-module failure**

Run: `python -m pytest -q tests/test_processed_spectral_screen.py`

Expected: collection fails because `scripts.screen_processed_spectra` does not exist.

- [x] **Step 3: Implement the deterministic dense screen**

Generate four seeded Hermitian fragments, normalize each to spectral norm `0.27`, and compare `kernel(1/r)^r` with `exp(-i sum(H_j))`. Match sorted eigenphases after unwrapping around the exact phases. Fit slopes by least squares in `log(r)` versus `log(error)` and emit all raw errors, not only fitted slopes.

- [x] **Step 4: Run the reduced test and full local screen**

Run:

```bash
python -m pytest -q tests/test_processed_spectral_screen.py
PYTHONPATH=src python -u -m scripts.screen_processed_spectra \
  --output docs/experiments/processed-spectral-screen/random-four-fragment.json
```

Expected: test passes; full screen completes locally in under ten minutes and preserves the qualitative split between low-order operator error and sixth-order phase error.

- [x] **Step 5: Document the diagnostic claim boundary**

The README must state that the screen checks stage order, phase matching, and candidate ranking only. It must not call the slopes rigorous, scalable, or a substitute for the `12 x 12` locality certificate.

- [ ] **Step 6: Commit the discovery gate**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/screen_processed_spectra.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_processed_spectral_screen.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processed-spectral-screen
git commit -m "test(issue128): add processed spectral discovery gate"
```

### Task 7: Local Promotion Decision and Regression Suite

**Files:**
- Create: `docs/experiments/processed-spectral-screen/promotion-decision.md`
- Modify: `docs/superpowers/plans/2026-07-31-processed-spectral-local-gate.md`

**Interfaces:**
- Consumes: all prior tasks.
- Produces: an explicit `promote`, `replace-s11-with-s8`, or `stop-before-hpc` decision with exact evidence and the next HPC plan boundary.

- [x] **Step 1: Run all focused tests**

Run:

```bash
python -m pytest -q \
  tests/test_processed_kernels.py \
  tests/test_rational_words.py \
  tests/test_effective_processor.py \
  tests/test_spectral_processing.py \
  tests/test_processed_kernel_audit.py \
  tests/test_processed_spectral_screen.py
```

Expected: all focused tests pass.

- [x] **Step 2: Run the full non-slow issue-128 suite**

Run: `python -m pytest -q`

Expected: no regressions relative to the pre-task baseline.

- [x] **Step 3: Write the promotion decision**

Record exact residual l1 values, processor coordinate l1 values, dense operator/phase slopes, stage budgets, runtime, peak memory, and one of these exact decisions:

- `promote s10+s11 to shared E7 design` when both pass;
- `promote s10+s8 to shared E7 design` when only `s11` fails the overlap/tail proxy;
- `stop before HPC` when neither pair preserves sixth-order phase behavior or the exact audit exceeds local feasibility.

- [x] **Step 4: Mark completed plan checkboxes**

Run: `rg -n '^- \[ \]' docs/superpowers/plans/2026-07-31-processed-spectral-local-gate.md`

Expected: only the intentionally deferred commit steps remain unchecked.

- [ ] **Step 5: Commit the decision**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processed-spectral-screen/promotion-decision.md \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/superpowers/plans/2026-07-31-processed-spectral-local-gate.md
git commit -m "docs(issue128): record processed kernel promotion gate"
```

The commit commands remain deliberately unexecuted because this is a shared,
already-dirty branch and the user did not request staging or commits.

## Self-Review

- Spec coverage: exact coefficient provenance, rational residuals, processor equations, spectral/QPE primitives, resource gates, canonical proof payload, dense discovery screen, and the pre-HPC decision are each assigned to a testable task.
- Deliberate exclusions: physical Pauli E7 evaluation, D6/D7 sidecars, a production finite-step locality theorem, Slurm scripts, manuscript edits, and any PRX/PRX Quantum claim are separate follow-on projects.
- Placeholder scan: the plan contains no `TBD`, deferred implementation marker, or unspecified test request.
- Type consistency: kernel names are strings; proof coefficients are `Fraction`; word polynomials are `dict[tuple[int, ...], Fraction]`; persistent rational values are `[numerator, denominator]` pairs.
