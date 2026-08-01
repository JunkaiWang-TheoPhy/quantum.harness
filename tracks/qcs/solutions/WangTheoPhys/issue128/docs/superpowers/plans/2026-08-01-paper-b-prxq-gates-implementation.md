# Paper B PRX Quantum Gates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decide with machine-checkable evidence whether gauge-aware spectral obstructions and correction-cost tradeoffs support a finite-step PRX Quantum paper.

**Architecture:** Build Paper B as a primal-dual certificate system.  Exact algebra produces either a processor certificate or a normalized commutant witness; exact E5/E7/E9 contractions and a rational E11-and-higher envelope determine whether the leading obstruction survives at finite step.  Promotion remains fail-closed and separate from Paper A.

**Tech Stack:** Python 3.12, exact `Fraction` and cubic-field arithmetic, NumPy/SymPy only for discovery and small-matrix checks, pytest, canonical JSON, SHA-256 manifests, Slurm arrays, LaTeX/latexmk.

## Global Constraints

- Use the gauge `Im(i ad_H) + span{I,H}` and handle degenerate energy blocks.
- Leading-order evidence must never be labeled finite-step evidence.
- Exact E9 alone never changes `finite_step_status` from `inconclusive`.
- A positive finite-step claim requires a strictly positive rational signed-margin interval.
- The standard comparison includes symplectic, symmetric, and composite correctors with explicit repeated-block costs.
- Preserve a positive correctable example; do not publish an indiscriminate no-go.
- E9 HPC requires an immutable source commit, new run root, complete manifests, recovery artifacts, and independent shard reruns.
- D8 remains outside this plan.

---

## File Structure

### New files

- `src/trottercert/spectral_duality.py`: gauge projection, exact small-matrix primal/dual certificates, and witness normalization.
- `tests/test_spectral_duality.py`: degenerate-block, gauge, weak/strong duality, and positive-example tests.
- `src/trottercert/trace_obstruction.py`: PF2/PF4 exact trace identities and certificate records.
- `tests/test_trace_obstruction.py`: symbolic coefficient and random-matrix identity checks.
- `src/trottercert/tfim_obstruction.py`: TFIM exact Pauli counts and scaling bound.
- `tests/test_tfim_obstruction.py`: finite-size counting and formula tests.
- `scripts/certify_tfim_obstruction.py`: canonical TFIM family certificate builder/verifier.
- `docs/experiments/processor-obstruction/tfim-family-obstruction.json`: frozen family certificate.
- `src/trottercert/corrector_cost.py`: symplectic/symmetric/composite block accounting and telescoping residuals.
- `tests/test_corrector_cost.py`: exact cancellation and linear-cost tests.
- `scripts/compare_corrector_architectures.py`: guarantee-aware comparison of baseline, processed, corrected, symmetric, and composite formulas.
- `tests/test_corrector_comparison.py`: metric-label, cost, and missing-evidence tests.
- `benchmarks/paper-b/corrector-comparison.json`: frozen fair-comparison artifact.
- `scripts/reduce_dual_e9_run.py`: clean-checkout E9 reduction and independent-rerun audit wrapper.
- `tests/test_reduce_dual_e9_run.py`: coverage, hash, duplicate, and missing-shard mutations.
- `src/trottercert/finite_step_obstruction.py`: combine exact pairings and tail interval into a signed decision.
- `tests/test_finite_step_obstruction.py`: positive, negative, crossing-zero, and digest-drift tests.
- `scripts/certify_finite_step_obstruction.py`: build/verify canonical finite-step certificate.
- `docs/experiments/processor-obstruction/finite-step-obstruction.json`: promoted or explicitly inconclusive result.
- `docs/paper-b/main.tex`: Paper B manuscript root.
- `docs/paper-b/sections/duality.tex`: B1.
- `docs/paper-b/sections/cost.tex`: B2.
- `docs/paper-b/sections/pf2-pf4.tex`: B3/B4.
- `docs/paper-b/sections/tfim.tex`: B5.
- `docs/paper-b/sections/finite-step.tex`: B6.
- `docs/paper-b/sections/sharpness.tex`: positive and negative examples.
- `docs/paper-b/validate_gates.py`: hash-bound PRX Quantum promotion gate.

### Existing files to modify

- `src/trottercert/spectral_gauge.py`: share established gauge helpers without changing frozen semantics.
- `src/trottercert/commutant_witness.py`: expose witness records needed by duality normalization.
- `src/trottercert/dual_log_pairing.py`: expose exact pairing records to the finite-step combiner.
- `src/trottercert/dual_log_tail.py`: expose the verified centered tail interval.
- `src/trottercert/dual_manifest_pairing.py`: retain exact reachability pruning and manifest-source checks.
- `scripts/certify_dual_e9_pairing.py`: add no new scientific semantics; only canonical reduction metadata if required.
- `docs/experiments/processor-obstruction/README.md`: distinguish leading-order and finite-step claims.

---

### Task 1: Gauge-Aware Spectral Duality Core

**Files:**
- Create: `src/trottercert/spectral_duality.py`
- Create: `tests/test_spectral_duality.py`
- Modify: `src/trottercert/spectral_gauge.py`
- Modify: `src/trottercert/commutant_witness.py`

**Interfaces:**
- Produces: `energy_blocks(h: np.ndarray, *, atol: float) -> tuple[np.ndarray, ...]` for discovery tests only.
- Produces: `gauge_residual(h: np.ndarray, e: np.ndarray, *, remove_identity: bool = True, remove_retiming: bool = True) -> np.ndarray`.
- Produces: `check_dual_witness(h: np.ndarray, e: np.ndarray, w: np.ndarray, *, atol: float) -> dict[str, float | bool]`.
- Exact production certificates use rational trace moments and do not rely on floating eigensolver tolerances.

- [ ] **Step 1: Write failing degenerate-block and gauge tests**

```python
def test_commutator_has_zero_energy_block_diagonal() -> None:
    h = np.diag([0.0, 0.0, 2.0])
    q = np.array([[0, 1j, 2j], [-1j, 0, 3j], [-2j, -3j, 0]], dtype=complex)
    comm = 1j * (q @ h - h @ q)
    residual = gauge_residual(h, comm)
    assert np.linalg.norm(residual, ord=2) < 1e-12


def test_retiming_and_phase_are_removed() -> None:
    h = np.diag([-1.0, 2.0])
    e = 3.0 * np.eye(2) - 0.25 * h
    assert np.linalg.norm(gauge_residual(h, e), ord=2) < 1e-12
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_spectral_duality.py`

Expected: missing module or interface.

- [ ] **Step 3: Implement block pinching and gauge projection for discovery tests**

Group equal eigenvalues before zeroing off-block entries.  Orthogonally project
the block-diagonal part away from `I` and `H` in the Hilbert--Schmidt inner
product, solving the 2 by 2 Gram system and handling a singular `span{I,H}`.

- [ ] **Step 4: Add a normalized dual-witness check**

```python
def check_dual_witness(h, e, w, *, atol):
    return {
        "commutes": bool(np.linalg.norm(w @ h - h @ w, ord=2) <= atol),
        "identity_orthogonal": bool(abs(np.trace(w)) <= atol),
        "retiming_orthogonal": bool(abs(np.trace(w @ h)) <= atol),
        "pairing": float(np.real(np.trace(w @ e))),
        "trace_norm": float(np.linalg.norm(w, ord="nuc")),
    }
```

- [ ] **Step 5: Test the existing Heisenberg witness and tuned qubit**

For the tuned qubit, assert the gauge residual is below `1e-12`.  For the
committed Heisenberg trace moments, assert exact nonzero pairing after identity
and retiming orthogonality.

- [ ] **Step 6: Run focused and existing gauge tests**

Run:

```bash
python -m pytest -q tests/test_spectral_duality.py \
  tests/test_spectral_gauge.py tests/test_commutant_witness.py \
  tests/test_gauge_aware_obstruction_audit.py
```

- [ ] **Step 7: Commit the duality core**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/spectral_duality.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/spectral_gauge.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/commutant_witness.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_spectral_duality.py
git commit -m "feat(issue128): add gauge-aware spectral duality core"
```

### Task 2: PF2 and PF4 Exact Trace Identities

**Files:**
- Create: `src/trottercert/trace_obstruction.py`
- Create: `tests/test_trace_obstruction.py`

**Interfaces:**
- Produces: `pf2_trace_obstruction(a: np.ndarray, b: np.ndarray) -> Fraction | float` for exact symbolic fixtures or numerical discovery.
- Produces: `pf4_trace_quadratic_form(trace_c2, trace_cd, trace_d2, gamma) -> Fraction`.
- Certificate records contain `formula`, `assumptions`, `exact_coefficients`, `exception_conditions`, and `identity_digest`.

- [ ] **Step 1: Write the failing PF4 coefficient test**

```python
def test_pf4_trace_quadratic_coefficients() -> None:
    value = pf4_trace_quadratic_form(
        trace_c2=Fraction(6), trace_cd=Fraction(5), trace_d2=Fraction(3), gamma=Fraction(7, 11)
    )
    expected = Fraction(7, 11) * (Fraction(6, 2) + Fraction(14, 3) * 5 + Fraction(4, 3) * 3)
    assert value == expected
```

- [ ] **Step 2: Run the test and verify failure**

Run: `python -m pytest -q tests/test_trace_obstruction.py`

- [ ] **Step 3: Implement the exact quadratic form**

Use the established coefficients from the derivation; preserve the mixed
`Tr(CD)` term rather than assuming symmetry until the symmetry test is applied.

- [ ] **Step 4: Add random-matrix identity tests**

Use seeded 2 by 2 and 3 by 3 Hermitian integer matrices.  Compute both the BCH
trace pairing and the closed quadratic form, and assert absolute agreement
below `1e-10`.  Include a commuting pair whose obstruction vanishes.

- [ ] **Step 5: Encode exception conditions**

Return explicit enum values `commuting`, `symmetry_protected_positive`,
`mixed_term_indefinite`, or `exactly_correctable`; never collapse the last two
into a universal no-go.

- [ ] **Step 6: Run tests and commit**

```bash
python -m pytest -q tests/test_trace_obstruction.py tests/test_lie_series.py
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/trace_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_trace_obstruction.py
git commit -m "feat(issue128): certify PF2 and PF4 trace obstructions"
```

### Task 3: TFIM Family Certificate

**Files:**
- Create: `src/trottercert/tfim_obstruction.py`
- Create: `tests/test_tfim_obstruction.py`
- Create: `scripts/certify_tfim_obstruction.py`
- Create: `docs/experiments/processor-obstruction/tfim-family-obstruction.json`

**Interfaces:**
- Produces: `tfim_trace_moments(length: int, h: Fraction, j: Fraction, periodic: bool) -> dict[str, Fraction]`.
- Produces exact keys `trace_c2_over_d`, `trace_d2_over_d`, `trace_cd_over_d`, `trace_h2_over_d`, `obstruction_over_d`, `operator_lower_bound`.

- [ ] **Step 1: Write failing finite-size formula tests**

```python
@pytest.mark.parametrize("length", (4, 6, 8))
def test_tfim_even_periodic_trace_formulas(length: int) -> None:
    moments = tfim_trace_moments(length, Fraction(2), Fraction(3), periodic=True)
    assert moments["trace_c2_over_d"] == 128 * length * Fraction(2) ** 4 * Fraction(3) ** 2
    assert moments["trace_d2_over_d"] == 128 * length * Fraction(2) ** 2 * Fraction(3) ** 4
    assert moments["trace_cd_over_d"] == 0
    assert moments["trace_h2_over_d"] == length * (Fraction(2) ** 2 + Fraction(3) ** 2)
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_tfim_obstruction.py`

- [ ] **Step 3: Implement exact Pauli counting**

Construct the finite periodic Pauli maps, compute exact trace moments by
orthogonality, and compare them to the closed formula.  Do not fit coefficients
from the checked lengths.

- [ ] **Step 4: Implement the family certificate builder and verifier mode**

The CLI accepts `--build OUTPUT` or `--verify INPUT`.  The certificate stores
the exact Suzuki coefficient field, checked lengths `(4, 6, 8)`, formula
strings, exact rational coefficients, assumptions, and source hashes.

- [ ] **Step 5: Add mutation tests**

Mutate `trace_cd_over_d`, the `alpha^2` coefficient, one checked length, and the
source hash.  Verify each mutation fails with a distinct error.

- [ ] **Step 6: Build and verify the frozen family artifact**

```bash
python -m pytest -q tests/test_tfim_obstruction.py
PYTHONPATH=src python scripts/certify_tfim_obstruction.py \
  --build docs/experiments/processor-obstruction/tfim-family-obstruction.json
PYTHONPATH=src python scripts/certify_tfim_obstruction.py \
  --verify docs/experiments/processor-obstruction/tfim-family-obstruction.json
```

- [ ] **Step 7: Commit the TFIM theorem artifact**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/certify_tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/tfim-family-obstruction.json
git commit -m "feat(issue128): certify TFIM spectral obstruction family"
```

### Task 4: Telescoping and Corrector-Cost Accounting

**Files:**
- Create: `src/trottercert/corrector_cost.py`
- Create: `tests/test_corrector_cost.py`

**Interfaces:**
- Produces: `CorrectorArchitecture(kind: Literal["symplectic","symmetric","composite"], kernel_blocks: int, endpoint_blocks: int, repeated_blocks_per_step: int)`.
- Produces: `total_blocks(architecture: CorrectorArchitecture, steps: int) -> int`.
- Produces: `telescoping_defect_order(kind: str) -> int | None`.

- [ ] **Step 1: Write failing exact cost tests**

```python
def test_symplectic_endpoint_cost_is_additive() -> None:
    arch = CorrectorArchitecture("symplectic", kernel_blocks=31, endpoint_blocks=2, repeated_blocks_per_step=0)
    assert total_blocks(arch, 95) == 31 * 95 + 2


def test_symmetric_correction_cost_is_linear() -> None:
    arch = CorrectorArchitecture("symmetric", kernel_blocks=31, endpoint_blocks=0, repeated_blocks_per_step=2)
    assert total_blocks(arch, 95) == 33 * 95
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_corrector_cost.py`

- [ ] **Step 3: Implement strict architecture validation**

Reject booleans, negative block counts, unknown architecture kinds, and zero
steps.  Record whether the architecture is a spectrum-preserving conjugation.

- [ ] **Step 4: Add a two-step symbolic telescoping test**

Use noncommuting SymPy symbols to verify
`(P*S*P^{-1})^2 = P*S^2*P^{-1}` as a word identity and verify that
`(F*S)^2 != F^2*S^2` unless the exact commutation condition is imposed.

- [ ] **Step 5: Run tests and commit**

```bash
python -m pytest -q tests/test_corrector_cost.py
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/corrector_cost.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_corrector_cost.py
git commit -m "feat(issue128): formalize corrector telescoping costs"
```

### Task 5: Guarantee-Aware Corrector Comparison

**Files:**
- Create: `scripts/compare_corrector_architectures.py`
- Create: `tests/test_corrector_comparison.py`
- Create: `benchmarks/paper-b/corrector-comparison.json`

**Interfaces:**
- Consumes: frozen standard-Suzuki resource data, processed-kernel audits, corrected-formula coefficient/cost records, and `CorrectorArchitecture`.
- Produces: `compare_architectures(records: tuple[MethodRecord, ...]) -> dict[str, object]`.
- Every `MethodRecord` has exact fields `name`, `guarantee`, `metric`, `formula_changed`, `kernel_blocks`, `endpoint_blocks`, `repeated_blocks_per_step`, `steps`, `total_blocks`, `source`.

- [ ] **Step 1: Write failing metric-separation tests**

```python
def test_comparison_rejects_mixed_metrics_in_one_ratio() -> None:
    records = (
        method_record("strict", metric="operator_norm", guarantee="finite_step_strict"),
        method_record("spectral", metric="eigenphase", guarantee="asymptotic"),
    )
    with pytest.raises(ValueError, match="ratios require the same metric and guarantee"):
        compare_architectures(records)


def test_comparison_recomputes_total_blocks() -> None:
    record = method_record("symmetric", repeated_blocks_per_step=2, steps=95, total_blocks=1)
    with pytest.raises(ValueError, match="total block count mismatch"):
        compare_architectures((record,))
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_corrector_comparison.py`

- [ ] **Step 3: Implement guarantee-aware records**

Allow ratios only within an identical `(metric, guarantee, benchmark,
normalization)` group.  For cross-guarantee rows, emit descriptive data without
a speedup ratio.  Recompute total blocks with `total_blocks()` and reject a
submitted count that differs.

- [ ] **Step 4: Freeze the comparison set**

Include the published control, the 95-step strict candidate, standard Suzuki
S4, every committed processed-kernel audit that passes its verifier, and the
corrected/symmetric/composite architectures used in B2.  Every external record
must contain a DOI or arXiv identifier and a verbatim guarantee label from the
source audit.

- [ ] **Step 5: Build and test the artifact**

```bash
python -m pytest -q tests/test_corrector_comparison.py tests/test_corrector_cost.py
PYTHONPATH=src python scripts/compare_corrector_architectures.py \
  --output benchmarks/paper-b/corrector-comparison.json
```

Expected: strict operator-norm ratios are separated from asymptotic eigenphase
comparisons, and every total block count is recomputed.

- [ ] **Step 6: Commit the comparison artifact**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/compare_corrector_architectures.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_corrector_comparison.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-b/corrector-comparison.json
git commit -m "feat(issue128): compare corrector architectures by guarantee"
```

### Task 6: Exact E9 Production Gate and Reduction Audit

**Files:**
- Create: `scripts/reduce_dual_e9_run.py`
- Create: `tests/test_reduce_dual_e9_run.py`
- Modify: `scripts/certify_dual_e9_pairing.py`
- Modify: `docs/experiments/processor-obstruction/README.md`

**Interfaces:**
- Consumes: immutable manifest index, exactly 64 shard files, source hashes, and two independent reruns.
- Produces: canonical `dual-e9-pairing.json` with complete coverage and exact cubic sums.

- [ ] **Step 1: Write failing coverage mutations**

```python
def test_reducer_rejects_missing_shard(tmp_path: Path) -> None:
    run = build_fake_run(tmp_path, shard_count=64)
    (run / "shards/shard-017.json").unlink()
    with pytest.raises(ValueError, match="missing shard 17"):
        audit_e9_run(run)


def test_reducer_rejects_duplicate_ordinal(tmp_path: Path) -> None:
    run = build_fake_run(tmp_path, shard_count=64)
    duplicate_first_ordinal(run / "shards/shard-001.json")
    with pytest.raises(ValueError, match="duplicate group ordinal"):
        audit_e9_run(run)
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_reduce_dual_e9_run.py`

- [ ] **Step 3: Implement the audit wrapper**

Require 64 source-bound shards, ordinal coverage `0..65535`, forward/reverse
sum equality, byte-identical index digest, and independent reruns for one median
and one heavy manifest.

- [ ] **Step 4: Run the local pre-HPC gate**

```bash
python -m pytest -q tests/test_dual_word_manifest.py tests/test_dual_manifest_pairing.py \
  tests/test_dual_e9_pairing.py tests/test_dual_e9_hpc_scripts.py \
  tests/test_reduce_dual_e9_run.py
bash -n ../../../../../hpc/issue128_e9_manifest.sbatch
bash -n ../../../../../hpc/issue128_e9_array.sbatch
bash -n ../../../../../hpc/issue128_e9_reduce.sbatch
```

- [ ] **Step 5: Commit the immutable production source**

Commit only E9 source, tests, Slurm scripts, and runbook.  Record the commit in
the run-root metadata before submission.

- [ ] **Step 6: Execute the documented Slurm run**

Use a new explicit `ISSUE128_E9_RUN_ROOT`.  Preserve scheduler IDs and accounting.
Do not edit manifests or shard outputs in place; missing work is rerun to a new
output path.

- [ ] **Step 7: Download and audit the complete result in a clean checkout**

Run:

```bash
python scripts/prepare_dual_log_manifests.py --verify "$ISSUE128_E9_RUN_ROOT/manifests/index.json"
python scripts/reduce_dual_e9_run.py --run-root "$ISSUE128_E9_RUN_ROOT" \
  --output "$ISSUE128_E9_RUN_ROOT/dual-e9-pairing.json"
python scripts/certify_dual_e9_pairing.py --index "$ISSUE128_E9_RUN_ROOT/manifests/index.json" \
  --verify "$ISSUE128_E9_RUN_ROOT/dual-e9-pairing.json"
```

- [ ] **Step 8: Commit only the reduced artifact and provenance**

Do not commit scheduler scratch or duplicated shard payloads to the main tree;
archive them externally and record their DOI/hash manifest.

### Task 7: Finite-Step Signed-Margin Certificate

**Files:**
- Create: `src/trottercert/finite_step_obstruction.py`
- Create: `tests/test_finite_step_obstruction.py`
- Create: `scripts/certify_finite_step_obstruction.py`
- Create: `docs/experiments/processor-obstruction/finite-step-obstruction.json`
- Modify: `src/trottercert/dual_log_pairing.py`
- Modify: `src/trottercert/dual_log_tail.py`

**Interfaces:**
- Produces: `decide_finite_step(q5: Cubic, q7: Cubic, q9: Cubic, root: RationalInterval, tail: RationalInterval, h: Fraction) -> FiniteStepDecision`; each cubic value is enclosed with `Cubic.enclose(root)` before interval combination.
- `FiniteStepDecision.status` is exactly `certified_obstruction`, `certified_no_margin`, or `inconclusive`.
- Produces exact fields `leading_interval`, `tail_interval`, `signed_margin_interval`, and source digests.

- [ ] **Step 1: Write failing decision tests**

```python
def test_positive_margin_certifies_obstruction() -> None:
    decision = decide_fixture(leading=(Fraction(10), Fraction(11)), tail=(Fraction(1), Fraction(2)))
    assert decision.status == "certified_obstruction"
    assert decision.signed_margin_interval.lower == 8


def test_crossing_zero_stays_inconclusive() -> None:
    decision = decide_fixture(leading=(Fraction(-1), Fraction(2)), tail=(Fraction(0), Fraction(1)))
    assert decision.status == "inconclusive"
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_finite_step_obstruction.py`

- [ ] **Step 3: Implement outward interval arithmetic**

Compute `h^4*q5 + h^6*q7 + h^8*q9` exactly in the cubic enclosure for the
effective local-log generator (the one-step logarithm divided by `h`), take a
sound absolute-value lower bound, subtract the tail upper bound, and promote
only when the margin lower endpoint is strictly positive.

- [ ] **Step 4: Bind every input digest and claim boundary**

The certificate must reject missing E9, mismatched E5/E7/E9 formula hashes,
wrong lattice normalization, boolean step counts, tail version drift, and an
attempt to change `status` without changing arithmetic.

- [ ] **Step 5: Add mutation tests**

Mutate each source digest, flip one exact cubic numerator, enlarge `h`, shrink
the tail by hand, and change `inconclusive` to `certified_obstruction`.  Every
mutation must fail verification.

- [ ] **Step 6: Build and verify the production decision**

```bash
python -m pytest -q tests/test_finite_step_obstruction.py \
  tests/test_dual_log_pairing.py tests/test_dual_log_tail.py \
  tests/test_dual_log_tail_certificate.py
PYTHONPATH=src python scripts/certify_finite_step_obstruction.py \
  --build docs/experiments/processor-obstruction/finite-step-obstruction.json
PYTHONPATH=src python scripts/certify_finite_step_obstruction.py \
  --verify docs/experiments/processor-obstruction/finite-step-obstruction.json
```

- [ ] **Step 7: Apply the promotion gate exactly once**

If `signed_margin_interval.lower > 0`, promote B6.  If the interval crosses
zero, retain `inconclusive`.  If one documented tail component dominates and a
precomputed best-case tightening can cross zero, authorize one focused
tightening; otherwise stop higher-degree expansion.

- [ ] **Step 8: Commit the finite-step certificate**

Use an explicit path list and commit:

`git commit -m "feat(issue128): certify finite-step spectral obstruction decision"`

### Task 8: Paper B Manuscript and Go/No-Go Artifact

**Files:**
- Create: `docs/paper-b/main.tex`
- Create: `docs/paper-b/sections/duality.tex`
- Create: `docs/paper-b/sections/cost.tex`
- Create: `docs/paper-b/sections/pf2-pf4.tex`
- Create: `docs/paper-b/sections/tfim.tex`
- Create: `docs/paper-b/sections/finite-step.tex`
- Create: `docs/paper-b/sections/sharpness.tex`
- Create: `docs/paper-b/validate_gates.py`
- Modify: `docs/experiments/processor-obstruction/README.md`
- Create: `docs/paper-b/prxq-gate.json`

**Interfaces:**
- Consumes only verified B1--B6 artifacts and corrector-cost data.
- Produces a gate record with booleans `b1_duality`, `b2_cost`, `family_theorem`, `b5_tfim`, `b6_finite_step`, `sharpness`, `exact_verifier`, `fair_comparison`.

- [ ] **Step 1: Write a gate validator**

The validator returns `PRX_QUANTUM` only when all eight booleans are true and
their evidence files exist with matching hashes.  Otherwise it returns
`QUANTUM_OR_PRR` and lists failed gates.

- [ ] **Step 2: Draft theorem sections from exact artifacts**

Each theorem states assumptions before conclusions and distinguishes
leading-order, finite-step, finite-size, and family-level claims.  The
finite-step section reads `finite-step-obstruction.json`; it never copies a
margin manually.

- [ ] **Step 3: Add the positive correctable example**

Include the tuned qubit exact processor, its zero quotient residual, and a
machine-readable certificate hash.  Contrast it with the TFIM/Heisenberg
witness rather than claiming universal failure.

- [ ] **Step 4: Build the manuscript and gate record**

```bash
python docs/paper-b/validate_gates.py
latexmk -pdf -interaction=nonstopmode -halt-on-error docs/paper-b/main.tex
rg -n "undefined|Overfull|Fatal error" docs/paper-b/*.log
```

- [ ] **Step 5: Run the Paper B full verification gate**

```bash
python -m pytest -q tests/test_spectral_duality.py tests/test_trace_obstruction.py \
  tests/test_tfim_obstruction.py tests/test_corrector_cost.py tests/test_corrector_comparison.py \
  tests/test_reduce_dual_e9_run.py tests/test_finite_step_obstruction.py
python scripts/certify_tfim_obstruction.py --verify \
  docs/experiments/processor-obstruction/tfim-family-obstruction.json
python scripts/certify_finite_step_obstruction.py --verify \
  docs/experiments/processor-obstruction/finite-step-obstruction.json
```

- [ ] **Step 6: Commit manuscript and gate decision**

If the PRXQ gate passes, commit with
`docs(issue128): assemble PRX Quantum spectral-obstruction manuscript`.
If it does not pass, commit the honest gate record with
`docs(issue128): record Paper B Quantum fallback scope`.

---

## Self-Review Checklist

- [ ] B1 handles degenerate energy blocks and removes both identity and retiming.
- [ ] PF2/PF4 code preserves exception conditions instead of forcing positivity.
- [ ] TFIM formulas are proved by exact counting rather than finite-size fitting.
- [ ] Corrector costs distinguish additive endpoints from repeated blocks.
- [ ] E9 requires all 64 shards and independent reruns.
- [ ] The finite-step decision remains inconclusive without verified E9 and tail inputs.
- [ ] Promotion uses a strictly positive lower endpoint, not a central estimate.
- [ ] A positive correctable example accompanies every no-go claim package.
- [ ] D8 is not introduced anywhere in this plan.
