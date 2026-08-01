# Affine Spectral Invariant Closure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade the exact finite-step dual local-log margin into a rigorous obstruction to the complete affine unitary orbit of the frozen Hamiltonian.

**Architecture:** Keep the dual signed-margin calculation as the first gate.  Add an exact rational centered-moment invariant whose linear term is the existing quadratic-witness pairing, bound every nonlinear term with the committed principal-log norm envelope, and expose a separate affine-spectral margin and status.  The production certificate regenerates all arithmetic from hash-bound E5, E7, E9, manifest, and tail inputs.

**Tech Stack:** Python 3.11+, `Fraction`, exact cubic-field intervals, canonical JSON, SHA-256, pytest.

## Global Constraints

- Never infer a finite-step spectral no-go from the tangent-gauge pairing alone.
- Use normalized trace moments `m2=54` and `m3=-27` for the frozen `12 x 12` instance.
- Convert the one-step log-defect cap to the effective-Hamiltonian cap by multiplying by exactly 97.
- Centering costs a factor of two: `epsilon=2*97*log_defect`.
- Promote only when both the dual signed-margin lower endpoint and the affine-invariant margin are strictly positive.
- Preserve exact rational endpoints in JSON; decimals are non-normative report text only.
- Any missing E9 manifest replay, source drift, lowered norm cap, or edited status fails closed.

---

### Task 1: Exact affine-invariant remainder kernel

**Files:**
- Modify: `src/trottercert/finite_step_obstruction.py`
- Modify: `tests/test_finite_step_obstruction.py`

**Interfaces:**
- Consumes: `dual_margin_lower: Fraction`, `cells: int`, `m2: Fraction`, `m3: Fraction`, `h_hs_cap: Fraction`, `effective_log_defect: Fraction`.
- Produces: `AffineSpectralDecision` with exact fields `centered_defect_cap`, `linear_invariant_lower`, `nonlinear_remainder_upper`, `invariant_margin`, and `status`.
- Produces: `decide_affine_spectral_obstruction(...) -> AffineSpectralDecision`.

- [ ] **Step 1: Write failing invariant and remainder tests**

```python
def test_affine_spectral_gate_uses_effective_log_cap() -> None:
    result = decide_affine_spectral_obstruction(
        dual_margin_lower=Fraction(1, 10**10),
        cells=36,
        m2=Fraction(54),
        m3=Fraction(-27),
        h_hs_cap=Fraction(15, 2),
        effective_log_defect=Fraction(1, 10**6),
    )
    assert result.centered_defect_cap == Fraction(1, 500000)
    assert result.nonlinear_remainder_upper > 0


def test_zero_affine_margin_stays_inconclusive() -> None:
    result = decide_affine_spectral_obstruction(
        dual_margin_lower=Fraction(0),
        cells=36,
        m2=Fraction(54),
        m3=Fraction(-27),
        h_hs_cap=Fraction(15, 2),
        effective_log_defect=Fraction(1, 10**6),
    )
    assert result.status == "inconclusive"
```

- [ ] **Step 2: Run the focused tests and verify failure**

Run: `PYTHONPATH=src python -m pytest -q tests/test_finite_step_obstruction.py`

Expected: import or attribute failure for the new interface.

- [ ] **Step 3: Implement exact remainder arithmetic**

Use

```python
epsilon = 2 * effective_log_defect
l2 = 2 * h_hs_cap * epsilon
n2 = epsilon**2
a2 = l2 + n2
l3 = 3 * m2 * epsilon
n3 = 3 * h_hs_cap * epsilon**2 + epsilon**3
a3 = l3 + n3
remainder = (
    m2**3 * (a3**2 + 2 * abs(m3) * n3)
    + m3**2 * (3 * m2**2 * n2 + 3 * m2 * a2**2 + a2**3)
)
linear = 6 * abs(m3) * m2**3 * cells * dual_margin_lower
margin = linear - remainder
```

Reject booleans, nonpositive cell counts, `m2<=0`, `m3=0`, a Hilbert--Schmidt cap below `sqrt(m2)` (check by squaring), negative margins, and negative defect caps.  Use status `certified_affine_spectral_obstruction` only when `margin>0`; otherwise use `inconclusive`.

- [ ] **Step 4: Add exact affine-orbit regression fixtures**

For rational diagonal eigenvalues, compute centered moments of `a+b*permutation(H)` and assert the cross-multiplied invariant is exactly zero.  Add a non-affine spectrum with nonzero invariant.  These tests document the global, not tangent-only, claim.

- [ ] **Step 5: Run and commit the kernel**

```bash
PYTHONPATH=src python -m pytest -q \
  tests/test_finite_step_obstruction.py tests/test_spectral_duality.py \
  tests/test_spectral_gauge.py
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/finite_step_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_finite_step_obstruction.py
git commit -m "feat(issue128): bound nonlinear spectral invariant"
```

### Task 2: Bind the nonlinear gate into the canonical certificate

**Files:**
- Modify: `scripts/certify_finite_step_obstruction.py`
- Modify: `tests/test_finite_step_obstruction_certificate.py`
- Modify: `artifacts/publication/paper-a-file-ownership.json`

**Interfaces:**
- Consumes: verified `envelope.log_defect`, `instance.steps`, `moments.h_hs_cap`, `moments.h_hs_squared`, and the positive dual margin.
- Produces: JSON sections `affine_spectral_invariant`, `effective_log_defect_bound`, `centered_defect_cap`, `linear_invariant_lower`, `nonlinear_remainder_upper`, and `invariant_margin`.
- Produces claim fields `local_log_status`, `finite_step_spectral_status`, and `promoted`.

- [ ] **Step 1: Write failing factor-97 and claim-forgery tests**

```python
def test_certificate_multiplies_one_step_log_defect_by_steps() -> None:
    payload = build_fixture()
    one_step = pair(payload["inputs"]["one_step_log_defect"])
    effective = pair(payload["affine_spectral_invariant"]["effective_log_defect_bound"])
    assert effective == 97 * one_step


def test_spectral_status_cannot_follow_local_log_status_without_margin() -> None:
    forged = build_fixture()
    forged["claim"]["finite_step_spectral_status"] = (
        "certified_affine_spectral_obstruction"
    )
    with pytest.raises(ValueError, match="regeneration mismatch"):
        verify_fixture(forged)
```

- [ ] **Step 2: Run and verify failure**

Run: `PYTHONPATH=src python -m pytest -q tests/test_finite_step_obstruction_certificate.py`

- [ ] **Step 3: Integrate exact tail geometry and the second decision**

Parse canonical rational pairs from the tail artifact, require `steps==97`,
`cells==36`, `h_hs_squared==54`, and `h_hs_cap==15/2`, compute
`effective_log_defect=steps*log_defect`, then call
`decide_affine_spectral_obstruction()` with the first gate's strict lower
endpoint.  Store every intermediate exact pair and regenerate the entire
payload during verification.

- [ ] **Step 4: Add mutation attacks**

Mutate the factor 97 to 1, halve the effective defect, lower the remainder,
change `m3`, flip only the public claim, and replace an E9 digest.  Require a
distinct verifier failure for each attack.

- [ ] **Step 5: Run the complete certificate suite and commit**

```bash
PYTHONPATH=src python -m pytest -q \
  tests/test_finite_step_obstruction.py \
  tests/test_finite_step_obstruction_certificate.py \
  tests/test_dual_e9_pairing.py tests/test_dual_log_tail_certificate.py
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/certify_finite_step_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_finite_step_obstruction_certificate.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-file-ownership.json
git commit -m "feat(issue128): gate finite-step affine spectrum"
```

### Task 3: Build production evidence and update the theorem boundary

**Files:**
- Create: `docs/experiments/processor-obstruction/dual-e9-pairing.json`
- Create: `docs/experiments/processor-obstruction/finite-step-obstruction.json`
- Create: `docs/report/finite-step-affine-spectral-ledger.md`
- Modify: `docs/experiments/processor-obstruction/README.md`
- Modify: `docs/report/dual-log-tail-ledger.md`
- Modify: `artifacts/publication/paper-a-file-ownership.json`

**Interfaces:**
- Consumes: fully audited 64-shard E9 archive, two independent reruns, and the frozen tail artifact.
- Produces: hash-bound reduced E9 and two-gate finite-step artifacts plus a human-readable theorem ledger.

- [ ] **Step 1: Verify and audit the complete E9 archive**

```bash
PYTHONPATH=src python scripts/prepare_dual_log_manifests.py \
  --verify "$RUN_ROOT/manifests/index.json"
PYTHONPATH=src python scripts/reduce_dual_e9_run.py \
  --run-root "$RUN_ROOT" --rerun "$RERUN_MEDIAN" --rerun "$RERUN_HEAVY" \
  --output "$RUN_ROOT/e9-run-audit.json"
PYTHONPATH=src python scripts/certify_dual_e9_pairing.py \
  --index "$RUN_ROOT/manifests/index.json" \
  --verify "$RUN_ROOT/dual-e9-pairing.json"
```

- [ ] **Step 2: Copy only canonical publication artifacts**

Copy the verified reduced E9 JSON into
`docs/experiments/processor-obstruction/dual-e9-pairing.json`.  Keep manifests,
workers, scheduler logs, and reruns in the external archive; record their
SHA-256 ledger and immutable source commit in the report.

- [ ] **Step 3: Build and verify the production two-gate certificate**

```bash
PYTHONPATH=src python scripts/certify_finite_step_obstruction.py \
  --build docs/experiments/processor-obstruction/finite-step-obstruction.json \
  --e9 docs/experiments/processor-obstruction/dual-e9-pairing.json \
  --index "$RUN_ROOT/manifests/index.json"
PYTHONPATH=src python scripts/certify_finite_step_obstruction.py \
  --verify docs/experiments/processor-obstruction/finite-step-obstruction.json \
  --e9 docs/experiments/processor-obstruction/dual-e9-pairing.json \
  --index "$RUN_ROOT/manifests/index.json"
```

- [ ] **Step 4: Write the theorem ledger from artifact fields**

State the affine-orbit invariant and remainder inequality before reporting the
two exact margins.  If the affine margin is not strictly positive, retain
`inconclusive` and describe only the local-log result.  Do not copy a decimal
margin into a claim-bearing file; the report cites JSON pointers and hashes.

- [ ] **Step 5: Run clean-checkout verification and commit**

```bash
PYTHONPATH=src python -m pytest -q \
  tests/test_reduce_dual_e9_run.py \
  tests/test_finite_step_obstruction.py \
  tests/test_finite_step_obstruction_certificate.py \
  tests/test_publication_scope.py
PYTHONPATH=src python scripts/certify_finite_step_obstruction.py \
  --verify docs/experiments/processor-obstruction/finite-step-obstruction.json \
  --e9 docs/experiments/processor-obstruction/dual-e9-pairing.json \
  --index "$RUN_ROOT/manifests/index.json"
git add tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/dual-e9-pairing.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/finite-step-obstruction.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/README.md \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/report/dual-log-tail-ledger.md \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/report/finite-step-affine-spectral-ledger.md \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-file-ownership.json
git commit -m "cert(issue128): close finite-step affine spectrum"
```

## Self-Review Checklist

- [ ] The invariant vanishes on the complete affine unitary orbit, not only its tangent space.
- [ ] The exact linear coefficient is `-162*54^3*tau(WD)`.
- [ ] The norm cap is multiplied by 97 before the nonlinear bound.
- [ ] Centering uses the explicit factor two.
- [ ] The spectral claim has its own strictly positive margin.
- [ ] E9 is verified against all manifests and two independent reruns.
- [ ] Local-log and affine-spectral statuses remain distinct in code, JSON, and prose.
