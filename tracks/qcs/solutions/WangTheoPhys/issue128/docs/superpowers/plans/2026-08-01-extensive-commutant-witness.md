# Extensive Commutant-Witness Theorem Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an exact multi-size certificate for the polynomial commutant witness, identify its alias-free extensive scaling, and emit the normalized finite-step remainder threshold without overstating the claim.

**Architecture:** Extend the existing sparse symplectic module with exact support-geometry and moment primitives. A dedicated CLI generates the expensive cubic E5 density once, evaluates several torus sizes, fits only exact stable formulas, and writes a deterministic artifact; a theorem note records which statements are verified identities versus proved counting relations.

**Tech Stack:** Python 3.12, `fractions.Fraction`, exact `Cubic` arithmetic, bit-packed symplectic Paulis, canonical JSON, SHA-256, pytest 8

## Global Constraints

- Work only under `tracks/qcs/solutions/WangTheoPhys/issue128` on branch `codex/issue128-hpc-six-hour`.
- Preserve unrelated dirty-worktree files and stage only paths listed below.
- Use normalized trace `tau(A)=Tr(A)/2^(L^2)` throughout.
- Reject, rather than fold, any density whose coordinates alias modulo the requested torus length.
- Label empirical stable-size formulas separately from analytic identities until a local counting proof is included.
- Keep `finite_step_status=inconclusive` unless a compatible logarithm branch and all-order remainder are supplied.
- Do not authorize or submit HPC work.

---

### Task 1: Exact Geometry and Witness Moments

**Files:**
- Modify: `src/trottercert/commutant_witness.py`
- Modify: `tests/test_commutant_witness.py`

**Interfaces:**
- Consumes: `CoordinateRegistry`, canonical cubic density, `SquareLattice`, exact rational Pauli maps.
- Produces: `density_aliases_on_torus`, `quadratic_witness_moments`, and `QuadraticWitnessMoments`.

- [x] **Step 1: Add failing support-geometry and moment tests**

Add tests that require:

```python
assert density_aliases_on_torus(registry, density, 4) is True
assert density_aliases_on_torus(registry, density, 6) is False

moments = quadratic_witness_moments(h, h2, tau_h_e5, tau_h2_e5)
assert moments.tau_h2 == Fraction(6)
assert moments.tau_h3 == Fraction(-3)
assert moments.h_coefficient == Fraction(1, 2)
assert moments.tau_w2 > 0
assert moments.tau_w_e5 != Cubic.zero()
assert moments.squared_normalized_pairing == moments.tau_w_e5**2 / moments.tau_w2
```

- [x] **Step 2: Run the focused tests and confirm missing-symbol failures**

Run:

```bash
PYTHONPATH=src python -m pytest -q tests/test_commutant_witness.py \
  -k 'alias or moments or normalized'
```

Expected: collection or assertion failure because the new interfaces do not exist.

- [x] **Step 3: Implement exact alias detection**

Implement:

```python
def density_aliases_on_torus(
    registry: CoordinateRegistry,
    density: Mapping[SymplecticPauli, object],
    length: int,
) -> bool:
    """Return whether any distinct coordinates in one Pauli term coincide mod L."""
```

Validate even `length >= 4`.  For each term, map every occupied coordinate to
`(x % length, y % length)` and return true if the number of images decreases.

- [x] **Step 4: Implement the exact moment dataclass**

Implement:

```python
@dataclass(frozen=True)
class QuadraticWitnessMoments:
    tau_h2: Fraction
    tau_h3: Fraction
    tau_h4: Fraction
    identity_coefficient: Fraction
    h_coefficient: Fraction
    tau_w2: Fraction
    tau_h_e5: Cubic
    tau_h2_e5: Cubic
    tau_w_e5: Cubic
    squared_normalized_pairing: Cubic

def quadratic_witness_moments(
    hamiltonian: Mapping[SymplecticPauli, Fraction],
    squared: Mapping[SymplecticPauli, Fraction],
    tau_h_e5: Cubic,
    tau_h2_e5: Cubic,
) -> QuadraticWitnessMoments:
    ...
```

Use Pauli orthogonality:

```text
tau_h2 = coefficient of I in H^2
tau_h3 = <H,H^2>
tau_h4 = <H^2,H^2>
c = -tau_h3/tau_h2
W = H^2 - tau_h2 I + c H
tau_w2 = tau_h4 - tau_h2^2 + 2 c tau_h3 + c^2 tau_h2
tau_w_e5 = tau_h2_e5 + c tau_h_e5
rho = tau_w_e5^2 / tau_w2
```

Reject nonpositive `tau_h2` or `tau_w2`.

- [x] **Step 5: Run focused tests**

Run:

```bash
PYTHONPATH=src python -m pytest -q tests/test_commutant_witness.py
```

Expected: all tests pass.

- [x] **Step 6: Commit Task 1**

```bash
git add src/trottercert/commutant_witness.py tests/test_commutant_witness.py
git commit -m "feat(issue128): add exact extensive witness moments"
```

---

### Task 2: Multi-Size Exact Scan and Artifact

**Files:**
- Create: `scripts/certify_extensive_commutant_witness.py`
- Create: `tests/test_extensive_commutant_witness.py`
- Create: `docs/experiments/processor-obstruction/extensive-commutant-witness.json`

**Interfaces:**
- Consumes: Task 1 primitives and `exact_log_e5_density`.
- Produces: `compute_size_record`, `build_extensive_payload`, `verify_extensive_payload`, and a canonical artifact for sizes `6,8,10,12` or the maximal alias-free subset.

- [x] **Step 1: Add failing payload tests**

Tests must require that the verifier:

```python
payload = build_extensive_payload(records)
verify_extensive_payload(payload)
assert payload["claim"]["leading_order_status"] == "no_go"
assert payload["claim"]["finite_step_status"] == "inconclusive"
assert payload["hpc_authorized"] is False
```

Mutating a size, moment, alias decision, stable ratio, source hash, or claim
must raise `ValueError`.

- [x] **Step 2: Confirm the tests fail before the script exists**

Run:

```bash
PYTHONPATH=src python -m pytest -q tests/test_extensive_commutant_witness.py
```

Expected: import failure for `scripts.certify_extensive_commutant_witness`.

- [x] **Step 3: Implement one-pass E5 multi-size evaluation**

Generate E5 once.  For each requested even length:

```python
lattice = SquareLattice(length)
h = heisenberg_symplectic_terms(lattice)
h2 = square_real_pauli_terms(h)
tau_h_e5 = pair_lifted_cubic_density(registry, e5, h, lattice)
tau_h2_e5 = pair_lifted_cubic_density(registry, e5, h2, lattice)
record = quadratic_witness_moments(h, h2, tau_h_e5, tau_h2_e5)
```

Reject aliased sizes before lifting.  Reconstruct `H` from all four matching
densities for every accepted size.

- [x] **Step 4: Implement canonical payload and verifier**

Store every rational as `[numerator, denominator]` and every cubic as three
rational pairs.  Include source hashes for `cubic_field.py`, `cubic_local.py`,
and `commutant_witness.py`.  Recompute all algebraic relationships and stable
per-cell comparisons in the lightweight verifier.  `--verify-full` must rerun
the expensive exact scan and compare the complete payload.

- [x] **Step 5: Generate and fully verify the artifact**

Run:

```bash
PYTHONPATH=src python -u -m scripts.certify_extensive_commutant_witness \
  --lengths 6 8 10 12
PYTHONPATH=src python -u -m scripts.certify_extensive_commutant_witness \
  --lengths 6 8 10 12 --verify-full
```

Expected: deterministic equality and a nonzero witness pairing for every
accepted size.  If a requested size aliases, record the rejection and rerun
with the alias-free subset; do not weaken alias validation.

- [x] **Step 6: Run focused tests and commit Task 2**

```bash
PYTHONPATH=src python -m pytest -q \
  tests/test_commutant_witness.py tests/test_extensive_commutant_witness.py
git add scripts/certify_extensive_commutant_witness.py \
  tests/test_extensive_commutant_witness.py \
  docs/experiments/processor-obstruction/extensive-commutant-witness.json
git commit -m "feat(issue128): certify extensive spectral obstruction"
```

---

### Task 3: Theorem Note and Finite-Step Gate

**Files:**
- Create: `docs/report/extensive-commutant-witness-theorem.md`
- Modify: `docs/experiments/processor-obstruction/README.md`
- Modify: `docs/experiments/processor-obstruction/gauge-aware-audit.json`
- Modify: `scripts/audit_gauge_aware_obstruction.py`
- Modify: `tests/test_gauge_aware_obstruction_audit.py`

**Interfaces:**
- Consumes: the verified multi-size artifact and exact normalized obstruction.
- Produces: a claim-level theorem note, a squared finite-step remainder threshold, and an audit bound to both witness artifacts.

- [x] **Step 1: Add failing audit tests for the extensive artifact**

Require the audit to include the extensive artifact SHA-256, accepted sizes,
stable per-cell pairing status, exact squared normalized obstruction, and:

```json
{
  "leading_order_status": "no_go",
  "finite_step_status": "inconclusive",
  "finite_step_gate": "dual_pairing_remainder"
}
```

Mutation or removal of the extensive artifact must fail.

- [x] **Step 2: Write the theorem note with explicit proof levels**

The note must separate:

1. algebraically proved identities (`W` commutes and is calibration-orthogonal);
2. exact verified finite-size identities;
3. stable formulas supported by all alias-free sizes but still awaiting a
   symbolic linked-cluster count; and
4. the sufficient finite-step remainder inequality.

It must not call item 3 a theorem until its counting derivation is included.

- [x] **Step 3: Update and regenerate the fail-closed audit**

Read and lightweight-verify both witness artifacts.  Bind their hashes and
copy only recomputed statuses.  Preserve `hpc_authorized=false`.

- [x] **Step 4: Run all verification**

```bash
PYTHONPATH=src python -m pytest -q \
  tests/test_commutant_witness.py \
  tests/test_extensive_commutant_witness.py \
  tests/test_gauge_aware_obstruction_audit.py \
  tests/test_spectral_gauge.py
PYTHONPATH=src python -m pytest -q -m 'not slow'
```

Expected: all focused and non-slow tests pass.

- [x] **Step 5: Commit Task 3**

```bash
git add docs/report/extensive-commutant-witness-theorem.md \
  docs/experiments/processor-obstruction/README.md \
  docs/experiments/processor-obstruction/gauge-aware-audit.json \
  scripts/audit_gauge_aware_obstruction.py \
  tests/test_gauge_aware_obstruction_audit.py
git commit -m "docs(issue128): state extensive spectral theorem gate"
```

## Self-Review

- Spec coverage: aliasing, moments, normalization, multi-size data, mutation
  resistance, proof-level separation, finite-step criterion, and HPC gate each
  map to an explicit task.
- Placeholder scan: no `TBD`, `TODO`, or unspecified implementation steps.
- Type consistency: Task 2 consumes exactly the dataclass and functions
  introduced in Task 1; Task 3 consumes the artifact produced by Task 2.
