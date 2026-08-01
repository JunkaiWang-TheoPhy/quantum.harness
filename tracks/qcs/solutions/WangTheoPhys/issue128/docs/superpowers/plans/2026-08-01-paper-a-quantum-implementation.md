# Paper A Quantum Submission Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a clean-clone-reproducible Quantum submission and DOI release for the 95-step proof-carrying finite-step product-formula certificate.

**Architecture:** Preserve the frozen 95-step numerical result while adding a publication ownership manifest, claim-to-artifact audit, general soundness statements, nondegenerate exact-error calibration, preregistered XXZ transfer data, and a held-out model.  All paper numbers are generated or checked from hash-bound JSON/CSV artifacts, and the standard-library reference verifier remains independent of `trottercert`.

**Tech Stack:** Python 3.12, `fractions.Fraction`, NumPy/SciPy for small exact-matrix calibration, pytest, Ruff, LaTeX/latexmk, JSON/CSV, SHA-256, GitHub Actions.

## Global Constraints

- The authoritative result is 95 steps, 2,851 merged groups, and exact published-control ratio `11791/2851`.
- The strengthened-control ratio must be labeled separately and never replace the published-control denominator.
- Do not integrate an unfrozen 94-step bundle, D8 output, or a fivefold claim.
- Do not stage unrelated files from the existing dirty worktree; every commit uses an explicit path allowlist.
- `scripts/reference_verify.py` must not import `trottercert`.
- Any successful mutation stops release until every affected artifact and manifest is regenerated.
- Paper A submission does not wait for Paper B finite-step closure.

---

## File Structure

### New files

- `scripts/audit_publication_scope.py`: classify tracked and untracked Issue 128 files and validate allowed ownership labels.
- `tests/test_publication_scope.py`: schema, coverage, and overlap tests for the ownership manifest.
- `artifacts/publication/paper-a-file-ownership.json`: reviewed file classification.
- `artifacts/publication/paper-a-claim-matrix.json`: machine-readable claim, source, verifier, and manuscript-location table.
- `scripts/audit_paper_a_claims.py`: fail-closed claim and artifact audit.
- `tests/test_paper_a_claims.py`: claim-denominator, hash, and manuscript-drift mutations.
- `scripts/run_ed_calibration.py`: deterministic nondegenerate exact-error sweep.
- `tests/test_ed_calibration.py`: Hamiltonian, circuit, and output-schema tests.
- `benchmarks/paper-a/ed-calibration.json`: frozen exact-error data.
- `scripts/run_xxz_transfer.py`: preregistered XXZ parameter-family certificate sweep.
- `tests/test_xxz_transfer.py`: parameter, normalization, and determinism tests.
- `benchmarks/paper-a/xxz-transfer.json`: frozen transfer results.
- `scripts/run_heldout_model.py`: held-out TFIM certificate/calibration runner.
- `tests/test_heldout_model.py`: held-out freeze and output tests.
- `benchmarks/paper-a/heldout-tfim.json`: frozen held-out results.
- `scripts/build_paper_a_ablation.py`: recompute each norming layer from frozen inputs.
- `tests/test_paper_a_ablation.py`: monotonicity and final-row tests.
- `benchmarks/paper-a/ablation.json`: frozen ablation table.
- `docs/manuscript/sections/soundness.tex`: A1--A4 theorem statements and proofs.
- `docs/manuscript/figures/certificate_gap.pdf`: actual/certified error figure.
- `docs/manuscript/figures/model_transfer.pdf`: XXZ/held-out transfer figure.
- `docs/manuscript/data/certificate_gap.csv`: plotted exact-error data.
- `docs/manuscript/data/model_transfer.csv`: plotted transfer data.
- `artifacts/publication/paper-a-reproduction-transcript.txt`: clean-clone evidence.

### Existing files to modify

- `src/trottercert/crosscheck.py`: expose deterministic model construction and operator-error calculation used by the calibration script.
- `src/trottercert/hamiltonian.py`: add explicit XXZ and TFIM constructors if no equivalent public functions exist.
- `src/trottercert/verify.py`: accept only schema changes required by the transfer certificates; do not alter the frozen Heisenberg arithmetic.
- `scripts/crosscheck_small.py`: delegate to the reusable calibration interface.
- `scripts/package_delivery.py`: include the publication claim matrix and benchmark datasets.
- `docs/manuscript/main.tex`: include the soundness section and new figures.
- `docs/manuscript/sections/abstract.tex`: lead with proof-carrying capability and label both controls.
- `docs/manuscript/sections/introduction.tex`: state A1--A4 and scope limitations.
- `docs/manuscript/sections/results.tex`: add calibration, transfer, and ablation.
- `docs/manuscript/sections/reproducibility.tex`: document clean-clone and claim audit.
- `docs/manuscript/sections/limitations.tex`: forbid physical-minimal-step and universal-speedup interpretations.
- `docs/manuscript/scripts/generate_figures.py`: read frozen CSV/JSON only.
- `docs/manuscript/scripts/validate_claims.py`: invoke the publication claim audit.
- `.github/workflows/issue128-certificate.yml`: run fast/reference verification, mutation tests, claim audit, and manuscript data validation.

---

### Task 1: Publication Ownership Manifest

**Files:**
- Create: `scripts/audit_publication_scope.py`
- Create: `tests/test_publication_scope.py`
- Create: `artifacts/publication/paper-a-file-ownership.json`

**Interfaces:**
- Consumes: repository-relative paths under `tracks/qcs/solutions/WangTheoPhys/issue128`.
- Produces: `load_ownership(path: Path) -> dict[str, str]` and CLI exit status zero only when every listed path has one of `paper-a`, `paper-b`, `shared`, `exploratory`, `obsolete`.

- [ ] **Step 1: Write the failing ownership-schema tests**

```python
def test_ownership_rejects_unknown_label(tmp_path: Path) -> None:
    path = tmp_path / "ownership.json"
    path.write_text('{"README.md":"misc"}')
    with pytest.raises(ValueError, match="unknown ownership label"):
        load_ownership(path)


def test_ownership_rejects_duplicate_release_roles(tmp_path: Path) -> None:
    data = {"certificates/issue128-d5-integrated-certificate.json": ["paper-a", "obsolete"]}
    path = tmp_path / "ownership.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="exactly one ownership label"):
        load_ownership(path)
```

- [ ] **Step 2: Run the focused test and confirm it fails**

Run: `python -m pytest -q tests/test_publication_scope.py`

Expected: collection fails because `scripts.audit_publication_scope` does not exist.

- [ ] **Step 3: Implement the manifest loader and CLI**

```python
ALLOWED_LABELS = {"paper-a", "paper-b", "shared", "exploratory", "obsolete"}


def load_ownership(path: Path) -> dict[str, str]:
    raw = json.loads(path.read_text())
    result: dict[str, str] = {}
    for name, label in raw.items():
        if not isinstance(label, str):
            raise ValueError(f"{name}: exactly one ownership label is required")
        if label not in ALLOWED_LABELS:
            raise ValueError(f"{name}: unknown ownership label {label!r}")
        result[name] = label
    return result
```

- [ ] **Step 4: Populate the reviewed manifest**

At minimum classify the frozen certificate, D4/D5 sidecars, verifier scripts,
mutation tests, manuscript, E9 artifacts, processed-kernel experiments, D6,
D8, and all current plan documents.  Generate the initial path list with:

`git status --short && rg --files tracks/qcs/solutions/WangTheoPhys/issue128`

- [ ] **Step 5: Run schema and repository coverage tests**

Run: `python -m pytest -q tests/test_publication_scope.py`

Expected: all tests pass and no Paper A release file is labeled `exploratory`.

- [ ] **Step 6: Commit only ownership files**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/audit_publication_scope.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_publication_scope.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-file-ownership.json
git commit -m "chore(issue128): classify dual-track publication files"
```

### Task 2: Claim-to-Artifact Matrix and Fail-Closed Audit

**Files:**
- Create: `scripts/audit_paper_a_claims.py`
- Create: `tests/test_paper_a_claims.py`
- Create: `artifacts/publication/paper-a-claim-matrix.json`
- Modify: `docs/manuscript/scripts/validate_claims.py`

**Interfaces:**
- Consumes: certificate JSON, SHA-256 manifest, verification transcript, manuscript sources.
- Produces: `audit_claims(matrix_path: Path, root: Path) -> list[str]`; an empty list means every claim has a source, semantic checker, and manuscript location.

- [ ] **Step 1: Write failing tests for denominator and manuscript drift**

```python
def test_claim_audit_rejects_wrong_published_ratio(tmp_path: Path) -> None:
    matrix = make_matrix(tmp_path, ratio=[11791, 2851])
    matrix["claims"]["published_ratio"]["value"] = [10591, 2851]
    write_matrix(tmp_path, matrix)
    assert "published ratio must equal 11791/2851" in audit_claims(tmp_path / "matrix.json", tmp_path)


def test_claim_audit_rejects_unlabeled_strengthened_control(tmp_path: Path) -> None:
    matrix = make_matrix(tmp_path, ratio=[11791, 2851])
    matrix["claims"]["strengthened_ratio"]["label"] = "improvement"
    write_matrix(tmp_path, matrix)
    assert "strengthened control must be labeled" in audit_claims(tmp_path / "matrix.json", tmp_path)
```

- [ ] **Step 2: Verify the tests fail**

Run: `python -m pytest -q tests/test_paper_a_claims.py`

Expected: import failure for the missing audit script.

- [ ] **Step 3: Implement strict claim records**

Each claim record must contain exactly:

```json
{
  "value": [11791, 2851],
  "source": "certificates/issue128-d5-integrated-certificate.json",
  "checker": "scripts/reference_verify.py",
  "manuscript_locations": ["docs/manuscript/sections/abstract.tex"]
}
```

Reject missing files, hash mismatches, unlabeled controls, and forbidden phrases
such as `physically minimal`, `universal 4.1357`, or `certified fivefold`.

- [ ] **Step 4: Populate the matrix from the frozen certificate**

Include benchmark normalization, lattice size, tolerance, accepted and adjacent
steps, published and candidate groups, both labeled ratios, D4/D5 counts, test
evidence, and trusted-computing-base statements.

- [ ] **Step 5: Make manuscript validation invoke the audit**

```python
errors = audit_claims(PAPER_ROOT / "artifacts/publication/paper-a-claim-matrix.json", PAPER_ROOT)
if errors:
    raise SystemExit("\n".join(errors))
```

- [ ] **Step 6: Run tests and current manuscript validation**

Run:

```bash
python -m pytest -q tests/test_paper_a_claims.py
python docs/manuscript/scripts/validate_claims.py
```

Expected: both commands pass.

- [ ] **Step 7: Commit the claim audit**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/audit_paper_a_claims.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_paper_a_claims.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-claim-matrix.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/manuscript/scripts/validate_claims.py
git commit -m "test(issue128): bind paper claims to frozen artifacts"
```

### Task 3: Nondegenerate Exact-Error Calibration

**Files:**
- Create: `scripts/run_ed_calibration.py`
- Create: `tests/test_ed_calibration.py`
- Modify: `src/trottercert/crosscheck.py`
- Modify: `scripts/crosscheck_small.py`
- Create: `benchmarks/paper-a/ed-calibration.json`

**Interfaces:**
- Produces: `calibrate_rectangles(sizes: tuple[tuple[int,int], ...], steps: tuple[int, ...]) -> dict[str, object]`.
- Output fields: `model`, `width`, `height`, `steps`, `actual_operator_error`, `certified_upper`, `published_upper`, `normalization`, `source_commit`.

- [ ] **Step 1: Write failing deterministic calibration tests**

```python
def test_calibration_is_deterministic() -> None:
    left = calibrate_rectangles(((2, 3),), (8, 12))
    right = calibrate_rectangles(((2, 3),), (8, 12))
    assert left == right


def test_certificate_dominates_actual_error() -> None:
    data = calibrate_rectangles(((2, 3),), (16,))
    row = data["rows"][0]
    assert Decimal(row["actual_operator_error"]) <= Decimal(row["certified_upper"])
```

- [ ] **Step 2: Run tests and observe failure**

Run: `python -m pytest -q tests/test_ed_calibration.py`

Expected: missing calibration interface.

- [ ] **Step 3: Extract reusable exact-matrix functions**

Add to `crosscheck.py`:

```python
def operator_error(hamiltonian: np.ndarray, product_step: np.ndarray, *, time: float, steps: int) -> float:
    exact = scipy.linalg.expm(-1j * time * hamiltonian)
    approximate = np.linalg.matrix_power(product_step, steps)
    return float(np.linalg.norm(exact - approximate, ord=2))
```

Use the existing four-matching construction for `product_step`; do not create a
second normalization convention.

- [ ] **Step 4: Implement the frozen sweep**

Use sizes `((2, 3), (2, 4), (3, 3))` and steps
`(8, 12, 16, 24, 32, 48, 64, 96)`.  Sort rows by `(width, height, steps)` and
serialize decimal strings with 17 significant digits.

- [ ] **Step 5: Run focused and existing crosscheck tests**

Run:

```bash
python -m pytest -q tests/test_ed_calibration.py tests/test_small_crosscheck.py
PYTHONPATH=src python scripts/run_ed_calibration.py --output benchmarks/paper-a/ed-calibration.json
```

Expected: tests pass and every row satisfies `actual_operator_error <= certified_upper`.

- [ ] **Step 6: Commit calibration code and frozen data**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/crosscheck.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/crosscheck_small.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/run_ed_calibration.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_ed_calibration.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/ed-calibration.json
git commit -m "feat(issue128): add nondegenerate exact-error calibration"
```

### Task 4: Preregistered XXZ Transfer Sweep

**Files:**
- Create: `scripts/run_xxz_transfer.py`
- Create: `tests/test_xxz_transfer.py`
- Modify: `src/trottercert/hamiltonian.py`
- Create: `benchmarks/paper-a/xxz-transfer.json`

**Interfaces:**
- Produces: `xxz_bond(delta: Fraction) -> dict[str, Fraction]` with keys `XX`, `YY`, `ZZ` and values `1/4`, `1/4`, `delta/4`.
- Produces: `run_xxz_transfer(deltas: tuple[Fraction, ...]) -> dict[str, object]`.

- [ ] **Step 1: Write failing normalization and freeze tests**

```python
def test_xxz_bond_reduces_to_heisenberg() -> None:
    assert xxz_bond(Fraction(1)) == {"XX": Fraction(1, 4), "YY": Fraction(1, 4), "ZZ": Fraction(1, 4)}


def test_preregistered_deltas_are_exact() -> None:
    assert PREREGISTERED_DELTAS == (
        Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(1),
        Fraction(3, 2), Fraction(2), Fraction(4),
    )
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_xxz_transfer.py`

- [ ] **Step 3: Implement exact model construction and output schema**

Each row must contain exact delta, D4/D5 counts, norm bounds, accepted steps if
the existing schema supports them, builder time, verifier time, and a status
from `certified`, `unsupported`, or `inconclusive`.  Never silently omit a delta.

- [ ] **Step 4: Add a small-fixture deterministic test**

Patch the expensive builder with a three-term exact fixture and assert byte-
identical canonical JSON on two runs.

- [ ] **Step 5: Run the preregistered sweep**

Run:

```bash
python -m pytest -q tests/test_xxz_transfer.py tests/test_lattice_hamiltonian.py
PYTHONPATH=src python scripts/run_xxz_transfer.py --output benchmarks/paper-a/xxz-transfer.json
```

Expected: seven sorted rows are present, including unsuccessful statuses.

- [ ] **Step 6: Commit the transfer artifact**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/hamiltonian.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/run_xxz_transfer.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_xxz_transfer.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/xxz-transfer.json
git commit -m "feat(issue128): add preregistered XXZ transfer sweep"
```

### Task 5: Held-Out TFIM Evaluation

**Files:**
- Create: `scripts/run_heldout_model.py`
- Create: `tests/test_heldout_model.py`
- Modify: `src/trottercert/hamiltonian.py`
- Create: `benchmarks/paper-a/heldout-tfim.json`

**Interfaces:**
- Produces: `tfim_terms(length: int, field: Fraction, coupling: Fraction, periodic: bool) -> dict[SymplecticPauli, Fraction]`.
- Produces a frozen result with `rules_frozen_at_commit` and refuses to run when the ownership/claim schemas are dirty relative to that commit.

- [ ] **Step 1: Write the failing held-out freeze test**

```python
def test_heldout_run_requires_freeze_commit() -> None:
    with pytest.raises(ValueError, match="rules_frozen_at_commit is required"):
        run_heldout(field=Fraction(1), coupling=Fraction(1), rules_frozen_at_commit="")
```

- [ ] **Step 2: Run the test and confirm failure**

Run: `python -m pytest -q tests/test_heldout_model.py`

- [ ] **Step 3: Implement the exact TFIM constructor**

Use `-coupling * X_i X_{i+1}` and `-field * Z_i`, explicitly record sign and
normalization, sort by canonical Pauli key, and handle periodic closure only
when requested.

- [ ] **Step 4: Implement the held-out runner**

Run `length in (4, 6, 8)` and `steps in (8, 12, 16, 24, 32, 48)`.  Preserve
unsupported or inconclusive outcomes as rows with an explanatory enum.

- [ ] **Step 5: Run tests and freeze the held-out artifact**

```bash
python -m pytest -q tests/test_heldout_model.py tests/test_lattice_hamiltonian.py
PYTHONPATH=src python scripts/run_heldout_model.py \
  --rules-frozen-at-commit "$(git rev-parse HEAD)" \
  --output benchmarks/paper-a/heldout-tfim.json
```

- [ ] **Step 6: Commit the held-out result**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/hamiltonian.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/run_heldout_model.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_heldout_model.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/heldout-tfim.json
git commit -m "feat(issue128): add held-out TFIM evaluation"
```

### Task 6: Computed Deferred-Norming Ablation

**Files:**
- Create: `scripts/build_paper_a_ablation.py`
- Create: `tests/test_paper_a_ablation.py`
- Create: `benchmarks/paper-a/ablation.json`

**Interfaces:**
- Produces: `build_ablation(certificate: Path) -> dict[str, object]` with ordered stages `published`, `complete_formula`, `pauli_aggregation`, `translation_aggregation`, `d4_anticommuting`, `finite_ledger`.

- [ ] **Step 1: Write failing monotonicity and endpoint tests**

```python
def test_ablation_finishes_at_frozen_candidate() -> None:
    data = build_ablation(CERTIFICATE)
    assert data["stages"][-1]["steps"] == 95
    assert data["stages"][-1]["groups"] == 2851


def test_ablation_bounds_do_not_increase_after_exact_aggregation() -> None:
    data = build_ablation(CERTIFICATE)
    bounds = [Fraction(*row["bound"]) for row in data["stages"]]
    assert all(right <= left for left, right in zip(bounds, bounds[1:]))
```

- [ ] **Step 2: Run tests and verify failure**

Run: `python -m pytest -q tests/test_paper_a_ablation.py`

- [ ] **Step 3: Implement each ablation stage from frozen inputs**

Do not infer intermediate values by interpolation.  Each row records the exact
input artifact and function used to recompute its bound.

- [ ] **Step 4: Generate and verify the artifact**

```bash
python -m pytest -q tests/test_paper_a_ablation.py
PYTHONPATH=src python scripts/build_paper_a_ablation.py \
  --certificate certificates/issue128-d5-integrated-certificate.json \
  --output benchmarks/paper-a/ablation.json
```

- [ ] **Step 5: Commit ablation code and data**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_paper_a_ablation.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_paper_a_ablation.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/ablation.json
git commit -m "feat(issue128): compute deferred-norming ablation"
```

### Task 7: Soundness Section and Publication Figures

**Files:**
- Create: `docs/manuscript/sections/soundness.tex`
- Modify: `docs/manuscript/main.tex`
- Modify: `docs/manuscript/sections/abstract.tex`
- Modify: `docs/manuscript/sections/introduction.tex`
- Modify: `docs/manuscript/sections/results.tex`
- Modify: `docs/manuscript/sections/limitations.tex`
- Modify: `docs/manuscript/scripts/generate_figures.py`
- Create: `docs/manuscript/data/certificate_gap.csv`
- Create: `docs/manuscript/data/model_transfer.csv`
- Create: `docs/manuscript/figures/certificate_gap.pdf`
- Create: `docs/manuscript/figures/model_transfer.pdf`

**Interfaces:**
- Consumes only frozen benchmark JSON and the claim matrix.
- Produces deterministic CSV and PDF figures plus A1--A4 theorem text.

- [ ] **Step 1: Add manuscript-data tests to claim validation**

Assert that every plotted row occurs in its source JSON, all exact ratios use
integer numerator/denominator pairs, and both baseline labels appear in the
abstract or results.

- [ ] **Step 2: Run claim validation and record the expected failure**

Run: `python docs/manuscript/scripts/validate_claims.py`

Expected: failure because the new figure data and soundness section are absent.

- [ ] **Step 3: Write A1--A4 with explicit assumptions**

Use theorem labels `thm:certificate-soundness`, `thm:deferred-norming`,
`thm:finite-lattice-transfer`, and `thm:fail-closed`.  Each proof cites the
exact schema field and verifier obligation it discharges.

- [ ] **Step 4: Generate figures from frozen artifacts**

No numerical arrays may be embedded in `generate_figures.py`; the script reads
the three Paper A benchmark JSON files and writes canonical CSV before plotting.

- [ ] **Step 5: Build and inspect the manuscript**

Run:

```bash
python docs/manuscript/scripts/generate_figures.py
python docs/manuscript/scripts/validate_claims.py
latexmk -pdf -interaction=nonstopmode -halt-on-error docs/manuscript/main.tex
rg -n "undefined|Overfull|Fatal error" docs/manuscript/*.log
```

Expected: validation and LaTeX pass; the final `rg` returns no blocking issue.

- [ ] **Step 6: Commit theorem prose and generated evidence**

Use an explicit `git add` list containing only the files in this task, then:

`git commit -m "docs(issue128): present proof-carrying soundness and transfer evidence"`

### Task 8: CI, Clean-Clone Reproduction, and Release Bundle

**Files:**
- Modify: `.github/workflows/issue128-certificate.yml`
- Modify: `scripts/package_delivery.py`
- Modify: `docs/manuscript/sections/reproducibility.tex`
- Create: `artifacts/publication/paper-a-reproduction-transcript.txt`

**Interfaces:**
- CI consumes the frozen certificate and benchmark artifacts.
- Release packaging produces a directory containing manuscript, source, certificates, sidecars, verifiers, mutation corpus, data, transcript, `CITATION.cff`, license, and `SHA256SUMS`.

- [ ] **Step 1: Add CI assertions**

The workflow must run:

```yaml
- run: python scripts/verify.py certificates/issue128-d5-integrated-certificate.json
- run: python scripts/reference_verify.py certificates/issue128-d5-integrated-certificate.json
- run: python -m pytest -q tests/test_certificate_mutations.py tests/test_reference_verify.py
- run: python scripts/audit_paper_a_claims.py artifacts/publication/paper-a-claim-matrix.json
- run: python docs/manuscript/scripts/validate_claims.py
```

- [ ] **Step 2: Run the full local release gate**

```bash
python -m pytest -q
python scripts/verify.py certificates/issue128-d5-integrated-certificate.json
python scripts/verify.py --deep certificates/issue128-d5-integrated-certificate.json
python scripts/reference_verify.py certificates/issue128-d5-integrated-certificate.json
python docs/manuscript/scripts/validate_claims.py
git diff --check
```

Expected: all commands pass with 95 accepted and 94 rejected.

- [ ] **Step 3: Reproduce from a clean clone**

Clone the exact release commit into a new temporary directory, install from
`requirements-reproducibility.txt`, rerun the local release gate, and record
commit, Python, dependency hashes, wall time, peak RSS, and artifact hashes in
`paper-a-reproduction-transcript.txt`.

- [ ] **Step 4: Build and verify the release package**

Run `scripts/package_delivery.py` with an output directory outside the source
tree.  Run `shasum -a 256 -c SHA256SUMS` from inside the package and compare the
packaged manuscript PDF hash to the source release hash.

- [ ] **Step 5: Commit release automation and transcript**

Use an explicit path list and commit:

`git commit -m "release(issue128): freeze Quantum submission artifact"`

- [ ] **Step 6: Tag and archive only after remote CI passes**

Create `issue128-paper-a-v1.0.0`, mint the DOI archive, update the claim matrix
with the DOI, rerun claim validation, and create a final metadata-only commit.

---

## Self-Review Checklist

- [ ] Every Paper A completion-gate item maps to a task above.
- [ ] No task changes D8 or integrates an unfrozen 94-step result.
- [ ] Both control ratios are separately labeled and tested.
- [ ] The held-out model is run only after rules are frozen.
- [ ] Actual-error calibration includes nondegenerate systems.
- [ ] Every generated figure reads frozen data rather than embedded arrays.
- [ ] Release requires clean-clone evidence and independent reference verification.
