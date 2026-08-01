# Strengthened Rational-Pair Control Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine with a proof-carrying artifact whether the informal 353-step, 10,591-group Heisenberg control is valid, and publish the exact derived result even if it differs.

**Architecture:** Refactor the existing fully rational center-17 paired-theorem calculation so it emits a deterministic witness sidecar, then verify that sidecar with an independent standard-library implementation.  A canonical summary derives the accepted and adjacent rejected step counts and the `G(r)=30r+1` resources; no code path is allowed to request or force `r=353`.

**Tech Stack:** Python 3.12, `Fraction`, rational intervals, deterministic gzip/JSON, SHA-256, pytest.

## Global Constraints

- Hamiltonian: `H=sum_<ij>(XX+YY+ZZ)/4` on the 12 by 12 periodic square lattice.
- Formula: five-copy fourth-order Suzuki formula with 31 merged stages.
- Target: `T=1`, operator-norm tolerance `epsilon=10^-6`.
- Control method name: `rational_pair_theorem_control`; never call it `complete_formula_log`.
- Candidate theorem parameters: center 17 and 24 decimal digits of outward coefficient intervals.
- Resource rule: `G(r)=30*r+1` merged schedule exponentials.
- Promotion requires exact `E_r <= epsilon` and exact `E_(r-1) > epsilon`.
- If the derived `r` is not 353, freeze the actual value and mark the 353 hypothesis refuted.
- Do not modify or regenerate the authoritative 95-step D4/D5 certificate.
- Use explicit file allowlists for every commit; the worktree contains unrelated concurrent changes.

---

## File Structure

### New files

- `src/trottercert/strengthened_control.py`: strict spec, summary schema, exact step/resource derivation, and semantic verification.
- `scripts/certify_strengthened_control.py`: generator CLI for witness and summary artifacts.
- `scripts/reference_verify_strengthened_control.py`: independent standard-library verifier that does not import `trottercert`.
- `tests/test_strengthened_control.py`: exact derivation, schema, adjacent-step, and mutation tests.
- `tests/test_reference_strengthened_control.py`: subprocess attacks against the independent verifier.
- `benchmarks/paper-a/strengthened-control.json`: canonical summary.
- `benchmarks/paper-a/strengthened-control-witness.json.gz`: deterministic raw-composition multiset, aggregated theorem records, and pairing witness.

### Existing files to modify after the scientific gate passes

- `src/trottercert/rigorous_fourth.py`: expose deterministic witness records without changing the existing aggregate result.
- `benchmarks/paper-a/ablation.json`: add a separately named `rational_pair_theorem_control` row; leave `complete_formula` as `missing_evidence`.
- `scripts/build_paper_a_ablation.py`: rebuild the new row only from the certified artifact.
- `artifacts/publication/paper-a-claim-matrix.json`: include the strengthened ratio only if the artifact is certified.
- `artifacts/publication/paper-a-file-ownership.json`: classify code, tests, summary, and witness.
- `scripts/audit_publication_scope.py`: add only the promoted Paper A paths to the release allowlist.

---

### Task 1: Freeze the Control Spec and Result Semantics

**Files:**
- Create: `src/trottercert/strengthened_control.py`
- Test: `tests/test_strengthened_control.py`

**Interfaces:**
- Produces: `StrengthenedControlSpec` dataclass.
- Produces: `derive_control_summary(spec, theorem_record, sources) -> dict[str, object]`.
- Produces: `verify_control_summary(payload, witness_path, root) -> list[str]`.

- [ ] **Step 1: Write strict spec tests**

```python
def test_frozen_control_spec_is_exact() -> None:
    spec = StrengthenedControlSpec.frozen()
    assert spec.length == 12
    assert spec.time == Fraction(1)
    assert spec.tolerance == Fraction(1, 10**6)
    assert spec.center == 17
    assert spec.decimal_digits == 24
    assert spec.stage_count == 31
    assert spec.method == "rational_pair_theorem_control"


@pytest.mark.parametrize("field,value", [("length", True), ("center", 17.0)])
def test_spec_rejects_python_numeric_aliases(field: str, value: object) -> None:
    values = asdict(StrengthenedControlSpec.frozen())
    values[field] = value
    with pytest.raises((TypeError, ValueError)):
        StrengthenedControlSpec(**values)
```

- [ ] **Step 2: Run the tests and verify the missing module failure**

Run: `python -m pytest -q tests/test_strengthened_control.py`

Expected: collection fails because `trottercert.strengthened_control` does not exist.

- [ ] **Step 3: Implement the frozen dataclass and strict rational-pair parser**

The dataclass fields are `model`, `normalization`, `length`, `boundary`, `time`, `tolerance`, `formula`, `stage_count`, `center`, `decimal_digits`, `method`, and `resource_rule`.  Use `type(value) is int` for integer fields.  For each JSON rational pair `[n,d]`, require `d > 0`, construct `value = Fraction(n, d)`, and require the submitted pair to equal `[value.numerator, value.denominator]`; this rejects non-lowest-term and negative-denominator encodings.

- [ ] **Step 4: Add hypothesis-neutral result states**

Permit exactly `certified_matches_353`, `certified_differs_from_353`, and `inconclusive`.  The factory selects the first two from a fully verified accepted step; callers cannot pass a status string.  A controlled `inconclusive` record is reachable only when the exact theorem construction is valid but a frozen reason enum (`resource_limit_before_search`, `interval_precision_not_separating_adjacent_steps`) leaves `accepted_step=null`; it carries no ratio, never promotes, and may pass semantic verification as an honest negative result.  Algebra/schema/source-verification failure produces no artifact at all.

- [ ] **Step 5: Run focused tests**

Run: `python -m pytest -q tests/test_strengthened_control.py`

Expected: all Task 1 tests pass.

- [ ] **Step 6: Commit Task 1**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/strengthened_control.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_strengthened_control.py
git commit -m "feat(issue128): freeze strengthened control semantics"
```

### Task 2: Emit a Deterministic Rational-Pair Witness

**Files:**
- Modify: `src/trottercert/rigorous_fourth.py`
- Modify: `tests/test_rigorous_fourth.py`
- Modify: `tests/test_strengthened_control.py`

**Interfaces:**
- Produces: `fourth_order_rational_pair_witness(center=17, decimal_digits=24) -> RationalPairWitness`.
- `RationalPairWitness.aggregate()` must equal the existing `FourthOrderRationalCertificate` exactly.
- `RationalPairWitness.raw_terms` retains every weak-composition contribution before aggregation; `RationalPairWitness.records` retains the exact per-`(outer,j)` sums used by the legacy certificate.

- [ ] **Step 1: Write an aggregate-preservation test**

```python
@pytest.mark.slow
def test_rational_pair_witness_reproduces_existing_aggregate() -> None:
    aggregate = fourth_order_rational_pair_certificate(center=17, decimal_digits=24)
    witness = fourth_order_rational_pair_witness(center=17, decimal_digits=24)
    assert witness.aggregate() == aggregate
    assert witness.raw_term_count == len(witness.raw_terms)
    assert witness.aggregated_record_count == len(witness.records)
    assert witness.raw_term_count == aggregate.theorem_terms
```

- [ ] **Step 2: Define immutable witness records**

Each raw record contains the fragment-index word, weak composition, outer fragment word, one-based theorem index `j`, and the exact raw scalar interval/upper bound.  Preserve the complete raw multiset, including repeated `(outer,j)` keys.  Each aggregated record contains that key's exact summed scalar plus the common coefficient denominator, sorted Pauli masks, lower/upper integer coefficient numerators, exact pair indices, singleton indices, and exact local bound.  The top-level witness also contains `raw_term_count`, `aggregated_record_count`, the 31 exact interval stages, cube-root interval, center, digit policy, and factorial denominator.

- [ ] **Step 3: Refactor the existing loop to collect records**

Keep `_greedy_pairs` discovery-only.  Immediately recheck every recorded pair with `_anticommutes`; compute each aggregated scalar exclusively from its raw-record multiset, and compute the certificate aggregate exclusively by summing the aggregated exact local bounds.  The legacy aggregate function calls the witness builder and returns `witness.aggregate()` so the two paths cannot drift.

- [ ] **Step 4: Add determinism and coverage tests**

Assert sorted record keys, unique Pauli indices, disjoint pair/singleton coverage, positive denominators, and byte-identical canonical serialization from two builds.

- [ ] **Step 5: Run the bounded slow calculation after setup confirmation**

Run: `python -m pytest -q -m slow tests/test_rigorous_fourth.py::test_rational_pair_witness_reproduces_existing_aggregate`

Expected: pass; record wall time and peak RSS.  If wall time is projected above 10 minutes or RSS above 16 GB, stop and move this task to the configured Slurm profile.

- [ ] **Step 6: Commit Task 2**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/rigorous_fourth.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_rigorous_fourth.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_strengthened_control.py
git commit -m "feat(issue128): expose rational-pair control witness"
```

### Task 3: Derive Adjacent-Step Evidence and Resources

**Files:**
- Modify: `src/trottercert/strengthened_control.py`
- Modify: `tests/test_strengthened_control.py`

**Interfaces:**
- Consumes: `FourthOrderRationalCertificate.site_density_upper`.
- Produces exact `global_error_upper(r) = 144*site_density_upper/r^4`.
- Produces `accepted_steps`, `previous_step`, both exact error pairs, `group_exponentials`, and exact comparison ratios.

- [ ] **Step 1: Write derivation tests without assuming 353**

```python
def test_summary_derives_the_minimal_step_from_exact_density() -> None:
    theorem = fake_theorem(site_density_upper=Fraction(5, 2))
    summary = derive_control_summary(FROZEN_SPEC, theorem, SOURCES)
    r = summary["accepted_steps"]
    assert Fraction(*summary["accepted_error_upper"]) <= FROZEN_SPEC.tolerance
    assert Fraction(*summary["previous_error_upper"]) > FROZEN_SPEC.tolerance
    assert summary["previous_step"] == r - 1
    assert summary["group_exponentials"] == 30 * r + 1
```

- [ ] **Step 2: Implement exact derivation**

Use `required_steps(144*site_density_upper, tolerance, 4, time)` and recompute both inequalities directly with `Fraction`.  Store `hypothesized_steps=353` and `hypothesized_group_exponentials=10591` as comparison metadata, never as inputs to the derivation.

- [ ] **Step 3: Add tamper tests**

Mutate accepted step, previous step, either error numerator, resource count, status, ratio, or a boolean masquerading as an integer; reseal the payload digest and require semantic verification to reject every mutation.

- [ ] **Step 4: Run tests**

Run: `python -m pytest -q tests/test_strengthened_control.py`

Expected: all Task 1 and Task 3 tests pass.

- [ ] **Step 5: Commit Task 3**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/strengthened_control.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_strengthened_control.py
git commit -m "feat(issue128): derive strengthened adjacent-step gate"
```

### Task 4: Build the Independent Reference Verifier

**Files:**
- Create: `scripts/reference_verify_strengthened_control.py`
- Create: `tests/test_reference_strengthened_control.py`

**Interfaces:**
- CLI: `python scripts/reference_verify_strengthened_control.py SUMMARY WITNESS`.
- Exit 0 and print canonical `{"valid":true,...}` only after full semantic reconstruction.

- [ ] **Step 1: Write subprocess mutation tests**

Create temporary summary/witness pairs and assert rejection after changing one pair endpoint, adding a duplicate Pauli, omitting a singleton, inserting a commuting pair, changing a stage interval, changing the root enclosure, or changing `r` while recomputing both SHA fields.

- [ ] **Step 2: Implement strict JSON and deterministic gzip readers**

Reject duplicate keys, non-ASCII canonical JSON drift, nonzero gzip timestamp, extra members, path traversal, floats, booleans in integer positions, and noncanonical rationals.

- [ ] **Step 3: Reconstruct the theorem sum independently**

The script implements the four checkerboard Heisenberg bond fragments, local-coordinate canonicalization, exact symplectic commutators, coefficient-interval reconstruction, symplectic anticommutation, pair/single coverage, exact rational pair bounds, theorem scalar summation, factorial division, site/global density, adjacent-step inequalities, `30*r+1`, payload digests, and source hashes.  It must not import any `trottercert` module.  It independently implements the exact 24-decimal outward cube-root rounding policy, reconstructs the five-copy fourth-order Suzuki formula and 31-stage merge, and requires the canonical root interval and every stage record to equal the witness exactly; root containment alone is insufficient.  For every theorem word it recomputes the full canonical Pauli coefficient map and requires equality with the witness; submitted coefficient endpoints are never trusted.

Independently enumerate the complete raw theorem multiset `(fragment-index word, weak composition, outer fragment word, j, raw scalar)` from the frozen 31 stages.  Require exact multiset equality, including multiplicities, with `raw_terms`; then independently sum the raw scalars per `(outer,j)`, reconstruct every canonical Pauli coefficient map from the frozen fragments and the corresponding outer/base partial sum, and require exact equality with each aggregated record.  A submitted `(outer,j)` key set alone is never sufficient.

- [ ] **Step 4: Run reference tests**

Run: `python -m pytest -q tests/test_reference_strengthened_control.py`

Expected: all valid fixtures accepted and all mutations rejected.

- [ ] **Step 5: Commit Task 4**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/reference_verify_strengthened_control.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_reference_strengthened_control.py
git commit -m "feat(issue128): add independent strengthened control verifier"
```

### Task 5: Freeze the Scientific Artifacts

**Files:**
- Create: `scripts/certify_strengthened_control.py`
- Create: `benchmarks/paper-a/strengthened-control.json`
- Create: `benchmarks/paper-a/strengthened-control-witness.json.gz`
- Modify: `tests/test_strengthened_control.py`

**Interfaces:**
- CLI: `PYTHONPATH=src python scripts/certify_strengthened_control.py --source-commit COMMIT --summary PATH --witness PATH`.
- Output is deterministic and contains no execution timings.

- [ ] **Step 1: Implement deterministic artifact writing**

Serialize the witness with sorted compact JSON plus newline and gzip `mtime=0` and an empty member filename.  The source closure explicitly roots at `scripts/certify_strengthened_control.py`, `scripts/reference_verify_strengthened_control.py`, `src/trottercert/strengthened_control.py`, `src/trottercert/rigorous_fourth.py`, `src/trottercert/resources.py`, `src/trottercert/higher_order.py`, `src/trottercert/local_commutators.py`, `src/trottercert/intervals.py`, and `src/trottercert/algebra.py`; recursively AST-walk every local import and reject any unlisted transitive dependency.  Hash the resulting canonical relative-path manifest, `pyproject.toml`, and `requirements-reproducibility.txt`.  Replace NumPy floating-score ordering in witness discovery with a pure-Python exact rational ordering key, so witness bytes do not depend on NumPy or platform sorting.  Include an explicit `source_commit`, and refuse it unless every closure path in the worktree byte-matches `git show SOURCE_COMMIT:path`.

- [ ] **Step 2: Commit the artifact builder before freezing provenance**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/certify_strengthened_control.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_strengthened_control.py
git commit -m "feat(issue128): add strengthened control artifact builder"
```

- [ ] **Step 3: Capture the source commit and run the certified build after setup confirmation**

Run:

```bash
SOURCE_COMMIT="$(git rev-parse HEAD)"
PYTHONPATH=src python scripts/certify_strengthened_control.py \
  --source-commit "$SOURCE_COMMIT" \
  --summary benchmarks/paper-a/strengthened-control.json \
  --witness benchmarks/paper-a/strengthened-control-witness.json.gz
python scripts/reference_verify_strengthened_control.py \
  benchmarks/paper-a/strengthened-control.json \
  benchmarks/paper-a/strengthened-control-witness.json.gz
```

Expected: reference verifier exits 0.  Read the derived step from the output; do not treat 353 as the expected command result.

- [ ] **Step 4: Add frozen byte-rebuild tests**

Create a clean detached worktree at the artifact's `source_commit`, regenerate both files with the same explicit commit, and assert byte equality and SHA-256 equality with the candidate artifacts.  Clean-clone provenance replay may read sources through `git show source_commit:path`.  The release/DOI verifier instead hashes the bundled canonical relative paths against the embedded closure manifest and does not require a `.git` directory.

- [ ] **Step 5: Run focused and regression tests**

Run:

```bash
python -m pytest -q tests/test_strengthened_control.py \
  tests/test_reference_strengthened_control.py tests/test_rigorous_fourth.py
```

- [ ] **Step 6: Independently review the result state**

The reviewer reports one of: derived 353/10591 and adjacent gate passes; derived a different exact step/resource count; or verifier inconclusive.  Review fails if the artifact simply echoes the hypothesis.

- [ ] **Step 7: Commit only the frozen artifacts**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/strengthened-control.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/strengthened-control-witness.json.gz
git commit -m "certify(issue128): freeze strengthened rational-pair control"
```

### Task 6: Promote or Refute the 3.71-Fold Comparison

**Files:**
- Modify: `scripts/build_paper_a_ablation.py`
- Modify: `tests/test_paper_a_ablation.py`
- Modify: `benchmarks/paper-a/ablation.json`
- Modify: `artifacts/publication/paper-a-claim-matrix.json`
- Modify: `artifacts/publication/paper-a-file-ownership.json`
- Modify: `artifacts/publication/paper-a-claim-SHA256SUMS`
- Modify: `scripts/audit_paper_a_claims.py`
- Modify: `scripts/audit_publication_scope.py`
- Modify: `docs/manuscript/sections/results.tex`
- Modify: `docs/manuscript/sections/limitations.tex`

**Interfaces:**
- Consumes only the independently accepted summary and witness hashes.
- Produces either an included exact strengthened-control claim or an explicit refutation/inconclusive record.

- [ ] **Step 1: Add a separately named ablation row**

Insert `rational_pair_theorem_control` after `published`.  Preserve `complete_formula` as `missing_evidence`; its semantics do not change.  Rebuild the row from the summary artifact, not from copied constants.

- [ ] **Step 2: Gate the claim matrix**

Extend the strict claim schema with one claim name only: `strengthened_control_result`.  Its typed fields include `status`, optional `accepted_step`, optional `resource_count`, optional `ratio`, both artifact hashes, and the independent verifier command.  Add that single name to required/source-checker allowlists, claim bindings, manuscript bindings, and the hash manifest; `strengthened_control_ratio` is not a separate claim name.  If and only if the derived state is `certified_matches_353`, the record contains `ratio=[10591,2851]`.  If the result differs, record the actual exact step/resource result and retain the 10591 token ban.  An inconclusive record has no ratio and cannot promote.  Refresh `paper-a-claim-SHA256SUMS` from the accepted canonical bytes.

- [ ] **Step 3: Update ownership and release policy**

Classify the summary, witness, generator, verifier, tests, checksum manifest, and the two manuscript sections as Paper A.  Include only the minimum reproduction closure in the release and verify transitive dependencies.  Add a manuscript result marker sourced from the claim matrix and a limitations sentence that distinguishes `rational_pair_theorem_control` from `complete_formula`.

- [ ] **Step 4: Run the complete Paper A claim gate**

Run:

```bash
python -m pytest -q tests/test_paper_a_ablation.py tests/test_paper_a_claims.py \
  tests/test_publication_scope.py
PYTHONPATH=src python scripts/audit_paper_a_claims.py \
  artifacts/publication/paper-a-claim-matrix.json --root .
PYTHONPATH=src python scripts/audit_publication_scope.py \
  artifacts/publication/paper-a-file-ownership.json --root .
```

Expected: all strengthened-control semantics pass.  Unrelated concurrent unclassified paths must be reported separately rather than added to this commit.

- [ ] **Step 5: Commit Task 6**

Stage only the explicitly listed integration files:

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_paper_a_ablation.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_paper_a_ablation.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/benchmarks/paper-a/ablation.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-claim-matrix.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-file-ownership.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-claim-SHA256SUMS \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/audit_paper_a_claims.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/audit_publication_scope.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/manuscript/sections/results.tex \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/manuscript/sections/limitations.tex
git commit -m "docs(issue128): gate strengthened control comparison"
```


## Self-Review Checklist

- [ ] The plan never assumes that 353 is true.
- [ ] The rational-pair method is never mislabeled as complete formal-log evaluation.
- [ ] The witness contains enough information for a verifier that does not import the generator.
- [ ] Both adjacent-step inequalities are exact.
- [ ] `complete_formula` remains missing evidence.
- [ ] A result different from 353 is a successful scientific outcome, not a test failure.
- [ ] The 95-step primary certificate is never modified.
