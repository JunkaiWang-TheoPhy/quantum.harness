# Direct Finite-Torus Grouped XXZ Compiler Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce the first positive, independently verified nonisotropic XXZ transfer certificate without reusing the isotropic D4/D5 sidecars.

**Architecture:** Build an exact weighted-Pauli compiler on a direct 4 by 4 periodic torus, where wraparound and aliasing are explicit.  Freeze the compiler after pilot validation at `Delta=1/2` and `Delta=2`, then evaluate the preregistered seven-point family and, only after the finite-torus gate passes, scale the same rules to the 12 by 12 publication geometry.

**Tech Stack:** Python 3.12, `Fraction`, exact Pauli/symplectic dictionaries, rational intervals, deterministic JSON/gzip, SHA-256, pytest, Slurm for projected runs above the local threshold.

## Global Constraints

- Hamiltonian: `H_Delta=sum_<uv>(X_uX_v+Y_uY_v+Delta*Z_uZ_v)/4`.
- Pilot lattice: 4 by 4 periodic square torus, 16 spins, 32 unique nearest-neighbor bonds.
- Four checkerboard fragments: horizontal even/odd and vertical even/odd; eight disjoint bonds per fragment.
- Formula: five-copy fourth-order Suzuki with 31 merged stages.
- Target: `T=1`, operator-norm tolerance `epsilon=10^-6`.
- Pilot anisotropies: exact `Delta=1/2` and `Delta=2`; full preregistered grid is `0,1/4,1/2,1,3/2,2,4`.
- `Delta=1` is an isotropic identity anchor and does not count as transfer.
- Primary resource metric: merged schedule-group exponentials `G(r)=30*r+1`.
- Do not claim a CNOT count until a separate general XXZ bond-synthesis certificate exists.
- Direct finite-torus evaluation must handle wraparound explicitly; unit-cell-density multiplication is forbidden on the 4 by 4 pilot.
- Each Delta gets its own term universe, grouping witness, baseline, finite-step bound, and adjacent-step check.
- Full scientific runs projected above 10 minutes or 16 GB run through the configured Slurm profile.
- Existing isotropic source, sidecars, and artifact remain immutable.

---

## File Structure

### New files

- `src/trottercert/grouped_xxz.py`: strict compile spec, direct finite-torus fragments, exact weighted Pauli backend, theorem ledger, grouping records, and semantic verifier.
- `scripts/compile_grouped_xxz.py`: deterministic per-Delta builder and bounded profiling mode.
- `scripts/reference_verify_grouped_xxz.py`: independent certificate verifier with no generator imports.
- `scripts/certify_grouped_xxz_scaling.py`: deterministic builder/primary verifier for the optional 12 by 12 geometry gate.
- `tests/test_grouped_xxz.py`: Hamiltonian, torus, schedule, ledger, grouping, finite-step, and schema tests.
- `tests/test_reference_grouped_xxz.py`: subprocess mutation suite.
- `tests/test_grouped_xxz_scaling.py`: support injectivity, orbit/stabilizer, alias, source-binding, and scaling-certificate tests.
- `benchmarks/paper-a/grouped-xxz-pilot.json`: pilot summary for `Delta=1/2,2`.
- `benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz`: exact pilot ledgers and grouping witnesses.
- `benchmarks/paper-a/grouped-xxz-family.json`: final seven-point result, including honest failures.
- `benchmarks/paper-a/grouped-xxz-family-witness.json.gz`: complete per-Delta ledgers, theorem-block groups, local-log records, and row digests for all seven points.
- `benchmarks/paper-a/grouped-xxz-12x12.json`: optional 12 by 12 result; absent or explicitly `future_non_claim` until its independent scaling certificate closes.
- `benchmarks/paper-a/grouped-xxz-12x12-witness.json.gz`: geometry-specific 12 by 12 ledger, growth/tail, baseline, adjacent-step, and grouping evidence when promoted.
- `artifacts/verification/grouped-xxz-scaling-12x12.json`: canonical geometry authorization, created only after every support/replay gate passes.

### Existing files to modify after pilot review

- `src/trottercert/hamiltonian.py`: expose exact XXZ bond and four-matching fragment constructors.
- `scripts/run_xxz_transfer.py`: consume the new compiler only after its spec commit is frozen.
- `tests/test_xxz_transfer.py`: replace the hardcoded nonisotropic prohibition with artifact-backed result validation.
- `benchmarks/paper-a/xxz-transfer.json`: regenerate the preregistered family.
- `artifacts/publication/paper-a-claim-matrix.json`: promote transfer only after a nonisotropic row passes.
- `artifacts/publication/paper-a-file-ownership.json`: classify new Paper A paths.
- `scripts/audit_publication_scope.py`: include the reviewed release closure.

---

### Task 1: Freeze the Exact XXZ Compile Spec

**Files:**
- Create: `src/trottercert/grouped_xxz.py`
- Create: `tests/test_grouped_xxz.py`

**Interfaces:**
- Produces: `XXZCompileSpec` dataclass.
- Produces: `XXZStatus = Literal["certified", "unsupported", "inconclusive"]`.
- Produces strict fraction-pair and canonical JSON helpers scoped to this artifact.

- [ ] **Step 1: Write exact setup tests**

```python
def test_pilot_spec_pins_the_physics() -> None:
    spec = XXZCompileSpec.pilot(Fraction(1, 2))
    assert spec.length == 4
    assert spec.boundary == "periodic"
    assert spec.time == Fraction(1)
    assert spec.tolerance == Fraction(1, 10**6)
    assert spec.normalization == "(XX+YY+delta*ZZ)/4"
    assert spec.primary_metric == "merged_group_exponentials"


@pytest.mark.parametrize("delta", [True, 0.5, "1/2"])
def test_spec_requires_exact_delta(delta: object) -> None:
    with pytest.raises(TypeError):
        XXZCompileSpec.pilot(delta)  # type: ignore[arg-type]
```

- [ ] **Step 2: Run the missing-module failure**

Run: `python -m pytest -q tests/test_grouped_xxz.py`

- [ ] **Step 3: Implement strict spec construction**

`XXZCompileSpec.pilot(delta)` requires `Delta` as `Fraction`, fixes `length=4`, literal periodic boundary, `T=1`, `epsilon=10^-6`, formula identifier, stage count 31, and primary metric.  A separate `XXZCompileSpec.periodic(delta, length, scaling_certificate)` accepts an even `length>=4` only after validating the geometry-specific scaling certificate; production callers cannot silently switch geometry.  Reject all extra fields, floats, booleans, noncanonical rational pairs, odd lengths, and lengths below four.

- [ ] **Step 4: Commit Task 1**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/grouped_xxz.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_grouped_xxz.py
git commit -m "feat(issue128): freeze grouped XXZ compile spec"
```

### Task 2: Implement Exact Bond and Four-Matching Fragments

**Files:**
- Modify: `src/trottercert/hamiltonian.py`
- Modify: `src/trottercert/grouped_xxz.py`
- Modify: `tests/test_grouped_xxz.py`
- Modify: `tests/test_lattice_hamiltonian.py`

**Interfaces:**
- Produces: `xxz_bond_operator(left, right, delta) -> PauliSum`.
- Produces: `four_matching_xxz_fragments(lattice, delta) -> tuple[PauliSum, PauliSum, PauliSum, PauliSum]`.
- Produces: `finite_torus_xxz_hamiltonian(lattice, delta) -> PauliSum`.

- [ ] **Step 1: Write bond normalization tests**

```python
def test_xxz_bond_operator_has_exact_axis_weights() -> None:
    bond = xxz_bond_operator(0, 5, Fraction(3, 2))
    assert bond.terms[PauliString({0: "X", 5: "X"})].real == Fraction(1, 4)
    assert bond.terms[PauliString({0: "Y", 5: "Y"})].real == Fraction(1, 4)
    assert bond.terms[PauliString({0: "Z", 5: "Z"})].real == Fraction(3, 8)
```

- [ ] **Step 2: Write 4 by 4 torus coverage tests**

Assert exactly four matching bond sets, eight unique bonds in each, pairwise disjoint bond sets, 32 total bonds, no self-loops, and equality between the sum of fragments and the independently enumerated periodic Hamiltonian for `Delta=0,1/2,1,2`.

- [ ] **Step 3: Implement constructors**

Use `SquareLattice.four_matchings()` as a source of bond coordinates but independently validate the cardinality and uniqueness invariants.  Build every term with exact `Fraction`; never route through the isotropic `heisenberg_bond` helper.

- [ ] **Step 4: Add isotropic bridge test**

At `Delta=1`, require exact equality of the normalized Pauli coefficient dictionaries with `four_matching_fragments(SquareLattice(4))`.

- [ ] **Step 5: Run focused tests**

Run: `python -m pytest -q tests/test_grouped_xxz.py tests/test_lattice_hamiltonian.py`

- [ ] **Step 6: Commit Task 2**

Stage only the four listed files and commit with:

`git commit -m "feat(issue128): construct exact finite-torus XXZ fragments"`

### Task 3: Freeze and Replay the 31-Stage Schedule

**Files:**
- Modify: `src/trottercert/grouped_xxz.py`
- Modify: `tests/test_grouped_xxz.py`
- Create: `scripts/compile_grouped_xxz.py` with a profiling-only skeleton; artifact build subcommands remain fail-closed until Task 8.

**Interfaces:**
- Produces: `xxz_suzuki_schedule() -> tuple[StageRecord, ...]`.
- Produces: `replay_schedule_resources(schedule, steps) -> int`.

- [ ] **Step 1: Write schedule tests**

Assert 31 stages, fragment indices in `0..3`, palindromic stage sequence, exact cubic-field coefficient coordinates, exact total coefficient one for each fragment, and `replay_schedule_resources(schedule,r)==30*r+1` for `r=1,2,95`.

- [ ] **Step 2: Implement explicit schedule records**

Convert `fourth_order_suzuki_cubic_stages(4)` into canonical records storing `fragment_index` and `[a0,a1,a2]` rational coordinates.  Verification reconstructs the sums rather than trusting `stage_count`.

- [ ] **Step 3: Add stage-order mutations**

Swap unequal adjacent stages, flip one cubic coordinate, duplicate a stage, or replace a fragment index; reseal the digest and require rejection.

- [ ] **Step 4: Run and commit**

Run: `python -m pytest -q tests/test_grouped_xxz.py`

Commit: `git commit -m "feat(issue128): bind grouped XXZ Suzuki schedule"`

### Task 4: Build the Direct Finite-Torus Theorem Ledger

**Files:**
- Modify: `src/trottercert/grouped_xxz.py`
- Modify: `tests/test_grouped_xxz.py`

**Interfaces:**
- Produces: `build_finite_xxz_ledger(spec, *, max_degree=5, progress=None) -> XXZLedger`.
- `XXZLedger` contains every raw theorem word, exact positive scalar interval, source fragment identity, immutable theorem norm-block key, and an unweighted canonical Pauli coefficient map for that key.

- [ ] **Step 1: Write reduced-instance algebra tests**

On a two-bond open fixture, compare every degree-three and degree-five nested commutator from the weighted symplectic backend with `PauliSum` commutator evaluation for `Delta=1/2` and `Delta=2`.

- [ ] **Step 2: Implement common-denominator weighted symplectic terms**

For `Delta=a/b`, use denominator `4b` and integer bond weights `(b,b,a)`.  Implement exact symplectic multiplication phases and canonical finite-site masks; no translation quotient is used in the 4 by 4 pilot.

- [ ] **Step 3: Enumerate theorem words from the frozen schedule**

Reuse only the general weak-composition theorem indexing from `fourth_order_published_triangle_certificate`, with theorem identity `published_high_order_triangle_v1`, order four, one-based center 20, factorial denominator `5!`, and the same left/right partial-sum convention.  A norm-block key `k` identifies one raw nested-commutator polynomial `C_k`; its canonical Pauli map never contains the theorem scalar.  Every raw weak-composition record stores its key and positive scalar upper weight, and `w_k` is defined once as the exact sum of all raw weights carrying key `k`, including multiplicity.  Include `print(..., flush=True)` progress every approximately five percent of records.

- [ ] **Step 4: Add full-coverage invariants**

Verify unique raw theorem identities, exact multiplicities, no missing fragment word or norm block, `w_k=sum(raw weights for k)`, sorted Pauli masks, nonzero coefficients only, and equality of each unweighted `C_k` Pauli polynomial with a separately accumulated reference on reduced fixtures.  Tests reject both a preweighted `C_k` map and a second multiplication of a submitted weight.  Never aggregate coefficient maps across different theorem norm blocks.

- [ ] **Step 5: Run a bounded profiling mode after setup confirmation**

Run:

```bash
PYTHONPATH=src python scripts/compile_grouped_xxz.py \
  --delta 1/2 --profile-only --pipeline direct-theorem --max-records 200
```

The command reports records/second, projected wall time, and projected peak storage without creating a scientific artifact.  If the full pilot projects above 10 minutes or 16 GB, package the same immutable commit for Slurm.

- [ ] **Step 6: Commit Task 4**

Commit: `git commit -m "feat(issue128): enumerate direct finite XXZ ledger"`

### Task 5: Certify Pairwise-Anticommuting Groups

**Files:**
- Modify: `src/trottercert/grouped_xxz.py`
- Modify: `tests/test_grouped_xxz.py`

**Interfaces:**
- Produces: `discover_xxz_groups(ledger) -> tuple[TheoremBlockGroupRecord, ...]`.
- Produces: `verify_xxz_groups(ledger, groups) -> Fraction`, where the return value is the positive weighted sum of independently bounded theorem blocks.

- [ ] **Step 1: Write exact group-verification tests**

Use a fixture containing `XI`, `ZI`, `IX`, and `IZ`.  Require exact coverage, reject duplicates and omissions, reject a commuting pair in one group, and reproduce the outward rational square-root upper bound.

- [ ] **Step 2: Implement discovery/trust separation**

Discovery may use deterministic greedy ordering by exact coefficient magnitude, but only within one immutable theorem norm block.  Verification ignores discovery scores, reconstructs all symplectic commutation relations, requires each Pauli term exactly once inside its own block, and recomputes every group norm with outward rational intervals.  Cross-block grouping or cancellation is a hard schema error.

- [ ] **Step 3: Add Delta-specific identity binding**

The group witness stores the exact Delta and ledger digest.  Reusing a `Delta=1/2` witness with a `Delta=2` ledger must fail even if Pauli masks coincide.

- [ ] **Step 4: Reconstruct the theorem sum without cross-block cancellation**

For each unweighted block polynomial `C_k`, aggregate identical Pauli masks and compute its certified grouped norm `U_k`.  Independently recompute `w_k` from the raw record multiset; neither `C_k` nor `U_k` contains `w_k`.  The final theorem constant is exactly `K = sum_k w_k*U_k/5!`, with each theorem weight `w_k>=0`.  Freeze the theorem identifier, one-based center 20, the Duhamel/factorial factor, and the finite-step error formula `E_r = K/r^4` for `T=1`.  The same-Delta baseline uses `sum_k w_k*||C_k||_Pauli-l1/5!`; it does not combine the `C_k` maps.

- [ ] **Step 5: Run and commit**

Run: `python -m pytest -q tests/test_grouped_xxz.py`

Commit: `git commit -m "feat(issue128): certify grouped XXZ norm witnesses"`

### Task 5A: Compress the Production Ledger by Exact D4 Orbits

**Files:**
- Create: `src/trottercert/grouped_xxz_compressed.py`
- Create: `tests/test_grouped_xxz_compressed.py`
- Modify: `scripts/compile_grouped_xxz.py`

**Interfaces:**
- Consumes: `_iter_raw_theorem_records()`, `_FiniteXXZSymplecticEvaluator`,
  `discover_anticommuting_groups()`, and the frozen actions in
  `grouped_xxz_orbits.py`.
- Produces: `build_compressed_xxz_ledger(spec, *, max_records=None, progress=None)`.
- Produces: `verify_compressed_xxz_ledger(ledger)` and
  `compile_compressed_grouped_xxz(spec, ledger=None)`.
- A bounded ledger has `complete=False` and cannot be passed to the compiler.

- [ ] **Step 1: Write failing raw-stream and orbit tests**

Require the first 200 raw records to have the same identities, exact weights,
actual-word coverage, and canonical stream digest as the uncompressed ledger.
For every active word, require one frozen symmetry that maps its canonical
representative to it.  Mutating a word weight, representative, symmetry id, or
stream digest must fail semantic replay after resealing outer digests.

- [ ] **Step 2: Implement the streaming theorem ledger**

Serialize each raw record with `canonical_json_bytes` as it is enumerated,
update one SHA-256 stream digest, and accumulate exact positive weights in a
dictionary keyed by the actual five-letter word.  Store only the count, digest,
and sorted per-word weights.  Do not retain the 61,677 raw dataclasses and do
not aggregate weights into representative records.

- [ ] **Step 3: Evaluate and verify canonical representative blocks**

For each active canonical representative, use the exact finite-torus evaluator
once, freeze its sorted `SymplecticCoefficient` map, and build
`deterministic_pair_only_bitset_v1` groups.  Verify the eight site/fragment
actions structurally and require exact transported-map equality on bounded
integration fixtures.  Store actual-word-to-representative symmetry metadata,
but never duplicate transported term maps.

- [ ] **Step 4: Reconstruct both theorem constants literally**

Loop over sorted actual words and add `w_k*U_rep(k)` or
`w_k*||C_rep(k)||_1` one word at a time before dividing by `5!`.  The verifier
independently repeats this loop.  It rejects a submitted representative-weight
sum even if the final scalar is numerically equal.

- [ ] **Step 5: Close full certificates and fail closed on prefixes**

For a complete ledger, compute the least positive candidate and baseline steps
using exact integer fourth-power comparisons, store both accepted errors and
strict predecessor errors, and replay `G(r)=30*r+1`.  Calling the certificate
compiler with `complete=False`, a missing word, or a partial orbit is an error.

- [ ] **Step 6: Run the 200-record end-to-end profile for both anisotropies**

Run:

```bash
python -m pytest -q -m slow tests/test_grouped_xxz_compressed.py \
  -k compressed_200_record_profile
```

Record raw count, active words, representative count, representative Pauli
terms, pair/singleton counts, grouped/triangle ratio, wall time, and peak RSS.
Authorize a full local run only if the measured projection is below 10 minutes
and 16 GB; otherwise package the same immutable commit for Slurm.

- [ ] **Step 7: Commit the compressed core**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/grouped_xxz_compressed.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_grouped_xxz_compressed.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/compile_grouped_xxz.py
git commit -m "perf(issue128): compress grouped XXZ theorem orbits"
```

- [ ] **Step 8: Conditional stronger-group fallback**

Only if the complete pair-only certificate fails the positive-transfer resource
gate, run the existing multi-member greedy discovery on canonical
representatives, bind a distinct algorithm identifier, and rerun the same
primary/reference mutation gates.  Do not change theorem weights, orbit rules,
baseline, or the preregistered Delta values after seeing the pair-only result.

### Task 6: Close the Finite-Step Bound and Same-Delta Baseline

**Files:**
- Modify: `src/trottercert/grouped_xxz.py`
- Modify: `tests/test_grouped_xxz.py`

**Interfaces:**
- Produces: `compile_grouped_xxz(spec) -> XXZCertificate`.
- Produces: exact candidate and published-theorem baseline steps for the same `Delta` and torus.

- [ ] **Step 1: Write exact adjacent-step tests**

For a small fake ledger, require candidate `E_r<=epsilon`, `E_(r-1)>epsilon`, baseline `B_s<=epsilon`, `B_(s-1)>epsilon`, and resources `30*r+1` and `30*s+1`.

- [ ] **Step 2: Implement the direct high-order theorem bound**

Use `K_grouped=sum_k w_k U_k/5!` as the candidate constant, with aggregation and anticommuting groups confined to each theorem block `C_k`.  Independently calculate `K_triangle=sum_k w_k ||C_k||_Pauli-l1/5!` as the same-Delta baseline.  For `T=1`, certify `E_r=K_grouped/r^4` and `B_s=K_triangle/s^4` by exact `Fraction`; do not reuse the isotropic 393-step value.  Store theorem identity, center 20, order four, `5!`, and the Duhamel convention in both certificate and witness.

- [ ] **Step 3: Add model-dependent growth metadata**

Record, but do not trust without recomputation, `bond_growth=max(1,(1+abs(Delta))/2)` and `cell_base=(2+abs(Delta))/2`.  Any future D6/D7/tail extension consumes these exact values rather than the isotropic constants.

- [ ] **Step 4: Keep claim scope explicit**

Mark the pilot method `direct_finite_high_order_theorem_grouped_norm`.  It is a rigorous finite-step transfer certificate, but it is not the flagship local-log D4-D7-tail pipeline.  The manuscript must compare methods by name.

- [ ] **Step 5: Run and commit**

Run: `python -m pytest -q tests/test_grouped_xxz.py`

Commit: `git commit -m "feat(issue128): close grouped XXZ finite-step bounds"`

### Task 7: Add an Independent Reference Verifier

**Files:**
- Create: `scripts/reference_verify_grouped_xxz.py`
- Create: `tests/test_reference_grouped_xxz.py`

**Interfaces:**
- CLI: `python scripts/reference_verify_grouped_xxz.py SUMMARY WITNESS`.
- The verifier uses only the standard library and explicit artifact data.

- [ ] **Step 1: Write end-to-end attacks**

Mutate Delta, axis weight, length, boundary, one periodic bond, stage order, theorem word, coefficient interval, Pauli mask, group coverage, anticommutation, candidate step, previous-step bound, baseline, resource count, source hash, and payload digest.  Recompute outer digests and require every semantic mutation to fail.

- [ ] **Step 2: Implement independent replay**

Reconstruct 16-site masks, all 32 bonds, four matching coverage, schedule identities, theorem scalar sums, ledger coverage, symplectic anticommutation, grouped and triangle norms, exact step inequalities, and resources.  Do not import `trottercert` or execute the generator.

- [ ] **Step 3: Run and commit**

Run: `python -m pytest -q tests/test_reference_grouped_xxz.py`

Commit: `git commit -m "feat(issue128): independently verify grouped XXZ"`

### Task 8: Run and Freeze the Two-Point Pilot

**Files:**
- Modify: `scripts/compile_grouped_xxz.py`
- Create: `benchmarks/paper-a/grouped-xxz-pilot.json`
- Create: `benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz`
- Modify: `tests/test_grouped_xxz.py`

**Interfaces:**
- Per-Delta CLI has mutually exclusive `--build` and `--profile-only` modes; accepts exact `--delta NUM/DEN`, `--pipeline direct-theorem|local-log`, explicit `--summary`/`--witness` output paths in build mode, and optional positive integer `--max-records` only in profile mode.  There is no implicit pipeline or output path in production mode.
- Aggregator CLI `--merge-pilot LEFT_SUMMARY LEFT_WITNESS RIGHT_SUMMARY RIGHT_WITNESS --summary OUT --witness OUT` verifies both inputs and builds one canonical pilot bundle.

- [ ] **Step 1: Freeze the source commit and compute-feasibility record**

Record local profile or Slurm job parameters before the full run.  Scientific artifacts contain no timing-dependent fields; timing and peak memory go to a separate transcript.

- [ ] **Step 2: Run `Delta=1/2` and `Delta=2` independently**

Use separate output roots and atomic writes.  Each process emits progress every approximately five percent and self-verifies before publishing.  Then merge them deterministically:

```bash
PYTHONPATH=src python scripts/compile_grouped_xxz.py --build --delta 1/2 --pipeline direct-theorem \
  --summary build/grouped-xxz/delta-1-2.json \
  --witness build/grouped-xxz/delta-1-2-witness.json.gz
PYTHONPATH=src python scripts/compile_grouped_xxz.py --build --delta 2/1 --pipeline direct-theorem \
  --summary build/grouped-xxz/delta-2-1.json \
  --witness build/grouped-xxz/delta-2-1-witness.json.gz
PYTHONPATH=src python scripts/compile_grouped_xxz.py --merge-pilot \
  build/grouped-xxz/delta-1-2.json build/grouped-xxz/delta-1-2-witness.json.gz \
  build/grouped-xxz/delta-2-1.json build/grouped-xxz/delta-2-1-witness.json.gz \
  --summary benchmarks/paper-a/grouped-xxz-pilot.json \
  --witness benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz
```

Every gzip member uses canonical compact sorted JSON plus newline, `mtime=0`, and an empty filename.  The merged witness contains the two complete per-Delta ledgers, not only their summaries.

- [ ] **Step 3: Run the independent verifier**

Both per-Delta results and the merged bundle must pass the reference verifier:

```bash
python scripts/reference_verify_grouped_xxz.py build/grouped-xxz/delta-1-2.json build/grouped-xxz/delta-1-2-witness.json.gz
python scripts/reference_verify_grouped_xxz.py build/grouped-xxz/delta-2-1.json build/grouped-xxz/delta-2-1-witness.json.gz
python scripts/reference_verify_grouped_xxz.py benchmarks/paper-a/grouped-xxz-pilot.json benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz
```

A certified easy-plane and easy-axis pair passes the positive pilot gate; a valid unsupported/inconclusive result remains publishable but does not pass it.

- [ ] **Step 4: Rebuild byte-identically**

Regenerate from the same source commit in a fresh temporary directory and compare summary and deterministic gzip bytes:

```bash
tmpdir="$(mktemp -d)"
PYTHONPATH=src python scripts/compile_grouped_xxz.py --build --delta 1/2 \
  --pipeline direct-theorem --summary "$tmpdir/delta-1-2.json" \
  --witness "$tmpdir/delta-1-2-witness.json.gz"
PYTHONPATH=src python scripts/compile_grouped_xxz.py --build --delta 2/1 \
  --pipeline direct-theorem --summary "$tmpdir/delta-2-1.json" \
  --witness "$tmpdir/delta-2-1-witness.json.gz"
cmp build/grouped-xxz/delta-1-2.json "$tmpdir/delta-1-2.json"
cmp build/grouped-xxz/delta-1-2-witness.json.gz "$tmpdir/delta-1-2-witness.json.gz"
cmp build/grouped-xxz/delta-2-1.json "$tmpdir/delta-2-1.json"
cmp build/grouped-xxz/delta-2-1-witness.json.gz "$tmpdir/delta-2-1-witness.json.gz"
shasum -a 256 "$tmpdir"/*
```

- [ ] **Step 5: Independent scientific review**

Review Hamiltonian signs, 32-bond torus coverage, alias handling, theorem scope, both baselines, adjacent steps, and resource semantics before commit.

- [ ] **Step 6: Commit Task 8**

Commit: `git commit -m "feat(issue128): certify grouped XXZ pilot"`

### Task 9: Generalize the Flagship Local-Log Ledger to Weighted XXZ

**Files:**
- Modify: `src/trottercert/refined_error.py`
- Modify: `src/trottercert/grouped_xxz.py`
- Modify: `scripts/reference_verify_grouped_xxz.py`
- Modify: `tests/test_grouped_xxz.py`
- Modify: `tests/test_refined_error.py`
- Modify: `tests/test_reference_grouped_xxz.py`

**Interfaces:**
- Produces: `HamiltonianGrowthProfile` with exact bond-axis l1, finite-torus Hamiltonian l1, all-order support-growth recurrence, convergence condition, and tail-base rationals.
- Produces: `compile_local_log_grouped_xxz(spec) -> XXZCertificate` with exact D4-D7 contributions and generator tail.

- [ ] **Step 1: Expose every isotropic hardcoded growth constant as exact data**

Write a regression fixture whose `Delta=1` profile reproduces the current bond growth, cell base, D6, D7, and tail values exactly.  The default production verifier rejects an omitted profile rather than silently selecting isotropic constants.

- [ ] **Step 2: Derive the finite-torus weighted profile**

For the 4 by 4 Hamiltonian, reconstruct the exact Pauli-l1 norm directly from all 32 bonds and the axis weights.  Compute low-degree support growth by exact bond adjacency on the torus, but do not extrapolate from the saturated 16-site support.

Freeze one finite proof-object grammar, `square_nn_linear_incidence_v1`; arbitrary submitted expressions are forbidden.  Its only exact rational/integer fields are `bond_pauli_l1`, `coordination=4`, `commutator_factor=2`, `max_new_sites_per_adjoint=1`, `initial_support`, `initial_l1`, `first_omitted_degree=8`, and the absolute stage-prefix list.  The verifier derives, rather than trusts,

```text
s_q = initial_support + q,
B_0 = initial_l1,
B_(q+1) = commutator_factor * bond_pauli_l1 * coordination * s_q * B_q.
```

Before using the recurrence, the reference verifier independently recomputes `bond_pauli_l1` from the exact Delta-weighted bond Hamiltonian and the absolute stage-prefix list from the frozen 31-stage schedule.  It reconstructs the complete first-omitted degree-eight seed ledger and requires submitted `initial_support` and `initial_l1` to be exact certified upper bounds for every seed term; neither base value is trusted.  Mutation tests lower each primitive/base/prefix field, reseal all hashes, and require rejection.  Only then does the verifier prove the overlap count by the square-lattice incidence bound `m(S)<=4|S|`, check each explicit D4-D7 support against `s_q`, divide `B_q` by the exact exponential Taylor and Duhamel factorials, and derive the tail terms for every stage prefix.  The certificate stores the first omitted term and a rational ratio majorant `rho(r)`; the verifier algebraically proves `a_(q+1)/a_q<=rho(r)<1` for every `q>=8` from this fixed recurrence, then recomputes the geometric sum `a_8/(1-rho)`.  Any unsupported recurrence kind, fitted coefficient, free-form expression, missing base/induction field, or downward-mutated primitive is rejected.

- [ ] **Step 3: Evaluate exact degree-five and degree-seven logarithms**

Use the frozen cubic-field 31-stage series with the weighted finite-site symplectic evaluator and freeze the identities

```text
D4 = 5 E5
D5 = 2 ad_H(E5)
D6 = 7 E7 + (2/3) ad_H^2(E5)
D7 = 3 ad_H(E7) + (1/6) ad_H^3(E5).
```

Aggregate identical Pauli masks before norms, certify D4/D5 grouping witnesses, and retain D6/D7 exact Pauli-l1 contributions.  For a finite torus, the degree-`d` contribution is exactly `||D_d||/((d+1)r^d)` at `T=1` (or `T^(d+1)||D_d||/((d+1)r^d)` in the general schema).  For a translation-density certificate it is `N*density_d/((d+1)r^d)` at `T=1`.  Add the degree-eight-and-higher generator tail only after the fixed `square_nn_linear_incidence_v1` recurrence and its ratio `<1` are verified.  Emit flushed progress by degree and word index.

- [ ] **Step 4: Require two finite-step closures**

The local-log certificate and the direct-theorem certificate must each independently satisfy their own exact adjacent-step gate.  `local_log <= direct_theorem` at a common step is a performance/promotion gate only: if it fails, keep both sound rows but do not promote local-log as the improved method.

- [ ] **Step 5: Add weighted-growth mutations**

Change one axis weight, bond growth, support-neighborhood count, tail recurrence coefficient, convergence witness, D4-D7 mapping coefficient, D6/D7 contribution, or final step and reseal all outer hashes; require both primary and reference verification to reject.  The reference verifier independently rebuilds the weighted growth profile, D4-D7 map, degree-eight tail recurrence, and adjacent-step inequalities.

- [ ] **Step 6: Run focused tests and the bounded profile**

Run:

```bash
python -m pytest -q tests/test_grouped_xxz.py tests/test_refined_error.py
PYTHONPATH=src python scripts/compile_grouped_xxz.py \
  --delta 1/2 --profile-only --pipeline local-log --max-records 200
```

Move the full local-log run to Slurm if it crosses the declared local threshold.

- [ ] **Step 7: Commit Task 9**

Commit: `git commit -m "feat(issue128): generalize local-log bounds to XXZ"`

### Task 10: Freeze the Seven-Point Family and Scale the Geometry

**Files:**
- Modify: `scripts/run_xxz_transfer.py`
- Modify: `tests/test_xxz_transfer.py`
- Modify: `benchmarks/paper-a/xxz-transfer.json`
- Create: `benchmarks/paper-a/grouped-xxz-family.json`
- Create: `benchmarks/paper-a/grouped-xxz-family-witness.json.gz`
- Create: `scripts/certify_grouped_xxz_scaling.py`
- Create: `tests/test_grouped_xxz_scaling.py`
- Create only after the scale gate passes: `artifacts/verification/grouped-xxz-scaling-12x12.json`
- Create only after the scale gate passes: `benchmarks/paper-a/grouped-xxz-12x12.json`
- Create only after the scale gate passes: `benchmarks/paper-a/grouped-xxz-12x12-witness.json.gz`

**Interfaces:**
- Consumes the frozen direct-theorem compiler from Task 8 and the frozen local-log compiler plus updated reference verifier from Task 9.
- Produces one canonical row for every preregistered Delta and a complete witness section/digest per row, with no post-hoc changes.  Report `direct_theorem` and `local_log` as separate method columns; never merge their constants.

- [ ] **Step 1: Remove the hardcoded `Delta!=1` prohibition**

Replace it with artifact-backed calls to the frozen compiler.  The schema directly records length, boundary, T, epsilon, formula, method, baseline, candidate, and witness hashes.

- [ ] **Step 2: Run the seven-point 4 by 4 family**

Run all points regardless of pilot outcome and retain certified, unsupported, and inconclusive rows.  Each Delta has a complete deterministic ledger/group/local-log witness and its own digest inside `grouped-xxz-family-witness.json.gz`; the summary row binds that digest.  Build and verify with:

```bash
PYTHONPATH=src python scripts/run_xxz_transfer.py --build-grouped-family \
  --deltas 0/1,1/4,1/2,1/1,3/2,2/1,4/1 --length 4 \
  --summary benchmarks/paper-a/grouped-xxz-family.json \
  --witness benchmarks/paper-a/grouped-xxz-family-witness.json.gz \
  --transfer benchmarks/paper-a/xxz-transfer.json
python scripts/reference_verify_grouped_xxz.py \
  benchmarks/paper-a/grouped-xxz-family.json \
  benchmarks/paper-a/grouped-xxz-family-witness.json.gz
```

The family generator uses atomic writes, canonical row ordering, canonical JSON, and deterministic gzip (`mtime=0`, empty filename).  A fresh same-commit rebuild must match all three output files byte-for-byte:

```bash
tmpdir="$(mktemp -d)"
PYTHONPATH=src python scripts/run_xxz_transfer.py --build-grouped-family \
  --deltas 0/1,1/4,1/2,1/1,3/2,2/1,4/1 --length 4 \
  --summary "$tmpdir/grouped-xxz-family.json" \
  --witness "$tmpdir/grouped-xxz-family-witness.json.gz" \
  --transfer "$tmpdir/xxz-transfer.json"
cmp benchmarks/paper-a/grouped-xxz-family.json "$tmpdir/grouped-xxz-family.json"
cmp benchmarks/paper-a/grouped-xxz-family-witness.json.gz "$tmpdir/grouped-xxz-family-witness.json.gz"
cmp benchmarks/paper-a/xxz-transfer.json "$tmpdir/xxz-transfer.json"
shasum -a 256 "$tmpdir"/*
```

- [ ] **Step 3: Add the 12 by 12 scaling gate**

Implement a translation-local backend only after its unit-cell ledger is proven equal to direct finite-torus results on every support for which the finite embedding is injective.  For every canonical support, explicitly enumerate integer lattice coordinates, reduce them modulo L, and require: distinct support sites remain distinct; translated Pauli masks have the computed stabilizer; orbit size times stabilizer equals `L^2`; and no two distinct infinite-lattice translates alias to the same finite mask with different coefficients.  The 4 by 4 pilot records and rejects aliased supports instead of serving as the scaling proof.  Validate nonaliased records on the smallest admissible larger torus.

`certify_grouped_xxz_scaling.py` has exclusive `--build --length 12 --output PATH` and `--verify PATH` modes.  Its strict artifact schema contains formula/stage hash, Delta grid, every canonical D4-D7 support in integer coordinates, minimal injective comparison length, L=12 reduced mask, stabilizer, orbit, coefficient-map equality, weighted-growth/tail profile hash, source closure, and payload digest.  The primary verifier and the standard-library grouped-XXZ reference verifier independently replay injectivity, orbit/stabilizer, coefficient equality, and the fixed tail schema; unknown fields, duplicate keys, altered supports, or missing grid rows fail closed.

By default, the 12 by 12 row is `future_non_claim`.  Promote it only after a geometry-specific build independently rederives all of: the L=12 D4-D7 ledger, weighted all-order growth/tail profile, same-Delta baseline, grouped witness, accepted and adjacent-step inequalities, resources, and canonical artifact hashes.  The 12 by 12 summary/witness use their own schema and are verified by:

```bash
PYTHONPATH=src python scripts/run_xxz_transfer.py --build-grouped-family \
  --deltas 0/1,1/4,1/2,1/1,3/2,2/1,4/1 --length 12 \
  --scaling-certificate artifacts/verification/grouped-xxz-scaling-12x12.json \
  --summary benchmarks/paper-a/grouped-xxz-12x12.json \
  --witness benchmarks/paper-a/grouped-xxz-12x12-witness.json.gz
PYTHONPATH=src python scripts/certify_grouped_xxz_scaling.py --verify \
  artifacts/verification/grouped-xxz-scaling-12x12.json
python scripts/reference_verify_grouped_xxz.py \
  benchmarks/paper-a/grouped-xxz-12x12.json \
  benchmarks/paper-a/grouped-xxz-12x12-witness.json.gz \
  --scaling-certificate artifacts/verification/grouped-xxz-scaling-12x12.json
```

Before that command, build the scaling artifact from the frozen post-Task-9 source and verify it independently:

```bash
PYTHONPATH=src python scripts/certify_grouped_xxz_scaling.py --build --length 12 \
  --output artifacts/verification/grouped-xxz-scaling-12x12.json
PYTHONPATH=src python scripts/certify_grouped_xxz_scaling.py --verify \
  artifacts/verification/grouped-xxz-scaling-12x12.json
python scripts/reference_verify_grouped_xxz.py --verify-scaling \
  artifacts/verification/grouped-xxz-scaling-12x12.json
```

The L=12 summary and witness embed the scaling-artifact SHA-256 and source-commit binding; the reference verifier requires the supplied artifact, recomputes its digest, and replays it rather than treating the hash as authority.

If projected resources exceed the local threshold, use the frozen six-hour CPU profile and record the executable submission/collection sequence in the computation transcript, for example:

```bash
SOURCE_COMMIT="$(git rev-parse HEAD)"
JOB_ID="$(sbatch --parsable --job-name=issue128-xxz12 --time=06:00:00 \
  --cpus-per-task=16 --mem=64G --output=build/slurm/xxz12-%j.out \
  --error=build/slurm/xxz12-%j.err --wrap="srun python scripts/run_xxz_transfer.py --build-grouped-family --deltas 0/1,1/4,1/2,1/1,3/2,2/1,4/1 --length 12 --scaling-certificate artifacts/verification/grouped-xxz-scaling-12x12.json --summary build/slurm/grouped-xxz-12x12.json --witness build/slurm/grouped-xxz-12x12-witness.json.gz")"
sacct -j "$JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode > "build/slurm/xxz12-$JOB_ID.receipt"
shasum -a 256 "build/slurm/xxz12-$JOB_ID.out" "build/slurm/xxz12-$JOB_ID.err" \
  "build/slurm/xxz12-$JOB_ID.receipt" build/slurm/grouped-xxz-12x12*
```

The transcript also records `SOURCE_COMMIT`, the exact cluster partition/account supplied by the configured environment, collection host/time, and the post-collection primary/reference verifier outputs.  Failure to close any L=12 item leaves 12 by 12 outside the claim without weakening the complete L=4 family result.

- [ ] **Step 4: Promote the transfer claim**

The claim matrix may say positive XXZ transfer only if at least one `Delta!=1` row is certified.  A family statement requires every preregistered row to be reported and at least one easy-plane and one easy-axis certification.

- [ ] **Step 5: Commit Task 10**

Commit: `git commit -m "feat(issue128): freeze preregistered XXZ transfer family"`

### Task 11: Publication Ownership and Regression Gate

**Files:**
- Modify: `artifacts/publication/paper-a-file-ownership.json`
- Modify: `artifacts/publication/paper-a-claim-matrix.json`
- Modify: `scripts/audit_publication_scope.py`
- Modify: `tests/test_publication_scope.py`
- Modify: `tests/test_paper_a_claims.py`

- [ ] **Step 1: Add exact ownership entries**

Classify the grouped-XXZ source, generator, verifier, tests, summary, witness, family artifact, scaling builder/test/artifact, optional L=12 artifacts, and computation transcript as Paper A; include only the reproduction closure.  The scaling artifact is mandatory in ownership and release closure whenever any L=12 row is promoted.

- [ ] **Step 2: Add claim and scope attacks**

Require rejection when a manuscript claims transfer from only `Delta=1`, claims CNOT costs without a synthesis artifact, omits a failed preregistered row, or cites a 4 by 4 result as 12 by 12.

- [ ] **Step 3: Run the publication gate**

Run:

```bash
python -m pytest -q tests/test_grouped_xxz.py tests/test_reference_grouped_xxz.py \
  tests/test_xxz_transfer.py tests/test_paper_a_claims.py tests/test_publication_scope.py
PYTHONPATH=src python scripts/audit_paper_a_claims.py \
  artifacts/publication/paper-a-claim-matrix.json --root .
PYTHONPATH=src python scripts/audit_publication_scope.py \
  artifacts/publication/paper-a-file-ownership.json --root .
```

- [ ] **Step 4: Commit Task 11**

Commit: `git commit -m "chore(issue128): release grouped XXZ evidence"`

## Self-Review Checklist

- [ ] The pilot is 4 by 4 direct finite-torus, not a relabeled 12 by 12 artifact.
- [ ] The 12 by 12 anchor is not reused as nonisotropic evidence.
- [ ] Every Delta has its own baseline and witness.
- [ ] Schedule groups and anticommuting certificate groups are never conflated.
- [ ] The plan does not promise CNOT counts.
- [ ] At least one easy-plane and one easy-axis certification are required for a family claim.
- [ ] Honest failures stay in the release.
- [ ] Scaling to 12 by 12 is gated by direct-versus-local equivalence and alias detection.
