# Dual E9 Manifest Contraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a hash-bound, pre-sharded degree-nine word manifest and use it to compute the exact dual pairing `q9=tau(W E9)/(L^2/4)` without constructing the full E9 operator or rebuilding the formal-log map in every worker.

**Architecture:** A preparation command generates the exact degree-nine free-log map once and publishes 64 canonical word manifests plus a digest index. A separate manifest-driven Pauli evaluator produces exact dual-pairing shard artifacts, and a fail-closed reducer verifies full coverage and emits `q9`. Existing frozen E5/E7 code and artifacts remain unchanged; the new path is independently regressed against them.

**Tech Stack:** Python 3.12, exact `fractions.Fraction`, the existing `Cubic` field and symplectic Pauli evaluator, canonical JSON, SHA-256, pytest, Ruff, Slurm.

## Global Constraints

- Do not modify `src/trottercert/dual_log_pairing.py` or any committed E7 shard; their hashes are frozen evidence.
- Do not construct or claim the full E9 operator.
- Do not promote `finite_step_status` from `inconclusive`; exact E9 does not certify the E11-and-higher tail.
- Every mathematical artifact uses canonical uncompressed JSON bytes terminated by one newline.
- Every rational is encoded canonically as `[numerator, positive_denominator]`.
- A partial or interrupted manifest/worker set is rejected, never interpreted as zero.
- The production degree is 9, torus length is 6, and production shard count is 64.
- The formula is the frozen five-copy fourth-order Suzuki formula over four matchings.
- Preserve all unrelated dirty-worktree changes and stage only files named by each task.
- Every commit uses the repository Lore body and `Co-authored-by: OmX <omx@oh-my-codex.dev>` trailer.

---

### Task 1: Canonical exact word manifests

**Files:**
- Create: `src/trottercert/dual_word_manifest.py`
- Create: `scripts/prepare_dual_log_manifests.py`
- Create: `tests/test_dual_word_manifest.py`

**Interfaces:**
- Consumes: `cubic_formula_log_series(stages, degree)[degree]` and `Sequence[CubicStage]`.
- Produces: `WordRecord`, `WordGroup`, `WordManifest`, `ManifestIndex`, `build_word_manifests(stages, degree, shard_count)`, `write_manifest_set(...)`, `load_word_manifest(path, index)`, and `verify_manifest_index(index, directory)`.

- [ ] **Step 1: Write failing model and assignment tests**

```python
def test_manifest_group_assignment_is_canonical_and_complete() -> None:
    stages = fourth_order_suzuki_cubic_stages(4)
    manifests = build_word_manifests(stages, degree=5, shard_count=3)
    groups = [group for manifest in manifests for group in manifest.groups]
    assert [group.ordinal for group in sorted(groups, key=lambda item: item.ordinal)] == list(
        range(manifests[0].total_groups)
    )
    assert all(
        group.ordinal % 3 == manifest.shard_index
        for manifest in manifests
        for group in manifest.groups
    )
    assert all(tuple(sorted(group.records, key=lambda item: item.word)) == group.records for group in groups)


def test_manifest_rejects_changed_word_or_coefficient(tmp_path: Path) -> None:
    index_path = write_manifest_set(
        tmp_path,
        fourth_order_suzuki_cubic_stages(4),
        degree=5,
        shard_count=2,
        implementation_sources={"source.py": "a" * 64},
    )
    index = load_manifest_index(index_path)
    shard_path = tmp_path / index.shards[0].path
    payload = json.loads(shard_path.read_text())
    payload["groups"][0]["records"][0]["word"] = "99999"
    shard_path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")
    with pytest.raises(ValueError, match="word|digest"):
        verify_manifest_index(index, tmp_path)
```

- [ ] **Step 2: Run the tests and confirm the missing-module failure**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_word_manifest.py
```

Expected: collection fails because `trottercert.dual_word_manifest` does not exist.

- [ ] **Step 3: Implement immutable manifest models and canonical codecs**

Implement these exact public records:

```python
@dataclass(frozen=True, slots=True)
class WordRecord:
    word: tuple[int, ...]
    coefficient: Cubic


@dataclass(frozen=True, slots=True)
class WordGroup:
    ordinal: int
    suffix: tuple[int, ...]
    records: tuple[WordRecord, ...]


@dataclass(frozen=True, slots=True)
class WordManifest:
    degree: int
    shard_index: int
    shard_count: int
    total_groups: int
    groups: tuple[WordGroup, ...]
    implementation_sources: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class ManifestShardRecord:
    shard_index: int
    path: str
    sha256: str
    group_count: int
    word_count: int


@dataclass(frozen=True, slots=True)
class ManifestIndex:
    degree: int
    shard_count: int
    total_groups: int
    total_words: int
    formula: str
    implementation_sources: tuple[tuple[str, str], ...]
    shards: tuple[ManifestShardRecord, ...]
```

Encode words as base-four digit strings and reject any digit outside `0..3` or
any word whose length differs from `degree`.  Reuse the canonical rational and
cubic rules from `certify_dual_e7_pairing.py`, but implement them locally so
the frozen E7 script is not modified.

- [ ] **Step 4: Implement deterministic group construction**

Use the exact ordering and assignment:

```python
word_map = cubic_formula_log_series(stages, degree)[degree]
ordered = sorted(word_map.items(), key=lambda item: (item[0][1:], item[0][0]))
raw_groups = groupby(ordered, key=lambda item: item[0][1:])
for ordinal, (suffix, items) in enumerate(raw_groups):
    shard_index = ordinal % shard_count
```

Require odd `degree >= 3`, `shard_count >= 1`, nonzero coefficients, group sizes
between one and four, a common suffix within each group, and globally unique
words.  `build_word_manifests` returns manifests ordered by shard index.

- [ ] **Step 5: Implement fail-closed two-phase publication**

`write_manifest_set` must:

1. create a temporary sibling directory with `tempfile.mkdtemp`;
2. write `manifest-000.json` through `manifest-063.json` canonically;
3. reread and validate every manifest;
4. compute each SHA-256 and construct `index.json`;
5. validate the complete index and all child digests; and
6. rename the temporary directory to the requested new output directory.

Reject an existing output directory rather than overwriting it.  This makes an
interrupted preparation visibly incomplete and prevents a mixed generation.

- [ ] **Step 6: Implement the preparation CLI**

The command line is:

```bash
PYTHONPATH=src:. python scripts/prepare_dual_log_manifests.py \
  --degree 9 \
  --shard-count 64 \
  --output-dir artifacts/dual-e9/manifests
```

Support a separate verification mode:

```bash
PYTHONPATH=src:. python scripts/prepare_dual_log_manifests.py \
  --verify artifacts/dual-e9/manifests/index.json
```

The source manifest must bind `cubic_field.py`, `cubic_local.py`, and the new
`dual_word_manifest.py`.

- [ ] **Step 7: Run focused tests and a degree-five CLI round trip**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_word_manifest.py
tmp_dir=$(mktemp -d)
PYTHONPATH=src:. python scripts/prepare_dual_log_manifests.py \
  --degree 5 --shard-count 4 --output-dir "$tmp_dir/manifests"
PYTHONPATH=src:. python scripts/prepare_dual_log_manifests.py \
  --verify "$tmp_dir/manifests/index.json"
```

Expected: tests pass; the CLI reports exact complete coverage and verification.

- [ ] **Step 8: Commit Task 1**

Stage only the three Task-1 files and commit with subject:

```text
feat(issue128): add canonical dual word manifests
```

The Lore body must record the degree-five round trip and state that no E9
production artifact has been generated.

---

### Task 2: Manifest-driven exact dual contraction

**Files:**
- Create: `src/trottercert/dual_manifest_pairing.py`
- Create: `tests/test_dual_manifest_pairing.py`
- Read for frozen reference only: `src/trottercert/dual_log_pairing.py`
- Read for frozen record only: `docs/experiments/processor-obstruction/dual-e7-pairing.json`

**Interfaces:**
- Consumes: one verified `WordManifest` and `length=6`.
- Produces: `ManifestPairingPartial` and `contract_word_manifest(manifest, *, length=6, progress=None)`.

- [ ] **Step 1: Write failing exact regression tests**

```python
@pytest.mark.slow
def test_degree_five_manifest_pairing_matches_frozen_record(tmp_path: Path) -> None:
    manifests = build_word_manifests(
        fourth_order_suzuki_cubic_stages(4), degree=5, shard_count=4
    )
    observed = [contract_word_manifest(manifest, length=6) for manifest in manifests]
    tau_h = sum((item.tau_h for item in observed), Cubic.zero())
    tau_h2 = sum((item.tau_h2 for item in observed), Cubic.zero())
    tau_w = sum((item.tau_w for item in observed), Cubic.zero())
    expected = next(
        record
        for record in json.loads(EXTENSIVE.read_text())["records"]
        if record["length"] == 6
    )
    assert tau_h == parse_cubic(expected["moments"]["tau_h_e5"])
    assert tau_h2 == parse_cubic(expected["moments"]["tau_h2_e5"])
    assert tau_w == parse_cubic(expected["moments"]["tau_w_e5"])


def test_manifest_pairing_rejects_wrong_degree_metadata() -> None:
    manifest = replace(valid_manifest(), degree=7)
    with pytest.raises(ValueError, match="word length|degree"):
        contract_word_manifest(manifest)
```

- [ ] **Step 2: Run the tests and confirm the missing-module failure**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_manifest_pairing.py
```

Expected: collection fails because `trottercert.dual_manifest_pairing` does not exist.

- [ ] **Step 3: Implement the independent manifest evaluator**

Define:

```python
@dataclass(frozen=True, slots=True)
class ManifestPairingPartial:
    degree: int
    length: int
    shard_index: int
    shard_count: int
    total_groups: int
    group_indices: tuple[int, ...]
    word_count: int
    nonzero_word_count: int
    retained_term_count: int
    tau_h: Cubic
    tau_h2: Cubic
    tau_w: Cubic
```

Copy no serializer or reducer logic.  Implement only the contraction using the
same trusted low-level objects as the frozen path:

```python
lattice = SquareLattice(length)
hamiltonian = heisenberg_symplectic_terms(lattice)
squared = square_real_pauli_terms(hamiltonian)
evaluator = SymplecticDyadicLocalDensityEvaluator(shared_coordinates=True)
denominator = degree * (1 << evaluator.denominator_exponent((0,) * degree))
```

For every manifest word, call `evaluator.evaluate(word)`, retain support at
most four, lift coordinates to `L=6` with explicit alias rejection, and apply
the exact cubic word coefficient only after rational Pauli contraction.  Clear
the cache after every suffix group.

- [ ] **Step 4: Add an in-memory path equivalence test on selected groups**

Build one degree-five manifest with at least two nonadjacent group ordinals.
Compare its result with the sum of two one-group manifests constructed from the
same records.  Require exact equality of all pairings and additive counters.
This catches accidental cache dependence across groups without rerunning the
full frozen E5 test.

- [ ] **Step 5: Run focused and slow exact regressions**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_manifest_pairing.py -m 'not slow'
PYTHONPATH=src:. pytest -q tests/test_dual_manifest_pairing.py
```

Expected: both pass; the slow test reproduces the frozen degree-five cubic
pairings exactly.

- [ ] **Step 6: Commit Task 2**

Stage only the two Task-2 files and commit with subject:

```text
feat(issue128): contract pre-sharded dual manifests
```

The Lore constraint must say that the frozen E7 core was not modified.

---

### Task 3: E9 worker artifacts and exact reducer

**Files:**
- Create: `scripts/certify_dual_e9_pairing.py`
- Create: `tests/test_dual_e9_pairing.py`
- Read for codec behavior: `scripts/certify_dual_e7_pairing.py`

**Interfaces:**
- Consumes worker mode `--index INDEX --manifest MANIFEST --output OUTPUT`.
- Consumes reducer mode `--index INDEX --reduce SHARD... --output OUTPUT`.
- Produces verified `issue128_dual_e9_pairing_shard` and `issue128_dual_e9_pairing` JSON artifacts.

- [ ] **Step 1: Write failing payload and reducer mutation tests**

```python
def test_e9_worker_payload_round_trips() -> None:
    payload = build_worker_payload(sample_partial(), sample_manifest_record(), SOURCES)
    verify_worker_payload(payload, sample_index())
    assert payload["degree"] == 9
    assert payload["claim"]["full_e9_operator"] == "not_computed"
    assert payload["claim"]["dual_e9_pairing"] == "exact"


def test_e9_reducer_rejects_manifest_and_group_attacks() -> None:
    shards, index = complete_sample_shards()
    forged = copy.deepcopy(shards)
    forged[0]["manifest_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="manifest"):
        reduce_worker_payloads(forged, index=index, parent_sha256=PARENT_HASHES)
    missing = shards[:-1]
    with pytest.raises(ValueError, match="complete|coverage"):
        reduce_worker_payloads(missing, index=index, parent_sha256=PARENT_HASHES[:-1])
```

- [ ] **Step 2: Run the tests and confirm the missing-script failure**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_e9_pairing.py
```

Expected: collection fails because `scripts.certify_dual_e9_pairing` does not exist.

- [ ] **Step 3: Implement worker payload validation**

The worker payload must include:

```text
schema_version=1
kind=issue128_dual_e9_pairing_shard
degree=9
length=6
shard_index, shard_count=64, total_groups=65536, group_indices
word_count, nonzero_word_count, retained_term_count
pairings={tau_h,tau_h2,tau_w}
manifest_index_sha256, manifest_sha256
implementation_sources, implementation_sources_digest
runtime={wall_seconds, peak_rss_bytes, scheduler}
claim={full_e9_operator:not_computed, dual_e9_pairing:exact,
       finite_step_status:inconclusive}
```

`runtime.scheduler` may be `null` locally or contain only `job_id`,
`array_task_id`, `cluster`, and `host`.  Exclude runtime metadata from the
mathematical payload digest and record both the full-file and mathematical
digests.

- [ ] **Step 4: Implement the reducer**

Require all 64 worker shards, each matching the index's shard digest and group
ordinals.  Verify the union is exactly `0..65535`, source manifests match, and
the sum is independent of forward/reverse shard order.  Emit:

```text
pairings.tau_h
pairings.tau_h2
pairings.tau_w
pairings.tau_w_per_cell = tau_w / 9
pairings.q9_over_97_pow_8
pairings.exact_e5_e7_e9_per_cell_at_r97
claim.finite_step_status = inconclusive
claim.missing = E11-and-higher dual tail
```

Load q5 from `extensive-commutant-witness.json` and q7 from
`dual-e7-pairing.json`, binding both parent SHA-256 values in the reduced
artifact.  Reject any attempt to set a finite-step pass in worker or reducer
payloads.

- [ ] **Step 5: Implement worker, reducer, and verifier CLI modes**

Worker:

```bash
PYTHONPATH=src:. python scripts/certify_dual_e9_pairing.py \
  --index artifacts/dual-e9/manifests/index.json \
  --manifest artifacts/dual-e9/manifests/manifest-000.json \
  --output artifacts/dual-e9/shards/shard-000.json
```

Reducer:

```bash
PYTHONPATH=src:. python scripts/certify_dual_e9_pairing.py \
  --index artifacts/dual-e9/manifests/index.json \
  --reduce artifacts/dual-e9/shards/shard-*.json \
  --output artifacts/dual-e9/dual-e9-pairing.json
```

Verifier:

```bash
PYTHONPATH=src:. python scripts/certify_dual_e9_pairing.py \
  --verify artifacts/dual-e9/dual-e9-pairing.json
```

- [ ] **Step 6: Run mutation tests and a local degree-nine one-group smoke**

The smoke manifest must be produced by the trusted degree-nine preparer and
contain exactly one suffix group.  Run the worker with a test-only index whose
coverage is exactly that group, then verify the worker artifact.  Do not call
the production reducer on the incomplete production set.

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_e9_pairing.py
```

Expected: all payload, mutation, one-group smoke, and fail-closed claim tests pass.

- [ ] **Step 7: Commit Task 3**

Stage only the Task-3 script and test and commit with subject:

```text
feat(issue128): add fail-closed dual E9 reducer
```

The Lore body must state that no complete E9 artifact or all-order tail is
claimed.

---

### Task 4: Slurm pipeline and local calibration

**Files:**
- Create: `hpc/issue128_e9_manifest.sbatch`
- Create: `hpc/issue128_e9_array.sbatch`
- Create: `hpc/issue128_e9_reduce.sbatch`
- Create: `docs/report/dual-e9-hpc-runbook.md`
- Create: `tests/test_dual_e9_hpc_scripts.py`

**Interfaces:**
- Consumes: `ISSUE128_ROOT`, `ISSUE128_E9_RUN_ROOT`, and standard Slurm array variables.
- Produces: manifest job, dependent `0-63` worker array, dependent reducer, and deterministic verification commands.

- [ ] **Step 1: Write failing static Slurm tests**

```python
@pytest.mark.parametrize(
    "name",
    [
        "issue128_e9_manifest.sbatch",
        "issue128_e9_array.sbatch",
        "issue128_e9_reduce.sbatch",
    ],
)
def test_e9_slurm_script_is_fail_closed(name: str) -> None:
    text = (REPO_ROOT / "hpc" / name).read_text()
    assert "set -euo pipefail" in text
    assert "PYTHONPATH=src:." in text
    assert "ISSUE128_E9_RUN_ROOT" in text
    assert "--verify" in text or "--index" in text


def test_array_has_exact_production_shape() -> None:
    text = (REPO_ROOT / "hpc/issue128_e9_array.sbatch").read_text()
    assert "#SBATCH --array=0-63" in text
    assert "#SBATCH --cpus-per-task=1" in text
    assert "#SBATCH --mem=4G" in text
    assert "#SBATCH --time=02:00:00" in text
```

- [ ] **Step 2: Run the static tests and confirm missing-file failures**

Run:

```bash
PYTHONPATH=src:. pytest -q tests/test_dual_e9_hpc_scripts.py
```

Expected: failures naming the three absent `.sbatch` files.

- [ ] **Step 3: Implement the three noninteractive Slurm scripts**

Every script must start with:

```bash
#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
set -euo pipefail
: "${ISSUE128_ROOT:?set ISSUE128_ROOT to the Issue-128 directory}"
: "${ISSUE128_E9_RUN_ROOT:?set ISSUE128_E9_RUN_ROOT to a new run directory}"
cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."
```

The manifest job requests 4 GiB and 30 minutes and refuses an existing
`$ISSUE128_E9_RUN_ROOT/manifests`.  The array requests `0-63`, 4 GiB, and two
hours and writes `shard-%03d.json` through a temporary file followed by an
atomic rename.  The reducer requests 8 GiB and one hour, verifies all 64 shard
files before reducing, writes its output atomically, and verifies it again.

- [ ] **Step 4: Write the exact submission and recovery runbook**

The runbook must give commands using parsed job IDs:

```bash
manifest_job=$(sbatch --parsable hpc/issue128_e9_manifest.sbatch)
array_job=$(sbatch --parsable --dependency="afterok:${manifest_job}" hpc/issue128_e9_array.sbatch)
reduce_job=$(sbatch --parsable --dependency="afterok:${array_job}" hpc/issue128_e9_reduce.sbatch)
```

It must also document `sacct` inspection, detection of missing shards, safe
resubmission of only failed array indices, one independent shard rerun into a
new directory, SHA comparison, artifact download, and local verifier commands.
Do not include a real account, partition, hostname, username, or secret.

- [ ] **Step 5: Run bash syntax and static tests**

Run:

```bash
bash -n hpc/issue128_e9_manifest.sbatch
bash -n hpc/issue128_e9_array.sbatch
bash -n hpc/issue128_e9_reduce.sbatch
PYTHONPATH=src:. pytest -q tests/test_dual_e9_hpc_scripts.py
```

Expected: all pass.

- [ ] **Step 6: Run local calibration before requesting submission**

Generate a degree-nine test manifest set with 64 shards, choose the smallest
nonempty shard prefix that contains at least one group, and time one complete
group contraction with `/usr/bin/time -l`.  Record measured wall time and peak
RSS in the runbook.  Scale only as a scheduler request estimate; do not claim
the extrapolation as completed production work.

- [ ] **Step 7: Run the complete verification gate**

Run:

```bash
python -m ruff check --select I,UP035 \
  src/trottercert/dual_word_manifest.py \
  src/trottercert/dual_manifest_pairing.py \
  scripts/prepare_dual_log_manifests.py \
  scripts/certify_dual_e9_pairing.py \
  tests/test_dual_word_manifest.py \
  tests/test_dual_manifest_pairing.py \
  tests/test_dual_e9_pairing.py \
  tests/test_dual_e9_hpc_scripts.py
PYTHONPATH=src:. pytest -q \
  tests/test_dual_word_manifest.py \
  tests/test_dual_manifest_pairing.py \
  tests/test_dual_e9_pairing.py \
  tests/test_dual_e9_hpc_scripts.py
PYTHONPATH=src:. pytest -q -m 'not slow'
git diff --check
```

Expected: all focused and complete non-slow tests pass; only known slow tests
are deselected.

- [ ] **Step 8: Commit Task 4**

Stage only the four Task-4 files and commit with subject:

```text
ops(issue128): prepare dual E9 Slurm pipeline
```

The Lore constraint must state that the scripts contain no site credentials
and that no job has been submitted.

---

## Plan self-review

- **Spec coverage:** Tasks 1--4 cover canonical manifests, atomic publication,
  exact manifest-driven contraction, complete reduction, mutation resistance,
  frozen E5/E7 regressions, local E9 smoke, Slurm execution, recovery, and
  independent rerun evidence.  The separate E11+ mathematical majorant is
  deliberately excluded and remains the next design/plan after exact E9 data.
- **Placeholder scan:** Every implementation and validation step names its
  concrete interface, behavior, command, and expected result.
- **Type consistency:** `WordManifest` is produced in Task 1, consumed in Task
  2, and represented by a digest-bound file/index pair in Task 3.  Worker and
  reducer CLI names match the Slurm scripts in Task 4.
- **Claim boundary:** Every task retains `finite_step_status=inconclusive`; only
  a later independently verified E11+ tail can promote it.
