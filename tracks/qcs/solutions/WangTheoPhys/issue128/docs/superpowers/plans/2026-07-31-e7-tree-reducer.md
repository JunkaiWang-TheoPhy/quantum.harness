# E7 Checkpointed Tree Reducer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and deploy a fail-closed two-level E7 reducer that preserves reusable partial checkpoints and emits the existing exact D6 artifact format.

**Architecture:** Slurm array cells merge disjoint contiguous source-shard ranges into canonical partial E7 checkpoints. A dependent final job validates complete coverage, merges the checkpoints, constructs exact D6, and emits hash-bound parent, payload, and summary artifacts.

**Tech Stack:** Python 3.11, exact `Fraction`/`Cubic` arithmetic, deterministic gzip/JSON, pytest, Slurm.

## Global Constraints

- Do not cancel or modify the running monolithic dual reducer.
- Do not accept missing, corrupt, overlapping, gapped, or mixed-source inputs.
- Write every checkpoint and manifest atomically.
- Do not promote a new multiplier without two different fanout layouts and exact coefficient-map equality.
- Preserve unrelated tracked and untracked workspace changes.

---

### Task 1: Partial checkpoint reducer

**Files:**
- Create: `scripts/reduce_e7_word_range.py`
- Test: `tests/test_e7_tree_reducer.py`

**Interfaces:**
- Consumes: source word-shard directory, total shard count, inclusive start index, exclusive stop index.
- Produces: `build_partial_checkpoint(...) -> dict[str, object]`, a canonical `issue128_exact_e7_partial` payload, and a complete manifest.

- [ ] **Step 1: Write failing successful-range and digest-corruption tests**

Create two synthetic source shards with contiguous word coverage and exact
cubic terms. Assert their sum and source digests in the partial payload, then
mutate a payload after manifest creation and require `ValueError`.

- [ ] **Step 2: Run the focused tests and observe the missing module failure**

Run: `pytest -q tests/test_e7_tree_reducer.py`
Expected: collection fails because `scripts.reduce_e7_word_range` is absent.

- [ ] **Step 3: Implement validated range reduction**

Validate `0 <= start < stop <= shard_count`, manifest status/digest, payload
kind/index/count, identical formula/source/stage/degree metadata, and contiguous
word coverage. Merge with exact `Cubic` addition and write the payload and
manifest through `write_shard_gzip` and `write_manifest_atomic`.

- [ ] **Step 4: Run focused tests**

Run: `pytest -q tests/test_e7_tree_reducer.py`
Expected: all Task 1 tests pass.

### Task 2: Final checkpoint reducer

**Files:**
- Create: `scripts/finalize_e7_partials.py`
- Modify: `tests/test_e7_tree_reducer.py`

**Interfaces:**
- Consumes: partial directory and expected partial/source-shard counts.
- Produces: the established exact D6 payload plus parent and summary manifests.

- [ ] **Step 1: Write failing finalization and coverage-gap tests**

Create two partials covering all synthetic source shards; assert final exact D6
metadata and parent provenance. Remove or alter one interval and require a
gap/overlap failure before any final D6 file exists.

- [ ] **Step 2: Implement finalization**

Validate partial payload digests and metadata, require ordered intervals to
cover `[0, source_shard_count)` and word intervals to be contiguous, merge exact
E7 terms, construct `D6=7E7+(2/3)ad_A^2(E5)`, and atomically write the same
degree-six payload fields used by `reduce_e7_word_shards.py`.

- [ ] **Step 3: Run focused and existing reducer tests**

Run: `pytest -q tests/test_e7_tree_reducer.py tests/test_compare_exact_degree_payloads.py`
Expected: all tests pass.

### Task 3: Slurm deployment wrappers

**Files:**
- Create: `hpc/issue128_e7_tree_partial.sbatch`
- Create: `hpc/issue128_e7_tree_finalize.sbatch`

**Interfaces:**
- Consumes: `ISSUE128_SOURCE_RUN_ID`, `ISSUE128_TREE_RUN_ID`,
  `ISSUE128_SOURCE_SHARD_COUNT`, `ISSUE128_PARTIAL_COUNT`.
- Produces: a partial array and a dependency-gated final D6 reducer.

- [ ] **Step 1: Add fail-closed array wrapper**

Map each array index `i` to
`start=floor(S*i/P), stop=floor(S*(i+1)/P)` and invoke
`reduce_e7_word_range.py` with deterministic environment variables.

- [ ] **Step 2: Add final wrapper**

Invoke `finalize_e7_partials.py` only after the array dependency succeeds.
Request enough memory for exact D6 construction but only one CPU because the
final Python merge is single-process.

- [ ] **Step 3: Validate shell syntax**

Run: `bash -n hpc/issue128_e7_tree_partial.sbatch hpc/issue128_e7_tree_finalize.sbatch`
Expected: exit code 0.

### Task 4: Deploy and harvest

**Files:**
- Read/write remotely: `/work/home/acamtw70yu/quantum-harness-issue128-hpc/results/`

**Interfaces:**
- Consumes: complete 544- and 1,088-shard producer lanes.
- Produces: one or two complete exact D6 outputs and an exact comparison audit.

- [ ] **Step 1: Sync only committed tree-reducer files and verify remote hashes**

- [ ] **Step 2: Smoke a small number of partial ranges against existing complete inputs**

- [ ] **Step 3: Submit production arrays and `afterok` finalizers**

- [ ] **Step 4: Validate summaries, hashes, exact-map equality, and certificate integration**

- [ ] **Step 5: Retain the existing 4.1357418449× claim unless all promotion gates pass**

## Self-Review

- The plan covers atomic partials, exact coverage, corruption rejection,
  deterministic finalization, Slurm deployment, and promotion gates.
- The two scripts have distinct responsibilities and stable artifact
  interfaces.
- No placeholder implementation or unrelated scientific route is included.

