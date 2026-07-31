# D6 Physical-Channel Bound Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the independently reconstructed exact D6 Pauli map into a physically interpretable, rigorously anticommuting-grouped norm bound and use it to raise the adjacent-integer certified resource ratio when the full ledger permits.

**Architecture:** Keep the exact D6 payload immutable. A new pure module decodes its coordinate Pauli strings, constructs deterministic symplectic masks, classifies spatial clusters, discovers scalable local anticommuting pairs, and re-certifies complete coverage with the existing exact anticommutation machinery. A separate builder writes a hash-bound grouping sidecar; the main certificate builder and verifier consume and regenerate that sidecar before using its tighter site bound.

**Tech Stack:** Python 3.11, exact `Fraction` arithmetic, cubic-field outward intervals, canonical JSON/gzip, pytest, existing `trottercert.anticommuting` and `trottercert.exact_series_certificate` modules.

## Global Constraints

- Preserve the frozen 12x12 periodic Heisenberg benchmark, `T=1`, tolerance `1e-6`, 31-stage five-copy fourth-order Suzuki formula, and `30*r+1` merged-exponential accounting.
- A multiplier update requires two independent fanout layouts with identical complete exact D6 coefficient maps.
- Every D6 term must appear exactly once in the grouped partition; missing and duplicate terms are fatal.
- Every nonsingleton group must pass exact symplectic pairwise anticommutation checks.
- Cubic coefficients must be enclosed outward at the payload-declared decimal precision.
- Discovery may fail to improve the l1 bound; it may never weaken or replace a tighter valid bound.
- Report an improved result only as adjacent integers `r PASS / (r-1) FAIL` after the verifier chain.
- Do not modify or overwrite unrelated dirty worktree files.

---

### Task 1: Coordinate-cluster and scalable anticommuting-pair core

**Files:**
- Create: `tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/d6_physical_channels.py`
- Create: `tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_d6_physical_channels.py`

**Interfaces:**
- Consumes: `CoordinateCubicTerms` from `trottercert.hpc_artifacts` and the interval precision from the exact D6 payload.
- Produces: `PhysicalChannelSummary`, `GroupedD6Bound`, `classify_coordinate_pauli(pauli)`, and `build_grouped_d6_bound(terms, decimal_digits, candidate_cap=128)`.

- [ ] **Step 1: Write failing spatial-classification tests**

```python
def test_classifies_connected_path_and_disconnected_support() -> None:
    path = ((0, 0, "X"), (1, 0, "Y"), (2, 0, "Z"))
    disconnected = ((0, 0, "X"), (2, 0, "Z"))
    assert classify_coordinate_pauli(path).shape == "path"
    assert classify_coordinate_pauli(path).bbox == (3, 1)
    assert classify_coordinate_pauli(path).manhattan_diameter == 2
    assert classify_coordinate_pauli(disconnected).component_count == 2
```

- [ ] **Step 2: Run the classification test and confirm import failure**

Run: `PYTHONPATH=src:. pytest -q tests/test_d6_physical_channels.py::test_classifies_connected_path_and_disconnected_support`

Expected: FAIL because `trottercert.d6_physical_channels` does not exist.

- [ ] **Step 3: Implement deterministic coordinate classification**

Implement immutable records:

```python
@dataclass(frozen=True)
class ChannelClass:
    support_size: int
    component_count: int
    edge_count: int
    bbox: tuple[int, int]
    manhattan_diameter: int
    shape: str
```

Build the nearest-neighbor induced graph on occupied coordinates. Use exact labels: `singleton`, `bond`, `path`, `cycle`, `branched_tree`, `connected_other`, and `disconnected`. Pauli axes do not change spatial shape.

- [ ] **Step 4: Write failing complete-partition and improvement tests**

```python
def test_grouped_bound_pairs_local_anticommuting_channels() -> None:
    terms = {
        ((0, 0, "X"),): Cubic.one(),
        ((0, 0, "Z"),): Cubic.one(),
        ((2, 0, "X"),): Cubic.one(),
    }
    result = build_grouped_d6_bound(terms, decimal_digits=12, candidate_cap=8)
    assert result.term_count == 3
    assert sorted(len(group) for group in result.groups) == [1, 2]
    assert result.grouped_cell_bound < result.l1_cell_bound
```

Also add tests that two runs produce identical groups and that `candidate_cap < 1` raises `ValueError`.

- [ ] **Step 5: Implement scalable local-pair discovery and exact certification**

Construct one deterministic global coordinate-to-bit registry, convert coordinate Paulis to symplectic masks, and enclose every cubic coefficient with `cube_root_four_interval(decimal_digits)`. Maintain an inverted index from occupied coordinate to sorted unpaired term indices. For each coefficient-descending term, examine at most `candidate_cap` deterministic local candidates that overlap its support, pair with the first exact anticommuting candidate, and otherwise emit a singleton. Pass the resulting full partition to `certify_anticommuting_partition`; never trust discovery alone.

Return:

```python
@dataclass(frozen=True)
class GroupedD6Bound:
    term_count: int
    groups: tuple[tuple[int, ...], ...]
    l1_cell_bound: Fraction
    grouped_cell_bound: Fraction
    grouped_site_bound: Fraction
    max_group_size: int
    channel_counts: tuple[tuple[str, int], ...]
    channel_l1_bounds: tuple[tuple[str, Fraction], ...]
```

- [ ] **Step 6: Run focused tests**

Run: `PYTHONPATH=src:. pytest -q tests/test_d6_physical_channels.py tests/test_anticommuting.py`

Expected: all tests pass.

- [ ] **Step 7: Commit the core**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/d6_physical_channels.py tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_d6_physical_channels.py
git commit -m "feat(issue128): certify D6 physical-channel pairs"
```

### Task 2: Canonical grouped-D6 sidecar builder

**Files:**
- Create: `tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_d6_grouped_sidecar.py`
- Create: `tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_d6_grouped_sidecar.py`

**Interfaces:**
- Consumes: a verified `issue128_exact_right_generator_degree` gzip and its SHA-256.
- Produces: canonical JSON `issue128-d6-groups.json` containing source payload digest, exact coverage indices, certified group bounds, physical-channel counts, l1 baseline, and grouped cell/site bounds.

- [ ] **Step 1: Write a failing deterministic-sidecar test**

Create a three-term degree-six fixture using `write_shard_gzip`, invoke `build_sidecar(input_path, output_path, candidate_cap=8)`, and assert:

```python
assert sidecar["kind"] == "issue128_d6_physical_channel_groups"
assert sidecar["source_payload_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
assert sidecar["term_count"] == 3
assert sidecar["grouped_site_bound"] < sidecar["site_pauli_l1_upper"]
assert sum(len(group["term_indices"]) for group in sidecar["groups"]) == 3
```

- [ ] **Step 2: Run the test and confirm import failure**

Run: `PYTHONPATH=src:. pytest -q tests/test_d6_grouped_sidecar.py`

Expected: FAIL because the builder does not exist.

- [ ] **Step 3: Implement the builder**

Use `read_portable_canonical_gzip` and `verify_exact_degree_payload` before grouping. Encode every `Fraction` as `[numerator, denominator]`, sort all dictionary keys, write through a temporary sibling followed by `Path.replace`, and include `source_commit`, `coefficient_interval_decimal_digits`, `candidate_cap`, group count, maximum group size, channel statistics, and both l1/grouped bounds.

- [ ] **Step 4: Add tamper and reproducibility tests**

Build twice and require byte-identical output. Change a source coefficient without fixing its l1 field and require rejection. Duplicate or remove a saved group index and verify the future verifier fixture rejects it.

- [ ] **Step 5: Run focused tests**

Run: `PYTHONPATH=src:. pytest -q tests/test_d6_grouped_sidecar.py tests/test_exact_series_certificate.py`

Expected: all tests pass.

- [ ] **Step 6: Commit the builder**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_d6_grouped_sidecar.py tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_d6_grouped_sidecar.py
git commit -m "feat(issue128): build grouped exact D6 sidecar"
```

### Task 3: Verifier regeneration and main-certificate integration

**Files:**
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/verify.py`
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_d6_integrated_certificate.py`
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_exact_series_certificate.py`
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_resources_verify.py`

**Interfaces:**
- Consumes: exact D6 gzip, parent artifact, and grouped-D6 JSON sidecar.
- Produces: `D6SidecarVerification.site_bound` equal to the regenerated grouped bound and a main certificate whose D6 metadata binds both artifacts by path and SHA-256.

- [ ] **Step 1: Extend the verifier fixture with a valid grouped sidecar**

Add `groups_path`, `groups_sha256`, `group_count`, `max_group_size`, `l1_site_norm_upper`, and grouped `site_norm_upper` to `candidate["d6_certificate"]`. Assert `_verify_d6_sidecar` regenerates the partition and returns the grouped site bound.

- [ ] **Step 2: Add failing corruption tests**

Require rejection for: escaping group path, wrong group digest, duplicate index, missing index, nonanticommuting pair, reduced stored group bound, source-payload digest mismatch, and grouped bound larger than the verified l1 baseline.

- [ ] **Step 3: Run the new verifier tests and observe failure**

Run: `PYTHONPATH=src:. pytest -q tests/test_exact_series_certificate.py tests/test_resources_verify.py -k 'd6 or grouped'`

Expected: FAIL because grouped-D6 metadata is not yet consumed.

- [ ] **Step 4: Implement fail-closed grouped-D6 regeneration**

Resolve the group sidecar inside the certificate directory, verify its digest and source-payload digest, reconstruct intervals and symplectic masks from the exact D6 terms, rebuild the listed partition with `certify_anticommuting_partition`, and require exact equality of coverage, group bounds, group count, maximum size, and cell/site total bounds.

- [ ] **Step 5: Integrate the grouped bound into integer search**

Modify `build_d6_integrated_certificate.py` to require the grouping sidecar, use `grouped.site_bound` as `d6_site_override`, preserve the exact-D6 l1 baseline in metadata, and record the physical-channel method in `proof_method`. Continue searching integers from 2 upward, then record both the accepted error and the previous-step error.

- [ ] **Step 6: Run focused integration tests**

Run: `PYTHONPATH=src:. pytest -q tests/test_d6_grouped_sidecar.py tests/test_exact_series_certificate.py tests/test_resources_verify.py tests/test_refined_error.py`

Expected: all tests pass, including every tamper rejection.

- [ ] **Step 7: Commit integration**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/verify.py tracks/qcs/solutions/WangTheoPhys/issue128/scripts/build_d6_integrated_certificate.py tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_exact_series_certificate.py tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_resources_verify.py
git commit -m "feat(issue128): verify grouped D6 certificate"
```

### Task 4: Production artifact, adjacent-integer audit, and delivery

**Files:**
- Create: `tracks/qcs/solutions/WangTheoPhys/issue128/certificates/issue128-d6-groups.json`
- Modify only on successful certification: `tracks/qcs/solutions/WangTheoPhys/issue128/certificates/issue128-d6-integrated-certificate.json`
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/d5-integrated/verification-transcript.txt`
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/docs/report/issue128-hpc-sprint-log.md`
- Modify: `tracks/qcs/solutions/WangTheoPhys/issue128/docs/presentation/issue128-10-minute-talk.md`

**Interfaces:**
- Consumes: the completed cross-lane audit and exact D6 artifact.
- Produces: the final verified resource ratio or an explicit negative/provisional result without altering the frozen certificate.

- [ ] **Step 1: Gate on two-lane equality**

Require the dual reducer to finish with both payloads, an exact comparison JSON reporting equality, and both SHA-256 values. If this gate fails, generate no production grouping sidecar and retain `4.1357418449...x`.

- [ ] **Step 2: Build the production grouped sidecar**

Run from the Issue #128 directory:

```bash
PYTHONPATH=src:. python -u scripts/build_d6_grouped_sidecar.py \
  --input certificates/issue128-d6-exact.json.gz \
  --output certificates/issue128-d6-groups.json \
  --candidate-cap 128
```

Record term count, channel counts, l1 site bound, grouped site bound, reduction fraction, group count, and maximum group size.

- [ ] **Step 3: Build and inspect the integer certificate**

Run:

```bash
PYTHONPATH=src:. python -u scripts/build_d6_integrated_certificate.py
PYTHONPATH=src:. python -u scripts/verify.py --fast
```

Require `candidate_error_upper <= 1e-6` and `previous_step_error_upper > 1e-6`. Compute the exact ratio as `Fraction(11791, 30*r+1)`.

- [ ] **Step 4: Run the verifier chain**

Run:

```bash
PYTHONPATH=src:. pytest -q -m 'not slow'
PYTHONPATH=src:. python -u scripts/verify.py --deep
```

Then copy the repository to a clean temporary worktree, copy only hash-bound production artifacts, and repeat fast verification. Any incomplete deep/clean gate leaves the result labeled candidate.

- [ ] **Step 5: Freeze evidence and explain the physical gain**

Append exact commands, test counts, artifact hashes, l1/grouped D6 bounds, all ledger contributions, adjacent integer results, and exact resource ratio to the transcript and sprint log. In the talk, describe the gain as a strict norm reduction from mutually incompatible multi-spin exchange-error channels, not as reducer acceleration.

- [ ] **Step 6: Update the PR only with the verified status**

Post one comment distinguishing one of three outcomes: certified improved ratio; candidate pending deep/clean reproduction; or no integer crossing with the frozen `4.1357418449...x` retained. Include both cross-lane hashes and the adjacent `PASS/FAIL` boundary for any certified claim.

- [ ] **Step 7: Commit production evidence**

Stage only the files changed by this feature, inspect `git diff --cached`, and commit with:

```bash
git commit -m "cert(issue128): tighten D6 physical-channel bound"
```

## Self-review

- Spec coverage: Tasks 1-4 cover physical classification, scalable discovery, exact certification, dual-lane trust, sidecar binding, ledger integration, adjacent-integer search, verifier/tamper gates, narrative, and PR delivery.
- Placeholder scan: no deferred implementation or unspecified error handling remains.
- Type consistency: `GroupedD6Bound.grouped_site_bound` is the value consumed by the sidecar, verifier, and ledger; all artifact fractions use `[numerator, denominator]` encoding.
