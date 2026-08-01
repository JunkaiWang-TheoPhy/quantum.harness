# Dual-Only E7 Pairing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a deterministic, exact, four-shard certificate for `tau(W E7)` without constructing the full E7 operator.

**Architecture:** A reusable contraction core groups free-log words by common suffix and contracts only target-compatible Pauli strings. A shard CLI emits exact partial sums; a reducer proves coverage and sums in both orders before a ledger compares E5 and E7 dual contributions.

**Tech Stack:** Python 3.12, exact `Fraction`/`Cubic`, symplectic Pauli maps, canonical JSON, SHA-256, pytest

## Global Constraints

- Preserve unrelated dirty files and stage only listed paths.
- Use an alias-free periodic `L=6` contraction torus and exact nine-cell normalization.
- Require a degree-five self-check against the frozen quadratic witness artifact.
- Keep each local shard below ten minutes and 4 GB.
- Do not claim a full E7 map, an all-order remainder, or HPC authorization.

### Task 1: Exact contraction core

**Files:**
- Create: `src/trottercert/dual_log_pairing.py`
- Create: `tests/test_dual_log_pairing.py`

**Interfaces:**
- Produces `DualPairingPartial` and `contract_log_degree_shard(stages, degree, shard_index, shard_count, length=6)`.

- [x] Add failing tests for shard validation, deterministic group assignment,
  exact cubic addition, and degree-five agreement with the frozen E5 pairings.
- [x] Run `PYTHONPATH=src python -m pytest -q tests/test_dual_log_pairing.py`
  and confirm missing-module failure.
- [x] Implement suffix grouping, one-translation target lookup, cache eviction,
  exact denominator handling, and partial counters.
- [x] Run the focused test and require exact E5 equality.
- [x] Commit with `feat(issue128): add exact dual log contraction`.

### Task 2: Shard artifacts and reducer

**Files:**
- Create: `scripts/certify_dual_e7_pairing.py`
- Create: `tests/test_dual_e7_pairing.py`
- Create: `docs/experiments/processor-obstruction/dual-e7-shards/shard-*.json`
- Create: `docs/experiments/processor-obstruction/dual-e7-pairing.json`

**Interfaces:**
- Shard mode consumes `--shard-index`, `--shard-count`, and `--output`.
- Reduce mode consumes `--reduce`, four shard paths, and `--output`.

- [x] Add failing mutation tests for overlap, missing groups, configuration
  mismatch, source-digest mismatch, and forged cubic totals.
- [x] Implement canonical shard serialization and lightweight verification.
- [x] Implement reducer coverage checks and forward/reverse exact equality.
- [x] Run four local shards concurrently with `shard_count=4`, then reduce.
- [x] Independently rerun one shard and compare its SHA-256.
- [x] Commit with `feat(issue128): certify exact dual E7 pairing`.

### Task 3: Dual finite-step ledger

**Files:**
- Create: `docs/report/dual-e7-remainder-ledger.md`
- Modify: `scripts/audit_gauge_aware_obstruction.py`
- Modify: `tests/test_gauge_aware_obstruction_audit.py`
- Modify: `docs/experiments/processor-obstruction/gauge-aware-audit.json`

- [x] Compute exact E5/r4 plus E7/r6 dual pairing at `r=95,96,97`, including
  interval-certified signs and margins.
- [x] Record the centered-stage path bounds showing which `r` values satisfy
  the elementary `< pi` logarithm-branch criterion.
- [x] Bind the E7 artifact in the fail-closed audit while retaining
  `finite_step_status=inconclusive` until the E9+ tail is supplied.
- [x] Run focused tests and the complete non-slow suite.
- [x] Commit with `docs(issue128): add dual E7 remainder ledger`.
