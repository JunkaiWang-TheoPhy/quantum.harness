# Exact D6 Harvest Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Execute this plan inline, task by task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote a complete E7 fanout reduction into an independently cross-checked exact-D6 certificate and report the resulting adjacent integer resource boundary.

**Architecture:** Reuse the completed 272-, 544-, and 1,088-way E7 producer lanes. Accept the first complete reducer only as a candidate, require a second differently partitioned reducer to produce the identical canonical D6 coefficient map, then bind the selected payload and parent manifest into a new certificate and run fast/deep verification plus focused corruption tests.

**Tech Stack:** Python 3.11, exact `Fraction` and cubic-field arithmetic, deterministic gzip/JSON proof artifacts, Slurm result manifests, pytest.

## Global Constraints

- Benchmark remains the periodic 12×12 isotropic Heisenberg model with N=144, T=1, and operator-norm tolerance 10⁻⁶.
- Formula and resource model remain the five-copy fourth-order Suzuki kernel with G(r)=30r+1.
- Preserve the certified 95-step / 4.1357418449× certificate; create separate D6 artifacts.
- Do not use incomplete shards, partial sums, floating-point discovery values, processor results, or 4×4-patch work.
- Record producer, reducer, and verifier provenance separately.
- Do not modify the report, presentation, or D5 transcript files already dirty in the shared checkout.

---

### Task 1: Harvest a complete reducer

**Files:**
- Read remotely: reducer summaries, parent manifests, and exact-D6 payloads for jobs `23044907`, `23044906`, and `23044816`
- Create locally when valid: `certificates/issue128-d6-exact.json.gz`
- Create locally when valid: `certificates/issue128-d6-parent.json`

**Interfaces:**
- Consumes: complete E7 word-shard directories from one fanout width.
- Produces: one canonical degree-six payload plus its parent provenance document.

- [ ] Check reducer job state and exit code; accept only `COMPLETED` with exit code `0:0`.
- [ ] Read the summary and require `status=complete` and `reduction_order_check=forward_equals_reverse`.
- [ ] Verify output and parent SHA-256 values before copying.
- [ ] Copy the payload and parent into the certificate directory without touching frozen D5 artifacts.

### Task 2: Cross-check a second fanout layout

**Files:**
- Read: first and second reducer outputs from different fanout widths
- Use: `scripts/compare_exact_degree_payloads.py`
- Test: `tests/test_compare_exact_degree_payloads.py`

**Interfaces:**
- Consumes: two complete exact degree-six payloads.
- Produces: an exact-map equality record; any coefficient mismatch aborts promotion.

- [ ] Run the focused comparison tests and require both agreement and corruption-rejection tests to pass.
- [ ] Compare the two remote/local payloads coefficient by coefficient, not only by term count or norm.
- [ ] Require degree, formula identifier, canonical Pauli keys, and exact cubic coefficients to match.
- [ ] Save the comparison record under `artifacts/d6-integrated/` only after exact equality.

### Task 3: Build the D6-integrated certificate

**Files:**
- Use: `scripts/build_d6_integrated_certificate.py`
- Create: `certificates/issue128-d6-integrated-certificate.json`
- Test: `tests/test_exact_series_certificate.py`
- Test: `tests/test_resources_verify.py`

**Interfaces:**
- Consumes: the selected exact-D6 payload, parent manifest, and frozen D5-integrated certificate.
- Produces: a separately named certificate with an automatically selected minimum passing step count.

- [ ] Verify the exact D6 site bound is no larger than the generic D6 majorant.
- [ ] Run the builder and record the candidate and previous-step exact rational error bounds.
- [ ] Recompute G(r)=30r+1 and the exact ratio 11791/G(r).
- [ ] Preserve the 95-step certificate unchanged.

### Task 4: Verify and package the result

**Files:**
- Create: `artifacts/d6-integrated/verification-transcript.txt`
- Use: `scripts/verify.py`
- Test: focused D6 binding, resource corruption, and exact payload corruption tests

**Interfaces:**
- Consumes: the D6-integrated certificate and its hash-bound artifacts.
- Produces: fast/deep verification evidence and a final claim boundary.

- [ ] Run fast verification and require `valid=true`.
- [ ] Run deep verification and require `deep_proof_regenerated=true`.
- [ ] Run focused corruption tests for coefficient, digest, parent, and resource-boundary mutations.
- [ ] Capture commands, stdout, stderr, exit codes, source commits, and artifact hashes in the transcript.
- [ ] Report the new multiplier only if every gate passes; otherwise retain 4.1357418449×.

## Self-Review

- Scope covers only E7 fanout harvesting, exact D6 comparison, certificate integration, and verification.
- No processor, patch, matching-order, or D8 restart appears in the plan.
- All promotion gates fail closed and preserve the conservative certificate.
