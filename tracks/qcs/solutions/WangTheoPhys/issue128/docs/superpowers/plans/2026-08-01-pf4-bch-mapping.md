# Exact PF4 BCH-to-Trace Mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Derive the true degree-five trace obstruction of the actual two-fragment five-copy Suzuki formula, independently verify it, and decide whether it yields a gauge-aware periodic-TFIM family theorem.

**Architecture:** Two exact implementations expand the same explicit 11-stage formula through degree five and reduce `Tr((A+B)L5)` in the cyclic free-word quotient without presupposing coefficients.  The derived identity replaces the legacy hardcoded form only after coefficient maps agree exactly; a separate commutant-witness layer decides fixed-target versus gauge-aware scope, and finite-step promotion remains independent.

**Tech Stack:** Python 3.12, exact `Fraction`, `Q(alpha)` with `alpha^3=4`, free associative words, cyclic trace reduction, exact Pauli algebra, canonical JSON, SHA-256, pytest.

## Global Constraints

- Formal formula: `S2(s)=exp(sA/2) exp(sB) exp(sA/2)`.
- Fourth-order formula: `S4(t)=S2(u*t)^2 S2(v*t) S2(u*t)^2`.
- Algebraic coefficients: `alpha^3=4`, `u=1/(4-alpha)`, `v=1-4u=-alpha*u`.
- The explicit merged product has 11 stages and is read left to right.
- Formal logarithm convention: `log S4(t)=t(A+B)+t^5*L5+O(t^7)`.
- The physical `A,B -> -iA,-iB` substitution and the sign relating `L5` to the Hermitian defect `E5` must be explicit.
- Define `C=[A,[A,B]]` and `D=[B,[B,A]]` exactly with this commutator orientation.
- Use either `Tr` throughout or normalized `tau=Tr/d` throughout each identity; record which one is used.
- Solve the cyclic coefficients without assuming the historical candidate `gamma*(1,-4,8/3)` or the legacy vector `(1/2,14/3,4/3)`.
- The historical candidate is a hypothesis; the legacy vector is a negative regression.
- TFIM family: `A=h*sum_i X_i`, `B=j*sum_i Z_i Z_(i+1)`, periodic even `L>=4`, exact rational `h,j`.
- A strict TFIM statement requires `h*j != 0`; the zero axes are explicit exceptions.
- Pairing with `H` proves only a fixed-target endpoint-conjugation obstruction because retiming can absorb the `H` direction.
- Gauge-aware promotion requires an exact commuting witness orthogonal to both `I` and `H`.
- Leading-order mapping never promotes a finite-step claim; the E7/E9/tail signed-margin gate remains separate.
- Do not modify the existing conditional TFIM artifact until both exact implementations pass.

---

## File Structure

### New files

- `src/trottercert/pf4_bch_mapping.py`: explicit 11-stage spec, cyclic free-word reduction, exact coefficient solving, physical-sign bridge, and canonical identity record.
- `scripts/derive_pf4_bch_mapping.py`: primary artifact builder.
- `scripts/reference_pf4_bch_mapping.py`: independent standard-library expansion with a handwritten stage table.
- `tests/test_pf4_bch_mapping.py`: order conditions, coefficient solution, conventions, matrices, and mutation tests.
- `tests/test_reference_pf4_bch_mapping.py`: subprocess cross-implementation tests.
- `docs/experiments/processor-obstruction/pf4-bch-mapping.json`: canonical exact mapping artifact.
- `src/trottercert/tfim_gauge_witness.py`: exact centered-polynomial commutant witnesses and TFIM trace-pairing records.
- `scripts/certify_tfim_gauge_witness.py`: canonical gauge-witness builder and primary semantic verifier.
- `scripts/reference_tfim_gauge_witness.py`: independent standard-library verifier for the gauge artifact.
- `tests/test_tfim_gauge_witness.py`: orthogonality, pairing, family-domain, and counterexample tests.
- `tests/test_reference_tfim_gauge_witness.py`: duplicate-key, source-closure, semantic-mutation, and subprocess tests.
- `docs/experiments/processor-obstruction/tfim-gauge-witness.json`: certified witness family or explicit no-promotion result.

### Existing files to modify after mapping review

- `src/trottercert/trace_obstruction.py`: replace the legacy conditional factory with the artifact-backed derived identity while preserving a legacy negative record.
- `tests/test_trace_obstruction.py`: add exact free-trace and legacy-rejection tests.
- `src/trottercert/tfim_obstruction.py`: contract the TFIM moments with the general derived cubic triple `(g_c2,g_cd,g_d2)` and distinguish fixed-target/gauge-aware status; optional gamma metadata is used only after a proved projective match.
- `scripts/certify_tfim_obstruction.py`: bind the mapping and gauge-witness artifact hashes.
- `tests/test_tfim_obstruction.py`: enforce the new gate and reject rehashed legacy promotion.
- `docs/experiments/processor-obstruction/tfim-family-obstruction.json`: regenerate only after all semantic checks.
- `artifacts/publication/paper-a-file-ownership.json`: retain all new files as Paper B and outside Paper A release.

---

### Task 1: Freeze the Explicit 11-Stage Formula

**Files:**
- Create: `src/trottercert/pf4_bch_mapping.py`
- Create: `tests/test_pf4_bch_mapping.py`

**Interfaces:**
- Produces: `pf4_two_fragment_stages() -> tuple[CubicStage, ...]`.
- Produces: `PF4ConventionRecord` containing formal, physical, trace, and commutator conventions.

- [ ] **Step 1: Write the exact stage-table test**

```python
def test_two_fragment_pf4_has_the_frozen_eleven_stage_table() -> None:
    u = Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
    v = Cubic.one() - 4 * u
    expected = (
        (0, u / 2), (1, u), (0, u), (1, u),
        (0, (u + v) / 2), (1, v), (0, (u + v) / 2),
        (1, u), (0, u), (1, u), (0, u / 2),
    )
    assert tuple((s.fragment_index, s.coefficient) for s in pf4_two_fragment_stages()) == expected
```

- [ ] **Step 2: Run the missing-module failure**

Run: `python -m pytest -q tests/test_pf4_bch_mapping.py`

- [ ] **Step 3: Implement an explicit factory and convention record**

The factory contains the eleven records directly and cross-checks them against `fourth_order_suzuki_cubic_stages(2)` only in tests.  The record pins left-to-right multiplication, commutator orientation, unnormalized `Tr` for the free identity, and the physical `-i` substitution.

Prove the exact parity lemma for the palindromic formula: `S4(-t)=S4(t)^(-1)`, hence the formal logarithm is odd.  Test this lemma algebraically and either explicitly verify the degree-six log map is zero or derive it from the lemma before recording `log S4=tH+t^5L5+O(t^7)`.

- [ ] **Step 4: Add convention mutation tests**

Mutate one noncentral stage coefficient, change the central `v` stage, swap one adjacent unequal A/B pair, swap A/B globally, change `D` to `[B,[A,B]]`, or change the physical sign; require canonical verification to reject each record even after digest resealing.  The 11-stage table is palindromic, so reversal must be proved identical and is not a negative test.  If multiplication-direction semantics need a rejection fixture, use a separate nonpalindromic auxiliary product.

- [ ] **Step 5: Commit Task 1**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/pf4_bch_mapping.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_pf4_bch_mapping.py
git commit -m "feat(issue128): freeze two-fragment PF4 conventions"
```

### Task 2: Derive L5 and Solve the Cyclic Trace Coefficients

**Files:**
- Modify: `src/trottercert/pf4_bch_mapping.py`
- Modify: `tests/test_pf4_bch_mapping.py`

**Interfaces:**
- Produces: `derive_pf4_cyclic_mapping() -> PF4CyclicMapping`.
- `PF4CyclicMapping` contains `g_c2`, `g_cd`, `g_d2`, a zero residual map, exact order-condition maps, and a payload digest.

- [ ] **Step 1: Write exact order-condition tests**

```python
def test_pf4_log_has_the_declared_order() -> None:
    mapping = derive_pf4_cyclic_mapping()
    assert mapping.log_terms[1] == {(0,): Cubic.one(), (1,): Cubic.one()}
    assert mapping.log_terms[2] == {}
    assert mapping.log_terms[3] == {}
    assert mapping.log_terms[4] == {}
    assert mapping.log_terms[5]
```

- [ ] **Step 2: Implement cyclic canonicalization**

Map each nonempty word to the lexicographically smallest rotation and sum exact `Cubic` coefficients.  Reject empty degree-six words and remove exact zeros only after cubic-field reduction.

- [ ] **Step 3: Construct the basis without coefficient assumptions**

Build free-word polynomials for `H=A+B`, `C=[A,[A,B]]`, and `D=[B,[B,A]]`; calculate cyclic maps for `H*L5`, `C*C`, `C*D`, and `D*D`.

- [ ] **Step 4: Solve for three cubic coefficients**

Treat each cubic coefficient as three rational coordinates.  Assemble exact rational linear systems from all cyclic-word coordinates, solve by fraction-preserving Gaussian elimination, and require a unique solution.  Reconstruct the full cyclic map and require an empty residual dictionary.

- [ ] **Step 5: Test candidate and legacy vectors projectively after solving**

Always retain the general exact triple `(g_c2,g_cd,g_d2)`.  Compare the solution projectively with candidate `(1,-4,8/3)` and legacy `(1/2,14/3,4/3)`: select a nonzero component, solve the exact candidate-specific scale, and require all three scaled coordinates to agree.  Store separate `candidate_projective_status` and `legacy_projective_status`; neither vector is an input to the solver.  Generic mixed-term fixtures must expose a nonzero residual for any nonproportional hypothesis.

- [ ] **Step 6: Prove coefficient signs and an optional projective prefactor**

Use `cube_root_four_interval` with increasing decimal digits to classify each of `g_c2`, `g_cd`, and `g_d2` as `negative`, `positive`, or `exact_zero`; never require a nonzero sign for a coefficient that vanishes algebraically.  Only if the solved triple is projectively proportional to the historical candidate may the record additionally contain an optional `gamma` coordinate triple and sign enclosure.  Verification rechecks `alpha^3=4`, every coefficient enclosure/status, projective residual, and optional gamma; absence of a single gamma never blocks a valid general solution.

- [ ] **Step 7: Run focused tests and commit**

Run: `python -m pytest -q tests/test_pf4_bch_mapping.py`

Commit with the Task 2 allowlist:

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/pf4_bch_mapping.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_pf4_bch_mapping.py
git commit -m "feat(issue128): derive exact PF4 cyclic mapping"
```

### Task 3: Implement an Independent Handwritten Expansion

**Files:**
- Create: `scripts/reference_pf4_bch_mapping.py`
- Create: `tests/test_reference_pf4_bch_mapping.py`

**Interfaces:**
- CLI: `python scripts/reference_pf4_bch_mapping.py ARTIFACT`.
- The script implements its own cubic arithmetic and a recursive BCH/Hall-Lyndon expansion, then maps the result into cyclic words and solves the exact trace system.  It does not reuse the primary product-series/log algorithm.

- [ ] **Step 1: Write a source-independence test**

Parse the reference script AST and reject imports from `trottercert`, NumPy, SymPy, or the primary generator.  Permit only standard-library modules.

- [ ] **Step 2: Handwrite the eleven stages**

Represent a cubic number as a three-`Fraction` tuple and implement multiplication modulo `alpha^3-4`.  The stage table is literal data in the reference script, not imported or parsed from the artifact.

- [ ] **Step 3: Implement an algorithmically independent logarithm through degree five**

Compose the literal stages with a recursive BCH expansion in a fixed Hall-Lyndon basis through degree five, then expand the Lie basis into free words only at the end.  Verify degrees one through four before reading the submitted degree-five coefficients.  A test compares the primary free-word result and the independent BCH result; code/source separation alone is not accepted as algorithmic independence.

- [ ] **Step 4: Independently solve the cyclic system and the physical bridge**

Rebuild H, C, D from commutator definitions and solve the coefficient system.  Compare the complete canonical cyclic map, general coefficient triple, projective comparison records, residual, root intervals, and physical-sign theorem with the artifact.  Reconstruct by homogeneity that `L5(-iA,-iB)=(-i)^5 L5(A,B)=-iL5(A,B)`; with `log U=-itH-it^5E5+O(t^7)`, require `E5=L5`.  A submitted sign string is never trusted.

- [ ] **Step 5: Add adversarial subprocess tests**

Change one stage, one cubic coordinate, one word coefficient, one solved coefficient, residual status, root interval, physical sign, or identity status; reseal the digest and require nonzero exit.

- [ ] **Step 6: Run and commit**

Run: `python -m pytest -q tests/test_reference_pf4_bch_mapping.py`

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/reference_pf4_bch_mapping.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_reference_pf4_bch_mapping.py
git commit -m "feat(issue128): independently derive PF4 mapping"
```

### Task 4: Add Generic Exact-Matrix Convention Checks

**Files:**
- Modify: `tests/test_pf4_bch_mapping.py`
- Modify: `tests/test_reference_pf4_bch_mapping.py`

**Interfaces:**
- Tests only; matrix substitutions are falsification checks and never replace the free-word proof.

- [ ] **Step 1: Add seeded 2 by 2 and 3 by 3 Hermitian integer matrices**

Evaluate L5 directly from the solved free-word polynomial and compare `Tr((A+B)L5)` with the solved `C2/CD/D2` combination using exact object arrays.

- [ ] **Step 2: Add discriminating mixed-term fixtures**

Require `Tr(CD)!=0`; reject fixtures whose mixed term vanishes.  Confirm that the solved identity agrees.  Compare each candidate projectively: fit its scale from one exact nonzero component, then require all remaining components and the exact matrix residual to agree.  At least one mixed fixture must reject every nonproportional legacy hypothesis.

- [ ] **Step 3: Add zero and commuting limits**

Check `A=0`, `B=0`, and diagonal commuting A/B all yield zero.  Check A/B swap against its separately derived convention.  Prove reversal of the palindromic PF4 stage table gives the same product; do not use reversal as a negative fixture.

- [ ] **Step 4: Run and commit**

Run: `python -m pytest -q tests/test_pf4_bch_mapping.py tests/test_reference_pf4_bch_mapping.py`

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_pf4_bch_mapping.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_reference_pf4_bch_mapping.py
git commit -m "test(issue128): cross-check PF4 trace conventions"
```

### Task 5: Freeze the Mapping Artifact and Replace the Legacy Record

**Files:**
- Create: `scripts/derive_pf4_bch_mapping.py`
- Create: `docs/experiments/processor-obstruction/pf4-bch-mapping.json`
- Modify: `src/trottercert/trace_obstruction.py`
- Modify: `tests/test_trace_obstruction.py`

**Interfaces:**
- Produces a canonical compact JSON artifact with source hashes and both implementation results.
- Produces: `pf4_derived_trace_identity_record()`.
- Preserves: `pf4_legacy_hypothesis_record()` with non-promoted status.
- Primary CLI uses exclusive `--build --output PATH` and `--verify PATH` modes; the reference CLI uses `--verify PATH`.

- [ ] **Step 1: Implement canonical artifact generation**

Store convention record, 11-stage map, degree-one through degree-five hashes, complete cyclic coefficient map, the general solved cubic coefficients `(g_c2,g_cd,g_d2)`, each coefficient enclosure/status, residual, candidate and legacy projective-comparison records, optional gamma only when proportionality is proved, physical bridge, primary/reference implementation hashes, and payload digest.  Include no timings.  Both CLIs reject duplicate JSON keys, floats, noncanonical cubic/rational encodings, unknown fields, noncanonical bytes, and source-hash drift.

The source closure roots at `src/trottercert/pf4_bch_mapping.py`, `scripts/derive_pf4_bch_mapping.py`, and `scripts/reference_pf4_bch_mapping.py`, and explicitly includes `src/trottercert/cubic_field.py`, `src/trottercert/intervals.py`, `src/trottercert/formulas.py`, `src/trottercert/series.py`, `src/trottercert/lie_series.py`, `src/trottercert/algebra.py`, `pyproject.toml`, and `requirements-reproducibility.txt`.  Recursively AST-walk local imports and reject a transitive path missing from the canonical manifest.  Resolve normalized relative paths within the issue root only, reject traversal/symlinks outside the root, hash every file, and require the reference verifier to validate the same closure before mathematical replay.

- [ ] **Step 2: Build and independently verify after setup confirmation**

Run:

```bash
PYTHONPATH=src python scripts/derive_pf4_bch_mapping.py --build \
  --output docs/experiments/processor-obstruction/pf4-bch-mapping.json
PYTHONPATH=src python scripts/derive_pf4_bch_mapping.py --verify \
  docs/experiments/processor-obstruction/pf4-bch-mapping.json
python scripts/reference_pf4_bch_mapping.py --verify \
  docs/experiments/processor-obstruction/pf4-bch-mapping.json
```

- [ ] **Step 3: Update record factories**

The derived factory reads or reconstructs only the canonical solved identity.  `verify_identity_record` must reject a legacy record relabeled as proved, a derived record with changed coefficients, and a digest-resealed promotion.

- [ ] **Step 4: Rebuild byte-identically and review**

Regenerate to a temporary file from the same immutable source commit, compare canonical bytes and SHA-256, run both implementations, and obtain an independent scientific review of stage order, solution uniqueness, residual zero, all three coefficient signs, optional projective prefactor, and physical bridge.

- [ ] **Step 5: Commit Task 5**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/derive_pf4_bch_mapping.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/pf4-bch-mapping.json \
  tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/trace_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_trace_obstruction.py
git commit -m "feat(issue128): certify exact PF4 BCH mapping"
```

### Task 6: Substitute the Periodic TFIM Family Exactly

**Files:**
- Modify: `src/trottercert/tfim_obstruction.py`
- Modify: `scripts/certify_tfim_obstruction.py`
- Modify: `tests/test_tfim_obstruction.py`
- Modify: `docs/experiments/processor-obstruction/tfim-family-obstruction.json`

**Interfaces:**
- Consumes the exact solved mapping artifact.
- Produces fixed-target normalized trace pairing and conditional operator lower bound from the general solved cubic coefficients.

- [ ] **Step 1: Write family substitution tests**

For even `L=4,6,8`, exact rational `h,j`, compute `g_c2*tau(C^2)+g_cd*tau(CD)+g_d2*tau(D^2)` directly; require zero on `h=0` or `j=0` and a strict outward nonzero interval when `h*j!=0` on the proved parameter domain.  TFIM may have `tau(CD)=0`, so it validates substitution but never distinguishes the mixed coefficient.

- [ ] **Step 2: Replace `gamma=1` conditional values**

Store the exact `Q(alpha)` pairing coordinates and first compute an outward interval `[p_lower,p_upper]`.  Derive a strict lower bound for `|tau(H E5)|` only when zero is excluded.  With normalized `tau`, use Hölder and the certified upper bound `M_H=L(|h|+|j|)>=||H||_infinity` to report `||E5||_infinity >= lower_abs(tau(H E5))/M_H`.  Handle `H=0` separately; on `h=0` or `j=0`, report the exact zero/exception state rather than divide.  Do not serialize a rational `operator_lower_bound` unless both the cubic absolute lower bound and the denominator direction have been certified outward.

- [ ] **Step 3: Bind the all-even counting argument**

Emit separate local Pauli-orbit multiplicity formulas for the wraparound case `L=4` and for generic even `L>=6`; do not infer an all-even formula from samples.  The verifier checks both closed families and separately evaluates L=4,6,8 as regression cases.

- [ ] **Step 4: Keep fixed-target scope explicit**

Set `fixed_target_leading_order_status` from the mapping.  Leave `gauge_aware_status=blocked_pending_orthogonal_commutant_witness` and `finite_step_no_go=not_claimed`.

- [ ] **Step 5: Run, rebuild, review, and commit**

Run: `python -m pytest -q tests/test_tfim_obstruction.py tests/test_trace_obstruction.py tests/test_pf4_bch_mapping.py`

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/certify_tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/tfim-family-obstruction.json
git commit -m "feat(issue128): bind TFIM family to derived PF4 mapping"
```

### Task 7: Search for and Certify a Gauge-Orthogonal Commutant Witness

**Files:**
- Create: `src/trottercert/tfim_gauge_witness.py`
- Create: `scripts/certify_tfim_gauge_witness.py`
- Create: `scripts/reference_tfim_gauge_witness.py`
- Create: `tests/test_tfim_gauge_witness.py`
- Create: `tests/test_reference_tfim_gauge_witness.py`
- Create: `docs/experiments/processor-obstruction/tfim-gauge-witness.json`

**Interfaces:**
- Produces: `centered_power_witness_moments(length,h,j,power) -> WitnessRecord` for powers 2, 3, and 4.
- Produces an exact promoted witness or an explicit `no_witness_in_tested_polynomial_span` result.
- Builder CLI uses exclusive `--build --output PATH` and `--verify PATH`; the independent standard-library CLI uses `--verify PATH`.

- [ ] **Step 1: Define exact projected polynomial witnesses**

For `K=H^power`, solve the exact two-by-two Gram system for `a,b` in `W=K-aI-bH` so that normalized `tau(W)=tau(W*H)=0`.  Reject singular Gram systems unless direct exact equations establish the projection.

- [ ] **Step 2: Write commutant and orthogonality tests**

Using exact Pauli sums for L=4 and L=6, require `[W,H]=0`, `tau(W)=0`, `tau(W*H)=0`, and exact normalization metadata.  Mutating `a`, `b`, h, j, or boundary must fail.

- [ ] **Step 3: Derive `tau(W*E5)` from free traces**

Extend the cyclic word calculation to `H^power*L5`, substitute exact TFIM Pauli orbit counts, and calculate the pairing without dense matrices.  Search powers in the frozen order 2,3,4; do not choose the order after seeing numerical magnitudes.

- [ ] **Step 4: Freeze success or failure honestly**

Promote the first power with an exact nonzero pairing over a stated nonempty parameter domain.  This proves qualitative gauge-aware non-membership.  Report a quantitative operator-norm lower bound only after independently certifying the normalized trace norm `tau(|W|)` (or an outward upper bound) and applying the exact duality denominator; commutation plus `W` orthogonal to `I,H` alone is not a quantitative norm certificate.  If all three powers vanish or change sign without a provable domain, freeze the negative result and retain fixed-target-only scope.

- [ ] **Step 5: Independently verify the witness artifact**

The builder emits canonical compact JSON with complete source closure and no timings.  Both verifiers reject duplicate keys, unknown fields, noncanonical rationals/bytes, source drift, and rehashed semantic mutations.  The reference verifier reconstructs the separate `L=4` and generic even `L>=6` Gram/orbit formulas, commutation, both orthogonality equations, cyclic pairing, parameter domain, and any trace-norm denominator used quantitatively.  Rehashed fake nonzero pairings must fail.

Build, verify, and byte-rebuild after setup confirmation:

```bash
PYTHONPATH=src python scripts/certify_tfim_gauge_witness.py --build \
  --output docs/experiments/processor-obstruction/tfim-gauge-witness.json
PYTHONPATH=src python scripts/certify_tfim_gauge_witness.py --verify \
  docs/experiments/processor-obstruction/tfim-gauge-witness.json
python scripts/reference_tfim_gauge_witness.py --verify \
  docs/experiments/processor-obstruction/tfim-gauge-witness.json
tmpdir="$(mktemp -d)"
PYTHONPATH=src python scripts/certify_tfim_gauge_witness.py --build \
  --output "$tmpdir/tfim-gauge-witness.json"
cmp docs/experiments/processor-obstruction/tfim-gauge-witness.json \
  "$tmpdir/tfim-gauge-witness.json"
shasum -a 256 docs/experiments/processor-obstruction/tfim-gauge-witness.json \
  "$tmpdir/tfim-gauge-witness.json"
```

- [ ] **Step 6: Commit Task 7**

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/src/trottercert/tfim_gauge_witness.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/certify_tfim_gauge_witness.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/scripts/reference_tfim_gauge_witness.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_tfim_gauge_witness.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_reference_tfim_gauge_witness.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/tfim-gauge-witness.json
git commit -m "feat(issue128): gate TFIM gauge-aware PF4 witness"
```

### Task 8: Enforce Leading-Order and Finite-Step Separation

**Files:**
- Modify: `scripts/certify_tfim_obstruction.py`
- Modify: `tests/test_tfim_obstruction.py`
- Modify: `docs/experiments/processor-obstruction/README.md`
- Modify: `artifacts/publication/paper-a-file-ownership.json`

- [ ] **Step 1: Add promotion-state tests**

Mapping-only evidence may set fixed-target leading-order proved.  Gauge-aware status requires the orthogonal witness artifact.  For the current two-fragment periodic TFIM, set `finite_step_no_go=not_claimed`: the existing finite-step artifact is a 31-stage four-matching 12 by 12 Heisenberg certificate and must be classified as incompatible.  A future finite-step status requires a separately verified positive signed-margin artifact and cannot be changed by either PF4 file.

- [ ] **Step 2: Add cross-artifact source closure**

Bind mapping, TFIM moments, gauge witness, and corrector-cost definitions by exact hashes while preserving their distinct statuses.  A future finite-step artifact is accepted only if it simultaneously matches the PF4 formula identifier and 11-stage table hash; the two-fragment Hamiltonian `A=h sum X`, `B=j sum ZZ`; periodic `L,h,j`; trace normalization; target witness; step size; and same-formula E7/E9/tail source closure.  Add a cross-formula mutation test that proves the current Heisenberg artifact is rejected rather than included as source closure.

- [ ] **Step 3: Update Paper B ownership only**

Classify all PF4 and gauge-witness files as Paper B with `release_included=false` for Paper A.  Do not add them to Paper A claims or release closure.

- [ ] **Step 4: Run the Paper B focused regression**

Run:

```bash
python -m pytest -q tests/test_pf4_bch_mapping.py \
  tests/test_reference_pf4_bch_mapping.py tests/test_trace_obstruction.py \
  tests/test_tfim_obstruction.py tests/test_tfim_gauge_witness.py \
  tests/test_corrector_cost.py tests/test_finite_step_obstruction.py
PYTHONPATH=src python scripts/derive_pf4_bch_mapping.py --verify \
  docs/experiments/processor-obstruction/pf4-bch-mapping.json
python scripts/reference_pf4_bch_mapping.py --verify \
  docs/experiments/processor-obstruction/pf4-bch-mapping.json
PYTHONPATH=src python scripts/certify_tfim_obstruction.py --verify \
  docs/experiments/processor-obstruction/tfim-family-obstruction.json
PYTHONPATH=src python scripts/certify_tfim_gauge_witness.py --verify \
  docs/experiments/processor-obstruction/tfim-gauge-witness.json
python scripts/reference_tfim_gauge_witness.py --verify \
  docs/experiments/processor-obstruction/tfim-gauge-witness.json
tmpdir="$(mktemp -d)"
PYTHONPATH=src python scripts/derive_pf4_bch_mapping.py --build \
  --output "$tmpdir/pf4-bch-mapping.json"
cmp docs/experiments/processor-obstruction/pf4-bch-mapping.json \
  "$tmpdir/pf4-bch-mapping.json"
shasum -a 256 docs/experiments/processor-obstruction/pf4-bch-mapping.json \
  "$tmpdir/pf4-bch-mapping.json"
```

- [ ] **Step 5: Independent two-stage review and commit**

One reviewer checks the exact free algebra and physical convention; a different reviewer checks claim scope, gauge orthogonality, finite-step separation, and ownership.  Commit only after both report no P1/P2.

```bash
git add tracks/qcs/solutions/WangTheoPhys/issue128/scripts/certify_tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/tests/test_tfim_obstruction.py \
  tracks/qcs/solutions/WangTheoPhys/issue128/docs/experiments/processor-obstruction/README.md \
  tracks/qcs/solutions/WangTheoPhys/issue128/artifacts/publication/paper-a-file-ownership.json
git commit -m "chore(issue128): gate PF4 family claims"
```

## Self-Review Checklist

- [ ] The coefficient solver has no expected vector as input.
- [ ] The old and historical candidate forms are hypotheses, not fixtures that define success.
- [ ] Generic mixed-term matrices distinguish the formulas; TFIM alone is not used for that purpose.
- [ ] The physical `-i` sign bridge is explicit and tested.
- [ ] Exact coefficient agreement is map-level, not digest-only.
- [ ] A pairing with H is never labeled gauge-aware.
- [ ] A gauge witness is exactly orthogonal to I and H.
- [ ] Mapping, gauge, and finite-step statuses remain separate.
- [ ] Failure to find a degree-2/3/4 polynomial witness is frozen as a negative result rather than hidden.
