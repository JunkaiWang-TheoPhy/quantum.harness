# Quantum Harness Issue 128 certificate

This directory contains a proof-carrying resource analysis for the periodic
`12 × 12` spin-1/2 isotropic Heisenberg benchmark at `T = 1` and operator-norm
tolerance `10^-6`.

The strongest completed certificate in this directory is:

- published-theorem control: 393 fourth-order Suzuki steps;
- certified candidate: 95 steps;
- merged group exponentials: 11,791 versus 2,851;
- exact resource ratio: `11791/2851 = 4.135741844966678...`;
- certified upper bound at 95:
  `9.904491402832450555... × 10^-7`;
- the same bound at 94:
  `1.046806116560370569... × 10^-6`.

The 94-step statement is adjacent minimality of this emitted upper-bound
function. It is not a lower bound on the true physical simulation error and is
not a claim that every 94-step circuit fails.

## Install

Python 3.12 is the frozen verification environment:

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

For exact reproduction of the recorded local package versions:

```bash
python -m pip install -r requirements-reproducibility.txt
python -m pip install -e . --no-deps
```

## Verification layers

Run the primary finite-step verifier:

```bash
python scripts/verify.py \
  certificates/issue128-d5-integrated-certificate.json
```

Fast mode still verifies the complete D4--D7 and tail ledger, the 95/94
boundary, both D4/D5 sidecar partitions, and all resource identities. It does
not regenerate every D4/D5 coefficient from the Suzuki formula.

Run the full algebraic replay:

```bash
python scripts/verify.py --deep \
  certificates/issue128-d5-integrated-certificate.json
```

Deep mode additionally regenerates all 31 published-theorem centers, selects
center 20, reconstructs the D4 and D5 coefficient maps, and requires exact
interval equality with the sidecars. On the recorded Python 3.12 environment,
the current-source run completed in 408.88 seconds with maximum RSS
1,538,768,896 bytes.

Run the structurally independent downstream checker:

```bash
python scripts/reference_verify.py \
  certificates/issue128-d5-integrated-certificate.json
```

`reference_verify.py` imports no `trottercert` code. It independently parses
the sidecars, verifies exact coverage, pairwise anticommutation, equal-support
D5 groups, rational square-root enclosures, every finite-step contribution,
the tail, the adjacent integer, and resource arithmetic. It deliberately does
not claim to regenerate the coefficient maps; that is the primary deep
verifier's role.

## Tests and active falsification

```bash
python -m pytest -q
python -m pytest -q tests/test_certificate_mutations.py
python -m pytest -q tests/test_reference_verify.py
```

The mutation suite includes attacks that previously distinguished internal
consistency from mathematical verification:

- set D6, D7, or tail to zero and recompute the submitted total;
- forge the adjacent-step value;
- alter an E7 formula constant;
- alter resources or the exact ratio;
- drift D4/D5 term counts, group counts, or maximum group sizes away from the
  verified sidecars;
- use booleans or strings where exact integers are required;
- repeat the tail attack against the D5-integrated certificate.

Every mutation must produce a verifier failure.

## Nondegenerate dense diagnostic

The original frozen package includes a degenerate periodic `2 × 2` algebra
check. A separate open `2 × 3` diagnostic has seven unique bonds and no
periodic edge duplication:

```bash
PYTHONPATH=src python scripts/crosscheck_small.py --width 2 --height 3
```

The observed 95-step operator-norm error is approximately
`1.8313133889571672 × 10^-10`, below the outward scaled periodic comparison
`4.129166666666667 × 10^-8`. This is an adversarial normalization and stage
ordering check, not a proof that the periodic certificate automatically
applies to every open-boundary model.

## Evidence layout

- `certificates/issue128-d5-integrated-certificate.json`: 95-step main record;
- `certificates/issue128-d4-groups.json`: grouped D4 witness;
- `certificates/issue128-d5-groups.json.gz`: grouped exact D5 witness;
- `artifacts/d5-integrated/verification-transcript.txt`: frozen run summary;
- `artifacts/d5-integrated/SHA256SUMS`: integrity manifest;
- `docs/report/issue128-d5-integrated-note.md`: reviewer-facing derivation;
- `docs/report/output/pdf/issue128-d5-integrated-note.pdf`: rendered note;
- `docs/manuscript/output/pdf/issue128-proof-carrying-trotter-paper.pdf`:
  the 31-page 95-step long-form manuscript;
- `docs/manuscript/`: complete manuscript source, figures, and claim validators;
- `docs/report/issue128-hpc-sprint-report.tex`: post-freeze D5/D6/D8 report.

The 97-step ten-file package is retained unchanged as a historical public
freeze. It is described as the predecessor inside the current manuscript. The
95-step bundle is a separately versioned improvement and must not be
represented as if the old manifest had silently changed.

## Claim boundary

This work certifies a resource reduction for one fixed model, formula,
normalization, tolerance, boundary condition, and compilation cost model. It
does not by itself establish universal product-formula optimality, a lower
bound on the true error, or state of the art over unrestricted quantum
simulation algorithms.
