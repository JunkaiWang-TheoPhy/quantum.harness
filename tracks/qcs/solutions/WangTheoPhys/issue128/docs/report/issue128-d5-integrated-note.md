# Issue 128 D5-integrated certificate: technical note

Date: 2026-07-31 (Asia/Shanghai)

## Certified result

For the unchanged periodic `12×12` isotropic Heisenberg benchmark at `T=1`
and global operator-norm tolerance `1e-6`, the post-freeze certificate
`certificates/issue128-d5-integrated-certificate.json` accepts 95 steps and
rejects 94. The five-copy fourth-order Suzuki circuit and compilation model
are unchanged.

| Resource | Published control | D5-integrated certificate |
|---|---:|---:|
| Trotter steps | 393 | 95 |
| Merged group exponentials | 11,791 | 2,851 |
| Bond propagators | 848,952 | 205,272 |
| CNOT upper bound | 2,546,856 | 615,816 |

The exact primary ratio is

```text
11791/2851 = 4.13574184496667835847...
```

This saves 643,680 bond propagators and 1,931,040 CNOTs relative to the pinned
control under the same upper-bound model.

## What changed

The original frozen schema-v3 certificate used an anticommuting D4 sidecar and
a generic D5 majorant. The new certificate keeps the same D4 proof and replaces
only the D5 majorant by the already generated exact grouped D5 proof:

- 605,832 canonical Pauli coefficients;
- 123,106 same-support pairwise-anticommuting groups;
- maximum group size 10;
- exact site-norm upper bound
  `44948270001027856175670154896253 / 4000000000000000000000000000000`.

The D5 sidecar is content-bound by SHA-256. The verifier reconstructs its
canonical gzip payload, checks coverage and support constraints, checks every
pairwise symplectic anticommutation relation, recomputes every outward group
bound, and verifies all metadata against the main certificate. The deep mode
then regenerates all 605,832 exact D5 coefficients from the 31-stage Suzuki
formula and requires exact interval equality with the sidecar.

## Exact finite-step ledger

At `r=95`, the outward global contributions are:

| Right-generator degree | Global upper bound |
|---:|---:|
| D4 | `5.7218800760132647e-7` |
| D5 grouped | `3.4853484615103630e-8` |
| D6 | `1.7893140877453180e-7` |
| D7 | `2.1921361303553973e-8` |
| D8 and above | `1.8255487798872918e-7` |
| **Total** | **`9.9044914028324506e-7`** |

The same exact-rational error function gives
`1.0468061165603706e-6` at `r=94`, so the adjacent integer boundary is closed.
Fast verification now regenerates every displayed contribution and the
adjacent-step value. In particular, setting D6, D7, or the tail to zero and
synchronously recomputing the submitted total is rejected.

## Verification evidence

Run the fast verifier:

```bash
PYTHONPATH=src python3 scripts/verify.py \
  certificates/issue128-d5-integrated-certificate.json
```

Regenerate the proof algebra:

```bash
PYTHONPATH=src python3 scripts/verify.py --deep \
  certificates/issue128-d5-integrated-certificate.json
```

Run the standard-library-only downstream checker:

```bash
python3 scripts/reference_verify.py \
  certificates/issue128-d5-integrated-certificate.json
```

Observed on 2026-07-31:

- fast verifier: `valid=true`, `verification_level=fast`,
  `finite_step_bound_recomputed=true`;
- deep verifier: `valid=true`, `deep_proof_regenerated=true`,
  `baseline_centers_scanned=31`, 408.88 s, maximum RSS 1,538,768,896 bytes;
- independent downstream checker: `valid=true`, 3.63 s, with no
  `trottercert` imports;
- normal test suite: 146 passed, 12 explicitly deselected slow tests;
- focused certificate-mutation suite: 10 passed;
- independent-checker suite: 8 passed.

The independent checker reimplements D4/D5 partition coverage, symplectic
anticommutation, square-root domination, the finite-step ledger, Suzuki tail,
95/94 comparison, resource arithmetic, and strict equality of the declared
D4/D5 term counts, group counts, and maximum group sizes using only the Python
standard library. It does not regenerate the D4/D5 coefficient maps; that
distinct obligation is discharged by deep mode.

A nondegenerate open `2×3` dense diagnostic contains seven unique bonds. At 95
steps its measured operator-norm error is
`1.8313133889571672e-10`, below the outward scaled periodic comparison
`4.129166666666667e-8`. This checks normalization and stage ordering but is not
presented as a proof for open boundaries.

The builder is deterministic from the frozen certificate and D5 sidecar:

```bash
PYTHONPATH=src python3 scripts/build_d5_integrated_certificate.py
```

## Provenance and claim boundary

The original ten-file delivery manifest remains frozen at 97 steps and
`11791/2911`. It is not silently overwritten. This post-freeze certificate is
a new, separately verifiable artifact that reuses the hash-bound D4/D5
sidecars and extends the verifier. A reviewer can therefore reproduce either
claim without ambiguity.

No fivefold claim is made. Fivefold resource arithmetic requires `r≤78`, but
the current D4 contribution alone evaluates to approximately `1.259084e-6` at
78 steps, already above the total tolerance. Exact-D8 HPC work remains useful
for diagnosis and future proof design, but cannot close fivefold without a
new D4 tightening.

The certificate is scoped to the stated Hamiltonian, normalization, periodic
boundary, formula, tolerance, and cost model. It does not prove that 94
physical steps fail, that the grouping heuristic is optimal, or that the
method is state-of-the-art over unrestricted quantum simulation algorithms.

Release acceptance requires all of the following:

- the 95-step primary and reference verifiers pass;
- deep mode reports 31 scanned centers and coefficient regeneration;
- all ordinary and mutation tests pass;
- the D5-integrated SHA-256 manifest is clean;
- the PDF is rebuilt from this TeX source and visually inspected.
