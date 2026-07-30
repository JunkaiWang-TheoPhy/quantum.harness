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

Observed on 2026-07-31:

- fast verifier: `valid=true`, `verification_level=fast`;
- deep verifier: `valid=true`, `deep_proof_regenerated=true`, 167.62 s,
  maximum RSS 1,533,149,184 bytes;
- normal test suite: 100 passed, 11 explicitly deselected slow tests;
- focused corruption/minimality suite: 3 passed.

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
