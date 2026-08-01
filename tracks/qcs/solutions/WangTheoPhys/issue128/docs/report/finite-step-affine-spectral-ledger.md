# Finite-step affine effective-spectrum obstruction ledger

## Outcome

For the frozen periodic `12 x 12` isotropic Heisenberg Hamiltonian, the
five-copy fourth-order Suzuki formula at step size `h=1/97` has a certified
nonzero dual local-log pairing after exact E5, E7, and E9 terms and a compatible
E11-and-higher tail are combined.  The exact lower endpoint of the signed
margin is strictly positive.  A separate centered-moment invariant then
promotes this to

```text
certified_affine_effective_spectral_obstruction
```

for the principal effective Hamiltonian.  This is the submission-level closure
of the earlier E7 and dual-tail ledgers.

## Fixed object and normalization

The instance is

```text
model: isotropic Heisenberg
Hamiltonian normalization: (XX + YY + ZZ)/4 per nearest-neighbor bond
lattice: 12 x 12 periodic square lattice
formula: five-copy fourth-order Suzuki, 31 merged stages
target time: 1
steps: 97
step size: 1/97
dual normalization: normalized trace per 2 x 2 cell
```

The commuting quadratic witness is

```text
W = H^2 - 54 I + H/2,
tau(W) = tau(WH) = 0,
[W,H] = 0.
```

The combined local-log series is

```text
h^4 q5 + h^6 q7 + h^8 q9 + R_{>=11}^{dual}.
```

All claim-bearing arithmetic below is stored as canonical exact rational pairs
or exact cubic-field coordinates.  Decimal values are diagnostics only.

## Complete degree-nine contraction

The E9 production manifest contains 64 deterministic shards and covers

| quantity | exact count |
|---|---:|
| suffix groups | 65,536 |
| degree-nine words | 262,140 |
| nonzero words | 196,608 |
| retained Pauli terms | 2,105,732,028 |

The manifest index SHA-256 is

```text
b61d2f96c117fd64f969d395293c00d96b8380ff2be400bd365554532b775e10
```

Forward and reverse reduction orders agree exactly.  The per-cell E9 pairing
in `Q(alpha)`, `alpha^3=4`, is

```text
q9 = -66319553467/248832000000000
     - (1446202237/10368000000000) alpha
     - (12151491869/82944000000000) alpha^2
   ~= -8.571081962296589e-4.
```

The reduced artifact is

```text
docs/experiments/processor-obstruction/dual-e9-pairing.json
file sha256:    796db9b2f2a3635176ea5f100abc5da6d07fac68c6b3f7cd5191f0d927e9e7b1
payload sha256: ea17c2504462d19a102590f5de5046ea24bbc7fb132bebc8e90322aa6b92e144
```

Its 64 parent records bind the exact worker-file hashes.  Deterministically
selected median-load shard 18 and maximum-load shard 36 were recomputed in an
independent Slurm array; their mathematical-payload SHA-256 values match the
production workers exactly.  The compact job, deployment, rerun, and audit hash
chain is stored in `dual-e9-provenance.json`.  The full production and rerun
trees remain in the ignored `results/hpc` archive.

The E9 artifact deliberately retains `finite_step_status = inconclusive`: it
is a coefficient-only artifact and does not own the higher-order tail.  The
separate tail artifact deliberately says that exact E9 is missing.  Only the
hash-bound combiner below is permitted to promote their joint result.

## Exact signed-margin gate

At `h=1/97`, the exact cubic interval for the E5+E7+E9 contribution is strictly
negative.  Its absolute interval and the compatible tail give

| quantity | outward decimal summary |
|---|---:|
| `|h^4 q5 + h^6 q7 + h^8 q9|` lower | `6.9236871012985846e-11` |
| E11-and-higher dual tail upper | `6.4950925459924388e-12` |
| signed margin lower | `6.2741778466993407e-11` |

The exact lower endpoint of the final row is a positive rational number stored
under `/intervals/signed_margin/lower`.  Therefore the local-log status is

```text
certified_local_log_obstruction.
```

## Nonlinear affine-spectrum promotion

A nonzero tangent-space pairing alone does not exclude a finite affine unitary
orbit.  The final gate therefore uses centered normalized moments

```text
C(A) = A - tau(A) I,
m_k(A) = tau(C(A)^k),
Phi_H(A) = 54^3 m_3(A)^2 - 27^2 m_2(A)^3.
```

For every

```text
A = a I + b U H U^dagger,
```

including either sign of `b`, `Phi_H(A)=0` exactly.  For the principal
effective Hamiltonian `A=H+D` at one `1/97` step, the certificate obtains

| gate term | exact-bound decimal summary |
|---|---:|
| effective-log defect bound | `1.4223209486507154e-6` |
| centered defect cap | `2.8446418973014308e-6` |
| linear invariant lower bound | `5.7617660431199422e-2` |
| nonlinear remainder upper bound | `3.5254820251395814e-2` |
| invariant margin | `2.2362840179803608e-2` |

The invariant margin is an exact strictly positive rational.  Hence the
principal effective Hamiltonian is not in the affine unitary orbit of `H`, and
the promoted status is

```text
certified_affine_effective_spectral_obstruction.
```

## Claim boundary

The proved statement is fixed-instance and fixed-step:

- it concerns the principal-branch effective Hamiltonian at `h=1/97`;
- it excludes `a I + b U H U^dagger`, covering global phase, energy/time
  scaling, and unitary conjugation;
- it uses the frozen `12 x 12` periodic Heisenberg instance and normalization.

It does not claim:

- a lower bound on the eigenphases of the total-time unitary modulo `2 pi`;
- a no-wrap or no-eigenvalue-relabeling theorem for 97 repeated steps;
- a finite-step theorem uniform in lattice size;
- a TFIM finite-step theorem; the PF4--TFIM result remains a separate
  leading-order family statement.

The total-time eigenphase statement is not a wording variant of the present
theorem.  It would require an additional spectral branch/no-relabeling gate and
must remain unclaimed until such a proof exists.

## HPC and independent-audit chain

The exact compute used the immutable source commit

```text
e16a4e21e0ff240eaea07fa52e4e26c4c99985f4
```

with these Slurm jobs:

| role | job | required terminal state |
|---|---:|---|
| smoke shard 2 | `23065469_2` | `COMPLETED`, `0:0` |
| production array | `23065497` | all 64 shards `COMPLETED`, `0:0` |
| reducer | `23065498` | `COMPLETED`, `0:0` |
| independent reruns 18 and 36 | `23065746` | both `COMPLETED`, `0:0` |

The independent records are

| shard | production file SHA-256 | rerun file SHA-256 | matching mathematical payload SHA-256 |
|---:|---|---|---|
| 18 | `d2487ee2483419ff50af37d70bddf251009c83347d8a5ebf94e5ba331e7a15e5` | `e4356275f80dd00ae2e4731e5672eaae3fc80fb8e3cb914b393e989d31f88fbc` | `5efe20b5cc174880a9d2f121545e8b55c2f1a31fdf891ba8f228a339abc36890` |
| 36 | `04058fffae7a48f54a3541dfe132f2d75c4584b518875acd1899d2227eebb442` | `fe873a12e9b27fdfe2d7dd22c65277091cbcf53c2a9db02c22dd0e8eb96a8e0c` | `991bbb51af3e383198fd371a8f2c017211a367cfb35860d309647aa3e8a894e0` |

Different whole-file hashes are expected because runtime metadata differs;
the independently verified mathematical payloads are identical.  The archived
full-run audit has SHA-256
`2335f206ac149c122606075084bcd13302ae33b046dc7cf4c33b777220d8ae3c`.

The reducer re-verifies every worker before publishing the reduced JSON.  The
post-run audit additionally verifies the manifest index and all manifests,
unique complete group coverage, every worker payload, the reduced parent hash
list, and the two independent mathematical-payload matches.

## Reproduction and verification

Let `$E9_RUN_ROOT` point to the archived production tree and
`$E9_RERUN_ROOT` to the archived independent-rerun tree.  From the Issue-128
directory:

```bash
export PYTHONPATH=src:.

python scripts/certify_dual_e9_pairing.py \
  --index "$E9_RUN_ROOT/manifests/index.json" \
  --verify docs/experiments/processor-obstruction/dual-e9-pairing.json

python scripts/reduce_dual_e9_run.py \
  --run-root "$E9_RUN_ROOT" \
  --rerun "$E9_RERUN_ROOT/shards/shard-018.json" \
  --rerun "$E9_RERUN_ROOT/shards/shard-036.json"

python scripts/certify_finite_step_obstruction.py \
  --verify docs/experiments/processor-obstruction/finite-step-obstruction.json \
  --e9 docs/experiments/processor-obstruction/dual-e9-pairing.json \
  --index "$E9_RUN_ROOT/manifests/index.json"
```

The combined certificate is

```text
docs/experiments/processor-obstruction/finite-step-obstruction.json
file sha256:    13805e2d560efce6051f142ef9b69ae53b1f718d37426fb6fbaa10bd3072a909
payload sha256: f972e5b31e8dfe459aa2250ef7d6bf19c57e93d75bbf4f7bd6578b908656e24e
```
