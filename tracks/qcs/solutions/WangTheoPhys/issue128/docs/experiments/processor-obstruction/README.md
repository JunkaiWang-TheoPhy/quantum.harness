# Exact processor obstruction audit

## Decision

**NO-GO for exact cancellation by an endpoint conjugation processor.**
**NO-GO for exact cancellation by every Pauli-support-at-most-four local
processor.**  Whether a longer-range processor can lower the certified
operator-norm bound enough while also repairing the high-order ledger is
**INCONCLUSIVE**.

This is an isolated experiment.  It does not modify the D4/D5/D6 sidecars,
the accepted certificate, or any Slurm artifact.

## Exact object

For the frozen 31-stage fourth-order formula,

```text
log S(t) = -i t H - i t^5 E5 + O(t^7),   D4 = 5 E5.
```

The script reconstructs the canonical two-by-two-cell density over the exact
cubic field `Q(alpha)`, `alpha^3 = 4`.  It finds 74,448 nonzero E5 Pauli
coefficients:

| Pauli support | Exact nonzero count |
|---:|---:|
| 2 | 144 |
| 4 | 8,448 |
| 6 | 65,856 |

No floating-point optimizer enters these statements.

## Obstruction 1: every endpoint processor

Let `P(t)=exp(-i t^4 Q)`.  Conjugating the kernel changes the leading defect
by a commutator, up to the irrelevant sign convention:

```text
E5 -> E5 + i[Q,H].
```

For the Hilbert--Schmidt inner product,

```text
<H, i[Q,H]> = 0
```

for every Hermitian `Q`, by cyclicity of trace.  The exact density
reconstruction instead gives

```text
<H,E5>_cell =
  3589/2400000
  + (30361/38400000) alpha
  + (2813/4800000) alpha^2
  in [0.004227229425845701..., 0.004227229425845703...].
```

The three rational coordinates are not all zero, so this is a nonzero element
of `Q(alpha)`.  Therefore `E5` is not in `image(ad_H)`: no endpoint
conjugation processor can cancel the complete leading defect exactly.

As a diagnostic (not an operator-norm certificate), the invariant projection
onto `H` forces at least `0.1183044` relative residual in the canonical Pauli
coefficient `L2` metric.

## Obstruction 2: short-support processors

Commutation with the two-local Heisenberg Hamiltonian increases Pauli support
by at most one: a disjoint bond commutes, while an overlapping bond adds at
most its other endpoint.  Consequently, if every Pauli term of `Q` has
support at most four, every term of `[Q,H]` has support at most five.

The 65,856 exact support-six coefficients of E5 are therefore invariant under
all such processors.  They carry

```text
0.295176947961207... of ||E5||_Pauli,L2^2,
```

so this restricted class has a relative Pauli-`L2` residual of at least
`sqrt(0.295176947961207...) = 0.543301894`.  The existing six-dimensional
color-chirality basis is a subset of this rejected class because its
generators have support three.

This `L2` result diagnoses the algebraic obstruction.  It is not substituted
for the grouped operator-norm D4 sidecar.

## Why this does not open fivefold

At `r=78`, the current D5-integrated ledger contains

| contribution | upper bound |
|---|---:|
| D4 | `1.259084e-6` |
| D5 | `9.34095e-8` |
| D6 | `5.84064e-7` |
| D7 | `8.71506e-8` |
| D8+ tail | `9.00764e-7` |

Even setting both D4 and D6 to zero leaves

```text
D5 + D7 + tail = 1.081324e-6 > 1e-6.
```

Thus a processor that merely reduces D4 cannot certify fivefold.  A viable
processed formula would have to change the high-order structure, explicitly
rebuild D5 and beyond, postpone the geometric tail, and pass a new exact
global ledger.

There is also a circuit-cost gate.  Fivefold permits at most
`floor(11791/5)=2358` merged exponentials.  The unprocessed `r=78` kernel uses
`30*78+1=2341`, leaving at most **17** additional merged exponentials for both
endpoint processors.  Any proposal exceeding that cost cannot retain 5x even
if its error bound closes.

## Reproduction

From `tracks/qcs/solutions/WangTheoPhys/issue128`:

```bash
PYTHONPATH=src:. python -u scripts/experimental_processor_obstruction.py \
  --output docs/experiments/processor-obstruction/exact-obstruction.json

PYTHONPATH=src:. pytest -q tests/test_processors.py tests/test_cubic_local.py \
  -m 'not slow'
```

The exact reconstruction takes about two minutes on the current local host.
The focused regression result observed for this audit was `3 passed, 1
deselected`.

## Next decision

Do not enlarge the present chirality basis: its support class is already
ruled out.  If processor work continues, the minimum credible candidate must
contain support-five terms, cost no more than 17 endpoint exponentials after
merging, and be assessed as a new processed formula through at least explicit
D8 with a delayed rigorous tail.  Given the exact universal obstruction,
the more promising route is a genuinely different composition or phase
cycle, not a claim of exact order lifting by endpoint conjugation.

## Calibration-aware refinement

The exact `<H,E5>` witness proves the endpoint no-go at fixed target time but
is absorbed by the allowed time/energy-calibration direction.  The independent
quadratic commutant witness

```text
W = H^2 - 54 I + H/2
```

closes that gap on the periodic 12-by-12 lattice.  Exact Pauli orthogonality
gives

```text
tau(H^2) = 54,        tau(H^3) = -27,
tau(W) = 0,           tau(WH) = 0,          [W,H] = 0,
tau(W E5) =
  -7807/100000
  - (66043/1600000) alpha
  - (6119/200000) alpha^2 != 0,
alpha^3 = 4.
```

Therefore, after quotienting by

```text
image(i ad_H) + span(I) + span(H),
```

the leading defect still has a nonzero spectral component.  The calibrated
**leading-order** status is now **NO-GO**.  This is not yet a finite-step
eigenphase lower bound at total time one.  The exact dual E7 refinement below
certifies a local-log branch at `r=97`; promotion there now requires an
E9-and-higher dual tail small enough to preserve the exact E5+E7 pairing.

Likewise, the support-six result remains a rigorous no-go for the stated
short-support processor class in operator space.  Support alone does not
lower-bound eigenphase error, so that result is not promoted to a spectral
no-go.

The fail-closed refinement is recorded in:

```text
quadratic-commutant-witness.json
sha256: d0892b7320cb22a8be1701a1fdb185d2ce9dbfbc20d191fc21bf4359cea5d391
gauge-aware-audit.json
sha256: 565df344d324f410df8f8d27df5e0527f000486cda0a9ba47997c1d8cea3ad03
source exact-obstruction sha256:
b993596dcacb714c20bbae7b3e38c254e639268e57ab529823079db732f61103
```

Reproduce and verify it from the Issue-128 directory:

```bash
PYTHONPATH=src python -u -m scripts.certify_quadratic_commutant_witness

PYTHONPATH=src python -u -m scripts.certify_quadratic_commutant_witness \
  --verify

PYTHONPATH=src python -u -m scripts.audit_gauge_aware_obstruction \
  --source docs/experiments/processor-obstruction/exact-obstruction.json \
  --output docs/experiments/processor-obstruction/gauge-aware-audit.json

PYTHONPATH=src python -u -m scripts.audit_gauge_aware_obstruction \
  --source docs/experiments/processor-obstruction/exact-obstruction.json \
  --output docs/experiments/processor-obstruction/gauge-aware-audit.json \
  --verify
```

The full quadratic-witness recomputation takes about four minutes on the
current local host and peaks below 200 MB.  The audit also contains three exact
reference problems: an off-diagonal defect removed by a commutator, an
`H`-parallel defect removed only after time calibration, and a genuine
commutant witness that survives phase and time calibration.  Both artifacts
explicitly record `hpc_authorized=false`.

## Extensive family refinement

The quadratic witness is not special to 12-by-12. For every even periodic
square lattice with `L >= 6`, `N=L^2`,

```text
W_L = H_L^2 - (3N/8) I + H_L/2,
tau(W_L^2) = 9N(2N-3)/64,
tau(W_L E5_L) = (N/4) q,
```

where

```text
q = -7807/3600000
    - (66043/57600000) alpha
    - (6119/7200000) alpha^2 != 0.
```

The normalized squared obstruction tends to `2q^2/9`, approximately
`8.35099162871e-6`, rather than disappearing with system size. Exact torus
records at `L=6,8,10,12` reproduce the same per-cell cubic coordinates; `L=4`
is rejected because periodic coordinate aliases invalidate the stable local
geometry.

The complete proof, including the fourth-moment graph count and the
finite-step dual-pairing gate, is in
`docs/report/extensive-commutant-witness-theorem.md`. The machine artifact is

```text
docs/experiments/processor-obstruction/extensive-commutant-witness.json
sha256: 5e5ba831312109df42897c147a85dc4bb77e6670fec4fc6d5befae24bdd8416d
```

Reproduce it with

```bash
PYTHONPATH=src python -u -m scripts.certify_extensive_commutant_witness \
  --lengths 4 6 8 10 12

PYTHONPATH=src python -u -m scripts.certify_extensive_commutant_witness \
  --lengths 4 6 8 10 12 --verify-full
```

The statement above records the historical gate at the time of that artifact.
The separate E9 finite-step production chain now addresses it; no conclusion
from that later chain is back-propagated into this leading-order artifact.

## Exact PF4--TFIM family theorem

The actual two-fragment five-copy Suzuki formula has the exact cyclic
free-trace identity

```text
Tr((A+B) E5) = gamma Tr(C^2 - 4 C D + (8/3) D^2),
C = [A,[A,B]],  D = [B,[B,A]],
gamma = 37/900000
      + (313/14400000) alpha
      + (29/1800000) alpha^2 > 0,
alpha^3 = 4, alpha > 0.
```

The verifier compares all ten cyclic degree-six word classes over the exact
cubic field.  Exact rational Hermitian matrix checks in dimensions two and
three are independent cross-checks, not substitutes for that free-word proof.

For the even periodic TFIM family

```text
A = h sum_i X_i,  B = j sum_i Z_i Z_{i+1},
Tr(C^2)/d = 128 L h^4 j^2,
Tr(CD)/d = 0,
Tr(D^2)/d = 128 L h^2 j^4,
Tr(H E5)/d = 128 L gamma h^2 j^2 (h^2 + 8 j^2/3).
```

Thus both nonzero couplings give a strict leading-order obstruction to
`E5=-i[Q,H]` for every endpoint-conjugation generator `Q`.  The exact
normalized-trace dual bound is also stored in the artifact.  Its scope is
fixed time and fixed Hamiltonian normalization.  Affine `a I + b H` gauge,
finite-step TFIM behavior, and total-time eigenphases modulo `2 pi` remain
explicitly unclaimed.

```text
artifact: tfim-family-obstruction.json
file sha256: a1f628104f816ea2cd35c7846b64c1205df76fae73eecf2b82cc69a1e9f789e4
payload sha256: ffd5bf98bf40138f221dde86b990eb43e97e4e101332bd898465bcbae53c09d3
```

Reproduce and verify from the Issue-128 directory:

```bash
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py \
  --verify docs/experiments/processor-obstruction/tfim-family-obstruction.json
PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py tests/test_tfim_obstruction.py
```

The earlier research-plan coefficients `1/2, 14/3, 4/3` were not a proved BCH
mapping and are superseded by this exact identity.

## Exact dual E7 refinement

A four-shard local contraction has now computed `tau(W E7)` exactly without
constructing the full E7 operator.  The per-cell result is

```text
18151321/8064000000
+ (1494729/1024000000) alpha
+ (96257773/110592000000) alpha^2
~= 0.00676126758644310.
```

It opposes the E5 pairing, but after the extra `1/r^2` suppression it changes
the calibrated signal by only `1.22e-4` at `r=95` and `1.17e-4` at `r=97`.
The exact artifact is

```text
docs/experiments/processor-obstruction/dual-e7-pairing.json
sha256: 9d68771af90096d0ad96bdf6a0861c10d8ec39a864b417d772e51e0c6a9d5d7b
```

The centered-stage path gives a simple principal-log certificate at `r=97`
but not at 95 or 96.  At 97 the only remaining finite-step gap is the
E9-and-higher dual tail.  See `docs/report/dual-e7-remainder-ledger.md` for the
exact margins, branch proof, shard coverage, and independent rerun evidence.

## Finite-step affine effective-spectrum closure

The E9 gap is now closed for the frozen periodic `12 x 12` isotropic
Heisenberg instance at `r=97`.  The manifest-bound 64-shard contraction covers
all 65,536 suffix groups and 262,140 degree-nine words.  The exact per-cell E9
pairing is

```text
-66319553467/248832000000000
- (1446202237/10368000000000) alpha
- (12151491869/82944000000000) alpha^2,
alpha^3 = 4,
```

approximately `-8.57108196229659e-4`.  Combining E5, E7, E9, and the compatible
E11-and-higher tail gives

```text
|q5/97^4 + q7/97^6 + q9/97^8| >= 6.92368710129858e-11,
dual tail                              <= 6.49509254599244e-12,
signed margin                          >= 6.27417784669934e-11 > 0.
```

All inequalities in the certificate use exact rational endpoints.  The
displayed decimals are summaries only.

The nonlinear promotion gate evaluates the centered-moment invariant

```text
Phi_H(A) = m2(H)^3 m3(A)^2 - m3(H)^2 m2(A)^3,
m2(H)=54, m3(H)=-27.
```

It certifies a linear lower bound `0.0576176604311994`, a nonlinear remainder
upper bound `0.0352548202513958`, and hence a strictly positive invariant
margin `0.0223628401798036`.  The resulting machine status is

```text
certified_affine_effective_spectral_obstruction
```

for the principal effective Hamiltonian at one `1/97` step.  Precisely, it
excludes the affine unitary orbit `A = a I + b U H U^dagger`.  It does **not**
claim a lower bound on total-time unitary eigenphases modulo `2 pi`; that would
require a separate no-wrap and no-relabeling theorem.

The tracked artifacts are

```text
dual-e9-pairing.json
file sha256:    796db9b2f2a3635176ea5f100abc5da6d07fac68c6b3f7cd5191f0d927e9e7b1
payload sha256: ea17c2504462d19a102590f5de5046ea24bbc7fb132bebc8e90322aa6b92e144

finite-step-obstruction.json
file sha256:    13805e2d560efce6051f142ef9b69ae53b1f718d37426fb6fbaa10bd3072a909
payload sha256: f972e5b31e8dfe459aa2250ef7d6bf19c57e93d75bbf4f7bd6578b908656e24e

dual-e9-provenance.json
file sha256:    4f9e10686ba2daad8a9b9ee855e12f60851848d16764f33a21681889a9267603
payload sha256: 9676ff11e366509957cbb02b2a6e6548c5ea322c8873cbdfc8943944df223935
```

The compact HPC and independent-rerun hash chain is recorded in
`dual-e9-provenance.json`; the proof and claim boundary are in
`docs/report/finite-step-affine-spectral-ledger.md`.  With the archived
manifest index available as `$E9_RUN_ROOT/manifests/index.json`, reproduce the
full gate from the Issue-128 directory with

```bash
PYTHONPATH=src:. python scripts/certify_dual_e9_pairing.py \
  --index "$E9_RUN_ROOT/manifests/index.json" \
  --verify docs/experiments/processor-obstruction/dual-e9-pairing.json

PYTHONPATH=src:. python scripts/certify_finite_step_obstruction.py \
  --verify docs/experiments/processor-obstruction/finite-step-obstruction.json \
  --e9 docs/experiments/processor-obstruction/dual-e9-pairing.json \
  --index "$E9_RUN_ROOT/manifests/index.json"
```
