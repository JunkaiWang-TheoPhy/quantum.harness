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
eigenphase lower bound at total time one: that promotion still requires a
certified local-log branch and an all-order remainder small enough to preserve
the leading pairing.

Likewise, the support-six result remains a rigorous no-go for the stated
short-support processor class in operator space.  Support alone does not
lower-bound eigenphase error, so that result is not promoted to a spectral
no-go.

The fail-closed refinement is recorded in:

```text
quadratic-commutant-witness.json
sha256: d0892b7320cb22a8be1701a1fdb185d2ce9dbfbc20d191fc21bf4359cea5d391
gauge-aware-audit.json
sha256: 6c25e184ad5815b67a03e1eaad7c61b124d50803aa5c5f22da6c493821fe8368
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

The finite-step status remains `inconclusive`: a compatible logarithm branch
and a Hilbert--Schmidt remainder satisfying
`tau(R_r^2) < rho_L/r^8` are still required. This refinement does not
authorize E7 or any shared HPC run.
