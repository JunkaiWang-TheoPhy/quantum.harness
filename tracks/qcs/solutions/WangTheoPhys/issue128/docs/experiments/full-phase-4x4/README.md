# Full-phase 4x4 D4 patch audit

This directory is an experimental, non-production feasibility audit. It does
not change the frozen Issue #128 certificate.

Run from the Issue #128 solution directory:

```bash
PYTHONPATH=src:. python scripts/experimental_full_phase_4x4.py
```

The experiment retains all four colored-cell phases in one open 4x4 patch.
It rejects periodic 4x4 aliasing, because wrapping a long-support Pauli term
on a small torus is not a valid cluster decomposition of the 12x12 operator.
The JSON output distinguishes exact rational upper bounds from the resulting
go/no-go interpretation.

The current experiment returns **NO-GO** for this patch size. A fixed-origin,
phase-preserving open 4x4 patch contains only 2,562 of 75,324 D4 terms
(7.80% of interval Pauli-l1 weight). Combining its strict anticommuting bound
with the unchanged outside partition gives a cell bound of 6.7389332, worse
than the frozen 6.4729265. More decisively, at r=78 the frozen D5+D7+tail
contributions already sum to 1.0813241e-6 even if both D4 and D6 vanished.
