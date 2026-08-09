# Processed finite-step pre-E7 budget

## Result

The exact rationalization residual and proof-only processor gates pass for all
four target integer points.  No point is yet a simulation certificate: the
physical processed-E7 density, all-higher-degree remainder, and a cited local
logarithm theorem are intentionally absent.

| candidate point | degree-3 total | degree-5 total | processor distance `rho` upper | status |
|---|---:|---:|---:|---|
| `s10, r=39` | `3.803026e-31` | `5.231400e-34` | `2.170113e-3` | awaiting proof |
| `s10, r=47` | `2.618562e-31` | `2.480190e-34` | `1.493351e-3` | awaiting proof |
| `s11, r=35` | `6.596204e-33` | `2.548526e-35` | `1.832720e-6` | awaiting proof |
| `s11, r=43` | `4.370119e-33` | `1.118633e-35` | `8.044412e-7` | awaiting proof |

The published-decimal rationalization is therefore negligible relative to
the `10^-6` simulation budget.  The conservative processor rotation also
passes the local `rho < 0.1` screen.  `s11` is especially favorable for a
state-overlap analysis because its `R2` component is nearly zero; this table
does not yet include the separate gap-dependent effective-Hamiltonian term.

## Proof contract

The conditional budget API uses, at total time one,

```text
N * C3_site / r^2
+ N * C5_site / r^4
+ N * C7_site / r^6
+ certified higher-order remainder.
```

It refuses to construct a claim unless `C7_site`, the higher-order remainder,
and a local-log theorem identifier are supplied.  The pre-E7 JSON therefore
fixes each of those fields to `null` and its status to
`awaiting_degree7_and_local_log_remainder`.  The verifier rejects a forged
`accepted` status or a changed processor bound.

## Artifact and reproduction

```text
pre-e7-local-budget.json
size: 11 KiB
sha256: f0ea58e594d02a76bf8e2acd4b8472a0f4c4176f5efa19063043daa7f889a223
```

From the issue-128 directory:

```bash
PYTHONPATH=src python -u -m scripts.report_processed_local_budget \
  --verify \
  --output docs/experiments/processed-finite-step-budget/pre-e7-local-budget.json
python -m pytest -q \
  tests/test_rational_lie_local.py \
  tests/test_processed_local_bounds.py \
  tests/test_processed_budget.py \
  tests/test_processed_local_budget_report.py
```

Observed verification:

```text
focused finite-step suite: 12 passed in 44.04 s
final full non-slow suite: 179 passed, 12 deselected in 182.80 s
```

## Consequence

Do not spend HPC time on coefficient refinement: its contribution is more
than twenty orders of magnitude below the error target.  The next useful work
is the all-order processed local-log remainder.  A shared physical E7 run is
requested only after that theorem produces a numerical envelope showing that
an exact `C7_site` can plausibly close `r=39` or `r=35`.
