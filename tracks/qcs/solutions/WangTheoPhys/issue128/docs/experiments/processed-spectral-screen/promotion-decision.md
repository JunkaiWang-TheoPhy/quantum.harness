# Processed kernel pre-HPC promotion decision

## Decision

`promote s10+s11 to shared E7 design`

Both primary candidates pass the exact-rational effective-order gate and the
independent dense spectral discovery gate.  This promotion authorizes the
shared physical-word E7 architecture design.  It does **not** authorize a new
HPC submission: the finite-step processed-defect remainder must first show
that an exact D6/D7 result can close either the current-beat or fivefold error
budget.

## Exact audit

| quantity | `s10` | `s11` |
|---|---:|---:|
| merged stages per step | 61 | 67 |
| current-beat step gate | 47 | 43 |
| fivefold step gate | 39 | 35 |
| `R2` free-word l1 | `6.101027913679052e-2` | `4.746167262185952e-32` |
| `R4` free-word l1 | `6.117422601654833e-2` | `4.224879765167434e-2` |
| processed degree-3 residual l1 | `8.777116841420768e-30` | `1.226093209398038e-31` |
| processed degree-5 residual l1 | `1.463642963908243e-29` | `4.625089116627788e-31` |
| processed degree-7 free-word l1 | `7.791094055418490e-1` | `5.604056490876308e-1` |

The degree-three and degree-five entries are nonzero because the published
decimal strings are treated as the exact rational numbers they denote.  They
are retained in the payload rather than rounded to zero or described as exact
algebraic order conditions.

Frozen payload:

```text
artifacts/processed-kernel-audit/s10-s11-local-audit.json
size: 18 MiB
sha256: 992b02d3fae0f16aad9a0791ff58b1f123cee8ef8794e39a825477776470b99b
```

## Dense discovery audit

| kernel | repeated operator slope | repeated phase slope |
|---|---:|---:|
| `s8` | `2.00151769` | `5.99184313` |
| `s10` | `2.00155973` | `5.99712310` |
| `s11` | `4.00258240` | `5.95571683` |

The result confirms that `s10` and `s11` must be presented as spectrally
sixth-order processed kernels, not conventional sixth-order operator
approximations.  The diagnostic is binary64 and non-rigorous; it fixes stage
ordering and qualitative behavior only.

## Regression status

```text
new focused suite: 20 passed in 21.29 s
final full non-slow suite: 179 passed, 12 deselected in 182.80 s
```

## Remaining HPC trigger

Before requesting new compute, the next local theorem prototype must provide:

1. a finite-step norm bound for
   `exp(R(h)) kernel(h) exp(-R(h)) - exp(-i h H)`;
2. an exact contribution from the recorded degree-three and degree-five
   rationalization residuals at `r=39` and `r=35`;
3. a rigorous all-higher-degree envelope small enough that exact physical E7
   has a plausible path to the `1e-6` total budget;
4. a processor near-identity bound suitable for a gap/overlap statement.

If this theorem gate passes, prepare one shared 16,380-word physical evaluator
with two exact weight accumulators and two independent reducer layouts.  If it
fails, stop before HPC and work on the remainder theorem rather than enumerating
E9 through E15 blindly.
