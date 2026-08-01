# Exact dual-E7 remainder ledger

## Outcome

The first omitted logarithm coefficient has been contracted exactly against
the extensive commutant witness without constructing the full E7 operator.
On one two-by-two cell,

```text
q5 = tau(W E5)/(N/4)
   = -7807/3600000
     - (66043/57600000) alpha
     - (6119/7200000) alpha^2
   ~= -0.00613020899555631,

q7 = tau(W E7)/(N/4)
   = 18151321/8064000000
     + (1494729/1024000000) alpha
     + (96257773/110592000000) alpha^2
   ~= 0.00676126758644310,

alpha^3 = 4.
```

The signs are opposite, but the total-time logarithm contains
`q5/r^4 + q7/r^6`.  Thus the relative E7 correction is only

```text
|q7/q5| / r^2,
```

which is `1.22e-4` at `r=95` and `1.17e-4` at `r=97`.  E7 cannot cancel or
reverse the leading calibrated spectral obstruction at the relevant step
counts.

## Exact computation

The seventh-degree free logarithm contains 16,380 words in 4,096 common
six-letter suffix groups.  Four deterministic local shards contracted only
Pauli strings of support at most four with `H` and `H^2` on an alias-free
`L=6` torus.  Even-translation invariance converts one target lookup into the
exact nine-cell pairing.

The reducer verified:

```text
group coverage:             4096/4096, no overlap
word coverage:              16380
nonzero words:              12288
retained Pauli occurrences: 30807360
reduction order:            forward equals reverse exactly
```

The reduced artifact is

```text
docs/experiments/processor-obstruction/dual-e7-pairing.json
sha256: 9d68771af90096d0ad96bdf6a0861c10d8ec39a864b417d772e51e0c6a9d5d7b
```

An independent rerun of shard 0 was byte-identical to the stored shard:

```text
sha256: d30e3fd9a8f5020b3e643d13ae31bf1e7d75901f99e2d5207910ba32acebbd5b
wall time: 519.30 s
peak memory: 237224704 bytes
```

The contraction core was separately validated at degree five, where it
reproduced all three frozen E5 cubic pairings exactly.

## Truncated logarithm margins

For `L=12`, `tau(W^2)=23085/4`.  The exact E5+E7 normalized dual magnitude is

| steps `r` | E5 only | E5+E7 | E7 fraction |
|---:|---:|---:|---:|
| 95 | `3.56655082e-11` | `3.56611495e-11` | `1.22210e-4` |
| 96 | `3.42025043e-11` | `3.41984110e-11` | `1.19677e-4` |
| 97 | `3.28137527e-11` | `3.28099062e-11` | `1.17222e-4` |

Every displayed sign is certified using a 50-digit rational interval for the
real root of `alpha^3-4`; the exact cubic coordinates are stored in the
gauge-aware audit.

## Elementary logarithm-branch gate

Each matching fragment contains 72 disjoint Heisenberg bonds and has spectrum
in `[-54,18]`.  Removing the allowed scalar center `-18 I` gives exact radius
36.  The interval-certified sum of absolute merged Suzuki coefficients is at
most

```text
8434723153831272114270621885831 / 10^30.
```

Consequently, the centered one-step unitary path has length at most

```text
36 * sum_j |a_j| / r.
```

Using the rational lower bound `pi > 314159/100000`:

| steps `r` | path upper | certified `< pi` |
|---:|---:|:---:|
| 95 | `3.19631614` | no |
| 96 | `3.16302118` | no |
| 97 | `3.13041272` | yes |

Thus `r=97` has an elementary global principal-log branch certificate after
the allowed scalar phase shift.  Failure at 95 and 96 means only that this
simple path-length proof is insufficient there, not that a branch does not
exist.

## Remaining gate

At `r=97`, the truncated E5+E7 dual obstruction and the elementary logarithm
branch are both certified.  The sole remaining coefficient-level gap is an
E9-and-higher dual tail smaller than `3.28099062e-11` in normalized pairing.

The generic E7 local-l1 majorant is roughly three orders of magnitude too
loose, while the exact dual contraction shows the true E7 correction is only
`1.17e-4` of the leading signal.  The next computation should therefore be a
dual-only E9 plus tail certificate, not a full E7 operator enumeration.

Degree nine has approximately sixteen times as many free words as degree
seven and larger commutator maps.  This is the first point at which a new
shared HPC allocation is justified; no job has been submitted by this work.
Until the E9+ tail is attached, `finite_step_status` remains `inconclusive`.
