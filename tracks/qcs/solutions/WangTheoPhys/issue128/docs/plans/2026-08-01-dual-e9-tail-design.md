# Dual E9 and all-order tail certificate design

**Date:** 2026-08-01

**Status:** Approved for implementation by the user's standing `做吧`
instruction after the exact dual-E7 review.

## Objective

Promote the calibrated leading-order obstruction at `r=97` toward a finite-step
eigenphase theorem by certifying the only unresolved coefficient-level gate:

```text
|sum_(odd d >= 9) q_d / 97^(d-1)|
  < |q_5 / 97^4 + q_7 / 97^6|,

q_d = tau(W_L E_d)/(L^2/4).
```

The computation must remain dual-only.  It must not construct the complete E9
operator, infer an all-order theorem from a finite list of coefficients, or
replace the required local-log Hilbert--Schmidt remainder with the existing
right-generator operator-norm error ledger.

## Current evidence and measured scale

The exact E7 contraction contains 16,380 words in 4,096 common-suffix groups.
One independent four-way shard took 519.30 seconds and 237,224,704 bytes.

A read-only local scale probe generated the complete degree-nine free-log word
map and measured:

```text
words:                 262140
nonzero coefficients:  262140
suffix groups:          65536
group size:             3 or 4 words
wall time:              344.54 s
maximum resident set:   907427840 bytes
peak memory footprint:  933069904 bytes
```

Thus E9 is exactly sixteen times the E7 word/group count.  Rebuilding this map
inside every array task would waste roughly 6.1 CPU-hours across 64 shards and
multiply shared-node startup memory pressure.

A production-shape calibration exposed a second constraint.  Lexical suffix
ordinals are base-four integers, so direct `ordinal % 64` assignment fixes the
three innermost fragment labels in each worker.  Shard zero therefore received
only suffixes ending in `000`, whose nested commutators vanish, while other
workers would carry the nonzero load.  The selected assignment now applies the
fixed SplitMix64 avalanche to the ordinal before reduction modulo the shard
count.  For all 65,536 groups this gives 961--1,086 groups per shard, and every
shard contains all 64 possible innermost three-letter patterns.

The balanced 64-file preparation completed in 482.94 seconds with a
933,577,904-byte peak footprint and produced 29 MB of canonical JSON.  Child
sizes range from 961 to 1,086 groups and 3,844 to 4,344 words.

## Approaches considered

### A. Rebuild the word map in every worker

Reuse `contract_log_degree_shard` unchanged with `degree=9`.  This minimizes
new code, but every worker repeats the measured 344-second, 0.93-GB formal-log
construction.  The repeated work also makes the worker artifact insufficient
to distinguish word-generation failures from Pauli-contraction failures.

### B. Pre-sharded canonical word manifests (selected)

Run one preparation job that generates the degree-nine map once, orders words
by `(suffix, first_letter)`, assigns suffix groups through a deterministic
SplitMix64 avalanche of the lexical ordinal, and
writes one canonical manifest per worker.  Every manifest records exact cubic
coefficients, group ordinals, the full manifest-set digest, and implementation
source hashes.  Workers consume only one manifest and perform the existing
local Pauli contraction.  A reducer verifies manifest and contraction coverage
before summing exact cubic pairings.

This separates the expensive deterministic formal-log stage from the
embarrassingly parallel contraction, permits a byte-identical manifest rerun,
and keeps each production task small enough for ordinary CPU partitions.

### C. Streaming ordinal generation

Replace the current formal-log dictionary with a new noncommutative streaming
recurrence that emits one suffix range at a time.  This could eventually avoid
the shared manifest, but it creates a second exact-series implementation and a
substantially larger proof surface.  It is deferred unless manifest I/O is
shown to dominate the actual E9 run.

## Architecture

### Generic exact word-manifest layer

Add a module whose trusted outputs are immutable records

```text
degree, shard_index, shard_count,
total_groups, group_ordinals,
[(word, exact_cubic_coefficient), ...],
formula configuration, source manifest, parent-set digest.
```

The preparation command must support any odd degree at least three.  It will
write all shard files to a temporary directory, validate their disjoint union,
compute every file digest, then publish a small index last.  Workers reject a
manifest unless its own digest, index digest, degree, shard count, formula, and
source manifest agree.

Words are encoded as compact base-four digit strings.  Cubic coefficients use
three canonical rational pairs.  JSON is retained for independent review;
gzip may be used only as a transport wrapper, with the digest binding the
canonical uncompressed bytes.

### Manifest-driven dual contraction

Refactor the contraction core so it can consume either a generated word map
or an iterable of manifest groups.  The scientific evaluator is unchanged:

1. evaluate one shared suffix tail once through the existing cache;
2. retain only final Pauli strings of support at most four;
3. lift one alias-free `L=6` density to nine two-by-two cells;
4. contract exactly with `H`, `H^2`, and
   `W=H^2-(3N/8)I+H/2`; and
5. clear the evaluator cache at group boundaries.

The manifest path adds an exact reachability prune that is absent from the
frozen E7 implementation.  At an intermediate suffix of length `m`, there are
`degree-m` outer commutators left and support can decrease by at most one per
commutator.  Terms with support greater than

```text
4 + degree - m
```

therefore cannot return to the support-at-most-four dual target and are deleted
immediately.  Full E5 and selected E7 comparisons require exact agreement with
the unpruned evaluator.  On the measured E9 heavy group, this reduced wall time
from 182.03 seconds to 5.71 seconds without changing any pairing; a 32-group
mixed sample took 21.57 seconds, or 0.674 seconds per group.

The worker artifact records word, nonzero-word, retained-term, and group
counts, exact `tau_h`, `tau_h2`, and `tau_w`, peak RSS, wall time, scheduler
metadata when present, manifest SHA-256, and implementation hashes.

### E9 reducer and decision gate

The reducer requires all shard indices exactly once, exact group coverage
`0..65535`, matching manifests and sources, canonical parent digests, and
forward/reverse equality of all cubic sums.  It emits `q9=tau(W E9)/(N/4)` and
the exact contribution `q9/97^8` alongside the frozen E5+E7 pairing.

E9 alone never changes `finite_step_status`.  It only selects the next tail
action:

- if an independently certified E11+ majorant fits inside the remaining
  signed margin, emit `tail_gate=pass`;
- otherwise emit the exact deficit and request a higher-degree refinement;
- never infer a geometric ratio solely from the observed q5, q7, and q9.

### E11-and-higher dual majorant

Build a separate nonnegative recurrence by moving nested commutators across
the normalized trace.  It starts from the Pauli terms of `H` and `H^2`, tracks
support size, connected-component count, overlap choices, and whether a path
can return to support at most four, and applies outward absolute cubic stage
coefficients.  Formal-log coefficient growth and Pauli-growth transitions are
kept as separate factors.

The recurrence must prove a degree-uniform two-step domination after an
explicit anchor degree.  In matrix form, with nonnegative rational transition
matrix `M`, it must provide a positive rational Collatz vector `v` and
`lambda` such that

```text
M v <= lambda v,       lambda / 97^2 < 1.
```

The tail is then summed as an outward rational geometric series.  Every state
aggregation must be justified by a worst-case transition inequality; empirical
counts or sampled supports may guide the partition but cannot enter the
trusted certificate.

If no sufficiently small domination is found, the artifact remains a failed
majorant and E11 exact contraction becomes the next explicit gate.  A loose
bound is evidence about the method, not evidence that the physical tail is
large.

## Data flow

```text
frozen Suzuki stages + exact cubic field
                |
                v
       one E9 manifest-prep job
                |
       canonical index + 64 manifests
                |
                v
       64 independent dual workers
                |
       exact shard artifacts + digests
                |
                v
       coverage-checking E9 reducer
                |
                +----> exact q9 and signed r=97 margin
                |
reverse-support transition proof
                |
                v
       rational E11+ tail certificate
                |
                v
       fail-closed finite-step audit
```

## Error handling and claim boundaries

- Reject even degrees, inconsistent formulas, noncanonical rationals, missing
  source files, malformed words, duplicate or missing groups, and changed
  parent digests.
- Publish the manifest index only after all shard files have been written and
  verified.
- Treat interrupted or partial workers as absent, never as zero contributions.
- Keep scheduler metadata outside the mathematical digest when necessary, but
  bind the mathematical payload and record the metadata policy explicitly.
- Require a byte-identical rerun of one manifest and one contraction shard.
- Keep `finite_step_status=inconclusive` unless the branch, exact E5+E7(+E9),
  and all-order E11+ tail are simultaneously verified.
- Do not use the old D8+ right-generator operator-norm tail as the local-log
  Hilbert--Schmidt remainder.  They bound different objects and normalizations.
- Do not claim an all-even-L finite-step theorem from a finite-torus numerical
  check; retain the finite-range/alias-free argument used by the extensive
  witness theorem.

## Testing and verification

1. Degree-five manifests must reproduce the frozen E5 `tau_h`, `tau_h2`, and
   `tau_w` exactly.
2. Degree-seven manifest-driven shards must reduce to the committed E7 artifact
   byte-for-byte at the mathematical payload level.
3. Mutation tests must reject changed words, coefficients, group ordinals,
   source hashes, index hashes, parent hashes, resource counters, and claims.
4. A small degree-nine shard must run locally and agree with the legacy
   in-memory path for the same group ordinals.
5. The tail recurrence must be checked both state-by-state and by an independent
   verifier that recomputes `M v <= lambda v` and the rational tail sum.
6. The focused suite and complete non-slow Issue-128 suite must pass before any
   promotion commit.

## HPC envelope

The initial production request is:

- one manifest-prep job: 1 CPU, 3 GiB, 30 minutes;
- 64 contraction array tasks: 1 CPU, 3 GiB, 2 hours each;
- one reducer: 1 CPU, 3 GiB, 1 hour;
- array concurrency capped according to shared-filesystem policy.

The final 3-GiB shape incorporates the target site's `DefMemPerCPU` preflight
gate.  It remains more than three times the measured 0.94-GB preparation peak;
the rejected 4/8-GiB site-neutral shape never entered the queue.

These are conservative envelopes derived from the measured E7 shard and E9
word-map probe.  A local one-group and one-percent-shard calibration must be
used to tighten time and memory before submission.  No HPC job is submitted by
this design document.

## Success criterion

The phase succeeds only when the repository contains:

1. an exact, hash-bound E9 dual pairing with complete shard coverage;
2. an independently verifiable E11+ tail certificate strictly inside the
   remaining `r=97` signed margin;
3. the already certified principal-log branch; and
4. a regenerated fail-closed audit that promotes the finite-step status only
   after all three inputs verify.

If item 2 fails, the exact E9 artifact remains useful evidence, but the
finite-step theorem remains open.
