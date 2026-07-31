# E7 Checkpointed Tree Reducer Design

## Goal

Turn a complete 272-, 544-, or 1,088-shard exact E7 producer run into a
recoverable exact D6 artifact while using multiple Slurm nodes for the expensive
merge. Preserve the existing fail-closed provenance and exact-arithmetic checks.

## Selected architecture

Use a two-level tree.

1. An array of partial reducers owns disjoint contiguous shard-index ranges.
   Each cell validates every source manifest and payload, merges its range, and
   atomically writes a canonical gzip checkpoint plus a manifest.
2. A final reducer validates that the checkpoints have identical source
   metadata and form an exact, gap-free cover of all source shards and words.
   It merges the much smaller checkpoint set, regenerates exact D6, and writes
   the existing parent, payload, and summary artifact shapes.

The partial checkpoint records its source shard interval, word interval,
source commit, formula, stage count, degree, input digests, and exact canonical
term map. The final parent records both checkpoint and original shard
provenance.

On the current `xhacnormalb` partition, partial cells request 8 CPUs with 16 GB
and final cells request 26 CPUs with 96 GB to satisfy the site's enforced
memory-per-CPU ratio. The Python merge remains single-process; the extra CPU
allocation is a scheduler constraint rather than a parallel-speedup claim.

## Alternatives considered

- Increasing CPU and memory for the monolithic reducer does not help enough:
  the Python merge is effectively single-core and emits no recoverable state.
- Running more independent monolithic reducers improves cross-checking but
  duplicates the same bottleneck.
- A general multi-level scheduler would support arbitrary tree depths, but a
  fixed two-level tree is sufficient for the present 272–1,088 inputs and can
  be deployed and audited within the remaining window.

## Failure behavior

Every output is atomic. A missing, corrupt, overlapping, gapped, or
metadata-inconsistent input aborts the affected stage. Completed partials remain
reusable after timeout. The final D6 artifact is never written unless the full
source range is covered exactly once.

## Verification

Unit tests use small synthetic canonical shards to prove successful reduction,
gap rejection, digest rejection, and equality with a direct exact merge. On
HPC, the partial array must complete before the final reducer is released.
Promotion still requires a second fanout layout to produce an identical exact
D6 coefficient map.
