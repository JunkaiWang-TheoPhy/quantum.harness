# Dual-E9 HPC deployment and recovery runbook

## Scope and claim boundary

This runbook executes the exact manifest-driven contraction

```text
q9 = tau(W E9)/(L^2/4)
```

for the frozen five-copy fourth-order Suzuki formula.  It does not construct
the full E9 operator and it does not certify the E11-and-higher tail.  A
successful reduction therefore retains

```text
finite_step_status = inconclusive.
```

No cluster hostname, account, partition, username, password, private key, or
other site credential is stored in these scripts.  No HPC job was submitted
while preparing this runbook.

## Measured local calibration

The production-shape manifest preparation and verifier were run locally on
2026-08-01:

| quantity | measured value |
|---|---:|
| formal-log degree | 9 |
| exact words | 262,140 |
| suffix groups | 65,536 |
| child manifests | 64 |
| canonical manifest size | 29 MB |
| preparation wall time | 482.94 s |
| preparation maximum RSS | 914,489,344 bytes |
| preparation peak footprint | 933,577,904 bytes |
| groups per balanced shard | 961--1,086 |
| words per balanced shard | 3,844--4,344 |

The SplitMix64 assignment gives every shard all 64 possible innermost
three-letter suffix patterns.  Direct lexical modulo assignment was rejected
because it pinned those three letters and created severe load imbalance.

The exact target-support reachability prune was benchmarked on E9 group 4818:

| evaluator | wall time | final retained terms |
|---|---:|---:|
| unpruned | 182.03 s | 112,299 |
| reachability-pruned | 5.71 s | 112,299 |

Both paths produced exactly

```text
tau(W E9)_group4818 =
  -8181997/4608000000000
  - (64514921/55296000000000) alpha
  - (39898969/55296000000000) alpha^2,
alpha^3 = 4.
```

A mixed 32-group sample contained 128 words, 84 target-nonzero words, and
580,293 retained terms.  It completed in 21.57 seconds, or 0.674 seconds per
group.  Linear scaling to the measured 1,039-group shard gives about 700
seconds.  This is only a scheduler estimate, not a completed production
result; the two-hour worker request leaves roughly a tenfold wall-time margin.

## Resource envelope

| job | shape | memory | wall time |
|---|---:|---:|---:|
| manifest preparation | 1 CPU | 3 GiB | 00:30:00 |
| contraction array | 64 tasks, 1 CPU each | 3 GiB/task | 02:00:00 |
| reducer | 1 CPU | 3 GiB | 01:00:00 |

The initial site-neutral request used 4 GiB for preparation/workers and 8 GiB
for reduction.  On 2026-08-01, `sbatch --test-only` on `xhacnormalb` rejected
that shape because its memory-per-CPU exceeded `DefMemPerCPU`; no job was
submitted.  The production scripts now request 3 GiB with one CPU.  This still
exceeds the measured 0.94-GB preparation peak by more than threefold, while
avoiding idle CPUs requested only to satisfy the scheduler ratio.

Use a site concurrency cap if the shared filesystem or account policy requires
one.  A starting cap of 16 simultaneous workers is conservative.

## Required environment

Work from one immutable checkout containing the manifest, contraction, and
reducer sources.  Define:

```bash
export ISSUE128_ROOT=/absolute/path/to/tracks/qcs/solutions/WangTheoPhys/issue128
export ISSUE128_E9_RUN_ROOT=/absolute/scratch/path/issue128-e9-YYYYMMDD-HHMMSS
```

`ISSUE128_E9_RUN_ROOT` must be a new path.  The scripts deliberately refuse to
overwrite manifests, worker shards, or the reduced artifact.  Export any
site-specific module or environment setup before `sbatch`, or invoke the jobs
through the site's standard wrapper.  Keep credentials outside the repository.

## Submission

From the Git repository root, set optional site flags as a shell array.  The
empty form is valid when the cluster supplies defaults:

```bash
site_args=()
```

Submit and parse job IDs without relying on human-formatted `sbatch` output:

```bash
manifest_submission=$(sbatch --parsable "${site_args[@]}" \
  hpc/issue128_e9_manifest.sbatch)
manifest_job=${manifest_submission%%;*}

array_submission=$(sbatch --parsable "${site_args[@]}" \
  --array=0-63%16 \
  --dependency="afterok:${manifest_job}" \
  hpc/issue128_e9_array.sbatch)
array_job=${array_submission%%;*}

reduce_submission=$(sbatch --parsable "${site_args[@]}" \
  --dependency="afterok:${array_job}" \
  hpc/issue128_e9_reduce.sbatch)
reduce_job=${reduce_submission%%;*}

printf '%s\n' \
  "manifest_job=${manifest_job}" \
  "array_job=${array_job}" \
  "reduce_job=${reduce_job}"
```

The command-line `--array` value intentionally overrides the uncapped
`#SBATCH --array=0-63` declaration while preserving the exact index set.

## Monitoring

Use scheduler state as operational evidence only:

```bash
squeue -j "${manifest_job},${array_job},${reduce_job}"

sacct -j "${manifest_job},${array_job},${reduce_job}" \
  --format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS,AllocCPUS \
  --units=M
```

Mathematical completion is established by the artifact verifiers, not by a
`COMPLETED` scheduler state.  After the manifest job, verify:

```bash
cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."
python scripts/prepare_dual_log_manifests.py \
  --verify "$ISSUE128_E9_RUN_ROOT/manifests/index.json"
```

After the array, require exactly 64 final shard files:

```bash
find "$ISSUE128_E9_RUN_ROOT/shards" \
  -maxdepth 1 -type f -name 'shard-*.json' | sort
```

The reducer independently verifies each child before summing it.

## Failed-array recovery

Do not remove successful shards and do not resubmit the complete array: worker
scripts reject existing final outputs by design.  Determine failed or missing
indices from `sacct` and the final filenames.  For example, if indices 7, 19,
and 42 are missing, submit only those indices:

```bash
recovery_submission=$(sbatch --parsable "${site_args[@]}" \
  --array=7,19,42 \
  hpc/issue128_e9_array.sbatch)
recovery_job=${recovery_submission%%;*}

reduce_submission=$(sbatch --parsable "${site_args[@]}" \
  --dependency="afterok:${recovery_job}" \
  hpc/issue128_e9_reduce.sbatch)
reduce_job=${reduce_submission%%;*}
```

If manifest preparation fails, use a new run root rather than editing or
reusing a partially published directory.  A failed worker may leave a file
whose name contains `.pending.`; it is not a parent artifact and the reducer's
`shard-*.json` glob ignores it.

## Independent worker rerun

Rerun one selected manifest into a new output directory without changing the
production shard:

```bash
rerun_root="$ISSUE128_E9_RUN_ROOT/independent-rerun"
mkdir -p "$rerun_root"
selected=17
manifest_path=$(printf \
  '%s/manifests/manifest-%03d.json' \
  "$ISSUE128_E9_RUN_ROOT" \
  "$selected")
rerun_output=$(printf '%s/shard-%03d.json' "$rerun_root" "$selected")

cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."
python scripts/certify_dual_e9_pairing.py \
  --index "$ISSUE128_E9_RUN_ROOT/manifests/index.json" \
  --manifest "$manifest_path" \
  --output "$rerun_output"
python scripts/certify_dual_e9_pairing.py \
  --index "$ISSUE128_E9_RUN_ROOT/manifests/index.json" \
  --verify "$rerun_output"
```

Runtime and scheduler metadata are intentionally outside the mathematical
digest.  Compare the mathematical payload digest, not the full-file SHA:

```bash
production_output=$(printf \
  '%s/shards/shard-%03d.json' \
  "$ISSUE128_E9_RUN_ROOT" \
  "$selected")
jq -r '.mathematical_payload_sha256' \
  "$production_output" "$rerun_output"
```

The two printed values must be identical.

## Reduction and local verification

The successful reducer publishes:

```text
$ISSUE128_E9_RUN_ROOT/dual-e9-pairing.json
```

Verify it on the cluster:

```bash
cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."
python scripts/certify_dual_e9_pairing.py \
  --index "$ISSUE128_E9_RUN_ROOT/manifests/index.json" \
  --verify "$ISSUE128_E9_RUN_ROOT/dual-e9-pairing.json"
```

Download the entire `manifests/`, `shards/`, and reduced JSON tree.  In a clean
checkout of the exact source commit, rerun both manifest and reduced-artifact
verification.  Record the source commit, cluster job IDs, `sacct` table,
manifest index SHA-256, 64 worker SHA-256 values, reduced artifact SHA-256, and
the independent mathematical digest comparison in the final evidence ledger.

Only after this verification may `q9` be used to design the separate
E11-and-higher dual-tail certificate.
