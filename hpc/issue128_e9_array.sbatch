#!/usr/bin/env bash
#SBATCH --job-name=issue128-e9-array
#SBATCH --array=0-63
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=02:00:00

set -euo pipefail

: "${ISSUE128_ROOT:?set ISSUE128_ROOT to the Issue-128 directory}"
: "${ISSUE128_E9_RUN_ROOT:?set ISSUE128_E9_RUN_ROOT to a new run directory}"
: "${SLURM_ARRAY_TASK_ID:?Slurm array task id is required}"

cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."

index_path="$ISSUE128_E9_RUN_ROOT/manifests/index.json"
manifest_path=$(printf \
  '%s/manifests/manifest-%03d.json' \
  "$ISSUE128_E9_RUN_ROOT" \
  "$SLURM_ARRAY_TASK_ID")
mkdir -p "$ISSUE128_E9_RUN_ROOT/shards"
final_output=$(printf \
  '%s/shards/shard-%03d.json' \
  "$ISSUE128_E9_RUN_ROOT" \
  "$SLURM_ARRAY_TASK_ID")
temporary_output="${final_output}.pending.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

[[ -f "$index_path" ]]
[[ -f "$manifest_path" ]]
[[ ! -e "$final_output" ]]
[[ ! -e "$temporary_output" ]]

python scripts/certify_dual_e9_pairing.py \
  --index "$index_path" \
  --manifest "$manifest_path" \
  --output "$temporary_output"

python scripts/certify_dual_e9_pairing.py \
  --index "$index_path" \
  --verify "$temporary_output"

mv "$temporary_output" "$final_output"
