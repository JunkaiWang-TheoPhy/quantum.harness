#!/usr/bin/env bash
#SBATCH --job-name=issue128-e9-manifest
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00

set -euo pipefail

: "${ISSUE128_ROOT:?set ISSUE128_ROOT to the Issue-128 directory}"
: "${ISSUE128_E9_RUN_ROOT:?set ISSUE128_E9_RUN_ROOT to a new run directory}"

cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."

mkdir -p "$ISSUE128_E9_RUN_ROOT"
final_output="$ISSUE128_E9_RUN_ROOT/manifests"
temporary_output="$ISSUE128_E9_RUN_ROOT/manifests.pending.${SLURM_JOB_ID}"

[[ ! -e "$final_output" ]]
[[ ! -e "$temporary_output" ]]

python scripts/prepare_dual_log_manifests.py \
  --degree 9 \
  --shard-count 64 \
  --output-dir "$temporary_output"

python scripts/prepare_dual_log_manifests.py \
  --verify "$temporary_output/index.json"

mv "$temporary_output" "$final_output"
