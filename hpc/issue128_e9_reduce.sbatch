#!/usr/bin/env bash
#SBATCH --job-name=issue128-e9-reduce
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00

set -euo pipefail
shopt -s nullglob

: "${ISSUE128_ROOT:?set ISSUE128_ROOT to the Issue-128 directory}"
: "${ISSUE128_E9_RUN_ROOT:?set ISSUE128_E9_RUN_ROOT to a new run directory}"

cd "$ISSUE128_ROOT"
export PYTHONPATH="src:."

index_path="$ISSUE128_E9_RUN_ROOT/manifests/index.json"
shards=("$ISSUE128_E9_RUN_ROOT"/shards/shard-*.json)
[[ "${#shards[@]}" -eq 64 ]]

for shard in "${shards[@]}"; do
  python scripts/certify_dual_e9_pairing.py \
    --index "$index_path" \
    --verify "$shard"
done

final_output="$ISSUE128_E9_RUN_ROOT/dual-e9-pairing.json"
temporary_output="${final_output}.pending.${SLURM_JOB_ID}"
[[ ! -e "$final_output" ]]
[[ ! -e "$temporary_output" ]]

python scripts/certify_dual_e9_pairing.py \
  --index "$index_path" \
  --reduce "${shards[@]}" \
  --output "$temporary_output"

python scripts/certify_dual_e9_pairing.py \
  --index "$index_path" \
  --verify "$temporary_output"

mv "$temporary_output" "$final_output"
