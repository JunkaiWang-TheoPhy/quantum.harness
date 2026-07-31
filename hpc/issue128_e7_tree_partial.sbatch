#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-tree-part
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:25:00
#SBATCH --output=results/slurm-%x-%A_%a.out
#SBATCH --error=results/slurm-%x-%A_%a.err

set -euo pipefail
: "${ISSUE128_SOURCE_RUN_ID:?}"
: "${ISSUE128_TREE_RUN_ID:?}"
: "${ISSUE128_SOURCE_SHARD_COUNT:?}"
: "${ISSUE128_PARTIAL_COUNT:?}"

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"

partial_index="${SLURM_ARRAY_TASK_ID:?}"
if (( partial_index < 0 || partial_index >= ISSUE128_PARTIAL_COUNT )); then
  echo "array index is outside ISSUE128_PARTIAL_COUNT" >&2
  exit 2
fi
shard_start=$(( ISSUE128_SOURCE_SHARD_COUNT * partial_index / ISSUE128_PARTIAL_COUNT ))
shard_stop=$(( ISSUE128_SOURCE_SHARD_COUNT * (partial_index + 1) / ISSUE128_PARTIAL_COUNT ))
if (( shard_start >= shard_stop )); then
  echo "partial range is empty; reduce ISSUE128_PARTIAL_COUNT" >&2
  exit 2
fi

OUT="$ROOT/results/$ISSUE128_TREE_RUN_ID/partials/$(printf 'partial-%04d' "$partial_index")"
mkdir -p "$OUT"
cd "$ISSUE"
"$ROOT/.venv/bin/python" -u scripts/reduce_e7_word_range.py \
  --input-root "$ROOT/results/$ISSUE128_SOURCE_RUN_ID/words" \
  --shard-count "$ISSUE128_SOURCE_SHARD_COUNT" \
  --shard-start "$shard_start" \
  --shard-stop "$shard_stop" \
  --output "$OUT/partial.json.gz" \
  --manifest "$OUT/manifest.json"
