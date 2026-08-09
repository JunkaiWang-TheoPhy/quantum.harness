#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-dual
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err

set -euo pipefail

: "${ISSUE128_LEFT_RUN_ID:?}"
: "${ISSUE128_LEFT_SHARD_COUNT:?}"
: "${ISSUE128_RIGHT_RUN_ID:?}"
: "${ISSUE128_RIGHT_SHARD_COUNT:?}"
: "${ISSUE128_HARVEST_RUN_ID:?}"

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
OUT="$ROOT/results/$ISSUE128_HARVEST_RUN_ID"
LEFT_OUT="$OUT/left"
RIGHT_OUT="$OUT/right"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$LEFT_OUT" "$RIGHT_OUT" "$OUT/comparison"
cd "$ISSUE"

reduce_lane() {
  local run_id="$1"
  local shard_count="$2"
  local output_root="$3"

  "$ROOT/.venv/bin/python" -u scripts/reduce_e7_word_shards.py \
    --input-root "$ROOT/results/$run_id/words" \
    --shard-count "$shard_count" \
    --output "$output_root/degree-6.json.gz" \
    --parent "$output_root/word-parent.json" \
    --summary "$output_root/summary.json"
}

reduce_lane \
  "$ISSUE128_LEFT_RUN_ID" \
  "$ISSUE128_LEFT_SHARD_COUNT" \
  "$LEFT_OUT" \
  >"$LEFT_OUT/reducer.log" 2>&1 &
left_pid=$!

reduce_lane \
  "$ISSUE128_RIGHT_RUN_ID" \
  "$ISSUE128_RIGHT_SHARD_COUNT" \
  "$RIGHT_OUT" \
  >"$RIGHT_OUT/reducer.log" 2>&1 &
right_pid=$!

while kill -0 "$left_pid" 2>/dev/null || kill -0 "$right_pid" 2>/dev/null; do
  date -u +"heartbeat=%Y-%m-%dT%H:%M:%SZ"
  ps -o pid=,etime=,rss=,pcpu=,command= -p "$left_pid","$right_pid" || true
  sleep 60
done

set +e
wait "$left_pid"
left_rc=$?
wait "$right_pid"
right_rc=$?
set -e

printf 'left_rc=%s right_rc=%s\n' "$left_rc" "$right_rc"
if [[ "$left_rc" -ne 0 || "$right_rc" -ne 0 ]]; then
  exit 1
fi

"$ROOT/.venv/bin/python" -u scripts/compare_exact_degree_payloads.py \
  --left "$LEFT_OUT/degree-6.json.gz" \
  --right "$RIGHT_OUT/degree-6.json.gz" \
  --degree 6 \
  --output "$OUT/comparison/comparison.json"
