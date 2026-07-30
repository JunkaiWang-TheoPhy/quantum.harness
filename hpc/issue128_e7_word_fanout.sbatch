#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-fanout
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --array=0-15%16
#SBATCH --output=results/slurm-%x-%A_%a.out
#SBATCH --error=results/slurm-%x-%A_%a.err
set -euo pipefail

: "${ISSUE128_RUN_ID:?}"
: "${ISSUE128_SOURCE_COMMIT:?}"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ISSUE128_WORKERS_PER_TASK="${ISSUE128_WORKERS_PER_TASK:-17}"
export ISSUE128_SHARD_COUNT="${ISSUE128_SHARD_COUNT:-272}"

if (( ISSUE128_WORKERS_PER_TASK * 16 != ISSUE128_SHARD_COUNT )); then
  printf 'expected shard count = 16 array cells * workers per cell\n' >&2
  exit 2
fi

ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
export ROOT ISSUE PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"

srun --kill-on-bad-exit=1 --exclusive \
  --ntasks="$ISSUE128_WORKERS_PER_TASK" --cpus-per-task=1 \
  bash -c '
    set -euo pipefail
    shard_index=$((SLURM_ARRAY_TASK_ID * ISSUE128_WORKERS_PER_TASK + SLURM_PROCID))
    out="$ROOT/results/$ISSUE128_RUN_ID/words/$(printf "word-%02d" "$shard_index")"
    mkdir -p "$out"
    cd "$ISSUE"
    "$ROOT/.venv/bin/python" -u scripts/build_e7_word_shard.py \
      --shard-index "$shard_index" \
      --shard-count "$ISSUE128_SHARD_COUNT" \
      --output "$out/shard.json.gz" \
      --manifest "$out/manifest.json"
  '
