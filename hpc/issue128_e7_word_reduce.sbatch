#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-reduce
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=96G
#SBATCH --time=00:35:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err
set -euo pipefail
: "${ISSUE128_RUN_ID:?}"
export ISSUE128_SHARD_COUNT="${ISSUE128_SHARD_COUNT:-16}"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"
OUT="$ROOT/results/$ISSUE128_RUN_ID/reduced"
mkdir -p "$OUT"
cd "$ISSUE"
"$ROOT/.venv/bin/python" -u scripts/reduce_e7_word_shards.py --input-root "$ROOT/results/$ISSUE128_RUN_ID/words" --shard-count "$ISSUE128_SHARD_COUNT" --output "$OUT/degree-6.json.gz" --parent "$OUT/word-parent.json" --summary "$OUT/summary.json"
