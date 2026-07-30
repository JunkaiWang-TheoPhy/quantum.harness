#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-word
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --array=0-15%16
#SBATCH --output=results/slurm-%x-%A_%a.out
#SBATCH --error=results/slurm-%x-%A_%a.err
set -euo pipefail
: "${ISSUE128_RUN_ID:?}"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"
OUT="$ROOT/results/$ISSUE128_RUN_ID/words/$(printf 'word-%02d' "$SLURM_ARRAY_TASK_ID")"
mkdir -p "$OUT"
cd "$ISSUE"
"$ROOT/.venv/bin/python" -u scripts/build_e7_word_shard.py --shard-index "$SLURM_ARRAY_TASK_ID" --shard-count 16 --output "$OUT/shard.json.gz" --manifest "$OUT/manifest.json"
