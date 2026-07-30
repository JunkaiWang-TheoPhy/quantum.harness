#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-compare
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=96G
#SBATCH --time=00:15:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err
set -euo pipefail

: "${ISSUE128_LEFT_RUN_ID:?}"
: "${ISSUE128_RIGHT_RUN_ID:?}"
: "${ISSUE128_COMPARE_RUN_ID:?}"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
OUT="$ROOT/results/$ISSUE128_COMPARE_RUN_ID"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$OUT"
cd "$ISSUE"

"$ROOT/.venv/bin/python" -u scripts/compare_exact_degree_payloads.py \
  --left "$ROOT/results/$ISSUE128_LEFT_RUN_ID/reduced/degree-6.json.gz" \
  --right "$ROOT/results/$ISSUE128_RIGHT_RUN_ID/reduced/degree-6.json.gz" \
  --degree 6 \
  --output "$OUT/comparison.json"
