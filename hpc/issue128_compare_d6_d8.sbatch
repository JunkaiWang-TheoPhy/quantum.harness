#!/usr/bin/env bash
#SBATCH --job-name=issue128-d6-d8-compare
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err

set -euo pipefail
: "${ISSUE128_D6_RUN_ID:?ISSUE128_D6_RUN_ID must name the order-6 run}"
: "${ISSUE128_D8_RUN_ID:?ISSUE128_D8_RUN_ID must name the order-8 run}"

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE_ROOT="$REPOSITORY_ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
export PYTHONPATH="$ISSUE_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
OUTPUT_ROOT="$REPOSITORY_ROOT/results/$ISSUE128_D6_RUN_ID/reduced"
mkdir -p "$OUTPUT_ROOT"

cd "$ISSUE_ROOT"
"$REPOSITORY_ROOT/.venv/bin/python" -u scripts/compare_exact_stage_degrees.py \
  --left-root "$REPOSITORY_ROOT/results/$ISSUE128_D6_RUN_ID/cells" \
  --right-root "$REPOSITORY_ROOT/results/$ISSUE128_D8_RUN_ID/cells" \
  --stages 25-30 \
  --degree 6 \
  --output "$OUTPUT_ROOT/d6-d8-cross-order-comparison.json"
