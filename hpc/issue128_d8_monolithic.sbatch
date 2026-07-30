#!/usr/bin/env bash
#SBATCH --job-name=issue128-d8-mono
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=63
#SBATCH --mem=240G
#SBATCH --time=04:15:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err

set -euo pipefail
: "${ISSUE128_RUN_ID:?ISSUE128_RUN_ID must name the result run}"

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE_ROOT="$REPOSITORY_ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
export PYTHONPATH="$ISSUE_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
OUT="$REPOSITORY_ROOT/results/$ISSUE128_RUN_ID/monolithic"
mkdir -p "$OUT"
cd "$ISSUE_ROOT"
"$REPOSITORY_ROOT/.venv/bin/python" -u scripts/build_d8_monolithic.py \
  --order 8 \
  --output "$OUT/full-series.json.gz" \
  --manifest "$OUT/manifest.json"
