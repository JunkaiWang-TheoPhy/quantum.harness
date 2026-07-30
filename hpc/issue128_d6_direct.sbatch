#!/usr/bin/env bash
#SBATCH --job-name=issue128-d6-direct
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=64G
#SBATCH --time=01:45:00
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
OUT="$REPOSITORY_ROOT/results/$ISSUE128_RUN_ID/direct"
mkdir -p "$OUT"

cd "$ISSUE_ROOT"
"$REPOSITORY_ROOT/.venv/bin/python" -u scripts/build_d6_direct.py \
  --output "$OUT/degree-6.json.gz" \
  --parent "$OUT/direct-parent.json" \
  --manifest "$OUT/manifest.json"
