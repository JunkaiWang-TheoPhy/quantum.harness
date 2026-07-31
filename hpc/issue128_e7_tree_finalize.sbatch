#!/usr/bin/env bash
#SBATCH --job-name=issue128-e7-tree-final
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err

set -euo pipefail
: "${ISSUE128_TREE_RUN_ID:?}"
: "${ISSUE128_SOURCE_SHARD_COUNT:?}"
: "${ISSUE128_PARTIAL_COUNT:?}"

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
OUT="$ROOT/results/$ISSUE128_TREE_RUN_ID/reduced"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$OUT"
cd "$ISSUE"
"$ROOT/.venv/bin/python" -u scripts/finalize_e7_partials.py \
  --input-root "$ROOT/results/$ISSUE128_TREE_RUN_ID/partials" \
  --partial-count "$ISSUE128_PARTIAL_COUNT" \
  --source-shard-count "$ISSUE128_SOURCE_SHARD_COUNT" \
  --output "$OUT/degree-6.json.gz" \
  --parent "$OUT/tree-parent.json" \
  --summary "$OUT/summary.json"
