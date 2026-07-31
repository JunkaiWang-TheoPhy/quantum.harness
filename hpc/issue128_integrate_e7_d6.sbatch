#!/usr/bin/env bash
#SBATCH --job-name=issue128-d6-integrate
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=64G
#SBATCH --time=00:12:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err
set -euo pipefail

: "${ISSUE128_INPUT_RUN_ID:?}"
: "${ISSUE128_INTEGRATE_RUN_ID:?}"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR:?}"
ISSUE="$ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
INPUT="$ROOT/results/$ISSUE128_INPUT_RUN_ID/reduced"
OUT="$ROOT/results/$ISSUE128_INTEGRATE_RUN_ID"
export PYTHONPATH="$ISSUE/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$OUT"

cp "$INPUT/degree-6.json.gz" "$ISSUE/certificates/issue128-d6-exact.json.gz"
cp "$INPUT/word-parent.json" "$ISSUE/certificates/issue128-d6-parent.json"
cp "$INPUT/summary.json" "$OUT/reduction-summary.json"
cd "$ISSUE"
"$ROOT/.venv/bin/python" -u scripts/build_d6_integrated_certificate.py \
  | tee "$OUT/build.json"
"$ROOT/.venv/bin/python" -u scripts/verify.py \
  certificates/issue128-d6-integrated-certificate.json \
  | tee "$OUT/fast-verification.json"
cp certificates/issue128-d6-integrated-certificate.json "$OUT/"
