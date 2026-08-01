#!/usr/bin/env bash
#SBATCH --job-name=issue128-xxz-verify
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err

set -euo pipefail

: "${ISSUE128_RUN_ID:?ISSUE128_RUN_ID must identify this immutable run}"
: "${ISSUE128_GROUPED_SOURCE_COMMIT:?missing frozen generator commit}"
: "${ISSUE128_SUBMISSION_COMMIT:?missing submission artifact commit}"
: "${ISSUE128_BUILD_SOURCE_ROOT:?missing staged frozen generator tree}"
: "${ISSUE128_SUBMISSION_SOURCE_ROOT:?missing staged submission tree}"

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

REMOTE_REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?}"
SUBMISSION_ROOT="$ISSUE128_SUBMISSION_SOURCE_ROOT"
SUBMISSION_ISSUE="$SUBMISSION_ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
OUT="$REMOTE_REPOSITORY_ROOT/results/$ISSUE128_RUN_ID/grouped-xxz"
BUILD_ROOT="$ISSUE128_BUILD_SOURCE_ROOT"
BUILD_ISSUE="$BUILD_ROOT/tracks/qcs/solutions/WangTheoPhys/issue128"
FRESH="$OUT/fresh"
PYTHON_BIN="${ISSUE128_PYTHON:-python3}"

test "$(tr -d '\n' < "$SUBMISSION_ROOT/.issue128-commit")" = \
  "$ISSUE128_SUBMISSION_COMMIT"
test "$(tr -d '\n' < "$BUILD_ROOT/.issue128-commit")" = \
  "$ISSUE128_GROUPED_SOURCE_COMMIT"
test -f "$SUBMISSION_ISSUE/benchmarks/paper-a/grouped-xxz-pilot.json"
test -f "$SUBMISSION_ISSUE/benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz"
test ! -e "$OUT"
mkdir -p "$OUT" "$FRESH"

set +e
/usr/bin/time -v "$PYTHON_BIN" -u \
  "$SUBMISSION_ISSUE/scripts/reference_verify_grouped_xxz.py" \
  "$SUBMISSION_ISSUE/benchmarks/paper-a/grouped-xxz-pilot.json" \
  "$SUBMISSION_ISSUE/benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz" \
  >"$OUT/reference.stdout" 2>"$OUT/reference.stderr"
REFERENCE_RC=$?

if [[ "$REFERENCE_RC" -eq 0 ]]; then
  /usr/bin/time -v "$PYTHON_BIN" -u \
    "$BUILD_ISSUE/scripts/compile_grouped_xxz.py" --build \
    --delta 1/2 --pipeline direct-theorem \
    --summary "$FRESH/delta-1-2.json" \
    --witness "$FRESH/delta-1-2-witness.json.gz" \
    >"$OUT/build-delta-1-2.stdout" 2>"$OUT/build-delta-1-2.stderr"
  DELTA_HALF_RC=$?
else
  DELTA_HALF_RC=125
fi

if [[ "$DELTA_HALF_RC" -eq 0 ]]; then
  /usr/bin/time -v "$PYTHON_BIN" -u \
    "$BUILD_ISSUE/scripts/compile_grouped_xxz.py" --build \
    --delta 2/1 --pipeline direct-theorem \
    --summary "$FRESH/delta-2-1.json" \
    --witness "$FRESH/delta-2-1-witness.json.gz" \
    >"$OUT/build-delta-2-1.stdout" 2>"$OUT/build-delta-2-1.stderr"
  DELTA_TWO_RC=$?
else
  DELTA_TWO_RC=125
fi

if [[ "$DELTA_TWO_RC" -eq 0 ]]; then
  /usr/bin/time -v "$PYTHON_BIN" -u \
    "$BUILD_ISSUE/scripts/compile_grouped_xxz.py" --merge-pilot \
    "$FRESH/delta-1-2.json" "$FRESH/delta-1-2-witness.json.gz" \
    "$FRESH/delta-2-1.json" "$FRESH/delta-2-1-witness.json.gz" \
    --summary "$FRESH/grouped-xxz-pilot.json" \
    --witness "$FRESH/grouped-xxz-pilot-witness.json.gz" \
    >"$OUT/merge.stdout" 2>"$OUT/merge.stderr"
  MERGE_RC=$?
else
  MERGE_RC=125
fi

if [[ "$MERGE_RC" -eq 0 ]]; then
  cmp -s "$FRESH/grouped-xxz-pilot.json" \
    "$SUBMISSION_ISSUE/benchmarks/paper-a/grouped-xxz-pilot.json"
  SUMMARY_CMP_RC=$?
  cmp -s "$FRESH/grouped-xxz-pilot-witness.json.gz" \
    "$SUBMISSION_ISSUE/benchmarks/paper-a/grouped-xxz-pilot-witness.json.gz"
  WITNESS_CMP_RC=$?
else
  SUMMARY_CMP_RC=125
  WITNESS_CMP_RC=125
fi
set -e

"$PYTHON_BIN" - "$OUT" "$ISSUE128_RUN_ID" \
  "$ISSUE128_GROUPED_SOURCE_COMMIT" "$ISSUE128_SUBMISSION_COMMIT" \
  "$REFERENCE_RC" "$DELTA_HALF_RC" "$DELTA_TWO_RC" "$MERGE_RC" \
  "$SUMMARY_CMP_RC" "$WITNESS_CMP_RC" <<'PY'
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

out = Path(sys.argv[1])
run_id = sys.argv[2]
source_commit = sys.argv[3]
submission_commit = sys.argv[4]
names = (
    "reference",
    "delta_half_build",
    "delta_two_build",
    "merge",
    "summary_byte_compare",
    "witness_byte_compare",
)
exit_codes = {name: int(value) for name, value in zip(names, sys.argv[5:], strict=True)}
files = {}
for path in sorted(out.rglob("*")):
    if path.is_file() and path.name != "manifest.json":
        files[path.relative_to(out).as_posix()] = {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
payload = {
    "schema": "issue128_grouped_xxz_hpc_verification_v1",
    "run_id": run_id,
    "source_commit": source_commit,
    "submission_commit": submission_commit,
    "status": "complete" if all(code == 0 for code in exit_codes.values()) else "failed",
    "exit_codes": exit_codes,
    "files": files,
}
raw = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()
temporary = out / "manifest.json.tmp"
temporary.write_bytes(raw)
temporary.replace(out / "manifest.json")
print(raw.decode(), end="", flush=True)
PY

test "$REFERENCE_RC" -eq 0
test "$DELTA_HALF_RC" -eq 0
test "$DELTA_TWO_RC" -eq 0
test "$MERGE_RC" -eq 0
test "$SUMMARY_CMP_RC" -eq 0
test "$WITNESS_CMP_RC" -eq 0
