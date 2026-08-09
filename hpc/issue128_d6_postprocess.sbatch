#!/usr/bin/env bash
#SBATCH --job-name=issue128-d6-post
#SBATCH --partition=xhacnormalb
#SBATCH --account=giggleliu
#SBATCH --qos=user_acamtw70yu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=96G
#SBATCH --time=01:15:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err

set -euo pipefail

: "${ISSUE128_HARVEST_RUN_ID:?}"
: "${ISSUE128_POSTPROCESS_RUN_ID:?}"
: "${ISSUE128_BUILDER_COMMIT:?}"
export ISSUE128_SELECTED_LANE="${ISSUE128_SELECTED_LANE:-left}"
if [[ "$ISSUE128_SELECTED_LANE" != "left" && "$ISSUE128_SELECTED_LANE" != "right" ]]; then
  echo "selected lane must be left or right" >&2
  exit 2
fi

export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR:?}"
HARVEST="$ROOT/results/$ISSUE128_HARVEST_RUN_ID"
OUT="$ROOT/results/$ISSUE128_POSTPROCESS_RUN_ID"
ISSUE="$OUT/source/issue128"
FINAL="$OUT/final"

test -d "$ISSUE/src"
mkdir -p "$ISSUE/certificates" "$FINAL"
cp "$HARVEST/$ISSUE128_SELECTED_LANE/degree-6.json.gz" \
  "$ISSUE/certificates/issue128-d6-exact.json.gz"
cp "$HARVEST/$ISSUE128_SELECTED_LANE/word-parent.json" \
  "$ISSUE/certificates/issue128-d6-parent.json"
cp "$HARVEST/comparison/comparison.json" \
  "$ISSUE/certificates/issue128-d6-comparison.json"

export PYTHONPATH="$ISSUE/src:$ISSUE${PYTHONPATH:+:$PYTHONPATH}"
cd "$ISSUE"

python3 -u scripts/build_d6_grouped_sidecar.py --candidate-cap 128 \
  | tee "$OUT/grouped-builder.log"
python3 -u scripts/build_d6_integrated_certificate.py \
  | tee "$OUT/certificate-builder.log"
python3 -u scripts/verify.py certificates/issue128-d6-integrated-certificate.json \
  | tee "$OUT/fast-verifier.log"

cp certificates/issue128-d6-exact.json.gz "$FINAL/"
cp certificates/issue128-d6-parent.json "$FINAL/"
cp certificates/issue128-d6-comparison.json "$FINAL/"
cp certificates/issue128-d6-groups.json "$FINAL/"
cp certificates/issue128-d6-integrated-certificate.json "$FINAL/"

python3 - "$FINAL" "$OUT/postprocess-summary.json" \
  "$ISSUE128_BUILDER_COMMIT" "$ISSUE128_SELECTED_LANE" <<'PY'
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

final = Path(sys.argv[1])
output = Path(sys.argv[2])
builder_commit = sys.argv[3]
selected_lane = sys.argv[4]
files = sorted(final.iterdir())
groups = json.loads((final / "issue128-d6-groups.json").read_text())
comparison = json.loads((final / "issue128-d6-comparison.json").read_text())

def fraction(pair):
    return pair[0] / pair[1]

l1_site = fraction(groups["site_pauli_l1_upper"])
grouped_site = fraction(groups["grouped_site_bound"])
payload = {
    "schema_version": 1,
    "kind": "issue128_d6_postprocess_summary",
    "status": "complete",
    "builder_commit": builder_commit,
    "selected_lane": selected_lane,
    "term_count": groups["term_count"],
    "group_count": groups["group_count"],
    "channel_counts": groups["channel_counts"],
    "site_pauli_l1_upper": groups["site_pauli_l1_upper"],
    "grouped_site_bound": groups["grouped_site_bound"],
    "grouped_over_l1_ratio": grouped_site / l1_site,
    "comparison_left_sha256": comparison["left"]["sha256"],
    "comparison_right_sha256": comparison["right"]["sha256"],
    "files": {
        path.name: {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for path in files
    },
}
temporary = output.with_suffix(output.suffix + ".tmp")
temporary.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")
temporary.replace(output)
print(json.dumps(payload, sort_keys=True), flush=True)
PY
