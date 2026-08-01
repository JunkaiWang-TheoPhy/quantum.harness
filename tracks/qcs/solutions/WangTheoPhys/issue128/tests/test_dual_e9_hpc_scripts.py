from __future__ import annotations

from pathlib import Path

import pytest

ISSUE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ISSUE_ROOT.parents[4]
HPC_ROOT = REPO_ROOT / "hpc"


@pytest.mark.parametrize(
    "name",
    [
        "issue128_e9_manifest.sbatch",
        "issue128_e9_array.sbatch",
        "issue128_e9_reduce.sbatch",
    ],
)
def test_e9_slurm_script_is_fail_closed(name: str) -> None:
    text = (HPC_ROOT / name).read_text()
    assert "set -euo pipefail" in text
    assert 'export PYTHONPATH="src:."' in text
    assert "ISSUE128_ROOT" in text
    assert "ISSUE128_E9_RUN_ROOT" in text
    assert "--verify" in text
    assert "mv \"$temporary_output\" \"$final_output\"" in text
    forbidden = ("--account=", "--partition=", "PRIVATE KEY", "password=")
    assert all(value not in text for value in forbidden)


def test_array_has_exact_production_shape() -> None:
    text = (HPC_ROOT / "issue128_e9_array.sbatch").read_text()
    assert "#SBATCH --array=0-63" in text
    assert "#SBATCH --cpus-per-task=1" in text
    assert "#SBATCH --mem=4G" in text
    assert "#SBATCH --time=02:00:00" in text
    assert "SLURM_ARRAY_TASK_ID" in text
    assert "manifest-%03d.json" in text
    assert "shard-%03d.json" in text


def test_manifest_and_reducer_resource_envelopes() -> None:
    manifest = (HPC_ROOT / "issue128_e9_manifest.sbatch").read_text()
    reducer = (HPC_ROOT / "issue128_e9_reduce.sbatch").read_text()
    assert "#SBATCH --mem=4G" in manifest
    assert "#SBATCH --time=00:30:00" in manifest
    assert "#SBATCH --mem=8G" in reducer
    assert "#SBATCH --time=01:00:00" in reducer
    assert '[[ "${#shards[@]}" -eq 64 ]]' in reducer
