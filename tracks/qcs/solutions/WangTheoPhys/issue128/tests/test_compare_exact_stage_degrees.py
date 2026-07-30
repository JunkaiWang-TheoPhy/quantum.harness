from __future__ import annotations

from pathlib import Path

import pytest

from scripts.compare_exact_stage_degrees import _parse_stages, compare
from trottercert.cubic_field import Cubic
from trottercert.hpc_artifacts import (
    coordinate_terms_to_json,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)


def _stage(root: Path, *, stage: int, order: int, coefficient: int) -> None:
    cell = root / f"stage-{stage:02d}"
    series = [[] for _ in range(order + 1)]
    series[1] = coordinate_terms_to_json(
        {((0, 0, "X"),): Cubic(coefficient, 0, 0)}
    )
    payload = {
        "schema_version": 1,
        "kind": "issue128_exact_right_generator_stage",
        "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
        "stage_index": stage,
        "stage_count": 31,
        "order": order,
        "series": series,
    }
    payload_path = cell / "shard.json.gz"
    write_shard_gzip(payload_path, payload)
    write_manifest_atomic(
        cell / "manifest.json",
        {
            "status": "complete",
            "git_commit": f"order-{order}",
            "output_sha256": sha256_file(payload_path),
        },
    )


def test_cross_order_exact_degree_comparison(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    _stage(left, stage=25, order=6, coefficient=3)
    _stage(right, stage=25, order=8, coefficient=3)
    output = tmp_path / "comparison.json"
    result = compare(left, right, stages=[25], degree=1, output=output)
    assert result["stages"][0]["term_count"] == 1
    assert output.is_file()


def test_cross_order_mismatch_and_stage_parser(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    _stage(left, stage=25, order=6, coefficient=3)
    _stage(right, stage=25, order=8, coefficient=4)
    with pytest.raises(ArithmeticError, match="maps differ"):
        compare(
            left,
            right,
            stages=[25],
            degree=1,
            output=tmp_path / "comparison.json",
        )
    assert _parse_stages("25-27,30") == [25, 26, 27, 30]
