#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from trottercert.hpc_artifacts import (
    coordinate_decode_terms,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
)


FORMULA_ID = "five_copy_suzuki_fourth_order_exact_cubic"


def _parse_stages(specification: str) -> list[int]:
    stages: set[int] = set()
    for part in specification.split(","):
        if "-" in part:
            first, last = (int(value) for value in part.split("-", 1))
            if first > last:
                raise ValueError("stage range is reversed")
            stages.update(range(first, last + 1))
        else:
            stages.add(int(part))
    if not stages or min(stages) < 0:
        raise ValueError("at least one nonnegative stage is required")
    return sorted(stages)


def _load_stage(root: Path, stage: int, degree: int) -> tuple[dict, dict]:
    cell = root / f"stage-{stage:02d}"
    manifest_path = cell / "manifest.json"
    payload_path = cell / "shard.json.gz"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "complete":
        raise ValueError(f"stage {stage} manifest is not complete")
    if manifest.get("output_sha256") != sha256_file(payload_path):
        raise ValueError(f"stage {stage} payload digest mismatch")
    payload = read_shard_gzip(payload_path)
    if payload.get("formula_id") != FORMULA_ID:
        raise ValueError(f"stage {stage} formula mismatch")
    if payload.get("stage_index") != stage:
        raise ValueError(f"stage {stage} index mismatch")
    order = payload.get("order")
    if isinstance(order, bool) or not isinstance(order, int) or order < degree:
        raise ValueError(f"stage {stage} does not cover degree {degree}")
    series = payload.get("series")
    if not isinstance(series, list) or len(series) != order + 1:
        raise ValueError(f"stage {stage} series coverage mismatch")
    return manifest, {
        "order": order,
        "payload_sha256": sha256_file(payload_path),
        "terms": coordinate_decode_terms(series[degree]),
    }


def compare(
    left_root: Path,
    right_root: Path,
    *,
    stages: list[int],
    degree: int,
    output: Path,
) -> dict[str, object]:
    records: list[dict[str, object]] = []
    for stage in stages:
        left_manifest, left = _load_stage(left_root, stage, degree)
        right_manifest, right = _load_stage(right_root, stage, degree)
        if left["terms"] != right["terms"]:
            raise ArithmeticError(f"stage {stage} degree-{degree} maps differ")
        records.append(
            {
                "stage_index": stage,
                "term_count": len(left["terms"]),
                "left_order": left["order"],
                "right_order": right["order"],
                "left_source_commit": left_manifest.get("git_commit"),
                "right_source_commit": right_manifest.get("git_commit"),
                "left_payload_sha256": left["payload_sha256"],
                "right_payload_sha256": right["payload_sha256"],
                "exact_degree_map_equal": True,
            }
        )
    result = {
        "schema_version": 1,
        "kind": "issue128_cross_order_stage_comparison",
        "status": "complete",
        "degree": degree,
        "stages": records,
    }
    write_manifest_atomic(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left-root", type=Path, required=True)
    parser.add_argument("--right-root", type=Path, required=True)
    parser.add_argument("--stages", required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = compare(
        arguments.left_root,
        arguments.right_root,
        stages=_parse_stages(arguments.stages),
        degree=arguments.degree,
        output=arguments.output,
    )
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
