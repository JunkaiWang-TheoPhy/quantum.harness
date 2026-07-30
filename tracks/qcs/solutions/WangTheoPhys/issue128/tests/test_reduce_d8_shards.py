from __future__ import annotations

import json
from pathlib import Path

from scripts.reduce_d8_shards import reduce_shards
from trottercert.cubic_field import Cubic
from trottercert.exact_series_certificate import verify_exact_degree_payload
from trottercert.hpc_artifacts import (
    coordinate_terms_to_json,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)


def test_reducer_emits_portable_exact_degree_artifact(tmp_path: Path) -> None:
    input_root = tmp_path / "cells"
    terms = {
        0: {((0, 0, "X"),): Cubic.one()},
        1: {((0, 0, "Y"),): Cubic(2, 0, 0)},
    }
    for stage in range(31):
        degree_terms = terms.get(stage, {})
        cell = input_root / f"stage-{stage:02d}"
        payload_path = cell / "shard.json.gz"
        manifest_path = cell / "manifest.json"
        payload = {
            "schema_version": 1,
            "kind": "issue128_exact_right_generator_stage",
            "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
            "stage_index": stage,
            "stage_count": 31,
            "order": 0,
            "series": [coordinate_terms_to_json(degree_terms)],
        }
        write_shard_gzip(payload_path, payload)
        write_manifest_atomic(
            manifest_path,
            {
                "status": "complete",
                "git_commit": "abc123",
                "degree_term_counts": [len(degree_terms)],
                "output_sha256": sha256_file(payload_path),
            },
        )

    output = tmp_path / "degree-0.json.gz"
    summary_path = tmp_path / "summary.json"
    summary = reduce_shards(
        input_root,
        expected_stages=31,
        order=0,
        output=output,
        summary_path=summary_path,
    )
    payload = read_shard_gzip(output)
    verified = verify_exact_degree_payload(
        payload,
        expected_degree=0,
        expected_source_commit="abc123",
    )
    parent = tmp_path / "reduction-parent.json"
    assert verified.term_count == 2
    assert verified.cell_l1_upper == 3
    assert payload["stage_count"] == 31
    assert verified.parent_sha256 == sha256_file(parent)
    assert summary["parent_sha256"] == sha256_file(parent)
    assert json.loads(parent.read_text())["source_stage_count"] == 31
