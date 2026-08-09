from __future__ import annotations

import json
from pathlib import Path

import pytest
from scripts.finalize_e7_partials import finalize_partials
from scripts.reduce_e7_word_range import reduce_word_range
from trottercert.cubic_field import Cubic
from trottercert.hpc_artifacts import (
    coordinate_decode_terms,
    coordinate_terms_to_json,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)


def _write_source(
    root: Path,
    *,
    index: int,
    count: int,
    first: int,
    last: int,
    key: tuple[tuple[int, int, str], ...],
    coefficient: Cubic,
) -> None:
    cell = root / f"word-{index:02d}"
    payload_path = cell / "shard.json.gz"
    manifest_path = cell / "manifest.json"
    payload = {
        "schema_version": 1,
        "kind": "issue128_exact_e7_word_shard",
        "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
        "source_commit": "fixture-commit",
        "stage_count": 31,
        "degree": 7,
        "shard_index": index,
        "shard_count": count,
        "word_first": first,
        "word_last": last,
        "word_total": count,
        "terms": coordinate_terms_to_json({key: coefficient}),
    }
    write_shard_gzip(payload_path, payload)
    write_manifest_atomic(
        manifest_path,
        {
            "schema_version": 1,
            "kind": "issue128_exact_e7_word_shard_manifest",
            "status": "complete",
            "source_commit": "fixture-commit",
            "shard_index": index,
            "shard_count": count,
            "word_first": first,
            "word_last": last,
            "output_sha256": sha256_file(payload_path),
        },
    )


def _sources(root: Path) -> None:
    _write_source(
        root,
        index=0,
        count=2,
        first=0,
        last=1,
        key=((0, 0, "X"),),
        coefficient=Cubic.one(),
    )
    _write_source(
        root,
        index=1,
        count=2,
        first=1,
        last=2,
        key=((0, 0, "X"),),
        coefficient=Cubic(2, 0, 0),
    )


def test_reduce_word_range_writes_exact_atomic_checkpoint(tmp_path: Path) -> None:
    source = tmp_path / "words"
    _sources(source)
    output = tmp_path / "partial.json.gz"
    manifest = tmp_path / "manifest.json"

    summary = reduce_word_range(
        input_root=source,
        shard_count=2,
        shard_start=0,
        shard_stop=2,
        output=output,
        manifest=manifest,
    )

    payload = read_shard_gzip(output)
    terms = coordinate_decode_terms(payload["terms"])
    assert terms[((0, 0, "X"),)] == Cubic(3, 0, 0)
    assert payload["shard_start"] == 0
    assert payload["shard_stop"] == 2
    assert len(payload["sources"]) == 2
    assert summary["status"] == "complete"
    assert json.loads(manifest.read_text())["output_sha256"] == sha256_file(output)


def test_reduce_word_range_rejects_corrupt_source_digest(tmp_path: Path) -> None:
    source = tmp_path / "words"
    _sources(source)
    payload = read_shard_gzip(source / "word-01" / "shard.json.gz")
    payload["terms"] = coordinate_terms_to_json(
        {((0, 0, "X"),): Cubic(4, 0, 0)}
    )
    write_shard_gzip(source / "word-01" / "shard.json.gz", payload)

    with pytest.raises(ValueError, match="incomplete or corrupt"):
        reduce_word_range(
            input_root=source,
            shard_count=2,
            shard_start=0,
            shard_stop=2,
            output=tmp_path / "partial.json.gz",
            manifest=tmp_path / "manifest.json",
        )


def test_finalize_partials_emits_exact_d6_shape_from_fixture(
    tmp_path: Path,
) -> None:
    source = tmp_path / "words"
    partials = tmp_path / "partials"
    _sources(source)
    for index in range(2):
        cell = partials / f"partial-{index:04d}"
        reduce_word_range(
            input_root=source,
            shard_count=2,
            shard_start=index,
            shard_stop=index + 1,
            output=cell / "partial.json.gz",
            manifest=cell / "manifest.json",
        )

    output = tmp_path / "degree-6.json.gz"
    parent = tmp_path / "parent.json"
    summary_path = tmp_path / "summary.json"
    summary = finalize_partials(
        input_root=partials,
        partial_count=2,
        source_shard_count=2,
        output=output,
        parent=parent,
        summary=summary_path,
        e5_override={((0, 0, "Z"),): Cubic.one()},
    )

    payload = read_shard_gzip(output)
    assert payload["kind"] == "issue128_exact_right_generator_degree"
    assert payload["degree"] == 6
    assert payload["term_count"] > 0
    assert payload["parent_sha256"] == sha256_file(parent)
    assert summary["status"] == "complete"
    assert summary["output_sha256"] == sha256_file(output)


def test_finalize_partials_rejects_gap_before_output(tmp_path: Path) -> None:
    source = tmp_path / "words"
    partials = tmp_path / "partials"
    _sources(source)
    for index, shard_index in enumerate((0, 0)):
        cell = partials / f"partial-{index:04d}"
        reduce_word_range(
            input_root=source,
            shard_count=2,
            shard_start=shard_index,
            shard_stop=shard_index + 1,
            output=cell / "partial.json.gz",
            manifest=cell / "manifest.json",
        )
    output = tmp_path / "degree-6.json.gz"

    with pytest.raises(ValueError, match="gap or overlap"):
        finalize_partials(
            input_root=partials,
            partial_count=2,
            source_shard_count=2,
            output=output,
            parent=tmp_path / "parent.json",
            summary=tmp_path / "summary.json",
        )
    assert not output.exists()
