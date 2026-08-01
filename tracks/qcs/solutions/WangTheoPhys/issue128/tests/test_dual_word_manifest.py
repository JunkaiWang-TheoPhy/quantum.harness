from __future__ import annotations

import json
from pathlib import Path

import pytest

from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.dual_word_manifest import (
    build_word_manifests,
    load_manifest_index,
    load_word_manifest,
    verify_manifest_index,
    write_manifest_set,
)


def test_manifest_group_assignment_is_canonical_and_complete() -> None:
    manifests = build_word_manifests(
        fourth_order_suzuki_cubic_stages(4),
        degree=5,
        shard_count=3,
        implementation_sources={"source.py": "a" * 64},
    )
    groups = [group for manifest in manifests for group in manifest.groups]
    ordered = sorted(groups, key=lambda item: item.ordinal)

    assert [group.ordinal for group in ordered] == list(
        range(manifests[0].total_groups)
    )
    assert all(
        group.ordinal % 3 == manifest.shard_index
        for manifest in manifests
        for group in manifest.groups
    )
    assert all(
        tuple(sorted(group.records, key=lambda item: item.word))
        == group.records
        for group in groups
    )
    assert sum(len(group.records) for group in groups) == 1020


def test_manifest_set_round_trips_and_binds_children(tmp_path: Path) -> None:
    output = tmp_path / "manifests"
    index_path = write_manifest_set(
        output,
        fourth_order_suzuki_cubic_stages(4),
        degree=5,
        shard_count=2,
        implementation_sources={"source.py": "a" * 64},
    )
    index = load_manifest_index(index_path)

    verify_manifest_index(index, output)
    assert index.total_groups == 256
    assert index.total_words == 1020
    assert [record.path for record in index.shards] == [
        "manifest-000.json",
        "manifest-001.json",
    ]

    shard_path = output / index.shards[0].path
    payload = json.loads(shard_path.read_text())
    payload["groups"][0]["records"][0]["word"] = "999"
    shard_path.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    )
    with pytest.raises(ValueError, match="digest"):
        load_word_manifest(shard_path, index)
    with pytest.raises(ValueError, match="digest"):
        verify_manifest_index(index, output)


def test_manifest_publication_refuses_existing_directory(tmp_path: Path) -> None:
    output = tmp_path / "manifests"
    output.mkdir()
    with pytest.raises(FileExistsError):
        write_manifest_set(
            output,
            fourth_order_suzuki_cubic_stages(4),
            degree=3,
            shard_count=2,
            implementation_sources={"source.py": "a" * 64},
        )


@pytest.mark.parametrize("degree", [0, 2, 4])
def test_manifest_rejects_invalid_degree(degree: int) -> None:
    with pytest.raises(ValueError, match="odd"):
        build_word_manifests(
            fourth_order_suzuki_cubic_stages(4),
            degree=degree,
            shard_count=1,
            implementation_sources={"source.py": "a" * 64},
        )
