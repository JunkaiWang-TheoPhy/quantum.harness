from __future__ import annotations

import json
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import pytest

from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.dual_log_pairing import contract_log_degree_shard
from trottercert.dual_manifest_pairing import contract_word_manifest
from trottercert.dual_word_manifest import build_word_manifests

ROOT = Path(__file__).resolve().parents[1]
EXTENSIVE = (
    ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
SOURCES = {"source.py": "a" * 64}


def _cubic(value: list[list[int]]) -> Cubic:
    return Cubic(*(Fraction(*coordinate) for coordinate in value))


def _degree_five_manifest():
    return build_word_manifests(
        fourth_order_suzuki_cubic_stages(4),
        degree=5,
        shard_count=1,
        implementation_sources=SOURCES,
    )[0]


def test_manifest_pairing_is_exactly_additive_across_groups() -> None:
    manifest = _degree_five_manifest()
    selected = replace(manifest, groups=manifest.groups[:2])
    together = contract_word_manifest(selected, length=6)
    separate = [
        contract_word_manifest(replace(manifest, groups=(group,)), length=6)
        for group in selected.groups
    ]

    assert together.tau_h == sum(
        (partial.tau_h for partial in separate), Cubic.zero()
    )
    assert together.tau_h2 == sum(
        (partial.tau_h2 for partial in separate), Cubic.zero()
    )
    assert together.tau_w == sum(
        (partial.tau_w for partial in separate), Cubic.zero()
    )
    assert together.word_count == sum(item.word_count for item in separate)
    assert together.retained_term_count == sum(
        item.retained_term_count for item in separate
    )


def test_manifest_pairing_rejects_word_length_degree_mismatch() -> None:
    manifest = _degree_five_manifest()
    with pytest.raises(ValueError, match="degree|length"):
        contract_word_manifest(replace(manifest, degree=7), length=6)


@pytest.mark.slow
def test_degree_five_manifest_pairing_matches_frozen_record() -> None:
    manifests = build_word_manifests(
        fourth_order_suzuki_cubic_stages(4),
        degree=5,
        shard_count=4,
        implementation_sources=SOURCES,
    )
    observed = [contract_word_manifest(manifest, length=6) for manifest in manifests]
    tau_h = sum((item.tau_h for item in observed), Cubic.zero())
    tau_h2 = sum((item.tau_h2 for item in observed), Cubic.zero())
    tau_w = sum((item.tau_w for item in observed), Cubic.zero())
    expected = next(
        record
        for record in json.loads(EXTENSIVE.read_text())["records"]
        if record["length"] == 6
    )

    assert tau_h == _cubic(expected["moments"]["tau_h_e5"])
    assert tau_h2 == _cubic(expected["moments"]["tau_h2_e5"])
    assert tau_w == _cubic(expected["moments"]["tau_w_e5"])


@pytest.mark.slow
def test_one_degree_seven_group_matches_frozen_in_memory_path() -> None:
    stages = fourth_order_suzuki_cubic_stages(4)
    shard_index = 123
    shard_count = 4096
    manifest = build_word_manifests(
        stages,
        degree=7,
        shard_count=shard_count,
        implementation_sources=SOURCES,
    )[shard_index]
    selected = replace(manifest, groups=manifest.groups[:1])
    group_ordinal = selected.groups[0].ordinal

    observed = contract_word_manifest(selected, length=6)
    expected = contract_log_degree_shard(
        stages,
        7,
        group_ordinal,
        shard_count,
        length=6,
    )

    assert observed.group_indices == expected.group_indices
    assert observed.word_count == expected.word_count
    assert observed.nonzero_word_count == expected.nonzero_word_count
    assert observed.retained_term_count == expected.retained_term_count
    assert observed.tau_h == expected.tau_h
    assert observed.tau_h2 == expected.tau_h2
    assert observed.tau_w == expected.tau_w
