from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import pytest

from scripts.build_d6_grouped_sidecar import (
    build_sidecar,
    verify_sidecar_payload,
)
from trottercert.cubic_field import Cubic
from trottercert.hpc_artifacts import coordinate_terms_to_json, write_shard_gzip


def _write_degree_six_fixture(path: Path) -> None:
    terms = {
        ((0, 0, "X"),): Cubic.one(),
        ((0, 0, "Z"),): Cubic.one(),
        ((2, 0, "X"),): Cubic.one(),
    }
    write_shard_gzip(
        path,
        {
            "schema_version": 1,
            "kind": "issue128_exact_right_generator_degree",
            "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
            "source_commit": "fixture-commit",
            "stage_count": 31,
            "degree": 6,
            "coefficient_interval_decimal_digits": 12,
            "term_count": 3,
            "terms": coordinate_terms_to_json(terms),
            "cell_pauli_l1_upper": [3, 1],
            "site_pauli_l1_upper": [3, 4],
            "parent_sha256": "0" * 64,
        },
    )


def test_builds_deterministic_hash_bound_sidecar(tmp_path: Path) -> None:
    source = tmp_path / "d6.json.gz"
    first_path = tmp_path / "first.json"
    second_path = tmp_path / "second.json"
    _write_degree_six_fixture(source)

    sidecar = build_sidecar(source, first_path, candidate_cap=8)
    build_sidecar(source, second_path, candidate_cap=8)

    assert sidecar["kind"] == "issue128_d6_physical_channel_groups"
    assert sidecar["source_payload_sha256"] == hashlib.sha256(
        source.read_bytes()
    ).hexdigest()
    assert sidecar["term_count"] == 3
    assert Fraction(*sidecar["grouped_site_bound"]) < Fraction(
        *sidecar["site_pauli_l1_upper"]
    )
    assert sum(len(group["term_indices"]) for group in sidecar["groups"]) == 3
    assert sum(
        (Fraction(*group["bound"]) for group in sidecar["groups"]),
        Fraction(),
    ) == Fraction(*sidecar["grouped_cell_bound"])
    assert first_path.read_bytes() == second_path.read_bytes()
    assert first_path.read_bytes() == (
        json.dumps(sidecar, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    assert not tuple(tmp_path.glob(".first.json.*.tmp"))


def test_rejects_tampered_source_coefficient(tmp_path: Path) -> None:
    source = tmp_path / "d6.json.gz"
    output = tmp_path / "groups.json"
    _write_degree_six_fixture(source)
    from trottercert.hpc_artifacts import read_shard_gzip

    payload = read_shard_gzip(source)
    payload["terms"][0][1] = [[2, 1], [0, 1], [0, 1]]
    # Deliberately retain the old l1 summary.
    write_shard_gzip(source, payload)
    with pytest.raises(ValueError, match="Pauli-l1 mismatch"):
        build_sidecar(source, output, candidate_cap=8)
    assert not output.exists()


@pytest.mark.parametrize("mutation", ["duplicate", "remove"])
def test_rejects_tampered_partition_coverage(
    tmp_path: Path,
    mutation: str,
) -> None:
    source = tmp_path / "d6.json.gz"
    output = tmp_path / "groups.json"
    _write_degree_six_fixture(source)
    sidecar = build_sidecar(source, output, candidate_cap=8)
    tampered = deepcopy(sidecar)
    groups = tampered["groups"]
    if mutation == "duplicate":
        groups[0]["term_indices"].append(groups[0]["term_indices"][0])
    else:
        groups[0]["term_indices"].pop()
    with pytest.raises(ValueError, match="differs from exact regeneration"):
        verify_sidecar_payload(tampered, source)


def test_rejects_source_payload_digest_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "d6.json.gz"
    output = tmp_path / "groups.json"
    _write_degree_six_fixture(source)
    sidecar = build_sidecar(source, output, candidate_cap=8)
    tampered = deepcopy(sidecar)
    tampered["source_payload_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="differs from exact regeneration"):
        verify_sidecar_payload(tampered, source)
