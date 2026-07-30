from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import pytest

from scripts.compare_exact_degree_payloads import compare
from trottercert.cubic_field import Cubic
from trottercert.hpc_artifacts import coordinate_terms_to_json, write_shard_gzip


def _payload(coefficient: Cubic) -> dict[str, object]:
    terms = {((0, 0, "X"),): coefficient}
    return {
        "schema_version": 1,
        "kind": "issue128_exact_right_generator_degree",
        "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
        "source_commit": "fixture",
        "stage_count": 31,
        "degree": 6,
        "term_count": 1,
        "terms": coordinate_terms_to_json(terms),
        "cell_pauli_l1_upper": [3, 2],
        "site_pauli_l1_upper": [3, 8],
        "parent_sha256": "0" * 64,
    }


def test_compare_exact_degree_payloads_records_exact_agreement(tmp_path: Path) -> None:
    left = tmp_path / "left.json.gz"
    right = tmp_path / "right.json.gz"
    write_shard_gzip(left, _payload(Cubic(Fraction(2, 3), 1, 0)))
    write_shard_gzip(right, _payload(Cubic(Fraction(2, 3), 1, 0)))

    result = compare(left, right, degree=6, output=tmp_path / "comparison.json")

    assert result["exact_coefficient_map_equal"] is True
    assert result["term_count"] == 1
    assert result["left"]["sha256"] == result["right"]["sha256"]


def test_compare_exact_degree_payloads_rejects_coefficient_mismatch(tmp_path: Path) -> None:
    left = tmp_path / "left.json.gz"
    right = tmp_path / "right.json.gz"
    write_shard_gzip(left, _payload(Cubic.one()))
    write_shard_gzip(right, _payload(Cubic(2, 0, 0)))

    with pytest.raises(ArithmeticError, match="coefficient maps differ"):
        compare(left, right, degree=6, output=tmp_path / "comparison.json")
