from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from scripts.extract_exact_degree import main
from trottercert.cubic_field import Cubic
from trottercert.hpc_artifacts import (
    coordinate_terms_to_json,
    read_shard_gzip,
    sha256_file,
    write_shard_gzip,
)


def test_extract_exact_degree_binds_parent_and_recomputes_l1(
    tmp_path: Path,
    monkeypatch,
) -> None:
    parent = tmp_path / "parent.json.gz"
    output = tmp_path / "d1.json.gz"
    manifest = tmp_path / "manifest.json"
    degree_zero = {((0, 0, "X"),): Cubic.one()}
    degree_one = {((0, 0, "Y"),): Cubic(Fraction(3, 5), 0, 0)}
    write_shard_gzip(
        parent,
        {
            "schema_version": 1,
            "kind": "issue128_exact_right_generator_monolithic",
            "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
            "source_commit": "abc123",
            "stage_count": 31,
            "order": 1,
            "series": [
                coordinate_terms_to_json(degree_zero),
                coordinate_terms_to_json(degree_one),
            ],
            "cell_pauli_l1_upper": [3, 5],
            "site_pauli_l1_upper": [3, 20],
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "extract_exact_degree.py",
            "--input",
            str(parent),
            "--degree",
            "1",
            "--decimal-digits",
            "12",
            "--output",
            str(output),
            "--manifest",
            str(manifest),
        ],
    )
    main()
    payload = read_shard_gzip(output)
    record = json.loads(manifest.read_text())
    assert payload["term_count"] == 1
    assert payload["cell_pauli_l1_upper"] == [3, 5]
    assert payload["parent_sha256"] == sha256_file(parent)
    assert record["status"] == "complete"
    assert record["output_sha256"] == sha256_file(output)
