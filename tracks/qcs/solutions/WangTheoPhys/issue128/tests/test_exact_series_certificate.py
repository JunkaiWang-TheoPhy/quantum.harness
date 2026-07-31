from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import pytest

from scripts.build_d6_grouped_sidecar import build_sidecar
from trottercert.cubic_field import Cubic
from trottercert.exact_series_certificate import (
    verify_exact_degree_payload,
    verify_exact_monolithic_series,
)
from trottercert.hpc_artifacts import coordinate_terms_to_json, write_shard_gzip
from trottercert.hpc_artifacts import write_manifest_atomic
from trottercert.intervals import cube_root_four_interval
from trottercert.verify import _verify_d6_sidecar


def _payload() -> dict[str, object]:
    degree_zero = {((0, 0, "X"),): Cubic.one()}
    degree_one = {
        ((0, 0, "Y"), (1, 0, "Z")): Cubic(Fraction(1, 3), -2, 1)
    }
    root = cube_root_four_interval(12)
    cell_l1 = next(iter(degree_one.values())).enclose(root).abs_upper()
    return {
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
        "cell_pauli_l1_upper": [cell_l1.numerator, cell_l1.denominator],
        "site_pauli_l1_upper": [
            (cell_l1 / 4).numerator,
            (cell_l1 / 4).denominator,
        ],
    }


def test_exact_monolithic_series_recomputes_highest_degree_l1() -> None:
    verified = verify_exact_monolithic_series(
        _payload(),
        expected_order=1,
        expected_source_commit="abc123",
        coefficient_interval_decimal_digits=12,
    )
    assert verified.degree_term_counts == (1, 1)
    assert verified.site_l1_upper == verified.cell_l1_upper / 4


def test_exact_monolithic_series_rejects_source_or_bound_corruption() -> None:
    payload = _payload()
    with pytest.raises(ValueError, match="source commit"):
        verify_exact_monolithic_series(
            payload,
            expected_order=1,
            expected_source_commit="wrong",
            coefficient_interval_decimal_digits=12,
        )
    payload["cell_pauli_l1_upper"] = [0, 1]
    with pytest.raises(ValueError, match="cell Pauli-l1"):
        verify_exact_monolithic_series(
            payload,
            expected_order=1,
            expected_source_commit="abc123",
            coefficient_interval_decimal_digits=12,
        )


def test_exact_degree_payload_recomputes_bound_and_rejects_corruption() -> None:
    monolithic = _payload()
    payload = {
        "schema_version": 1,
        "kind": "issue128_exact_right_generator_degree",
        "formula_id": monolithic["formula_id"],
        "source_commit": monolithic["source_commit"],
        "stage_count": monolithic["stage_count"],
        "degree": 1,
        "coefficient_interval_decimal_digits": 12,
        "term_count": 1,
        "terms": monolithic["series"][1],
        "cell_pauli_l1_upper": monolithic["cell_pauli_l1_upper"],
        "site_pauli_l1_upper": monolithic["site_pauli_l1_upper"],
        "parent_sha256": "a" * 64,
    }
    verified = verify_exact_degree_payload(
        payload,
        expected_degree=1,
        expected_source_commit="abc123",
    )
    assert verified.term_count == 1
    payload["term_count"] = 2
    with pytest.raises(ValueError, match="term count"):
        verify_exact_degree_payload(
            payload,
            expected_degree=1,
            expected_source_commit="abc123",
        )
    payload["term_count"] = 1
    payload["cell_pauli_l1_upper"] = [0, 1]
    with pytest.raises(ValueError, match="cell Pauli-l1"):
        verify_exact_degree_payload(
            payload,
            expected_degree=1,
            expected_source_commit="abc123",
        )


def test_d6_main_sidecar_binding_rejects_digest_corruption(
    tmp_path: Path,
) -> None:
    monolithic = _payload()
    parent = tmp_path / "parent.json"
    parent.write_text('{"kind":"test-parent"}\n')
    parent_sha256 = hashlib.sha256(parent.read_bytes()).hexdigest()
    payload = {
        "schema_version": 1,
        "kind": "issue128_exact_right_generator_degree",
        "formula_id": monolithic["formula_id"],
        "source_commit": monolithic["source_commit"],
        "stage_count": monolithic["stage_count"],
        "degree": 6,
        "coefficient_interval_decimal_digits": 12,
        "term_count": 1,
        "terms": monolithic["series"][1],
        "cell_pauli_l1_upper": monolithic["cell_pauli_l1_upper"],
        "site_pauli_l1_upper": monolithic["site_pauli_l1_upper"],
        "parent_sha256": parent_sha256,
    }
    sidecar = tmp_path / "d6.json.gz"
    main = tmp_path / "main.json"
    write_shard_gzip(sidecar, payload)
    main.write_text("{}")
    raw = sidecar.read_bytes()
    candidate = {
        "d6_certificate": {
            "path": sidecar.name,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "parent_sha256": parent_sha256,
            "parent_path": parent.name,
            "source_commit": "abc123",
            "coefficient_interval_decimal_digits": 12,
            "term_count": 1,
            "cell_norm_upper": monolithic["cell_pauli_l1_upper"],
            "site_norm_upper": monolithic["site_pauli_l1_upper"],
        }
    }
    verified = _verify_d6_sidecar(main, candidate)
    assert verified.artifact.term_count == 1
    parent.write_text('{"kind":"tampered"}\n')
    with pytest.raises(ValueError, match="D6 parent artifact digest"):
        _verify_d6_sidecar(main, candidate)
    candidate["d6_certificate"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="D6 sidecar digest"):
        _verify_d6_sidecar(main, candidate)


def test_d6_grouped_sidecar_regenerates_bound_and_rejects_tampering(
    tmp_path: Path,
) -> None:
    terms = {
        ((0, 0, "X"),): Cubic.one(),
        ((0, 0, "Z"),): Cubic.one(),
        ((2, 0, "X"),): Cubic.one(),
    }
    root = cube_root_four_interval(12)
    cell_l1 = sum(
        (coefficient.enclose(root).abs_upper() for coefficient in terms.values()),
        Fraction(),
    )
    parent = tmp_path / "parent.json"
    parent.write_text('{"kind":"test-parent"}\n')
    parent_sha256 = hashlib.sha256(parent.read_bytes()).hexdigest()
    payload = {
        "schema_version": 1,
        "kind": "issue128_exact_right_generator_degree",
        "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
        "source_commit": "abc123",
        "stage_count": 31,
        "degree": 6,
        "coefficient_interval_decimal_digits": 12,
        "term_count": len(terms),
        "terms": coordinate_terms_to_json(terms),
        "cell_pauli_l1_upper": [cell_l1.numerator, cell_l1.denominator],
        "site_pauli_l1_upper": [
            (cell_l1 / 4).numerator,
            (cell_l1 / 4).denominator,
        ],
        "parent_sha256": parent_sha256,
    }
    sidecar = tmp_path / "d6.json.gz"
    groups = tmp_path / "d6-groups.json"
    comparison = tmp_path / "d6-comparison.json"
    main = tmp_path / "main.json"
    write_shard_gzip(sidecar, payload)
    grouped = build_sidecar(sidecar, groups, candidate_cap=8)
    sidecar_sha256 = hashlib.sha256(sidecar.read_bytes()).hexdigest()
    write_manifest_atomic(
        comparison,
        {
            "schema_version": 1,
            "kind": "issue128_exact_degree_payload_comparison",
            "status": "complete",
            "degree": 6,
            "term_count": len(terms),
            "exact_coefficient_map_equal": True,
            "cell_pauli_l1_upper": payload["cell_pauli_l1_upper"],
            "site_pauli_l1_upper": payload["site_pauli_l1_upper"],
            "left": {
                "path": "left.json.gz",
                "sha256": sidecar_sha256,
                "source_commit": "abc123",
                "parent_sha256": parent_sha256,
            },
            "right": {
                "path": "right.json.gz",
                "sha256": "f" * 64,
                "source_commit": "abc123",
                "parent_sha256": "e" * 64,
            },
        },
    )
    main.write_text("{}")
    candidate = {
        "d6_certificate": {
            "path": sidecar.name,
            "sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
            "parent_sha256": parent_sha256,
            "parent_path": parent.name,
            "source_commit": "abc123",
            "coefficient_interval_decimal_digits": 12,
            "term_count": len(terms),
            "groups_path": groups.name,
            "groups_sha256": hashlib.sha256(groups.read_bytes()).hexdigest(),
            "comparison_path": comparison.name,
            "comparison_sha256": hashlib.sha256(
                comparison.read_bytes()
            ).hexdigest(),
            "group_count": grouped["group_count"],
            "max_group_size": grouped["max_group_size"],
            "l1_site_norm_upper": grouped["site_pauli_l1_upper"],
            "cell_norm_upper": grouped["grouped_cell_bound"],
            "site_norm_upper": grouped["grouped_site_bound"],
        }
    }
    verified = _verify_d6_sidecar(main, candidate)
    assert verified.site_bound == Fraction(*grouped["grouped_site_bound"])
    assert verified.site_bound < cell_l1 / 4
    assert verified.group_count == grouped["group_count"]

    tampered = json.loads(groups.read_text())
    tampered["groups"][0]["bound"] = [0, 1]
    groups.write_text(
        json.dumps(tampered, sort_keys=True, separators=(",", ":")) + "\n"
    )
    candidate["d6_certificate"]["groups_sha256"] = hashlib.sha256(
        groups.read_bytes()
    ).hexdigest()
    with pytest.raises(ValueError, match="group bound mismatch"):
        _verify_d6_sidecar(main, candidate)

    build_sidecar(sidecar, groups, candidate_cap=8)
    candidate["d6_certificate"]["groups_sha256"] = hashlib.sha256(
        groups.read_bytes()
    ).hexdigest()
    comparison_tampered = json.loads(comparison.read_text())
    comparison_tampered["exact_coefficient_map_equal"] = False
    comparison.write_text(
        json.dumps(
            comparison_tampered,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    )
    candidate["d6_certificate"]["comparison_sha256"] = hashlib.sha256(
        comparison.read_bytes()
    ).hexdigest()
    with pytest.raises(ValueError, match="comparison audit status"):
        _verify_d6_sidecar(main, candidate)
