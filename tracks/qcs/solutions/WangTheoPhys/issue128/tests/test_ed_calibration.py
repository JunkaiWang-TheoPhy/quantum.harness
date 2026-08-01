from __future__ import annotations

import hashlib
import json
import math
import re
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

from scripts.run_ed_calibration import (
    FROZEN_SIZES,
    FROZEN_STEPS,
    SCIENTIFIC_DEPENDENCY_PATHS,
    calibrate_rectangles,
)
from trottercert.crosscheck import (
    fourth_order_product_step,
    open_rectangle_dense_fragments,
    open_rectangle_published_triangle_certificate,
    operator_error,
)


def test_calibration_is_deterministic() -> None:
    left = calibrate_rectangles(((2, 3),), (8, 12))
    right = calibrate_rectangles(((2, 3),), (8, 12))
    assert left == right


def test_certificate_dominates_actual_error() -> None:
    data = calibrate_rectangles(((2, 3),), (16,))
    row = data["rows"][0]
    assert Decimal(row["actual_operator_error"]) <= Decimal(
        row["certified_upper"]
    )


def test_frozen_grid_and_row_order() -> None:
    assert FROZEN_SIZES == ((2, 3), (2, 4), (3, 3))
    assert FROZEN_STEPS == (8, 12, 16, 24, 32, 48, 64, 96)
    data = calibrate_rectangles(((2, 3),), (12, 8))
    assert [row["steps"] for row in data["rows"]] == [8, 12]


def test_certified_decimal_fields_are_outward_and_actual_is_diagnostic() -> None:
    row = calibrate_rectangles(((2, 3),), (8,))["rows"][0]
    assert re.fullmatch(r"\d\.\d{8}e[+-]\d{2,3}", row["actual_operator_error"])
    assert "diagnostic_not_a_certified_interval" in row[
        "actual_operator_error_status"
    ]
    for field in ("certified_upper", "published_upper"):
        exact = Fraction(*row[f"{field}_exact"])
        serialized = Fraction(Decimal(row[field]))
        assert serialized >= exact, (field, serialized, exact)


def test_calibration_records_four_matching_normalization() -> None:
    row = calibrate_rectangles(((2, 3),), (8,))["rows"][0]
    assert row["model"] == "open_rectangular_spin_half_isotropic_heisenberg"
    assert row["normalization"] == "(XX+YY+ZZ)/4"
    assert row["fragment_count"] == 4
    assert row["unique_bond_count"] == 7
    assert row["boundary"] == "open"


def test_provenance_has_no_commit_identity_and_binds_scientific_inputs() -> None:
    data = calibrate_rectangles(((2, 3),), (8,))
    assert "source_commit" not in data
    assert all("source_commit" not in row for row in data["rows"])
    provenance = data["provenance"]
    assert provenance["scientific_dependency_sha256"] == {
        relative: hashlib.sha256(
            (Path(__file__).resolve().parents[1] / relative).read_bytes()
        ).hexdigest()
        for relative in SCIENTIFIC_DEPENDENCY_PATHS
    }
    reference = provenance["reference_certificate"]
    assert reference["verification"]["valid"] is True
    assert reference["verification"]["coefficient_generation_replayed"] is False
    assert "not transferred to open boundaries" in reference["role"]
    assert set(provenance["runtime_environment"]) >= {
        "python",
        "numpy",
        "scipy",
    }


def test_reusable_operator_error_matches_calibration() -> None:
    _groups, fragments, hamiltonian = open_rectangle_dense_fragments(2, 3)
    product_step = fourth_order_product_step(fragments, time=1.0, steps=8)
    error = operator_error(hamiltonian, product_step, time=1.0, steps=8)
    row = calibrate_rectangles(((2, 3),), (8,))["rows"][0]
    assert math.isclose(
        error,
        float(row["actual_operator_error"]),
        rel_tol=5e-9,
        abs_tol=5e-16,
    )


def test_open_certificate_is_rebuilt_from_finite_pauli_commutators() -> None:
    data = calibrate_rectangles(((2, 3),), (8,))
    record = data["open_rectangle_certificates"][0]
    rebuilt = open_rectangle_published_triangle_certificate(2, 3)
    assert Fraction(*record["constant_upper_exact"]) == rebuilt.constant_upper
    assert Fraction(Decimal(record["constant_upper"])) >= rebuilt.constant_upper
    assert record["theorem_terms"] == rebuilt.theorem_terms == 8_316
    assert record["expanded_commutator_keys"] == 1_024
    assert record["norm_method"] == "direct_finite_open_exact_Pauli_axis_l1"
    assert data["rows"][0]["open_triangle_certificate_digest"] == record[
        "identity_digest"
    ]
    unsigned = {
        key: value for key, value in record.items() if key != "identity_digest"
    }
    encoded = json.dumps(
        unsigned, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    assert hashlib.sha256(encoded).hexdigest() == record["identity_digest"]


def test_frozen_artifact_covers_grid_and_is_dominated() -> None:
    root = Path(__file__).resolve().parents[1]
    payload = json.loads(
        (root / "benchmarks/paper-a/ed-calibration.json").read_text(
            encoding="utf-8"
        )
    )
    rows = payload["rows"]
    assert len(rows) == len(FROZEN_SIZES) * len(FROZEN_STEPS) == 24
    assert [
        (row["width"], row["height"], row["steps"]) for row in rows
    ] == sorted(
        (width, height, steps)
        for width, height in FROZEN_SIZES
        for steps in FROZEN_STEPS
    )
    assert all(
        Decimal(row["actual_operator_error"])
        <= Decimal(row["certified_upper"])
        for row in rows
    )
    assert all(
        row["candidate_grouped_open_boundary_transfer_claimed"] is False
        for row in rows
    )
    assert all(
        row["certified_upper_method"].startswith("direct_finite_open_")
        for row in rows
    )
    for record in payload["open_rectangle_certificates"]:
        rebuilt = open_rectangle_published_triangle_certificate(
            record["width"], record["height"]
        )
        assert Fraction(*record["constant_upper_exact"]) == (
            rebuilt.constant_upper
        )
    for row in rows:
        record = next(
            certificate
            for certificate in payload["open_rectangle_certificates"]
            if certificate["identity_digest"]
            == row["open_triangle_certificate_digest"]
        )
        expected_published = (
            Fraction(*record["constant_upper_exact"]) / row["steps"] ** 4
        )
        assert Fraction(*row["published_upper_exact"]) == expected_published
        assert Fraction(*row["certified_upper_exact"]) == min(
            Fraction(2), expected_published
        )
        for field in ("certified_upper", "published_upper"):
            assert Fraction(Decimal(row[field])) >= Fraction(
                *row[f"{field}_exact"]
            )
    assert payload["provenance"]["scientific_dependency_sha256"] == {
        relative: hashlib.sha256((root / relative).read_bytes()).hexdigest()
        for relative in SCIENTIFIC_DEPENDENCY_PATHS
    }
    assert payload["provenance"]["reference_certificate"]["verification"][
        "valid"
    ] is True
    reference = payload["provenance"]["reference_certificate"]
    assert hashlib.sha256((root / reference["path"]).read_bytes()).hexdigest() == (
        reference["sha256"]
    )
    provenance_encoded = json.dumps(
        payload["provenance"],
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    assert hashlib.sha256(provenance_encoded).hexdigest() == payload[
        "provenance_digest"
    ]
