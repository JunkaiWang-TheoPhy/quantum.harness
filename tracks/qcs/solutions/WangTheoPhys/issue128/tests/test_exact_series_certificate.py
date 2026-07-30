from __future__ import annotations

from fractions import Fraction

import pytest

from trottercert.cubic_field import Cubic
from trottercert.exact_series_certificate import verify_exact_monolithic_series
from trottercert.hpc_artifacts import coordinate_terms_to_json
from trottercert.intervals import cube_root_four_interval


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
