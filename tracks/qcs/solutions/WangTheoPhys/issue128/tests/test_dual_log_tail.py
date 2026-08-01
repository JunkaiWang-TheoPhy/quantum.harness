from __future__ import annotations

from dataclasses import replace
from fractions import Fraction

import pytest

from trottercert.dual_log_tail import (
    average_stage_tail,
    build_dual_tail_envelope,
    issue128_dual_tail_geometry,
    pointwise_stage_tail,
    verify_analytic_caps,
)
from trottercert.refined_error import build_refined_fourth_order_constants
from trottercert.rigorous_fourth import (
    fourth_order_suzuki_interval_stages,
)


def test_geometric_stage_tail_matches_long_finite_sum() -> None:
    stages, _ = fourth_order_suzuki_interval_stages(
        4, decimal_digits=30
    )
    steps = 97
    first = 10
    prefix = Fraction()
    average_reference = Fraction()
    pointwise_reference = Fraction()
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        average_reference += coefficient * sum(
            (ratio**degree for degree in range(first, 160)),
            Fraction(),
        )
        pointwise_reference += coefficient * sum(
            (
                (degree + 1) * ratio**degree
                for degree in range(first, 160)
            ),
            Fraction(),
        )
        prefix += coefficient
    average = average_stage_tail(stages, steps, first)
    pointwise = pointwise_stage_tail(stages, steps, first)
    assert average > average_reference
    assert pointwise > pointwise_reference
    assert average - average_reference < Fraction(1, 10**100)
    assert pointwise - pointwise_reference < Fraction(1, 10**98)


def test_stage_tail_rejects_invalid_convergence_region() -> None:
    stages, _ = fourth_order_suzuki_interval_stages(
        4, decimal_digits=30
    )
    with pytest.raises(ValueError, match="convergence"):
        average_stage_tail(stages, 1, 10)
    with pytest.raises(ValueError, match="omitted degree"):
        pointwise_stage_tail(stages, 97, -1)


def test_selected_analytic_caps_are_proved_by_rational_gates() -> None:
    verify_analytic_caps(
        log_radius_cap=Fraction(3, 2),
        dexp_cap=Fraction(8, 5),
        multiplier_difference_cap=Fraction(1),
        even_series_cap=Fraction(10, 33),
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("log_radius_cap", Fraction(151, 100), "log radius"),
        ("dexp_cap", Fraction(159, 100), "dexp cap"),
        (
            "multiplier_difference_cap",
            Fraction(99, 100),
            "multiplier",
        ),
        ("even_series_cap", Fraction(3, 10), "even-series"),
    ],
)
def test_analytic_caps_reject_unsafe_claims(
    field: str,
    value: Fraction,
    message: str,
) -> None:
    arguments = {
        "log_radius_cap": Fraction(3, 2),
        "dexp_cap": Fraction(8, 5),
        "multiplier_difference_cap": Fraction(1),
        "even_series_cap": Fraction(10, 33),
    }
    arguments[field] = value
    with pytest.raises(ValueError, match=message):
        verify_analytic_caps(**arguments)


def test_issue128_geometry_has_exact_moment_caps() -> None:
    geometry = issue128_dual_tail_geometry()
    assert geometry.n_sites == 144
    assert geometry.cells == 36
    assert geometry.steps == 97
    assert geometry.h_hs_squared == 54
    assert geometry.w_hs_squared == Fraction(23085, 4)
    assert geometry.h_hs_cap**2 >= geometry.h_hs_squared
    assert (
        geometry.w_hs_per_cell_cap**2
        >= geometry.w_hs_squared / geometry.cells**2
    )
    assert geometry.max_w_pauli_coefficient == Fraction(1, 4)


def test_fixed_envelope_regression_and_branch_gate() -> None:
    constants = build_refined_fourth_order_constants(
        decimal_digits=30,
        quantization_digits=24,
    )
    result = build_dual_tail_envelope(
        constants, issue128_dual_tail_geometry()
    )
    assert result.centered_exact_phase_radius < Fraction(3, 2)
    assert result.centered_log_radius_bound < Fraction(3, 2)
    assert Fraction(2146, 10**9) < result.average_generator_defect
    assert result.average_generator_defect < Fraction(2148, 10**9)
    assert Fraction(1187, 10**8) < result.pointwise_generator_defect
    assert result.pointwise_generator_defect < Fraction(1188, 10**8)
    assert Fraction(354, 10**10) < result.log_defect
    assert result.log_defect < Fraction(355, 10**10)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("n_sites", 143, "site count"),
        ("cells", 35, "cell count"),
        ("steps", 96, "step count"),
        ("h_hs_squared", Fraction(53), "H moment"),
        ("w_hs_squared", Fraction(23084, 4), "W moment"),
        (
            "max_w_pauli_coefficient",
            Fraction(1, 5),
            "W Pauli",
        ),
    ],
)
def test_geometry_mutations_fail_closed(
    field: str,
    value: object,
    message: str,
) -> None:
    constants = build_refined_fourth_order_constants(
        decimal_digits=30,
        quantization_digits=24,
    )
    forged = replace(issue128_dual_tail_geometry(), **{field: value})
    with pytest.raises(ValueError, match=message):
        build_dual_tail_envelope(constants, forged)
