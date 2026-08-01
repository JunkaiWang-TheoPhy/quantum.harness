from __future__ import annotations

from fractions import Fraction

import pytest

from trottercert.dual_log_tail import (
    average_stage_tail,
    pointwise_stage_tail,
    verify_analytic_caps,
)
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
