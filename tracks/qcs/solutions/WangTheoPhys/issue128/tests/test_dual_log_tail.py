from __future__ import annotations

from dataclasses import replace
from fractions import Fraction

import pytest

from trottercert.dual_log_tail import (
    average_stage_tail,
    build_dual_tail_envelope,
    certify_issue128_dual_log_tail,
    issue128_dual_tail_geometry,
    issue128_generator_constants,
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
    constants = issue128_generator_constants()
    result = build_dual_tail_envelope(
        constants, issue128_dual_tail_geometry()
    )
    assert result.centered_exact_phase_radius < Fraction(3, 2)
    assert result.centered_log_radius_bound < Fraction(3, 2)
    assert Fraction(8889, 10**10) < result.average_generator_defect
    assert result.average_generator_defect < Fraction(8891, 10**10)
    assert Fraction(5478, 10**9) < result.pointwise_generator_defect
    assert result.pointwise_generator_defect < Fraction(5480, 10**9)
    assert Fraction(1466, 10**11) < result.log_defect
    assert result.log_defect < Fraction(1467, 10**11)


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
    constants = issue128_generator_constants()
    forged = replace(issue128_dual_tail_geometry(), **{field: value})
    with pytest.raises(ValueError, match=message):
        build_dual_tail_envelope(constants, forged)


def test_complete_e11_tail_is_compatible_and_below_preaudit_gate() -> None:
    result = certify_issue128_dual_log_tail()
    assert Fraction(6448, 10**15) < result.direct_even_generator_tail
    assert result.direct_even_generator_tail < Fraction(6450, 10**15)
    assert Fraction(46, 10**15) < result.dexp_correction_tail
    assert result.dexp_correction_tail < Fraction(47, 10**15)
    assert result.total_log_tail == (
        result.direct_even_generator_tail + result.dexp_correction_tail
    )
    assert Fraction(6495, 10**15) < result.total_log_tail
    assert result.total_log_tail < Fraction(6496, 10**15)
    assert result.first_omitted_log_degree == 11
    assert result.claim_scope == "fixed_12x12_r97_dual_local_log"


def test_direct_tail_has_the_exact_dual_normalization() -> None:
    result = certify_issue128_dual_log_tail()
    constants = issue128_generator_constants()
    raw_cell = Fraction(3, 2) * average_stage_tail(
        constants.stages,
        97,
        10,
    )
    assert result.direct_even_generator_tail == raw_cell * Fraction(1, 4)


def test_envelope_decreases_at_larger_step_counts() -> None:
    constants = issue128_generator_constants()
    base = issue128_dual_tail_geometry()
    at_97 = build_dual_tail_envelope(constants, base)
    at_98 = build_dual_tail_envelope(
        constants,
        replace(base, steps=98),
    )
    assert at_98.average_generator_defect < at_97.average_generator_defect
    assert at_98.log_defect < at_97.log_defect
