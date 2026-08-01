from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

from .rigorous_fourth import IntervalStage

LOG_RADIUS_CAP = Fraction(3, 2)
DEXP_CAP = Fraction(8, 5)
MULTIPLIER_DIFFERENCE_CAP = Fraction(1)
EVEN_SERIES_CAP = Fraction(10, 33)


def _validate_tail_inputs(
    steps: int,
    first_omitted_degree: int,
) -> None:
    if not isinstance(steps, int) or isinstance(steps, bool) or steps < 1:
        raise ValueError("steps must be a positive integer")
    if (
        not isinstance(first_omitted_degree, int)
        or isinstance(first_omitted_degree, bool)
        or first_omitted_degree < 0
    ):
        raise ValueError("first omitted degree must be nonnegative")


def average_stage_tail(
    stages: Sequence[IntervalStage],
    steps: int,
    first_omitted_degree: int,
) -> Fraction:
    """Average-in-time absolute conjugation tail without its density factor."""

    _validate_tail_inputs(steps, first_omitted_degree)
    prefix = Fraction()
    total = Fraction()
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("stage tail is outside its convergence region")
        total += coefficient * ratio**first_omitted_degree / (1 - ratio)
        prefix += coefficient
    return total


def pointwise_stage_tail(
    stages: Sequence[IntervalStage],
    steps: int,
    first_omitted_degree: int,
) -> Fraction:
    """Endpoint absolute conjugation tail without its density factor."""

    _validate_tail_inputs(steps, first_omitted_degree)
    prefix = Fraction()
    total = Fraction()
    degree = first_omitted_degree
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("stage tail is outside its convergence region")
        total += (
            coefficient
            * ratio**degree
            * ((degree + 1) - degree * ratio)
            / (1 - ratio) ** 2
        )
        prefix += coefficient
    return total


def verify_analytic_caps(
    *,
    log_radius_cap: Fraction,
    dexp_cap: Fraction,
    multiplier_difference_cap: Fraction,
    even_series_cap: Fraction,
) -> None:
    """Verify the rational scalar lemmas used by the ``dexp`` conversion."""

    radius = Fraction(log_radius_cap)
    if radius < 0 or radius > LOG_RADIUS_CAP:
        raise ValueError("log radius exceeds the proved interval")

    # sin(x)/x >= 1-x^2/6 on [0, 3/2].
    required_dexp = 1 / (1 - radius**2 / 6)
    if Fraction(dexp_cap) < required_dexp:
        raise ValueError("dexp cap is below the rational sine bound")

    # For y=2x and |x|<=3/2, the cotangent remainder obeys
    # |x*cot(x)-1|/|x| <= x/(3*(1-x^2/6)) <= 4/5.  Therefore
    # |(f(iy)-1)/(iy)|^2 <= 1/4+(2/5)^2 < 1.
    if Fraction(multiplier_difference_cap) < MULTIPLIER_DIFFERENCE_CAP:
        raise ValueError(
            "multiplier-difference cap is below the proved cap"
        )

    # The k=1 term is 1/3! and every later ratio is at most 9/20.
    required_even_series = Fraction(1, 6) / (1 - Fraction(9, 20))
    if Fraction(even_series_cap) < required_even_series:
        raise ValueError("even-series cap is below the ratio majorant")
