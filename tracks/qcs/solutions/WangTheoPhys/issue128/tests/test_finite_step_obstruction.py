from __future__ import annotations

from fractions import Fraction

import pytest

from trottercert.cubic_field import Cubic
from trottercert.finite_step_obstruction import decide_finite_step
from trottercert.intervals import RationalInterval, cube_root_four_interval


ROOT = cube_root_four_interval(36)


def _point(value: int | Fraction) -> RationalInterval:
    return RationalInterval.point(value)


def test_positive_margin_certifies_obstruction() -> None:
    decision = decide_finite_step(
        Cubic(10, 0, 0),
        Cubic.zero(),
        Cubic.zero(),
        ROOT,
        _point(1),
        Fraction(1),
    )

    assert decision.leading_interval == _point(10)
    assert decision.signed_margin_interval == _point(9)
    assert decision.status == "certified_obstruction"


def test_negative_leading_interval_uses_sound_absolute_lower_bound() -> None:
    decision = decide_finite_step(
        Cubic(-10, 0, 0),
        Cubic.zero(),
        Cubic.zero(),
        ROOT,
        _point(3),
        Fraction(1),
    )

    assert decision.absolute_leading_interval == _point(10)
    assert decision.signed_margin_interval == _point(7)
    assert decision.status == "certified_obstruction"


def test_crossing_zero_stays_inconclusive() -> None:
    decision = decide_finite_step(
        Cubic(1, -1, 0),
        Cubic.zero(),
        Cubic.zero(),
        RationalInterval(Fraction(1), Fraction(2)),
        _point(Fraction(1, 4)),
        Fraction(1),
    )

    assert decision.leading_interval == RationalInterval(Fraction(-1), Fraction(0))
    assert decision.signed_margin_interval == RationalInterval(
        Fraction(-1, 4), Fraction(3, 4)
    )
    assert decision.status == "inconclusive"


def test_strictly_negative_margin_records_failed_bound() -> None:
    decision = decide_finite_step(
        Cubic(1, 0, 0),
        Cubic.zero(),
        Cubic.zero(),
        ROOT,
        _point(2),
        Fraction(1),
    )

    assert decision.signed_margin_interval == _point(-1)
    assert decision.status == "certified_no_margin"


def test_degree_scaling_matches_e9_reducer_convention() -> None:
    decision = decide_finite_step(
        Cubic(1, 0, 0),
        Cubic(1, 0, 0),
        Cubic(1, 0, 0),
        ROOT,
        _point(0),
        Fraction(1, 97),
    )

    expected = Fraction(1, 97**4) + Fraction(1, 97**6) + Fraction(1, 97**8)
    assert decision.leading_interval == _point(expected)


@pytest.mark.parametrize(
    ("root", "tail", "h", "message"),
    (
        (
            RationalInterval(Fraction(1), Fraction(3, 2)),
            _point(0),
            Fraction(1),
            "cube root",
        ),
        (ROOT, RationalInterval(Fraction(-1), Fraction(1)), Fraction(1), "nonnegative"),
        (ROOT, _point(0), Fraction(0), "positive"),
    ),
)
def test_invalid_intervals_and_step_are_rejected(
    root: RationalInterval,
    tail: RationalInterval,
    h: Fraction,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        decide_finite_step(Cubic.zero(), Cubic.zero(), Cubic.zero(), root, tail, h)


def test_boolean_step_is_rejected() -> None:
    with pytest.raises(TypeError, match="positive rational"):
        decide_finite_step(
            Cubic.zero(), Cubic.zero(), Cubic.zero(), ROOT, _point(0), True
        )
