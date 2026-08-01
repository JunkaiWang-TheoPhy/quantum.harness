from __future__ import annotations

from fractions import Fraction

import pytest
import sympy as sp

from trottercert.cubic_field import Cubic
from trottercert.finite_step_obstruction import (
    decide_affine_spectral_obstruction,
    decide_finite_step,
)
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


def _centered_moments(values: tuple[Fraction, ...]) -> tuple[Fraction, Fraction]:
    mean = sum(values, Fraction()) / len(values)
    centered = tuple(value - mean for value in values)
    return (
        sum((value**2 for value in centered), Fraction()) / len(values),
        sum((value**3 for value in centered), Fraction()) / len(values),
    )


def _affine_invariant(
    reference: tuple[Fraction, ...], candidate: tuple[Fraction, ...]
) -> Fraction:
    reference_m2, reference_m3 = _centered_moments(reference)
    candidate_m2, candidate_m3 = _centered_moments(candidate)
    return (
        reference_m2**3 * candidate_m3**2
        - reference_m3**2 * candidate_m2**3
    )


def test_centered_moment_invariant_vanishes_on_affine_spectra() -> None:
    reference = (Fraction(-3), Fraction(0), Fraction(1), Fraction(2))
    candidate = tuple(Fraction(7) - 5 * value for value in reversed(reference))

    assert _affine_invariant(reference, candidate) == 0
    assert _affine_invariant(reference, candidate[:-1] + (Fraction(100),)) != 0


def test_affine_spectral_gate_uses_effective_log_cap() -> None:
    result = decide_affine_spectral_obstruction(
        dual_margin_lower=Fraction(1, 10**10),
        cells=36,
        m2=Fraction(54),
        m3=Fraction(-27),
        h_hs_cap=Fraction(15, 2),
        effective_log_defect=Fraction(1, 10**6),
    )

    assert result.centered_defect_cap == Fraction(1, 500000)
    assert result.nonlinear_remainder_upper > 0
    assert result.invariant_margin == (
        result.linear_invariant_lower - result.nonlinear_remainder_upper
    )


def test_zero_affine_margin_stays_inconclusive() -> None:
    result = decide_affine_spectral_obstruction(
        dual_margin_lower=Fraction(0),
        cells=36,
        m2=Fraction(54),
        m3=Fraction(-27),
        h_hs_cap=Fraction(15, 2),
        effective_log_defect=Fraction(1, 10**6),
    )

    assert result.status == "inconclusive"
    assert result.invariant_margin < 0


def test_issue128_scale_fixture_has_positive_affine_margin() -> None:
    result = decide_affine_spectral_obstruction(
        dual_margin_lower=Fraction(6, 10**11),
        cells=36,
        m2=Fraction(54),
        m3=Fraction(-27),
        h_hs_cap=Fraction(15, 2),
        effective_log_defect=Fraction(3, 2_000_000),
    )

    assert result.status == "certified_affine_effective_spectral_obstruction"
    assert result.invariant_margin > 0


def test_affine_remainder_cap_bounds_exact_noncommuting_fixture() -> None:
    hamiltonian = sp.diag(-2, 0, 3)
    defect = sp.Matrix(
        [
            [sp.Rational(1, 100), sp.Rational(1, 200), 0],
            [sp.Rational(1, 200), sp.Rational(-1, 150), sp.Rational(1, 300)],
            [0, sp.Rational(1, 300), sp.Rational(1, 400)],
        ]
    )
    dimension = hamiltonian.rows

    def tau(matrix: sp.MatrixBase) -> sp.Expr:
        return sp.trace(matrix) / dimension

    identity = sp.eye(dimension)
    centered_h = hamiltonian - tau(hamiltonian) * identity
    centered_d = defect - tau(defect) * identity
    candidate = centered_h + centered_d
    moment2 = tau(centered_h**2)
    moment3 = tau(centered_h**3)
    candidate_m2 = tau(candidate**2)
    candidate_m3 = tau(candidate**3)
    witness = centered_h**2 - moment2 * identity - (moment3 / moment2) * centered_h
    pairing = tau(witness * defect)
    invariant = moment2**3 * candidate_m3**2 - moment3**2 * candidate_m2**3
    linear = 6 * moment3 * moment2**3 * pairing
    defect_row_sum = max(
        sum(abs(defect[row, column]) for column in range(dimension))
        for row in range(dimension)
    )
    result = decide_affine_spectral_obstruction(
        dual_margin_lower=Fraction(abs(pairing)),
        cells=1,
        m2=Fraction(moment2),
        m3=Fraction(moment3),
        h_hs_cap=Fraction(3),
        effective_log_defect=Fraction(defect_row_sum),
    )

    assert abs(Fraction(invariant - linear)) <= result.nonlinear_remainder_upper


@pytest.mark.parametrize(
    ("changes", "error"),
    (
        ({"dual_margin_lower": Fraction(-1)}, "nonnegative"),
        ({"cells": True}, "positive integer"),
        ({"m2": Fraction(0)}, "m2 must be positive"),
        ({"m3": Fraction(0)}, "m3 must be nonzero"),
        ({"h_hs_cap": Fraction(7)}, "Hilbert--Schmidt cap"),
        ({"effective_log_defect": Fraction(-1)}, "nonnegative"),
    ),
)
def test_affine_spectral_gate_rejects_invalid_inputs(
    changes: dict[str, object], error: str
) -> None:
    arguments: dict[str, object] = {
        "dual_margin_lower": Fraction(1),
        "cells": 36,
        "m2": Fraction(54),
        "m3": Fraction(-27),
        "h_hs_cap": Fraction(15, 2),
        "effective_log_defect": Fraction(1, 10**6),
    }
    arguments.update(changes)
    with pytest.raises((TypeError, ValueError), match=error):
        decide_affine_spectral_obstruction(**arguments)
