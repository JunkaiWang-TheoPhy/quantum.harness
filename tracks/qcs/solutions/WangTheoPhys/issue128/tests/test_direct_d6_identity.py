from __future__ import annotations

from fractions import Fraction

from trottercert.cubic_field import Cubic
from trottercert.cubic_local import _add_cubic_operator


def test_direct_d6_identity_scalars_are_exact() -> None:
    e7 = {(1, 0): Cubic(Fraction(2, 5), 0, 0)}
    second_adjoint = {(0, 1): Cubic(Fraction(3, 7), 0, 0)}
    result = {}
    _add_cubic_operator(result, e7, 7)
    _add_cubic_operator(result, second_adjoint, Fraction(2, 3))
    assert result[(1, 0)] == Cubic(Fraction(14, 5), 0, 0)
    assert result[(0, 1)] == Cubic(Fraction(2, 7), 0, 0)
