from fractions import Fraction

import pytest

from trottercert.rational_lie_local import (
    local_terms_l1,
    project_lie_polynomial_to_cell,
)


def _add_terms(left, right):
    result = dict(left)
    for pauli, coefficient in right.items():
        updated = result.get(pauli, Fraction()) + coefficient
        if updated:
            result[pauli] = updated
        else:
            result.pop(pauli, None)
    return result


def test_zero_and_commutator_projection() -> None:
    assert not project_lie_polynomial_to_cell({})
    projected = project_lie_polynomial_to_cell(
        {(0, 1): Fraction(1), (1, 0): Fraction(-1)}
    )
    assert projected
    assert local_terms_l1(projected) > 0


def test_projection_is_exactly_linear() -> None:
    left = {(0, 1): Fraction(1), (1, 0): Fraction(-1)}
    right = {(2, 3): Fraction(2), (3, 2): Fraction(-2)}
    combined = {**left, **right}
    assert project_lie_polynomial_to_cell(combined) == _add_terms(
        project_lie_polynomial_to_cell(left),
        project_lie_polynomial_to_cell(right),
    )


def test_projection_rejects_mixed_degrees() -> None:
    with pytest.raises(ValueError, match="homogeneous"):
        project_lie_polynomial_to_cell(
            {(0, 1): Fraction(1), (0, 1, 2): Fraction(1)}
        )
