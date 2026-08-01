from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction

import numpy as np
import pytest

from trottercert.cubic_field import Cubic
from trottercert.pf4_bch_mapping import (
    CubicWordPolynomial,
    Word,
    cyclic_trace_classes,
    pf4_suzuki_gamma,
    pf4_suzuki_trace_polynomials,
    verify_pf4_suzuki_trace_identity,
)


def _add(
    left: Mapping[Word, Fraction], right: Mapping[Word, Fraction]
) -> dict[Word, Fraction]:
    result = dict(left)
    for word, coefficient in right.items():
        updated = result.get(word, Fraction()) + coefficient
        if updated:
            result[word] = updated
        else:
            result.pop(word, None)
    return result


def _scale(
    polynomial: Mapping[Word, Fraction], scalar: Fraction | int
) -> dict[Word, Fraction]:
    factor = Fraction(scalar)
    return {
        word: coefficient * factor
        for word, coefficient in polynomial.items()
        if coefficient * factor
    }


def _multiply(
    left: Mapping[Word, Fraction], right: Mapping[Word, Fraction]
) -> dict[Word, Fraction]:
    result: dict[Word, Fraction] = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            result = _add(
                result,
                {word: left_coefficient * right_coefficient},
            )
    return result


def _commutator(
    left: Mapping[Word, Fraction], right: Mapping[Word, Fraction]
) -> dict[Word, Fraction]:
    return _add(_multiply(left, right), _scale(_multiply(right, left), -1))


def _mutated_quadratic(
    c2: Fraction, cd: Fraction, d2: Fraction
) -> CubicWordPolynomial:
    a = {(0,): Fraction(1)}
    b = {(1,): Fraction(1)}
    c = _commutator(a, _commutator(a, b))
    d = _commutator(b, _commutator(b, a))
    rational = _add(
        _scale(_multiply(c, c), c2),
        _add(
            _scale(_multiply(c, d), cd),
            _scale(_multiply(d, d), d2),
        ),
    )
    gamma = pf4_suzuki_gamma()
    return {word: gamma * coefficient for word, coefficient in rational.items()}


def _evaluate_trace(
    polynomial: Mapping[Word, Cubic],
    a: np.ndarray,
    b: np.ndarray,
) -> Cubic:
    dimension = a.shape[0]
    identity = np.eye(dimension, dtype=object)
    generators = (a, b)
    result = Cubic.zero()
    for word, coefficient in polynomial.items():
        term = identity
        for letter in word:
            term = term @ generators[letter]
        trace = sum((Fraction(term[i, i]) for i in range(dimension)), Fraction())
        result += coefficient * trace
    return result


def test_pf4_suzuki_trace_identity_matches_frozen_cyclic_classes() -> None:
    left, right = pf4_suzuki_trace_polynomials()
    left_classes = cyclic_trace_classes(left)
    right_classes = cyclic_trace_classes(right)
    gamma = pf4_suzuki_gamma()
    expected_multiples = {
        (0, 0, 0, 0, 1, 1): Fraction(2),
        (0, 0, 0, 1, 0, 1): Fraction(-8),
        (0, 0, 0, 1, 1, 1): Fraction(-8),
        (0, 0, 1, 0, 0, 1): Fraction(6),
        (0, 0, 1, 0, 1, 1): Fraction(12),
        (0, 0, 1, 1, 0, 1): Fraction(12),
        (0, 0, 1, 1, 1, 1): Fraction(16, 3),
        (0, 1, 0, 1, 0, 1): Fraction(-16),
        (0, 1, 0, 1, 1, 1): Fraction(-64, 3),
        (0, 1, 1, 0, 1, 1): Fraction(16),
    }

    assert pf4_suzuki_gamma() == Cubic(
        Fraction(37, 900000),
        Fraction(313, 14400000),
        Fraction(29, 1800000),
    )
    assert left_classes == {
        word: gamma * multiple for word, multiple in expected_multiples.items()
    }
    assert right_classes == left_classes
    verify_pf4_suzuki_trace_identity()


@pytest.mark.parametrize(
    ("c2", "cd", "d2"),
    (
        (Fraction(2), Fraction(-4), Fraction(8, 3)),
        (Fraction(1), Fraction(-3), Fraction(8, 3)),
        (Fraction(1), Fraction(-4), Fraction(7, 3)),
    ),
)
def test_pf4_trace_identity_rejects_mutated_quadratic_coefficients(
    c2: Fraction,
    cd: Fraction,
    d2: Fraction,
) -> None:
    left, _ = pf4_suzuki_trace_polynomials()

    assert cyclic_trace_classes(_mutated_quadratic(c2, cd, d2)) != (
        cyclic_trace_classes(left)
    )


@pytest.mark.parametrize("dimension", (2, 3))
def test_pf4_trace_identity_on_exact_symmetric_matrices(dimension: int) -> None:
    rng = np.random.default_rng(12840 + dimension)
    raw_a = rng.integers(-3, 4, size=(dimension, dimension))
    raw_b = rng.integers(-3, 4, size=(dimension, dimension))
    a = np.asarray(raw_a + raw_a.T, dtype=object)
    b = np.asarray(raw_b + raw_b.T, dtype=object)
    left, right = pf4_suzuki_trace_polynomials()

    assert _evaluate_trace(left, a, b) == _evaluate_trace(right, a, b)
