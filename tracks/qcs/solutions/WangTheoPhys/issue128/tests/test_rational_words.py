from fractions import Fraction

from trottercert.processed_kernels import RationalStage
from trottercert.rational_words import (
    commutator_polynomial,
    lyndon_words,
    rational_formula_log_series,
    standard_lyndon_bracket,
    word_l1,
)


def test_commutator_and_lyndon_dimensions() -> None:
    left = {(0,): Fraction(1)}
    right = {(1,): Fraction(1)}
    assert commutator_polynomial(left, right) == {
        (0, 1): Fraction(1),
        (1, 0): Fraction(-1),
    }
    assert len(lyndon_words(4, 2)) == 6
    degree_four = lyndon_words(4, 4)
    assert len(degree_four) == 60
    assert all(standard_lyndon_bracket(word) for word in degree_four)


def test_exact_strang_log_has_no_degree_two() -> None:
    stages = (
        RationalStage(0, Fraction(1, 2)),
        RationalStage(1, Fraction(1)),
        RationalStage(0, Fraction(1, 2)),
    )
    logarithm = rational_formula_log_series(stages, 3)
    assert logarithm[1] == {(0,): Fraction(1), (1,): Fraction(1)}
    assert word_l1(logarithm[2]) == 0
    assert word_l1(logarithm[3]) > 0


def test_exact_log_series_drops_zero_coefficients() -> None:
    stages = (
        RationalStage(0, Fraction(1)),
        RationalStage(0, Fraction(-1)),
    )
    logarithm = rational_formula_log_series(stages, 4)
    assert all(not degree for degree in logarithm)
