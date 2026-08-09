from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Sequence

import sympy as sp

from .processed_kernels import (
    ProcessedKernel,
    alternating_kernel_stages,
)
from .rational_words import (
    Word,
    WordPolynomial,
    add_polynomials,
    commutator_polynomial,
    lyndon_words,
    rational_formula_log_series,
    scale_polynomial,
    standard_lyndon_bracket,
    word_l1,
)


@dataclass(frozen=True)
class ProcessorCoordinates:
    basis_words: tuple[Word, ...]
    coordinates: tuple[Fraction, ...]
    polynomial: WordPolynomial


@dataclass(frozen=True)
class LeastSquaresLieImage:
    coordinates: tuple[Fraction, ...]
    approximation: WordPolynomial
    residual: WordPolynomial


@dataclass(frozen=True)
class EffectiveOrderAudit:
    kernel_name: str
    r2: ProcessorCoordinates
    r4: ProcessorCoordinates
    processed_degree_three: WordPolynomial
    processed_degree_five: WordPolynomial
    processed_degree_seven: WordPolynomial

    @property
    def processed_degree_three_l1(self) -> Fraction:
        return word_l1(self.processed_degree_three)

    @property
    def processed_degree_five_l1(self) -> Fraction:
        return word_l1(self.processed_degree_five)

    @property
    def processed_degree_seven_l1(self) -> Fraction:
        return word_l1(self.processed_degree_seven)


def _linear_combination(
    basis: Sequence[WordPolynomial],
    coordinates: Sequence[Fraction],
) -> WordPolynomial:
    if len(basis) != len(coordinates):
        raise ValueError("basis and coordinate lengths disagree")
    result: WordPolynomial = {}
    for polynomial, coordinate in zip(basis, coordinates):
        result = add_polynomials(
            result,
            scale_polynomial(polynomial, coordinate),
        )
    return result


def _as_sympy(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


def _as_fraction(value: sp.Expr) -> Fraction:
    rational = sp.cancel(value)
    numerator, denominator = rational.as_numer_denom()
    return Fraction(int(numerator), int(denominator))


def least_squares_lie_image(
    target: WordPolynomial,
    image_basis: Sequence[WordPolynomial],
) -> LeastSquaresLieImage:
    """Project a rational word polynomial onto a rational Lie-image span."""

    column_count = len(image_basis)
    if not column_count:
        return LeastSquaresLieImage((), {}, dict(target))

    rows: dict[Word, dict[int, Fraction]] = {}
    for column, polynomial in enumerate(image_basis):
        for word, coefficient in polynomial.items():
            rows.setdefault(word, {})[column] = coefficient

    gram = [
        [Fraction() for _ in range(column_count)]
        for _ in range(column_count)
    ]
    rhs = [Fraction() for _ in range(column_count)]
    for word, entries in rows.items():
        target_value = target.get(word, Fraction())
        ordered_entries = tuple(entries.items())
        for left, left_value in ordered_entries:
            rhs[left] += left_value * target_value
            for right, right_value in ordered_entries:
                gram[left][right] += left_value * right_value

    matrix = sp.Matrix(
        [[_as_sympy(value) for value in row] for row in gram]
    )
    vector = sp.Matrix([_as_sympy(value) for value in rhs])
    solution_set = sp.linsolve((matrix, vector))
    raw_solution = next(iter(solution_set), ())
    if not raw_solution:
        raise ArithmeticError("Lie-image normal equations have no solution")
    substitutions = {
        symbol: sp.Integer(0)
        for expression in raw_solution
        for symbol in expression.free_symbols
    }
    coordinates = tuple(
        _as_fraction(expression.subs(substitutions))
        for expression in raw_solution
    )
    approximation = _linear_combination(image_basis, coordinates)
    residual = add_polynomials(
        target,
        scale_polynomial(approximation, -1),
    )
    return LeastSquaresLieImage(
        coordinates=coordinates,
        approximation=approximation,
        residual=residual,
    )


def _sum_polynomials(*polynomials: WordPolynomial) -> WordPolynomial:
    result: WordPolynomial = {}
    for polynomial in polynomials:
        result = add_polynomials(result, polynomial)
    return result


@lru_cache(maxsize=8)
def audit_effective_order_six(
    kernel: ProcessedKernel,
) -> EffectiveOrderAudit:
    logarithm = rational_formula_log_series(
        alternating_kernel_stages(kernel),
        7,
    )
    a = logarithm[1]
    b3 = logarithm[3]
    b5 = logarithm[5]
    b7 = logarithm[7]

    r2_words = lyndon_words(4, 2)
    r2_basis = tuple(standard_lyndon_bracket(word) for word in r2_words)
    r2_images = tuple(
        commutator_polynomial(polynomial, a)
        for polynomial in r2_basis
    )
    r2_fit = least_squares_lie_image(
        scale_polynomial(b3, -1),
        r2_images,
    )
    r2_polynomial = _linear_combination(r2_basis, r2_fit.coordinates)
    r2 = ProcessorCoordinates(
        basis_words=r2_words,
        coordinates=r2_fit.coordinates,
        polynomial=r2_polynomial,
    )
    r2_a = commutator_polynomial(r2_polynomial, a)
    processed_degree_three = _sum_polynomials(b3, r2_a)

    degree_five_without_r4 = _sum_polynomials(
        b5,
        commutator_polynomial(r2_polynomial, b3),
        scale_polynomial(
            commutator_polynomial(r2_polynomial, r2_a),
            Fraction(1, 2),
        ),
    )
    r4_words = lyndon_words(4, 4)
    r4_basis = tuple(standard_lyndon_bracket(word) for word in r4_words)
    r4_images = tuple(
        commutator_polynomial(polynomial, a)
        for polynomial in r4_basis
    )
    r4_fit = least_squares_lie_image(
        scale_polynomial(degree_five_without_r4, -1),
        r4_images,
    )
    r4_polynomial = _linear_combination(r4_basis, r4_fit.coordinates)
    r4 = ProcessorCoordinates(
        basis_words=r4_words,
        coordinates=r4_fit.coordinates,
        polynomial=r4_polynomial,
    )
    r4_a = commutator_polynomial(r4_polynomial, a)
    processed_degree_five = _sum_polynomials(
        degree_five_without_r4,
        r4_a,
    )

    r2_b3 = commutator_polynomial(r2_polynomial, b3)
    processed_degree_seven = _sum_polynomials(
        b7,
        commutator_polynomial(r2_polynomial, b5),
        commutator_polynomial(r4_polynomial, b3),
        scale_polynomial(
            commutator_polynomial(r2_polynomial, r2_b3),
            Fraction(1, 2),
        ),
        scale_polynomial(
            commutator_polynomial(r2_polynomial, r4_a),
            Fraction(1, 2),
        ),
        scale_polynomial(
            commutator_polynomial(r4_polynomial, r2_a),
            Fraction(1, 2),
        ),
        scale_polynomial(
            commutator_polynomial(
                r2_polynomial,
                commutator_polynomial(r2_polynomial, r2_a),
            ),
            Fraction(1, 6),
        ),
    )
    return EffectiveOrderAudit(
        kernel_name=kernel.name,
        r2=r2,
        r4=r4,
        processed_degree_three=processed_degree_three,
        processed_degree_five=processed_degree_five,
        processed_degree_seven=processed_degree_seven,
    )
