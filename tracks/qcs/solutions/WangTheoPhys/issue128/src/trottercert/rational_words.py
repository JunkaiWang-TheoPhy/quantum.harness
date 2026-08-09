from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import factorial
from typing import Sequence

from .processed_kernels import RationalStage


Word = tuple[int, ...]
WordPolynomial = dict[Word, Fraction]
WordSeries = list[WordPolynomial]


def add_polynomials(
    left: WordPolynomial,
    right: WordPolynomial,
) -> WordPolynomial:
    result = dict(left)
    for word, coefficient in right.items():
        updated = result.get(word, Fraction()) + coefficient
        if updated:
            result[word] = updated
        else:
            result.pop(word, None)
    return result


def scale_polynomial(
    polynomial: WordPolynomial,
    scalar: Fraction | int,
) -> WordPolynomial:
    factor = Fraction(scalar)
    if not factor:
        return {}
    return {
        word: coefficient * factor
        for word, coefficient in polynomial.items()
        if coefficient * factor
    }


def multiply_polynomials(
    left: WordPolynomial,
    right: WordPolynomial,
) -> WordPolynomial:
    result: WordPolynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            updated = result.get(word, Fraction()) + (
                left_coefficient * right_coefficient
            )
            if updated:
                result[word] = updated
            else:
                result.pop(word, None)
    return result


def commutator_polynomial(
    left: WordPolynomial,
    right: WordPolynomial,
) -> WordPolynomial:
    return add_polynomials(
        multiply_polynomials(left, right),
        scale_polynomial(multiply_polynomials(right, left), -1),
    )


def word_l1(polynomial: WordPolynomial) -> Fraction:
    return sum((abs(value) for value in polynomial.values()), Fraction())


def _empty_series(order: int) -> WordSeries:
    return [{} for _ in range(order + 1)]


def _identity_series(order: int) -> WordSeries:
    result = _empty_series(order)
    result[0][()] = Fraction(1)
    return result


def _add_series(left: WordSeries, right: WordSeries) -> WordSeries:
    if len(left) != len(right):
        raise ValueError("word-series truncation orders disagree")
    return [
        add_polynomials(left_degree, right_degree)
        for left_degree, right_degree in zip(left, right)
    ]


def _scale_series(series: WordSeries, scalar: Fraction | int) -> WordSeries:
    return [scale_polynomial(degree, scalar) for degree in series]


def _multiply_series(left: WordSeries, right: WordSeries) -> WordSeries:
    if len(left) != len(right):
        raise ValueError("word-series truncation orders disagree")
    order = len(left) - 1
    result = _empty_series(order)
    for total_degree in range(order + 1):
        degree_result: WordPolynomial = {}
        for left_degree in range(total_degree + 1):
            product_degree = multiply_polynomials(
                left[left_degree],
                right[total_degree - left_degree],
            )
            degree_result = add_polynomials(degree_result, product_degree)
        result[total_degree] = degree_result
    return result


def rational_formula_product_series(
    stages: Sequence[RationalStage],
    order: int,
) -> WordSeries:
    if order < 0:
        raise ValueError("series order must be nonnegative")
    result = _identity_series(order)
    for stage in stages:
        updated = _empty_series(order)
        for source_degree, source in enumerate(result):
            for power in range(order - source_degree + 1):
                scalar = stage.coefficient**power / factorial(power)
                suffix = (stage.fragment_index,) * power
                target = updated[source_degree + power]
                for word, coefficient in source.items():
                    expanded_word = word + suffix
                    value = target.get(expanded_word, Fraction()) + (
                        coefficient * scalar
                    )
                    if value:
                        target[expanded_word] = value
                    else:
                        target.pop(expanded_word, None)
        result = updated
    return result


def rational_formula_log_series(
    stages: Sequence[RationalStage],
    order: int,
) -> WordSeries:
    product_series = rational_formula_product_series(stages, order)
    delta = [dict(degree) for degree in product_series]
    identity_value = delta[0].get((), Fraction()) - 1
    if identity_value:
        delta[0][()] = identity_value
    else:
        delta[0].pop((), None)

    logarithm = _empty_series(order)
    power = _identity_series(order)
    for exponent in range(1, order + 1):
        power = _multiply_series(power, delta)
        signed_inverse = Fraction(1 if exponent % 2 else -1, exponent)
        logarithm = _add_series(
            logarithm,
            _scale_series(power, signed_inverse),
        )
    return logarithm


def _is_lyndon(word: Word) -> bool:
    if not word:
        return False
    return all(
        word < word[shift:] + word[:shift]
        for shift in range(1, len(word))
    )


@lru_cache(maxsize=None)
def lyndon_words(alphabet_size: int, degree: int) -> tuple[Word, ...]:
    if alphabet_size < 1:
        raise ValueError("alphabet size must be positive")
    if degree < 1:
        raise ValueError("Lyndon degree must be positive")
    return tuple(
        word
        for word in product(range(alphabet_size), repeat=degree)
        if _is_lyndon(word)
    )


@lru_cache(maxsize=None)
def standard_lyndon_bracket(word: Word) -> WordPolynomial:
    if not _is_lyndon(word):
        raise ValueError("standard bracketing requires a Lyndon word")
    if len(word) == 1:
        return {word: Fraction(1)}
    for split in range(1, len(word)):
        suffix = word[split:]
        if not _is_lyndon(suffix):
            continue
        prefix = word[:split]
        if not _is_lyndon(prefix):
            raise ArithmeticError("invalid standard Lyndon factorization")
        return commutator_polynomial(
            standard_lyndon_bracket(prefix),
            standard_lyndon_bracket(suffix),
        )
    raise ArithmeticError("Lyndon word has no standard factorization")
