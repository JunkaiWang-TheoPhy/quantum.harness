"""Exact cyclic free-trace identity for the five-copy Suzuki PF4 formula."""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction

from .cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from .cubic_local import cubic_formula_log_series

Word = tuple[int, ...]
CubicWordPolynomial = dict[Word, Cubic]


def pf4_suzuki_gamma() -> Cubic:
    """Return the exact positive prefactor in the PF4 trace identity."""

    return Cubic(
        Fraction(37, 900000),
        Fraction(313, 14400000),
        Fraction(29, 1800000),
    )


def _add(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    result = dict(left)
    for word, coefficient in right.items():
        updated = result.get(word, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            result.pop(word, None)
        else:
            result[word] = updated
    return result


def _scale(
    polynomial: Mapping[Word, Cubic],
    scalar: Cubic | Fraction | int,
) -> CubicWordPolynomial:
    factor = Cubic.coerce(scalar)
    return {
        word: scaled
        for word, coefficient in polynomial.items()
        if (scaled := coefficient * factor) != Cubic.zero()
    }


def _multiply(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    result: CubicWordPolynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            updated = result.get(word, Cubic.zero()) + (
                left_coefficient * right_coefficient
            )
            if updated == Cubic.zero():
                result.pop(word, None)
            else:
                result[word] = updated
    return result


def _commutator(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    return _add(_multiply(left, right), _scale(_multiply(right, left), -1))


def cyclic_trace_classes(
    polynomial: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    """Collect a word polynomial modulo exact cyclic trace equivalence."""

    classes: CubicWordPolynomial = {}
    for word, coefficient in polynomial.items():
        representative = (
            min(word[offset:] + word[:offset] for offset in range(len(word)))
            if word
            else ()
        )
        updated = classes.get(representative, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            classes.pop(representative, None)
        else:
            classes[representative] = updated
    return classes


def pf4_suzuki_trace_polynomials() -> tuple[
    CubicWordPolynomial,
    CubicWordPolynomial,
]:
    """Return the raw free-word polynomials on the two sides of the identity."""

    a = {(0,): Cubic.one()}
    b = {(1,): Cubic.one()}
    hamiltonian = _add(a, b)
    logarithm = cubic_formula_log_series(
        fourth_order_suzuki_cubic_stages(2),
        5,
    )
    if logarithm[1] != hamiltonian or any(
        logarithm[degree] for degree in (2, 3, 4)
    ):
        raise ArithmeticError("exact Suzuki stages are not fourth order")

    c = _commutator(a, _commutator(a, b))
    d = _commutator(b, _commutator(b, a))
    quadratic = _add(
        _multiply(c, c),
        _add(
            _scale(_multiply(c, d), -4),
            _scale(_multiply(d, d), Fraction(8, 3)),
        ),
    )
    return (
        _multiply(hamiltonian, logarithm[5]),
        _scale(quadratic, pf4_suzuki_gamma()),
    )


def verify_pf4_suzuki_trace_identity() -> None:
    """Raise unless all ten exact cyclic word classes agree."""

    left, right = pf4_suzuki_trace_polynomials()
    left_classes = cyclic_trace_classes(left)
    right_classes = cyclic_trace_classes(right)
    if len(left_classes) != 10 or left_classes != right_classes:
        raise ArithmeticError("exact PF4 cyclic free-trace identity failed")
