from __future__ import annotations

from fractions import Fraction

from .local_commutators import (
    CoordinateRegistry,
    SymplecticDyadicLocalDensityEvaluator,
    SymplecticPauli,
)
from .rational_words import WordPolynomial
from .refined_error import canonicalize_symplectic_unit_cell


CoordinatePauli = tuple[tuple[int, int, str], ...]
RationalLocalLieTerms = dict[CoordinatePauli, Fraction]


def _coordinate_pauli(
    registry: CoordinateRegistry,
    pauli: SymplecticPauli,
) -> CoordinatePauli:
    x_mask, z_mask = pauli
    sites = x_mask | z_mask
    result: list[tuple[int, int, str]] = []
    while sites:
        bit = sites & -sites
        site = bit.bit_length() - 1
        has_x = bool(x_mask & bit)
        has_z = bool(z_mask & bit)
        operator = "Y" if has_x and has_z else ("X" if has_x else "Z")
        x, y = registry.coordinate(site)
        result.append((x, y, operator))
        sites ^= bit
    return tuple(sorted(result))


def local_terms_l1(terms: RationalLocalLieTerms) -> Fraction:
    return sum((abs(value) for value in terms.values()), Fraction())


def project_lie_polynomial_to_cell(
    polynomial: WordPolynomial,
) -> RationalLocalLieTerms:
    """Project a homogeneous free-Lie element to the shared `2 x 2` cell.

    The input uses its associative-word representation.  The
    Dynkin--Specht--Wever projection maps each degree-`d` word to its
    right-nested commutator divided by `d`.  The local evaluator supplies the
    exact matching-fragment commutator with a common dyadic denominator.
    """

    nonzero = {
        word: Fraction(coefficient)
        for word, coefficient in polynomial.items()
        if coefficient
    }
    if not nonzero:
        return {}
    degrees = {len(word) for word in nonzero}
    if len(degrees) != 1 or 0 in degrees:
        raise ValueError("local Lie projection requires a homogeneous degree")
    degree = next(iter(degrees))
    if any(any(not 0 <= letter < 4 for letter in word) for word in nonzero):
        raise ValueError("local Lie projection requires four fragment letters")

    evaluator = SymplecticDyadicLocalDensityEvaluator(
        shared_coordinates=True
    )
    registry = evaluator.registries[0]
    denominator = degree * (
        1 << evaluator.denominator_exponent((0,) * degree)
    )
    result: RationalLocalLieTerms = {}
    for word, word_coefficient in nonzero.items():
        for raw_pauli, numerator in evaluator.evaluate(word).items():
            canonical = canonicalize_symplectic_unit_cell(
                registry,
                raw_pauli,
            )
            pauli = _coordinate_pauli(registry, canonical)
            updated = result.get(pauli, Fraction()) + (
                word_coefficient * Fraction(numerator, denominator)
            )
            if updated:
                result[pauli] = updated
            else:
                result.pop(pauli, None)
    return result
