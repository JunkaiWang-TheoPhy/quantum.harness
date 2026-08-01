from __future__ import annotations

from fractions import Fraction

import pytest

from trottercert.cubic_field import Cubic
from trottercert.pf4_bch_mapping import pf4_suzuki_gamma
from trottercert.tfim_gauge_witness import (
    CENTERED_POLYNOMIAL_SEARCH_ORDER,
    build_tfim_centered_polynomial_witnesses,
    count_tfim_pf4_pauli_orbits,
    tfim_gauge_orbit_formula,
    verify_tfim_gauge_orbit_formula,
)


def _cubic_over_common_denominator(
    numerators: tuple[int, int, int], denominator: int
) -> Cubic:
    return Cubic(*(Fraction(value, denominator) for value in numerators))


def test_centered_polynomial_search_order_is_frozen() -> None:
    assert CENTERED_POLYNOMIAL_SEARCH_ORDER == (2, 3, 4)


@pytest.mark.parametrize(("length", "regime"), [(4, "l4_wraparound"), (6, "generic_even_l_ge_6")])
def test_parameter_free_pauli_orbit_polynomial_certificate(
    length: int, regime: str
) -> None:
    certificate = count_tfim_pf4_pauli_orbits(length)
    assert certificate.regime == regime
    assert certificate.matches_closed_formula
    assert certificate.tau_h2_coefficients == (
        ((0, 2), Fraction(length)),
        ((2, 0), Fraction(length)),
    )


def test_l4_operator_witnesses_are_centered_commuting_and_exact() -> None:
    result = build_tfim_centered_polynomial_witnesses(
        4, Fraction(1), Fraction(1)
    )
    assert tuple(record.power for record in result.witnesses) == (2, 3, 4)
    assert all(record.tau_w == 0 for record in result.witnesses)
    assert all(record.tau_wh == 0 for record in result.witnesses)
    assert all(record.commutator_is_zero for record in result.witnesses)
    assert result.witnesses[0].tau_w_e5 == Cubic.zero()
    assert result.witnesses[1].hamiltonian_coefficient == 21
    assert result.witnesses[1].tau_w_e5 == _cubic_over_common_denominator(
        (4144, 2191, 1624), 84375
    )
    assert result.witnesses[2].tau_w_e5 == Cubic.zero()
    assert result.first_nonzero_power == 3
    verify_tfim_gauge_orbit_formula(result)


@pytest.mark.parametrize(
    ("length", "b", "numerators", "denominator"),
    [
        (6, 30, (592, 313, 232), 5625),
        (8, 42, (2368, 1252, 928), 16875),
        (10, 54, (592, 313, 232), 3375),
    ],
)
def test_generic_even_operator_pairing_matches_orbit_formula(
    length: int,
    b: int,
    numerators: tuple[int, int, int],
    denominator: int,
) -> None:
    result = build_tfim_centered_polynomial_witnesses(
        length, Fraction(1), Fraction(1)
    )
    p2, p3, p4 = result.witnesses
    assert p2.tau_w_e5 == Cubic.zero()
    assert p3.hamiltonian_coefficient == b == 6 * length - 6
    assert p3.tau_w_e5 == _cubic_over_common_denominator(
        numerators, denominator
    )
    assert p3.tau_w_e5 == pf4_suzuki_gamma() * Fraction(1280, 3) * length
    assert p4.tau_w_e5 == Cubic.zero()
    assert all(record.commutator_is_zero for record in result.witnesses)
    verify_tfim_gauge_orbit_formula(result)


def test_generic_orbit_formula_is_exact_away_from_h_equals_j() -> None:
    h = Fraction(2)
    j = Fraction(3)
    result = build_tfim_centered_polynomial_witnesses(6, h, j)
    formula = tfim_gauge_orbit_formula(6, h, j)
    p3 = result.witnesses[1]
    assert p3.hamiltonian_coefficient == formula.p3_hamiltonian_coefficient
    assert p3.tau_w_e5 == pf4_suzuki_gamma() * formula.tau_w3_e5_over_gamma
    assert formula.tau_w3_e5_over_gamma > 0
    verify_tfim_gauge_orbit_formula(result)


def test_l4_wraparound_exception_is_not_promoted_to_generic_nonzero_claim() -> None:
    h = Fraction(4)
    j = Fraction(3)
    formula = tfim_gauge_orbit_formula(4, h, j)
    assert formula.regime == "l4_wraparound"
    assert formula.tau_w3_e5_over_gamma == 0
    assert formula.nonzero_condition == "h*j!=0 and 9*h^2!=16*j^2"

    result = build_tfim_centered_polynomial_witnesses(4, h, j)
    assert result.witnesses[1].tau_w_e5 == Cubic.zero()
    assert result.first_nonzero_power is None
    verify_tfim_gauge_orbit_formula(result)


@pytest.mark.parametrize(
    "arguments",
    [
        (True, Fraction(1), Fraction(1)),
        (5, Fraction(1), Fraction(1)),
        (2, Fraction(1), Fraction(1)),
        (4, 1, Fraction(1)),
        (4, Fraction(1), 1.0),
        (4, Fraction(0), Fraction(0)),
    ],
)
def test_witness_builder_rejects_out_of_scope_inputs(
    arguments: tuple[object, object, object],
) -> None:
    with pytest.raises((TypeError, ValueError)):
        build_tfim_centered_polynomial_witnesses(*arguments)  # type: ignore[arg-type]
