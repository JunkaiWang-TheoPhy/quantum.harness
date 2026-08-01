"""Exact Pauli counting for the periodic transverse-field Ising family.

Only the Pauli trace moments and the resulting *conditional* PF4 algebraic
quantities are certified here.  Turning those quantities into an operator
obstruction still requires the independently missing PF4 BCH mapping.
"""

from __future__ import annotations

from fractions import Fraction

from .algebra import PauliString, PauliSum, QComplex, commutator
from .trace_obstruction import pf4_trace_quadratic_form

MOMENT_KEYS = (
    "trace_c2_over_d",
    "trace_d2_over_d",
    "trace_cd_over_d",
    "trace_h2_over_d",
    "obstruction_over_d",
    "operator_lower_bound",
)


def _validate_family_inputs(
    length: int,
    h: Fraction,
    j: Fraction,
    periodic: bool,
) -> None:
    if isinstance(length, bool) or not isinstance(length, int):
        raise TypeError("length must be an integer, not bool")
    if not isinstance(periodic, bool):
        raise TypeError("periodic must be bool")
    if not isinstance(h, Fraction) or not isinstance(j, Fraction):
        raise TypeError("h and j must be Fractions")
    if periodic and length % 2:
        raise ValueError("the certified periodic TFIM family requires even length")
    minimum = 4 if periodic else 2
    if length < minimum:
        raise ValueError(f"length must be at least {minimum}")


def _tfim_operators(
    length: int,
    h: Fraction,
    j: Fraction,
    periodic: bool,
) -> tuple[PauliSum, PauliSum]:
    field = PauliSum.zero()
    interaction = PauliSum.zero()
    for site in range(length):
        field += PauliSum.term(PauliString({site: "X"}), h)
    bond_count = length if periodic else length - 1
    for site in range(bond_count):
        interaction += PauliSum.term(
            PauliString({site: "Z", (site + 1) % length: "Z"}),
            j,
        )
    return field, interaction


def _normalized_hilbert_schmidt_square(operator: PauliSum) -> Fraction:
    result = Fraction()
    for coefficient in operator.terms.values():
        product = coefficient.conjugate() * coefficient
        if product.imag:
            raise ArithmeticError("Pauli Hilbert--Schmidt square is not real")
        result += product.real
    return result


def _normalized_trace_pairing(left: PauliSum, right: PauliSum) -> Fraction:
    result = QComplex()
    for pauli, coefficient in left.terms.items():
        result += coefficient.conjugate() * right.terms.get(pauli, QComplex())
    if result.imag:
        raise ArithmeticError("Pauli trace pairing is not real")
    return result.real


def _closed_periodic_trace_moments(
    length: int,
    h: Fraction,
    j: Fraction,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Return formulas derived from local Pauli multiplicities, not a fit."""

    return (
        128 * length * h**4 * j**2,
        128 * length * h**2 * j**4,
        Fraction(),
        length * (h**2 + j**2),
    )


def tfim_trace_moments(
    length: int,
    h: Fraction,
    j: Fraction,
    periodic: bool,
) -> dict[str, Fraction]:
    """Count exact normalized TFIM trace moments by Pauli orthogonality.

    With ``A = h sum_i X_i`` and ``B = j sum_i Z_i Z_{i+1}``, this uses
    ``C = [A,[A,B]]`` and ``D = [B,[B,A]]``.  The returned
    ``obstruction_over_d`` is only the plan-stated rational quadratic core

    ``Tr(C^2)/(2d) + 14 Tr(CD)/(3d) + 4 Tr(D^2)/(3d)``.

    ``operator_lower_bound`` divides its absolute value by the Pauli l1 upper
    bound on ``||A+B||`` when that bound is nonzero, and is defined to be zero
    at ``h=j=0``.  Both quantities are conditional algebraic values, not
    certified PF4 operator obstructions until the BCH mapping is proved.
    """

    _validate_family_inputs(length, h, j, periodic)
    field, interaction = _tfim_operators(length, h, j, periodic)
    hamiltonian = field + interaction
    c_operator = commutator(field, commutator(field, interaction))
    d_operator = commutator(interaction, commutator(interaction, field))

    c2 = _normalized_hilbert_schmidt_square(c_operator)
    d2 = _normalized_hilbert_schmidt_square(d_operator)
    cd = _normalized_trace_pairing(c_operator, d_operator)
    h2 = _normalized_hilbert_schmidt_square(hamiltonian)
    if periodic:
        expected = _closed_periodic_trace_moments(length, h, j)
        if (c2, d2, cd, h2) != expected:
            raise ArithmeticError("exact Pauli count disagrees with closed TFIM formula")

    quadratic_core = pf4_trace_quadratic_form(c2, cd, d2, Fraction(1))
    hamiltonian_norm_upper = hamiltonian.exact_real_l1()
    if hamiltonian_norm_upper:
        conditional_lower_bound = abs(quadratic_core) / hamiltonian_norm_upper
    else:
        if quadratic_core:
            raise ArithmeticError("zero Hamiltonian norm with nonzero trace form")
        conditional_lower_bound = Fraction()
    return {
        "trace_c2_over_d": c2,
        "trace_d2_over_d": d2,
        "trace_cd_over_d": cd,
        "trace_h2_over_d": h2,
        "obstruction_over_d": quadratic_core,
        "operator_lower_bound": conditional_lower_bound,
    }
