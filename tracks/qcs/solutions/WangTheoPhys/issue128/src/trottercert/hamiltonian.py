from __future__ import annotations

from fractions import Fraction

from .algebra import PauliString, PauliSum
from .lattice import SquareLattice

SymplecticPauli = tuple[int, int]


def tfim_terms(
    length: int,
    field: Fraction,
    coupling: Fraction,
    periodic: bool,
) -> dict[SymplecticPauli, Fraction]:
    """Return canonical exact terms for ``-J XX - h Z`` on a 1D chain.

    Keys are phase-free symplectic masks ``(x_mask, z_mask)``.  Open nearest
    neighbours are always included; the last-to-first bond is included only
    when ``periodic`` is true.  Zero coefficients are omitted and insertion
    order is the canonical tuple ordering of the masks.
    """

    if isinstance(length, bool) or not isinstance(length, int) or length < 2:
        raise ValueError("TFIM length must be an integer at least two")
    if not isinstance(field, Fraction) or not isinstance(coupling, Fraction):
        raise TypeError("TFIM field and coupling must be exact Fractions")
    if not isinstance(periodic, bool):
        raise TypeError("TFIM periodic flag must be bool")

    terms: dict[SymplecticPauli, Fraction] = {}
    if coupling:
        bonds = {(site, site + 1) for site in range(length - 1)}
        if periodic:
            bonds.add((0, length - 1))
        for left, right in bonds:
            key = ((1 << left) | (1 << right), 0)
            terms[key] = terms.get(key, Fraction()) - coupling
    if field:
        for site in range(length):
            terms[(0, 1 << site)] = -field
    return {key: terms[key] for key in sorted(terms) if terms[key]}


def xxz_bond(delta: Fraction) -> dict[str, Fraction]:
    """Return exact coefficients for ``(XX + YY + delta ZZ) / 4``."""

    if not isinstance(delta, Fraction):
        raise TypeError("delta must be a Fraction")
    return {
        "XX": Fraction(1, 4),
        "YY": Fraction(1, 4),
        "ZZ": delta / 4,
    }


def heisenberg_bond(u: int, v: int) -> PauliSum:
    result = PauliSum.zero()
    for op in ("X", "Y", "Z"):
        result += PauliSum.term(PauliString({u: op, v: op}), Fraction(1, 4))
    return result


def fragment_from_bonds(bonds: tuple[tuple[int, int], ...]) -> PauliSum:
    result = PauliSum.zero()
    for u, v in bonds:
        result += heisenberg_bond(u, v)
    return result


def four_matching_fragments(lattice: SquareLattice) -> tuple[PauliSum, ...]:
    return tuple(fragment_from_bonds(group) for group in lattice.four_matchings())


def full_heisenberg_hamiltonian(lattice: SquareLattice) -> PauliSum:
    return fragment_from_bonds(lattice.bonds())
