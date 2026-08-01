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


def xxz_bond_operator(left: int, right: int, delta: Fraction) -> PauliSum:
    """Return the exact operator ``(XX + YY + delta*ZZ) / 4`` on one bond."""

    for name, site in (("left", left), ("right", right)):
        if isinstance(site, bool) or not isinstance(site, int):
            raise TypeError(f"{name} site must be an integer")
        if site < 0:
            raise ValueError(f"{name} site must be nonnegative")
    if left == right:
        raise ValueError("XXZ self bonds are not allowed")
    if not isinstance(delta, Fraction):
        raise TypeError("delta must be an exact Fraction")

    result = PauliSum.zero()
    for operator, coefficient in (
        ("X", Fraction(1, 4)),
        ("Y", Fraction(1, 4)),
        ("Z", delta / 4),
    ):
        result += PauliSum.term(
            PauliString({left: operator, right: operator}), coefficient
        )
    return result


def _xxz_fragment_from_bonds(
    bonds: tuple[tuple[int, int], ...],
    delta: Fraction,
) -> PauliSum:
    result = PauliSum.zero()
    for left, right in bonds:
        result += xxz_bond_operator(left, right, delta)
    return result


def _validated_four_matchings(
    lattice: SquareLattice,
) -> tuple[tuple[tuple[int, int], ...], ...]:
    if not isinstance(lattice, SquareLattice):
        raise TypeError("lattice must be a SquareLattice")
    matchings = lattice.four_matchings()
    expected_per_matching = lattice.n_sites // 2
    if len(matchings) != 4:
        raise ValueError("square-lattice XXZ decomposition requires four matchings")
    if any(len(group) != expected_per_matching for group in matchings):
        raise ValueError("each XXZ matching must cover every site exactly once")
    if any(left == right for group in matchings for left, right in group):
        raise ValueError("XXZ matching contains a self bond")
    if any(len(set(group)) != len(group) for group in matchings):
        raise ValueError("XXZ matching contains a duplicate bond")
    if any(
        not set(matchings[left]).isdisjoint(matchings[right])
        for left in range(4)
        for right in range(left + 1, 4)
    ):
        raise ValueError("XXZ matching bond sets must be pairwise disjoint")
    flattened_sites = [
        [site for bond in group for site in bond] for group in matchings
    ]
    if any(len(sites) != len(set(sites)) for sites in flattened_sites):
        raise ValueError("a matching may not use one site twice")
    if set().union(*(set(group) for group in matchings)) != set(lattice.bonds()):
        raise ValueError("four XXZ matchings must cover the periodic bond set exactly")
    return matchings


def four_matching_xxz_fragments(
    lattice: SquareLattice,
    delta: Fraction,
) -> tuple[PauliSum, PauliSum, PauliSum, PauliSum]:
    """Return exact horizontal/vertical even/odd XXZ matching fragments."""

    if not isinstance(delta, Fraction):
        raise TypeError("delta must be an exact Fraction")
    matchings = _validated_four_matchings(lattice)
    fragments = tuple(_xxz_fragment_from_bonds(group, delta) for group in matchings)
    return fragments[0], fragments[1], fragments[2], fragments[3]


def finite_torus_xxz_hamiltonian(
    lattice: SquareLattice,
    delta: Fraction,
) -> PauliSum:
    """Independently enumerate the exact periodic finite-torus XXZ Hamiltonian."""

    if not isinstance(lattice, SquareLattice):
        raise TypeError("lattice must be a SquareLattice")
    if not isinstance(delta, Fraction):
        raise TypeError("delta must be an exact Fraction")
    return _xxz_fragment_from_bonds(lattice.bonds(), delta)


def fragment_from_bonds(bonds: tuple[tuple[int, int], ...]) -> PauliSum:
    result = PauliSum.zero()
    for u, v in bonds:
        result += heisenberg_bond(u, v)
    return result


def four_matching_fragments(lattice: SquareLattice) -> tuple[PauliSum, ...]:
    return tuple(fragment_from_bonds(group) for group in lattice.four_matchings())


def full_heisenberg_hamiltonian(lattice: SquareLattice) -> PauliSum:
    return fragment_from_bonds(lattice.bonds())
