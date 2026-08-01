"""Strict spatial-orbit primitives for the periodic 4 by 4 grouped XXZ model.

The eight frozen transformations are the D4 point symmetries about lattice
site ``(0, 0)``.  Each record binds both its permutation of the 16 finite
sites and the induced permutation of the four matching fragments.  Word
orbits are classification metadata only: this module deliberately has no API
that combines theorem weights across distinct raw words or norm blocks.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction


SITE_COUNT = 16
FRAGMENT_COUNT = 4
WORD_LENGTH = 5
SymplecticPauli = tuple[int, int]
GaussianRational = tuple[Fraction, Fraction]
FiveLetterWord = tuple[int, int, int, int, int]


@dataclass(frozen=True, slots=True)
class XXZSpatialSymmetry:
    """One exact source-to-target D4 action on sites and fragments."""

    name: str
    site_permutation: tuple[int, ...]
    fragment_permutation: tuple[int, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise TypeError("spatial symmetry name must be a nonempty string")
        _validate_permutation(
            self.site_permutation, SITE_COUNT, "site permutation"
        )
        _validate_permutation(
            self.fragment_permutation,
            FRAGMENT_COUNT,
            "fragment permutation",
        )


@dataclass(frozen=True, slots=True)
class FiveLetterWordOrbit:
    """A D4 word orbit with no coefficient or theorem-weight aggregation."""

    representative: FiveLetterWord
    members: tuple[FiveLetterWord, ...]
    stabilizer_size: int
    raw_word_count: int


@dataclass(frozen=True, slots=True)
class CommutatorOrbitInvariantReport:
    """Exact invariants established for one transported commutator map."""

    source_word: FiveLetterWord
    target_word: FiveLetterWord
    symmetry_name: str
    term_count: int
    source_axis_l1: Fraction
    target_axis_l1: Fraction
    source_hilbert_schmidt_squared: Fraction
    target_hilbert_schmidt_squared: Fraction
    exact_map_transport: bool
    axis_l1_invariant: bool
    hilbert_schmidt_squared_invariant: bool
    commutation_graph_invariant: bool


def _validate_permutation(
    permutation: object,
    size: int,
    label: str,
) -> tuple[int, ...]:
    if not isinstance(permutation, tuple) or len(permutation) != size:
        raise TypeError(f"{label} must be a tuple of length {size}")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in permutation):
        raise TypeError(f"{label} entries must be integers")
    if tuple(sorted(permutation)) != tuple(range(size)):
        raise ValueError(f"{label} must contain each target exactly once")
    return permutation


XXZ_4X4_SPATIAL_SYMMETRIES = (
    XXZSpatialSymmetry(
        "identity",
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
        (0, 1, 2, 3),
    ),
    XXZSpatialSymmetry(
        "rotate_90",
        (0, 4, 8, 12, 3, 7, 11, 15, 2, 6, 10, 14, 1, 5, 9, 13),
        (2, 3, 1, 0),
    ),
    XXZSpatialSymmetry(
        "rotate_180",
        (0, 3, 2, 1, 12, 15, 14, 13, 8, 11, 10, 9, 4, 7, 6, 5),
        (1, 0, 3, 2),
    ),
    XXZSpatialSymmetry(
        "rotate_270",
        (0, 12, 8, 4, 1, 13, 9, 5, 2, 14, 10, 6, 3, 15, 11, 7),
        (3, 2, 0, 1),
    ),
    XXZSpatialSymmetry(
        "reflect_x_axis",
        (0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7),
        (0, 1, 3, 2),
    ),
    XXZSpatialSymmetry(
        "reflect_y_axis",
        (0, 3, 2, 1, 4, 7, 6, 5, 8, 11, 10, 9, 12, 15, 14, 13),
        (1, 0, 2, 3),
    ),
    XXZSpatialSymmetry(
        "reflect_main_diagonal",
        (0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15),
        (2, 3, 0, 1),
    ),
    XXZSpatialSymmetry(
        "reflect_anti_diagonal",
        (0, 12, 8, 4, 3, 15, 11, 7, 2, 14, 10, 6, 1, 13, 9, 5),
        (3, 2, 1, 0),
    ),
)


def _compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def _site(x: int, y: int) -> int:
    return (y % 4) * 4 + (x % 4)


def _canonical_bond(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _frozen_matchings() -> tuple[frozenset[tuple[int, int]], ...]:
    groups: list[set[tuple[int, int]]] = [set(), set(), set(), set()]
    for y in range(4):
        for x in range(0, 4, 2):
            groups[0].add(_canonical_bond(_site(x, y), _site(x + 1, y)))
            groups[1].add(_canonical_bond(_site(x + 1, y), _site(x + 2, y)))
    for x in range(4):
        for y in range(0, 4, 2):
            groups[2].add(_canonical_bond(_site(x, y), _site(x, y + 1)))
            groups[3].add(_canonical_bond(_site(x, y + 1), _site(x, y + 2)))
    return tuple(frozenset(group) for group in groups)


def verify_frozen_spatial_symmetries() -> None:
    """Verify D4 closure and the site-induced matching action."""

    if len(XXZ_4X4_SPATIAL_SYMMETRIES) != 8:
        raise ArithmeticError("the frozen spatial action must contain eight elements")
    if len({symmetry.name for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES}) != 8:
        raise ArithmeticError("spatial symmetry names must be unique")
    by_sites = {
        symmetry.site_permutation: symmetry
        for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
    }
    if len(by_sites) != 8:
        raise ArithmeticError("spatial site permutations must be unique")
    matchings = _frozen_matchings()
    for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES:
        for source, target in enumerate(symmetry.fragment_permutation):
            transported = frozenset(
                _canonical_bond(
                    symmetry.site_permutation[left],
                    symmetry.site_permutation[right],
                )
                for left, right in matchings[source]
            )
            if transported != matchings[target]:
                raise ArithmeticError(
                    f"{symmetry.name} fragment permutation is not site-induced"
                )
    for left in XXZ_4X4_SPATIAL_SYMMETRIES:
        for right in XXZ_4X4_SPATIAL_SYMMETRIES:
            composed_sites = _compose(
                left.site_permutation, right.site_permutation
            )
            composed_fragments = _compose(
                left.fragment_permutation, right.fragment_permutation
            )
            composed = by_sites.get(composed_sites)
            if composed is None or composed.fragment_permutation != composed_fragments:
                raise ArithmeticError("site and fragment D4 actions are not a homomorphism")


def _validate_symmetry(symmetry: object) -> XXZSpatialSymmetry:
    if not isinstance(symmetry, XXZSpatialSymmetry):
        raise TypeError("symmetry must be an XXZSpatialSymmetry")
    if symmetry not in XXZ_4X4_SPATIAL_SYMMETRIES:
        raise ValueError("symmetry is not one of the frozen 4 by 4 records")
    return symmetry


def _validate_word(word: object) -> FiveLetterWord:
    if not isinstance(word, tuple) or len(word) != WORD_LENGTH:
        raise TypeError("five-letter word must be a tuple of length five")
    if any(isinstance(index, bool) or not isinstance(index, int) for index in word):
        raise TypeError("five-letter word entries must be integer fragment indices")
    if any(not 0 <= index < FRAGMENT_COUNT for index in word):
        raise ValueError("five-letter word contains an invalid fragment index")
    return word  # type: ignore[return-value]


def transport_five_letter_word(
    word: FiveLetterWord,
    symmetry: XXZSpatialSymmetry,
) -> FiveLetterWord:
    """Apply one frozen fragment permutation without changing multiplicity."""

    checked_word = _validate_word(word)
    checked_symmetry = _validate_symmetry(symmetry)
    return tuple(
        checked_symmetry.fragment_permutation[index] for index in checked_word
    )  # type: ignore[return-value]


def five_letter_word_orbit(word: FiveLetterWord) -> FiveLetterWordOrbit:
    """Return the unique D4 orbit; no weights or coefficient maps are merged."""

    checked = _validate_word(word)
    actions = tuple(
        transport_five_letter_word(checked, symmetry)
        for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
    )
    members = tuple(sorted(set(actions)))
    stabilizer = sum(transformed == checked for transformed in actions)
    if len(members) * stabilizer != 8:
        raise ArithmeticError("five-letter word orbit-stabilizer identity failed")
    return FiveLetterWordOrbit(
        representative=members[0],
        members=members,
        stabilizer_size=stabilizer,
        raw_word_count=len(members),
    )


def canonical_five_letter_word(word: FiveLetterWord) -> FiveLetterWord:
    return five_letter_word_orbit(word).representative


def _validate_pauli_mask(mask: object) -> SymplecticPauli:
    if not isinstance(mask, tuple) or len(mask) != 2:
        raise TypeError("Pauli mask must be an (x_mask,z_mask) tuple")
    x_mask, z_mask = mask
    if any(isinstance(value, bool) or not isinstance(value, int) for value in mask):
        raise TypeError("Pauli masks must be integers")
    if x_mask < 0 or z_mask < 0 or x_mask >= 1 << 16 or z_mask >= 1 << 16:
        raise ValueError("Pauli mask must be supported on exactly the 16-site universe")
    return x_mask, z_mask


def _transport_bit_mask(mask: int, permutation: tuple[int, ...]) -> int:
    result = 0
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        source = bit.bit_length() - 1
        result |= 1 << permutation[source]
        remaining ^= bit
    return result


def transport_pauli_mask(
    mask: SymplecticPauli,
    symmetry: XXZSpatialSymmetry,
) -> SymplecticPauli:
    checked_mask = _validate_pauli_mask(mask)
    checked_symmetry = _validate_symmetry(symmetry)
    return (
        _transport_bit_mask(checked_mask[0], checked_symmetry.site_permutation),
        _transport_bit_mask(checked_mask[1], checked_symmetry.site_permutation),
    )


def pauli_masks_anticommute(
    left: SymplecticPauli,
    right: SymplecticPauli,
) -> bool:
    lhs = _validate_pauli_mask(left)
    rhs = _validate_pauli_mask(right)
    return bool(
        (((lhs[0] & rhs[1]).bit_count() + (lhs[1] & rhs[0]).bit_count()) & 1)
    )


def _validate_coefficient(value: object) -> GaussianRational:
    if not isinstance(value, tuple) or len(value) != 2:
        raise TypeError("commutator coefficient must be an exact (real,imag) tuple")
    if any(not isinstance(component, Fraction) for component in value):
        raise TypeError("commutator coefficient components must be Fractions")
    if value == (Fraction(), Fraction()):
        raise ValueError("zero commutator coefficients must be omitted")
    return value


def _canonical_commutator_map(
    coefficients: object,
) -> dict[SymplecticPauli, GaussianRational]:
    if not isinstance(coefficients, Mapping):
        raise TypeError("commutator map must be a mapping")
    result: dict[SymplecticPauli, GaussianRational] = {}
    for mask, coefficient in coefficients.items():
        checked_mask = _validate_pauli_mask(mask)
        checked_coefficient = _validate_coefficient(coefficient)
        result[checked_mask] = checked_coefficient
    return {mask: result[mask] for mask in sorted(result)}


def transport_commutator_map(
    coefficients: Mapping[SymplecticPauli, GaussianRational],
    symmetry: XXZSpatialSymmetry,
) -> dict[SymplecticPauli, GaussianRational]:
    """Transport one unweighted word map without combining different words."""

    source = _canonical_commutator_map(coefficients)
    checked_symmetry = _validate_symmetry(symmetry)
    result: dict[SymplecticPauli, GaussianRational] = {}
    for mask, coefficient in source.items():
        transported = transport_pauli_mask(mask, checked_symmetry)
        if transported in result:
            raise ArithmeticError("site permutation unexpectedly collided Pauli masks")
        result[transported] = coefficient
    return {mask: result[mask] for mask in sorted(result)}


def _axis_l1(coefficients: Mapping[SymplecticPauli, GaussianRational]) -> Fraction:
    return sum(
        (abs(real) + abs(imag) for real, imag in coefficients.values()),
        Fraction(),
    )


def _hilbert_schmidt_squared(
    coefficients: Mapping[SymplecticPauli, GaussianRational],
) -> Fraction:
    return sum(
        (real**2 + imag**2 for real, imag in coefficients.values()),
        Fraction(),
    )


def commutation_graph(
    masks: Sequence[SymplecticPauli],
) -> set[tuple[SymplecticPauli, SymplecticPauli]]:
    """Return exact anticommutation edges of a finite Pauli mask sequence."""

    if not isinstance(masks, Sequence):
        raise TypeError("commutation graph masks must be a sequence")
    checked = tuple(_validate_pauli_mask(mask) for mask in masks)
    if len(set(checked)) != len(checked):
        raise ValueError("commutation graph masks must be unique")
    return {
        tuple(sorted((left, right)))  # type: ignore[misc]
        for index, left in enumerate(checked)
        for right in checked[index + 1 :]
        if pauli_masks_anticommute(left, right)
    }


def verify_commutation_graph_invariant(
    masks: Sequence[SymplecticPauli],
    symmetry: XXZSpatialSymmetry,
) -> bool:
    """Explicitly verify graph isomorphism for a supplied finite mask set."""

    checked_symmetry = _validate_symmetry(symmetry)
    checked_masks = tuple(_validate_pauli_mask(mask) for mask in masks)
    source_graph = commutation_graph(checked_masks)
    transported_masks = tuple(
        transport_pauli_mask(mask, checked_symmetry) for mask in checked_masks
    )
    target_graph = commutation_graph(transported_masks)
    expected = {
        tuple(
            sorted(
                (
                    transport_pauli_mask(left, checked_symmetry),
                    transport_pauli_mask(right, checked_symmetry),
                )
            )
        )  # type: ignore[misc]
        for left, right in source_graph
    }
    if target_graph != expected:
        raise ArithmeticError("Pauli commutation graph is not symmetry-invariant")
    return True


def verify_commutator_orbit_invariants(
    source_word: FiveLetterWord,
    source_coefficients: Mapping[SymplecticPauli, GaussianRational],
    target_word: FiveLetterWord,
    target_coefficients: Mapping[SymplecticPauli, GaussianRational],
    symmetry: XXZSpatialSymmetry,
) -> CommutatorOrbitInvariantReport:
    """Verify exact map transport and its coefficient/graph invariants.

    Exact equality with the transported map proves commutation-graph
    invariance for arbitrarily large maps because a common site permutation
    preserves the binary symplectic form.  The separate explicit graph API is
    available for small audit fixtures without imposing quadratic cost on a
    production theorem block.
    """

    checked_symmetry = _validate_symmetry(symmetry)
    checked_source_word = _validate_word(source_word)
    checked_target_word = _validate_word(target_word)
    expected_word = transport_five_letter_word(
        checked_source_word, checked_symmetry
    )
    if checked_target_word != expected_word:
        raise ValueError("target word is not the frozen spatial transport")
    source = _canonical_commutator_map(source_coefficients)
    target = _canonical_commutator_map(target_coefficients)
    expected_map = transport_commutator_map(source, checked_symmetry)
    if target != expected_map:
        raise ValueError("target commutator map is not the exact spatial transport")
    source_l1 = _axis_l1(source)
    target_l1 = _axis_l1(target)
    source_hs = _hilbert_schmidt_squared(source)
    target_hs = _hilbert_schmidt_squared(target)
    l1_invariant = source_l1 == target_l1
    hs_invariant = source_hs == target_hs
    if not l1_invariant or not hs_invariant:
        raise ArithmeticError("commutator coefficient norm changed under transport")
    # A bijective common permutation of x/z bit positions preserves
    # x.left*z.right + z.left*x.right modulo two for every mask pair.
    graph_invariant = tuple(sorted(checked_symmetry.site_permutation)) == tuple(
        range(SITE_COUNT)
    )
    if not graph_invariant:
        raise ArithmeticError("site action does not preserve the commutation graph")
    return CommutatorOrbitInvariantReport(
        source_word=checked_source_word,
        target_word=checked_target_word,
        symmetry_name=checked_symmetry.name,
        term_count=len(source),
        source_axis_l1=source_l1,
        target_axis_l1=target_l1,
        source_hilbert_schmidt_squared=source_hs,
        target_hilbert_schmidt_squared=target_hs,
        exact_map_transport=True,
        axis_l1_invariant=l1_invariant,
        hilbert_schmidt_squared_invariant=hs_invariant,
        commutation_graph_invariant=graph_invariant,
    )


verify_frozen_spatial_symmetries()
