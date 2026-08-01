from __future__ import annotations

from fractions import Fraction
from itertools import product

import pytest

from trottercert.grouped_xxz_orbits import (
    XXZ_4X4_SPATIAL_SYMMETRIES,
    canonical_five_letter_word,
    commutation_graph,
    five_letter_word_orbit,
    pauli_masks_anticommute,
    transport_commutator_map,
    transport_five_letter_word,
    transport_pauli_mask,
    verify_commutation_graph_invariant,
    verify_commutator_orbit_invariants,
    verify_frozen_spatial_symmetries,
)
from trottercert.lattice import SquareLattice


EXPECTED_FRAGMENT_PERMUTATIONS = (
    (0, 1, 2, 3),
    (2, 3, 1, 0),
    (1, 0, 3, 2),
    (3, 2, 0, 1),
    (0, 1, 3, 2),
    (1, 0, 2, 3),
    (2, 3, 0, 1),
    (3, 2, 1, 0),
)


def _compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Return left-after-right for source-to-target permutations."""

    return tuple(left[right[index]] for index in range(len(left)))


def test_eight_spatial_symmetries_have_frozen_permutation_tables() -> None:
    assert tuple(symmetry.name for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES) == (
        "identity",
        "rotate_90",
        "rotate_180",
        "rotate_270",
        "reflect_x_axis",
        "reflect_y_axis",
        "reflect_main_diagonal",
        "reflect_anti_diagonal",
    )
    assert tuple(
        symmetry.fragment_permutation for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
    ) == EXPECTED_FRAGMENT_PERMUTATIONS
    assert all(
        tuple(sorted(symmetry.site_permutation)) == tuple(range(16))
        for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
    )
    assert len(
        {symmetry.site_permutation for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES}
    ) == 8
    verify_frozen_spatial_symmetries()


def test_site_permutations_induce_the_frozen_fragment_permutations() -> None:
    lattice = SquareLattice(4)
    matchings = tuple(set(group) for group in lattice.four_matchings())
    for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES:
        for source_index, target_index in enumerate(symmetry.fragment_permutation):
            transported = {
                lattice.canonical_bond(
                    symmetry.site_permutation[left],
                    symmetry.site_permutation[right],
                )
                for left, right in matchings[source_index]
            }
            assert transported == matchings[target_index]


def test_spatial_tables_are_closed_under_composition() -> None:
    sites = {symmetry.site_permutation for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES}
    fragments = {
        symmetry.fragment_permutation for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
    }
    for left in XXZ_4X4_SPATIAL_SYMMETRIES:
        for right in XXZ_4X4_SPATIAL_SYMMETRIES:
            assert _compose(left.site_permutation, right.site_permutation) in sites
            assert (
                _compose(left.fragment_permutation, right.fragment_permutation)
                in fragments
            )


def test_all_five_letter_words_have_canonical_orbit_representatives() -> None:
    for word in product(range(4), repeat=5):
        orbit = five_letter_word_orbit(word)
        representative = canonical_five_letter_word(word)
        assert representative == orbit.representative == min(orbit.members)
        assert word in orbit.members
        assert len(orbit.members) * orbit.stabilizer_size == 8
        assert orbit.raw_word_count == len(orbit.members)
        for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES:
            transformed = transport_five_letter_word(word, symmetry)
            assert transformed in orbit.members
            assert canonical_five_letter_word(transformed) == representative


def test_word_orbit_keeps_distinct_words_and_has_no_weight_aggregation() -> None:
    word = (0, 0, 1, 2, 3)
    orbit = five_letter_word_orbit(word)
    assert len(orbit.members) > 1
    assert len(set(orbit.members)) == len(orbit.members)
    assert not hasattr(orbit, "weight")
    assert not hasattr(orbit, "combined_weight")


@pytest.mark.parametrize(
    "word",
    [(), (0, 1, 2, 3), (0, 1, 2, 3, 4), (0, 1, 2, 3, True)],
)
def test_five_letter_word_validation_is_strict(word: tuple[object, ...]) -> None:
    with pytest.raises((TypeError, ValueError), match="word"):
        canonical_five_letter_word(word)  # type: ignore[arg-type]


def test_pauli_mask_transport_round_trips_and_preserves_symplectic_form() -> None:
    masks = (
        (1 << 0, 0),
        (0, 1 << 0),
        ((1 << 1) | (1 << 6), 1 << 6),
        ((1 << 3) | (1 << 12), (1 << 7) | (1 << 12)),
    )
    for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES:
        inverse = next(
            candidate
            for candidate in XXZ_4X4_SPATIAL_SYMMETRIES
            if _compose(candidate.site_permutation, symmetry.site_permutation)
            == tuple(range(16))
        )
        transported = tuple(transport_pauli_mask(mask, symmetry) for mask in masks)
        assert tuple(
            transport_pauli_mask(mask, inverse) for mask in transported
        ) == masks
        for left, right in product(range(len(masks)), repeat=2):
            assert pauli_masks_anticommute(masks[left], masks[right]) == (
                pauli_masks_anticommute(transported[left], transported[right])
            )


def test_commutator_map_norm_and_commutation_graph_are_invariant() -> None:
    symmetry = XXZ_4X4_SPATIAL_SYMMETRIES[1]
    word = (0, 2, 1, 3, 0)
    coefficients = {
        (1 << 0, 0): (Fraction(1, 2), Fraction()),
        (0, 1 << 0): (Fraction(-3, 4), Fraction()),
        ((1 << 1) | (1 << 5), 1 << 5): (
            Fraction(),
            Fraction(5, 7),
        ),
        ((1 << 2) | (1 << 6), (1 << 2) | (1 << 6)): (
            Fraction(2, 3),
            Fraction(),
        ),
    }
    transported = transport_commutator_map(coefficients, symmetry)
    report = verify_commutator_orbit_invariants(
        word,
        coefficients,
        transport_five_letter_word(word, symmetry),
        transported,
        symmetry,
    )
    assert report.term_count == 4
    assert report.exact_map_transport
    assert report.axis_l1_invariant
    assert report.hilbert_schmidt_squared_invariant
    assert report.commutation_graph_invariant
    assert report.source_axis_l1 == sum(
        abs(real) + abs(imag) for real, imag in coefficients.values()
    )
    assert report.source_hilbert_schmidt_squared == sum(
        real**2 + imag**2 for real, imag in coefficients.values()
    )
    assert verify_commutation_graph_invariant(tuple(coefficients), symmetry)
    assert commutation_graph(tuple(transported)) == {
        tuple(sorted((transport_pauli_mask(left, symmetry), transport_pauli_mask(right, symmetry))))
        for left, right in commutation_graph(tuple(coefficients))
    }


@pytest.mark.parametrize("mutation", ["word", "coefficient", "mask", "missing"])
def test_commutator_orbit_verifier_rejects_mutations(mutation: str) -> None:
    symmetry = XXZ_4X4_SPATIAL_SYMMETRIES[6]
    word = (0, 1, 2, 3, 0)
    source = {
        (1 << 0, 0): (Fraction(1, 2), Fraction()),
        (0, 1 << 0): (Fraction(-3, 4), Fraction()),
    }
    target_word = transport_five_letter_word(word, symmetry)
    target = transport_commutator_map(source, symmetry)
    if mutation == "word":
        target_word = (target_word[0],) * 5
    elif mutation == "coefficient":
        key = next(iter(target))
        target[key] = (target[key][0] + 1, target[key][1])
    elif mutation == "mask":
        key = next(iter(target))
        value = target.pop(key)
        target[(key[0] ^ (1 << 15), key[1])] = value
    elif mutation == "missing":
        target.pop(next(iter(target)))
    else:
        raise AssertionError(mutation)
    with pytest.raises(ValueError, match="word|map"):
        verify_commutator_orbit_invariants(
            word, source, target_word, target, symmetry
        )
