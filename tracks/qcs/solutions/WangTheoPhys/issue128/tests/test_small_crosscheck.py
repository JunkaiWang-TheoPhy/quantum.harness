from fractions import Fraction

from trottercert.crosscheck import (
    open_rectangular_matchings,
    small_exact_crosscheck,
    small_open_exact_crosscheck,
)


def test_small_exact_evolution_is_below_certificate() -> None:
    result = small_exact_crosscheck(2, Fraction(1, 1000))
    assert result["bound_dominates_empirical_error"]
    assert result["bound_meets_tolerance"]
    assert result["steps"] == 97


def test_non_degenerate_open_evolution_uses_direct_finite_certificate() -> None:
    groups = open_rectangular_matchings(2, 3)
    bonds = tuple(bond for group in groups for bond in group)
    assert len(bonds) == len(set(bonds)) == 7

    result = small_open_exact_crosscheck(2, 3, Fraction(1, 10**6))
    assert result["empirical_error_is_diagnostic_only"] is True
    assert result["steps"] == 95
    assert result["unique_bond_count"] == 7
    assert result["nonempty_fragment_count"] == 3
    assert result["certified_upper_method"] == (
        "direct_finite_open_published_triangle_exact_pauli_l1"
    )
    assert result["bound_dominates_empirical_error"]
    assert result["bound_meets_tolerance"]
