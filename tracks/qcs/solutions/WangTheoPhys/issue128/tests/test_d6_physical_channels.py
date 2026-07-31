from __future__ import annotations

from fractions import Fraction

import pytest

from trottercert.cubic_field import Cubic
from trottercert.d6_physical_channels import (
    build_grouped_d6_bound,
    classify_coordinate_pauli,
)


def test_classifies_connected_path_and_disconnected_support() -> None:
    path = ((0, 0, "X"), (1, 0, "Y"), (2, 0, "Z"))
    disconnected = ((0, 0, "X"), (2, 0, "Z"))
    assert classify_coordinate_pauli(path).shape == "path"
    assert classify_coordinate_pauli(path).bbox == (3, 1)
    assert classify_coordinate_pauli(path).manhattan_diameter == 2
    assert classify_coordinate_pauli(disconnected).component_count == 2
    assert classify_coordinate_pauli(disconnected).shape == "disconnected"


def test_classifies_cycle_and_branched_tree() -> None:
    cycle = (
        (0, 0, "X"),
        (0, 1, "X"),
        (1, 0, "X"),
        (1, 1, "X"),
    )
    branch = (
        (0, 0, "X"),
        (-1, 0, "X"),
        (1, 0, "X"),
        (0, 1, "X"),
    )
    assert classify_coordinate_pauli(cycle).shape == "cycle"
    assert classify_coordinate_pauli(cycle).edge_count == 4
    assert classify_coordinate_pauli(branch).shape == "branched_tree"


def test_grouped_bound_pairs_local_anticommuting_channels() -> None:
    terms = {
        ((0, 0, "X"),): Cubic.one(),
        ((0, 0, "Z"),): Cubic.one(),
        ((2, 0, "X"),): Cubic.one(),
    }
    result = build_grouped_d6_bound(
        terms,
        decimal_digits=12,
        candidate_cap=8,
    )
    assert result.term_count == 3
    assert sorted(len(group) for group in result.groups) == [1, 2]
    assert result.grouped_cell_bound < result.l1_cell_bound
    assert result.grouped_site_bound == result.grouped_cell_bound / 4
    assert result.channel_counts == (("singleton", 3),)
    assert result.channel_l1_bounds == (("singleton", Fraction(3)),)


def test_grouped_bound_is_deterministic_and_covers_sorted_term_indices() -> None:
    terms = {
        ((1, 0, "Y"),): Cubic(Fraction(2, 3), 0, 0),
        ((0, 0, "X"),): Cubic.one(),
        ((0, 0, "Z"),): Cubic.one(),
    }
    first = build_grouped_d6_bound(terms, decimal_digits=12, candidate_cap=8)
    second = build_grouped_d6_bound(terms, decimal_digits=12, candidate_cap=8)
    assert first == second
    assert sorted(index for group in first.groups for index in group) == [0, 1, 2]
    assert first.max_group_size == 2


def test_grouped_bound_rejects_invalid_candidate_cap() -> None:
    with pytest.raises(ValueError, match="candidate cap"):
        build_grouped_d6_bound({}, decimal_digits=12, candidate_cap=0)
