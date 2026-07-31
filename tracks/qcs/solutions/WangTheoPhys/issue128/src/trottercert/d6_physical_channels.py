from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from heapq import heappop, heappush

from .anticommuting import (
    certify_anticommuting_partition,
    symplectic_anticommutes,
)
from .cubic_field import Cubic
from .hpc_artifacts import CoordinateCubicTerms, CoordinatePauli
from .intervals import RationalInterval, cube_root_four_interval
from .local_commutators import SymplecticPauli


Coordinate = tuple[int, int]


@dataclass(frozen=True)
class ChannelClass:
    support_size: int
    component_count: int
    edge_count: int
    bbox: tuple[int, int]
    manhattan_diameter: int
    shape: str


@dataclass(frozen=True)
class PhysicalChannelSummary:
    shape: str
    term_count: int
    l1_bound: Fraction


@dataclass(frozen=True)
class GroupedD6Bound:
    term_count: int
    groups: tuple[tuple[int, ...], ...]
    group_bounds: tuple[Fraction, ...]
    l1_cell_bound: Fraction
    grouped_cell_bound: Fraction
    grouped_site_bound: Fraction
    max_group_size: int
    channel_counts: tuple[tuple[str, int], ...]
    channel_l1_bounds: tuple[tuple[str, Fraction], ...]


def _occupied_coordinates(pauli: CoordinatePauli) -> tuple[Coordinate, ...]:
    coordinates = tuple(sorted((x, y) for x, y, _ in pauli))
    if len(coordinates) != len(set(coordinates)):
        raise ValueError("coordinate Pauli contains duplicate occupied sites")
    return coordinates


def classify_coordinate_pauli(pauli: CoordinatePauli) -> ChannelClass:
    """Classify the nearest-neighbor geometry of a coordinate Pauli string."""

    coordinates = _occupied_coordinates(pauli)
    if not coordinates:
        raise ValueError("coordinate Pauli support must be nonempty")
    occupied = set(coordinates)
    adjacency: dict[Coordinate, tuple[Coordinate, ...]] = {}
    for x, y in coordinates:
        adjacency[(x, y)] = tuple(
            sorted(
                neighbor
                for neighbor in (
                    (x - 1, y),
                    (x + 1, y),
                    (x, y - 1),
                    (x, y + 1),
                )
                if neighbor in occupied
            )
        )

    remaining = set(coordinates)
    component_count = 0
    while remaining:
        component_count += 1
        stack = [min(remaining)]
        remaining.remove(stack[0])
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)

    support_size = len(coordinates)
    edge_count = sum(len(neighbors) for neighbors in adjacency.values()) // 2
    width = max(x for x, _ in coordinates) - min(x for x, _ in coordinates) + 1
    height = max(y for _, y in coordinates) - min(y for _, y in coordinates) + 1
    manhattan_diameter = max(
        abs(left[0] - right[0]) + abs(left[1] - right[1])
        for left in coordinates
        for right in coordinates
    )
    degrees = tuple(len(adjacency[coordinate]) for coordinate in coordinates)

    if component_count > 1:
        shape = "disconnected"
    elif support_size == 1:
        shape = "singleton"
    elif support_size == 2 and edge_count == 1:
        shape = "bond"
    elif edge_count == support_size - 1 and max(degrees) <= 2:
        shape = "path"
    elif edge_count == support_size and all(degree == 2 for degree in degrees):
        shape = "cycle"
    elif edge_count == support_size - 1 and max(degrees) >= 3:
        shape = "branched_tree"
    else:
        shape = "connected_other"
    return ChannelClass(
        support_size=support_size,
        component_count=component_count,
        edge_count=edge_count,
        bbox=(width, height),
        manhattan_diameter=manhattan_diameter,
        shape=shape,
    )


def _coordinate_masks(
    paulis: tuple[CoordinatePauli, ...],
) -> tuple[SymplecticPauli, ...]:
    coordinates = sorted(
        {(x, y) for pauli in paulis for x, y, _ in pauli}
    )
    coordinate_to_bit = {
        coordinate: 1 << index for index, coordinate in enumerate(coordinates)
    }
    masks: list[SymplecticPauli] = []
    for pauli in paulis:
        _occupied_coordinates(pauli)
        x_mask = z_mask = 0
        for x, y, axis in pauli:
            bit = coordinate_to_bit[(x, y)]
            if axis in {"X", "Y"}:
                x_mask |= bit
            if axis in {"Z", "Y"}:
                z_mask |= bit
        masks.append((x_mask, z_mask))
    if len(masks) != len(set(masks)):
        raise ValueError("coordinate Paulis do not map injectively to masks")
    return tuple(masks)


def _bounded_local_candidates(
    support: frozenset[Coordinate],
    inverted: Mapping[Coordinate, list[int]],
    unpaired: set[int],
    cap: int,
) -> tuple[int, ...]:
    """Read a bounded sorted union from lazily cleaned posting heaps."""

    frontier: list[tuple[int, Coordinate]] = []
    removed: dict[Coordinate, list[int]] = {}

    def expose_next(coordinate: Coordinate) -> None:
        posting = inverted[coordinate]
        while posting and posting[0] not in unpaired:
            heappop(posting)
        if posting:
            candidate = heappop(posting)
            removed.setdefault(coordinate, []).append(candidate)
            heappush(frontier, (candidate, coordinate))

    for coordinate in sorted(support):
        expose_next(coordinate)

    result: list[int] = []
    seen: set[int] = set()
    while frontier and len(result) < cap:
        candidate, coordinate = heappop(frontier)
        if candidate not in seen:
            seen.add(candidate)
            result.append(candidate)
        expose_next(coordinate)

    for coordinate, candidates in removed.items():
        posting = inverted[coordinate]
        for candidate in candidates:
            heappush(posting, candidate)
    return tuple(result)


def build_grouped_d6_bound(
    terms: CoordinateCubicTerms,
    decimal_digits: int,
    candidate_cap: int = 128,
) -> GroupedD6Bound:
    """Discover local anticommuting pairs and certify the complete partition.

    Discovery is deterministic and deliberately bounded per term.  The final
    result is recomputed by the existing exact partition certifier, so neither
    the candidate ordering nor the cap belongs to the trusted base.
    """

    if candidate_cap < 1:
        raise ValueError("candidate cap must be positive")
    if decimal_digits < 1:
        raise ValueError("decimal digits must be positive")

    paulis = tuple(sorted(terms))
    masks = _coordinate_masks(paulis)
    root = cube_root_four_interval(decimal_digits)
    intervals = tuple(terms[pauli].enclose(root) for pauli in paulis)
    coefficients: dict[SymplecticPauli, RationalInterval] = {
        mask: interval for mask, interval in zip(masks, intervals)
    }
    l1_cell_bound = sum(
        (interval.abs_upper() for interval in intervals),
        Fraction(),
    )

    inverted: dict[Coordinate, list[int]] = {}
    for index, pauli in enumerate(paulis):
        for coordinate in _occupied_coordinates(pauli):
            inverted.setdefault(coordinate, []).append(index)

    order = tuple(
        sorted(
            range(len(paulis)),
            key=lambda index: (-intervals[index].abs_upper(), paulis[index]),
        )
    )
    unpaired = set(range(len(paulis)))
    groups: list[tuple[int, ...]] = []
    for index in order:
        if index not in unpaired:
            continue
        unpaired.remove(index)
        support = frozenset(_occupied_coordinates(paulis[index]))
        candidates = _bounded_local_candidates(
            support,
            inverted,
            unpaired,
            candidate_cap,
        )
        partner = next(
            (
                candidate
                for candidate in candidates
                if symplectic_anticommutes(masks[index], masks[candidate])
            ),
            None,
        )
        if partner is None:
            groups.append((index,))
        else:
            unpaired.remove(partner)
            groups.append((index, partner))

    mask_groups = tuple(
        tuple(masks[index] for index in group) for group in groups
    )
    certificate = certify_anticommuting_partition(coefficients, mask_groups)
    if certificate.bound > l1_cell_bound:
        raise ArithmeticError("grouped D6 bound exceeds its Pauli-l1 baseline")

    channel_count_map: Counter[str] = Counter()
    channel_l1_map: dict[str, Fraction] = {}
    for pauli, interval in zip(paulis, intervals):
        shape = classify_coordinate_pauli(pauli).shape
        channel_count_map[shape] += 1
        channel_l1_map[shape] = (
            channel_l1_map.get(shape, Fraction()) + interval.abs_upper()
        )

    return GroupedD6Bound(
        term_count=len(paulis),
        groups=tuple(groups),
        group_bounds=tuple(group.bound for group in certificate.groups),
        l1_cell_bound=l1_cell_bound,
        grouped_cell_bound=certificate.bound,
        grouped_site_bound=certificate.bound / 4,
        max_group_size=max((len(group) for group in groups), default=0),
        channel_counts=tuple(sorted(channel_count_map.items())),
        channel_l1_bounds=tuple(sorted(channel_l1_map.items())),
    )
