#!/usr/bin/env python3
"""Exact leading-defect obstruction audit for endpoint processors.

This experiment is deliberately disconnected from certificate generation.  It
reconstructs the exact cubic-field E5 density of the frozen fourth-order
Suzuki formula and emits algebraic invariants that a processor commutator
cannot change.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import exact_log_e5_density
from trottercert.intervals import cube_root_four_interval


def _fraction(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _cubic(value: Cubic) -> list[list[int]]:
    return [_fraction(value.a0), _fraction(value.a1), _fraction(value.a2)]


def _add(table: dict[int, Cubic], key: int, value: Cubic) -> None:
    table[key] = table.get(key, Cubic.zero()) + value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    registry, e5 = exact_log_e5_density(fourth_order_suzuki_cubic_stages(4))
    counts: Counter[int] = Counter()
    squared_by_weight: dict[int, Cubic] = {}
    nearest_neighbor_sum = Cubic.zero()
    nearest_neighbor_terms = 0

    for (x_mask, z_mask), coefficient in e5.items():
        support = x_mask | z_mask
        weight = support.bit_count()
        counts[weight] += 1
        _add(squared_by_weight, weight, coefficient * coefficient)
        if weight != 2:
            continue
        sites = [
            bit
            for bit in range(support.bit_length())
            if support & (1 << bit)
        ]
        first = registry.coordinate(sites[0])
        second = registry.coordinate(sites[1])
        if abs(first[0] - second[0]) + abs(first[1] - second[1]) == 1:
            nearest_neighbor_sum += coefficient
            nearest_neighbor_terms += 1

    total_squared = sum(squared_by_weight.values(), Cubic.zero())
    # Every nearest-neighbor Heisenberg Pauli has coefficient 1/4.  The
    # canonical 2x2 density contains eight bonds and three axes per bond.
    e5_hs_overlap_h = nearest_neighbor_sum / 4
    root = cube_root_four_interval(40)
    squared_intervals = {
        str(weight): {
            "exact_cubic": _cubic(value),
            "interval": [
                str(value.enclose(root).lower),
                str(value.enclose(root).upper),
            ],
        }
        for weight, value in sorted(squared_by_weight.items())
    }
    weight_six_fraction = (
        squared_by_weight[6].enclose(root) / total_squared.enclose(root)
    )

    payload = {
        "schema_version": 1,
        "formula": "five-copy fourth-order Suzuki, four matchings",
        "coefficient_field": "Q(alpha), alpha^3=4",
        "e5_term_count": len(e5),
        "support_counts": {str(k): v for k, v in sorted(counts.items())},
        "squared_pauli_coefficient_norm_by_weight": squared_intervals,
        "squared_pauli_coefficient_norm_total": {
            "exact_cubic": _cubic(total_squared),
            "interval": [
                str(total_squared.enclose(root).lower),
                str(total_squared.enclose(root).upper),
            ],
        },
        "weight_six_squared_fraction_interval": [
            str(weight_six_fraction.lower),
            str(weight_six_fraction.upper),
        ],
        "nearest_neighbor_pauli_terms": nearest_neighbor_terms,
        "e5_hilbert_schmidt_overlap_with_h_per_cell": {
            "exact_cubic": _cubic(e5_hs_overlap_h),
            "interval": [
                str(e5_hs_overlap_h.enclose(root).lower),
                str(e5_hs_overlap_h.enclose(root).upper),
            ],
            "nonzero": e5_hs_overlap_h != Cubic.zero(),
        },
        "exact_consequences": {
            "all_endpoint_processors": (
                "E5 is not in image(ad_H), because its Hilbert-Schmidt "
                "overlap with H is nonzero while <H,[Q,H]>=0 for every Q"
            ),
            "support_at_most_four": (
                "[Q,H] has Pauli support at most five, but E5 has 65856 "
                "nonzero support-six terms"
            ),
        },
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)


if __name__ == "__main__":
    main()
