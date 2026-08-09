#!/usr/bin/env python3
"""Experimental, non-production 4x4 full-phase D4 feasibility audit.

This script does not modify a certificate.  It checks which canonical D4
density terms admit a genuine open 4x4 cluster placement while retaining
their (2,2)-cell phase, builds the corresponding interval-valued patch
operator, and reports a rigorous hybrid upper bound.  Terms that do not fit
remain in the frozen sidecar's anticommuting partition.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from trottercert.anticommuting import (
    certify_anticommuting_partition,
    discover_anticommuting_partition,
    sqrt_fraction_upper,
)
from trottercert.intervals import RationalInterval
from trottercert.local_commutators import SymplecticDyadicLocalDensityEvaluator
from trottercert.refined_error import interval_formula_log_series
from trottercert.rigorous_fourth import fourth_order_suzuki_interval_stages


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SIDECAR = ROOT / "certificates" / "issue128-d4-groups.json"
DEFAULT_OUTPUT = ROOT / "docs" / "experiments" / "full-phase-4x4" / "result.json"


def pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def reconstruct_registry() -> SymplecticDyadicLocalDensityEvaluator:
    """Reproduce the deterministic site registry used by the D4 sidecar."""

    stages, _ = fourth_order_suzuki_interval_stages(4, decimal_digits=12)
    logarithm = interval_formula_log_series(stages, 5)
    evaluator = SymplecticDyadicLocalDensityEvaluator(shared_coordinates=True)
    for word in logarithm[5]:
        evaluator.evaluate(word)
        evaluator.cache.pop(word, None)
    return evaluator


def coordinates_for_pauli(registry, pauli: tuple[int, int]):
    x_mask, z_mask = pauli
    sites = x_mask | z_mask
    result = []
    while sites:
        bit = sites & -sites
        site = bit.bit_length() - 1
        has_x = bool(x_mask & bit)
        has_z = bool(z_mask & bit)
        op = "Y" if has_x and has_z else ("X" if has_x else "Z")
        x, y = registry.coordinate(site)
        result.append((x, y, op))
        sites ^= bit
    return tuple(result)


def embed(coords, dx: int, dy: int) -> tuple[int, int]:
    x_mask = z_mask = 0
    for x, y, op in coords:
        bit = 1 << ((y + dy) * 4 + x + dx)
        if op in {"X", "Y"}:
            x_mask |= bit
        if op in {"Z", "Y"}:
            z_mask |= bit
    return x_mask, z_mask


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    payload = json.loads(args.sidecar.read_text())
    coefficient_denominator = int(payload["coefficient_denominator"])
    sqrt_denominator = int(payload["sqrt_denominator"])
    evaluator = reconstruct_registry()
    registry = evaluator.registries[0]

    coefficients: list[RationalInterval] = []
    coordinates = []
    fits = []
    bbox_counts: dict[str, int] = {}
    for x_mask, z_mask, lower, upper in payload["terms"]:
        coords = coordinates_for_pauli(registry, (int(x_mask), int(z_mask)))
        interval = RationalInterval(
            Fraction(int(lower), coefficient_denominator),
            Fraction(int(upper), coefficient_denominator),
        )
        min_x = min(x for x, _, _ in coords)
        min_y = min(y for _, y, _ in coords)
        max_x = max(x for x, _, _ in coords)
        max_y = max(y for _, y, _ in coords)
        width = max_x - min_x + 1
        height = max_y - min_y + 1
        key = f"{width}x{height}"
        bbox_counts[key] = bbox_counts.get(key, 0) + 1
        # Retaining a colored-cell phase means coordinates may not be shifted
        # by one site merely to make them fit.  Genuine patch placements use
        # nonnegative translations in strides of two.
        fit = min_x >= 0 and min_y >= 0 and max_x < 4 and max_y < 4
        coefficients.append(interval)
        coordinates.append(coords)
        fits.append(fit)

    patch_coefficients: dict[tuple[int, int], RationalInterval] = {}
    for coords, interval, fit in zip(coordinates, coefficients, fits):
        if not fit:
            continue
        max_x = max(x for x, _, _ in coords)
        max_y = max(y for _, y, _ in coords)
        offsets = tuple(
            (dx, dy)
            for dy in range(0, 4 - max_y, 2)
            for dx in range(0, 4 - max_x, 2)
        )
        weight = Fraction(1, len(offsets))
        for dx, dy in offsets:
            pauli = embed(coords, dx, dy)
            patch_coefficients[pauli] = patch_coefficients.get(
                pauli, RationalInterval.point(0)
            ) + interval * weight

    patch_groups = discover_anticommuting_partition(
        patch_coefficients, max_group_size=10
    )
    patch_certificate = certify_anticommuting_partition(
        patch_coefficients, patch_groups
    )

    current_cell_bound = Fraction(*payload["cell_bound"])
    outside_bound = Fraction()
    all_fit_groups = mixed_groups = outside_groups = 0
    for indices, _ in payload["groups"]:
        outside_squared = sum(
            (
                coefficients[int(index)].abs_upper() ** 2
                for index in indices
                if not fits[int(index)]
            ),
            Fraction(),
        )
        outside_bound += sqrt_fraction_upper(outside_squared)
        flags = {fits[int(index)] for index in indices}
        if flags == {True}:
            all_fit_groups += 1
        elif flags == {False}:
            outside_groups += 1
        else:
            mixed_groups += 1

    hybrid_cell_bound = outside_bound + patch_certificate.bound
    total_l1 = sum((value.abs_upper() for value in coefficients), Fraction())
    fit_l1 = sum(
        (value.abs_upper() for value, fit in zip(coefficients, fits) if fit),
        Fraction(),
    )
    patch_l1 = sum(
        (value.abs_upper() for value in patch_coefficients.values()), Fraction()
    )

    r = 78
    d4_contribution = Fraction(144) * current_cell_bound / (4 * 5 * r**4)
    hybrid_d4_contribution = Fraction(144) * hybrid_cell_bound / (4 * 5 * r**4)
    d4_alone_target_cell = Fraction(1, 10**6) * (4 * 5 * r**4) / 144
    # Frozen D5-integrated ledger, supplied by the independent budget audit.
    non_d4_r78 = (
        Fraction(934095, 10**13)
        + Fraction(584064, 10**12)
        + Fraction(871506, 10**13)
        + Fraction(900764, 10**12)
    )
    non_d4_d6_r78 = (
        Fraction(934095, 10**13)
        + Fraction(871506, 10**13)
        + Fraction(900764, 10**12)
    )

    result = {
        "status": "NO-GO",
        "scope": "experimental_nonproduction_full_phase_4x4_D4_audit",
        "reason": (
            "A phase-preserving open 4x4 patch covers only a small subset of "
            "the D4 density, and splitting that subset from the frozen "
            "partition worsens the rigorous upper bound. Independently, the "
            "frozen non-D4 r=78 ledger already exceeds the total tolerance."
        ),
        "sidecar": str(args.sidecar.relative_to(ROOT)),
        "term_count": len(coefficients),
        "phase_preserving_4x4_fit_terms": sum(fits),
        "outside_terms": len(fits) - sum(fits),
        "fit_fraction": pair(Fraction(sum(fits), len(fits))),
        "bbox_counts": dict(sorted(bbox_counts.items())),
        "current_cell_bound": pair(current_cell_bound),
        "current_site_bound": pair(current_cell_bound / 4),
        "current_r78_d4_contribution": pair(d4_contribution),
        "d4_alone_r78_target_cell_bound": pair(d4_alone_target_cell),
        "fit_term_l1": pair(fit_l1),
        "total_term_l1": pair(total_l1),
        "fit_l1_fraction": pair(fit_l1 / total_l1),
        "patch_merged_term_count": len(patch_coefficients),
        "patch_l1_bound": pair(patch_l1),
        "patch_anticommuting_group_count": len(patch_certificate.groups),
        "patch_anticommuting_bound": pair(patch_certificate.bound),
        "frozen_partition_outside_bound": pair(outside_bound),
        "hybrid_cell_bound": pair(hybrid_cell_bound),
        "hybrid_r78_d4_contribution": pair(hybrid_d4_contribution),
        "hybrid_relative_to_current": pair(hybrid_cell_bound / current_cell_bound),
        "partition_group_classes": {
            "all_fit": all_fit_groups,
            "mixed": mixed_groups,
            "all_outside": outside_groups,
        },
        "non_d4_r78_contribution_from_frozen_ledger": pair(non_d4_r78),
        "non_d4_d6_r78_contribution_from_frozen_ledger": pair(non_d4_d6_r78),
        "tolerance": [1, 10**6],
        "rigor_notes": [
            "All coefficient and bound arithmetic is Fraction-based.",
            "Every reported patch group is rechecked for exact pairwise anticommutation.",
            "The patch result is an upper bound, not a spectral-norm estimate.",
            "Periodic 4x4 aliasing is excluded because it is not an open-cluster decomposition of the 12x12 operator.",
            "This experiment is not a production certificate and does not alter the frozen result.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    def decimal(value: Fraction) -> str:
        return f"{float(value):.12g}"

    print("status=NO-GO")
    print(f"fit_terms={sum(fits)}/{len(fits)}")
    print(f"fit_l1_fraction={decimal(fit_l1 / total_l1)}")
    print(f"current_cell_bound={decimal(current_cell_bound)}")
    print(f"patch_bound={decimal(patch_certificate.bound)}")
    print(f"outside_bound={decimal(outside_bound)}")
    print(f"hybrid_cell_bound={decimal(hybrid_cell_bound)}")
    print(f"hybrid/current={decimal(hybrid_cell_bound / current_cell_bound)}")
    print(f"r78_non_d4={decimal(non_d4_r78)}")
    print(f"r78_non_d4_d6={decimal(non_d4_d6_r78)}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
