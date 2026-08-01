#!/usr/bin/env python3
"""Generate the frozen nondegenerate dense-error calibration grid."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from decimal import ROUND_CEILING, Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath
import numpy as np
import scipy
from scipy.linalg import expm

from trottercert.crosscheck import (
    fourth_order_product_step,
    open_rectangle_dense_fragments,
    open_rectangle_published_triangle_certificate,
)

if __package__:
    from scripts.reference_verify import verify as reference_verify
else:
    from reference_verify import verify as reference_verify

PAPER_ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = PAPER_ROOT / "certificates/issue128-d5-integrated-certificate.json"
FROZEN_SIZES = ((2, 3), (2, 4), (3, 3))
FROZEN_STEPS = (8, 12, 16, 24, 32, 48, 64, 96)
NORMALIZATION = "(XX+YY+ZZ)/4"
SCIENTIFIC_DEPENDENCY_PATHS = (
    "scripts/run_ed_calibration.py",
    "scripts/reference_verify.py",
    "src/trottercert/algebra.py",
    "src/trottercert/crosscheck.py",
    "src/trottercert/hamiltonian.py",
    "src/trottercert/higher_order.py",
    "src/trottercert/intervals.py",
    "src/trottercert/lattice.py",
    "src/trottercert/local_commutators.py",
    "src/trottercert/rigorous_fourth.py",
)
BOUND_SIGNIFICANT_DIGITS = 16
ACTUAL_DISPLAY_SIGNIFICANT_DIGITS = 9


def _diagnostic_decimal(value: float) -> str:
    return format(value, f".{ACTUAL_DISPLAY_SIGNIFICANT_DIGITS - 1}e")


def _outward_decimal(value: Fraction) -> str:
    """Serialize a nonnegative exact rational toward positive infinity."""

    if value < 0:
        raise ValueError("upper-bound serialization requires a nonnegative value")
    with localcontext() as context:
        context.prec = BOUND_SIGNIFICANT_DIGITS
        context.rounding = ROUND_CEILING
        rounded = Decimal(value.numerator) / Decimal(value.denominator)
        return format(rounded, f".{BOUND_SIGNIFICANT_DIGITS - 1}E").lower()


def _fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _scientific_dependency_hashes() -> dict[str, str]:
    return {
        relative: hashlib.sha256((PAPER_ROOT / relative).read_bytes()).hexdigest()
        for relative in SCIENTIFIC_DEPENDENCY_PATHS
    }


def _canonical_digest(payload: object) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


@lru_cache(maxsize=1)
def _reference_certificate_record() -> dict[str, object]:
    verification = reference_verify(CERTIFICATE_PATH)
    if verification.get("valid") is not True:
        raise ArithmeticError("reference certificate did not verify")
    return {
        "path": str(CERTIFICATE_PATH.relative_to(PAPER_ROOT)),
        "sha256": hashlib.sha256(CERTIFICATE_PATH.read_bytes()).hexdigest(),
        "role": (
            "verified periodic reference only; not transferred to open boundaries"
        ),
        "reference_verifier": "scripts/reference_verify.py",
        "verification": verification,
    }


def _runtime_environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
    }


def _open_certificate_record(width: int, height: int) -> dict[str, object]:
    certificate = open_rectangle_published_triangle_certificate(width, height)
    unsigned: dict[str, object] = {
        "model": "open_rectangular_spin_half_isotropic_heisenberg",
        "normalization": NORMALIZATION,
        "width": width,
        "height": height,
        "boundary": "open",
        "theorem": "published_high_order_commutator_triangle_bound",
        "norm_method": "direct_finite_open_exact_Pauli_axis_l1",
        "theorem_center": certificate.center,
        "coefficient_interval_decimal_digits": (
            certificate.coefficient_interval_decimal_digits
        ),
        "root_interval": {
            "lower": _fraction_pair(certificate.root_interval.lower),
            "upper": _fraction_pair(certificate.root_interval.upper),
        },
        "constant_upper_exact": _fraction_pair(certificate.constant_upper),
        "constant_upper": _outward_decimal(certificate.constant_upper),
        "theorem_terms": certificate.theorem_terms,
        "expanded_commutator_keys": certificate.expanded_commutator_keys,
    }
    return {**unsigned, "identity_digest": _canonical_digest(unsigned)}


def calibrate_rectangles(
    sizes: tuple[tuple[int, int], ...],
    steps: tuple[int, ...],
    *,
    progress: bool = False,
) -> dict[str, object]:
    """Evaluate deterministic dense S4 errors on sorted open rectangles.

    ``certified_upper`` is rebuilt directly for each finite open rectangle by
    exact Pauli commutators and outward rational Suzuki coefficients, then
    capped by the universal distance-two bound for unitaries.  No periodic
    site density is transferred to the open boundary.
    """

    normalized_sizes = tuple(sorted(set(sizes)))
    normalized_steps = tuple(sorted(set(steps)))
    if not normalized_sizes or not normalized_steps:
        raise ValueError("at least one size and step count are required")
    if any(width < 2 or height < 2 for width, height in normalized_sizes):
        raise ValueError("calibration rectangles require width,height >= 2")
    if any(isinstance(step, bool) or not isinstance(step, int) or step < 1 for step in normalized_steps):
        raise ValueError("calibration step counts must be positive integers")

    rows: list[dict[str, object]] = []
    certificate_records: list[dict[str, object]] = []
    for width, height in normalized_sizes:
        open_record = _open_certificate_record(width, height)
        certificate_records.append(open_record)
        constant_upper = Fraction(*open_record["constant_upper_exact"])
        groups, fragment_matrices, hamiltonian = open_rectangle_dense_fragments(
            width, height
        )
        n_sites = width * height
        exact = expm(-1j * hamiltonian)
        bonds = tuple(bond for group in groups for bond in group)
        for step_count in normalized_steps:
            product_step = fourth_order_product_step(
                fragment_matrices,
                time=1.0,
                steps=step_count,
            )
            approximate = np.linalg.matrix_power(product_step, step_count)
            actual = float(np.linalg.norm(exact - approximate, ord=2))
            if not math.isfinite(actual):
                raise ArithmeticError("dense operator error is not finite")
            published = constant_upper / step_count**4
            certified = min(Fraction(2), published)
            if actual > float(certified):
                raise ArithmeticError(
                    f"dense error exceeds certified upper at {width}x{height}, "
                    f"r={step_count}"
                )
            rows.append(
                {
                    "model": "open_rectangular_spin_half_isotropic_heisenberg",
                    "width": width,
                    "height": height,
                    "n_sites": n_sites,
                    "boundary": "open",
                    "normalization": NORMALIZATION,
                    "fragment_count": 4,
                    "nonempty_fragment_count": sum(bool(group) for group in groups),
                    "unique_bond_count": len(set(bonds)),
                    "steps": step_count,
                    "actual_operator_error": _diagnostic_decimal(actual),
                    "actual_operator_error_status": (
                        "dense_float64_diagnostic_not_a_certified_interval"
                    ),
                    "certified_upper": _outward_decimal(certified),
                    "certified_upper_exact": _fraction_pair(certified),
                    "published_upper": _outward_decimal(published),
                    "published_upper_exact": _fraction_pair(published),
                    "certified_upper_method": (
                        "direct_finite_open_published_triangle_exact_pauli_l1_"
                        "capped_by_unitary_distance"
                    ),
                    "open_triangle_certificate_digest": open_record[
                        "identity_digest"
                    ],
                    "candidate_grouped_open_boundary_transfer_claimed": False,
                }
            )
            if progress:
                print(
                    f"calibrated={width}x{height} steps={step_count} "
                    f"actual={_diagnostic_decimal(actual)}",
                    flush=True,
                )
    provenance = {
        "scientific_dependency_sha256": _scientific_dependency_hashes(),
        "reference_certificate": _reference_certificate_record(),
        "runtime_environment": _runtime_environment(),
    }
    return {
        "schema_version": 2,
        "kind": "paper_a_nondegenerate_ed_calibration",
        "normalization": NORMALIZATION,
        "time": "1",
        "decimal_policy": {
            "certified_bounds": (
                "16_significant_digits_Decimal_ROUND_CEILING"
            ),
            "actual_operator_error": (
                "9_significant_digits_float64_diagnostic_display"
            ),
        },
        "actual_error_status": (
            "dense_float64_diagnostic_not_a_certified_interval"
        ),
        "candidate_grouped_open_boundary_transfer_claimed": False,
        "open_rectangle_certificates": certificate_records,
        "provenance": provenance,
        "provenance_digest": _canonical_digest(provenance),
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=PAPER_ROOT / "benchmarks/paper-a/ed-calibration.json",
    )
    arguments = parser.parse_args()
    payload = calibrate_rectangles(FROZEN_SIZES, FROZEN_STEPS, progress=True)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    print(f"output={arguments.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
