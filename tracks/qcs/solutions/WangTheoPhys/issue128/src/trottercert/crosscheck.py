from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial

import numpy as np
from scipy.linalg import expm

from .algebra import to_dense
from .hamiltonian import four_matching_fragments, fragment_from_bonds
from .higher_order import (
    _multinomial,
    fourth_order_suzuki_stages,
    nested_commutator,
    weak_compositions,
)
from .intervals import RationalInterval
from .lattice import SquareLattice
from .rigorous_fourth import fourth_order_suzuki_interval_stages

# A short outward-rounded summary of the schema-v3 bound at N=144, r=97.
# The exact rational value remains in certificates/issue128-certificate.json.
REFINED_BOUND_N144 = Fraction(996, 10**9)
REFINED_STEPS = 97
# The D5-integrated candidate step count is retained for dense diagnostics;
# its periodic bound is not transferred to finite open rectangles.
REFINED_D5_STEPS = 95


@dataclass(frozen=True, slots=True)
class OpenRectanglePublishedTriangleCertificate:
    """Exact finite-open instantiation of the published S4 triangle bound."""

    width: int
    height: int
    center: int
    coefficient_interval_decimal_digits: int
    root_interval: RationalInterval
    constant_upper: Fraction
    theorem_terms: int
    expanded_commutator_keys: int


def open_rectangular_matchings(
    width: int,
    height: int,
) -> tuple[tuple[tuple[int, int], ...], ...]:
    """Four nonoverlapping color classes for an open rectangular grid."""

    if width < 2 or height < 2:
        raise ValueError("open rectangular cross-check requires both dimensions >= 2")

    def site(x: int, y: int) -> int:
        return y * width + x

    groups: list[list[tuple[int, int]]] = [[], [], [], []]
    for y in range(height):
        for x in range(width - 1):
            groups[x % 2].append((site(x, y), site(x + 1, y)))
    for y in range(height - 1):
        for x in range(width):
            groups[2 + y % 2].append((site(x, y), site(x, y + 1)))
    return tuple(tuple(group) for group in groups)


def open_rectangle_dense_fragments(
    width: int,
    height: int,
) -> tuple[
    tuple[tuple[tuple[int, int], ...], ...],
    tuple[np.ndarray, ...],
    np.ndarray,
]:
    """Build the existing normalized four-matching open-grid matrices."""

    groups = open_rectangular_matchings(width, height)
    n_sites = width * height
    fragments = tuple(fragment_from_bonds(group) for group in groups)
    matrices = tuple(to_dense(fragment, n_sites) for fragment in fragments)
    hamiltonian = sum(matrices, np.zeros_like(matrices[0]))
    return groups, matrices, hamiltonian


@lru_cache(maxsize=32)
def open_rectangle_published_triangle_certificate(
    width: int,
    height: int,
    *,
    center: int = 20,
    decimal_digits: int = 18,
) -> OpenRectanglePublishedTriangleCertificate:
    """Rebuild a strict published-theorem constant on one open rectangle.

    Unlike the periodic certificate, this finite-volume proof never transfers
    a per-site density.  It builds the four open-boundary matching fragments,
    evaluates every depth-five Pauli commutator exactly, and combines their
    exact axis-l1 norms with outward interval weights for the Suzuki
    coefficients.  The resulting ``constant_upper`` gives the global bound
    ``constant_upper * time**5 / steps**4`` at total evolution time ``time``.
    """

    groups = open_rectangular_matchings(width, height)
    fragments = tuple(fragment_from_bonds(group) for group in groups)
    stages_left, root = fourth_order_suzuki_interval_stages(
        4,
        decimal_digits=decimal_digits,
    )
    stages = tuple(reversed(stages_left))
    count = len(stages)
    if not 1 <= center <= count:
        raise ValueError("center must use one-based stage indexing")

    order = 4
    weights: dict[tuple[int, ...], Fraction] = {}
    theorem_terms = 0

    def collect(
        j: int,
        indices: tuple[int, ...],
        composition: tuple[int, ...],
    ) -> None:
        nonlocal theorem_terms
        theorem_terms += 1
        outer: list[int] = []
        scalar = RationalInterval.point(_multinomial(order, composition))
        for stage_index, power in zip(indices, composition):
            stage = stages[stage_index - 1]
            scalar *= stage.coefficient**power
            outer.extend([stage.fragment_index] * power)
        scalar_upper = scalar.abs_upper()
        for base_index in range(1, j):
            base_stage = stages[base_index - 1]
            key = tuple(outer) + (base_stage.fragment_index,)
            weights[key] = weights.get(key, Fraction()) + (
                scalar_upper * base_stage.coefficient.abs_upper()
            )

    for j in range(2, center + 1):
        indices = tuple(range(center, j - 1, -1))
        for composition in weak_compositions(order, len(indices)):
            if composition[-1]:
                collect(j, indices, composition)
    for j in range(center + 1, count + 1):
        indices = tuple(range(center + 1, j + 1))
        for composition in weak_compositions(order, len(indices)):
            if composition[-1]:
                collect(j, indices, composition)

    commutator_cache = {}
    total = Fraction()
    for key, weight in weights.items():
        operator = nested_commutator(fragments, key, commutator_cache)
        total += weight * operator.exact_axis_l1()

    return OpenRectanglePublishedTriangleCertificate(
        width=width,
        height=height,
        center=center,
        coefficient_interval_decimal_digits=decimal_digits,
        root_interval=root,
        constant_upper=total / factorial(order + 1),
        theorem_terms=theorem_terms,
        expanded_commutator_keys=len(weights),
    )


def fourth_order_product_step(
    fragment_matrices: tuple[np.ndarray, ...],
    *,
    time: float,
    steps: int,
) -> np.ndarray:
    """Return one S4 macro-step using the repository's four-fragment order."""

    if len(fragment_matrices) != 4:
        raise ValueError("four fragment matrices are required")
    if steps < 1 or time <= 0:
        raise ValueError("time and steps must be positive")
    shape = fragment_matrices[0].shape
    if len(shape) != 2 or shape[0] != shape[1]:
        raise ValueError("fragment matrices must be square")
    if any(matrix.shape != shape for matrix in fragment_matrices):
        raise ValueError("fragment matrices must have a common shape")

    product_step = np.eye(shape[0], dtype=np.complex128)
    exponential_cache: dict[tuple[int, str], np.ndarray] = {}
    for stage in fourth_order_suzuki_stages(4):
        coefficient_key = str(stage.coefficient)
        key = stage.fragment_index, coefficient_key
        exponential = exponential_cache.get(key)
        if exponential is None:
            exponential = expm(
                -1j
                * float(stage.coefficient)
                * (time / steps)
                * fragment_matrices[stage.fragment_index]
            )
            exponential_cache[key] = exponential
        product_step = product_step @ exponential
    return product_step


def operator_error(
    hamiltonian: np.ndarray,
    product_step: np.ndarray,
    *,
    time: float,
    steps: int,
) -> float:
    """Dense spectral-norm error for a fixed one-step product matrix."""

    if steps < 1 or time <= 0:
        raise ValueError("time and steps must be positive")
    if hamiltonian.shape != product_step.shape:
        raise ValueError("Hamiltonian and product step must have the same shape")
    exact = expm(-1j * time * hamiltonian)
    approximate = np.linalg.matrix_power(product_step, steps)
    return float(np.linalg.norm(exact - approximate, ord=2))


def small_exact_crosscheck(
    length: int,
    tolerance: Fraction,
) -> dict[str, object]:
    if length == 2:
        groups = (
            ((0, 1), (2, 3)),
            ((0, 1), (2, 3)),
            ((0, 2), (1, 3)),
            ((0, 2), (1, 3)),
        )
        fragments = tuple(fragment_from_bonds(group) for group in groups)
        n_sites = 4
    else:
        lattice = SquareLattice(length)
        fragments = four_matching_fragments(lattice)
        n_sites = lattice.n_sites
    matrices = [to_dense(fragment, n_sites) for fragment in fragments]
    hamiltonian = sum(matrices, np.zeros_like(matrices[0]))
    steps = REFINED_STEPS
    delta = 1 / steps
    one_step = np.eye(1 << n_sites, dtype=np.complex128)
    for stage in fourth_order_suzuki_stages(4):
        one_step = one_step @ expm(
            -1j * float(stage.coefficient) * delta * matrices[stage.fragment_index]
        )
    approximation = np.linalg.matrix_power(one_step, steps)
    exact = expm(-1j * hamiltonian)
    error = float(np.linalg.norm(approximation - exact, ord=2))
    bound = float(REFINED_BOUND_N144 * Fraction(n_sites, 144))
    return {
        "length": length,
        "n_sites": n_sites,
        "model_note": (
            "degenerate 2x2 periodic algebra sanity check for schema-v3 r=97"
            if length == 2
            else "periodic square lattice"
        ),
        "steps": steps,
        "empirical_operator_norm_error": error,
        "certified_upper_bound": bound,
        "requested_tolerance": float(tolerance),
        "bound_dominates_empirical_error": bound >= error,
        "bound_meets_tolerance": bound <= float(tolerance),
    }


def small_open_exact_crosscheck(
    width: int,
    height: int,
    tolerance: Fraction,
) -> dict[str, object]:
    """Dense diagnostic checked against a directly rebuilt open bound."""

    groups, matrices, hamiltonian = open_rectangle_dense_fragments(width, height)
    n_sites = width * height
    steps = REFINED_D5_STEPS
    one_step = fourth_order_product_step(matrices, time=1.0, steps=steps)
    error = operator_error(hamiltonian, one_step, time=1.0, steps=steps)
    open_certificate = open_rectangle_published_triangle_certificate(
        width,
        height,
    )
    certified_bound = open_certificate.constant_upper / steps**4
    bonds = tuple(bond for group in groups for bond in group)
    return {
        "width": width,
        "height": height,
        "n_sites": n_sites,
        "boundary": "open",
        "model_note": "nondegenerate open-grid diagnostic with finite-open proof",
        "steps": steps,
        "bond_count": len(bonds),
        "unique_bond_count": len(set(bonds)),
        "nonempty_fragment_count": sum(bool(group) for group in groups),
        "empirical_operator_norm_error": error,
        "certified_upper_bound": float(certified_bound),
        "certified_upper_exact": [
            certified_bound.numerator,
            certified_bound.denominator,
        ],
        "certified_upper_method": (
            "direct_finite_open_published_triangle_exact_pauli_l1"
        ),
        "requested_tolerance": float(tolerance),
        "bound_dominates_empirical_error": float(certified_bound) >= error,
        "bound_meets_tolerance": certified_bound <= tolerance,
        "empirical_error_is_diagnostic_only": True,
    }
