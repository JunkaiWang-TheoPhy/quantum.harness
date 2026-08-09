from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial

from .effective_processor import audit_effective_order_six
from .processed_kernels import ProcessedKernel
from .rational_lie_local import (
    local_terms_l1,
    project_lie_polynomial_to_cell,
)


@dataclass(frozen=True)
class ProcessedLocalBounds:
    kernel_name: str
    r2_cell_l1: Fraction
    r2_site_l1: Fraction
    r4_cell_l1: Fraction
    r4_site_l1: Fraction
    degree3_residual_cell_l1: Fraction
    degree3_residual_site_l1: Fraction
    degree5_residual_cell_l1: Fraction
    degree5_residual_site_l1: Fraction
    r2_term_count: int
    r4_term_count: int
    degree3_residual_term_count: int
    degree5_residual_term_count: int


@lru_cache(maxsize=8)
def build_processed_local_bounds(
    kernel: ProcessedKernel,
) -> ProcessedLocalBounds:
    audit = audit_effective_order_six(kernel)
    r2_terms = project_lie_polynomial_to_cell(audit.r2.polynomial)
    r4_terms = project_lie_polynomial_to_cell(audit.r4.polynomial)
    degree3_terms = project_lie_polynomial_to_cell(
        audit.processed_degree_three
    )
    degree5_terms = project_lie_polynomial_to_cell(
        audit.processed_degree_five
    )
    r2_cell = local_terms_l1(r2_terms)
    r4_cell = local_terms_l1(r4_terms)
    degree3_cell = local_terms_l1(degree3_terms)
    degree5_cell = local_terms_l1(degree5_terms)
    return ProcessedLocalBounds(
        kernel_name=kernel.name,
        r2_cell_l1=r2_cell,
        r2_site_l1=r2_cell / 4,
        r4_cell_l1=r4_cell,
        r4_site_l1=r4_cell / 4,
        degree3_residual_cell_l1=degree3_cell,
        degree3_residual_site_l1=degree3_cell / 4,
        degree5_residual_cell_l1=degree5_cell,
        degree5_residual_site_l1=degree5_cell / 4,
        r2_term_count=len(r2_terms),
        r4_term_count=len(r4_terms),
        degree3_residual_term_count=len(degree3_terms),
        degree5_residual_term_count=len(degree5_terms),
    )


def rational_exp_upper(
    value: Fraction | int,
    *,
    terms: int = 32,
    quantization_digits: int = 60,
) -> Fraction:
    """Return an exact rational upper enclosure of ``exp(value)``."""

    if isinstance(value, bool):
        raise ValueError("exponential input must be nonnegative")
    x = Fraction(value)
    if x < 0:
        raise ValueError("exponential input must be nonnegative")
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 0:
        raise ValueError("exponential term count must be nonnegative")
    if (
        isinstance(quantization_digits, bool)
        or not isinstance(quantization_digits, int)
        or quantization_digits < 1
    ):
        raise ValueError("exponential quantization digits must be positive")
    if not x:
        return Fraction(1)
    grid = 10**quantization_digits

    def upward_grid(rational: Fraction) -> Fraction:
        numerator = (
            rational.numerator * grid + rational.denominator - 1
        ) // rational.denominator
        return Fraction(numerator, grid)

    # Exact processor coordinates can have very large denominators.  An
    # outward decimal grid keeps the proof object portable without weakening
    # rigor: exp is monotone on the nonnegative real axis.
    x_upper = upward_grid(x)
    ratio = x_upper / (terms + 2)
    if ratio >= 1:
        raise ValueError("exponential tail ratio must remain below one")
    partial = sum(
        (
            x_upper**power / factorial(power)
            for power in range(terms + 1)
        ),
        Fraction(),
    )
    next_term = x_upper ** (terms + 1) / factorial(terms + 1)
    return upward_grid(partial + next_term / (1 - ratio))


def processor_generator_norm_upper(
    bounds: ProcessedLocalBounds,
    n_sites: int,
    steps: int,
) -> Fraction:
    if isinstance(n_sites, bool) or not isinstance(n_sites, int) or n_sites < 1:
        raise ValueError("site count must be a positive integer")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
        raise ValueError("step count must be a positive integer")
    return Fraction(n_sites) * (
        bounds.r2_site_l1 / steps**2
        + bounds.r4_site_l1 / steps**4
    )


def processor_distance_upper(
    bounds: ProcessedLocalBounds,
    n_sites: int,
    steps: int,
) -> Fraction:
    generator_norm = processor_generator_norm_upper(
        bounds,
        n_sites,
        steps,
    )
    return rational_exp_upper(generator_norm) - 1
