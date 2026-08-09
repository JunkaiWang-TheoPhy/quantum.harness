from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .processed_kernels import published_effective_order_six_kernel
from .processed_local_bounds import (
    build_processed_local_bounds,
    processor_distance_upper,
)


@dataclass(frozen=True)
class ProcessedFiniteStepInputs:
    kernel_name: str
    steps: int
    n_sites: int
    tolerance: Fraction
    degree7_site_l1: Fraction | None
    higher_order_remainder: Fraction | None
    local_log_theorem_id: str | None

    def __post_init__(self) -> None:
        published_effective_order_six_kernel(self.kernel_name)
        if (
            isinstance(self.steps, bool)
            or not isinstance(self.steps, int)
            or self.steps < 1
        ):
            raise ValueError("step count must be a positive integer")
        if (
            isinstance(self.n_sites, bool)
            or not isinstance(self.n_sites, int)
            or self.n_sites < 1
        ):
            raise ValueError("site count must be a positive integer")
        if isinstance(self.tolerance, bool) or not isinstance(
            self.tolerance,
            (Fraction, int),
        ):
            raise ValueError("tolerance must be an exact rational value")
        tolerance = Fraction(self.tolerance)
        if tolerance <= 0:
            raise ValueError("tolerance must be positive")
        object.__setattr__(self, "tolerance", tolerance)
        if self.degree7_site_l1 is None:
            raise ValueError("physical degree-seven site bound is required")
        if self.higher_order_remainder is None:
            raise ValueError("certified higher-order remainder is required")
        if any(
            isinstance(value, bool) or not isinstance(value, (Fraction, int))
            for value in (
                self.degree7_site_l1,
                self.higher_order_remainder,
            )
        ):
            raise ValueError("proof inputs must be exact rational values")
        degree7 = Fraction(self.degree7_site_l1)
        remainder = Fraction(self.higher_order_remainder)
        if degree7 < 0 or remainder < 0:
            raise ValueError("proof inputs must be nonnegative")
        object.__setattr__(self, "degree7_site_l1", degree7)
        object.__setattr__(self, "higher_order_remainder", remainder)
        if not isinstance(self.local_log_theorem_id, str) or not (
            self.local_log_theorem_id.strip()
        ):
            raise ValueError("a local-log theorem identifier is required")


@dataclass(frozen=True)
class ProcessedFiniteStepBudget:
    kernel_name: str
    steps: int
    degree3_residual: Fraction
    degree5_residual: Fraction
    degree7: Fraction
    higher_order_remainder: Fraction
    total: Fraction
    tolerance: Fraction
    remaining_budget: Fraction
    processor_distance: Fraction
    local_log_theorem_id: str
    accepted: bool


def evaluate_processed_budget(
    inputs: ProcessedFiniteStepInputs,
) -> ProcessedFiniteStepBudget:
    """Evaluate a conditional repeated-local-log error ledger.

    This arithmetic does not prove the supplied local-log theorem or the
    higher-order remainder.  It only combines those required proof inputs
    with exact rationalization-residual and physical degree-seven densities.
    """

    kernel = published_effective_order_six_kernel(inputs.kernel_name)
    bounds = build_processed_local_bounds(kernel)
    n_sites = Fraction(inputs.n_sites)
    degree3 = (
        n_sites * bounds.degree3_residual_site_l1 / inputs.steps**2
    )
    degree5 = (
        n_sites * bounds.degree5_residual_site_l1 / inputs.steps**4
    )
    degree7 = n_sites * inputs.degree7_site_l1 / inputs.steps**6
    total = degree3 + degree5 + degree7 + inputs.higher_order_remainder
    return ProcessedFiniteStepBudget(
        kernel_name=inputs.kernel_name,
        steps=inputs.steps,
        degree3_residual=degree3,
        degree5_residual=degree5,
        degree7=degree7,
        higher_order_remainder=inputs.higher_order_remainder,
        total=total,
        tolerance=inputs.tolerance,
        remaining_budget=inputs.tolerance - total,
        processor_distance=processor_distance_upper(
            bounds,
            inputs.n_sites,
            inputs.steps,
        ),
        local_log_theorem_id=inputs.local_log_theorem_id,
        accepted=total <= inputs.tolerance,
    )
