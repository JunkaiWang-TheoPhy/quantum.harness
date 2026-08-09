from fractions import Fraction

import pytest

from trottercert.processed_kernels import (
    alternating_kernel_stages,
    groups_for_repetitions,
    maximum_steps_for_group_budget,
    merged_groups_per_step,
    published_effective_order_six_kernel,
)


@pytest.mark.parametrize(
    ("name", "stage_count", "beat_steps", "fivefold_steps"),
    [("s8", 49, 59, 49), ("s10", 61, 47, 39), ("s11", 67, 43, 35)],
)
def test_published_kernel_resource_gates(
    name: str,
    stage_count: int,
    beat_steps: int,
    fivefold_steps: int,
) -> None:
    kernel = published_effective_order_six_kernel(name)
    stages = alternating_kernel_stages(kernel)
    assert len(stages) == stage_count
    assert sum((stage.coefficient for stage in stages), Fraction()) == 4
    assert merged_groups_per_step(kernel) == stage_count
    assert maximum_steps_for_group_budget(kernel, 2850) == beat_steps
    assert maximum_steps_for_group_budget(kernel, 2358) == fivefold_steps


def test_registry_preserves_exact_decimal_provenance() -> None:
    kernel = published_effective_order_six_kernel("s10")
    assert kernel.seed_coefficients[0] == Fraction(
        "0.100838384835000970361478569216"
    )
    assert sum(kernel.composition_coefficients, Fraction()) == 1
    assert kernel.source_doi == "10.5281/zenodo.10814897"
    assert kernel.source_file == "metemprar.py"
    assert kernel.source_md5 == "35a37af90096f371e7c82366cc7fad8b"


def test_registry_rejects_unknown_kernel() -> None:
    with pytest.raises(ValueError, match="unknown processed kernel"):
        published_effective_order_six_kernel("s9")


def test_resource_arithmetic_rejects_nonintegral_inputs() -> None:
    kernel = published_effective_order_six_kernel("s10")
    with pytest.raises(ValueError, match="positive integer"):
        groups_for_repetitions(kernel, 39.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="positive integer"):
        maximum_steps_for_group_budget(kernel, 2358.0)  # type: ignore[arg-type]
