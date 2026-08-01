from __future__ import annotations

import pytest
from sympy.combinatorics.free_groups import free_group

from trottercert.corrector_cost import (
    CorrectorArchitecture,
    telescoping_defect_order,
    total_blocks,
)


def test_symplectic_endpoint_cost_is_additive() -> None:
    architecture = CorrectorArchitecture(
        "symplectic",
        kernel_blocks=31,
        endpoint_blocks=2,
        repeated_blocks_per_step=0,
    )

    assert architecture.spectrum_preserving_conjugation
    assert total_blocks(architecture, 95) == 31 * 95 + 2
    assert telescoping_defect_order(architecture.kind) is None


def test_symmetric_correction_cost_is_linear() -> None:
    architecture = CorrectorArchitecture(
        "symmetric",
        kernel_blocks=31,
        endpoint_blocks=0,
        repeated_blocks_per_step=2,
    )

    assert not architecture.spectrum_preserving_conjugation
    assert total_blocks(architecture, 95) == (31 + 2) * 95
    assert telescoping_defect_order(architecture.kind) == 2


def test_composite_cost_can_include_endpoint_and_repeated_blocks() -> None:
    architecture = CorrectorArchitecture(
        "composite",
        kernel_blocks=7,
        endpoint_blocks=2,
        repeated_blocks_per_step=3,
    )

    assert total_blocks(architecture, 4) == 42
    assert telescoping_defect_order(architecture.kind) == 2


@pytest.mark.parametrize(
    ("field", "value", "exception"),
    (
        ("kernel_blocks", True, TypeError),
        ("endpoint_blocks", False, TypeError),
        ("repeated_blocks_per_step", -1, ValueError),
    ),
)
def test_architecture_rejects_invalid_counts(
    field: str, value: int, exception: type[Exception]
) -> None:
    values = {
        "kernel_blocks": 1,
        "endpoint_blocks": 0,
        "repeated_blocks_per_step": 0,
    }
    values[field] = value

    with pytest.raises(exception):
        CorrectorArchitecture("symplectic", **values)


def test_architecture_rejects_unknown_kind() -> None:
    with pytest.raises(ValueError, match="unknown corrector architecture"):
        CorrectorArchitecture("unknown", 1, 0, 0)  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="unknown corrector architecture"):
        telescoping_defect_order("unknown")


@pytest.mark.parametrize(("steps", "exception"), ((0, ValueError), (-1, ValueError), (True, TypeError)))
def test_total_blocks_rejects_invalid_steps(
    steps: int, exception: type[Exception]
) -> None:
    architecture = CorrectorArchitecture("symplectic", 1, 0, 0)
    with pytest.raises(exception):
        total_blocks(architecture, steps)


def test_two_step_conjugation_telescopes_as_a_free_word() -> None:
    _, processor, step = free_group("processor, step")
    processed_step = processor * step * processor**-1

    assert processed_step**2 == processor * step**2 * processor**-1


def test_two_step_nonconjugate_correction_does_not_generically_telescope() -> None:
    _, correction, step = free_group("correction, step")

    assert (correction * step) ** 2 != correction**2 * step**2
