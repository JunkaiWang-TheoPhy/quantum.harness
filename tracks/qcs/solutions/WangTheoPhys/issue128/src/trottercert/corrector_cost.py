"""Exact block accounting for processor and corrector architectures.

This module deliberately separates word-level telescoping from approximation
order.  ``telescoping_defect_order`` reports the first number of composed
steps at which a generic word fails to telescope; it is not an order in the
step size.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

CorrectorKind = Literal["symplectic", "symmetric", "composite"]
_KINDS = frozenset(("symplectic", "symmetric", "composite"))


def _validate_count(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer, not a boolean")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")


def _validate_kind(kind: str) -> CorrectorKind:
    if kind not in _KINDS:
        allowed = ", ".join(sorted(_KINDS))
        raise ValueError(f"unknown corrector architecture {kind!r}; expected {allowed}")
    return kind  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class CorrectorArchitecture:
    """Block-cost data for one repeated product-formula architecture.

    ``endpoint_blocks`` are paid once for the complete repeated evolution.
    ``repeated_blocks_per_step`` are paid alongside every kernel application.
    A symplectic architecture denotes an exact endpoint conjugation and is
    therefore spectrum preserving at the kernel level.
    """

    kind: CorrectorKind
    kernel_blocks: int
    endpoint_blocks: int
    repeated_blocks_per_step: int

    def __post_init__(self) -> None:
        _validate_kind(self.kind)
        _validate_count("kernel_blocks", self.kernel_blocks)
        _validate_count("endpoint_blocks", self.endpoint_blocks)
        _validate_count("repeated_blocks_per_step", self.repeated_blocks_per_step)

    @property
    def spectrum_preserving_conjugation(self) -> bool:
        """Whether the architecture is an exact conjugation of its kernel."""

        return self.kind == "symplectic"


def total_blocks(architecture: CorrectorArchitecture, steps: int) -> int:
    """Return the exact number of primitive blocks for ``steps`` repetitions."""

    if not isinstance(architecture, CorrectorArchitecture):
        raise TypeError("architecture must be a CorrectorArchitecture")
    if isinstance(steps, bool) or not isinstance(steps, int):
        raise TypeError("steps must be an integer, not a boolean")
    if steps <= 0:
        raise ValueError("steps must be positive")
    return architecture.endpoint_blocks + steps * (
        architecture.kernel_blocks + architecture.repeated_blocks_per_step
    )


def telescoping_defect_order(kind: str) -> int | None:
    """Return the first generic composition length with a word defect.

    Exact endpoint conjugations telescope for every number of steps and return
    ``None``.  Generic symmetric and composite nonconjugate corrections fail
    the corresponding word identity already at two steps and return ``2``.
    This value is a composition length, not a power of the step size.
    """

    checked_kind = _validate_kind(kind)
    if checked_kind == "symplectic":
        return None
    return 2
