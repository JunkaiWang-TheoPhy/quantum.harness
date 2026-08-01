"""Frozen setup primitives for the grouped finite-torus XXZ compiler.

This module intentionally contains only the setup, schedule, and canonical
serialization surface established by implementation-plan Tasks 1--3.  The
scientific ledger and artifact build paths remain unavailable until their
later tasks are implemented and independently verified.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal, Mapping

from .cubic_field import Cubic, fourth_order_suzuki_cubic_stages


XXZStatus = Literal["certified", "unsupported", "inconclusive"]
Boundary = Literal["periodic"]

FORMULA_IDENTIFIER = "five_copy_fourth_order_suzuki_four_matchings"
NORMALIZATION = "(XX+YY+delta*ZZ)/4"
PRIMARY_METRIC = "merged_group_exponentials"
PILOT_LENGTH = 4
STAGE_COUNT = 31


def fraction_pair(value: Fraction) -> list[int]:
    """Return the unique positive-denominator pair for an exact fraction."""

    if not isinstance(value, Fraction):
        raise TypeError("value must be an exact Fraction")
    return [value.numerator, value.denominator]


def strict_fraction_pair(value: object, *, field: str) -> Fraction:
    """Decode a canonical JSON fraction pair without implicit coercions."""

    if not isinstance(field, str) or not field:
        raise TypeError("field must be a nonempty string")
    if not isinstance(value, list) or len(value) != 2:
        raise TypeError(f"{field} must be a two-element JSON array")
    numerator, denominator = value
    if (
        isinstance(numerator, bool)
        or isinstance(denominator, bool)
        or not isinstance(numerator, int)
        or not isinstance(denominator, int)
    ):
        raise TypeError(f"{field} numerator and denominator must be integers")
    if denominator <= 0:
        raise ValueError(f"{field} denominator must be positive")
    result = Fraction(numerator, denominator)
    if value != fraction_pair(result):
        raise ValueError(f"{field} must be a canonical reduced fraction pair")
    return result


def canonical_json_bytes(payload: object) -> bytes:
    """Encode compact, sorted, UTF-8 JSON with one trailing newline."""

    return (
        json.dumps(
            payload,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


@dataclass(frozen=True, slots=True)
class XXZCompileSpec:
    """Immutable physical and resource setup for one grouped XXZ compile."""

    delta: Fraction
    length: int
    boundary: Boundary
    time: Fraction
    tolerance: Fraction
    normalization: str
    formula_identifier: str
    stage_count: int
    primary_metric: str
    scaling_certificate_sha256: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.delta, Fraction):
            raise TypeError("delta must be an exact Fraction")
        if isinstance(self.length, bool) or not isinstance(self.length, int):
            raise TypeError("length must be an integer")
        if self.length != PILOT_LENGTH:
            raise ValueError(
                "only the frozen 4 by 4 pilot length is enabled before the "
                "Task 10 scaling verifier"
            )
        if self.boundary != "periodic":
            raise ValueError("boundary must be literal 'periodic'")
        if not isinstance(self.time, Fraction) or self.time != 1:
            raise ValueError("time must be the exact Fraction 1")
        if (
            not isinstance(self.tolerance, Fraction)
            or self.tolerance != Fraction(1, 10**6)
        ):
            raise ValueError("tolerance must be the exact Fraction 1/1000000")
        if self.normalization != NORMALIZATION:
            raise ValueError("XXZ normalization is frozen")
        if self.formula_identifier != FORMULA_IDENTIFIER:
            raise ValueError("formula identifier is frozen")
        if isinstance(self.stage_count, bool) or self.stage_count != STAGE_COUNT:
            raise ValueError("stage count must be exactly 31")
        if self.primary_metric != PRIMARY_METRIC:
            raise ValueError("primary resource metric is frozen")
        if self.scaling_certificate_sha256 is not None:
            raise ValueError("the 4 by 4 pilot must not cite a scaling certificate")

    @classmethod
    def pilot(cls, delta: Fraction) -> XXZCompileSpec:
        if not isinstance(delta, Fraction):
            raise TypeError("delta must be an exact Fraction")
        return cls(
            delta=delta,
            length=PILOT_LENGTH,
            boundary="periodic",
            time=Fraction(1),
            tolerance=Fraction(1, 10**6),
            normalization=NORMALIZATION,
            formula_identifier=FORMULA_IDENTIFIER,
            stage_count=STAGE_COUNT,
            primary_metric=PRIMARY_METRIC,
        )

    @classmethod
    def periodic(
        cls,
        delta: Fraction,
        length: int,
        scaling_certificate: Mapping[str, object],
    ) -> XXZCompileSpec:
        """Reserve nonpilot construction and fail closed until Task 10.

        Merely presenting a mapping is not geometry authorization.  Task 10
        replaces this guard with semantic scaling-certificate verification.
        """

        if not isinstance(delta, Fraction):
            raise TypeError("delta must be an exact Fraction")
        if isinstance(length, bool) or not isinstance(length, int):
            raise TypeError("length must be an integer")
        if length < PILOT_LENGTH or length % 2:
            raise ValueError("periodic XXZ length must be even and at least four")
        if not isinstance(scaling_certificate, Mapping):
            raise TypeError("scaling_certificate must be a verified mapping")
        raise NotImplementedError(
            "nonpilot periodic construction remains fail-closed until Task 10"
        )


@dataclass(frozen=True, slots=True)
class StageRecord:
    """Canonical exact record for one merged Suzuki stage."""

    fragment_index: int
    coefficient_coordinates: tuple[Fraction, Fraction, Fraction]

    def __post_init__(self) -> None:
        if (
            isinstance(self.fragment_index, bool)
            or not isinstance(self.fragment_index, int)
        ):
            raise TypeError("stage fragment index must be an integer")
        if not 0 <= self.fragment_index < 4:
            raise ValueError("stage fragment index must lie in 0..3")
        coordinates = self.coefficient_coordinates
        if not isinstance(coordinates, tuple) or len(coordinates) != 3:
            raise TypeError("stage coefficient coordinates must be a triple")
        if any(not isinstance(value, Fraction) for value in coordinates):
            raise TypeError("stage coefficient coordinates must be exact Fractions")

    @property
    def coefficient(self) -> Cubic:
        return Cubic(*self.coefficient_coordinates)


def xxz_suzuki_schedule() -> tuple[StageRecord, ...]:
    """Return the frozen exact 31-stage fourth-order four-fragment schedule."""

    return tuple(
        StageRecord(
            fragment_index=stage.fragment_index,
            coefficient_coordinates=(
                stage.coefficient.a0,
                stage.coefficient.a1,
                stage.coefficient.a2,
            ),
        )
        for stage in fourth_order_suzuki_cubic_stages(4)
    )


def verify_xxz_suzuki_schedule(schedule: tuple[StageRecord, ...]) -> None:
    """Reject any deviation from the frozen exact schedule and its identities."""

    if not isinstance(schedule, tuple) or any(
        not isinstance(stage, StageRecord) for stage in schedule
    ):
        raise TypeError("XXZ schedule must be a tuple of StageRecord values")
    if len(schedule) != STAGE_COUNT:
        raise ValueError("XXZ schedule must contain exactly 31 stages")
    if schedule != tuple(reversed(schedule)):
        raise ValueError("XXZ schedule must be palindromic")
    for fragment in range(4):
        total = sum(
            (
                stage.coefficient
                for stage in schedule
                if stage.fragment_index == fragment
            ),
            Cubic.zero(),
        )
        if total != Cubic.one():
            raise ValueError("XXZ schedule coefficient sum is not one")
    if schedule != xxz_suzuki_schedule():
        raise ValueError("XXZ schedule differs from the frozen cubic coordinates")


def replay_schedule_resources(
    schedule: tuple[StageRecord, ...],
    steps: int,
) -> int:
    """Count merged fragment exponentials by replaying stage adjacency."""

    verify_xxz_suzuki_schedule(schedule)
    if isinstance(steps, bool) or not isinstance(steps, int):
        raise TypeError("steps must be an integer")
    if steps < 1:
        raise ValueError("steps must be positive")

    within_step_changes = sum(
        left.fragment_index != right.fragment_index
        for left, right in zip(schedule, schedule[1:])
    )
    boundary_changes = schedule[-1].fragment_index != schedule[0].fragment_index
    return 1 + steps * within_step_changes + (steps - 1) * boundary_changes
