"""Exact finite-torus primitives for the grouped XXZ certificate.

The setup and Suzuki schedule are frozen exactly.  The direct-theorem ledger
keeps the published triangle expansion as a raw multiset and keeps every
unweighted nested-commutator polynomial in its own immutable norm block.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import pairwise
from math import comb, factorial, isqrt
from typing import Literal

from .cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from .higher_order import _multinomial, weak_compositions
from .intervals import RationalInterval
from .lattice import SquareLattice
from .rigorous_fourth import fourth_order_suzuki_interval_stages

XXZStatus = Literal["certified", "unsupported", "inconclusive"]
Boundary = Literal["periodic"]

FORMULA_IDENTIFIER = "five_copy_fourth_order_suzuki_four_matchings"
NORMALIZATION = "(XX+YY+delta*ZZ)/4"
PRIMARY_METRIC = "merged_group_exponentials"
PILOT_LENGTH = 4
STAGE_COUNT = 31
THEOREM_IDENTIFIER = "published_high_order_triangle_v1"
THEOREM_ORDER = 4
THEOREM_CENTER = 20
THEOREM_FACTORIAL_DENOMINATOR = factorial(THEOREM_ORDER + 1)
THEOREM_DECIMAL_DIGITS = 18
THEOREM_DUHAMEL_CONVENTION = "published_triangle_duhamel_1_over_5_factorial"
FINITE_STEP_ERROR_FORMULA = "E_r=K/r^4_for_T=1"
GROUPING_ALGORITHM_IDENTIFIER = "deterministic_pair_only_bitset_v1"


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
        for left, right in pairwise(schedule)
    )
    boundary_changes = schedule[-1].fragment_index != schedule[0].fragment_index
    return 1 + steps * within_step_changes + (steps - 1) * boundary_changes


SymplecticPauli = tuple[int, int]
GaussianInteger = tuple[int, int]
FragmentBonds = tuple[tuple[int, int], ...]


@dataclass(frozen=True, order=True, slots=True)
class IntegerSymplecticTerm:
    """One phase-free Pauli with an integer numerator."""

    x_mask: int
    z_mask: int
    numerator: int

    def __post_init__(self) -> None:
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in (self.x_mask, self.z_mask, self.numerator)
        ):
            raise TypeError("symplectic masks and numerator must be integers")
        if self.x_mask < 0 or self.z_mask < 0:
            raise ValueError("symplectic masks must be nonnegative")
        if self.numerator == 0:
            raise ValueError("zero integer Pauli terms must be omitted")

    @property
    def mask(self) -> SymplecticPauli:
        return self.x_mask, self.z_mask


@dataclass(frozen=True, slots=True)
class WeightedXXZFragment:
    """An XXZ fragment over one exact common denominator."""

    common_denominator: int
    axis_numerators: tuple[int, int, int]
    terms: tuple[IntegerSymplecticTerm, ...]


@dataclass(frozen=True, order=True, slots=True)
class SymplecticCoefficient:
    """Canonical Gaussian-rational coefficient of a phase-free Pauli."""

    x_mask: int
    z_mask: int
    real: Fraction
    imag: Fraction

    def __post_init__(self) -> None:
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in (self.x_mask, self.z_mask)
        ):
            raise TypeError("symplectic masks must be integers")
        if self.x_mask < 0 or self.z_mask < 0:
            raise ValueError("symplectic masks must be nonnegative")
        if not isinstance(self.real, Fraction) or not isinstance(self.imag, Fraction):
            raise TypeError("Pauli coefficients must be exact Fractions")
        if self.real == 0 and self.imag == 0:
            raise ValueError("zero Pauli coefficients must be omitted")

    @property
    def mask(self) -> SymplecticPauli:
        return self.x_mask, self.z_mask

    @property
    def is_zero(self) -> bool:
        return self.real == 0 and self.imag == 0


def _canonical_bonds(bonds: Sequence[tuple[int, int]]) -> FragmentBonds:
    canonical: list[tuple[int, int]] = []
    for bond in bonds:
        if not isinstance(bond, tuple) or len(bond) != 2:
            raise TypeError("each fragment bond must be an integer pair")
        left, right = bond
        if any(isinstance(site, bool) or not isinstance(site, int) for site in bond):
            raise TypeError("bond sites must be integers")
        if left < 0 or right < 0:
            raise ValueError("bond sites must be nonnegative")
        if left == right:
            raise ValueError("XXZ self bonds are not allowed")
        canonical.append((left, right) if left < right else (right, left))
    result = tuple(sorted(canonical))
    if len(result) != len(set(result)):
        raise ValueError("fragment bonds must be unique")
    return result


def weighted_xxz_fragment(
    bonds: Sequence[tuple[int, int]],
    delta: Fraction,
) -> WeightedXXZFragment:
    """Return exact integer XXZ terms over denominator ``4*delta.denominator``."""

    if not isinstance(delta, Fraction):
        raise TypeError("delta must be an exact Fraction")
    canonical_bonds = _canonical_bonds(bonds)
    denominator = 4 * delta.denominator
    axis_numerators = (
        delta.denominator,
        delta.denominator,
        delta.numerator,
    )
    accumulated: dict[SymplecticPauli, int] = {}
    for left, right in canonical_bonds:
        sites = (1 << left) | (1 << right)
        for pauli, numerator in zip(
            ((sites, 0), (sites, sites), (0, sites)),
            axis_numerators,
        ):
            if numerator:
                accumulated[pauli] = accumulated.get(pauli, 0) + numerator
    terms = tuple(
        IntegerSymplecticTerm(x_mask, z_mask, accumulated[(x_mask, z_mask)])
        for x_mask, z_mask in sorted(accumulated)
        if accumulated[(x_mask, z_mask)]
    )
    return WeightedXXZFragment(denominator, axis_numerators, terms)


def _symplectic_anticommutes(
    left: SymplecticPauli,
    right: SymplecticPauli,
) -> bool:
    return bool(
        (
            (left[0] & right[1]).bit_count()
            + (left[1] & right[0]).bit_count()
        )
        & 1
    )


def _symplectic_product_phase(
    left: SymplecticPauli,
    right: SymplecticPauli,
) -> tuple[int, SymplecticPauli]:
    """Return ``e, mask`` such that ``P(left)P(right)=i**e P(mask)``."""

    result = (left[0] ^ right[0], left[1] ^ right[1])
    exponent = (
        (left[0] & left[1]).bit_count()
        + (right[0] & right[1]).bit_count()
        + 2 * (left[1] & right[0]).bit_count()
        - (result[0] & result[1]).bit_count()
    ) % 4
    return exponent, result


def _times_i_power(value: GaussianInteger, exponent: int) -> GaussianInteger:
    real, imag = value
    return (
        (real, imag),
        (-imag, real),
        (-real, -imag),
        (imag, -real),
    )[exponent % 4]


class _FiniteXXZSymplecticEvaluator:
    """Exact sparse evaluator with a common denominator factored by degree."""

    def __init__(
        self,
        fragment_bonds: Sequence[Sequence[tuple[int, int]]],
        delta: Fraction,
    ) -> None:
        if not isinstance(fragment_bonds, Sequence) or not fragment_bonds:
            raise TypeError("fragment_bonds must be a nonempty sequence")
        self.fragments = tuple(
            weighted_xxz_fragment(bonds, delta) for bonds in fragment_bonds
        )
        denominators = {fragment.common_denominator for fragment in self.fragments}
        if len(denominators) != 1:
            raise ArithmeticError("XXZ fragments do not share one denominator")
        self.denominator = denominators.pop()
        self._integer_terms = tuple(
            {term.mask: term.numerator for term in fragment.terms}
            for fragment in self.fragments
        )
        self._cache: dict[tuple[int, ...], dict[SymplecticPauli, GaussianInteger]] = {}

    def integer_numerators(
        self,
        key: tuple[int, ...],
    ) -> dict[SymplecticPauli, GaussianInteger]:
        if not isinstance(key, tuple) or not key:
            raise TypeError("commutator key must be a nonempty tuple")
        if any(
            isinstance(index, bool)
            or not isinstance(index, int)
            or not 0 <= index < len(self.fragments)
            for index in key
        ):
            raise ValueError("commutator key contains an invalid fragment index")
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        if len(key) == 1:
            result = {
                pauli: (numerator, 0)
                for pauli, numerator in self._integer_terms[key[0]].items()
            }
            self._cache[key] = result
            return result

        inner = self.integer_numerators(key[1:])
        result: dict[SymplecticPauli, GaussianInteger] = {}
        for left, left_numerator in self._integer_terms[key[0]].items():
            for right, right_numerator in inner.items():
                if not _symplectic_anticommutes(left, right):
                    continue
                exponent, product = _symplectic_product_phase(left, right)
                phased = _times_i_power(right_numerator, exponent)
                contribution = (
                    2 * left_numerator * phased[0],
                    2 * left_numerator * phased[1],
                )
                previous = result.get(product, (0, 0))
                updated = (
                    previous[0] + contribution[0],
                    previous[1] + contribution[1],
                )
                if updated == (0, 0):
                    result.pop(product, None)
                else:
                    result[product] = updated
        self._cache[key] = result
        return result

    def evaluate(self, key: tuple[int, ...]) -> tuple[SymplecticCoefficient, ...]:
        denominator = self.denominator ** len(key)
        return tuple(
            SymplecticCoefficient(
                x_mask,
                z_mask,
                Fraction(real, denominator),
                Fraction(imag, denominator),
            )
            for (x_mask, z_mask), (real, imag) in sorted(
                self.integer_numerators(key).items()
            )
            if real or imag
        )


def weighted_symplectic_nested_commutator(
    fragment_bonds: Sequence[Sequence[tuple[int, int]]],
    delta: Fraction,
    key: tuple[int, ...],
) -> tuple[SymplecticCoefficient, ...]:
    """Evaluate one finite-site nested commutator without floating point."""

    return _FiniteXXZSymplecticEvaluator(fragment_bonds, delta).evaluate(key)


@dataclass(frozen=True, slots=True)
class RawTheoremRecord:
    """One unaggregated triangle-inequality theorem contribution."""

    record_id: str
    side: Literal["left", "right"]
    partial_sum_index: int
    adjoint_stage_indices: tuple[int, ...]
    composition: tuple[int, ...]
    base_stage_index: int
    fragment_word: tuple[int, ...]
    block_key: tuple[int, ...]
    weight_interval: RationalInterval


@dataclass(frozen=True, slots=True)
class TheoremNormBlock:
    """One unweighted nested-commutator polynomial and its raw weight sum."""

    block_key: tuple[int, ...]
    terms: tuple[SymplecticCoefficient, ...]
    raw_weight_upper: Fraction


@dataclass(frozen=True, slots=True)
class XXZLedger:
    """Finite theorem ledger; bounded prefixes are explicitly non-complete."""

    spec: XXZCompileSpec
    theorem_identifier: str
    order: int
    center: int
    factorial_denominator: int
    duhamel_convention: str
    finite_step_error_formula: str
    max_degree: int
    projected_record_count: int
    record_limit: int | None
    complete: bool
    raw_records: tuple[RawTheoremRecord, ...]
    blocks: tuple[TheoremNormBlock, ...]
    ledger_digest: str


def _positive_absolute_interval(interval: RationalInterval) -> RationalInterval:
    return RationalInterval(interval.abs_lower(), interval.abs_upper())


def _theorem_interval_stages() -> tuple[object, ...]:
    stages_left, _ = fourth_order_suzuki_interval_stages(
        4,
        decimal_digits=THEOREM_DECIMAL_DIGITS,
    )
    stages = tuple(reversed(stages_left))
    if len(stages) != STAGE_COUNT:
        raise ArithmeticError("published theorem schedule is not 31 stages")
    return stages


def _projected_raw_record_count() -> int:
    left = sum(
        (j - 1)
        * comb(
            (THEOREM_CENTER - j + 1) + THEOREM_ORDER - 2,
            THEOREM_ORDER - 1,
        )
        for j in range(2, THEOREM_CENTER + 1)
    )
    right = sum(
        (j - 1)
        * comb(
            (j - THEOREM_CENTER) + THEOREM_ORDER - 2,
            THEOREM_ORDER - 1,
        )
        for j in range(THEOREM_CENTER + 1, STAGE_COUNT + 1)
    )
    return left + right


def _iter_raw_theorem_records() -> Iterator[RawTheoremRecord]:
    stages = _theorem_interval_stages()

    def records_for(
        side: Literal["left", "right"],
        j: int,
        indices: tuple[int, ...],
        composition: tuple[int, ...],
    ) -> Iterator[RawTheoremRecord]:
        outer: list[int] = []
        scalar = RationalInterval.point(_multinomial(THEOREM_ORDER, composition))
        for stage_index, power in zip(indices, composition):
            stage = stages[stage_index - 1]
            scalar *= _positive_absolute_interval(stage.coefficient) ** power
            outer.extend([stage.fragment_index] * power)
        for base_index in range(1, j):
            base_stage = stages[base_index - 1]
            weight = scalar * _positive_absolute_interval(base_stage.coefficient)
            word = tuple(outer) + (base_stage.fragment_index,)
            composition_text = ",".join(str(value) for value in composition)
            yield RawTheoremRecord(
                record_id=f"{side}:{j}:{composition_text}:{base_index}",
                side=side,
                partial_sum_index=j,
                adjoint_stage_indices=indices,
                composition=composition,
                base_stage_index=base_index,
                fragment_word=word,
                block_key=word,
                weight_interval=weight,
            )

    for j in range(2, THEOREM_CENTER + 1):
        indices = tuple(range(THEOREM_CENTER, j - 1, -1))
        for composition in weak_compositions(THEOREM_ORDER, len(indices)):
            if composition[-1]:
                yield from records_for("left", j, indices, composition)
    for j in range(THEOREM_CENTER + 1, STAGE_COUNT + 1):
        indices = tuple(range(THEOREM_CENTER + 1, j + 1))
        for composition in weak_compositions(THEOREM_ORDER, len(indices)):
            if composition[-1]:
                yield from records_for("right", j, indices, composition)


def _fraction_payload(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _interval_payload(value: RationalInterval) -> list[list[int]]:
    return [_fraction_payload(value.lower), _fraction_payload(value.upper)]


def _ledger_payload(ledger: XXZLedger) -> dict[str, object]:
    return {
        "spec": {
            "delta": _fraction_payload(ledger.spec.delta),
            "length": ledger.spec.length,
            "boundary": ledger.spec.boundary,
            "time": _fraction_payload(ledger.spec.time),
            "tolerance": _fraction_payload(ledger.spec.tolerance),
            "normalization": ledger.spec.normalization,
            "formula_identifier": ledger.spec.formula_identifier,
            "stage_count": ledger.spec.stage_count,
            "primary_metric": ledger.spec.primary_metric,
        },
        "theorem_identifier": ledger.theorem_identifier,
        "order": ledger.order,
        "center": ledger.center,
        "factorial_denominator": ledger.factorial_denominator,
        "duhamel_convention": ledger.duhamel_convention,
        "finite_step_error_formula": ledger.finite_step_error_formula,
        "max_degree": ledger.max_degree,
        "projected_record_count": ledger.projected_record_count,
        "record_limit": ledger.record_limit,
        "complete": ledger.complete,
        "raw_records": [
            {
                "record_id": record.record_id,
                "side": record.side,
                "partial_sum_index": record.partial_sum_index,
                "adjoint_stage_indices": list(record.adjoint_stage_indices),
                "composition": list(record.composition),
                "base_stage_index": record.base_stage_index,
                "fragment_word": list(record.fragment_word),
                "block_key": list(record.block_key),
                "weight_interval": _interval_payload(record.weight_interval),
            }
            for record in ledger.raw_records
        ],
        "blocks": [
            {
                "block_key": list(block.block_key),
                "raw_weight_upper": _fraction_payload(block.raw_weight_upper),
                "terms": [
                    [
                        term.x_mask,
                        term.z_mask,
                        _fraction_payload(term.real),
                        _fraction_payload(term.imag),
                    ]
                    for term in block.terms
                ],
            }
            for block in ledger.blocks
        ],
    }


def _ledger_digest(ledger: XXZLedger) -> str:
    return sha256(canonical_json_bytes(_ledger_payload(ledger))).hexdigest()


def _finite_fragment_bonds(spec: XXZCompileSpec) -> tuple[FragmentBonds, ...]:
    lattice = SquareLattice(spec.length)
    return tuple(_canonical_bonds(bonds) for bonds in lattice.four_matchings())


def build_finite_xxz_ledger(
    spec: XXZCompileSpec,
    *,
    max_degree: int = 5,
    progress: Callable[[int, int], None] | None = None,
    max_records: int | None = None,
) -> XXZLedger:
    """Build the exact direct-theorem ledger or a deterministic bounded prefix."""

    if not isinstance(spec, XXZCompileSpec):
        raise TypeError("spec must be an XXZCompileSpec")
    if isinstance(max_degree, bool) or not isinstance(max_degree, int):
        raise TypeError("max_degree must be an integer")
    if max_degree < THEOREM_ORDER + 1:
        raise ValueError("max_degree must permit the degree-five theorem words")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable")
    if max_records is not None and (
        isinstance(max_records, bool)
        or not isinstance(max_records, int)
        or max_records < 1
    ):
        raise ValueError("max_records must be a positive integer")

    projected = _projected_raw_record_count()
    limit = projected if max_records is None else min(max_records, projected)
    records: list[RawTheoremRecord] = []
    next_progress = 1
    for record in _iter_raw_theorem_records():
        if len(records) == limit:
            break
        records.append(record)
        while next_progress <= 20 and len(records) * 20 >= next_progress * limit:
            if progress is not None:
                progress(len(records), projected)
            next_progress += 1

    raw_records = tuple(records)
    if len(raw_records) != limit:
        raise ArithmeticError("raw theorem enumeration ended before its projection")
    weights: dict[tuple[int, ...], Fraction] = {}
    for record in raw_records:
        weights[record.block_key] = (
            weights.get(record.block_key, Fraction()) + record.weight_interval.upper
        )

    evaluator = _FiniteXXZSymplecticEvaluator(
        _finite_fragment_bonds(spec),
        spec.delta,
    )
    blocks = tuple(
        TheoremNormBlock(
            block_key=key,
            terms=evaluator.evaluate(key),
            raw_weight_upper=weights[key],
        )
        for key in sorted(weights)
    )
    unsealed = XXZLedger(
        spec=spec,
        theorem_identifier=THEOREM_IDENTIFIER,
        order=THEOREM_ORDER,
        center=THEOREM_CENTER,
        factorial_denominator=THEOREM_FACTORIAL_DENOMINATOR,
        duhamel_convention=THEOREM_DUHAMEL_CONVENTION,
        finite_step_error_formula=FINITE_STEP_ERROR_FORMULA,
        max_degree=max_degree,
        projected_record_count=projected,
        record_limit=max_records,
        complete=limit == projected,
        raw_records=raw_records,
        blocks=blocks,
        ledger_digest="",
    )
    return XXZLedger(
        spec=unsealed.spec,
        theorem_identifier=unsealed.theorem_identifier,
        order=unsealed.order,
        center=unsealed.center,
        factorial_denominator=unsealed.factorial_denominator,
        duhamel_convention=unsealed.duhamel_convention,
        finite_step_error_formula=unsealed.finite_step_error_formula,
        max_degree=unsealed.max_degree,
        projected_record_count=unsealed.projected_record_count,
        record_limit=unsealed.record_limit,
        complete=unsealed.complete,
        raw_records=unsealed.raw_records,
        blocks=unsealed.blocks,
        ledger_digest=_ledger_digest(unsealed),
    )


def verify_finite_xxz_ledger(ledger: XXZLedger) -> None:
    """Re-enumerate raw coverage and every unweighted finite Pauli block."""

    if not isinstance(ledger, XXZLedger):
        raise TypeError("ledger must be an XXZLedger")
    if (
        ledger.theorem_identifier != THEOREM_IDENTIFIER
        or ledger.order != THEOREM_ORDER
        or ledger.center != THEOREM_CENTER
        or ledger.factorial_denominator != THEOREM_FACTORIAL_DENOMINATOR
        or ledger.duhamel_convention != THEOREM_DUHAMEL_CONVENTION
        or ledger.finite_step_error_formula != FINITE_STEP_ERROR_FORMULA
    ):
        raise ValueError("ledger theorem metadata differs from the frozen theorem")
    if ledger.max_degree < THEOREM_ORDER + 1:
        raise ValueError("ledger max_degree does not cover degree five")
    projected = _projected_raw_record_count()
    if ledger.projected_record_count != projected:
        raise ValueError("ledger record projection is incorrect")
    if ledger.record_limit is not None and (
        isinstance(ledger.record_limit, bool)
        or not isinstance(ledger.record_limit, int)
        or ledger.record_limit < 1
    ):
        raise ValueError("ledger record limit is invalid")
    expected_count = (
        projected
        if ledger.record_limit is None
        else min(ledger.record_limit, projected)
    )
    if len(ledger.raw_records) != expected_count:
        raise ValueError("ledger raw record multiplicity is incorrect")
    expected_records: list[RawTheoremRecord] = []
    for record in _iter_raw_theorem_records():
        if len(expected_records) == expected_count:
            break
        expected_records.append(record)
    if ledger.raw_records != tuple(expected_records):
        raise ValueError("ledger raw theorem identities or weights are incorrect")
    if len({record.record_id for record in ledger.raw_records}) != expected_count:
        raise ValueError("ledger raw theorem identities are not unique")
    if any(
        record.fragment_word != record.block_key
        or len(record.block_key) != THEOREM_ORDER + 1
        or record.weight_interval.lower <= 0
        for record in ledger.raw_records
    ):
        raise ValueError("ledger raw theorem word or positive interval is invalid")

    weights: dict[tuple[int, ...], Fraction] = {}
    for record in ledger.raw_records:
        weights[record.block_key] = (
            weights.get(record.block_key, Fraction()) + record.weight_interval.upper
        )
    if tuple(block.block_key for block in ledger.blocks) != tuple(sorted(weights)):
        raise ValueError("ledger norm-block coverage is incorrect")
    evaluator = _FiniteXXZSymplecticEvaluator(
        _finite_fragment_bonds(ledger.spec),
        ledger.spec.delta,
    )
    for block in ledger.blocks:
        if block.raw_weight_upper != weights[block.block_key]:
            raise ValueError("ledger block weight is not the exact raw-weight sum")
        if block.terms != tuple(sorted(block.terms, key=lambda term: term.mask)):
            raise ValueError("ledger Pauli masks are not canonically sorted")
        if any(term.is_zero for term in block.terms):
            raise ValueError("ledger contains a zero Pauli coefficient")
        if block.terms != evaluator.evaluate(block.block_key):
            raise ValueError("ledger C_k map is not the unweighted finite commutator")
    if ledger.complete != (expected_count == projected):
        raise ValueError("ledger completeness flag is incorrect")
    if ledger.ledger_digest != _ledger_digest(ledger):
        raise ValueError("ledger digest is incorrect")


@lru_cache(maxsize=8192)
def sqrt_fraction_interval(
    value: Fraction,
    *,
    decimal_places: int = 30,
) -> RationalInterval:
    """Return a deterministic outward rational enclosure of ``sqrt(value)``."""

    if not isinstance(value, Fraction):
        raise TypeError("square-root radicand must be an exact Fraction")
    if value < 0:
        raise ValueError("cannot bound the square root of a negative value")
    if isinstance(decimal_places, bool) or not isinstance(decimal_places, int):
        raise TypeError("decimal_places must be an integer")
    if decimal_places < 0:
        raise ValueError("decimal_places must be nonnegative")
    if value == 0:
        return RationalInterval.point(0)
    scale = 10**decimal_places
    quotient = value.numerator * scale * scale // value.denominator
    lower_integer = isqrt(quotient)
    lower = Fraction(lower_integer, scale)
    upper_integer = lower_integer
    if lower * lower < value:
        upper_integer += 1
    upper = Fraction(upper_integer, scale)
    if lower * lower > value or upper * upper < value:
        raise ArithmeticError("outward square-root interval failed")
    return RationalInterval(lower, upper)


@dataclass(frozen=True, slots=True)
class AnticommutingGroupRecord:
    terms: tuple[SymplecticCoefficient, ...]
    squared_norm: Fraction
    norm_interval: RationalInterval


@dataclass(frozen=True, slots=True)
class TheoremBlockGroupRecord:
    theorem_identifier: str
    order: int
    center: int
    factorial_denominator: int
    duhamel_convention: str
    finite_step_error_formula: str
    delta: Fraction
    ledger_digest: str
    grouping_algorithm: str
    block_key: tuple[int, ...]
    groups: tuple[AnticommutingGroupRecord, ...]


def _group_record(
    terms: Sequence[SymplecticCoefficient],
) -> AnticommutingGroupRecord:
    squared = sum((term.real * term.real for term in terms), Fraction())
    return AnticommutingGroupRecord(
        terms=tuple(terms),
        squared_norm=squared,
        norm_interval=sqrt_fraction_interval(squared),
    )


def discover_anticommuting_groups(
    terms: Sequence[SymplecticCoefficient],
) -> tuple[AnticommutingGroupRecord, ...]:
    """Deterministically match each term with its earliest anticommuting partner.

    Terms are ordered by decreasing exact coefficient magnitude and then by
    their symplectic mask.  Two bitsets per occupied site coordinate encode,
    for every still-unmatched candidate, the exact symplectic linear
    functional (32 bitsets on the frozen 16-site pilot).  The least
    significant candidate bit is therefore the earliest available
    anticommuting partner in the frozen order.  Every output group is a
    singleton or a pair.
    """

    if any(not isinstance(term, SymplecticCoefficient) for term in terms):
        raise TypeError("group discovery terms must be SymplecticCoefficient values")
    if any(term.imag for term in terms):
        raise ValueError("anticommuting norm groups require Hermitian real coefficients")
    if len({term.mask for term in terms}) != len(terms):
        raise ValueError("coefficient map contains duplicate Pauli masks")
    ordered = sorted(terms, key=lambda term: (-abs(term.real), term.mask))
    if not ordered:
        return ()

    coordinate_width = max(
        (term.x_mask | term.z_mask).bit_length() for term in ordered
    )
    coordinate_bitsets = [0] * (2 * coordinate_width)
    for index, term in enumerate(ordered):
        candidate_bit = 1 << index
        z_mask = term.z_mask
        while z_mask:
            low_bit = z_mask & -z_mask
            coordinate_bitsets[low_bit.bit_length() - 1] |= candidate_bit
            z_mask -= low_bit
        x_mask = term.x_mask
        while x_mask:
            low_bit = x_mask & -x_mask
            coordinate_bitsets[
                coordinate_width + low_bit.bit_length() - 1
            ] |= candidate_bit
            x_mask -= low_bit

    unmatched = (1 << len(ordered)) - 1
    groups: list[AnticommutingGroupRecord] = []
    for index, term in enumerate(ordered):
        term_bit = 1 << index
        if not unmatched & term_bit:
            continue
        unmatched ^= term_bit

        anticommuting = 0
        functional = term.x_mask | (term.z_mask << coordinate_width)
        while functional:
            low_bit = functional & -functional
            anticommuting ^= coordinate_bitsets[low_bit.bit_length() - 1]
            functional -= low_bit
        candidates = anticommuting & unmatched
        if candidates:
            partner_bit = candidates & -candidates
            partner_index = partner_bit.bit_length() - 1
            unmatched ^= partner_bit
            groups.append(_group_record((term, ordered[partner_index])))
        else:
            groups.append(_group_record((term,)))
    if unmatched:
        raise ArithmeticError("pair-only group discovery left unmatched terms")
    return tuple(groups)


def verify_anticommuting_groups(
    terms: Sequence[SymplecticCoefficient],
    groups: Sequence[AnticommutingGroupRecord],
) -> Fraction:
    """Verify exact coverage, coefficients, commutation, and outward norms."""

    expected = {term.mask: term for term in terms}
    if len(expected) != len(terms):
        raise ValueError("coefficient map contains duplicate Pauli masks")
    if any(term.imag for term in terms):
        raise ValueError("anticommuting norm groups require real coefficients")
    if any(not isinstance(group, AnticommutingGroupRecord) for group in groups):
        raise TypeError("groups must contain AnticommutingGroupRecord values")
    flattened = tuple(term for group in groups for term in group.terms)
    masks = tuple(term.mask for term in flattened)
    if len(masks) != len(set(masks)):
        raise ValueError("group coverage contains duplicate Pauli terms")
    if set(masks) != set(expected):
        raise ValueError("group coverage differs from the coefficient map")
    if any(term != expected[term.mask] for term in flattened):
        raise ValueError("group term coefficient differs from the unweighted block")

    total = Fraction()
    for group in groups:
        if not group.terms:
            raise ValueError("anticommuting groups must be nonempty")
        if len(group.terms) > 2:
            raise ValueError("pair-only groups must contain a singleton or pair")
        for index, left in enumerate(group.terms):
            for right in group.terms[index + 1 :]:
                if not _symplectic_anticommutes(left.mask, right.mask):
                    raise ValueError("group members do not anticommute")
        squared = sum(
            (term.real * term.real for term in group.terms),
            Fraction(),
        )
        interval = sqrt_fraction_interval(squared)
        if group.squared_norm != squared or group.norm_interval != interval:
            raise ValueError("submitted group norm interval is incorrect")
        total += interval.upper
    return total


def discover_xxz_groups(
    ledger: XXZLedger,
) -> tuple[TheoremBlockGroupRecord, ...]:
    """Discover deterministic groups independently inside every theorem block."""

    verify_finite_xxz_ledger(ledger)
    return tuple(
        TheoremBlockGroupRecord(
            theorem_identifier=ledger.theorem_identifier,
            order=ledger.order,
            center=ledger.center,
            factorial_denominator=ledger.factorial_denominator,
            duhamel_convention=ledger.duhamel_convention,
            finite_step_error_formula=ledger.finite_step_error_formula,
            delta=ledger.spec.delta,
            ledger_digest=ledger.ledger_digest,
            grouping_algorithm=GROUPING_ALGORITHM_IDENTIFIER,
            block_key=block.block_key,
            groups=discover_anticommuting_groups(block.terms),
        )
        for block in ledger.blocks
    )


def _raw_upper_weights(ledger: XXZLedger) -> dict[tuple[int, ...], Fraction]:
    weights: dict[tuple[int, ...], Fraction] = {}
    for record in ledger.raw_records:
        if record.weight_interval.upper < 0:
            raise ValueError("theorem weights must be nonnegative")
        weights[record.block_key] = (
            weights.get(record.block_key, Fraction()) + record.weight_interval.upper
        )
    return weights


def verify_xxz_groups(
    ledger: XXZLedger,
    groups: tuple[TheoremBlockGroupRecord, ...],
) -> Fraction:
    """Return ``sum_k w_k U_k/5!`` after replaying every block witness."""

    verify_finite_xxz_ledger(ledger)
    if not isinstance(groups, tuple) or any(
        not isinstance(record, TheoremBlockGroupRecord) for record in groups
    ):
        raise TypeError("theorem group witness must be a tuple of block records")
    by_key: dict[tuple[int, ...], TheoremBlockGroupRecord] = {}
    for record in groups:
        if record.block_key in by_key:
            raise ValueError("theorem group witness contains duplicate blocks")
        by_key[record.block_key] = record
    blocks = {block.block_key: block for block in ledger.blocks}
    if set(by_key) != set(blocks):
        raise ValueError("theorem group witness block coverage is incorrect")
    weights = _raw_upper_weights(ledger)
    total = Fraction()
    for key in sorted(blocks):
        record = by_key[key]
        if (
            record.theorem_identifier != ledger.theorem_identifier
            or record.order != ledger.order
            or record.center != ledger.center
            or record.factorial_denominator != ledger.factorial_denominator
            or record.duhamel_convention != ledger.duhamel_convention
            or record.finite_step_error_formula != ledger.finite_step_error_formula
        ):
            raise ValueError("group witness theorem metadata is incorrect")
        if record.delta != ledger.spec.delta:
            raise ValueError("group witness Delta does not match the ledger")
        if record.ledger_digest != ledger.ledger_digest:
            raise ValueError("group witness ledger digest does not match")
        if record.grouping_algorithm != GROUPING_ALGORITHM_IDENTIFIER:
            raise ValueError("group witness algorithm identifier is incorrect")
        block_bound = verify_anticommuting_groups(
            blocks[key].terms,
            record.groups,
        )
        expected_groups = discover_anticommuting_groups(blocks[key].terms)
        if record.groups != expected_groups:
            raise ValueError("group witness differs from deterministic pair discovery")
        total += weights[key] * block_bound
    return total / ledger.factorial_denominator


def xxz_triangle_baseline(ledger: XXZLedger) -> Fraction:
    """Return the same-Delta blockwise Pauli-l1 theorem constant."""

    verify_finite_xxz_ledger(ledger)
    weights = _raw_upper_weights(ledger)
    total = Fraction()
    for block in ledger.blocks:
        pauli_l1 = sum((abs(term.real) for term in block.terms), Fraction())
        if any(term.imag for term in block.terms):
            raise ValueError("degree-five theorem block is not Hermitian")
        total += weights[block.block_key] * pauli_l1
    return total / ledger.factorial_denominator


XXZ_METHOD = "direct_finite_high_order_theorem_grouped_norm"


@dataclass(frozen=True, slots=True)
class XXZCertificate:
    """Exact finite-step closure for one complete same-Delta theorem ledger."""

    status: XXZStatus
    method: str
    spec: XXZCompileSpec
    ledger_digest: str
    raw_record_count: int
    theorem_identifier: str
    order: int
    center: int
    factorial_denominator: int
    duhamel_convention: str
    finite_step_error_formula: str
    grouping_algorithm: str
    grouped_constant: Fraction
    triangle_constant: Fraction
    candidate_steps: int
    candidate_error: Fraction
    candidate_previous_error: Fraction | None
    baseline_steps: int
    baseline_error: Fraction
    baseline_previous_error: Fraction | None
    candidate_resources: int
    baseline_resources: int
    bond_growth: Fraction
    cell_base: Fraction


def _minimal_fourth_order_steps(constant: Fraction, tolerance: Fraction) -> int:
    """Return the least positive ``r`` satisfying ``constant/r**4 <= tolerance``."""

    if not isinstance(constant, Fraction) or not isinstance(tolerance, Fraction):
        raise TypeError("finite-step constants and tolerance must be exact Fractions")
    if constant < 0:
        raise ValueError("finite-step constant must be nonnegative")
    if tolerance <= 0:
        raise ValueError("finite-step tolerance must be positive")
    if constant <= tolerance:
        return 1
    lower = 1
    upper = 2
    while constant > tolerance * upper**THEOREM_ORDER:
        lower = upper
        upper *= 2
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if constant <= tolerance * middle**THEOREM_ORDER:
            upper = middle
        else:
            lower = middle
    return upper


def _fourth_order_error(constant: Fraction, steps: int) -> Fraction:
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
        raise ValueError("finite-step count must be a positive integer")
    return constant / steps**THEOREM_ORDER


def _close_xxz_certificate(
    ledger: XXZLedger,
    *,
    grouped_constant: Fraction,
    triangle_constant: Fraction,
) -> XXZCertificate:
    """Close exact arithmetic after a caller has verified a complete ledger.

    This helper deliberately does not establish ledger authenticity.  The
    production compiler and :func:`verify_xxz_certificate` do that before
    accepting the returned finite-step fields.
    """

    if not isinstance(ledger, XXZLedger):
        raise TypeError("ledger must be an XXZLedger")
    if not ledger.complete:
        raise ValueError("finite-step closure requires a complete ledger")
    if not isinstance(grouped_constant, Fraction) or not isinstance(
        triangle_constant,
        Fraction,
    ):
        raise TypeError("theorem constants must be exact Fractions")
    if grouped_constant < 0 or triangle_constant < 0:
        raise ValueError("theorem constants must be nonnegative")

    tolerance = ledger.spec.tolerance
    candidate_steps = _minimal_fourth_order_steps(grouped_constant, tolerance)
    baseline_steps = _minimal_fourth_order_steps(triangle_constant, tolerance)
    schedule = xxz_suzuki_schedule()
    delta_absolute = abs(ledger.spec.delta)
    certificate = XXZCertificate(
        status="certified",
        method=XXZ_METHOD,
        spec=ledger.spec,
        ledger_digest=ledger.ledger_digest,
        raw_record_count=len(ledger.raw_records),
        theorem_identifier=ledger.theorem_identifier,
        order=ledger.order,
        center=ledger.center,
        factorial_denominator=ledger.factorial_denominator,
        duhamel_convention=ledger.duhamel_convention,
        finite_step_error_formula=ledger.finite_step_error_formula,
        grouping_algorithm=GROUPING_ALGORITHM_IDENTIFIER,
        grouped_constant=grouped_constant,
        triangle_constant=triangle_constant,
        candidate_steps=candidate_steps,
        candidate_error=_fourth_order_error(grouped_constant, candidate_steps),
        candidate_previous_error=(
            None
            if candidate_steps == 1
            else _fourth_order_error(grouped_constant, candidate_steps - 1)
        ),
        baseline_steps=baseline_steps,
        baseline_error=_fourth_order_error(triangle_constant, baseline_steps),
        baseline_previous_error=(
            None
            if baseline_steps == 1
            else _fourth_order_error(triangle_constant, baseline_steps - 1)
        ),
        candidate_resources=replay_schedule_resources(schedule, candidate_steps),
        baseline_resources=replay_schedule_resources(schedule, baseline_steps),
        bond_growth=max(Fraction(1), (1 + delta_absolute) / 2),
        cell_base=(2 + delta_absolute) / 2,
    )
    verify_xxz_finite_step_fields(certificate)
    return certificate


def verify_xxz_finite_step_fields(certificate: XXZCertificate) -> None:
    """Verify exact minimal-step, adjacent-step, resource, and growth fields."""

    if not isinstance(certificate, XXZCertificate):
        raise TypeError("certificate must be an XXZCertificate")
    if certificate.status != "certified" or certificate.method != XXZ_METHOD:
        raise ValueError("XXZ certificate status or method is incorrect")
    if certificate.grouping_algorithm != GROUPING_ALGORITHM_IDENTIFIER:
        raise ValueError("XXZ certificate grouping algorithm is incorrect")
    if (
        certificate.theorem_identifier != THEOREM_IDENTIFIER
        or certificate.order != THEOREM_ORDER
        or certificate.center != THEOREM_CENTER
        or certificate.factorial_denominator != THEOREM_FACTORIAL_DENOMINATOR
        or certificate.duhamel_convention != THEOREM_DUHAMEL_CONVENTION
        or certificate.finite_step_error_formula != FINITE_STEP_ERROR_FORMULA
    ):
        raise ValueError("XXZ certificate theorem metadata is incorrect")
    for name, constant in (
        ("grouped", certificate.grouped_constant),
        ("triangle", certificate.triangle_constant),
    ):
        if not isinstance(constant, Fraction) or constant < 0:
            raise ValueError(f"{name} theorem constant must be an exact nonnegative Fraction")

    tolerance = certificate.spec.tolerance

    def verify_row(
        name: str,
        constant: Fraction,
        steps: int,
        accepted_error: Fraction,
        previous_error: Fraction | None,
        resources: int,
    ) -> None:
        expected_steps = _minimal_fourth_order_steps(constant, tolerance)
        if steps != expected_steps:
            raise ValueError(f"{name} step is not the exact minimal accepted step")
        expected_error = _fourth_order_error(constant, steps)
        if accepted_error != expected_error or accepted_error > tolerance:
            raise ValueError(f"{name} accepted finite-step error is incorrect")
        expected_previous = (
            None
            if steps == 1
            else _fourth_order_error(constant, steps - 1)
        )
        if previous_error != expected_previous:
            raise ValueError(f"{name} previous-step error is incorrect")
        if previous_error is not None and previous_error <= tolerance:
            raise ValueError(f"{name} previous step does not prove minimality")
        expected_resources = replay_schedule_resources(xxz_suzuki_schedule(), steps)
        if resources != expected_resources or resources != 30 * steps + 1:
            raise ValueError(f"{name} merged-group resource count is incorrect")

    verify_row(
        "candidate",
        certificate.grouped_constant,
        certificate.candidate_steps,
        certificate.candidate_error,
        certificate.candidate_previous_error,
        certificate.candidate_resources,
    )
    verify_row(
        "baseline",
        certificate.triangle_constant,
        certificate.baseline_steps,
        certificate.baseline_error,
        certificate.baseline_previous_error,
        certificate.baseline_resources,
    )
    delta_absolute = abs(certificate.spec.delta)
    if certificate.bond_growth != max(Fraction(1), (1 + delta_absolute) / 2):
        raise ValueError("XXZ bond-growth metadata is incorrect")
    if certificate.cell_base != (2 + delta_absolute) / 2:
        raise ValueError("XXZ cell-base metadata is incorrect")


def verify_xxz_certificate(
    certificate: XXZCertificate,
    ledger: XXZLedger,
    groups: tuple[TheoremBlockGroupRecord, ...],
) -> None:
    """Recompute a production certificate from its complete ledger and groups."""

    if not isinstance(certificate, XXZCertificate):
        raise TypeError("certificate must be an XXZCertificate")
    if certificate.grouping_algorithm != GROUPING_ALGORITHM_IDENTIFIER:
        raise ValueError("XXZ certificate grouping algorithm is incorrect")
    if not isinstance(ledger, XXZLedger) or not ledger.complete:
        raise ValueError("XXZ certificate verification requires a complete ledger")
    verify_finite_xxz_ledger(ledger)
    grouped_constant = verify_xxz_groups(ledger, groups)
    triangle_constant = xxz_triangle_baseline(ledger)
    expected = _close_xxz_certificate(
        ledger,
        grouped_constant=grouped_constant,
        triangle_constant=triangle_constant,
    )
    if certificate != expected:
        raise ValueError("XXZ certificate differs from complete-ledger replay")


def compile_grouped_xxz(
    spec: XXZCompileSpec,
    *,
    ledger: XXZLedger | None = None,
    groups: tuple[TheoremBlockGroupRecord, ...] | None = None,
) -> XXZCertificate:
    """Compile the direct finite theorem, requiring full raw-ledger coverage."""

    if not isinstance(spec, XXZCompileSpec):
        raise TypeError("spec must be an XXZCompileSpec")
    selected_ledger = ledger if ledger is not None else build_finite_xxz_ledger(spec)
    if selected_ledger.spec != spec:
        raise ValueError("supplied ledger does not match the requested XXZ spec")
    # Fail before any expensive replay or grouping when a profiling prefix is
    # accidentally supplied to the scientific compile path.
    if not selected_ledger.complete:
        raise ValueError("compile_grouped_xxz requires a complete ledger")
    verify_finite_xxz_ledger(selected_ledger)
    selected_groups = (
        discover_xxz_groups(selected_ledger) if groups is None else groups
    )
    grouped_constant = verify_xxz_groups(selected_ledger, selected_groups)
    triangle_constant = xxz_triangle_baseline(selected_ledger)
    certificate = _close_xxz_certificate(
        selected_ledger,
        grouped_constant=grouped_constant,
        triangle_constant=triangle_constant,
    )
    verify_xxz_certificate(certificate, selected_ledger, selected_groups)
    return certificate
