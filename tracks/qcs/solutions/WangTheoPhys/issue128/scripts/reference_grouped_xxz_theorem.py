"""Independent standard-library oracle for the grouped-XXZ raw theorem stream.

This module deliberately has no imports from ``trottercert``, ``scripts``, or
third-party packages.  It reconstructs the algebraic Suzuki schedule and the
rational interval upper weights used by ``published_high_order_triangle_v1``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

STAGE_COUNT = 31
THEOREM_IDENTIFIER = "published_high_order_triangle_v1"
THEOREM_ORDER = 4
THEOREM_CENTER = 20
THEOREM_DECIMAL_DIGITS = 18
PROJECTED_RECORD_COUNT = 61_677


@dataclass(frozen=True, slots=True)
class Cubic:
    """An element of Q(alpha), represented modulo alpha**3 - 4."""

    a0: Fraction
    a1: Fraction
    a2: Fraction

    def __init__(
        self,
        a0: int | Fraction,
        a1: int | Fraction,
        a2: int | Fraction,
    ) -> None:
        object.__setattr__(self, "a0", Fraction(a0))
        object.__setattr__(self, "a1", Fraction(a1))
        object.__setattr__(self, "a2", Fraction(a2))

    @classmethod
    def one(cls) -> Cubic:
        return cls(1, 0, 0)

    @classmethod
    def coerce(cls, value: Cubic | int | Fraction) -> Cubic:
        return value if isinstance(value, cls) else cls(value, 0, 0)

    def __add__(self, other: Cubic | int | Fraction) -> Cubic:
        rhs = self.coerce(other)
        return Cubic(self.a0 + rhs.a0, self.a1 + rhs.a1, self.a2 + rhs.a2)

    __radd__ = __add__

    def __neg__(self) -> Cubic:
        return Cubic(-self.a0, -self.a1, -self.a2)

    def __sub__(self, other: Cubic | int | Fraction) -> Cubic:
        return self + (-self.coerce(other))

    def __rsub__(self, other: Cubic | int | Fraction) -> Cubic:
        return self.coerce(other) - self

    def __mul__(self, other: Cubic | int | Fraction) -> Cubic:
        rhs = self.coerce(other)
        return Cubic(
            self.a0 * rhs.a0
            + 4 * (self.a1 * rhs.a2 + self.a2 * rhs.a1),
            self.a0 * rhs.a1
            + self.a1 * rhs.a0
            + 4 * self.a2 * rhs.a2,
            self.a0 * rhs.a2 + self.a1 * rhs.a1 + self.a2 * rhs.a0,
        )

    __rmul__ = __mul__

    def __truediv__(self, denominator: int | Fraction) -> Cubic:
        divisor = Fraction(denominator)
        if divisor == 0:
            raise ZeroDivisionError("division by zero")
        return Cubic(self.a0 / divisor, self.a1 / divisor, self.a2 / divisor)

    def enclose(self, root: Interval) -> Interval:
        return Interval.point(self.a0) + self.a1 * root + self.a2 * root**2


@dataclass(frozen=True, slots=True)
class Interval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if self.lower > self.upper:
            raise ValueError("interval endpoints are reversed")

    @classmethod
    def point(cls, value: int | Fraction) -> Interval:
        exact = Fraction(value)
        return cls(exact, exact)

    @classmethod
    def coerce(cls, value: Interval | int | Fraction) -> Interval:
        return value if isinstance(value, cls) else cls.point(value)

    def __add__(self, other: Interval | int | Fraction) -> Interval:
        rhs = self.coerce(other)
        return Interval(self.lower + rhs.lower, self.upper + rhs.upper)

    __radd__ = __add__

    def __neg__(self) -> Interval:
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other: Interval | int | Fraction) -> Interval:
        return self + (-self.coerce(other))

    def __rsub__(self, other: Interval | int | Fraction) -> Interval:
        return self.coerce(other) - self

    def __mul__(self, other: Interval | int | Fraction) -> Interval:
        rhs = self.coerce(other)
        products = (
            self.lower * rhs.lower,
            self.lower * rhs.upper,
            self.upper * rhs.lower,
            self.upper * rhs.upper,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self) -> Interval:
        if self.lower <= 0 <= self.upper:
            raise ZeroDivisionError("interval contains zero")
        return Interval(1 / self.upper, 1 / self.lower)

    def __truediv__(self, other: Interval | int | Fraction) -> Interval:
        return self * self.coerce(other).reciprocal()

    def __pow__(self, exponent: int) -> Interval:
        if isinstance(exponent, bool) or not isinstance(exponent, int):
            raise TypeError("interval exponent must be an integer")
        if exponent < 0:
            return self.reciprocal() ** (-exponent)
        result = self.point(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result

    def absolute(self) -> Interval:
        lower = Fraction() if self.lower <= 0 <= self.upper else min(
            abs(self.lower), abs(self.upper)
        )
        return Interval(lower, max(abs(self.lower), abs(self.upper)))


@dataclass(frozen=True, slots=True)
class Stage:
    fragment_index: int
    algebraic_coefficient: Cubic
    coefficient_interval: Interval


@dataclass(frozen=True, slots=True)
class RawTheoremRecord:
    record_id: str
    side: str
    partial_sum_index: int
    adjoint_stage_indices: tuple[int, ...]
    composition: tuple[int, ...]
    base_stage_index: int
    fragment_word: tuple[int, ...]
    weight_upper: Fraction


@dataclass(frozen=True, slots=True)
class RawTheoremSummary:
    record_count: int
    projected_record_count: int
    complete: bool
    stream_sha256: str
    word_weights: tuple[tuple[tuple[int, ...], Fraction], ...]


def _root_four_interval(decimal_digits: int) -> Interval:
    if isinstance(decimal_digits, bool) or not isinstance(decimal_digits, int):
        raise TypeError("decimal_digits must be an integer")
    if decimal_digits < 1:
        raise ValueError("decimal_digits must be positive")
    denominator = 10**decimal_digits
    low = denominator
    high = 4 * denominator
    target = 4 * denominator**3
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**3 <= target:
            low = middle
        else:
            high = middle
    return Interval(Fraction(low, denominator), Fraction(high, denominator))


def _outward_quantize(interval: Interval, denominator: int) -> Interval:
    lower_scaled = interval.lower * denominator
    upper_scaled = interval.upper * denominator
    lower = lower_scaled.numerator // lower_scaled.denominator
    upper = -((-upper_scaled.numerator) // upper_scaled.denominator)
    return Interval(Fraction(lower, denominator), Fraction(upper, denominator))


def _second_order_stages(
    algebraic_scale: Cubic,
    interval_scale: Interval,
) -> list[Stage]:
    algebraic_half = algebraic_scale / 2
    interval_half = interval_scale / 2
    stages = [Stage(index, algebraic_half, interval_half) for index in range(3)]
    stages.append(Stage(3, algebraic_scale, interval_scale))
    stages.extend(
        Stage(index, algebraic_half, interval_half)
        for index in reversed(range(3))
    )
    return stages


def _merge_stages(stages: Sequence[Stage]) -> tuple[Stage, ...]:
    merged: list[Stage] = []
    for stage in stages:
        if merged and merged[-1].fragment_index == stage.fragment_index:
            previous = merged.pop()
            merged.append(
                Stage(
                    stage.fragment_index,
                    previous.algebraic_coefficient + stage.algebraic_coefficient,
                    previous.coefficient_interval + stage.coefficient_interval,
                )
            )
        else:
            merged.append(stage)
    return tuple(merged)


def theorem_stages() -> tuple[Stage, ...]:
    """Return the independently reconstructed one-based theorem schedule."""

    root = _root_four_interval(THEOREM_DECIMAL_DIGITS)
    grid = 10**THEOREM_DECIMAL_DIGITS
    algebraic_u = Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
    if (4 - Cubic(0, 1, 0)) * algebraic_u != Cubic.one():
        raise ArithmeticError("algebraic Suzuki inverse identity failed")
    interval_u = _outward_quantize(Interval.point(1) / (4 - root), grid)
    algebraic_scales = (
        algebraic_u,
        algebraic_u,
        Cubic.one() - 4 * algebraic_u,
        algebraic_u,
        algebraic_u,
    )
    interval_scales = (
        interval_u,
        interval_u,
        1 - 4 * interval_u,
        interval_u,
        interval_u,
    )
    unmerged: list[Stage] = []
    for algebraic_scale, interval_scale in zip(
        algebraic_scales, interval_scales, strict=True
    ):
        unmerged.extend(_second_order_stages(algebraic_scale, interval_scale))
    result = tuple(reversed(_merge_stages(unmerged)))
    if len(result) != STAGE_COUNT:
        raise ArithmeticError("Suzuki schedule does not have 31 stages")
    for stage in result:
        exact_enclosure = stage.algebraic_coefficient.enclose(root)
        interval = stage.coefficient_interval
        if interval.lower > exact_enclosure.lower or interval.upper < exact_enclosure.upper:
            raise ArithmeticError("rational schedule does not enclose Q(alpha) stage")
    return result


def weak_compositions(total: int, length: int) -> Iterator[tuple[int, ...]]:
    if any(isinstance(value, bool) or not isinstance(value, int) for value in (total, length)):
        raise TypeError("composition parameters must be integers")
    if total < 0 or length < 1:
        raise ValueError("composition parameters are outside their domain")
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first,) + tail


def _multinomial(composition: Sequence[int]) -> int:
    result = factorial(sum(composition))
    for value in composition:
        result //= factorial(value)
    return result


def projected_record_count() -> int:
    left = sum(
        (j - 1) * comb(THEOREM_CENTER - j + THEOREM_ORDER - 1, THEOREM_ORDER - 1)
        for j in range(2, THEOREM_CENTER + 1)
    )
    right = sum(
        (j - 1) * comb(j - THEOREM_CENTER + THEOREM_ORDER - 2, THEOREM_ORDER - 1)
        for j in range(THEOREM_CENTER + 1, STAGE_COUNT + 1)
    )
    if left + right != PROJECTED_RECORD_COUNT:
        raise ArithmeticError("raw theorem count differs from the frozen projection")
    return left + right


def iter_raw_theorem_records(
    max_records: int | None = None,
) -> Iterator[RawTheoremRecord]:
    """Yield the theorem multiset in its canonical primary enumeration order."""

    if max_records is not None and (
        isinstance(max_records, bool)
        or not isinstance(max_records, int)
        or max_records < 1
    ):
        raise ValueError("max_records must be a positive integer or None")
    stages = theorem_stages()
    emitted = 0

    def records_for(
        side: str,
        j: int,
        indices: tuple[int, ...],
        composition: tuple[int, ...],
    ) -> Iterator[RawTheoremRecord]:
        nonlocal emitted
        outer: list[int] = []
        scalar_upper = Fraction(_multinomial(composition))
        for stage_index, power in zip(indices, composition, strict=True):
            stage = stages[stage_index - 1]
            scalar_upper *= stage.coefficient_interval.absolute().upper**power
            outer.extend([stage.fragment_index] * power)
        weight_cache: dict[Fraction, Fraction] = {}
        for base_index in range(1, j):
            if max_records is not None and emitted >= max_records:
                return
            base_stage = stages[base_index - 1]
            base_upper = base_stage.coefficient_interval.absolute().upper
            weight_upper = weight_cache.get(base_upper)
            if weight_upper is None:
                weight_upper = scalar_upper * base_upper
                weight_cache[base_upper] = weight_upper
            composition_text = ",".join(str(value) for value in composition)
            emitted += 1
            yield RawTheoremRecord(
                record_id=f"{side}:{j}:{composition_text}:{base_index}",
                side=side,
                partial_sum_index=j,
                adjoint_stage_indices=indices,
                composition=composition,
                base_stage_index=base_index,
                fragment_word=tuple(outer) + (base_stage.fragment_index,),
                weight_upper=weight_upper,
            )

    for j in range(2, THEOREM_CENTER + 1):
        indices = tuple(range(THEOREM_CENTER, j - 1, -1))
        for composition in weak_compositions(THEOREM_ORDER, len(indices)):
            if composition[-1]:
                yield from records_for("left", j, indices, composition)
                if max_records is not None and emitted >= max_records:
                    return
    for j in range(THEOREM_CENTER + 1, STAGE_COUNT + 1):
        indices = tuple(range(THEOREM_CENTER + 1, j + 1))
        for composition in weak_compositions(THEOREM_ORDER, len(indices)):
            if composition[-1]:
                yield from records_for("right", j, indices, composition)
                if max_records is not None and emitted >= max_records:
                    return


def fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def record_payload(record: RawTheoremRecord) -> dict[str, object]:
    if not isinstance(record, RawTheoremRecord):
        raise TypeError("record must be a RawTheoremRecord")
    return {
        "record_id": record.record_id,
        "side": record.side,
        "partial_sum_index": record.partial_sum_index,
        "adjoint_stage_indices": list(record.adjoint_stage_indices),
        "composition": list(record.composition),
        "base_stage_index": record.base_stage_index,
        "fragment_word": list(record.fragment_word),
        "weight_upper": fraction_pair(record.weight_upper),
    }


def canonical_json_bytes(payload: object) -> bytes:
    return (
        json.dumps(
            payload,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def _reject_numeric_aliases(value: object, path: str = "payload") -> None:
    if isinstance(value, (bool, float)):
        raise TypeError(f"{path}: bool and float numeric aliases are forbidden")
    if isinstance(value, list):
        for index, item in enumerate(value):
            _reject_numeric_aliases(item, f"{path}[{index}]")
    elif isinstance(value, Mapping):
        for key, item in value.items():
            _reject_numeric_aliases(item, f"{path}.{key}")


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_canonical_json(raw: bytes) -> object:
    """Load canonical JSON while rejecting duplicate keys and numeric aliases."""

    if not isinstance(raw, bytes):
        raise TypeError("canonical JSON input must be bytes")
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=_unique_object,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ValueError(f"invalid JSON constant: {token}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("invalid canonical JSON") from error
    _reject_numeric_aliases(value)
    if raw != canonical_json_bytes(value):
        raise ValueError("JSON is not in canonical byte form")
    return value


def strict_fraction_pair(value: object, field: str) -> Fraction:
    if not isinstance(value, list) or len(value) != 2:
        raise TypeError(f"{field} must be a two-element array")
    numerator, denominator = value
    if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
        raise TypeError(f"{field} entries must be integers")
    if denominator <= 0:
        raise ValueError(f"{field} denominator must be positive")
    result = Fraction(numerator, denominator)
    if value != fraction_pair(result):
        raise ValueError(f"{field} must be a reduced canonical fraction")
    return result


def _strict_int(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{field} must be an integer")
    return value


def _strict_int_tuple(value: object, field: str) -> tuple[int, ...]:
    if not isinstance(value, list):
        raise TypeError(f"{field} must be an array")
    return tuple(_strict_int(item, f"{field}[{index}]") for index, item in enumerate(value))


def decode_record_payload(payload: object) -> RawTheoremRecord:
    """Strictly decode and semantically replay one canonical record payload."""

    fields = {
        "record_id",
        "side",
        "partial_sum_index",
        "adjoint_stage_indices",
        "composition",
        "base_stage_index",
        "fragment_word",
        "weight_upper",
    }
    if not isinstance(payload, Mapping) or set(payload) != fields:
        raise ValueError("raw record field set mismatch")
    _reject_numeric_aliases(payload)
    record_id = payload["record_id"]
    side = payload["side"]
    if not isinstance(record_id, str) or not isinstance(side, str):
        raise TypeError("record identity and side must be strings")
    j = _strict_int(payload["partial_sum_index"], "partial_sum_index")
    base_index = _strict_int(payload["base_stage_index"], "base_stage_index")
    indices = _strict_int_tuple(payload["adjoint_stage_indices"], "adjoint_stage_indices")
    composition = _strict_int_tuple(payload["composition"], "composition")
    word = _strict_int_tuple(payload["fragment_word"], "fragment_word")
    weight_upper = strict_fraction_pair(payload["weight_upper"], "weight_upper")
    if side == "left" and 2 <= j <= THEOREM_CENTER:
        expected_indices = tuple(range(THEOREM_CENTER, j - 1, -1))
    elif side == "right" and THEOREM_CENTER < j <= STAGE_COUNT:
        expected_indices = tuple(range(THEOREM_CENTER + 1, j + 1))
    else:
        raise ValueError("record side or partial-sum index is invalid")
    if indices != expected_indices:
        raise ValueError("record adjoint-stage indices are invalid")
    if (
        len(composition) != len(indices)
        or any(value < 0 for value in composition)
        or sum(composition) != THEOREM_ORDER
        or not composition[-1]
    ):
        raise ValueError("record composition is invalid")
    if not 1 <= base_index < j:
        raise ValueError("record base-stage index is invalid")
    expected_id = f"{side}:{j}:{','.join(map(str, composition))}:{base_index}"
    if record_id != expected_id:
        raise ValueError("record identifier is not canonical")
    stages = theorem_stages()
    expected_word: list[int] = []
    scalar_upper = Fraction(_multinomial(composition))
    for stage_index, power in zip(indices, composition, strict=True):
        stage = stages[stage_index - 1]
        expected_word.extend([stage.fragment_index] * power)
        scalar_upper *= stage.coefficient_interval.absolute().upper**power
    base_stage = stages[base_index - 1]
    expected_word.append(base_stage.fragment_index)
    expected_weight = scalar_upper * base_stage.coefficient_interval.absolute().upper
    if word != tuple(expected_word) or len(word) != THEOREM_ORDER + 1:
        raise ValueError("record actual fragment word is invalid")
    if weight_upper != expected_weight or weight_upper <= 0:
        raise ValueError("record upper weight is invalid")
    return RawTheoremRecord(
        record_id=record_id,
        side=side,
        partial_sum_index=j,
        adjoint_stage_indices=indices,
        composition=composition,
        base_stage_index=base_index,
        fragment_word=word,
        weight_upper=weight_upper,
    )


def summarize_raw_theorem(max_records: int | None = None) -> RawTheoremSummary:
    digest = hashlib.sha256()
    weights: dict[tuple[int, ...], Fraction] = {}
    count = 0
    for record in iter_raw_theorem_records(max_records):
        digest.update(canonical_json_bytes(record_payload(record)))
        weights[record.fragment_word] = (
            weights.get(record.fragment_word, Fraction()) + record.weight_upper
        )
        count += 1
    projected = projected_record_count()
    expected = projected if max_records is None else min(max_records, projected)
    if count != expected:
        raise ArithmeticError("raw theorem stream ended at the wrong count")
    return RawTheoremSummary(
        record_count=count,
        projected_record_count=projected,
        complete=count == projected,
        stream_sha256=digest.hexdigest(),
        word_weights=tuple(sorted(weights.items())),
    )


def summary_payload(summary: RawTheoremSummary) -> dict[str, object]:
    return {
        "schema_version": 1,
        "theorem_identifier": THEOREM_IDENTIFIER,
        "order": THEOREM_ORDER,
        "center": THEOREM_CENTER,
        "stage_count": STAGE_COUNT,
        "record_count": summary.record_count,
        "projected_record_count": summary.projected_record_count,
        "coverage_status": "complete" if summary.complete else "bounded_prefix",
        "stream_sha256": summary.stream_sha256,
        "word_weights": [
            {"fragment_word": list(word), "weight_upper": fraction_pair(weight)}
            for word, weight in summary.word_weights
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-records", type=int)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    summary = summarize_raw_theorem(arguments.max_records)
    encoded = canonical_json_bytes(summary_payload(summary))
    if arguments.output is None:
        print(encoded.decode("ascii"), end="")
    else:
        if arguments.output.exists():
            parser.error(f"refusing to overwrite existing path: {arguments.output}")
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
