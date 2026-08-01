"""Independent exact verifier for the two-fragment Suzuki PF4 BCH mapping.

This file is deliberately self contained.  It uses a recursive BCH equation in
the Hall--Lyndon free-Lie basis and expands into associative words only for the
final cyclic-trace comparison.  It does not import the primary implementation.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from itertools import product
from math import factorial
from pathlib import Path

MAX_DEGREE = 5
SCHEMA_VERSION = 1
ARTIFACT_KIND = "pf4_two_fragment_exact_cyclic_mapping"
IDENTITY_STATUS = "proved_exact_suzuki_pf4_free_trace_identity"
ISSUE_ROOT = Path(__file__).resolve().parents[1]
WRAPPER_SCHEMA_VERSION = 2
WRAPPER_KIND = "pf4_bch_mapping_source_closed_artifact"
REFERENCE_ALGORITHM = "recursive_bch_hall_lyndon_v1"
SOURCE_ROOTS = (
    "src/trottercert/pf4_bch_mapping.py",
    "src/trottercert/cubic_field.py",
    "src/trottercert/intervals.py",
    "scripts/derive_pf4_bch_mapping.py",
    "scripts/reference_pf4_bch_mapping.py",
    "pyproject.toml",
    "requirements-reproducibility.txt",
)
WRAPPER_FIELDS = {
    "schema_version",
    "kind",
    "mapping",
    "implementation_sources",
    "reference_algorithm",
    "payload_sha256",
}

Word = tuple[int, ...]
RationalWordPolynomial = dict[Word, Fraction]


@dataclass(frozen=True, slots=True)
class Cubic:
    """An exact element of Q(alpha), alpha**3 = 4."""

    a0: Fraction
    a1: Fraction
    a2: Fraction

    def __init__(
        self,
        a0: int | Fraction = 0,
        a1: int | Fraction = 0,
        a2: int | Fraction = 0,
    ) -> None:
        object.__setattr__(self, "a0", Fraction(a0))
        object.__setattr__(self, "a1", Fraction(a1))
        object.__setattr__(self, "a2", Fraction(a2))

    @classmethod
    def coerce(cls, value: Cubic | int | Fraction) -> Cubic:
        return value if isinstance(value, cls) else cls(Fraction(value))

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
            self.a0 * rhs.a0 + 4 * (self.a1 * rhs.a2 + self.a2 * rhs.a1),
            self.a0 * rhs.a1 + self.a1 * rhs.a0 + 4 * self.a2 * rhs.a2,
            self.a0 * rhs.a2 + self.a1 * rhs.a1 + self.a2 * rhs.a0,
        )

    __rmul__ = __mul__

    def __truediv__(self, denominator: int | Fraction) -> Cubic:
        divisor = Fraction(denominator)
        if not divisor:
            raise ZeroDivisionError("division by zero")
        return Cubic(self.a0 / divisor, self.a1 / divisor, self.a2 / divisor)

    def is_zero(self) -> bool:
        return not (self.a0 or self.a1 or self.a2)

    def coordinates(self) -> tuple[Fraction, Fraction, Fraction]:
        return self.a0, self.a1, self.a2


ZERO = Cubic()
ONE = Cubic(1)


def _add_term(target: dict[Word, Cubic], key: Word, value: Cubic) -> None:
    if value.is_zero():
        return
    updated = target.get(key, ZERO) + value
    if updated.is_zero():
        target.pop(key, None)
    else:
        target[key] = updated


# A Lie polynomial is stored only in the standard Hall--Lyndon basis.  Its
# keys are Lyndon words, whose standard bracketings identify the basis vectors.
LiePolynomial = dict[Word, Cubic]
AuxiliaryPolynomial = dict[int, LiePolynomial]


def _is_lyndon(word: Word) -> bool:
    return bool(word) and all(word < word[offset:] for offset in range(1, len(word)))


@cache
def _lyndon_words(degree: int) -> tuple[Word, ...]:
    return tuple(word for word in product((0, 1), repeat=degree) if _is_lyndon(word))


@cache
def _standard_factorization(word: Word) -> tuple[Word, Word]:
    if len(word) < 2 or not _is_lyndon(word):
        raise ValueError("standard factorization requires a nontrivial Lyndon word")
    for offset in range(1, len(word)):
        suffix = word[offset:]
        if _is_lyndon(suffix):
            prefix = word[:offset]
            if _is_lyndon(prefix):
                return prefix, suffix
    raise ArithmeticError("Lyndon standard factorization was not found")


def _rational_add_term(
    target: RationalWordPolynomial,
    word: Word,
    coefficient: Fraction,
) -> None:
    updated = target.get(word, Fraction()) + coefficient
    if updated:
        target[word] = updated
    else:
        target.pop(word, None)


def _rational_multiply(
    left: Mapping[Word, Fraction],
    right: Mapping[Word, Fraction],
) -> RationalWordPolynomial:
    result: RationalWordPolynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            _rational_add_term(
                result,
                left_word + right_word,
                left_coefficient * right_coefficient,
            )
    return result


def _rational_commutator(
    left: Mapping[Word, Fraction],
    right: Mapping[Word, Fraction],
) -> RationalWordPolynomial:
    result = _rational_multiply(left, right)
    for word, coefficient in _rational_multiply(right, left).items():
        _rational_add_term(result, word, -coefficient)
    return result


@cache
def _hall_word_expansion(word: Word) -> tuple[tuple[Word, Fraction], ...]:
    if len(word) == 1:
        return ((word, Fraction(1)),)
    prefix, suffix = _standard_factorization(word)
    expansion = _rational_commutator(
        dict(_hall_word_expansion(prefix)),
        dict(_hall_word_expansion(suffix)),
    )
    if expansion.get(word) != 1:
        raise ArithmeticError("nonunit leading coefficient in Hall--Lyndon basis")
    return tuple(sorted(expansion.items()))


def _decompose_lie_words(polynomial: Mapping[Word, Fraction]) -> dict[Word, Fraction]:
    if not polynomial:
        return {}
    degrees = {len(word) for word in polynomial}
    if len(degrees) != 1:
        raise ValueError("Hall decomposition requires a homogeneous polynomial")
    degree = next(iter(degrees))
    residual = dict(polynomial)
    result: dict[Word, Fraction] = {}
    for basis_word in _lyndon_words(degree):
        coefficient = residual.get(basis_word, Fraction())
        if not coefficient:
            continue
        result[basis_word] = coefficient
        for word, value in _hall_word_expansion(basis_word):
            _rational_add_term(residual, word, -coefficient * value)
    if residual:
        raise ArithmeticError("commutator did not decompose in the Hall--Lyndon basis")
    return result


@cache
def _basis_bracket(left: Word, right: Word) -> tuple[tuple[Word, Fraction], ...]:
    degree = len(left) + len(right)
    if degree > MAX_DEGREE:
        return ()
    words = _rational_commutator(
        dict(_hall_word_expansion(left)),
        dict(_hall_word_expansion(right)),
    )
    return tuple(sorted(_decompose_lie_words(words).items()))


def _lie_add(left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]) -> LiePolynomial:
    result = dict(left)
    for word, coefficient in right.items():
        _add_term(result, word, coefficient)
    return result


def _lie_scale(
    polynomial: Mapping[Word, Cubic],
    scalar: Cubic | int | Fraction,
) -> LiePolynomial:
    factor = Cubic.coerce(scalar)
    return {
        word: value
        for word, coefficient in polynomial.items()
        if not (value := coefficient * factor).is_zero()
    }


def _lie_bracket(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> LiePolynomial:
    result: LiePolynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            if len(left_word) + len(right_word) > MAX_DEGREE:
                continue
            scalar = left_coefficient * right_coefficient
            for basis_word, structure in _basis_bracket(left_word, right_word):
                _add_term(result, basis_word, scalar * structure)
    return result


def _aux_add(
    left: Mapping[int, LiePolynomial],
    right: Mapping[int, LiePolynomial],
) -> AuxiliaryPolynomial:
    result = {power: dict(value) for power, value in left.items()}
    for power, value in right.items():
        updated = _lie_add(result.get(power, {}), value)
        if updated:
            result[power] = updated
        else:
            result.pop(power, None)
    return result


def _aux_scale(
    polynomial: Mapping[int, LiePolynomial], scalar: int | Fraction
) -> AuxiliaryPolynomial:
    return {
        power: scaled
        for power, value in polynomial.items()
        if (scaled := _lie_scale(value, scalar))
    }


def _aux_bracket(
    left: Mapping[int, LiePolynomial],
    right: Mapping[int, LiePolynomial],
    *,
    max_power: int,
) -> AuxiliaryPolynomial:
    result: AuxiliaryPolynomial = {}
    for left_power, left_value in left.items():
        for right_power, right_value in right.items():
            power = left_power + right_power
            if power > max_power:
                continue
            bracket = _lie_bracket(left_value, right_value)
            if bracket:
                result[power] = _lie_add(result.get(power, {}), bracket)
    return result


def _exp_ad(left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]) -> LiePolynomial:
    result = dict(right)
    term = dict(right)
    for exponent in range(1, MAX_DEGREE):
        term = _lie_bracket(left, term)
        if not term:
            break
        result = _lie_add(result, _lie_scale(term, Fraction(1, factorial(exponent))))
    return result


def _bch(left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]) -> LiePolynomial:
    """Return log(exp(left) exp(right)) through degree five.

    For F(s)=exp(left)exp(s*right)=exp(Z(s)), its right logarithmic
    derivative is W=exp(ad_left)right and

        Z' = W - [Z,W]/2 + [Z,[Z,W]]/12
             - [Z,[Z,[Z,[Z,W]]]]/720.

    Solving the auxiliary-s coefficients recursively keeps the calculation in
    the Hall--Lyndon Lie basis throughout.
    """

    if not left:
        return dict(right)
    if not right:
        return dict(left)
    w = _exp_ad(left, right)
    w_aux = {0: w}
    z_aux: AuxiliaryPolynomial = {0: dict(left)}
    for power in range(1, MAX_DEGREE + 1):
        max_aux = power - 1
        ad1 = _aux_bracket(z_aux, w_aux, max_power=max_aux)
        ad2 = _aux_bracket(z_aux, ad1, max_power=max_aux)
        ad3 = _aux_bracket(z_aux, ad2, max_power=max_aux)
        ad4 = _aux_bracket(z_aux, ad3, max_power=max_aux)
        rhs = _aux_add(
            w_aux,
            _aux_add(
                _aux_scale(ad1, Fraction(-1, 2)),
                _aux_add(
                    _aux_scale(ad2, Fraction(1, 12)),
                    _aux_scale(ad4, Fraction(-1, 720)),
                ),
            ),
        )
        coefficient = _lie_scale(rhs.get(max_aux, {}), Fraction(1, power))
        if coefficient:
            z_aux[power] = coefficient
    result: LiePolynomial = {}
    for coefficient in z_aux.values():
        result = _lie_add(result, coefficient)
    return result


def _literal_stages() -> tuple[tuple[int, Cubic], ...]:
    u = Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
    v = ONE - 4 * u
    return (
        (0, u / 2),
        (1, u),
        (0, u),
        (1, u),
        (0, (u + v) / 2),
        (1, v),
        (0, (u + v) / 2),
        (1, u),
        (0, u),
        (1, u),
        (0, u / 2),
    )


def derive_log_hall() -> LiePolynomial:
    logarithm: LiePolynomial = {}
    for fragment, coefficient in _literal_stages():
        logarithm = _bch(logarithm, {(fragment,): coefficient})
    expected_linear = {(0,): ONE, (1,): ONE}
    if {
        word: value for word, value in logarithm.items() if len(word) == 1
    } != expected_linear:
        raise ArithmeticError("PF4 linear order condition failed")
    for degree in (2, 3, 4):
        if any(len(word) == degree for word in logarithm):
            raise ArithmeticError(f"PF4 degree-{degree} order condition failed")
    if not any(len(word) == 5 for word in logarithm):
        raise ArithmeticError("PF4 degree-five logarithm vanished unexpectedly")
    return logarithm


def _expand_lie(polynomial: Mapping[Word, Cubic]) -> dict[Word, Cubic]:
    result: dict[Word, Cubic] = {}
    for basis_word, coefficient in polynomial.items():
        for word, rational in _hall_word_expansion(basis_word):
            _add_term(result, word, coefficient * rational)
    return result


def _word_multiply(
    left: Mapping[Word, Cubic], right: Mapping[Word, Cubic]
) -> dict[Word, Cubic]:
    result: dict[Word, Cubic] = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            _add_term(
                result,
                left_word + right_word,
                left_coefficient * right_coefficient,
            )
    return result


def _cyclic_classes(polynomial: Mapping[Word, Cubic]) -> dict[Word, Cubic]:
    result: dict[Word, Cubic] = {}
    for word, coefficient in polynomial.items():
        representative = min(
            word[offset:] + word[:offset] for offset in range(len(word))
        )
        _add_term(result, representative, coefficient)
    return result


def _solve_overdetermined(
    matrix: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    if len(matrix) != len(rhs) or not matrix:
        raise ValueError("invalid linear system dimensions")
    columns = len(matrix[0])
    rows = [list(row) + [value] for row, value in zip(matrix, rhs)]
    if any(len(row) != columns + 1 for row in rows):
        raise ValueError("ragged linear system")
    pivot_row = 0
    pivots: list[int] = []
    for column in range(columns):
        pivot = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        divisor = rows[pivot_row][column]
        rows[pivot_row] = [value / divisor for value in rows[pivot_row]]
        for index, row in enumerate(rows):
            if index == pivot_row or not row[column]:
                continue
            factor = row[column]
            rows[index] = [
                value - factor * pivot_value
                for value, pivot_value in zip(row, rows[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
    if any(not any(row[:columns]) and row[columns] for row in rows):
        raise ArithmeticError("cyclic coefficient system is inconsistent")
    if len(pivots) != columns:
        raise ArithmeticError("cyclic coefficient solution is not unique")
    solution = [Fraction() for _ in range(columns)]
    for row_index, column in enumerate(pivots):
        solution[column] = rows[row_index][columns]
    return tuple(solution)


@dataclass(frozen=True, slots=True)
class ProjectiveComparison:
    hypothesis: tuple[Fraction, Fraction, Fraction]
    status: str
    scale: Cubic
    residual: tuple[Cubic, Cubic, Cubic]


def _projective_comparison(
    solved: Sequence[Cubic], candidate: Sequence[Fraction]
) -> ProjectiveComparison:
    pivot = next((index for index, value in enumerate(candidate) if value), None)
    if pivot is None:
        raise ValueError("projective candidate is zero")
    scale = solved[pivot] / candidate[pivot]
    residual = tuple(
        value - scale * coefficient for value, coefficient in zip(solved, candidate)
    )
    status = (
        "proportional"
        if all(value.is_zero() for value in residual)
        else "not_proportional"
    )
    return ProjectiveComparison(
        hypothesis=(candidate[0], candidate[1], candidate[2]),
        status=status,
        scale=scale,
        residual=(residual[0], residual[1], residual[2]),
    )


@dataclass(frozen=True, slots=True)
class DerivedMapping:
    logarithm: LiePolynomial
    target_cyclic: dict[Word, Cubic]
    c2_cyclic: dict[Word, Cubic]
    cd_cyclic: dict[Word, Cubic]
    d2_cyclic: dict[Word, Cubic]
    reconstructed_cyclic: dict[Word, Cubic]
    residual: dict[Word, Cubic]
    coefficients: tuple[Cubic, Cubic, Cubic]
    candidate_comparison: ProjectiveComparison
    legacy_comparison: ProjectiveComparison


def derive_mapping() -> DerivedMapping:
    logarithm = derive_log_hall()
    a = {(0,): ONE}
    b = {(1,): ONE}
    h = _lie_add(a, b)
    c = _lie_bracket(a, _lie_bracket(a, b))
    d = _lie_bracket(b, _lie_bracket(b, a))
    l5 = {word: value for word, value in logarithm.items() if len(word) == 5}
    target = _cyclic_classes(_word_multiply(_expand_lie(h), _expand_lie(l5)))
    basis = (
        _cyclic_classes(_word_multiply(_expand_lie(c), _expand_lie(c))),
        _cyclic_classes(_word_multiply(_expand_lie(c), _expand_lie(d))),
        _cyclic_classes(_word_multiply(_expand_lie(d), _expand_lie(d))),
    )
    words = sorted(set(target).union(*(set(polynomial) for polynomial in basis)))
    matrix = [[polynomial.get(word, ZERO).a0 for polynomial in basis] for word in words]
    coordinate_solutions = []
    for coordinate in range(3):
        rhs = [target.get(word, ZERO).coordinates()[coordinate] for word in words]
        coordinate_solutions.append(_solve_overdetermined(matrix, rhs))
    solved = tuple(
        Cubic(*(coordinate_solutions[coordinate][index] for coordinate in range(3)))
        for index in range(3)
    )
    reconstructed: dict[Word, Cubic] = {}
    for coefficient, polynomial in zip(solved, basis):
        for word, value in polynomial.items():
            _add_term(reconstructed, word, coefficient * value)
    if reconstructed != target:
        raise ArithmeticError("nonzero residual after cyclic coefficient solve")
    candidate_comparison = _projective_comparison(
        solved,
        (Fraction(1), Fraction(-4), Fraction(8, 3)),
    )
    legacy_comparison = _projective_comparison(
        solved,
        (Fraction(1, 2), Fraction(14, 3), Fraction(4, 3)),
    )
    return DerivedMapping(
        logarithm=logarithm,
        target_cyclic=target,
        c2_cyclic=basis[0],
        cd_cyclic=basis[1],
        d2_cyclic=basis[2],
        reconstructed_cyclic=reconstructed,
        residual={},
        coefficients=(solved[0], solved[1], solved[2]),
        candidate_comparison=candidate_comparison,
        legacy_comparison=legacy_comparison,
    )


def _floor_cube_root(value: int) -> int:
    if value < 0:
        raise ValueError("cube-root input must be nonnegative")
    low, high = 0, 1
    while high**3 <= value:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**3 <= value:
            low = middle
        else:
            high = middle
    return low


def _root_interval(digits: int = 12) -> tuple[Fraction, Fraction]:
    scale = 10**digits
    lower_numerator = _floor_cube_root(4 * scale**3)
    lower = Fraction(lower_numerator, scale)
    upper = Fraction(lower_numerator + 1, scale)
    if not (lower**3 <= 4 <= upper**3):
        raise ArithmeticError("cube-root interval construction failed")
    return lower, upper


def _interval_multiply(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    values = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(values), max(values)


def _cubic_interval(
    value: Cubic, root: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    root_squared = _interval_multiply(root, root)
    result = (value.a0, value.a0)
    for coefficient, interval in ((value.a1, root), (value.a2, root_squared)):
        scaled = (
            min(coefficient * interval[0], coefficient * interval[1]),
            max(coefficient * interval[0], coefficient * interval[1]),
        )
        result = result[0] + scaled[0], result[1] + scaled[1]
    return result


def _classify_sign(value: Cubic, root: tuple[Fraction, Fraction]) -> str:
    if value.is_zero():
        return "exact_zero"
    lower, upper = _cubic_interval(value, root)
    if lower > 0:
        return "positive"
    if upper < 0:
        return "negative"
    return "unresolved"


def _fraction_json(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [_fraction_json(coordinate) for coordinate in value.coordinates()]


def _word_records(polynomial: Mapping[Word, Cubic]) -> list[dict[str, object]]:
    return [
        {"word": list(word), "coefficient": _cubic_json(value)}
        for word, value in sorted(polynomial.items())
    ]


def _log_word_maps(logarithm: Mapping[Word, Cubic]) -> list[dict[str, object]]:
    expanded = _expand_lie(logarithm)
    return [
        {
            "degree": degree,
            "terms": _word_records(
                {word: value for word, value in expanded.items() if len(word) == degree}
            ),
        }
        for degree in range(MAX_DEGREE + 1)
    ]


def _stage_digest_json(stage: tuple[int, Cubic]) -> list[object]:
    return [stage[0], _cubic_json(stage[1])]


def _unsigned_conventions() -> dict[str, object]:
    return {
        "stage_table": [_stage_digest_json(stage) for stage in _literal_stages()],
        "multiplication_order": "left_to_right",
        "formal_s2": "S2(s)=exp(s*A/2)exp(s*B)exp(s*A/2)",
        "formal_s4": "S4(t)=S2(u*t)^2S2(v*t)S2(u*t)^2",
        "formal_log": "log S4(t)=t(A+B)+t^5*L5+O(t^7)",
        "parity_identity": "S4(-t)=S4(t)^(-1)",
        "trace_normalization": "unnormalized_Tr",
        "commutator_c": "[A,[A,B]]",
        "commutator_d": "[B,[B,A]]",
        "physical_substitution": "A,B -> -i*A,-i*B",
        "physical_log": "log U4(t)=-i*t*H-i*t^5*E5+O(t^7)",
        "physical_bridge": "E5=L5",
    }


def _convention_digest(conventions: Mapping[str, object]) -> str:
    encoded = json.dumps(
        conventions,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _comparison_json(comparison: ProjectiveComparison) -> dict[str, object]:
    return {
        "hypothesis": [_fraction_json(value) for value in comparison.hypothesis],
        "status": comparison.status,
        "scale": _cubic_json(comparison.scale),
        "residual": [_cubic_json(value) for value in comparison.residual],
    }


def _interval_json(interval: tuple[Fraction, Fraction]) -> list[list[int]]:
    return [_fraction_json(interval[0]), _fraction_json(interval[1])]


def _sign_entry(
    value: Cubic,
    root: tuple[Fraction, Fraction],
) -> dict[str, object] | None:
    if value.is_zero():
        zero = _fraction_json(Fraction())
        return {"status": "exact_zero", "interval": [zero, zero]}
    interval = _cubic_interval(value, root)
    status = _classify_sign(value, root)
    if status == "unresolved":
        return None
    return {"status": status, "interval": _interval_json(interval)}


def _sign_certificate(
    coefficients: Mapping[str, Cubic],
) -> dict[str, object]:
    for decimal_digits in (12, 24, 48, 96):
        root = _root_interval(decimal_digits)
        entries = {
            name: _sign_entry(value, root) for name, value in coefficients.items()
        }
        if all(entry is not None for entry in entries.values()):
            return {
                "minimal_polynomial": "alpha^3-4",
                "decimal_digits": decimal_digits,
                "alpha_interval": _interval_json(root),
                "coefficients": entries,
            }
    raise ArithmeticError("could not certify every nonzero PF4 coefficient sign")


def _canonical_json_bytes(payload: object) -> bytes:
    return (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _payload_digest(payload: Mapping[str, object]) -> str:
    unsigned = {key: value for key, value in payload.items() if key != "payload_sha256"}
    return hashlib.sha256(_canonical_json_bytes(unsigned)).hexdigest()


def _checked_source_path(relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ValueError(f"source path is not a relative path: {relative!r}")
    candidate = ISSUE_ROOT / relative
    if candidate.is_symlink():
        raise ValueError(f"source path must not be a symlink: {relative}")
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(ISSUE_ROOT.resolve(strict=True))
    except ValueError as error:
        raise ValueError(f"source path escapes issue root: {relative}") from error
    if not resolved.is_file():
        raise ValueError(f"source path is not a file: {relative}")
    return resolved


def _module_source(module: str) -> str | None:
    if not module or module.startswith("."):
        return None
    parts = module.split(".")
    for base in (ISSUE_ROOT / "src", ISSUE_ROOT):
        module_file = base.joinpath(*parts).with_suffix(".py")
        package_file = base.joinpath(*parts, "__init__.py")
        for candidate in (module_file, package_file):
            if candidate.is_file():
                return candidate.relative_to(ISSUE_ROOT).as_posix()
    return None


def _relative_import_module(source: str, node: ast.ImportFrom) -> str:
    source_path = Path(source)
    if source_path.parts[:1] != ("src",):
        raise ValueError(f"relative import outside src package: {source}")
    module_parts = list(source_path.with_suffix("").parts[1:-1])
    ascend = node.level - 1
    if ascend > len(module_parts):
        raise ValueError(f"relative import escapes package: {source}")
    if ascend:
        module_parts = module_parts[:-ascend]
    if node.module:
        module_parts.extend(node.module.split("."))
    return ".".join(module_parts)


def _local_source_imports(source: str) -> set[str]:
    path = _checked_source_path(source)
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=source)
    except (SyntaxError, UnicodeDecodeError) as error:
        raise ValueError(f"cannot parse source dependency {source}") from error
    imports: set[str] = set()
    for node in ast.walk(tree):
        modules: list[str] = []
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            module = (
                _relative_import_module(source, node)
                if node.level
                else (node.module or "")
            )
            modules.append(module)
        for module in modules:
            local = _module_source(module)
            if local is not None:
                imports.add(local)
            elif module.split(".", 1)[0] in {"trottercert", "scripts"}:
                raise ValueError(
                    f"unresolved local import {module!r} from source {source}"
                )
    return imports


def _source_closure() -> tuple[str, ...]:
    pending = list(SOURCE_ROOTS)
    closure: set[str] = set()
    while pending:
        relative = Path(pending.pop()).as_posix()
        _checked_source_path(relative)
        if relative in closure:
            continue
        closure.add(relative)
        if relative.endswith(".py"):
            pending.extend(sorted(_local_source_imports(relative) - closure))
    return tuple(sorted(closure))


def _implementation_sources() -> dict[str, str]:
    return {
        relative: hashlib.sha256(
            _checked_source_path(relative).read_bytes()
        ).hexdigest()
        for relative in _source_closure()
    }


def build_reference_artifact() -> dict[str, object]:
    mapping = derive_mapping()
    g_c2, g_cd, g_d2 = mapping.coefficients
    candidate = mapping.candidate_comparison
    gamma = candidate.scale if candidate.status == "proportional" else None
    unsigned_conventions = _unsigned_conventions()
    conventions = dict(unsigned_conventions)
    conventions.pop("stage_table")
    conventions["convention_digest"] = _convention_digest(unsigned_conventions)
    sign_values = {"g_c2": g_c2, "g_cd": g_cd, "g_d2": g_d2}
    if gamma is not None:
        sign_values["optional_gamma"] = gamma
    payload: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "kind": ARTIFACT_KIND,
        "conventions": conventions,
        "stages": [
            {"fragment_index": fragment, "coefficient": _cubic_json(coefficient)}
            for fragment, coefficient in _literal_stages()
        ],
        "log_terms": _log_word_maps(mapping.logarithm),
        "cyclic_map": {
            "target_h_l5": _word_records(mapping.target_cyclic),
            "basis_c2": _word_records(mapping.c2_cyclic),
            "basis_cd": _word_records(mapping.cd_cyclic),
            "basis_d2": _word_records(mapping.d2_cyclic),
            "reconstructed": _word_records(mapping.reconstructed_cyclic),
        },
        "solved_coefficients": {
            "g_c2": _cubic_json(g_c2),
            "g_cd": _cubic_json(g_cd),
            "g_d2": _cubic_json(g_d2),
        },
        "residual": _word_records(mapping.residual),
        "projective_comparisons": {
            "historical_candidate": _comparison_json(candidate),
            "legacy": _comparison_json(mapping.legacy_comparison),
        },
        "optional_gamma": None if gamma is None else _cubic_json(gamma),
        "sign_certificate": _sign_certificate(sign_values),
        "physical_bridge": {
            "substitution": "A,B -> -i*A,-i*B",
            "degree_one_factor": [0, -1],
            "degree_five_factor": [0, -1],
            "physical_log": "log U4(t)=-i*t*H-i*t^5*E5+O(t^7)",
            "relation": "E5=L5",
        },
        "identity_status": IDENTITY_STATUS,
        "payload_sha256": "",
    }
    payload["payload_sha256"] = _payload_digest(payload)
    return payload


def build_reference_wrapper() -> dict[str, object]:
    """Build v2 independently from the mathematical replay and source closure."""

    payload: dict[str, object] = {
        "schema_version": WRAPPER_SCHEMA_VERSION,
        "kind": WRAPPER_KIND,
        "mapping": build_reference_artifact(),
        "implementation_sources": _implementation_sources(),
        "reference_algorithm": REFERENCE_ALGORITHM,
        "payload_sha256": "",
    }
    payload["payload_sha256"] = _payload_digest(payload)
    return payload


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _load_artifact(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON: {error}") from error
    if not isinstance(payload, dict):
        raise TypeError("artifact must be a JSON object")
    if raw != _canonical_json_bytes(payload):
        raise ValueError("artifact bytes are not canonical JSON")
    return payload


def _reject_json_numeric_aliases(value: object, path: str = "artifact") -> None:
    if isinstance(value, float):
        raise TypeError(f"{path}: floating-point values are forbidden")
    if isinstance(value, bool):
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _reject_json_numeric_aliases(item, f"{path}[{index}]")
    elif isinstance(value, dict):
        for key, item in value.items():
            _reject_json_numeric_aliases(item, f"{path}.{key}")


def _verify_digest(payload: Mapping[str, object], label: str) -> None:
    digest = payload.get("payload_sha256")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise ValueError(f"{label} payload_sha256 must be a 64-hex string")
    if digest != _payload_digest(payload):
        raise ValueError(f"{label} payload digest mismatch")


def _verify_v1_mapping(payload: Mapping[str, object]) -> None:
    expected = build_reference_artifact()
    if set(payload) != set(expected):
        missing = sorted(set(expected) - set(payload))
        extra = sorted(set(payload) - set(expected))
        raise ValueError(f"top-level schema mismatch missing={missing} extra={extra}")
    _verify_digest(payload, "mapping")
    if dict(payload) != expected:
        for key in expected:
            if payload.get(key) != expected[key]:
                raise ValueError(f"semantic mismatch in {key}")
        raise ValueError("artifact differs from independent derivation")


def _verify_v2_wrapper(payload: Mapping[str, object]) -> None:
    if set(payload) != WRAPPER_FIELDS:
        raise ValueError("PF4 wrapper field set mismatch")
    if payload.get("schema_version") != WRAPPER_SCHEMA_VERSION:
        raise ValueError("PF4 wrapper schema version mismatch")
    if payload.get("kind") != WRAPPER_KIND:
        raise ValueError("PF4 wrapper kind mismatch")
    if payload.get("reference_algorithm") != REFERENCE_ALGORITHM:
        raise ValueError("PF4 reference algorithm mismatch")
    _verify_digest(payload, "wrapper")
    mapping = payload.get("mapping")
    if not isinstance(mapping, Mapping):
        raise TypeError("PF4 mathematical mapping is missing")
    _verify_v1_mapping(mapping)
    sources = payload.get("implementation_sources")
    if not isinstance(sources, Mapping) or any(
        not isinstance(path, str)
        or not isinstance(source_digest, str)
        or len(source_digest) != 64
        or any(character not in "0123456789abcdef" for character in source_digest)
        for path, source_digest in sources.items()
    ):
        raise ValueError("PF4 implementation source manifest is malformed")
    if dict(sources) != _implementation_sources():
        raise ValueError("PF4 implementation source closure mismatch")
    expected = build_reference_wrapper()
    if dict(payload) != expected:
        raise ValueError("PF4 wrapper differs from independent rebuild")


def verify_artifact_payload(payload: object) -> list[str]:
    try:
        if not isinstance(payload, dict):
            raise TypeError("artifact must be an object")
        _reject_json_numeric_aliases(payload)
        if payload.get("schema_version") == WRAPPER_SCHEMA_VERSION:
            _verify_v2_wrapper(payload)
        else:
            _verify_v1_mapping(payload)
    except (ArithmeticError, OSError, TypeError, ValueError) as error:
        return [str(error)]
    return []


def verify_artifact(path: Path) -> list[str]:
    try:
        payload = _load_artifact(path)
    except (OSError, TypeError, ValueError) as error:
        return [str(error)]
    return verify_artifact_payload(payload)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path, required=True, metavar="ARTIFACT")
    arguments = parser.parse_args(argv)
    errors = verify_artifact(arguments.verify)
    if errors:
        print("\n".join(errors))
        return 1
    submitted = _load_artifact(arguments.verify)
    mapping = (
        submitted["mapping"]
        if submitted.get("schema_version") == WRAPPER_SCHEMA_VERSION
        else submitted
    )
    print(
        json.dumps(
            {
                "valid": True,
                "algorithm": REFERENCE_ALGORITHM,
                "identity_status": mapping["identity_status"],
                "payload_sha256": submitted["payload_sha256"],
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
