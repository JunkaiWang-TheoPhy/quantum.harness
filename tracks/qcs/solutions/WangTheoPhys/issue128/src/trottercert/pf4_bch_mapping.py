"""Exact primary BCH-to-trace mapping for the two-fragment Suzuki PF4.

The authoritative calculation in this module starts from the explicit merged
eleven-stage product.  It expands the product and its formal logarithm in the
free associative algebra over ``Q(alpha)``, reduces ``Tr((A+B)L5)`` modulo
cyclic word rotations, and solves for the coefficients of ``Tr(C^2)``,
``Tr(CD)``, and ``Tr(D^2)`` without supplying an expected coefficient vector.

Historical helpers such as :func:`pf4_suzuki_gamma` remain compatibility
wrappers.  They are derived from the general solved triple and fail closed if
the historical projective form is not what the exact solve produces.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import factorial

from .cubic_field import Cubic, CubicStage
from .intervals import RationalInterval, cube_root_four_interval

Word = tuple[int, ...]
CubicWordPolynomial = dict[Word, Cubic]
CubicWordSeries = tuple[CubicWordPolynomial, ...]

_CANDIDATE_VECTOR = (Fraction(1), Fraction(-4), Fraction(8, 3))
_LEGACY_VECTOR = (Fraction(1, 2), Fraction(14, 3), Fraction(4, 3))


@dataclass(frozen=True, slots=True)
class PF4ConventionRecord:
    """Digest-bound formal, physical, trace, and commutator conventions."""

    stage_table: tuple[CubicStage, ...]
    multiplication_order: str
    formal_s2: str
    formal_s4: str
    formal_log: str
    parity_identity: str
    trace_normalization: str
    commutator_c: str
    commutator_d: str
    physical_substitution: str
    physical_log: str
    physical_bridge: str
    convention_digest: str


@dataclass(frozen=True, slots=True)
class ProjectiveComparison:
    """Exact comparison of a solved cubic triple with a rational hypothesis."""

    hypothesis: tuple[Fraction, Fraction, Fraction]
    status: str
    scale: Cubic
    residual: tuple[Cubic, Cubic, Cubic]


@dataclass(frozen=True, slots=True)
class PF4CyclicMapping:
    """Complete primary exact solution of the degree-five cyclic mapping."""

    convention: PF4ConventionRecord
    log_terms: CubicWordSeries
    target_cyclic: CubicWordPolynomial
    c2_cyclic: CubicWordPolynomial
    cd_cyclic: CubicWordPolynomial
    d2_cyclic: CubicWordPolynomial
    g_c2: Cubic
    g_cd: Cubic
    g_d2: Cubic
    reconstructed_cyclic: CubicWordPolynomial
    residual: CubicWordPolynomial
    candidate_comparison: ProjectiveComparison
    legacy_comparison: ProjectiveComparison
    payload_sha256: str

    @property
    def coefficients(self) -> tuple[Cubic, Cubic, Cubic]:
        """Return the authoritative general coefficient triple."""

        return self.g_c2, self.g_cd, self.g_d2


def pf4_two_fragment_stages() -> tuple[CubicStage, ...]:
    """Return the literal left-to-right eleven-stage two-fragment PF4 table."""

    alpha = Cubic(0, 1, 0)
    u = Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
    v = Cubic.one() - 4 * u
    if (Cubic(4, 0, 0) - alpha) * u != Cubic.one() or v != -alpha * u:
        raise ArithmeticError("exact Suzuki coefficient identities failed")
    return (
        CubicStage(0, u / 2),
        CubicStage(1, u),
        CubicStage(0, u),
        CubicStage(1, u),
        CubicStage(0, (u + v) / 2),
        CubicStage(1, v),
        CubicStage(0, (u + v) / 2),
        CubicStage(1, u),
        CubicStage(0, u),
        CubicStage(1, u),
        CubicStage(0, u / 2),
    )


def physical_substitution_factor(degree: int) -> tuple[int, int]:
    """Return ``(-i)**degree`` as the exact Gaussian pair ``(real, imag)``."""

    if isinstance(degree, bool) or not isinstance(degree, int):
        raise TypeError("degree must be an integer, not bool")
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    return ((1, 0), (0, -1), (-1, 0), (0, 1))[degree % 4]


def _fraction_payload(value: Fraction) -> list[int]:
    exact = Fraction(value)
    return [exact.numerator, exact.denominator]


def _cubic_payload(value: Cubic) -> list[list[int]]:
    return [
        _fraction_payload(value.a0),
        _fraction_payload(value.a1),
        _fraction_payload(value.a2),
    ]


def _stage_payload(stage: CubicStage) -> list[object]:
    return [stage.fragment_index, _cubic_payload(stage.coefficient)]


def _convention_payload(record: PF4ConventionRecord) -> dict[str, object]:
    return {
        "stage_table": [_stage_payload(stage) for stage in record.stage_table],
        "multiplication_order": record.multiplication_order,
        "formal_s2": record.formal_s2,
        "formal_s4": record.formal_s4,
        "formal_log": record.formal_log,
        "parity_identity": record.parity_identity,
        "trace_normalization": record.trace_normalization,
        "commutator_c": record.commutator_c,
        "commutator_d": record.commutator_d,
        "physical_substitution": record.physical_substitution,
        "physical_log": record.physical_log,
        "physical_bridge": record.physical_bridge,
    }


def convention_record_digest(record: PF4ConventionRecord) -> str:
    """Hash every semantic convention field, excluding the supplied digest."""

    encoded = json.dumps(
        _convention_payload(record),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return sha256(encoded).hexdigest()


@lru_cache(maxsize=1)
def pf4_convention_record() -> PF4ConventionRecord:
    """Return the canonical convention record for the primary derivation."""

    stages = pf4_two_fragment_stages()
    if stages != tuple(reversed(stages)):
        raise ArithmeticError("the frozen PF4 table is not palindromic")
    if physical_substitution_factor(1) != physical_substitution_factor(5):
        raise ArithmeticError("formal-to-physical degree-five sign bridge failed")
    unsigned = PF4ConventionRecord(
        stage_table=stages,
        multiplication_order="left_to_right",
        formal_s2="S2(s)=exp(s*A/2)exp(s*B)exp(s*A/2)",
        formal_s4="S4(t)=S2(u*t)^2S2(v*t)S2(u*t)^2",
        formal_log="log S4(t)=t(A+B)+t^5*L5+O(t^7)",
        parity_identity="S4(-t)=S4(t)^(-1)",
        trace_normalization="unnormalized_Tr",
        commutator_c="[A,[A,B]]",
        commutator_d="[B,[B,A]]",
        physical_substitution="A,B -> -i*A,-i*B",
        physical_log="log U4(t)=-i*t*H-i*t^5*E5+O(t^7)",
        physical_bridge="E5=L5",
        convention_digest="",
    )
    return replace(unsigned, convention_digest=convention_record_digest(unsigned))


def verify_pf4_convention_record(record: PF4ConventionRecord) -> bool:
    """Verify the digest and exact canonical semantics, not the digest alone."""

    if not isinstance(record, PF4ConventionRecord):
        return False
    canonical = pf4_convention_record()
    return (
        record == canonical
        and convention_record_digest(record) == record.convention_digest
    )


def _add(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    result = dict(left)
    for word, coefficient in right.items():
        updated = result.get(word, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            result.pop(word, None)
        else:
            result[word] = updated
    return result


def _scale(
    polynomial: Mapping[Word, Cubic],
    scalar: Cubic | Fraction | int,
) -> CubicWordPolynomial:
    factor = Cubic.coerce(scalar)
    return {
        word: scaled
        for word, coefficient in polynomial.items()
        if (scaled := coefficient * factor) != Cubic.zero()
    }


def _multiply(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    result: CubicWordPolynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            updated = result.get(word, Cubic.zero()) + (
                left_coefficient * right_coefficient
            )
            if updated == Cubic.zero():
                result.pop(word, None)
            else:
                result[word] = updated
    return result


def _commutator(
    left: Mapping[Word, Cubic],
    right: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    return _add(_multiply(left, right), _scale(_multiply(right, left), -1))


def _empty_series(order: int) -> list[CubicWordPolynomial]:
    return [{} for _ in range(order + 1)]


def _identity_series(order: int) -> list[CubicWordPolynomial]:
    result = _empty_series(order)
    result[0][()] = Cubic.one()
    return result


def _series_add(
    left: Sequence[Mapping[Word, Cubic]],
    right: Sequence[Mapping[Word, Cubic]],
) -> list[CubicWordPolynomial]:
    if len(left) != len(right):
        raise ValueError("word-series truncation orders disagree")
    return [_add(lhs, rhs) for lhs, rhs in zip(left, right)]


def _series_scale(
    series: Sequence[Mapping[Word, Cubic]],
    scalar: Fraction | int,
) -> list[CubicWordPolynomial]:
    return [_scale(degree, scalar) for degree in series]


def _series_multiply(
    left: Sequence[Mapping[Word, Cubic]],
    right: Sequence[Mapping[Word, Cubic]],
) -> list[CubicWordPolynomial]:
    if len(left) != len(right):
        raise ValueError("word-series truncation orders disagree")
    order = len(left) - 1
    result = _empty_series(order)
    for total_degree in range(order + 1):
        degree_result: CubicWordPolynomial = {}
        for left_degree in range(total_degree + 1):
            degree_result = _add(
                degree_result,
                _multiply(
                    left[left_degree],
                    right[total_degree - left_degree],
                ),
            )
        result[total_degree] = degree_result
    return result


def _stage_exponential(stage: CubicStage, order: int) -> list[CubicWordPolynomial]:
    result = _empty_series(order)
    for degree in range(order + 1):
        result[degree][(stage.fragment_index,) * degree] = (
            stage.coefficient**degree / factorial(degree)
        )
    return result


def _formula_log_series(
    stages: Sequence[CubicStage],
    order: int,
) -> CubicWordSeries:
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("order must be an integer, not bool")
    if order < 0:
        raise ValueError("order must be nonnegative")
    product = _identity_series(order)
    for stage in stages:
        product = _series_multiply(product, _stage_exponential(stage, order))
    delta = [dict(degree) for degree in product]
    delta[0] = _add(delta[0], {(): -Cubic.one()})

    logarithm = _empty_series(order)
    power = _identity_series(order)
    for exponent in range(1, order + 1):
        power = _series_multiply(power, delta)
        logarithm = _series_add(
            logarithm,
            _series_scale(
                power,
                Fraction(1 if exponent % 2 else -1, exponent),
            ),
        )
    return tuple(logarithm)


def cyclic_trace_classes(
    polynomial: Mapping[Word, Cubic],
) -> CubicWordPolynomial:
    """Collect a word polynomial modulo exact cyclic trace equivalence."""

    classes: CubicWordPolynomial = {}
    for word, coefficient in polynomial.items():
        representative = (
            min(word[offset:] + word[:offset] for offset in range(len(word)))
            if word
            else ()
        )
        updated = classes.get(representative, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            classes.pop(representative, None)
        else:
            classes[representative] = updated
    return classes


def _free_trace_objects() -> tuple[
    CubicWordPolynomial,
    CubicWordPolynomial,
    CubicWordPolynomial,
    CubicWordPolynomial,
    CubicWordPolynomial,
]:
    a = {(0,): Cubic.one()}
    b = {(1,): Cubic.one()}
    hamiltonian = _add(a, b)
    c = _commutator(a, _commutator(a, b))
    d = _commutator(b, _commutator(b, a))
    return hamiltonian, c, d, a, b


def _rational_basis_coefficient(value: Cubic) -> Fraction:
    if value.a1 or value.a2:
        raise ArithmeticError("cyclic basis unexpectedly left the rational field")
    return value.a0


def _solve_cyclic_coefficients(
    target: Mapping[Word, Cubic],
    basis: Sequence[Mapping[Word, Cubic]],
) -> tuple[Cubic, Cubic, Cubic]:
    if len(basis) != 3:
        raise ValueError("the PF4 cyclic solve requires exactly three basis maps")
    words = sorted(set(target).union(*(set(polynomial) for polynomial in basis)))
    rows: list[tuple[list[Fraction], Cubic]] = []
    for word in words:
        coefficients = [
            _rational_basis_coefficient(polynomial.get(word, Cubic.zero()))
            for polynomial in basis
        ]
        rows.append((coefficients, target.get(word, Cubic.zero())))

    pivot_row = 0
    pivots: list[int] = []
    for column in range(3):
        selected = next(
            (
                row_index
                for row_index in range(pivot_row, len(rows))
                if rows[row_index][0][column]
            ),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        coefficients, rhs = rows[pivot_row]
        pivot = coefficients[column]
        coefficients = [value / pivot for value in coefficients]
        rhs = rhs / pivot
        rows[pivot_row] = (coefficients, rhs)

        for row_index, (other_coefficients, other_rhs) in enumerate(rows):
            if row_index == pivot_row:
                continue
            factor = other_coefficients[column]
            if not factor:
                continue
            rows[row_index] = (
                [
                    value - factor * pivot_value
                    for value, pivot_value in zip(other_coefficients, coefficients)
                ],
                other_rhs - factor * rhs,
            )
        pivots.append(column)
        pivot_row += 1

    for coefficients, rhs in rows:
        if not any(coefficients) and rhs != Cubic.zero():
            raise ArithmeticError("PF4 cyclic coefficient system is inconsistent")
    if pivots != [0, 1, 2]:
        raise ArithmeticError("PF4 cyclic coefficient system is not uniquely solvable")

    solution = [Cubic.zero(), Cubic.zero(), Cubic.zero()]
    for coefficients, rhs in rows:
        nonzero = [index for index, value in enumerate(coefficients) if value]
        if len(nonzero) == 1 and coefficients[nonzero[0]] == 1:
            solution[nonzero[0]] = rhs
    return solution[0], solution[1], solution[2]


def compare_projectively(
    solved: Sequence[Cubic],
    hypothesis: Sequence[Fraction | int],
) -> ProjectiveComparison:
    """Compare a solved triple to a rational direction after the exact solve."""

    if len(solved) != 3 or len(hypothesis) != 3:
        raise ValueError("projective comparisons require two triples")
    solved_triple = tuple(Cubic.coerce(value) for value in solved)
    hypothesis_triple = tuple(Fraction(value) for value in hypothesis)
    if not any(hypothesis_triple):
        raise ValueError("projective hypothesis must be nonzero")
    index = next(i for i, value in enumerate(hypothesis_triple) if value)
    scale = solved_triple[index] / hypothesis_triple[index]
    residual = tuple(
        value - scale * direction
        for value, direction in zip(solved_triple, hypothesis_triple)
    )
    proportional = all(value == Cubic.zero() for value in residual) and any(
        value != Cubic.zero() for value in solved_triple
    )
    return ProjectiveComparison(
        hypothesis=(
            hypothesis_triple[0],
            hypothesis_triple[1],
            hypothesis_triple[2],
        ),
        status="proportional" if proportional else "not_proportional",
        scale=scale,
        residual=(residual[0], residual[1], residual[2]),
    )


def _word_map_payload(polynomial: Mapping[Word, Cubic]) -> list[list[object]]:
    return [
        [list(word), _cubic_payload(coefficient)]
        for word, coefficient in sorted(polynomial.items())
    ]


def _comparison_payload(comparison: ProjectiveComparison) -> dict[str, object]:
    return {
        "hypothesis": [_fraction_payload(value) for value in comparison.hypothesis],
        "status": comparison.status,
        "scale": _cubic_payload(comparison.scale),
        "residual": [_cubic_payload(value) for value in comparison.residual],
    }


def _mapping_digest_payload(mapping: PF4CyclicMapping) -> dict[str, object]:
    return {
        "convention_digest": mapping.convention.convention_digest,
        "log_terms": [
            _word_map_payload(polynomial) for polynomial in mapping.log_terms
        ],
        "target_cyclic": _word_map_payload(mapping.target_cyclic),
        "c2_cyclic": _word_map_payload(mapping.c2_cyclic),
        "cd_cyclic": _word_map_payload(mapping.cd_cyclic),
        "d2_cyclic": _word_map_payload(mapping.d2_cyclic),
        "coefficients": [_cubic_payload(value) for value in mapping.coefficients],
        "reconstructed_cyclic": _word_map_payload(mapping.reconstructed_cyclic),
        "residual": _word_map_payload(mapping.residual),
        "candidate_comparison": _comparison_payload(mapping.candidate_comparison),
        "legacy_comparison": _comparison_payload(mapping.legacy_comparison),
    }


def _mapping_digest(mapping: PF4CyclicMapping) -> str:
    encoded = json.dumps(
        _mapping_digest_payload(mapping),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return sha256(encoded).hexdigest()


def _map_records(polynomial: Mapping[Word, Cubic]) -> list[dict[str, object]]:
    return [
        {"word": list(word), "coefficient": _cubic_payload(coefficient)}
        for word, coefficient in sorted(polynomial.items())
    ]


def _interval_payload(interval: RationalInterval) -> list[list[int]]:
    return [_fraction_payload(interval.lower), _fraction_payload(interval.upper)]


def _sign_entry(value: Cubic, root: RationalInterval) -> dict[str, object] | None:
    if value == Cubic.zero():
        return {
            "status": "exact_zero",
            "interval": [_fraction_payload(Fraction()), _fraction_payload(Fraction())],
        }
    enclosure = value.enclose(root)
    if enclosure.lower > 0:
        status = "positive"
    elif enclosure.upper < 0:
        status = "negative"
    else:
        return None
    return {"status": status, "interval": _interval_payload(enclosure)}


def _sign_certificate(mapping: PF4CyclicMapping) -> dict[str, object]:
    values = {
        "g_c2": mapping.g_c2,
        "g_cd": mapping.g_cd,
        "g_d2": mapping.g_d2,
    }
    if mapping.candidate_comparison.status == "proportional":
        values["optional_gamma"] = mapping.candidate_comparison.scale
    for decimal_digits in (12, 24, 48, 96):
        root = cube_root_four_interval(decimal_digits)
        entries = {name: _sign_entry(value, root) for name, value in values.items()}
        if all(entry is not None for entry in entries.values()):
            return {
                "minimal_polynomial": "alpha^3-4",
                "decimal_digits": decimal_digits,
                "alpha_interval": _interval_payload(root),
                "coefficients": entries,
            }
    raise ArithmeticError("could not certify every nonzero PF4 coefficient sign")


def canonical_payload_bytes(payload: Mapping[str, object]) -> bytes:
    """Serialize a mathematical mapping payload in canonical compact JSON."""

    try:
        encoded = json.dumps(
            dict(payload),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise ValueError("PF4 mapping payload is not canonical JSON data") from error
    return (encoded + "\n").encode("ascii")


def payload_digest(payload: Mapping[str, object]) -> str:
    """Return the SHA-256 of the canonical mathematical payload without its digest."""

    unsigned = {key: value for key, value in payload.items() if key != "payload_sha256"}
    return sha256(canonical_payload_bytes(unsigned)).hexdigest()


def to_payload(mapping: PF4CyclicMapping | None = None) -> dict[str, object]:
    """Return the reference-compatible mathematical payload for ``mapping``.

    Source closure is intentionally deferred to the Task-5 artifact wrapper;
    every mathematical field needed by an independent implementation is
    present here.
    """

    solved = derive_pf4_cyclic_mapping() if mapping is None else mapping
    convention_payload = _convention_payload(solved.convention)
    convention_payload.pop("stage_table")
    convention_payload["convention_digest"] = solved.convention.convention_digest
    optional_gamma = (
        _cubic_payload(solved.candidate_comparison.scale)
        if solved.candidate_comparison.status == "proportional"
        else None
    )
    payload: dict[str, object] = {
        "schema_version": 1,
        "kind": "pf4_two_fragment_exact_cyclic_mapping",
        "conventions": convention_payload,
        "stages": [
            {
                "fragment_index": stage.fragment_index,
                "coefficient": _cubic_payload(stage.coefficient),
            }
            for stage in solved.convention.stage_table
        ],
        "log_terms": [
            {"degree": degree, "terms": _map_records(polynomial)}
            for degree, polynomial in enumerate(solved.log_terms)
        ],
        "cyclic_map": {
            "target_h_l5": _map_records(solved.target_cyclic),
            "basis_c2": _map_records(solved.c2_cyclic),
            "basis_cd": _map_records(solved.cd_cyclic),
            "basis_d2": _map_records(solved.d2_cyclic),
            "reconstructed": _map_records(solved.reconstructed_cyclic),
        },
        "solved_coefficients": {
            "g_c2": _cubic_payload(solved.g_c2),
            "g_cd": _cubic_payload(solved.g_cd),
            "g_d2": _cubic_payload(solved.g_d2),
        },
        "residual": _map_records(solved.residual),
        "projective_comparisons": {
            "historical_candidate": _comparison_payload(
                solved.candidate_comparison
            ),
            "legacy": _comparison_payload(solved.legacy_comparison),
        },
        "optional_gamma": optional_gamma,
        "sign_certificate": _sign_certificate(solved),
        "physical_bridge": {
            "substitution": solved.convention.physical_substitution,
            "degree_one_factor": list(physical_substitution_factor(1)),
            "degree_five_factor": list(physical_substitution_factor(5)),
            "physical_log": solved.convention.physical_log,
            "relation": solved.convention.physical_bridge,
        },
        "identity_status": "proved_exact_suzuki_pf4_free_trace_identity",
        "payload_sha256": "",
    }
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def parse_payload(payload: Mapping[str, object]) -> PF4CyclicMapping:
    """Strictly parse and authoritatively verify a mathematical mapping payload."""

    if not isinstance(payload, Mapping):
        raise TypeError("PF4 mapping payload must be a mapping")
    expected_fields = {
        "schema_version",
        "kind",
        "conventions",
        "stages",
        "log_terms",
        "cyclic_map",
        "solved_coefficients",
        "residual",
        "projective_comparisons",
        "optional_gamma",
        "sign_certificate",
        "physical_bridge",
        "identity_status",
        "payload_sha256",
    }
    if set(payload) != expected_fields:
        raise ValueError("PF4 mapping payload field set mismatch")
    digest = payload.get("payload_sha256")
    if not isinstance(digest, str) or digest != payload_digest(payload):
        raise ValueError("PF4 mapping payload digest mismatch")
    expected = to_payload()
    if canonical_payload_bytes(payload) != canonical_payload_bytes(expected):
        raise ValueError("PF4 mapping payload differs from authoritative rebuild")
    return derive_pf4_cyclic_mapping()


@lru_cache(maxsize=1)
def derive_pf4_cyclic_mapping() -> PF4CyclicMapping:
    """Derive and solve the exact PF4 cyclic mapping with no expected vector."""

    convention = pf4_convention_record()
    if not verify_pf4_convention_record(convention):
        raise ArithmeticError("canonical PF4 convention record is invalid")
    logarithm = _formula_log_series(convention.stage_table, 5)
    hamiltonian, c, d, _, _ = _free_trace_objects()
    if logarithm[1] != hamiltonian or any(
        logarithm[degree] for degree in (2, 3, 4)
    ):
        raise ArithmeticError("exact Suzuki stages are not fourth order")
    if not logarithm[5]:
        raise ArithmeticError("PF4 degree-five logarithm unexpectedly vanished")

    target = cyclic_trace_classes(_multiply(hamiltonian, logarithm[5]))
    c2 = cyclic_trace_classes(_multiply(c, c))
    cd = cyclic_trace_classes(_multiply(c, d))
    d2 = cyclic_trace_classes(_multiply(d, d))
    coefficients = _solve_cyclic_coefficients(target, (c2, cd, d2))
    reconstructed = _add(
        _scale(c2, coefficients[0]),
        _add(_scale(cd, coefficients[1]), _scale(d2, coefficients[2])),
    )
    residual = _add(target, _scale(reconstructed, -1))
    if residual:
        raise ArithmeticError("PF4 cyclic reconstruction has a nonzero residual")

    candidate = compare_projectively(coefficients, _CANDIDATE_VECTOR)
    legacy = compare_projectively(coefficients, _LEGACY_VECTOR)
    unsigned = PF4CyclicMapping(
        convention=convention,
        log_terms=logarithm,
        target_cyclic=target,
        c2_cyclic=c2,
        cd_cyclic=cd,
        d2_cyclic=d2,
        g_c2=coefficients[0],
        g_cd=coefficients[1],
        g_d2=coefficients[2],
        reconstructed_cyclic=reconstructed,
        residual=residual,
        candidate_comparison=candidate,
        legacy_comparison=legacy,
        payload_sha256="",
    )
    return replace(unsigned, payload_sha256=_mapping_digest(unsigned))


def pf4_suzuki_gamma() -> Cubic:
    """Return the historical prefactor, derived from the general exact solve."""

    comparison = derive_pf4_cyclic_mapping().candidate_comparison
    if comparison.status != "proportional":
        raise ArithmeticError(
            "the exact PF4 mapping is not proportional to the historical candidate"
        )
    return comparison.scale


def pf4_suzuki_trace_polynomials() -> tuple[
    CubicWordPolynomial,
    CubicWordPolynomial,
]:
    """Return raw free-word polynomials for the derived trace identity."""

    mapping = derive_pf4_cyclic_mapping()
    hamiltonian, c, d, _, _ = _free_trace_objects()
    right = _add(
        _scale(_multiply(c, c), mapping.g_c2),
        _add(
            _scale(_multiply(c, d), mapping.g_cd),
            _scale(_multiply(d, d), mapping.g_d2),
        ),
    )
    return _multiply(hamiltonian, mapping.log_terms[5]), right


def verify_pf4_suzuki_trace_identity() -> None:
    """Raise unless the general solve and compatibility form both verify."""

    mapping = derive_pf4_cyclic_mapping()
    if mapping.residual or mapping.target_cyclic != mapping.reconstructed_cyclic:
        raise ArithmeticError("exact PF4 general cyclic mapping failed")
    if mapping.candidate_comparison.status != "proportional":
        raise ArithmeticError("historical PF4 compatibility form is not derived")
    left, right = pf4_suzuki_trace_polynomials()
    left_classes = cyclic_trace_classes(left)
    right_classes = cyclic_trace_classes(right)
    if len(left_classes) != 10 or left_classes != right_classes:
        raise ArithmeticError("exact PF4 cyclic free-trace identity failed")
