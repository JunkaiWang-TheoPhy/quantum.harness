"""Exact trace obstructions for symmetric product formulas.

The PF2 identity in this module is an exact free-trace identity.  The PF4
quadratic form is deliberately narrower: it evaluates the coefficients stated
in the research plan exactly, but does not certify an as-yet absent derivation
connecting those moments to a particular fourth-order BCH defect.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from hashlib import sha256
from numbers import Integral
from typing import Any

import numpy as np


class ObstructionException(str, Enum):
    """Logically distinct ways a trace-obstruction argument can branch."""

    COMMUTING = "commuting"
    SYMMETRY_PROTECTED_POSITIVE = "symmetry_protected_positive"
    MIXED_TERM_INDEFINITE = "mixed_term_indefinite"
    EXACTLY_CORRECTABLE = "exactly_correctable"


@dataclass(frozen=True)
class TraceIdentityRecord:
    """Digest-bound statement of what an implemented identity establishes."""

    identity_kind: str
    formula: str
    assumptions: tuple[str, ...]
    exact_coefficients: tuple[tuple[str, Fraction], ...]
    exception_conditions: tuple[tuple[str, str], ...]
    identity_status: str
    identity_digest: str


def _exact_fraction(value: Any, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact rational, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    raise TypeError(f"{name} must be an exact rational (Fraction or integer)")


def _is_exact_rational_matrix(matrix: np.ndarray) -> bool:
    return all(
        not isinstance(entry, (bool, np.bool_))
        and isinstance(entry, (Fraction, Integral))
        for entry in matrix.flat
    )


def _validate_matrix_pair(
    a: np.ndarray,
    b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, bool]:
    left = np.asarray(a)
    right = np.asarray(b)
    if left.ndim != 2 or right.ndim != 2:
        raise ValueError("A and B must be two-dimensional matrices")
    if left.shape != right.shape or left.shape[0] != left.shape[1]:
        raise ValueError("A and B must be square matrices of the same shape")
    if any(
        isinstance(entry, (bool, np.bool_))
        for matrix in (left, right)
        for entry in matrix.flat
    ):
        raise TypeError("boolean matrix entries are not permitted")

    exact = _is_exact_rational_matrix(left) and _is_exact_rational_matrix(right)
    if exact:
        exact_left = np.asarray(
            [[Fraction(entry) for entry in row] for row in left], dtype=object
        )
        exact_right = np.asarray(
            [[Fraction(entry) for entry in row] for row in right], dtype=object
        )
        if not np.array_equal(exact_left, exact_left.T) or not np.array_equal(
            exact_right, exact_right.T
        ):
            raise ValueError("A and B must be Hermitian")
        return exact_left, exact_right, True

    try:
        numeric_left = np.asarray(left, dtype=np.complex128)
        numeric_right = np.asarray(right, dtype=np.complex128)
    except (TypeError, ValueError) as error:
        raise TypeError("matrix entries must be numeric") from error
    if not np.isfinite(numeric_left).all() or not np.isfinite(numeric_right).all():
        raise ValueError("matrix entries must be finite")
    if not np.allclose(numeric_left, numeric_left.conj().T) or not np.allclose(
        numeric_right, numeric_right.conj().T
    ):
        raise ValueError("A and B must be Hermitian")
    return numeric_left, numeric_right, False


def pf2_trace_obstruction(a: np.ndarray, b: np.ndarray) -> Fraction | float:
    """Return the exact leading Strang trace obstruction.

    For the convention

    ``log(exp(t A/2) exp(t B) exp(t A/2)) = t(A+B) + t^3 L3 + ...``,

    cyclicity of trace gives

    ``Tr((A+B)L3) = Tr([A,B]^dagger [A,B]) / 24``.

    Thus the result is nonnegative and is zero exactly when the finite
    dimensional Hermitian matrices commute.  Integer/Fraction matrices use
    exact arithmetic; other numeric Hermitian matrices return a float.
    """

    left, right, exact = _validate_matrix_pair(a, b)
    commutator = left @ right - right @ left
    if exact:
        norm_square = sum(
            (entry * entry for entry in commutator.flat), Fraction()
        )
        return norm_square / 24
    norm_square = float(np.vdot(commutator, commutator).real)
    return norm_square / 24.0


def pf4_trace_quadratic_form(
    trace_c2: Fraction | int,
    trace_cd: Fraction | int,
    trace_d2: Fraction | int,
    gamma: Fraction | int,
) -> Fraction:
    """Evaluate the plan-stated PF4 moment form with exact coefficients.

    This function is an exact algebraic evaluator.  It intentionally makes no
    claim that ``C`` and ``D`` have been connected to a fourth-order BCH defect;
    that missing derivation is recorded by :func:`pf4_trace_identity_record`.
    """

    c2 = _exact_fraction(trace_c2, "trace_c2")
    cd = _exact_fraction(trace_cd, "trace_cd")
    d2 = _exact_fraction(trace_d2, "trace_d2")
    scale = _exact_fraction(gamma, "gamma")
    return scale * (
        Fraction(1, 2) * c2
        + Fraction(14, 3) * cd
        + Fraction(4, 3) * d2
    )


def classify_pf4_exception(
    trace_c2: Fraction | int,
    trace_cd: Fraction | int,
    trace_d2: Fraction | int,
    gamma: Fraction | int,
    *,
    commuting: bool = False,
    symmetry_forces_mixed_zero: bool = False,
    exactly_correctable: bool = False,
) -> ObstructionException:
    """Classify a PF4 form without inferring exact correction from cancellation.

    ``exactly_correctable`` is an externally established premise.  A vanishing
    scalar form alone is never promoted to exact correctability.
    """

    flags = (commuting, symmetry_forces_mixed_zero, exactly_correctable)
    if not all(isinstance(flag, bool) for flag in flags):
        raise TypeError("classification flags must be bool")
    if sum(flags) > 1:
        raise ValueError("classification premises are mutually exclusive")
    c2 = _exact_fraction(trace_c2, "trace_c2")
    cd = _exact_fraction(trace_cd, "trace_cd")
    d2 = _exact_fraction(trace_d2, "trace_d2")
    scale = _exact_fraction(gamma, "gamma")
    value = pf4_trace_quadratic_form(c2, cd, d2, scale)

    if commuting:
        if any(moment != 0 for moment in (c2, cd, d2)):
            raise ValueError("commuting classification requires zero moments")
        return ObstructionException.COMMUTING
    if exactly_correctable:
        if value != 0:
            raise ValueError("a nonzero trace form cannot be exactly correctable")
        return ObstructionException.EXACTLY_CORRECTABLE
    if symmetry_forces_mixed_zero:
        if cd != 0:
            raise ValueError("the mixed moment is nonzero despite the symmetry premise")
        if scale <= 0 or c2 < 0 or d2 < 0 or value <= 0:
            raise ValueError("the supplied moments do not establish positivity")
        return ObstructionException.SYMMETRY_PROTECTED_POSITIVE
    return ObstructionException.MIXED_TERM_INDEFINITE


def _record_digest_payload(record: TraceIdentityRecord) -> bytes:
    payload = {
        "identity_kind": record.identity_kind,
        "formula": record.formula,
        "assumptions": list(record.assumptions),
        "exact_coefficients": [
            [name, str(value)] for name, value in record.exact_coefficients
        ],
        "exception_conditions": [
            [name, condition] for name, condition in record.exception_conditions
        ],
        "identity_status": record.identity_status,
    }
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def _build_record(
    *,
    identity_kind: str,
    formula: str,
    assumptions: tuple[str, ...],
    exact_coefficients: tuple[tuple[str, Fraction], ...],
    identity_status: str,
) -> TraceIdentityRecord:
    conditions = (
        (ObstructionException.COMMUTING.value, "all defining commutators vanish"),
        (
            ObstructionException.SYMMETRY_PROTECTED_POSITIVE.value,
            "symmetry kills the mixed moment and the diagonal form is positive",
        ),
        (
            ObstructionException.MIXED_TERM_INDEFINITE.value,
            "the mixed moment can cancel the diagonal moments",
        ),
        (
            ObstructionException.EXACTLY_CORRECTABLE.value,
            "an external operator identity proves exact correction",
        ),
    )
    unsigned = TraceIdentityRecord(
        identity_kind=identity_kind,
        formula=formula,
        assumptions=assumptions,
        exact_coefficients=exact_coefficients,
        exception_conditions=conditions,
        identity_status=identity_status,
        identity_digest="",
    )
    return TraceIdentityRecord(
        identity_kind=unsigned.identity_kind,
        formula=unsigned.formula,
        assumptions=unsigned.assumptions,
        exact_coefficients=unsigned.exact_coefficients,
        exception_conditions=unsigned.exception_conditions,
        identity_status=unsigned.identity_status,
        identity_digest=sha256(_record_digest_payload(unsigned)).hexdigest(),
    )


def pf2_trace_identity_record() -> TraceIdentityRecord:
    """Return the digest-bound record for the proved PF2 identity."""

    return _build_record(
        identity_kind="pf2",
        formula=(
            "Tr((A+B)L3) = (1/24) "
            "Tr([A,B]^dagger [A,B])"
        ),
        assumptions=(
            "A and B are finite-dimensional Hermitian matrices",
            (
                "L3 is the degree-three term in "
                "log(exp(tA/2) exp(tB) exp(tA/2))"
            ),
            "trace is cyclic",
        ),
        exact_coefficients=(("commutator_norm_square", Fraction(1, 24)),),
        identity_status="proved_exact_free_trace_identity",
    )


def pf4_trace_identity_record() -> TraceIdentityRecord:
    """Return a record that explicitly limits the status of the PF4 form."""

    return _build_record(
        identity_kind="pf4",
        formula=(
            "gamma*((1/2)Tr(C^2) + (14/3)Tr(CD) "
            "+ (4/3)Tr(D^2))"
        ),
        assumptions=(
            "trace moments are supplied independently as exact rationals",
            "the plan-stated coefficients were not independently derived from a BCH mapping",
            "no PF4 impossibility conclusion follows from this algebraic record alone",
        ),
        exact_coefficients=(
            ("trace_c2", Fraction(1, 2)),
            ("trace_cd", Fraction(14, 3)),
            ("trace_d2", Fraction(4, 3)),
        ),
        identity_status="algebraic_form_only_unverified_bch_mapping",
    )


def verify_identity_record(
    record: TraceIdentityRecord,
    expected_kind: str,
) -> bool:
    """Match a record to the pinned canonical statement for ``expected_kind``.

    Recomputing a digest after changing the formula or status is insufficient:
    every semantic field must equal the corresponding factory-owned record.
    """

    if not isinstance(record, TraceIdentityRecord):
        return False
    if not isinstance(expected_kind, str):
        return False
    factories = {
        "pf2": pf2_trace_identity_record,
        "pf4": pf4_trace_identity_record,
    }
    factory = factories.get(expected_kind)
    if factory is None or record.identity_kind != expected_kind:
        return False
    canonical = factory()
    if record != canonical:
        return False
    return sha256(_record_digest_payload(record)).hexdigest() == record.identity_digest


def verify_identity_digest(record: TraceIdentityRecord) -> bool:
    """Verify both the digest and the canonical record selected by its kind."""

    if not isinstance(record, TraceIdentityRecord):
        return False
    if not isinstance(record.identity_kind, str):
        return False
    return verify_identity_record(record, record.identity_kind)
