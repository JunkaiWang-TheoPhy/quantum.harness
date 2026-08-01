from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from functools import reduce
from math import gcd
from typing import Sequence, TypeAlias

import sympy as sp


ExactScalar: TypeAlias = int | Fraction | sp.Expr


class GaugeScope(str, Enum):
    FULL_COMMUTATOR_IMAGE = "full_commutator_image"
    RESTRICTED_PROCESSOR_SPAN = "restricted_processor_span"


@dataclass(frozen=True)
class GaugeDecomposition:
    target: tuple[sp.Expr, ...]
    generators: tuple[tuple[sp.Expr, ...], ...]
    independent_generators: tuple[tuple[sp.Expr, ...], ...]
    pivot_indices: tuple[int, ...]
    metric: tuple[sp.Expr, ...]
    coefficients: tuple[sp.Expr, ...]
    projection: tuple[sp.Expr, ...]
    residual: tuple[sp.Expr, ...]
    primitive_witness: tuple[int, ...] | None
    witness_target_pairing: sp.Expr
    scope: GaugeScope
    completeness_id: str | None
    status: str


def _exact_rational(value: object, *, field: str) -> sp.Rational:
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError(f"{field} must contain exact rational values")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    try:
        expression = sp.sympify(value)
    except (sp.SympifyError, TypeError) as exc:
        raise ValueError(f"{field} must contain exact rational values") from exc
    if expression.has(sp.Float) or expression.is_Rational is not True:
        if expression.is_real is False:
            raise ValueError(f"{field} values must be real")
        raise ValueError(f"{field} must contain exact rational values")
    if expression.is_real is not True:
        raise ValueError(f"{field} values must be real")
    return sp.Rational(expression)


def _exact_vector(values: Sequence[object], *, field: str) -> tuple[sp.Expr, ...]:
    return tuple(_exact_rational(value, field=field) for value in values)


def primitive_integer_vector(
    values: Sequence[sp.Expr],
) -> tuple[int, ...] | None:
    rationals = tuple(_exact_rational(value, field="witness") for value in values)
    if all(value == 0 for value in rationals):
        return None
    denominator_lcm = sp.ilcm(*(int(value.q) for value in rationals))
    integers = [int(value * denominator_lcm) for value in rationals]
    common = reduce(gcd, (abs(value) for value in integers if value), 0)
    normalized = [value // common for value in integers]
    first = next(value for value in normalized if value)
    if first < 0:
        normalized = [-value for value in normalized]
    return tuple(normalized)


def _as_tuple(column: sp.MatrixBase) -> tuple[sp.Expr, ...]:
    return tuple(sp.cancel(column[row, 0]) for row in range(column.rows))


def decompose_gauge(
    target: Sequence[ExactScalar],
    generators: Sequence[Sequence[ExactScalar]],
    *,
    metric: Sequence[ExactScalar] | None = None,
    scope: GaugeScope,
    completeness_id: str | None = None,
) -> GaugeDecomposition:
    """Project a target onto an exact gauge span.

    A nonzero residual is a spectral obstruction only when ``scope`` declares
    that the supplied generators span the full commutator image and the caller
    supplies a nonempty completeness identifier.  Otherwise it is only a
    restricted-processor obstruction.
    """

    if not target:
        raise ValueError("target must be nonempty")
    exact_target = _exact_vector(target, field="target")
    dimension = len(exact_target)
    exact_generators: list[tuple[sp.Expr, ...]] = []
    for generator in generators:
        exact_generator = _exact_vector(generator, field="generator")
        if len(exact_generator) != dimension:
            raise ValueError("generator dimension must match target dimension")
        exact_generators.append(exact_generator)

    if metric is None:
        exact_metric = tuple(sp.Integer(1) for _ in range(dimension))
    else:
        exact_metric = _exact_vector(metric, field="metric")
        if len(exact_metric) != dimension:
            raise ValueError("metric dimension must match target dimension")
        if any(value <= 0 for value in exact_metric):
            raise ValueError("metric weights must be positive")

    try:
        exact_scope = GaugeScope(scope)
    except ValueError as exc:
        raise ValueError("unknown gauge scope") from exc
    if exact_scope is GaugeScope.FULL_COMMUTATOR_IMAGE and not completeness_id:
        raise ValueError("full commutator image requires completeness evidence")
    if completeness_id is not None and (
        not isinstance(completeness_id, str) or not completeness_id.strip()
    ):
        raise ValueError("completeness identifier must be a nonempty string")

    target_column = sp.Matrix(exact_target)
    metric_matrix = sp.diag(*exact_metric)
    if exact_generators:
        generator_matrix = sp.Matrix.hstack(
            *(sp.Matrix(generator) for generator in exact_generators)
        )
        pivot_indices = tuple(int(index) for index in generator_matrix.rref()[1])
        independent_columns = tuple(
            generator_matrix[:, index] for index in pivot_indices
        )
    else:
        pivot_indices = ()
        independent_columns = ()

    if independent_columns:
        basis = sp.Matrix.hstack(*independent_columns)
        gram = basis.T * metric_matrix * basis
        rhs = basis.T * metric_matrix * target_column
        coefficient_column = gram.inv() * rhs
        projection_column = basis * coefficient_column
    else:
        coefficient_column = sp.zeros(0, 1)
        projection_column = sp.zeros(dimension, 1)
    residual_column = (target_column - projection_column).applyfunc(sp.cancel)

    projection = _as_tuple(projection_column)
    residual = _as_tuple(residual_column)
    primitive_witness = primitive_integer_vector(residual)
    if primitive_witness is None:
        status = "removable"
        witness_target_pairing = sp.Integer(0)
    else:
        status = (
            "spectral_obstruction"
            if exact_scope is GaugeScope.FULL_COMMUTATOR_IMAGE
            else "restricted_obstruction"
        )
        witness_column = sp.Matrix(primitive_witness)
        witness_target_pairing = sp.cancel(
            (witness_column.T * metric_matrix * target_column)[0, 0]
        )

    return GaugeDecomposition(
        target=exact_target,
        generators=tuple(exact_generators),
        independent_generators=tuple(_as_tuple(column) for column in independent_columns),
        pivot_indices=pivot_indices,
        metric=exact_metric,
        coefficients=_as_tuple(coefficient_column),
        projection=projection,
        residual=residual,
        primitive_witness=primitive_witness,
        witness_target_pairing=witness_target_pairing,
        scope=exact_scope,
        completeness_id=completeness_id,
        status=status,
    )


def verify_gauge_decomposition(result: GaugeDecomposition) -> None:
    """Recompute an exact decomposition and reject any changed field."""

    expected = decompose_gauge(
        result.target,
        result.generators,
        metric=result.metric,
        scope=result.scope,
        completeness_id=result.completeness_id,
    )
    if result != expected:
        raise ValueError("gauge decomposition status or data mismatch")
    if result.status == "removable":
        if result.primitive_witness is not None or any(result.residual):
            raise ValueError("removable status carries a nonzero witness")
        return

    if result.primitive_witness is None or result.witness_target_pairing == 0:
        raise ValueError("obstruction status requires a nonzero witness pairing")
    metric_matrix = sp.diag(*result.metric)
    witness = sp.Matrix(result.primitive_witness)
    for generator in result.independent_generators:
        pairing = (witness.T * metric_matrix * sp.Matrix(generator))[0, 0]
        if sp.cancel(pairing) != 0:
            raise ValueError("witness is not orthogonal to the gauge span")
    if (
        result.status == "spectral_obstruction"
        and result.scope is not GaugeScope.FULL_COMMUTATOR_IMAGE
    ):
        raise ValueError("spectral obstruction requires the full commutator image")
