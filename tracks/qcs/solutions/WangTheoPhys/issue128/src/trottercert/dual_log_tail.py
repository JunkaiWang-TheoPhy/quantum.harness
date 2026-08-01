from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any

from .intervals import (
    RationalInterval,
    cube_root_four_interval,
    outward_quantize,
)

LOG_RADIUS_CAP = Fraction(3, 2)
DEXP_CAP = Fraction(8, 5)
MULTIPLIER_DIFFERENCE_CAP = Fraction(1)
EVEN_SERIES_CAP = Fraction(10, 33)
REFERENCE_CERTIFICATE = (
    Path(__file__).resolve().parents[2]
    / "certificates/issue128-d5-integrated-certificate.json"
)


@dataclass(frozen=True, slots=True)
class DualIntervalStage:
    fragment_index: int
    coefficient: RationalInterval


@dataclass(frozen=True, slots=True)
class DualGeneratorConstants:
    stages: tuple[DualIntervalStage, ...]
    d4_site: Fraction
    d5_site: Fraction
    d6_site: Fraction
    d7_site: Fraction
    source_steps: int
    coefficient_interval_decimal_digits: int
    source_certificate: Path


@dataclass(frozen=True, slots=True)
class DualTailGeometry:
    n_sites: int
    cells: int
    steps: int
    centered_h_op_cap: Fraction
    h_hs_squared: Fraction
    h_hs_cap: Fraction
    w_hs_squared: Fraction
    w_hs_per_cell_cap: Fraction
    max_w_pauli_coefficient: Fraction
    log_radius_cap: Fraction
    dexp_cap: Fraction
    multiplier_difference_cap: Fraction
    even_series_cap: Fraction


@dataclass(frozen=True, slots=True)
class DualTailEnvelope:
    average_generator_defect: Fraction
    one_step_unitary_defect: Fraction
    pointwise_generator_defect: Fraction
    relative_log_defect: Fraction
    centered_exact_phase_radius: Fraction
    centered_log_radius_bound: Fraction
    log_defect: Fraction
    log_derivative_defect_hs: Fraction
    log_commutator_defect_hs: Fraction


@dataclass(frozen=True, slots=True)
class DualLogTailBound:
    geometry: DualTailGeometry
    envelope: DualTailEnvelope
    first_omitted_generator_degree: int
    first_omitted_log_degree: int
    direct_even_generator_tail: Fraction
    dexp_correction_tail: Fraction
    total_log_tail: Fraction
    claim_scope: str


def issue128_dual_tail_geometry() -> DualTailGeometry:
    return DualTailGeometry(
        n_sites=144,
        cells=36,
        steps=97,
        centered_h_op_cap=Fraction(144),
        h_hs_squared=Fraction(54),
        h_hs_cap=Fraction(15, 2),
        w_hs_squared=Fraction(23085, 4),
        w_hs_per_cell_cap=Fraction(17, 8),
        max_w_pauli_coefficient=Fraction(1, 4),
        log_radius_cap=LOG_RADIUS_CAP,
        dexp_cap=DEXP_CAP,
        multiplier_difference_cap=MULTIPLIER_DIFFERENCE_CAP,
        even_series_cap=EVEN_SERIES_CAP,
    )


def _second_order_interval_stages(
    scale: RationalInterval,
) -> list[DualIntervalStage]:
    half = scale / 2
    stages = [DualIntervalStage(index, half) for index in range(3)]
    stages.append(DualIntervalStage(3, scale))
    stages.extend(
        DualIntervalStage(index, half) for index in reversed(range(3))
    )
    return stages


def _merge_interval_stages(
    stages: Sequence[DualIntervalStage],
) -> tuple[DualIntervalStage, ...]:
    merged: list[DualIntervalStage] = []
    for stage in stages:
        if merged and merged[-1].fragment_index == stage.fragment_index:
            previous = merged.pop()
            merged.append(
                DualIntervalStage(
                    stage.fragment_index,
                    previous.coefficient + stage.coefficient,
                )
            )
        else:
            merged.append(stage)
    return tuple(merged)


def issue128_interval_stages(
    decimal_digits: int,
) -> tuple[DualIntervalStage, ...]:
    """Independently rebuild the frozen four-fragment Suzuki stages."""

    if (
        not isinstance(decimal_digits, int)
        or isinstance(decimal_digits, bool)
        or decimal_digits < 1
    ):
        raise ValueError("stage decimal digits must be a positive integer")
    root = cube_root_four_interval(decimal_digits)
    grid = 10**decimal_digits
    u = outward_quantize(RationalInterval.point(1) / (4 - root), grid)
    scales = (u, u, 1 - 4 * u, u, u)
    stages: list[DualIntervalStage] = []
    for scale in scales:
        stages.extend(_second_order_interval_stages(scale))
    result = _merge_interval_stages(stages)
    if len(result) != 31:
        raise ArithmeticError("unexpected merged Suzuki stage count")
    return result


def _parse_pair(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(
            not isinstance(entry, int) or isinstance(entry, bool)
            for entry in value
        )
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if [result.numerator, result.denominator] != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def generator_constants_from_certificate(
    payload: Mapping[str, Any],
    *,
    source_certificate: Path,
) -> DualGeneratorConstants:
    """Recover valid D4--D7 density caps from a verified main certificate."""

    if payload.get("schema_version") != 3:
        raise ValueError("reference certificate schema mismatch")
    candidate = payload.get("candidate")
    if not isinstance(candidate, Mapping):
        raise ValueError("reference certificate candidate is missing")
    steps = candidate.get("steps")
    decimal_digits = candidate.get("coefficient_interval_decimal_digits")
    if steps != 95:
        raise ValueError("reference certificate step count mismatch")
    if decimal_digits != 12:
        raise ValueError("reference coefficient precision mismatch")
    contributions = candidate.get("contributions")
    if not isinstance(contributions, Mapping):
        raise ValueError("reference contributions are missing")
    densities: dict[int, Fraction] = {}
    for degree in range(4, 8):
        contribution = _parse_pair(
            contributions.get(f"degree{degree}"),
            f"degree-{degree} contribution",
        )
        if contribution <= 0:
            raise ValueError("reference contribution must be positive")
        densities[degree] = (
            contribution * (degree + 1) * steps**degree / 144
        )
    return DualGeneratorConstants(
        stages=issue128_interval_stages(decimal_digits),
        d4_site=densities[4],
        d5_site=densities[5],
        d6_site=densities[6],
        d7_site=densities[7],
        source_steps=steps,
        coefficient_interval_decimal_digits=decimal_digits,
        source_certificate=source_certificate,
    )


@lru_cache(maxsize=1)
def issue128_generator_constants() -> DualGeneratorConstants:
    try:
        payload = json.loads(REFERENCE_CERTIFICATE.read_bytes())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("cannot read the reference certificate") from exc
    if not isinstance(payload, dict):
        raise ValueError("reference certificate must contain a JSON object")
    return generator_constants_from_certificate(
        payload,
        source_certificate=REFERENCE_CERTIFICATE,
    )


def _validate_tail_inputs(
    steps: int,
    first_omitted_degree: int,
) -> None:
    if not isinstance(steps, int) or isinstance(steps, bool) or steps < 1:
        raise ValueError("steps must be a positive integer")
    if (
        not isinstance(first_omitted_degree, int)
        or isinstance(first_omitted_degree, bool)
        or first_omitted_degree < 0
    ):
        raise ValueError("first omitted degree must be nonnegative")


def average_stage_tail(
    stages: Sequence[DualIntervalStage],
    steps: int,
    first_omitted_degree: int,
) -> Fraction:
    """Average-in-time absolute conjugation tail without its density factor."""

    _validate_tail_inputs(steps, first_omitted_degree)
    prefix = Fraction()
    total = Fraction()
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("stage tail is outside its convergence region")
        total += coefficient * ratio**first_omitted_degree / (1 - ratio)
        prefix += coefficient
    return total


def pointwise_stage_tail(
    stages: Sequence[DualIntervalStage],
    steps: int,
    first_omitted_degree: int,
) -> Fraction:
    """Endpoint absolute conjugation tail without its density factor."""

    _validate_tail_inputs(steps, first_omitted_degree)
    prefix = Fraction()
    total = Fraction()
    degree = first_omitted_degree
    for stage in stages:
        coefficient = stage.coefficient.abs_upper()
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("stage tail is outside its convergence region")
        total += (
            coefficient
            * ratio**degree
            * ((degree + 1) - degree * ratio)
            / (1 - ratio) ** 2
        )
        prefix += coefficient
    return total


def verify_analytic_caps(
    *,
    log_radius_cap: Fraction,
    dexp_cap: Fraction,
    multiplier_difference_cap: Fraction,
    even_series_cap: Fraction,
) -> None:
    """Verify the rational scalar lemmas used by the ``dexp`` conversion."""

    radius = Fraction(log_radius_cap)
    if radius < 0 or radius > LOG_RADIUS_CAP:
        raise ValueError("log radius exceeds the proved interval")

    # sin(x)/x >= 1-x^2/6 on [0, 3/2].
    required_dexp = 1 / (1 - radius**2 / 6)
    if Fraction(dexp_cap) < required_dexp:
        raise ValueError("dexp cap is below the rational sine bound")

    # For y=2x and |x|<=3/2, the cotangent remainder obeys
    # |x*cot(x)-1|/|x| <= x/(3*(1-x^2/6)) <= 4/5.  Therefore
    # |(f(iy)-1)/(iy)|^2 <= 1/4+(2/5)^2 < 1.
    if Fraction(multiplier_difference_cap) < MULTIPLIER_DIFFERENCE_CAP:
        raise ValueError(
            "multiplier-difference cap is below the proved cap"
        )

    # The k=1 term is 1/3! and every later ratio is at most 9/20.
    required_even_series = Fraction(1, 6) / (1 - Fraction(9, 20))
    if Fraction(even_series_cap) < required_even_series:
        raise ValueError("even-series cap is below the ratio majorant")


def _verify_mathematical_geometry(geometry: DualTailGeometry) -> None:
    if geometry.n_sites != 144:
        raise ValueError("fixed-instance site count mismatch")
    if geometry.cells != 36:
        raise ValueError("fixed-instance cell count mismatch")
    if geometry.steps < 97:
        raise ValueError("step count lies below the proved convergence range")
    if geometry.centered_h_op_cap != 144:
        raise ValueError("centered Hamiltonian cap mismatch")
    if geometry.h_hs_squared != 54:
        raise ValueError("H moment mismatch")
    if geometry.h_hs_cap != Fraction(15, 2):
        raise ValueError("H Hilbert--Schmidt cap mismatch")
    if geometry.h_hs_cap**2 < geometry.h_hs_squared:
        raise ValueError("H Hilbert--Schmidt cap is invalid")
    if geometry.w_hs_squared != Fraction(23085, 4):
        raise ValueError("W moment mismatch")
    if geometry.w_hs_per_cell_cap != Fraction(17, 8):
        raise ValueError("W Hilbert--Schmidt cap mismatch")
    if (
        geometry.w_hs_per_cell_cap**2
        < geometry.w_hs_squared / geometry.cells**2
    ):
        raise ValueError("W Hilbert--Schmidt cap is invalid")
    if geometry.max_w_pauli_coefficient != Fraction(1, 4):
        raise ValueError("W Pauli coefficient cap mismatch")
    if geometry.log_radius_cap != LOG_RADIUS_CAP:
        raise ValueError("log-radius cap mismatch")
    if geometry.dexp_cap != DEXP_CAP:
        raise ValueError("dexp cap mismatch")
    if geometry.multiplier_difference_cap != MULTIPLIER_DIFFERENCE_CAP:
        raise ValueError("multiplier-difference cap mismatch")
    if geometry.even_series_cap != EVEN_SERIES_CAP:
        raise ValueError("even-series cap mismatch")


def _verify_fixed_geometry(geometry: DualTailGeometry) -> None:
    _verify_mathematical_geometry(geometry)
    if geometry.steps != 97:
        raise ValueError("fixed-instance step count mismatch")


def build_dual_tail_envelope(
    constants: DualGeneratorConstants,
    geometry: DualTailGeometry,
) -> DualTailEnvelope:
    """Build the fixed-instance right-generator-to-log envelope."""

    _verify_mathematical_geometry(geometry)
    verify_analytic_caps(
        log_radius_cap=geometry.log_radius_cap,
        dexp_cap=geometry.dexp_cap,
        multiplier_difference_cap=geometry.multiplier_difference_cap,
        even_series_cap=geometry.even_series_cap,
    )
    t = Fraction(1, geometry.steps)
    n_sites = Fraction(geometry.n_sites)
    direct_average = (
        constants.d4_site * t**4 / 5
        + constants.d5_site * t**5 / 6
        + constants.d6_site * t**6 / 7
        + constants.d7_site * t**7 / 8
    )
    average_defect = n_sites * (
        direct_average
        + Fraction(3, 8)
        * average_stage_tail(constants.stages, geometry.steps, 8)
    )
    pointwise_defect = n_sites * (
        constants.d4_site * t**4
        + constants.d5_site * t**5
        + constants.d6_site * t**6
        + constants.d7_site * t**7
        + Fraction(3, 8)
        * pointwise_stage_tail(constants.stages, geometry.steps, 8)
    )
    one_step = t * average_defect
    if one_step >= 1:
        raise ArithmeticError("unitary defect is outside the logarithm gate")
    relative_log = one_step / (1 - one_step)
    exact_radius = t * geometry.centered_h_op_cap
    centered_radius = exact_radius + relative_log
    if centered_radius > geometry.log_radius_cap:
        raise ArithmeticError("centered principal-log radius gate failed")
    log_defect = geometry.dexp_cap * relative_log
    log_derivative_defect = (
        geometry.dexp_cap * pointwise_defect
        + 2
        * geometry.multiplier_difference_cap
        * log_defect
        * geometry.h_hs_cap
    )
    commutator_defect = 2 * (
        exact_radius * log_derivative_defect
        + log_defect * geometry.h_hs_cap
        + log_defect * log_derivative_defect
    )
    return DualTailEnvelope(
        average_generator_defect=average_defect,
        one_step_unitary_defect=one_step,
        pointwise_generator_defect=pointwise_defect,
        relative_log_defect=relative_log,
        centered_exact_phase_radius=exact_radius,
        centered_log_radius_bound=centered_radius,
        log_defect=log_defect,
        log_derivative_defect_hs=log_derivative_defect,
        log_commutator_defect_hs=commutator_defect,
    )


def certify_issue128_dual_log_tail(
    constants: DualGeneratorConstants | None = None,
) -> DualLogTailBound:
    """Certify the compatible E11-and-higher dual local-log tail."""

    if constants is None:
        constants = issue128_generator_constants()
    geometry = issue128_dual_tail_geometry()
    _verify_fixed_geometry(geometry)
    envelope = build_dual_tail_envelope(constants, geometry)
    direct = Fraction(3, 8) * average_stage_tail(
        constants.stages,
        geometry.steps,
        first_omitted_degree=10,
    )
    correction = (
        geometry.w_hs_per_cell_cap
        * 2
        * envelope.log_defect
        * envelope.log_commutator_defect_hs
        * geometry.even_series_cap
        / 11
    )
    total = direct + correction
    if total <= 0:
        raise ArithmeticError("dual local-log tail must be positive")
    return DualLogTailBound(
        geometry=geometry,
        envelope=envelope,
        first_omitted_generator_degree=10,
        first_omitted_log_degree=11,
        direct_even_generator_tail=direct,
        dexp_correction_tail=correction,
        total_log_tail=total,
        claim_scope="fixed_12x12_r97_dual_local_log",
    )
