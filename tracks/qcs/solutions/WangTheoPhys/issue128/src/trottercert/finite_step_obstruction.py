"""Exact signed-margin arithmetic for the Issue-128 finite-step obstruction.

The trusted input is a dual local-log expansion per lattice cell,

    h**4 q5 + h**6 q7 + h**8 q9 + R_{>= 11},

where the three coefficients lie in ``Q[cuberoot(4)]`` and the final input is
an outward rational enclosure of the absolute tail bound.  This module makes
no claim about where those inputs came from; artifact verifiers are
responsible for establishing their provenance before calling the decision
routine.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .cubic_field import Cubic
from .intervals import RationalInterval


FINITE_STEP_STATUSES = (
    "certified_obstruction",
    "certified_no_margin",
    "inconclusive",
)
AFFINE_SPECTRAL_STATUSES = (
    "certified_affine_effective_spectral_obstruction",
    "inconclusive",
)


@dataclass(frozen=True, slots=True)
class FiniteStepDecision:
    """An exact interval decision for the finite-step dual pairing."""

    h: Fraction
    root_interval: RationalInterval
    leading_interval: RationalInterval
    absolute_leading_interval: RationalInterval
    tail_interval: RationalInterval
    signed_margin_interval: RationalInterval
    status: str

    def __post_init__(self) -> None:
        if self.status not in FINITE_STEP_STATUSES:
            raise ValueError("unknown finite-step decision status")


@dataclass(frozen=True, slots=True)
class AffineSpectralDecision:
    """A nonlinear centered-moment obstruction to the affine unitary orbit."""

    centered_defect_cap: Fraction
    linear_invariant_lower: Fraction
    nonlinear_remainder_upper: Fraction
    invariant_margin: Fraction
    status: str

    def __post_init__(self) -> None:
        if self.status not in AFFINE_SPECTRAL_STATUSES:
            raise ValueError("unknown affine-spectral decision status")


def _validate_positive_fraction(value: Fraction, *, field: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{field} must be a positive rational")
    result = Fraction(value)
    if result <= 0:
        raise ValueError(f"{field} must be positive")
    return result


def _absolute_interval(value: RationalInterval) -> RationalInterval:
    return RationalInterval(value.abs_lower(), value.abs_upper())


def decide_finite_step(
    q5: Cubic,
    q7: Cubic,
    q9: Cubic,
    root: RationalInterval,
    tail: RationalInterval,
    h: Fraction,
) -> FiniteStepDecision:
    """Combine exact pairings and a tail bound without a floating conversion.

    ``tail`` encloses the certified nonnegative upper bound on
    ``|R_{>=11}^{dual}|``.  A claim is promoted only when every allowed root
    and tail value gives a strictly positive margin.  A strictly negative
    interval records that this particular bound cannot establish a margin;
    it is not a proof that the physical obstruction vanishes.
    """

    if not all(isinstance(value, Cubic) for value in (q5, q7, q9)):
        raise TypeError("q5, q7, and q9 must be exact Cubic values")
    step = _validate_positive_fraction(h, field="h")
    if root.lower <= 0 or root.upper**3 < 4 or root.lower**3 > 4:
        raise ValueError("root interval must enclose the positive cube root of four")
    if tail.lower < 0:
        raise ValueError("tail interval must be nonnegative")

    leading = (
        q5.enclose(root) * step**4
        + q7.enclose(root) * step**6
        + q9.enclose(root) * step**8
    )
    absolute_leading = _absolute_interval(leading)
    margin = absolute_leading - tail
    if margin.lower > 0:
        status = "certified_obstruction"
    elif margin.upper < 0:
        status = "certified_no_margin"
    else:
        status = "inconclusive"
    return FiniteStepDecision(
        h=step,
        root_interval=root,
        leading_interval=leading,
        absolute_leading_interval=absolute_leading,
        tail_interval=tail,
        signed_margin_interval=margin,
        status=status,
    )


def decide_affine_spectral_obstruction(
    *,
    dual_margin_lower: Fraction,
    cells: int,
    m2: Fraction,
    m3: Fraction,
    h_hs_cap: Fraction,
    effective_log_defect: Fraction,
) -> AffineSpectralDecision:
    """Bound a centered-moment invariant away from the affine unitary orbit.

    The invariant is

    ``m2(H)^3 m3(A)^2 - m3(H)^2 m2(A)^3``.

    It vanishes whenever ``A = a I + b U H U^dagger``.  Its exact linear term
    at ``A=H+D`` is ``6*m3*m2^3*tau(WD)`` for the Issue-128 witness convention.
    The returned remainder uses only normalized Schatten inequalities and an
    operator-norm cap on the effective logarithm defect.
    """

    if isinstance(cells, bool) or not isinstance(cells, int) or cells <= 0:
        raise ValueError("cells must be a positive integer")
    values = {
        "dual margin": dual_margin_lower,
        "m2": m2,
        "m3": m3,
        "Hamiltonian Hilbert--Schmidt cap": h_hs_cap,
        "effective log defect": effective_log_defect,
    }
    exact: dict[str, Fraction] = {}
    for field, value in values.items():
        if isinstance(value, bool):
            raise TypeError(f"{field} must be rational")
        exact[field] = Fraction(value)
    margin = exact["dual margin"]
    moment2 = exact["m2"]
    moment3 = exact["m3"]
    hs_cap = exact["Hamiltonian Hilbert--Schmidt cap"]
    log_defect = exact["effective log defect"]
    if margin < 0:
        raise ValueError("dual margin lower bound must be nonnegative")
    if moment2 <= 0:
        raise ValueError("m2 must be positive")
    if moment3 == 0:
        raise ValueError("m3 must be nonzero")
    if hs_cap <= 0 or hs_cap**2 < moment2:
        raise ValueError("Hamiltonian Hilbert--Schmidt cap is invalid")
    if log_defect < 0:
        raise ValueError("effective log defect must be nonnegative")

    epsilon = 2 * log_defect
    linear_d2 = 2 * hs_cap * epsilon
    nonlinear_d2 = epsilon**2
    absolute_d2 = linear_d2 + nonlinear_d2
    linear_d3 = 3 * moment2 * epsilon
    nonlinear_d3 = 3 * hs_cap * epsilon**2 + epsilon**3
    absolute_d3 = linear_d3 + nonlinear_d3
    remainder = moment2**3 * (
        absolute_d3**2 + 2 * abs(moment3) * nonlinear_d3
    ) + moment3**2 * (
        3 * moment2**2 * nonlinear_d2
        + 3 * moment2 * absolute_d2**2
        + absolute_d2**3
    )
    linear_lower = (
        6 * abs(moment3) * moment2**3 * cells * margin
    )
    invariant_margin = linear_lower - remainder
    status = (
        "certified_affine_effective_spectral_obstruction"
        if invariant_margin > 0
        else "inconclusive"
    )
    return AffineSpectralDecision(
        centered_defect_cap=epsilon,
        linear_invariant_lower=linear_lower,
        nonlinear_remainder_upper=remainder,
        invariant_margin=invariant_margin,
        status=status,
    )
