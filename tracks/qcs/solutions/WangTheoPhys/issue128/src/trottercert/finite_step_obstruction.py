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
