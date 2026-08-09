from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import asin


@dataclass(frozen=True)
class SpectralPhaseCertificate:
    one_step_error: Fraction
    repetitions: int
    repeated_error: Fraction
    phase_radius: float


def spectral_phase_certificate(
    one_step_error: Fraction | int,
    repetitions: int,
) -> SpectralPhaseCertificate:
    """Convert a processed one-step norm bound to a repeated phase radius.

    If ``V = P Sigma P^dagger`` and ``||V-U|| <= delta``, telescoping gives
    ``||V^r-U^r|| <= r delta``.  Unitary spectral variation then places each
    eigenphase within ``2 asin(r delta / 2)`` of the other phase set.  The
    conversion is rejected when the norm bound reaches the vacuous value two.
    """

    if isinstance(one_step_error, bool) or not isinstance(
        one_step_error,
        (Fraction, int),
    ):
        raise ValueError("one-step error must be an exact rational value")
    error = Fraction(one_step_error)
    if error < 0:
        raise ValueError("one-step error must be nonnegative")
    if isinstance(repetitions, bool) or not isinstance(repetitions, int):
        raise ValueError("repetitions must be a positive integer")
    if repetitions < 1:
        raise ValueError("repetitions must be a positive integer")
    repeated = repetitions * error
    if repeated >= 2:
        raise ValueError("repeated unitary error must remain below two")
    return SpectralPhaseCertificate(
        one_step_error=error,
        repetitions=repetitions,
        repeated_error=repeated,
        phase_radius=2 * asin(float(repeated) / 2),
    )


def target_projector_loss(
    processor_distance: Fraction | int,
    effective_hamiltonian_error: Fraction | int,
    gap: Fraction | int,
) -> Fraction:
    """Conservative loss budget for a processed target eigenspace.

    The two terms track the near-identity processor rotation and a standard
    gap-dependent spectral-projector perturbation bound.  Values at or above
    one are deliberately retained so a caller can identify a vacuous claim.
    """

    if any(
        isinstance(value, bool) or not isinstance(value, (Fraction, int))
        for value in (processor_distance, effective_hamiltonian_error, gap)
    ):
        raise ValueError("projector inputs must be exact rational values")
    distance = Fraction(processor_distance)
    hamiltonian_error = Fraction(effective_hamiltonian_error)
    spectral_gap = Fraction(gap)
    if distance < 0 or hamiltonian_error < 0:
        raise ValueError("projector error inputs must be nonnegative")
    if spectral_gap <= 0:
        raise ValueError("target projector bound requires a positive gap")
    return 2 * distance + 2 * hamiltonian_error / spectral_gap
