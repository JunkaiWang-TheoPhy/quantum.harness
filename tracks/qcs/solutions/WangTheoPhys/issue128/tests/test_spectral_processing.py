from fractions import Fraction
from math import asin

import pytest

from trottercert.spectral_processing import (
    spectral_phase_certificate,
    target_projector_loss,
)


def test_repeated_spectral_phase_certificate() -> None:
    result = spectral_phase_certificate(Fraction(1, 1000), 35)
    assert result.repeated_error == Fraction(35, 1000)
    assert result.phase_radius == pytest.approx(2 * asin(0.035 / 2))


def test_phase_certificate_and_projector_loss_fail_closed() -> None:
    with pytest.raises(ValueError, match="below two"):
        spectral_phase_certificate(Fraction(1, 10), 20)
    assert target_projector_loss(
        Fraction(1, 100),
        Fraction(1, 1000),
        Fraction(1, 5),
    ) == Fraction(3, 100)
    with pytest.raises(ValueError, match="positive gap"):
        target_projector_loss(Fraction(), Fraction(), Fraction())


def test_spectral_bookkeeping_rejects_negative_or_boolean_inputs() -> None:
    with pytest.raises(ValueError, match="nonnegative"):
        spectral_phase_certificate(Fraction(-1, 100), 2)
    with pytest.raises(ValueError, match="positive integer"):
        spectral_phase_certificate(Fraction(), True)
    with pytest.raises(ValueError, match="nonnegative"):
        target_projector_loss(Fraction(-1), Fraction(), Fraction(1))
    with pytest.raises(ValueError, match="exact rational"):
        spectral_phase_certificate(0.001, 35)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="exact rational"):
        target_projector_loss(0.01, Fraction(), Fraction(1))  # type: ignore[arg-type]
