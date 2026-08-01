import copy
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from trottercert.commutant_witness import (
    quadratic_witness_record,
    verify_quadratic_witness_payload,
)
from trottercert.cubic_field import Cubic
from trottercert.spectral_duality import (
    check_dual_witness,
    energy_blocks,
    gauge_residual,
)
from trottercert.spectral_gauge import check_exact_dual_witness

ROOT = Path(__file__).resolve().parents[1]
STORED_WITNESS = (
    ROOT
    / "docs/experiments/processor-obstruction/quadratic-commutant-witness.json"
)


def test_energy_blocks_preserve_degenerate_eigenspaces() -> None:
    h = np.diag([0.0, 0.0, 2.0])
    blocks = energy_blocks(h, atol=1e-12)
    assert tuple(block.shape for block in blocks) == ((3, 2), (3, 1))
    pinching_identity = sum(block @ block.conj().T for block in blocks)
    np.testing.assert_allclose(pinching_identity, np.eye(3), atol=1e-14)


def test_commutator_has_zero_energy_block_diagonal() -> None:
    h = np.diag([0.0, 0.0, 2.0])
    q = np.array(
        [[0, 1j, 2j], [-1j, 0, 3j], [-2j, -3j, 0]],
        dtype=complex,
    )
    comm = 1j * (q @ h - h @ q)
    residual = gauge_residual(h, comm)
    assert np.linalg.norm(residual, ord=2) < 1e-12


def test_degenerate_block_component_survives_commutator_quotient() -> None:
    h = np.diag([0.0, 0.0, 2.0])
    defect = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    np.testing.assert_allclose(gauge_residual(h, defect), defect, atol=1e-12)


def test_retiming_and_phase_are_removed() -> None:
    h = np.diag([-1.0, 2.0])
    e = 3.0 * np.eye(2) - 0.25 * h
    assert np.linalg.norm(gauge_residual(h, e), ord=2) < 1e-12


def test_singular_identity_retiming_span_is_handled() -> None:
    h = 2.0 * np.eye(3)
    e = -7.0 * np.eye(3)
    assert np.linalg.norm(gauge_residual(h, e), ord=2) < 1e-12


def test_nearly_proportional_identity_and_retiming_are_projected_stably() -> None:
    h = np.eye(2) + 1e-8 * np.diag([-1.0, 1.0])
    assert np.linalg.norm(gauge_residual(h, h), ord=2) < 1e-12


def test_tuned_qubit_is_gauge_correctable() -> None:
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.diag([1.0, -1.0]).astype(complex)
    processor = y / 2
    defect = 1j * (processor @ z - z @ processor) + 3 * np.eye(2) - z / 4
    np.testing.assert_allclose(defect, -x + 3 * np.eye(2) - z / 4, atol=1e-14)
    assert np.linalg.norm(gauge_residual(z, defect), ord=2) < 1e-12


def test_dual_witness_check_and_weak_duality() -> None:
    h = np.diag([-1.0, 0.0, 2.0])
    e = np.diag([2.0, -3.0, 1.0])
    w = np.diag([2.0, -3.0, 1.0])
    # Make the witness orthogonal to both I and H without changing its
    # commutant status.
    gram = np.array(
        [[np.trace(np.eye(3)), np.trace(h)], [np.trace(h), np.trace(h @ h)]],
        dtype=float,
    )
    rhs = np.array([np.trace(w), np.trace(w @ h)], dtype=float)
    phase, retiming = np.linalg.solve(gram, rhs)
    w = w - phase * np.eye(3) - retiming * h

    result = check_dual_witness(h, e, w, atol=1e-12)
    assert result["commutes"] is True
    assert result["identity_orthogonal"] is True
    assert result["retiming_orthogonal"] is True
    assert result["nonzero"] is True
    assert result["valid"] is True
    assert result["trace_norm"] > 0
    assert abs(result["pairing"]) / result["trace_norm"] <= (
        np.linalg.norm(gauge_residual(h, e), ord=2) + 1e-12
    )


def test_simple_diagonal_witness_saturates_spectral_primal() -> None:
    h = np.diag([-1.0, 0.0, 1.0])
    e = np.diag([1.0, -2.0, 1.0])
    result = check_dual_witness(h, e, e, atol=1e-12)
    dual_value = abs(result["pairing"]) / result["trace_norm"]
    # The optimal phase shift is +I/2, giving diagonal entries +/-3/2.
    primal_value = np.linalg.norm(e + 0.5 * np.eye(3), ord=2)
    assert abs(dual_value - primal_value) < 1e-12


def test_zero_witness_is_not_valid_or_normalizable() -> None:
    h = np.diag([-1.0, 1.0])
    result = check_dual_witness(h, h, np.zeros((2, 2)), atol=1e-12)
    assert result["nonzero"] is False
    assert result["valid"] is False
    assert result["rigorous"] is False
    assert result["trace_norm"] == 0
    assert "lower_bound" not in result
    assert "normalized_pairing" not in result


def test_large_gauge_only_scale_does_not_create_a_quotient_pairing() -> None:
    h = np.diag([-1.0, 0.0, 1.0])
    gauge_only = 1e12 * np.eye(3) - 3e11 * h
    # This witness passes loose discovery tolerances but is not exactly
    # identity-orthogonal.  Pairing against raw E would spuriously be O(1).
    witness = np.diag([1.0, -2.0, 1.0]) + 2.5e-13 * np.eye(3)
    result = check_dual_witness(h, gauge_only, witness, atol=1e-12)
    assert result["valid"] is True
    assert result["rigorous"] is False
    assert abs(np.trace(witness @ gauge_only)) > 0.1
    assert abs(result["pairing"]) < 1e-3
    assert "lower_bound" not in result


def test_exact_dual_witness_check_handles_degenerate_blocks() -> None:
    h = sp.diag(0, 0, 2)
    e = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    result = check_exact_dual_witness(h, e, e)
    assert result.commutes is True
    assert result.identity_orthogonal is True
    assert result.retiming_orthogonal is True
    assert result.pairing == 2


def test_heisenberg_witness_record_has_exact_nonzero_calibrated_pairing() -> None:
    payload = json.loads(STORED_WITNESS.read_text())
    record = quadratic_witness_record(payload)
    assert record.identity_pairing == 0
    assert record.hamiltonian_pairing == 0
    assert record.defect_pairing != Cubic.zero()
    assert record.h2_coefficient == 1
    assert record.exact_nonzero_obstruction is True
    assert record.normalization_status == "exact_norm_certificate_required"
    assert record.defect_pairing == Cubic(
        Fraction(-7807, 100000),
        Fraction(-66043, 1600000),
        Fraction(-6119, 200000),
    )


def test_quadratic_witness_rejects_mutated_h2_coefficient() -> None:
    payload = json.loads(STORED_WITNESS.read_text())
    forged = copy.deepcopy(payload)
    forged["witness"]["h2_coefficient"] = [2, 1]
    with pytest.raises(ValueError, match="H-squared coefficient"):
        verify_quadratic_witness_payload(forged)
    with pytest.raises(ValueError, match="H-squared coefficient"):
        quadratic_witness_record(forged)
