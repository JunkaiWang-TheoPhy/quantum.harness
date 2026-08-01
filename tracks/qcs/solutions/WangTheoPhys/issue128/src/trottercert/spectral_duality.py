"""Numerical discovery helpers for gauge-aware spectral duality.

The routines in this module deliberately operate on dense floating-point
matrices.  They are useful for small-matrix discovery and regression tests,
but they are not part of the trusted path for production certificates.  Exact
certificates use :mod:`trottercert.spectral_gauge` and rational trace moments.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np


def _hermitian_matrix(value: np.ndarray, *, field: str) -> np.ndarray:
    matrix = np.asarray(value, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not matrix.size:
        raise ValueError(f"{field} must be a nonempty square matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{field} must contain finite entries")
    scale = max(1.0, float(np.linalg.norm(matrix, ord=2)))
    tolerance = 64 * np.finfo(float).eps * scale
    if np.linalg.norm(matrix - matrix.conj().T, ord=2) > tolerance:
        raise ValueError(f"{field} must be Hermitian")
    return (matrix + matrix.conj().T) / 2


def _block_tolerance(hamiltonian: np.ndarray) -> float:
    scale = max(1.0, float(np.linalg.norm(hamiltonian, ord=2)))
    return 64 * np.finfo(float).eps * scale


def energy_blocks(h: np.ndarray, *, atol: float) -> tuple[np.ndarray, ...]:
    """Return orthonormal eigenvector isometries for degenerate energy blocks.

    Each returned array has shape ``(dimension, block_multiplicity)``.  Thus
    ``block @ block.conj().T`` is the spectral projector for that block.  The
    grouping is intended only for numerical discovery; exact claims must not
    depend on ``atol``.
    """

    hamiltonian = _hermitian_matrix(h, field="hamiltonian")
    if isinstance(atol, bool) or not np.isscalar(atol):
        raise ValueError("atol must be a finite nonnegative scalar")
    tolerance = float(atol)
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("atol must be a finite nonnegative scalar")

    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    starts = [0]
    # Compare with the first eigenvalue in a block, rather than chaining
    # adjacent near-equalities, so a long tolerance chain cannot merge energy
    # levels whose endpoints are farther apart than ``atol``.
    reference = float(eigenvalues[0])
    for index in range(1, len(eigenvalues)):
        if abs(float(eigenvalues[index]) - reference) > tolerance:
            starts.append(index)
            reference = float(eigenvalues[index])
    starts.append(len(eigenvalues))
    return tuple(
        eigenvectors[:, first:last].copy()
        for first, last in pairwise(starts)
    )


def _pinch_energy_blocks(hamiltonian: np.ndarray, defect: np.ndarray) -> np.ndarray:
    result = np.zeros_like(defect, dtype=complex)
    for block in energy_blocks(hamiltonian, atol=_block_tolerance(hamiltonian)):
        result += block @ (block.conj().T @ defect @ block) @ block.conj().T
    return (result + result.conj().T) / 2


def _project_column_space(
    target: np.ndarray,
    generators: list[np.ndarray],
) -> np.ndarray:
    """Remove the HS column space without forming a normal-equation Gram matrix."""

    if not generators:
        return target.copy()
    columns = np.column_stack([generator.reshape(-1) for generator in generators])
    left, singular_values, _ = np.linalg.svd(columns, full_matrices=False)
    if not len(singular_values) or singular_values[0] == 0:
        return target.copy()
    rank_tolerance = (
        max(columns.shape) * np.finfo(float).eps * float(singular_values[0])
    )
    rank = int(np.count_nonzero(singular_values > rank_tolerance))
    if rank == 0:
        return target.copy()
    basis = left[:, :rank]
    flattened = target.reshape(-1)
    return (flattened - basis @ (basis.conj().T @ flattened)).reshape(target.shape)


def gauge_residual(
    h: np.ndarray,
    e: np.ndarray,
    *,
    remove_identity: bool = True,
    remove_retiming: bool = True,
) -> np.ndarray:
    """Return an HS-orthogonal representative modulo the spectral gauge.

    The image of ``i ad_H`` is removed by pinching ``e`` into the degenerate
    energy blocks of ``h``.  The block-diagonal result is then projected away
    from the requested subset of ``span{I, H}``.  The projection uses an SVD
    of the generator columns directly; it does not form the normal-equation
    Gram matrix, whose condition number would be squared when ``H`` is nearly
    proportional to the identity.
    """

    if not isinstance(remove_identity, bool) or not isinstance(remove_retiming, bool):
        raise TypeError("gauge-removal flags must be boolean")
    hamiltonian = _hermitian_matrix(h, field="hamiltonian")
    defect = _hermitian_matrix(e, field="defect")
    if hamiltonian.shape != defect.shape:
        raise ValueError("hamiltonian and defect must have the same square dimension")

    residual = _pinch_energy_blocks(hamiltonian, defect)
    generators: list[np.ndarray] = []
    if remove_identity:
        generators.append(np.eye(hamiltonian.shape[0], dtype=complex))
    if remove_retiming:
        generators.append(hamiltonian)
    residual = _project_column_space(residual, generators)
    return (residual + residual.conj().T) / 2


def check_dual_witness(
    h: np.ndarray,
    e: np.ndarray,
    w: np.ndarray,
    *,
    atol: float,
) -> dict[str, float | bool]:
    """Check a numerical commutant witness for discovery only.

    ``valid`` means only that the nonzero witness passes the requested floating
    tolerance tests.  It is not a rigorous certificate and this function does
    not output a normalized lower bound.  ``pairing`` is evaluated against the
    numerically projected quotient representative, preventing a large
    gauge-only component from amplifying small orthogonality residuals.
    """

    hamiltonian = _hermitian_matrix(h, field="hamiltonian")
    defect = _hermitian_matrix(e, field="defect")
    witness = _hermitian_matrix(w, field="witness")
    if hamiltonian.shape != defect.shape or hamiltonian.shape != witness.shape:
        raise ValueError("hamiltonian, defect, and witness must have the same dimension")
    if isinstance(atol, bool) or not np.isscalar(atol):
        raise ValueError("atol must be a finite nonnegative scalar")
    tolerance = float(atol)
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("atol must be a finite nonnegative scalar")

    commutes = bool(
        np.linalg.norm(witness @ hamiltonian - hamiltonian @ witness, ord=2)
        <= tolerance
    )
    identity_orthogonal = bool(abs(np.trace(witness)) <= tolerance)
    retiming_orthogonal = bool(
        abs(np.trace(witness @ hamiltonian)) <= tolerance
    )
    trace_norm = float(np.linalg.norm(witness, ord="nuc"))
    nonzero = bool(trace_norm > tolerance)
    quotient = gauge_residual(hamiltonian, defect)
    return {
        "commutes": commutes,
        "identity_orthogonal": identity_orthogonal,
        "retiming_orthogonal": retiming_orthogonal,
        "nonzero": nonzero,
        "valid": bool(
            commutes and identity_orthogonal and retiming_orthogonal and nonzero
        ),
        "rigorous": False,
        "pairing": float(np.real(np.trace(witness @ quotient))),
        "trace_norm": trace_norm,
    }
