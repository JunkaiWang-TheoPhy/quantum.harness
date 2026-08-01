from fractions import Fraction

import pytest
import sympy as sp

from trottercert.spectral_gauge import (
    GaugeScope,
    decompose_gauge,
    decompose_matrix_spectral_gauge,
    hermitian_basis,
    hermitian_coordinates,
    matrix_from_hermitian_coordinates,
    verify_gauge_decomposition,
)


def test_exact_primal_reconstruction() -> None:
    result = decompose_gauge(
        (Fraction(2), Fraction(3)),
        ((1, 0), (0, 1)),
        scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
        completeness_id="test-full-span",
    )
    assert result.status == "removable"
    assert result.residual == (0, 0)
    assert result.primitive_witness is None
    verify_gauge_decomposition(result)


def test_full_and_restricted_dual_claims_differ() -> None:
    full = decompose_gauge(
        (1, 2, 3),
        ((1, 0, 0),),
        scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
        completeness_id="test-full-span",
    )
    restricted = decompose_gauge(
        (1, 2, 3),
        ((1, 0, 0),),
        scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
    )
    assert full.status == "spectral_obstruction"
    assert restricted.status == "restricted_obstruction"
    assert full.primitive_witness == restricted.primitive_witness == (0, 2, 3)
    assert full.witness_target_pairing == 13
    verify_gauge_decomposition(full)
    verify_gauge_decomposition(restricted)


def test_weighted_metric_changes_projection() -> None:
    result = decompose_gauge(
        (1, 0),
        ((1, 1),),
        metric=(1, 2),
        scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
    )
    assert result.projection == (Fraction(1, 3), Fraction(1, 3))
    assert result.residual == (Fraction(2, 3), Fraction(-1, 3))
    assert result.primitive_witness == (2, -1)
    assert result.witness_target_pairing == 2
    verify_gauge_decomposition(result)


def test_dependent_generators_are_reduced_exactly() -> None:
    result = decompose_gauge(
        (2, 5),
        ((1, 0), (2, 0)),
        scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
    )
    assert result.pivot_indices == (0,)
    assert result.independent_generators == ((1, 0),)
    assert result.projection == (2, 0)
    assert result.primitive_witness == (0, 1)
    verify_gauge_decomposition(result)


@pytest.mark.parametrize(
    ("target", "generators", "metric", "match"),
    [
        ((), (), None, "nonempty"),
        ((0.1, 1), (), None, "exact"),
        ((sp.Float("0.1"), 1), (), None, "exact"),
        ((True, 1), (), None, "exact"),
        ((1, sp.I), (), None, "real"),
        ((1, 2), ((1,),), None, "dimension"),
        ((1, 2), (), (1,), "metric"),
        ((1, 2), (), (1, 0), "positive"),
    ],
)
def test_invalid_exact_vector_inputs_fail_closed(
    target: tuple[object, ...],
    generators: tuple[tuple[object, ...], ...],
    metric: tuple[object, ...] | None,
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        decompose_gauge(
            target,
            generators,
            metric=metric,
            scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
        )


def test_full_scope_requires_completeness_evidence() -> None:
    with pytest.raises(ValueError, match="completeness"):
        decompose_gauge(
            (1, 0),
            (),
            scope=GaugeScope.FULL_COMMUTATOR_IMAGE,
        )


def test_verifier_rejects_mutated_status() -> None:
    result = decompose_gauge(
        (0, 1),
        ((1, 0),),
        scope=GaugeScope.RESTRICTED_PROCESSOR_SPAN,
    )
    forged = result.__class__(**{**result.__dict__, "status": "spectral_obstruction"})
    with pytest.raises(ValueError, match="status"):
        verify_gauge_decomposition(forged)


def test_hermitian_coordinate_round_trip_and_metric_order() -> None:
    matrix = sp.Matrix(
        [
            [2, 3 + 5 * sp.I, 7 - 11 * sp.I],
            [3 - 5 * sp.I, 13, 17 + 19 * sp.I],
            [7 + 11 * sp.I, 17 - 19 * sp.I, 23],
        ]
    )
    coordinates = hermitian_coordinates(matrix)
    assert coordinates == (2, 13, 23, 3, 7, 17, 5, -11, 19)
    assert matrix_from_hermitian_coordinates(coordinates, 3) == matrix
    basis = hermitian_basis(3)
    assert len(basis) == 9
    assert tuple(sp.trace(item * item) for item in basis) == (
        1,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
    )


def test_off_diagonal_defect_is_commutator_removable() -> None:
    hamiltonian = sp.diag(-1, 1)
    defect = sp.Matrix([[0, 1], [1, 0]])
    result = decompose_matrix_spectral_gauge(hamiltonian, defect)
    assert result.gauge.status == "removable"
    assert result.witness_matrix is None
    assert result.residual_matrix == sp.zeros(2)


def test_time_calibration_absorbs_hamiltonian_parallel_defect() -> None:
    hamiltonian = sp.diag(-1, 1)
    fixed_time = decompose_matrix_spectral_gauge(
        hamiltonian,
        hamiltonian,
        include_time_calibration=False,
    )
    calibrated = decompose_matrix_spectral_gauge(
        hamiltonian,
        hamiltonian,
        include_time_calibration=True,
    )
    assert fixed_time.gauge.status == "spectral_obstruction"
    assert calibrated.gauge.status == "removable"
    assert fixed_time.witness_matrix == -hamiltonian
    assert calibrated.witness_matrix is None


def test_commutant_witness_survives_phase_and_time_gauge() -> None:
    hamiltonian = sp.diag(-1, 0, 2)
    defect = sp.diag(2, -3, 1)
    result = decompose_matrix_spectral_gauge(hamiltonian, defect)
    assert result.gauge.status == "spectral_obstruction"
    assert result.witness_matrix is not None
    commutator = result.witness_matrix * hamiltonian - hamiltonian * result.witness_matrix
    assert commutator.applyfunc(sp.simplify) == sp.zeros(3)
    assert sp.trace(result.witness_matrix) == 0
    assert sp.trace(result.witness_matrix * hamiltonian) == 0
    assert sp.trace(result.witness_matrix * defect) != 0


def test_degenerate_block_defect_is_a_spectral_obstruction() -> None:
    hamiltonian = sp.diag(0, 0, 2)
    defect = sp.diag(1, -1, 0)
    result = decompose_matrix_spectral_gauge(hamiltonian, defect)
    assert result.gauge.status == "spectral_obstruction"
    assert result.witness_matrix == defect


def test_matrix_inputs_reject_floats_and_nonhermitian_values() -> None:
    with pytest.raises(ValueError, match="exact"):
        decompose_matrix_spectral_gauge(sp.diag(0.1, 1), sp.eye(2))
    with pytest.raises(ValueError, match="Hermitian"):
        decompose_matrix_spectral_gauge(
            sp.Matrix([[0, 1], [0, 0]]),
            sp.eye(2),
        )
    with pytest.raises(ValueError, match="same square dimension"):
        decompose_matrix_spectral_gauge(sp.eye(2), sp.eye(3))
