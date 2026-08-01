from fractions import Fraction

import pytest
import sympy as sp

from trottercert.spectral_gauge import (
    GaugeScope,
    decompose_gauge,
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
