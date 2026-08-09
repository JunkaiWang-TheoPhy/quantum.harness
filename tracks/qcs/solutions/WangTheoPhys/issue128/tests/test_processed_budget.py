from fractions import Fraction

import pytest

from trottercert.processed_budget import (
    ProcessedFiniteStepInputs,
    evaluate_processed_budget,
)


def test_budget_requires_degree7_and_remainder_proof_inputs() -> None:
    with pytest.raises(ValueError, match="degree-seven"):
        ProcessedFiniteStepInputs(
            kernel_name="s10",
            steps=39,
            n_sites=144,
            tolerance=Fraction(1, 10**6),
            degree7_site_l1=None,
            higher_order_remainder=None,
            local_log_theorem_id=None,
        )
    with pytest.raises(ValueError, match="higher-order"):
        ProcessedFiniteStepInputs(
            kernel_name="s10",
            steps=39,
            n_sites=144,
            tolerance=Fraction(1, 10**6),
            degree7_site_l1=Fraction(),
            higher_order_remainder=None,
            local_log_theorem_id="test-only",
        )


def test_budget_sums_exact_log_contributions() -> None:
    inputs = ProcessedFiniteStepInputs(
        kernel_name="s10",
        steps=39,
        n_sites=144,
        tolerance=Fraction(1, 10**6),
        degree7_site_l1=Fraction(1, 10**6),
        higher_order_remainder=Fraction(1, 10**8),
        local_log_theorem_id="test-only",
    )
    result = evaluate_processed_budget(inputs)
    assert result.total == (
        result.degree3_residual
        + result.degree5_residual
        + result.degree7
        + result.higher_order_remainder
    )
    assert result.accepted
    assert result.remaining_budget == inputs.tolerance - result.total


def test_budget_rejects_unknown_kernel_and_negative_proof_values() -> None:
    with pytest.raises(ValueError, match="unknown processed kernel"):
        ProcessedFiniteStepInputs(
            kernel_name="unknown",
            steps=39,
            n_sites=144,
            tolerance=Fraction(1, 10**6),
            degree7_site_l1=Fraction(),
            higher_order_remainder=Fraction(),
            local_log_theorem_id="test-only",
        )
    with pytest.raises(ValueError, match="positive integer"):
        ProcessedFiniteStepInputs(
            kernel_name="s10",
            steps=39.0,  # type: ignore[arg-type]
            n_sites=144,
            tolerance=Fraction(1, 10**6),
            degree7_site_l1=Fraction(),
            higher_order_remainder=Fraction(),
            local_log_theorem_id="test-only",
        )
    with pytest.raises(ValueError, match="exact rational"):
        ProcessedFiniteStepInputs(
            kernel_name="s10",
            steps=39,
            n_sites=144,
            tolerance=Fraction(1, 10**6),
            degree7_site_l1=0.0,  # type: ignore[arg-type]
            higher_order_remainder=Fraction(),
            local_log_theorem_id="test-only",
        )
    with pytest.raises(ValueError, match="nonnegative"):
        ProcessedFiniteStepInputs(
            kernel_name="s10",
            steps=39,
            n_sites=144,
            tolerance=Fraction(1, 10**6),
            degree7_site_l1=Fraction(-1),
            higher_order_remainder=Fraction(),
            local_log_theorem_id="test-only",
        )
