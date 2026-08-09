from fractions import Fraction

import pytest

from trottercert.effective_processor import audit_effective_order_six
from trottercert.processed_kernels import published_effective_order_six_kernel
from trottercert.rational_words import word_l1


@pytest.mark.parametrize(
    (
        "name",
        "maximum_degree_three_l1",
        "maximum_degree_five_l1",
        "minimum_degree_seven_l1",
        "maximum_degree_seven_l1",
    ),
    [
        ("s8", Fraction(1, 10**20), Fraction(1, 10**20), Fraction(1), Fraction(3)),
        (
            "s10",
            Fraction(1, 10**24),
            Fraction(1, 10**24),
            Fraction(1, 2),
            Fraction(1),
        ),
        (
            "s11",
            Fraction(1, 10**26),
            Fraction(1, 10**26),
            Fraction(1, 3),
            Fraction(3, 4),
        ),
    ],
)
def test_rationalized_kernels_have_recorded_small_residuals(
    name: str,
    maximum_degree_three_l1: Fraction,
    maximum_degree_five_l1: Fraction,
    minimum_degree_seven_l1: Fraction,
    maximum_degree_seven_l1: Fraction,
) -> None:
    audit = audit_effective_order_six(
        published_effective_order_six_kernel(name)
    )
    assert audit.processed_degree_three_l1 > 0
    assert audit.processed_degree_five_l1 > 0
    assert audit.processed_degree_three_l1 < maximum_degree_three_l1
    assert audit.processed_degree_five_l1 < maximum_degree_five_l1
    assert len(audit.r2.coordinates) == 6
    assert len(audit.r4.coordinates) == 60
    assert minimum_degree_seven_l1 < audit.processed_degree_seven_l1
    assert audit.processed_degree_seven_l1 < maximum_degree_seven_l1


def test_s11_processor_has_negligible_degree_two_generator() -> None:
    audit = audit_effective_order_six(
        published_effective_order_six_kernel("s11")
    )
    assert word_l1(audit.r2.polynomial) < Fraction(1, 10**25)
    assert Fraction(1, 100) < word_l1(audit.r4.polynomial) < Fraction(1, 10)
