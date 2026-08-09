from fractions import Fraction

from trottercert.processed_kernels import published_effective_order_six_kernel
from trottercert.processed_local_bounds import (
    build_processed_local_bounds,
    processor_distance_upper,
    processor_generator_norm_upper,
    rational_exp_upper,
)


def test_rational_exp_upper_encloses_simple_values() -> None:
    assert rational_exp_upper(Fraction()) == 1
    assert rational_exp_upper(Fraction(1, 10)) > Fraction(11, 10)
    assert rational_exp_upper(Fraction(1, 10)) < Fraction(111, 100)


def test_s10_and_s11_local_bounds_are_nonzero_and_ordered() -> None:
    s10 = build_processed_local_bounds(
        published_effective_order_six_kernel("s10")
    )
    s11 = build_processed_local_bounds(
        published_effective_order_six_kernel("s11")
    )
    assert 0 < s10.degree3_residual_site_l1 < Fraction(1, 10**25)
    assert 0 < s11.degree3_residual_site_l1 < s10.degree3_residual_site_l1
    assert 0 < s10.degree5_residual_site_l1 < Fraction(1, 10**25)
    assert 0 < s11.degree5_residual_site_l1 < s10.degree5_residual_site_l1
    assert s11.r2_site_l1 < s10.r2_site_l1
    assert processor_distance_upper(s10, 144, 39) < Fraction(1, 10)
    assert processor_distance_upper(s11, 144, 35) < Fraction(1, 1000)


def test_processor_distance_is_monotone_in_steps() -> None:
    bounds = build_processed_local_bounds(
        published_effective_order_six_kernel("s10")
    )
    assert processor_generator_norm_upper(bounds, 144, 40) < (
        processor_generator_norm_upper(bounds, 144, 39)
    )
    assert processor_distance_upper(bounds, 144, 40) < (
        processor_distance_upper(bounds, 144, 39)
    )
