from scripts.screen_processed_spectra import run_random_fragment_screen


def test_processed_spectral_screen_separates_operator_and_phase_order() -> None:
    result = run_random_fragment_screen(
        kernels=("s10", "s11"),
        dimension=6,
        repetitions=(2, 3, 4, 5, 6),
        seed=128,
    )
    assert 1.7 < result["s10"]["operator_slope"] < 2.3
    assert 5.3 < result["s10"]["phase_slope"] < 6.7
    assert 3.5 < result["s11"]["operator_slope"] < 4.5
    assert 5.3 < result["s11"]["phase_slope"] < 6.7


def test_processed_spectral_screen_is_deterministic() -> None:
    arguments = {
        "kernels": ("s10",),
        "dimension": 4,
        "repetitions": (2, 3, 4, 5),
        "seed": 128,
    }
    assert run_random_fragment_screen(**arguments) == run_random_fragment_screen(
        **arguments
    )
