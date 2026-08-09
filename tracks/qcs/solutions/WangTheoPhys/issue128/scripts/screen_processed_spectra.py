#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import expm

from trottercert.processed_kernels import (
    alternating_kernel_stages,
    published_effective_order_six_kernel,
)


def _random_fragments(
    dimension: int,
    seed: int,
) -> tuple[np.ndarray, ...]:
    if dimension < 2:
        raise ValueError("dense screen dimension must be at least two")
    generator = np.random.default_rng(seed)
    fragments = []
    for _ in range(4):
        raw = generator.normal(size=(dimension, dimension)) + 1j * (
            generator.normal(size=(dimension, dimension))
        )
        hermitian = (raw + raw.conj().T) / 2
        fragments.append(hermitian * (0.27 / np.linalg.norm(hermitian, 2)))
    return tuple(fragments)


def _kernel_matrix(
    name: str,
    fragments: Sequence[np.ndarray],
    step_size: float,
) -> np.ndarray:
    stages = alternating_kernel_stages(
        published_effective_order_six_kernel(name)
    )
    dimension = fragments[0].shape[0]
    unitary = np.eye(dimension, dtype=complex)
    for stage in stages:
        unitary = unitary @ expm(
            -1j
            * float(stage.coefficient)
            * step_size
            * fragments[stage.fragment_index]
        )
    return unitary


def _phase_set_error(left: np.ndarray, right: np.ndarray) -> float:
    left_phases = np.sort(np.angle(np.linalg.eigvals(left)))
    right_phases = np.sort(np.angle(np.linalg.eigvals(right)))
    return float(np.max(np.abs(left_phases - right_phases)))


def _slope(repetitions: Sequence[int], errors: Sequence[float]) -> float:
    if any(error <= 0 for error in errors):
        raise ArithmeticError("dense screen encountered a nonpositive error")
    # Do not let the binary64 eigensolver floor flatten the phase fit.  The
    # raw tail values remain in the emitted diagnostic for inspection.
    floor = 64 * np.finfo(float).eps
    selected = [
        (repetition, error)
        for repetition, error in zip(repetitions, errors)
        if error > floor
    ]
    if len(selected) < 3:
        selected = list(zip(repetitions, errors))[:3]
    fit = np.polyfit(
        np.log(np.asarray([item[0] for item in selected], dtype=float)),
        np.log(np.asarray([item[1] for item in selected], dtype=float)),
        1,
    )
    return float(-fit[0])


def run_random_fragment_screen(
    *,
    kernels: Sequence[str],
    dimension: int,
    repetitions: Sequence[int],
    seed: int,
) -> dict[str, dict[str, object]]:
    if not repetitions or any(
        isinstance(value, bool) or value < 1 for value in repetitions
    ):
        raise ValueError("repetitions must contain positive integers")
    fragments = _random_fragments(dimension, seed)
    exact = expm(-1j * sum(fragments))
    result: dict[str, dict[str, object]] = {}
    for name in kernels:
        operator_errors = []
        phase_errors = []
        for repetition in repetitions:
            step = _kernel_matrix(name, fragments, 1 / repetition)
            repeated = np.linalg.matrix_power(step, repetition)
            operator_errors.append(
                float(np.linalg.norm(repeated - exact, 2))
            )
            phase_errors.append(_phase_set_error(repeated, exact))
        result[name] = {
            "operator_slope": _slope(repetitions, operator_errors),
            "phase_slope": _slope(repetitions, phase_errors),
            "operator_errors": operator_errors,
            "phase_errors": phase_errors,
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernels", nargs="+", default=("s8", "s10", "s11"))
    parser.add_argument("--dimension", type=int, default=8)
    parser.add_argument(
        "--repetitions",
        nargs="+",
        type=int,
        default=(2, 3, 4, 5, 6, 8, 10, 12),
    )
    parser.add_argument("--seed", type=int, default=128)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(
        f"dense processed-spectral screen dimension={args.dimension} "
        f"kernels={','.join(args.kernels)}",
        flush=True,
    )
    results = run_random_fragment_screen(
        kernels=args.kernels,
        dimension=args.dimension,
        repetitions=args.repetitions,
        seed=args.seed,
    )
    payload = {
        "schema_version": 1,
        "kind": "issue128_processed_spectral_dense_diagnostic",
        "claim_boundary": (
            "nonrigorous discovery diagnostic; not a certificate bound"
        ),
        "dimension": args.dimension,
        "fragment_count": 4,
        "fragment_spectral_norm": 0.27,
        "seed": args.seed,
        "repetitions": args.repetitions,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for name, values in results.items():
        print(
            f"kernel={name} operator_slope={values['operator_slope']:.8f} "
            f"phase_slope={values['phase_slope']:.8f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
