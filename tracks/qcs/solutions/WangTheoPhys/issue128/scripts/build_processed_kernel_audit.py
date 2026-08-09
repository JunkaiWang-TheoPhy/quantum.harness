#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from trottercert.effective_processor import audit_effective_order_six
from trottercert.processed_kernels import (
    alternating_kernel_stages,
    maximum_steps_for_group_budget,
    published_effective_order_six_kernel,
)
from trottercert.rational_words import WordPolynomial, word_l1


ISSUE_ROOT = Path(__file__).resolve().parents[1]
_IMPLEMENTATION_FILES = (
    ISSUE_ROOT / "src/trottercert/processed_kernels.py",
    ISSUE_ROOT / "src/trottercert/rational_words.py",
    ISSUE_ROOT / "src/trottercert/effective_processor.py",
    ISSUE_ROOT / "src/trottercert/spectral_processing.py",
    ISSUE_ROOT / "scripts/build_processed_kernel_audit.py",
    ISSUE_ROOT / "scripts/verify_processed_kernel_audit.py",
)


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _word_text(word: tuple[int, ...]) -> str:
    if any(not 0 <= letter <= 9 for letter in word):
        raise ValueError("word alphabet cannot be encoded as digits")
    return "".join(str(letter) for letter in word)


def encode_polynomial(polynomial: WordPolynomial) -> list[list[object]]:
    return [
        [_word_text(word), coefficient.numerator, coefficient.denominator]
        for word, coefficient in sorted(polynomial.items())
    ]


def _encode_coordinates(
    basis_words: Sequence[tuple[int, ...]],
    coordinates: Sequence[Fraction],
) -> list[list[object]]:
    if len(basis_words) != len(coordinates):
        raise ValueError("processor basis and coordinate lengths disagree")
    return [
        [_word_text(word), coordinate.numerator, coordinate.denominator]
        for word, coordinate in zip(basis_words, coordinates)
    ]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _implementation_hashes() -> dict[str, str]:
    return {path.name: _sha256(path) for path in _IMPLEMENTATION_FILES}


def _kernel_payload(name: str) -> dict[str, Any]:
    kernel = published_effective_order_six_kernel(name)
    stages = alternating_kernel_stages(kernel)
    audit = audit_effective_order_six(kernel)
    r2_l1 = word_l1(audit.r2.polynomial)
    r4_l1 = word_l1(audit.r4.polynomial)
    return {
        "name": name,
        "source_strings": list(kernel.source_strings),
        "composition_coefficients": [
            _pair(value) for value in kernel.composition_coefficients
        ],
        "composition_stage_count": kernel.composition_stage_count,
        "stage_count": len(stages),
        "resource_gates": {
            "beat_current_group_budget": 2850,
            "beat_current_steps": maximum_steps_for_group_budget(kernel, 2850),
            "fivefold_group_budget": 2358,
            "fivefold_steps": maximum_steps_for_group_budget(kernel, 2358),
        },
        "processor": {
            "r2_lyndon": _encode_coordinates(
                audit.r2.basis_words,
                audit.r2.coordinates,
            ),
            "r4_lyndon": _encode_coordinates(
                audit.r4.basis_words,
                audit.r4.coordinates,
            ),
            "r6_policy": "zero",
            "r2_word_l1": _pair(r2_l1),
            "r4_word_l1": _pair(r4_l1),
        },
        "residuals": {
            "degree3": encode_polynomial(audit.processed_degree_three),
            "degree5": encode_polynomial(audit.processed_degree_five),
        },
        "processed_degree7": encode_polynomial(
            audit.processed_degree_seven
        ),
        "word_l1": {
            "processed_degree3_residual": _pair(
                audit.processed_degree_three_l1
            ),
            "processed_degree5_residual": _pair(
                audit.processed_degree_five_l1
            ),
            "processed_degree7": _pair(
                audit.processed_degree_seven_l1
            ),
        },
    }


def build_payload(
    kernel_names: Sequence[str],
    *,
    progress: bool = False,
) -> dict[str, Any]:
    names = tuple(kernel_names)
    if not names:
        raise ValueError("at least one processed kernel is required")
    if len(set(names)) != len(names):
        raise ValueError("processed kernel names must be unique")
    kernels = []
    for name in names:
        if progress:
            print(f"auditing processed kernel {name}", flush=True)
        kernels.append(_kernel_payload(name))
    return {
        "schema_version": 1,
        "kind": "issue128_processed_kernel_local_audit",
        "rationalization_policy": (
            "rationalized_from_published_decimal_strings"
        ),
        "generator_count": 4,
        "source": {
            "paper_arxiv": "2404.04340",
            "doi": "10.5281/zenodo.10814897",
            "file": "metemprar.py",
            "md5": "35a37af90096f371e7c82366cc7fad8b",
        },
        "processor_convention": {
            "generator": "R(h)=h^2*R2+h^4*R4",
            "conjugation": "exp(ad_R) log(kernel)",
            "r6_policy": "zero",
        },
        "implementation_sha256": _implementation_hashes(),
        "kernels": kernels,
    }


def canonical_json_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def write_payload(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json_bytes(payload))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernels", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = build_payload(args.kernels, progress=True)
    write_payload(args.output, payload)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": _sha256(args.output),
                "kernels": args.kernels,
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
