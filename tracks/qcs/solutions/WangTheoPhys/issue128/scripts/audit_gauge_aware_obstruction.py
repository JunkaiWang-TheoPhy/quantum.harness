#!/usr/bin/env python3
"""Emit a fail-closed calibration-aware audit of the endpoint no-go claim."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from trottercert.commutant_witness import verify_quadratic_witness_payload
from trottercert.spectral_gauge import decompose_matrix_spectral_gauge

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/exact-obstruction.json"
)
DEFAULT_WITNESS = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/quadratic-commutant-witness.json"
)


def _source_label(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ISSUE_ROOT))
    except ValueError:
        return str(resolved)


def _read_source(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        encoded = path.read_bytes()
        payload = json.loads(encoded)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read source obstruction artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("source obstruction artifact must contain a JSON object")
    return payload, encoded


def _read_witness(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        encoded = path.read_bytes()
        payload = json.loads(encoded)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read quadratic witness artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("quadratic witness artifact must contain a JSON object")
    try:
        verify_quadratic_witness_payload(payload)
    except ValueError as exc:
        raise ValueError(f"quadratic witness artifact is invalid: {exc}") from exc
    return payload, encoded


def _matrix_entries(matrix: sp.MatrixBase | None) -> list[list[str]] | None:
    if matrix is None:
        return None
    return [
        [str(sp.simplify(matrix[row, column])) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _reference_examples() -> dict[str, Any]:
    two_level_h = sp.diag(-1, 1)
    off_diagonal = decompose_matrix_spectral_gauge(
        two_level_h,
        sp.Matrix([[0, 1], [1, 0]]),
    )
    h_fixed = decompose_matrix_spectral_gauge(
        two_level_h,
        two_level_h,
        include_time_calibration=False,
    )
    h_calibrated = decompose_matrix_spectral_gauge(two_level_h, two_level_h)

    three_level_h = sp.diag(-1, 0, 2)
    independent_defect = sp.diag(2, -3, 1)
    independent = decompose_matrix_spectral_gauge(
        three_level_h,
        independent_defect,
    )
    witness = independent.witness_matrix
    if witness is None:
        raise ArithmeticError("reference commutant example did not emit a witness")

    return {
        "off_diagonal": {
            "hamiltonian": _matrix_entries(two_level_h),
            "defect": _matrix_entries(sp.Matrix([[0, 1], [1, 0]])),
            "status": off_diagonal.gauge.status,
        },
        "hamiltonian_parallel": {
            "hamiltonian": _matrix_entries(two_level_h),
            "defect": _matrix_entries(two_level_h),
            "fixed_time_status": h_fixed.gauge.status,
            "fixed_time_witness": _matrix_entries(h_fixed.witness_matrix),
            "calibrated_status": h_calibrated.gauge.status,
        },
        "independent_commutant": {
            "hamiltonian": _matrix_entries(three_level_h),
            "defect": _matrix_entries(independent_defect),
            "status": independent.gauge.status,
            "witness": _matrix_entries(witness),
            "orthogonal_to_identity": sp.trace(witness) == 0,
            "orthogonal_to_hamiltonian": (
                sp.trace(witness * three_level_h) == 0
            ),
            "commutes_with_hamiltonian": (
                witness * three_level_h - three_level_h * witness
                == sp.zeros(3)
            ),
        },
    }


def build_payload(
    source: Path,
    witness: Path = DEFAULT_WITNESS,
) -> dict[str, Any]:
    source_payload, encoded = _read_source(source)
    witness_payload, witness_encoded = _read_witness(witness)
    overlap = source_payload.get(
        "e5_hilbert_schmidt_overlap_with_h_per_cell"
    )
    if not isinstance(overlap, dict) or overlap.get("nonzero") is not True:
        raise ValueError("source must prove a nonzero exact overlap with H")
    exact_cubic = overlap.get("exact_cubic")
    if not isinstance(exact_cubic, list) or not exact_cubic:
        raise ValueError("source nonzero overlap must include exact cubic coordinates")

    return {
        "schema_version": 1,
        "kind": "issue128_gauge_aware_processor_obstruction_audit",
        "source": {
            "path": _source_label(source),
            "sha256": hashlib.sha256(encoded).hexdigest(),
        },
        "allowed_calibration_gauge": [
            "image(i ad_H)",
            "span(I)",
            "span(H)",
        ],
        "fixed_time_endpoint_processor": {
            "status": "no_go",
            "witness": "H",
            "source_nonzero": True,
            "exact_cubic_overlap": exact_cubic,
            "claim": "E5 is not in image(i ad_H) at fixed target time",
        },
        "calibrated_spectral_obstruction": {
            "leading_order_status": "no_go",
            "finite_step_status": "inconclusive",
            "witness": "H^2 - 54 I + H/2",
            "witness_artifact": {
                "path": _source_label(witness),
                "sha256": hashlib.sha256(witness_encoded).hexdigest(),
            },
            "exact_cubic_pairing": witness_payload["pairings"]["tau_w_e5"],
            "claim": "E5 is not in image(i ad_H) + span(I,H)",
            "finite_step_missing": (
                "a certified local-log branch and all-order remainder small "
                "enough to preserve the leading pairing"
            ),
        },
        "restricted_support_processor": {
            "status": "operator_no_go_only",
            "spectral_implication": "not_established",
            "reason": (
                "support-six residuals cannot be canceled by support-at-most-four "
                "processors, but support alone does not lower-bound eigenphase error"
            ),
        },
        "reference_examples": _reference_examples(),
        "hpc_authorized": False,
        "next_gate": (
            "close the processed local-log branch and all-order remainder "
            "before claiming a finite-step eigenphase lower bound or running E7"
        ),
    }


def verify_payload(
    payload: dict[str, Any],
    source: Path,
    witness: Path = DEFAULT_WITNESS,
) -> None:
    expected = build_payload(source, witness)
    submitted_source = payload.get("source")
    if not isinstance(submitted_source, dict):
        raise ValueError("source metadata is missing")
    if submitted_source.get("sha256") != expected["source"]["sha256"]:
        raise ValueError("source digest mismatch")
    if payload.get("calibrated_spectral_obstruction") != expected[
        "calibrated_spectral_obstruction"
    ]:
        raise ValueError("calibrated spectral status mismatch")
    if payload != expected:
        raise ValueError("gauge-aware obstruction payload mismatch")


def _encoded(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--witness", type=Path, default=DEFAULT_WITNESS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    if args.verify:
        try:
            payload = json.loads(args.output.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"cannot read audit output: {args.output}") from exc
        verify_payload(payload, args.source, args.witness)
        print(
            "gauge-aware audit valid: "
            f"source_sha256={payload['source']['sha256']}",
            flush=True,
        )
        return

    payload = build_payload(args.source, args.witness)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(_encoded(payload))
    print(
        f"wrote gauge-aware audit: {args.output} "
        f"source_sha256={payload['source']['sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
