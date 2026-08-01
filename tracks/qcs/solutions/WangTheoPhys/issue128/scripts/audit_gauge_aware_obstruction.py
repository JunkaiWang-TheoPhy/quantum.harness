#!/usr/bin/env python3
"""Emit a fail-closed calibration-aware audit of the endpoint no-go claim."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_dual_e7_pairing import verify_reduced_payload
from scripts.certify_extensive_commutant_witness import verify_extensive_payload
from trottercert.commutant_witness import verify_quadratic_witness_payload
from trottercert.cubic_field import Cubic
from trottercert.intervals import cube_root_four_interval
from trottercert.rigorous_fourth import fourth_order_suzuki_interval_stages
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
DEFAULT_EXTENSIVE_WITNESS = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
DEFAULT_DUAL_E7 = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
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


def _read_extensive_witness(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        encoded = path.read_bytes()
        payload = json.loads(encoded)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read extensive witness artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("extensive witness artifact must contain a JSON object")
    try:
        verify_extensive_payload(payload)
    except ValueError as exc:
        raise ValueError(f"extensive witness artifact is invalid: {exc}") from exc
    return payload, encoded


def _read_dual_e7(path: Path) -> tuple[dict[str, Any], bytes]:
    try:
        encoded = path.read_bytes()
        payload = json.loads(encoded)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read dual E7 artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("dual E7 artifact must contain a JSON object")
    try:
        verify_reduced_payload(payload)
    except ValueError as exc:
        raise ValueError(f"dual E7 artifact is invalid: {exc}") from exc
    return payload, encoded


def _fraction(value: list[int]) -> Fraction:
    return Fraction(value[0], value[1])


def _cubic(value: list[list[int]]) -> Cubic:
    return Cubic(*(_fraction(coordinate) for coordinate in value))


def _fraction_json(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [
        _fraction_json(value.a0),
        _fraction_json(value.a1),
        _fraction_json(value.a2),
    ]


def _truncated_log_gate(
    extensive_payload: dict[str, Any],
    dual_e7_payload: dict[str, Any],
) -> dict[str, Any]:
    q5 = _cubic(extensive_payload["stable_relations"]["tau_w_e5_per_cell"])
    q7 = _cubic(dual_e7_payload["pairings"]["tau_w_per_cell"])
    root = cube_root_four_interval(50)
    stages, _ = fourth_order_suzuki_interval_stages(4, decimal_digits=30)
    stage_abs_upper = sum(
        (stage.coefficient.abs_upper() for stage in stages),
        Fraction(),
    )
    centered_fragment_radius = 36
    pi_lower = Fraction(314159, 100000)
    cells = 36
    tau_w2 = Fraction(23085, 4)
    result: dict[str, Any] = {}
    for steps in (95, 96, 97):
        combined = q5 / steps**4 + q7 / steps**6
        enclosure = combined.enclose(root)
        if enclosure.upper >= 0:
            raise ArithmeticError("E5+E7 truncated dual pairing is not negative")
        path_upper = centered_fragment_radius * stage_abs_upper / steps
        result[str(steps)] = {
            "per_cell_exact_cubic": _cubic_json(combined),
            "interval": {
                "lower": _fraction_json(enclosure.lower),
                "upper": _fraction_json(enclosure.upper),
            },
            "sign": "negative",
            "squared_normalized_l12": _cubic_json(
                (cells * combined) ** 2 / tau_w2
            ),
            "centered_stage_path_upper": _fraction_json(path_upper),
            "pi_lower": _fraction_json(pi_lower),
            "branch_path_lt_pi": path_upper < pi_lower,
        }
    return result


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
    extensive_witness: Path = DEFAULT_EXTENSIVE_WITNESS,
    dual_e7: Path = DEFAULT_DUAL_E7,
) -> dict[str, Any]:
    source_payload, encoded = _read_source(source)
    witness_payload, witness_encoded = _read_witness(witness)
    extensive_payload, extensive_encoded = _read_extensive_witness(
        extensive_witness
    )
    dual_e7_payload, dual_e7_encoded = _read_dual_e7(dual_e7)
    overlap = source_payload.get(
        "e5_hilbert_schmidt_overlap_with_h_per_cell"
    )
    if not isinstance(overlap, dict) or overlap.get("nonzero") is not True:
        raise ValueError("source must prove a nonzero exact overlap with H")
    exact_cubic = overlap.get("exact_cubic")
    if not isinstance(exact_cubic, list) or not exact_cubic:
        raise ValueError("source nonzero overlap must include exact cubic coordinates")
    l12_records = [
        record
        for record in extensive_payload["records"]
        if record.get("length") == 12
    ]
    if len(l12_records) != 1:
        raise ValueError("extensive witness must contain exactly one L=12 record")

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
            "extensive_witness_artifact": {
                "path": _source_label(extensive_witness),
                "sha256": hashlib.sha256(extensive_encoded).hexdigest(),
            },
            "dual_e7_artifact": {
                "path": _source_label(dual_e7),
                "sha256": hashlib.sha256(dual_e7_encoded).hexdigest(),
            },
            "dual_e7_pairing_status": "exact",
            "exact_tau_w_e7_per_cell": dual_e7_payload["pairings"][
                "tau_w_per_cell"
            ],
            "verified_lengths": extensive_payload["accepted_lengths"],
            "rejected_alias_lengths": extensive_payload[
                "rejected_alias_lengths"
            ],
            "stable_tau_w_e5_per_cell": extensive_payload[
                "stable_relations"
            ]["tau_w_e5_per_cell"],
            "squared_normalized_pairing_l12": l12_records[0]["moments"][
                "squared_normalized_pairing"
            ],
            "exact_cubic_pairing": witness_payload["pairings"]["tau_w_e5"],
            "claim": "E5 is not in image(i ad_H) + span(I,H)",
            "finite_step_gate": "dual_pairing_remainder",
            "truncated_log_gate": _truncated_log_gate(
                extensive_payload,
                dual_e7_payload,
            ),
            "finite_step_missing": (
                "at r=97, an E9-and-higher dual tail small enough to preserve "
                "the exact E5+E7 pairing"
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
            "certify the E9-and-higher dual tail at r=97 before claiming a "
            "finite-step eigenphase lower bound"
        ),
    }


def verify_payload(
    payload: dict[str, Any],
    source: Path,
    witness: Path = DEFAULT_WITNESS,
    extensive_witness: Path = DEFAULT_EXTENSIVE_WITNESS,
    dual_e7: Path = DEFAULT_DUAL_E7,
) -> None:
    expected = build_payload(source, witness, extensive_witness, dual_e7)
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
    parser.add_argument(
        "--extensive-witness",
        type=Path,
        default=DEFAULT_EXTENSIVE_WITNESS,
    )
    parser.add_argument("--dual-e7", type=Path, default=DEFAULT_DUAL_E7)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    if args.verify:
        try:
            payload = json.loads(args.output.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"cannot read audit output: {args.output}") from exc
        verify_payload(
            payload,
            args.source,
            args.witness,
            args.extensive_witness,
            args.dual_e7,
        )
        print(
            "gauge-aware audit valid: "
            f"source_sha256={payload['source']['sha256']}",
            flush=True,
        )
        return

    payload = build_payload(
        args.source,
        args.witness,
        args.extensive_witness,
        args.dual_e7,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(_encoded(payload))
    print(
        f"wrote gauge-aware audit: {args.output} "
        f"source_sha256={payload['source']['sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
