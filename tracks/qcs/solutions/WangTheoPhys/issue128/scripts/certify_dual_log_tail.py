from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

from trottercert.dual_log_tail import (
    REFERENCE_CERTIFICATE,
    certify_issue128_dual_log_tail,
    issue128_generator_constants,
)

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/dual-log-tail.json"
)
E5_PATH = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
E7_PATH = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
)
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/dual_log_tail.py",
    ISSUE_ROOT / "src/trottercert/intervals.py",
    Path(__file__).resolve(),
)


def _pair(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def _source_label(path: Path) -> str:
    return str(path.resolve().relative_to(ISSUE_ROOT.resolve()))


def build_payload() -> dict[str, Any]:
    constants = issue128_generator_constants()
    result = certify_issue128_dual_log_tail(constants)
    geometry = result.geometry
    envelope = result.envelope
    return {
        "schema_version": 1,
        "kind": "issue128_dual_log_tail",
        "formula": "five_copy_fourth_order_suzuki_four_matchings",
        "claim_scope": result.claim_scope,
        "instance": {
            "boundary": "periodic",
            "length_x": 12,
            "length_y": 12,
            "n_sites": geometry.n_sites,
            "cells": geometry.cells,
            "steps": geometry.steps,
        },
        "moments": {
            "h_hs_squared": _pair(geometry.h_hs_squared),
            "h_hs_cap": _pair(geometry.h_hs_cap),
            "w_hs_squared": _pair(geometry.w_hs_squared),
            "w_hs_per_cell_cap": _pair(
                geometry.w_hs_per_cell_cap
            ),
            "max_w_pauli_coefficient": _pair(
                geometry.max_w_pauli_coefficient
            ),
            "centered_h_operator_cap": _pair(
                geometry.centered_h_op_cap
            ),
        },
        "analytic_caps": {
            "log_radius": _pair(geometry.log_radius_cap),
            "dexp": _pair(geometry.dexp_cap),
            "multiplier_difference": _pair(
                geometry.multiplier_difference_cap
            ),
            "even_series": _pair(geometry.even_series_cap),
        },
        "right_generator_inputs": {
            "source_certificate_steps": constants.source_steps,
            "coefficient_interval_decimal_digits": (
                constants.coefficient_interval_decimal_digits
            ),
            "stage_count": len(constants.stages),
            "d4_site": _pair(constants.d4_site),
            "d5_site": _pair(constants.d5_site),
            "d6_site": _pair(constants.d6_site),
            "d7_site": _pair(constants.d7_site),
        },
        "envelope": {
            "average_generator_defect": _pair(
                envelope.average_generator_defect
            ),
            "one_step_unitary_defect": _pair(
                envelope.one_step_unitary_defect
            ),
            "pointwise_generator_defect": _pair(
                envelope.pointwise_generator_defect
            ),
            "relative_log_defect": _pair(
                envelope.relative_log_defect
            ),
            "centered_exact_phase_radius": _pair(
                envelope.centered_exact_phase_radius
            ),
            "centered_log_radius_bound": _pair(
                envelope.centered_log_radius_bound
            ),
            "log_defect": _pair(envelope.log_defect),
            "log_derivative_defect_hs": _pair(
                envelope.log_derivative_defect_hs
            ),
            "log_commutator_defect_hs": _pair(
                envelope.log_commutator_defect_hs
            ),
        },
        "bounds": {
            "first_omitted_generator_degree": (
                result.first_omitted_generator_degree
            ),
            "first_omitted_log_degree": result.first_omitted_log_degree,
            "direct_even_generator_tail": _pair(
                result.direct_even_generator_tail
            ),
            "dexp_correction_tail": _pair(
                result.dexp_correction_tail
            ),
            "total_log_tail": _pair(result.total_log_tail),
        },
        "inputs": {
            "right_generator_certificate": {
                "path": _source_label(REFERENCE_CERTIFICATE),
                "sha256": _sha(REFERENCE_CERTIFICATE),
            },
            "e5_extensive_witness": {
                "path": _source_label(E5_PATH),
                "sha256": _sha(E5_PATH),
            },
            "e7_dual_pairing": {
                "path": _source_label(E7_PATH),
                "sha256": _sha(E7_PATH),
            },
        },
        "implementation_sources": {
            _source_label(path): _sha(path) for path in SOURCE_PATHS
        },
        "claim": {
            "dual_e11_plus_tail": "certified",
            "finite_step_status": "inconclusive",
            "missing": "exact dual E9 pairing",
        },
    }


def verify_payload(payload: Mapping[str, object]) -> None:
    if not isinstance(payload, Mapping):
        raise ValueError("dual-log tail payload must be a mapping")
    if dict(payload) != build_payload():
        raise ValueError("payload regeneration mismatch")


def _load_payload(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_bytes())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read dual-log tail artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("dual-log tail artifact must contain a JSON object")
    return payload


def _write_payload(path: Path, *, replace: bool) -> None:
    if replace and path.resolve() != DEFAULT_OUTPUT.resolve():
        raise ValueError("--replace is allowed only for the default artifact")
    if path.exists() and not replace:
        raise FileExistsError(f"refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.pending.{os.getpid()}")
    if temporary.exists():
        raise FileExistsError(f"temporary artifact already exists: {temporary}")
    try:
        temporary.write_bytes(_canonical_bytes(build_payload()))
        verify_payload(_load_payload(temporary))
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build or verify the compatible Issue-128 dual-log tail"
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--verify", type=Path)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    if args.verify is not None:
        if args.replace:
            parser.error("--replace cannot be combined with --verify")
        verify_payload(_load_payload(args.verify))
        print(f"dual-log tail artifact valid: {args.verify}", flush=True)
        return

    _write_payload(args.output, replace=args.replace)
    print(f"wrote dual-log tail artifact: {args.output}", flush=True)


if __name__ == "__main__":
    main()
