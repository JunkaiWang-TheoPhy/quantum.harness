#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

from trottercert.processed_kernels import (
    alternating_kernel_stages,
    groups_for_repetitions,
    published_effective_order_six_kernel,
)
from trottercert.processed_local_bounds import (
    build_processed_local_bounds,
    processor_distance_upper,
    processor_generator_norm_upper,
)


ISSUE_ROOT = Path(__file__).resolve().parents[1]
LOCAL_AUDIT = (
    ISSUE_ROOT
    / "artifacts/processed-kernel-audit/s10-s11-local-audit.json"
)
POINTS = (("s10", 39), ("s10", 47), ("s11", 35), ("s11", 43))
STATUS = "awaiting_degree7_and_local_log_remainder"
_HASHED_IMPLEMENTATION = (
    ISSUE_ROOT / "src/trottercert/rational_lie_local.py",
    ISSUE_ROOT / "src/trottercert/processed_local_bounds.py",
    ISSUE_ROOT / "src/trottercert/processed_budget.py",
    ISSUE_ROOT / "scripts/report_processed_local_budget.py",
)


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _point_payload(name: str, steps: int) -> dict[str, Any]:
    kernel = published_effective_order_six_kernel(name)
    bounds = build_processed_local_bounds(kernel)
    n_sites = 144
    degree3_total = (
        Fraction(n_sites) * bounds.degree3_residual_site_l1 / steps**2
    )
    degree5_total = (
        Fraction(n_sites) * bounds.degree5_residual_site_l1 / steps**4
    )
    return {
        "kernel": name,
        "steps": steps,
        "n_sites": n_sites,
        "groups": groups_for_repetitions(kernel, steps),
        "stage_count": len(alternating_kernel_stages(kernel)),
        "degree3_residual_cell_l1": _pair(
            bounds.degree3_residual_cell_l1
        ),
        "degree3_residual_site_l1": _pair(
            bounds.degree3_residual_site_l1
        ),
        "degree3_residual_total": _pair(degree3_total),
        "degree5_residual_cell_l1": _pair(
            bounds.degree5_residual_cell_l1
        ),
        "degree5_residual_site_l1": _pair(
            bounds.degree5_residual_site_l1
        ),
        "degree5_residual_total": _pair(degree5_total),
        "processor_generator_norm_upper": _pair(
            processor_generator_norm_upper(bounds, n_sites, steps)
        ),
        "processor_distance_upper": _pair(
            processor_distance_upper(bounds, n_sites, steps)
        ),
        "local_term_counts": {
            "r2": bounds.r2_term_count,
            "r4": bounds.r4_term_count,
            "degree3_residual": bounds.degree3_residual_term_count,
            "degree5_residual": bounds.degree5_residual_term_count,
        },
        "degree7_site_l1": None,
        "higher_order_remainder": None,
        "local_log_theorem": None,
        "status": STATUS,
    }


def build_report(*, progress: bool = False) -> dict[str, Any]:
    if not LOCAL_AUDIT.is_file():
        raise FileNotFoundError(f"missing processed kernel audit: {LOCAL_AUDIT}")
    points = []
    for name, steps in POINTS:
        if progress:
            print(f"building local budget kernel={name} steps={steps}", flush=True)
        points.append(_point_payload(name, steps))
    return {
        "schema_version": 1,
        "kind": "issue128_processed_pre_e7_local_budget",
        "benchmark": {
            "model": "periodic_12x12_isotropic_heisenberg",
            "n_sites": 144,
            "time": [1, 1],
            "tolerance": [1, 10**6],
            "fragment_count": 4,
        },
        "claim_boundary": (
            "residual and processor bounds only; physical degree7, "
            "all-order remainder, and local-log theorem are missing"
        ),
        "processed_kernel_audit": {
            "path": str(LOCAL_AUDIT.relative_to(ISSUE_ROOT)),
            "sha256": _sha256(LOCAL_AUDIT),
        },
        "implementation_sha256": {
            path.name: _sha256(path) for path in _HASHED_IMPLEMENTATION
        },
        "points": points,
    }


def _canonical_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def write_report(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_canonical_bytes(payload))


def _load_report(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("processed local budget is not canonical JSON") from error
    if not isinstance(payload, dict) or raw != _canonical_bytes(payload):
        raise ValueError("processed local budget is not canonical JSON")
    return payload


def verify_report(path: Path) -> tuple[tuple[str, int], ...]:
    payload = _load_report(path)
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported processed local budget schema")
    if payload.get("kind") != "issue128_processed_pre_e7_local_budget":
        raise ValueError("unexpected processed local budget kind")
    expected = build_report()
    actual_points = payload.get("points")
    if not isinstance(actual_points, list) or len(actual_points) != len(POINTS):
        raise ValueError("processed local budget point coverage mismatch")
    for actual, rebuilt in zip(actual_points, expected["points"]):
        if actual.get("status") != STATUS:
            raise ValueError("processed local budget status must await proof")
        if actual.get("processor_distance_upper") != rebuilt.get(
            "processor_distance_upper"
        ):
            raise ValueError("processed local budget processor distance mismatch")
        if actual != rebuilt:
            raise ValueError("processed local budget point mismatch")
    if payload != expected:
        raise ValueError("processed local budget metadata mismatch")
    return POINTS


def _decimal(pair: list[int]) -> str:
    return f"{float(Fraction(*pair)):.12e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_report(progress=True)
    write_report(args.output, payload)
    if args.verify:
        verify_report(args.output)
    for point in payload["points"]:
        print(
            " ".join(
                (
                    f"kernel={point['kernel']}",
                    f"steps={point['steps']}",
                    f"d3={_decimal(point['degree3_residual_total'])}",
                    f"d5={_decimal(point['degree5_residual_total'])}",
                    f"rho={_decimal(point['processor_distance_upper'])}",
                    f"status={point['status']}",
                )
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
