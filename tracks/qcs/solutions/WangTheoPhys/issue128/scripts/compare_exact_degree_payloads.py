#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from trottercert.hpc_artifacts import (
    coordinate_decode_terms,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
)


KIND = "issue128_exact_right_generator_degree"
FORMULA_ID = "five_copy_suzuki_fourth_order_exact_cubic"


def _load(path: Path, degree: int) -> tuple[dict[str, object], dict]:
    payload = read_shard_gzip(path)
    if payload.get("kind") != KIND:
        raise ValueError(f"{path}: exact-degree kind mismatch")
    if payload.get("formula_id") != FORMULA_ID:
        raise ValueError(f"{path}: formula mismatch")
    if payload.get("stage_count") != 31 or payload.get("degree") != degree:
        raise ValueError(f"{path}: stage or degree mismatch")
    terms = coordinate_decode_terms(payload.get("terms"))
    if payload.get("term_count") != len(terms):
        raise ValueError(f"{path}: declared term count mismatch")
    return payload, terms


def compare(left: Path, right: Path, *, degree: int, output: Path) -> dict[str, object]:
    left_payload, left_terms = _load(left, degree)
    right_payload, right_terms = _load(right, degree)
    if left_terms != right_terms:
        raise ArithmeticError(f"degree-{degree} exact coefficient maps differ")
    for field in ("cell_pauli_l1_upper", "site_pauli_l1_upper"):
        if left_payload.get(field) != right_payload.get(field):
            raise ArithmeticError(f"degree-{degree} {field} differs")
    result: dict[str, object] = {
        "schema_version": 1,
        "kind": "issue128_exact_degree_payload_comparison",
        "status": "complete",
        "degree": degree,
        "term_count": len(left_terms),
        "exact_coefficient_map_equal": True,
        "cell_pauli_l1_upper": left_payload["cell_pauli_l1_upper"],
        "site_pauli_l1_upper": left_payload["site_pauli_l1_upper"],
        "left": {
            "path": str(left),
            "sha256": sha256_file(left),
            "source_commit": left_payload.get("source_commit"),
            "parent_sha256": left_payload.get("parent_sha256"),
        },
        "right": {
            "path": str(right),
            "sha256": sha256_file(right),
            "source_commit": right_payload.get("source_commit"),
            "parent_sha256": right_payload.get("parent_sha256"),
        },
    }
    write_manifest_atomic(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.left, args.right, degree=args.degree, output=args.output)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
