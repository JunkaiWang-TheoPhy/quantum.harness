#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from scripts.build_processed_kernel_audit import (
    build_payload,
    canonical_json_bytes,
)


def _load_canonical(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("processed kernel audit is not canonical JSON") from error
    if not isinstance(payload, dict):
        raise ValueError("processed kernel audit root must be an object")
    if raw != canonical_json_bytes(payload):
        raise ValueError("processed kernel audit is not canonical JSON")
    return payload


def _fraction_text(pair: object) -> str:
    if (
        not isinstance(pair, list)
        or len(pair) != 2
        or isinstance(pair[0], bool)
        or isinstance(pair[1], bool)
        or not isinstance(pair[0], int)
        or not isinstance(pair[1], int)
        or pair[1] <= 0
    ):
        raise ValueError("word-l1 value is not a rational pair")
    value = Fraction(pair[0], pair[1])
    return f"{float(value):.12e}"


def verify_payload(path: Path) -> tuple[str, ...]:
    payload = _load_canonical(path)
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported processed kernel audit schema")
    if payload.get("kind") != "issue128_processed_kernel_local_audit":
        raise ValueError("unexpected processed kernel audit kind")
    raw_kernels = payload.get("kernels")
    if not isinstance(raw_kernels, list) or not raw_kernels:
        raise ValueError("processed kernel audit has no kernels")
    names = tuple(
        kernel.get("name") if isinstance(kernel, dict) else None
        for kernel in raw_kernels
    )
    if any(not isinstance(name, str) for name in names):
        raise ValueError("processed kernel name is invalid")
    expected = build_payload(names)
    if payload.get("implementation_sha256") != expected.get(
        "implementation_sha256"
    ):
        raise ValueError("processed kernel implementation hash mismatch")
    if payload.get("source") != expected.get("source"):
        raise ValueError("processed kernel source provenance mismatch")
    expected_kernels = expected["kernels"]
    for actual, rebuilt in zip(raw_kernels, expected_kernels):
        if actual.get("stage_count") != rebuilt.get("stage_count"):
            raise ValueError(f"{actual.get('name')} stage count mismatch")
        if actual.get("resource_gates") != rebuilt.get("resource_gates"):
            raise ValueError(f"{actual.get('name')} resource gate mismatch")
        if actual.get("word_l1") != rebuilt.get("word_l1"):
            raise ValueError(f"{actual.get('name')} word-l1 mismatch")
        if actual != rebuilt:
            raise ValueError(f"{actual.get('name')} exact audit mismatch")
    if payload != expected:
        raise ValueError("processed kernel audit metadata mismatch")
    return names


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    payload = _load_canonical(args.path)
    names = verify_payload(args.path)
    for kernel in payload["kernels"]:
        l1 = kernel["word_l1"]
        gates = kernel["resource_gates"]
        print(
            " ".join(
                (
                    f"kernel={kernel['name']}",
                    f"d3_residual_l1={_fraction_text(l1['processed_degree3_residual'])}",
                    f"d5_residual_l1={_fraction_text(l1['processed_degree5_residual'])}",
                    f"d7_l1={_fraction_text(l1['processed_degree7'])}",
                    f"beat_steps={gates['beat_current_steps']}",
                    f"fivefold_steps={gates['fivefold_steps']}",
                )
            ),
            flush=True,
        )
    print(f"verified kernels={','.join(names)}", flush=True)


if __name__ == "__main__":
    main()
