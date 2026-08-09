#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import exact_d6_density, exact_log_e5_density
from trottercert.hpc_artifacts import (
    CoordinateCubicTerms,
    coordinate_decode_terms,
    coordinate_encode_terms,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)
from trottercert.intervals import cube_root_four_interval
from trottercert.local_commutators import CoordinateRegistry


def _merge(target: CoordinateCubicTerms, raw_terms: object) -> None:
    for key, coefficient in coordinate_decode_terms(raw_terms).items():
        updated = target.get(key, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            target.pop(key, None)
        else:
            target[key] = updated


def _symplectic(
    registry: CoordinateRegistry, terms: CoordinateCubicTerms
) -> dict[tuple[int, int], Cubic]:
    result: dict[tuple[int, int], Cubic] = {}
    for key, coefficient in terms.items():
        x_mask = 0
        z_mask = 0
        for x, y, operator in key:
            bit = 1 << registry.site((x, y))
            if operator in {"X", "Y"}:
                x_mask |= bit
            if operator in {"Y", "Z"}:
                z_mask |= bit
        result[(x_mask, z_mask)] = coefficient
    return result


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def finalize_partials(
    *,
    input_root: Path,
    partial_count: int,
    source_shard_count: int,
    output: Path,
    parent: Path,
    summary: Path,
    e5_override: CoordinateCubicTerms | None = None,
) -> dict[str, object]:
    if partial_count <= 0 or source_shard_count <= 0:
        raise ValueError("partial and source shard counts must be positive")
    forward: CoordinateCubicTerms = {}
    provenance: list[dict[str, object]] = []
    expected_shard_start = 0
    expected_word_first: int | None = None
    metadata: tuple[str, str, int, int, int] | None = None
    for index in range(partial_count):
        cell = input_root / f"partial-{index:04d}"
        payload_path = cell / "partial.json.gz"
        manifest_path = cell / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest.get("status") != "complete"
            or manifest.get("output_sha256") != sha256_file(payload_path)
        ):
            raise ValueError(f"partial {index} is incomplete or corrupt")
        payload = read_shard_gzip(payload_path)
        if (
            payload.get("kind") != "issue128_exact_e7_partial"
            or payload.get("shard_count") != source_shard_count
        ):
            raise ValueError(f"partial {index} metadata mismatch")
        if int(payload.get("shard_start", -1)) != expected_shard_start:
            raise ValueError("partial shard coverage has a gap or overlap")
        expected_shard_start = int(payload.get("shard_stop", -1))
        current_metadata = (
            str(payload.get("formula_id")),
            str(payload.get("source_commit")),
            int(payload.get("stage_count", -1)),
            int(payload.get("degree", -1)),
            int(payload.get("word_total", -1)),
        )
        if metadata is None:
            metadata = current_metadata
            expected_word_first = int(payload["word_first"])
        elif current_metadata != metadata:
            raise ValueError("partial source metadata mismatch")
        if int(payload["word_first"]) != expected_word_first:
            raise ValueError("partial word coverage has a gap or overlap")
        expected_word_first = int(payload["word_last"])
        _merge(forward, payload["terms"])
        provenance.append(
            {
                "partial_index": index,
                "shard_start": int(payload["shard_start"]),
                "shard_stop": int(payload["shard_stop"]),
                "payload_sha256": sha256_file(payload_path),
                "manifest_sha256": sha256_file(manifest_path),
                "sources": payload["sources"],
            }
        )
    if expected_shard_start != source_shard_count:
        raise ValueError("partial shard coverage has a gap or overlap")
    assert metadata is not None
    assert expected_word_first is not None
    if expected_word_first != metadata[4]:
        raise ValueError("partial word coverage is incomplete")

    registry = CoordinateRegistry()
    e7 = _symplectic(registry, forward)
    if e5_override is None:
        e5_registry, e5_raw = exact_log_e5_density(
            fourth_order_suzuki_cubic_stages()
        )
        e5_coordinates = coordinate_decode_terms(
            coordinate_encode_terms(e5_registry, e5_raw)
        )
    else:
        e5_coordinates = e5_override
    e5 = _symplectic(registry, e5_coordinates)
    d6 = exact_d6_density(registry, e5, e7)
    root = cube_root_four_interval(30)
    cell_l1 = sum(
        (coefficient.enclose(root).abs_upper() for coefficient in d6.values()),
        Fraction(),
    )
    parent_payload = {
        "schema_version": 1,
        "kind": "issue128_e7_tree_reduction_parent",
        "source_commit": metadata[1],
        "source_shard_count": source_shard_count,
        "partial_count": partial_count,
        "reducer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "partials": provenance,
        "e7_term_count": len(forward),
        "identity": "D6=7*E7+(2/3)*ad_A^2(E5)",
    }
    write_manifest_atomic(parent, parent_payload)
    degree_payload = {
        "schema_version": 1,
        "kind": "issue128_exact_right_generator_degree",
        "formula_id": metadata[0],
        "source_commit": metadata[1],
        "stage_count": metadata[2],
        "degree": 6,
        "coefficient_interval_decimal_digits": 30,
        "term_count": len(d6),
        "terms": coordinate_encode_terms(registry, d6),
        "cell_pauli_l1_upper": _pair(cell_l1),
        "site_pauli_l1_upper": _pair(cell_l1 / 4),
        "parent_sha256": sha256_file(parent),
    }
    write_shard_gzip(output, degree_payload)
    result = {
        "schema_version": 1,
        "kind": "issue128_e7_tree_reduction_summary",
        "status": "complete",
        "source_commit": metadata[1],
        "source_shard_count": source_shard_count,
        "partial_count": partial_count,
        "e7_term_count": len(forward),
        "d6_term_count": len(d6),
        "cell_pauli_l1_upper": _pair(cell_l1),
        "site_pauli_l1_upper": _pair(cell_l1 / 4),
        "parent_sha256": sha256_file(parent),
        "output_sha256": sha256_file(output),
        "coverage_check": "complete_gap_free",
    }
    write_manifest_atomic(summary, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--partial-count", type=int, required=True)
    parser.add_argument("--source-shard-count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    result = finalize_partials(
        input_root=args.input_root,
        partial_count=args.partial_count,
        source_shard_count=args.source_shard_count,
        output=args.output,
        parent=args.parent,
        summary=args.summary,
    )
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
