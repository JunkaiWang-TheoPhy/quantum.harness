#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import exact_d6_density, exact_log_e5_density
from trottercert.hpc_artifacts import CoordinateCubicTerms, coordinate_decode_terms, coordinate_encode_terms, coordinate_terms_to_json, read_shard_gzip, sha256_file, write_manifest_atomic, write_shard_gzip
from trottercert.intervals import cube_root_four_interval
from trottercert.local_commutators import CoordinateRegistry


def _merge(target: CoordinateCubicTerms, raw_terms: object) -> None:
    for key, coefficient in coordinate_decode_terms(raw_terms).items():
        updated = target.get(key, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            target.pop(key, None)
        else:
            target[key] = updated


def _symplectic(registry: CoordinateRegistry, terms: CoordinateCubicTerms):
    result = {}
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    paths = []
    provenance = []
    commits = set()
    forward: CoordinateCubicTerms = {}
    expected_first = 0
    for index in range(args.shard_count):
        cell = args.input_root / f"word-{index:02d}"
        payload_path = cell / "shard.json.gz"
        manifest_path = cell / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("status") != "complete" or manifest.get("output_sha256") != sha256_file(payload_path):
            raise ValueError(f"word shard {index} is incomplete or corrupt")
        payload = read_shard_gzip(payload_path)
        if payload.get("kind") != "issue128_exact_e7_word_shard" or payload.get("shard_index") != index or payload.get("shard_count") != args.shard_count:
            raise ValueError(f"word shard {index} metadata mismatch")
        if payload.get("word_first") != expected_first:
            raise ValueError("word-shard coverage has a gap or overlap")
        expected_first = int(payload["word_last"])
        _merge(forward, payload["terms"])
        commits.add(str(payload["source_commit"]))
        paths.append(payload_path)
        provenance.append({"shard_index": index, "payload_sha256": sha256_file(payload_path), "manifest_sha256": sha256_file(manifest_path)})
    if len(commits) != 1 or expected_first != int(read_shard_gzip(paths[-1])["word_total"]):
        raise ValueError("word-shard source or coverage mismatch")
    reverse: CoordinateCubicTerms = {}
    for path in reversed(paths):
        _merge(reverse, read_shard_gzip(path)["terms"])
    if forward != reverse:
        raise ArithmeticError("forward and reverse E7 reductions differ")
    registry = CoordinateRegistry()
    e7 = _symplectic(registry, forward)
    e5_registry, e5_raw = exact_log_e5_density(fourth_order_suzuki_cubic_stages())
    e5_coordinates = coordinate_decode_terms(coordinate_encode_terms(e5_registry, e5_raw))
    e5 = _symplectic(registry, e5_coordinates)
    d6 = exact_d6_density(registry, e5, e7)
    root = cube_root_four_interval(30)
    cell_l1 = sum((coefficient.enclose(root).abs_upper() for coefficient in d6.values()), Fraction())
    parent_payload = {"schema_version": 1, "kind": "issue128_e7_word_reduction_parent", "source_commit": next(iter(commits)), "shard_count": args.shard_count, "reducer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "shards": provenance, "e7_term_count": len(forward), "identity": "D6=7*E7+(2/3)*ad_A^2(E5)"}
    write_manifest_atomic(args.parent, parent_payload)
    payload = {"schema_version": 1, "kind": "issue128_exact_right_generator_degree", "formula_id": "five_copy_suzuki_fourth_order_exact_cubic", "source_commit": next(iter(commits)), "stage_count": 31, "degree": 6, "coefficient_interval_decimal_digits": 30, "term_count": len(d6), "terms": coordinate_encode_terms(registry, d6), "cell_pauli_l1_upper": _pair(cell_l1), "site_pauli_l1_upper": _pair(cell_l1 / 4), "parent_sha256": sha256_file(args.parent)}
    write_shard_gzip(args.output, payload)
    summary = {"schema_version": 1, "kind": "issue128_e7_word_reduction_summary", "status": "complete", "source_commit": next(iter(commits)), "e7_term_count": len(forward), "d6_term_count": len(d6), "cell_pauli_l1_upper": _pair(cell_l1), "site_pauli_l1_upper": _pair(cell_l1 / 4), "parent_sha256": sha256_file(args.parent), "output_sha256": sha256_file(args.output), "reduction_order_check": "forward_equals_reverse"}
    write_manifest_atomic(args.summary, summary)
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
