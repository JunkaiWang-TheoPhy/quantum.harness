#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections.abc import Mapping
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from trottercert.d6_physical_channels import build_grouped_d6_bound
from trottercert.exact_series_certificate import (
    read_portable_canonical_gzip,
    verify_exact_degree_payload,
)
from trottercert.hpc_artifacts import write_manifest_atomic


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "issue128-d6-exact.json.gz"
DEFAULT_OUTPUT = ROOT / "certificates" / "issue128-d6-groups.json"


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _fraction_pair(value: object, *, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or isinstance(value[0], bool)
        or isinstance(value[1], bool)
        or not isinstance(value[0], int)
        or not isinstance(value[1], int)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a rational pair")
    return Fraction(value[0], value[1])


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _verified_source(path: Path):
    raw = path.read_bytes()
    payload = read_portable_canonical_gzip(path)
    source_commit = payload.get("source_commit")
    if not isinstance(source_commit, str) or not source_commit:
        raise ValueError("exact D6 source commit is missing")
    verified = verify_exact_degree_payload(
        payload,
        expected_degree=6,
        expected_source_commit=source_commit,
    )
    return raw, payload, verified


def make_sidecar_payload(
    source_path: str | Path,
    *,
    candidate_cap: int = 128,
) -> dict[str, object]:
    path = Path(source_path)
    raw, source, verified = _verified_source(path)
    decimal_digits = int(source["coefficient_interval_decimal_digits"])
    grouped = build_grouped_d6_bound(
        verified.terms,
        decimal_digits,
        candidate_cap=candidate_cap,
    )
    group_bounds = grouped.group_bounds
    if sum(group_bounds, Fraction()) != grouped.grouped_cell_bound:
        raise ValueError("grouped D6 per-group bounds do not sum to total")
    if grouped.l1_cell_bound != verified.cell_l1_upper:
        raise ValueError("grouped D6 l1 baseline differs from exact source")
    if grouped.grouped_cell_bound > grouped.l1_cell_bound:
        raise ValueError("grouped D6 bound is weaker than its l1 baseline")

    return {
        "schema_version": 1,
        "kind": "issue128_d6_physical_channel_groups",
        "source_payload_sha256": _sha256(raw),
        "source_commit": verified.source_commit,
        "degree": 6,
        "coefficient_interval_decimal_digits": decimal_digits,
        "candidate_cap": candidate_cap,
        "term_count": grouped.term_count,
        "group_count": len(grouped.groups),
        "max_group_size": grouped.max_group_size,
        "grouping_method": "local_overlap_pairs_exact_symplectic_recertification",
        "groups": [
            {
                "term_indices": list(group),
                "bound": _pair(bound),
            }
            for group, bound in zip(grouped.groups, group_bounds)
        ],
        "channel_counts": {
            name: count for name, count in grouped.channel_counts
        },
        "channel_l1_bounds": {
            name: _pair(bound) for name, bound in grouped.channel_l1_bounds
        },
        "cell_pauli_l1_upper": _pair(grouped.l1_cell_bound),
        "site_pauli_l1_upper": _pair(grouped.l1_cell_bound / 4),
        "grouped_cell_bound": _pair(grouped.grouped_cell_bound),
        "grouped_site_bound": _pair(grouped.grouped_site_bound),
    }


def verify_sidecar_payload(
    payload: Mapping[str, object],
    source_path: str | Path,
) -> None:
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported grouped D6 sidecar schema")
    if payload.get("kind") != "issue128_d6_physical_channel_groups":
        raise ValueError("unexpected grouped D6 sidecar kind")
    candidate_cap = payload.get("candidate_cap")
    if (
        isinstance(candidate_cap, bool)
        or not isinstance(candidate_cap, int)
        or candidate_cap < 1
    ):
        raise ValueError("grouped D6 candidate cap is invalid")
    expected = make_sidecar_payload(source_path, candidate_cap=candidate_cap)
    if dict(payload) != expected:
        raise ValueError("grouped D6 sidecar differs from exact regeneration")

    # Give malformed coverage a specific fail-closed check even though exact
    # regeneration above is already sufficient to reject it.
    groups = payload.get("groups")
    if not isinstance(groups, list):
        raise ValueError("grouped D6 groups must be a list")
    indices: list[int] = []
    for group in groups:
        if not isinstance(group, dict) or not isinstance(
            group.get("term_indices"), list
        ):
            raise ValueError("grouped D6 group is malformed")
        indices.extend(group["term_indices"])
        _fraction_pair(group.get("bound"), field="grouped D6 group bound")
    term_count = payload.get("term_count")
    if not isinstance(term_count, int) or sorted(indices) != list(range(term_count)):
        raise ValueError("grouped D6 partition coverage mismatch")


def build_sidecar(
    input_path: str | Path,
    output_path: str | Path,
    *,
    candidate_cap: int = 128,
) -> dict[str, object]:
    payload = make_sidecar_payload(input_path, candidate_cap=candidate_cap)
    verify_sidecar_payload(payload, input_path)
    # write_manifest_atomic emits sorted, compact canonical JSON and replaces
    # a temporary sibling atomically after fsync.
    write_manifest_atomic(output_path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--candidate-cap", type=int, default=128)
    arguments = parser.parse_args()
    payload = build_sidecar(
        arguments.input,
        arguments.output,
        candidate_cap=arguments.candidate_cap,
    )
    print(
        json.dumps(
            {
                "output": str(arguments.output),
                "source_payload_sha256": payload["source_payload_sha256"],
                "term_count": payload["term_count"],
                "group_count": payload["group_count"],
                "max_group_size": payload["max_group_size"],
                "site_pauli_l1_upper": payload["site_pauli_l1_upper"],
                "grouped_site_bound": payload["grouped_site_bound"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
