#!/usr/bin/env python3
"""Generate and reduce exact dual-only E7 pairing shards."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.dual_log_pairing import (
    DualPairingPartial,
    contract_log_degree_shard,
    suffix_group_ordinals,
)

ISSUE_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/cubic_local.py",
    ISSUE_ROOT / "src/trottercert/local_commutators.py",
    ISSUE_ROOT / "src/trottercert/commutant_witness.py",
    ISSUE_ROOT / "src/trottercert/dual_log_pairing.py",
)


def _rational_json(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [_rational_json(value.a0), _rational_json(value.a1), _rational_json(value.a2)]


def _parse_rational(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(not isinstance(entry, int) or isinstance(entry, bool) for entry in value)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if _rational_json(result) != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _parse_cubic(value: object, field: str) -> Cubic:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must contain three cubic coordinates")
    return Cubic(*(_parse_rational(entry, field) for entry in value))


def _canonical_digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(ISSUE_ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in SOURCE_PATHS
    }


def build_shard_payload(
    partial: DualPairingPartial,
    implementation_sources: Mapping[str, str],
) -> dict[str, Any]:
    sources = dict(sorted(implementation_sources.items()))
    payload = {
        "schema_version": 1,
        "kind": "issue128_dual_log_pairing_shard",
        "degree": partial.degree,
        "length": partial.length,
        "shard_index": partial.shard_index,
        "shard_count": partial.shard_count,
        "total_groups": partial.total_groups,
        "group_indices": list(partial.group_indices),
        "word_count": partial.word_count,
        "nonzero_word_count": partial.nonzero_word_count,
        "retained_term_count": partial.retained_term_count,
        "pairings": {
            "tau_h": _cubic_json(partial.tau_h),
            "tau_h2": _cubic_json(partial.tau_h2),
            "tau_w": _cubic_json(partial.tau_w),
        },
        "implementation_sources": sources,
        "implementation_sources_digest": _canonical_digest(sources),
        "hpc_used": False,
    }
    verify_shard_payload(payload)
    return payload


def _verify_sources(payload: Mapping[str, object]) -> Mapping[str, str]:
    sources = payload.get("implementation_sources")
    if not isinstance(sources, Mapping) or not sources:
        raise ValueError("source manifest is missing")
    if any(
        not isinstance(path, str)
        or not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
        for path, digest in sources.items()
    ):
        raise ValueError("source manifest is malformed")
    if payload.get("implementation_sources_digest") != _canonical_digest(sources):
        raise ValueError("source digest mismatch")
    return sources


def verify_shard_payload(payload: Mapping[str, object]) -> None:
    if payload.get("kind") != "issue128_dual_log_pairing_shard":
        raise ValueError("unexpected dual-pairing shard kind")
    degree = payload.get("degree")
    length = payload.get("length")
    shard_index = payload.get("shard_index")
    shard_count = payload.get("shard_count")
    total_groups = payload.get("total_groups")
    if not isinstance(degree, int) or degree < 3 or degree % 2 == 0:
        raise ValueError("shard degree is invalid")
    if not isinstance(length, int) or length < 6 or length % 2:
        raise ValueError("shard torus length is invalid")
    if not all(isinstance(value, int) for value in (shard_index, shard_count, total_groups)):
        raise ValueError("shard metadata is malformed")
    expected_groups = suffix_group_ordinals(total_groups, shard_index, shard_count)
    if payload.get("group_indices") != list(expected_groups):
        raise ValueError("shard group assignment mismatch")
    for field in ("word_count", "nonzero_word_count", "retained_term_count"):
        value = payload.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"shard {field} is invalid")
    pairings = payload.get("pairings")
    if not isinstance(pairings, Mapping):
        raise ValueError("shard pairings are missing")
    tau_h = _parse_cubic(pairings.get("tau_h"), "tau_h pairing")
    tau_h2 = _parse_cubic(pairings.get("tau_h2"), "tau_h2 pairing")
    tau_w = _parse_cubic(pairings.get("tau_w"), "tau_w pairing")
    if tau_w != tau_h2 + tau_h / 2:
        raise ValueError("witness pairing mismatch")
    _verify_sources(payload)
    if payload.get("hpc_used") is not False:
        raise ValueError("dual-pairing shard must record hpc_used=false")


def reduce_shard_payloads(
    shards: Sequence[Mapping[str, object]],
    *,
    parent_sha256: Sequence[str],
) -> dict[str, Any]:
    if len(shards) != len(parent_sha256) or not shards:
        raise ValueError("one parent digest is required per shard")
    for shard in shards:
        verify_shard_payload(shard)
    first = shards[0]
    shard_count = first["shard_count"]
    if len(shards) != shard_count:
        raise ValueError("complete shard set is required")
    indices = [shard["shard_index"] for shard in shards]
    if sorted(indices) != list(range(shard_count)):
        raise ValueError("shard indices must be unique and complete")
    common_fields = (
        "degree",
        "length",
        "shard_count",
        "total_groups",
        "implementation_sources",
        "implementation_sources_digest",
    )
    for shard in shards[1:]:
        for field in common_fields:
            if shard[field] != first[field]:
                label = "source" if "source" in field else "shard"
                raise ValueError(f"{label} configuration mismatch")

    ordered = sorted(shards, key=lambda shard: shard["shard_index"])
    groups = [group for shard in ordered for group in shard["group_indices"]]
    if len(groups) != len(set(groups)) or set(groups) != set(range(first["total_groups"])):
        raise ValueError("shard group coverage is incomplete or overlapping")

    def sum_pairing(name: str, sequence) -> Cubic:
        total = Cubic.zero()
        for shard in sequence:
            total += _parse_cubic(shard["pairings"][name], name)
        return total

    forward = {
        name: sum_pairing(name, ordered) for name in ("tau_h", "tau_h2", "tau_w")
    }
    reverse = {
        name: sum_pairing(name, reversed(ordered))
        for name in ("tau_h", "tau_h2", "tau_w")
    }
    if forward != reverse:
        raise ArithmeticError("forward and reverse shard reductions differ")
    cells = first["length"] ** 2 // 4
    payload = {
        "schema_version": 1,
        "kind": "issue128_dual_e7_pairing",
        "degree": first["degree"],
        "length": first["length"],
        "shard_count": shard_count,
        "parents": [
            {
                "shard_index": shard["shard_index"],
                "sha256": digest,
            }
            for shard, digest in zip(ordered, parent_sha256)
        ],
        "coverage": {
            "group_count": len(groups),
            "word_count": sum(shard["word_count"] for shard in ordered),
            "nonzero_word_count": sum(
                shard["nonzero_word_count"] for shard in ordered
            ),
            "retained_term_count": sum(
                shard["retained_term_count"] for shard in ordered
            ),
            "reduction_order_check": "forward_equals_reverse",
        },
        "pairings": {
            **{name: _cubic_json(value) for name, value in forward.items()},
            "tau_w_per_cell": _cubic_json(forward["tau_w"] / cells),
            "tau_w_nonzero": forward["tau_w"] != Cubic.zero(),
        },
        "implementation_sources": first["implementation_sources"],
        "implementation_sources_digest": first["implementation_sources_digest"],
        "claim": {
            "full_e7_operator": "not_computed",
            "dual_e7_pairing": "exact",
            "finite_step_status": "inconclusive",
            "missing": "E9-and-higher dual remainder and logarithm branch",
        },
        "hpc_authorized": False,
    }
    verify_reduced_payload(payload)
    return payload


def verify_reduced_payload(payload: Mapping[str, object]) -> None:
    if payload.get("kind") != "issue128_dual_e7_pairing":
        raise ValueError("unexpected reduced dual-pairing kind")
    coverage = payload.get("coverage")
    if not isinstance(coverage, Mapping) or coverage.get("reduction_order_check") != "forward_equals_reverse":
        raise ValueError("reduction coverage proof is missing")
    parents = payload.get("parents")
    shard_count = payload.get("shard_count")
    if not isinstance(parents, list) or len(parents) != shard_count:
        raise ValueError("reduced parent shard list mismatch")
    if [parent.get("shard_index") for parent in parents] != list(range(shard_count)):
        raise ValueError("reduced parent shard indices mismatch")
    if any(
        not isinstance(parent.get("sha256"), str) or len(parent["sha256"]) != 64
        for parent in parents
    ):
        raise ValueError("reduced parent digest is malformed")
    pairings = payload.get("pairings")
    if not isinstance(pairings, Mapping):
        raise ValueError("reduced pairings are missing")
    tau_h = _parse_cubic(pairings.get("tau_h"), "reduced tau_h")
    tau_h2 = _parse_cubic(pairings.get("tau_h2"), "reduced tau_h2")
    tau_w = _parse_cubic(pairings.get("tau_w"), "reduced tau_w")
    if tau_w != tau_h2 + tau_h / 2:
        raise ValueError("reduced witness pairing mismatch")
    cells = payload.get("length") ** 2 // 4
    if _parse_cubic(pairings.get("tau_w_per_cell"), "reduced per-cell pairing") != tau_w / cells:
        raise ValueError("reduced per-cell pairing mismatch")
    if pairings.get("tau_w_nonzero") != (tau_w != Cubic.zero()):
        raise ValueError("reduced nonzero pairing status mismatch")
    _verify_sources(payload)
    claim = payload.get("claim")
    if not isinstance(claim, Mapping) or claim.get("full_e7_operator") != "not_computed":
        raise ValueError("reduced full-operator claim mismatch")
    if claim.get("finite_step_status") != "inconclusive":
        raise ValueError("reduced finite-step claim exceeds evidence")
    if payload.get("hpc_authorized") is not False:
        raise ValueError("reduced artifact does not authorize HPC")


def _encoded(payload: Mapping[str, object]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _load(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"JSON artifact must be an object: {path}")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--shard-count", type=int)
    parser.add_argument("--reduce", type=Path, nargs="+")
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    modes = sum(
        (args.shard_index is not None, args.reduce is not None, args.verify is not None)
    )
    if modes != 1:
        raise SystemExit("choose exactly one of shard, reduce, or verify mode")
    if args.verify is not None:
        payload = _load(args.verify)
        if payload.get("kind") == "issue128_dual_log_pairing_shard":
            verify_shard_payload(payload)
        else:
            verify_reduced_payload(payload)
        print(f"dual E7 artifact valid: {args.verify}", flush=True)
        return
    if args.output is None:
        raise SystemExit("--output is required")

    if args.reduce is not None:
        shards = [_load(path) for path in args.reduce]
        digests = [hashlib.sha256(path.read_bytes()).hexdigest() for path in args.reduce]
        payload = reduce_shard_payloads(shards, parent_sha256=digests)
    else:
        if args.shard_count is None:
            raise SystemExit("--shard-count is required in shard mode")

        def progress(completed: int, total: int, partial: DualPairingPartial) -> None:
            if completed % 100 == 0 or completed == total:
                print(
                    f"shard {partial.shard_index}/{partial.shard_count}: "
                    f"groups={completed}/{total} words={partial.word_count} "
                    f"retained={partial.retained_term_count}",
                    flush=True,
                )

        partial = contract_log_degree_shard(
            fourth_order_suzuki_cubic_stages(4),
            7,
            args.shard_index,
            args.shard_count,
            progress=progress,
        )
        payload = build_shard_payload(partial, _source_hashes())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(_encoded(payload))
    print(f"wrote dual E7 artifact: {args.output}", flush=True)


if __name__ == "__main__":
    main()
