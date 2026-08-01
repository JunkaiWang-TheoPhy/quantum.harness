#!/usr/bin/env python3
"""Audit a complete Issue-128 E9 run and its independent worker reruns."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from scripts.certify_dual_e9_pairing import (
    verify_reduced_payload,
    verify_worker_payload,
)
from trottercert.dual_word_manifest import (
    ManifestIndex,
    load_manifest_index,
    load_word_manifest,
    verify_manifest_index,
)


def _canonical_bytes(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path, field: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_bytes())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {field}: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{field} must contain a JSON object")
    return payload


def recommended_rerun_shards(index: ManifestIndex) -> tuple[int, int]:
    """Return deterministic median-load and maximum-load shard indices."""

    ordered = sorted(
        index.shards,
        key=lambda record: (record.word_count, record.group_count, record.shard_index),
    )
    median = ordered[(len(ordered) - 1) // 2].shard_index
    heavy = max(
        index.shards,
        key=lambda record: (record.word_count, record.group_count, -record.shard_index),
    ).shard_index
    if median == heavy:
        alternatives = [record.shard_index for record in reversed(ordered) if record.shard_index != heavy]
        if not alternatives:
            raise ValueError("independent rerun audit requires at least two shards")
        median = alternatives[0]
    return tuple(sorted((median, heavy)))


def _worker_paths(run_root: Path, shard_count: int) -> tuple[Path, ...]:
    shard_root = run_root / "shards"
    expected = tuple(shard_root / f"shard-{index:03d}.json" for index in range(shard_count))
    missing = [index for index, path in enumerate(expected) if not path.is_file()]
    if missing:
        raise ValueError(f"missing shard {missing[0]}")
    actual = sorted(shard_root.glob("shard-*.json"))
    if actual != list(expected):
        raise ValueError("worker shard directory contains a noncanonical or duplicate shard set")
    return expected


def audit_e9_run(
    run_root: Path,
    *,
    rerun_paths: Sequence[Path],
    expected_shards: int = 64,
) -> dict[str, Any]:
    """Verify manifests, workers, reduction, parent hashes, and rerun equality."""

    run_root = run_root.resolve()
    index_path = run_root / "manifests/index.json"
    index = load_manifest_index(index_path)
    if index.degree != 9 or index.shard_count != expected_shards:
        raise ValueError("E9 run has the wrong degree or shard count")
    verify_manifest_index(index, index_path.parent)
    manifests = tuple(
        load_word_manifest(index_path.parent / record.path, index)
        for record in index.shards
    )
    worker_paths = _worker_paths(run_root, index.shard_count)
    workers: list[dict[str, Any]] = []
    groups: list[int] = []
    for shard_index, (path, manifest) in enumerate(zip(worker_paths, manifests)):
        worker = _load(path, f"worker shard {shard_index}")
        verify_worker_payload(worker, index, manifest)
        if worker.get("shard_index") != shard_index:
            raise ValueError("worker filename and shard index mismatch")
        raw_groups = worker.get("group_indices")
        if not isinstance(raw_groups, list):
            raise ValueError("worker group coverage is missing")
        groups.extend(raw_groups)
        workers.append(worker)
    if len(groups) != len(set(groups)):
        raise ValueError("duplicate group ordinal")
    if sorted(groups) != list(range(index.total_groups)):
        raise ValueError("worker group coverage is incomplete")

    reduced_path = run_root / "dual-e9-pairing.json"
    reduced = _load(reduced_path, "reduced E9 artifact")
    verify_reduced_payload(reduced, index=index, manifests=manifests)
    parents = reduced.get("parents")
    if not isinstance(parents, list) or len(parents) != index.shard_count:
        raise ValueError("reduced worker parent list is incomplete")
    worker_sha256 = tuple(_sha(path) for path in worker_paths)
    if [parent.get("sha256") for parent in parents] != list(worker_sha256):
        raise ValueError("reduced worker parent digest does not match archived shards")

    required_reruns = recommended_rerun_shards(index)
    reruns: dict[int, dict[str, str]] = {}
    for path in rerun_paths:
        rerun = _load(path, "independent E9 rerun")
        shard_index = rerun.get("shard_index")
        if not isinstance(shard_index, int) or isinstance(shard_index, bool):
            raise ValueError("independent rerun shard index is invalid")
        if shard_index in reruns:
            raise ValueError("independent rerun shard index is duplicated")
        if not 0 <= shard_index < index.shard_count:
            raise ValueError("independent rerun shard index is out of range")
        verify_worker_payload(rerun, index, manifests[shard_index])
        production_digest = workers[shard_index].get("mathematical_payload_sha256")
        rerun_digest = rerun.get("mathematical_payload_sha256")
        if rerun_digest != production_digest:
            raise ValueError("independent rerun mathematical payload mismatch")
        reruns[shard_index] = {
            "path": str(path.resolve()),
            "sha256": _sha(path),
            "mathematical_payload_sha256": str(rerun_digest),
        }
    if tuple(sorted(reruns)) != required_reruns:
        raise ValueError(
            "independent reruns must cover the deterministic median and heavy shards "
            f"{required_reruns}"
        )

    return {
        "schema_version": 1,
        "kind": "issue128_dual_e9_run_audit",
        "run_root": str(run_root),
        "degree": index.degree,
        "shard_count": index.shard_count,
        "coverage": {
            "group_count": len(groups),
            "word_count": sum(worker["word_count"] for worker in workers),
            "group_ordinals": "complete_unique_0_through_total_minus_1",
        },
        "manifest_index": {
            "path": str(index_path),
            "sha256": _sha(index_path),
            "word_set_digest": index.word_set_digest,
        },
        "workers": [
            {
                "shard_index": shard_index,
                "path": str(path),
                "sha256": worker_sha256[shard_index],
                "mathematical_payload_sha256": workers[shard_index][
                    "mathematical_payload_sha256"
                ],
            }
            for shard_index, path in enumerate(worker_paths)
        ],
        "reduced": {
            "path": str(reduced_path),
            "sha256": _sha(reduced_path),
            "mathematical_payload_sha256": reduced["mathematical_payload_sha256"],
        },
        "independent_reruns": [
            {"shard_index": shard_index, **reruns[shard_index]}
            for shard_index in required_reruns
        ],
        "claim": {
            "complete_manifest_bound_coverage": True,
            "reduction_verified": True,
            "independent_reruns_match": True,
            "finite_step_status": "inconclusive",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--rerun", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit_e9_run(args.run_root, rerun_paths=args.rerun)
    encoded = _canonical_bytes(payload)
    if args.output is None:
        print(encoded.decode(), end="")
        return
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing audit: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(encoded)
    print(f"wrote E9 run audit: {args.output}", flush=True)


if __name__ == "__main__":
    main()
