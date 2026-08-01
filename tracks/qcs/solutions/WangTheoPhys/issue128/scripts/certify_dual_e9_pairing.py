#!/usr/bin/env python3
"""Generate, reduce, and verify exact manifest-driven dual E9 pairings."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import sys
import time
from collections.abc import Mapping, Sequence
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path
from typing import Any

from scripts.certify_dual_e7_pairing import (
    verify_reduced_payload as verify_dual_e7_payload,
)
from scripts.certify_extensive_commutant_witness import verify_extensive_payload
from trottercert.cubic_field import Cubic
from trottercert.dual_manifest_pairing import (
    ManifestPairingPartial,
    contract_word_manifest,
)
from trottercert.dual_word_manifest import (
    ManifestIndex,
    ManifestShardRecord,
    WordManifest,
    load_manifest_index,
    load_word_manifest,
    verify_manifest_index,
)

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXTENSIVE = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
DEFAULT_DUAL_E7 = (
    ISSUE_ROOT / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
)
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/cubic_local.py",
    ISSUE_ROOT / "src/trottercert/local_commutators.py",
    ISSUE_ROOT / "src/trottercert/commutant_witness.py",
    ISSUE_ROOT / "src/trottercert/dual_word_manifest.py",
    ISSUE_ROOT / "src/trottercert/dual_manifest_pairing.py",
    ISSUE_ROOT / "scripts/certify_dual_e9_pairing.py",
)
E9_CONTRACTION_LENGTH = 12
E9_CONTRACTION_CELLS = E9_CONTRACTION_LENGTH**2 // 4


def _canonical_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def _digest(payload: object) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_digest(value: object, field: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{field} must be a lowercase SHA-256 digest")
    return value


def _rational_json(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _parse_rational(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(
            not isinstance(entry, int) or isinstance(entry, bool)
            for entry in value
        )
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if _rational_json(result) != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [
        _rational_json(value.a0),
        _rational_json(value.a1),
        _rational_json(value.a2),
    ]


def _parse_cubic(value: object, field: str) -> Cubic:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must contain three cubic coordinates")
    return Cubic(*(_parse_rational(entry, field) for entry in value))


def _source_tuple(
    implementation_sources: Mapping[str, str],
) -> tuple[tuple[str, str], ...]:
    if not implementation_sources:
        raise ValueError("implementation source manifest is empty")
    result = tuple(sorted(implementation_sources.items()))
    for path, digest in result:
        if not isinstance(path, str) or not path:
            raise ValueError("implementation source path is invalid")
        _valid_digest(digest, f"implementation source {path}")
    return result


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(ISSUE_ROOT)): _sha256_file(path)
        for path in SOURCE_PATHS
    }


def _index_payload(index: ManifestIndex) -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "issue128_dual_word_manifest_index",
        "degree": index.degree,
        "shard_count": index.shard_count,
        "total_groups": index.total_groups,
        "total_words": index.total_words,
        "word_set_digest": index.word_set_digest,
        "formula": index.formula,
        "implementation_sources": dict(index.implementation_sources),
        "shards": [
            {
                "shard_index": record.shard_index,
                "path": record.path,
                "sha256": record.sha256,
                "group_count": record.group_count,
                "word_count": record.word_count,
            }
            for record in index.shards
        ],
    }


def _index_digest(index: ManifestIndex) -> str:
    return _digest(_index_payload(index))


def _index_record(index: ManifestIndex, shard_index: int) -> ManifestShardRecord:
    matching = tuple(
        record for record in index.shards if record.shard_index == shard_index
    )
    if len(matching) != 1:
        raise ValueError("manifest index has no unique record for worker shard")
    return matching[0]


def _runtime_json(runtime: Mapping[str, object]) -> dict[str, object]:
    wall = runtime.get("wall_seconds")
    peak = runtime.get("peak_rss_bytes")
    scheduler = runtime.get("scheduler")
    if not isinstance(wall, str):
        raise ValueError("runtime wall_seconds must be a decimal string")
    try:
        parsed_wall = Decimal(wall)
    except InvalidOperation as exc:
        raise ValueError("runtime wall_seconds must be a decimal string") from exc
    if not parsed_wall.is_finite() or parsed_wall < 0:
        raise ValueError("runtime wall_seconds must be finite and nonnegative")
    if not isinstance(peak, int) or isinstance(peak, bool) or peak < 0:
        raise ValueError("runtime peak_rss_bytes must be nonnegative")
    if scheduler is not None:
        if not isinstance(scheduler, Mapping):
            raise ValueError("runtime scheduler metadata must be an object or null")
        allowed = {"job_id", "array_task_id", "cluster", "host"}
        if set(scheduler) - allowed or any(
            not isinstance(value, str) for value in scheduler.values()
        ):
            raise ValueError("runtime scheduler metadata is malformed")
        scheduler = dict(sorted(scheduler.items()))
    return {
        "wall_seconds": wall,
        "peak_rss_bytes": peak,
        "scheduler": scheduler,
    }


def _worker_mathematical_view(payload: Mapping[str, object]) -> dict[str, object]:
    fields = (
        "schema_version",
        "kind",
        "degree",
        "length",
        "shard_index",
        "shard_count",
        "total_groups",
        "group_indices",
        "word_count",
        "nonzero_word_count",
        "retained_term_count",
        "pairings",
        "manifest_index_sha256",
        "manifest_sha256",
        "word_set_digest",
        "implementation_sources",
        "implementation_sources_digest",
        "claim",
    )
    return {field: payload.get(field) for field in fields}


def build_worker_payload(
    partial: ManifestPairingPartial,
    index: ManifestIndex,
    manifest: WordManifest,
    manifest_sha256: str,
    *,
    implementation_sources: Mapping[str, str],
    runtime: Mapping[str, object],
) -> dict[str, object]:
    if partial.degree != 9 or index.degree != 9 or manifest.degree != 9:
        raise ValueError("dual E9 worker requires degree nine")
    if partial.length != E9_CONTRACTION_LENGTH:
        raise ValueError("dual E9 worker requires the alias-free 12x12 torus")
    if (
        partial.shard_index != manifest.shard_index
        or partial.shard_index >= index.shard_count
        or partial.shard_count != index.shard_count
        or manifest.shard_count != index.shard_count
        or partial.total_groups != index.total_groups
        or manifest.total_groups != index.total_groups
    ):
        raise ValueError("worker and manifest shard configuration mismatch")
    record = _index_record(index, partial.shard_index)
    manifest_sha256 = _valid_digest(manifest_sha256, "manifest digest")
    if manifest_sha256 != record.sha256:
        raise ValueError("worker manifest digest does not match index")
    expected_groups = tuple(group.ordinal for group in manifest.groups)
    if partial.group_indices != expected_groups:
        raise ValueError("worker group coverage does not match manifest")
    expected_words = sum(len(group.records) for group in manifest.groups)
    if partial.word_count != expected_words or expected_words != record.word_count:
        raise ValueError("worker word count does not match manifest")
    sources = _source_tuple(implementation_sources)
    payload: dict[str, object] = {
        "schema_version": 1,
        "kind": "issue128_dual_e9_pairing_shard",
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
        "manifest_index_sha256": _index_digest(index),
        "manifest_sha256": manifest_sha256,
        "word_set_digest": index.word_set_digest,
        "implementation_sources": dict(sources),
        "implementation_sources_digest": _digest(dict(sources)),
        "claim": {
            "full_e9_operator": "not_computed",
            "dual_e9_pairing": "exact",
            "finite_step_status": "inconclusive",
        },
        "runtime": _runtime_json(runtime),
    }
    payload["mathematical_payload_sha256"] = _digest(
        _worker_mathematical_view(payload)
    )
    verify_worker_payload(payload, index, manifest)
    return payload


def verify_worker_payload(
    payload: Mapping[str, object],
    index: ManifestIndex,
    manifest: WordManifest,
) -> None:
    if payload.get("schema_version") != 1:
        raise ValueError("worker schema version mismatch")
    if payload.get("kind") != "issue128_dual_e9_pairing_shard":
        raise ValueError("unexpected dual E9 worker artifact kind")
    if payload.get("degree") != 9 or manifest.degree != 9 or index.degree != 9:
        raise ValueError("dual E9 worker degree mismatch")
    shard_index = payload.get("shard_index")
    if not isinstance(shard_index, int) or isinstance(shard_index, bool):
        raise ValueError("worker shard index is invalid")
    record = _index_record(index, shard_index)
    if (
        payload.get("length") != E9_CONTRACTION_LENGTH
        or payload.get("shard_count") != index.shard_count
        or payload.get("total_groups") != index.total_groups
        or manifest.shard_index != shard_index
        or manifest.shard_count != index.shard_count
        or manifest.total_groups != index.total_groups
        or manifest.total_words != index.total_words
        or manifest.word_set_digest != index.word_set_digest
        or manifest.formula != index.formula
        or manifest.implementation_sources != index.implementation_sources
    ):
        raise ValueError("worker shard configuration mismatch")
    expected_groups = [group.ordinal for group in manifest.groups]
    if payload.get("group_indices") != expected_groups:
        raise ValueError("worker group coverage does not match manifest")
    expected_words = sum(len(group.records) for group in manifest.groups)
    if payload.get("word_count") != expected_words or record.word_count != expected_words:
        raise ValueError("worker word count does not match manifest")
    for field in ("nonzero_word_count", "retained_term_count"):
        value = payload.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"worker {field} is invalid")
    if payload["nonzero_word_count"] > payload["word_count"]:
        raise ValueError("worker nonzero word count exceeds total words")
    pairings = payload.get("pairings")
    if not isinstance(pairings, Mapping):
        raise ValueError("worker pairings are missing")
    tau_h = _parse_cubic(pairings.get("tau_h"), "worker tau_h pairing")
    tau_h2 = _parse_cubic(pairings.get("tau_h2"), "worker tau_h2 pairing")
    tau_w = _parse_cubic(pairings.get("tau_w"), "worker tau_w pairing")
    if tau_w != tau_h2 + tau_h / 2:
        raise ValueError("worker witness pairing mismatch")
    if payload.get("manifest_index_sha256") != _index_digest(index):
        raise ValueError("worker manifest index digest mismatch")
    if payload.get("manifest_sha256") != record.sha256:
        raise ValueError("worker manifest digest mismatch")
    if payload.get("word_set_digest") != index.word_set_digest:
        raise ValueError("worker word-set digest mismatch")
    raw_sources = payload.get("implementation_sources")
    if not isinstance(raw_sources, Mapping):
        raise ValueError("worker implementation sources are missing")
    sources = _source_tuple(raw_sources)
    if payload.get("implementation_sources_digest") != _digest(dict(sources)):
        raise ValueError("worker implementation source digest mismatch")
    claim = payload.get("claim")
    if not isinstance(claim, Mapping):
        raise ValueError("worker claim is missing")
    if (
        claim.get("full_e9_operator") != "not_computed"
        or claim.get("dual_e9_pairing") != "exact"
    ):
        raise ValueError("worker dual E9 claim mismatch")
    if claim.get("finite_step_status") != "inconclusive":
        raise ValueError("worker finite-step claim exceeds evidence")
    runtime = payload.get("runtime")
    if not isinstance(runtime, Mapping):
        raise ValueError("worker runtime metadata is missing")
    _runtime_json(runtime)
    if payload.get("mathematical_payload_sha256") != _digest(
        _worker_mathematical_view(payload)
    ):
        raise ValueError("worker mathematical payload digest mismatch")


def _load_json(path: Path, field: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_bytes())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {field}: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{field} must contain a JSON object")
    return payload


def _sum_pairing(
    workers: Sequence[Mapping[str, object]], name: str
) -> Cubic:
    total = Cubic.zero()
    for worker in workers:
        pairings = worker["pairings"]
        total += _parse_cubic(pairings[name], f"worker {name}")
    return total


def _reduced_mathematical_view(payload: Mapping[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in payload.items()
        if key != "mathematical_payload_sha256"
    }


def reduce_worker_payloads(
    workers: Sequence[Mapping[str, object]],
    *,
    index: ManifestIndex,
    manifests: Sequence[WordManifest],
    parent_sha256: Sequence[str],
    extensive_path: Path = DEFAULT_EXTENSIVE,
    dual_e7_path: Path = DEFAULT_DUAL_E7,
) -> dict[str, object]:
    if index.degree != 9:
        raise ValueError("dual E9 reducer requires a degree-nine index")
    if (
        len(workers) != index.shard_count
        or len(manifests) != index.shard_count
        or len(parent_sha256) != index.shard_count
    ):
        raise ValueError("complete worker, manifest, and parent shard sets are required")
    manifest_by_index = {manifest.shard_index: manifest for manifest in manifests}
    if set(manifest_by_index) != set(range(index.shard_count)):
        raise ValueError("manifest shard set is incomplete or duplicated")
    worker_by_index: dict[int, Mapping[str, object]] = {}
    parent_by_index: dict[int, str] = {}
    for worker, parent_digest in zip(workers, parent_sha256):
        shard_index = worker.get("shard_index")
        if not isinstance(shard_index, int) or shard_index in worker_by_index:
            raise ValueError("worker shard set is malformed or duplicated")
        manifest = manifest_by_index.get(shard_index)
        if manifest is None:
            raise ValueError("worker has no matching manifest shard")
        verify_worker_payload(worker, index, manifest)
        worker_by_index[shard_index] = worker
        parent_by_index[shard_index] = _valid_digest(
            parent_digest, "worker parent digest"
        )
    if set(worker_by_index) != set(range(index.shard_count)):
        raise ValueError("worker shard set is incomplete")
    ordered_workers = tuple(worker_by_index[index] for index in range(index.shard_count))
    groups = [
        group
        for worker in ordered_workers
        for group in worker["group_indices"]
    ]
    if groups != list(range(index.total_groups)):
        if len(groups) != len(set(groups)) or set(groups) != set(range(index.total_groups)):
            raise ValueError("worker group coverage is incomplete or overlapping")
        groups = sorted(groups)
    forward = {
        name: _sum_pairing(ordered_workers, name)
        for name in ("tau_h", "tau_h2", "tau_w")
    }
    reverse = {
        name: _sum_pairing(tuple(reversed(ordered_workers)), name)
        for name in ("tau_h", "tau_h2", "tau_w")
    }
    if forward != reverse:
        raise ArithmeticError("forward and reverse E9 reductions differ")
    extensive = _load_json(extensive_path, "extensive witness artifact")
    dual_e7 = _load_json(dual_e7_path, "dual E7 artifact")
    verify_extensive_payload(extensive)
    verify_dual_e7_payload(dual_e7)
    q5 = _parse_cubic(
        extensive["stable_relations"]["tau_w_e5_per_cell"],
        "q5 per-cell pairing",
    )
    q7 = _parse_cubic(
        dual_e7["pairings"]["tau_w_per_cell"],
        "q7 per-cell pairing",
    )
    q9 = forward["tau_w"] / E9_CONTRACTION_CELLS
    q9_scaled = q9 / 97**8
    combined = q5 / 97**4 + q7 / 97**6 + q9_scaled
    sources = ordered_workers[0]["implementation_sources"]
    if any(worker["implementation_sources"] != sources for worker in ordered_workers):
        raise ValueError("worker implementation source manifests differ")
    parents = []
    for shard_index in range(index.shard_count):
        parents.append(
            {
                "shard_index": shard_index,
                "sha256": parent_by_index[shard_index],
                "manifest_sha256": index.shards[shard_index].sha256,
            }
        )
    payload: dict[str, object] = {
        "schema_version": 1,
        "kind": "issue128_dual_e9_pairing",
        "degree": 9,
        "length": E9_CONTRACTION_LENGTH,
        "shard_count": index.shard_count,
        "manifest_index_sha256": _index_digest(index),
        "word_set_digest": index.word_set_digest,
        "parents": parents,
        "inputs": {
            "extensive_witness_sha256": _sha256_file(extensive_path),
            "dual_e7_sha256": _sha256_file(dual_e7_path),
            "q5_per_cell": _cubic_json(q5),
            "q7_per_cell": _cubic_json(q7),
        },
        "coverage": {
            "group_count": len(groups),
            "word_count": sum(worker["word_count"] for worker in ordered_workers),
            "nonzero_word_count": sum(
                worker["nonzero_word_count"] for worker in ordered_workers
            ),
            "retained_term_count": sum(
                worker["retained_term_count"] for worker in ordered_workers
            ),
            "reduction_order_check": "forward_equals_reverse",
        },
        "pairings": {
            **{name: _cubic_json(value) for name, value in forward.items()},
            "tau_w_per_cell": _cubic_json(q9),
            "q9_over_97_pow_8": _cubic_json(q9_scaled),
            "exact_e5_e7_e9_per_cell_at_r97": _cubic_json(combined),
            "tau_w_nonzero": forward["tau_w"] != Cubic.zero(),
        },
        "implementation_sources": sources,
        "implementation_sources_digest": ordered_workers[0][
            "implementation_sources_digest"
        ],
        "claim": {
            "full_e9_operator": "not_computed",
            "dual_e9_pairing": "exact",
            "finite_step_status": "inconclusive",
            "missing": "E11-and-higher dual tail",
        },
    }
    payload["mathematical_payload_sha256"] = _digest(
        _reduced_mathematical_view(payload)
    )
    verify_reduced_payload(payload, index=index, manifests=manifests)
    return payload


def verify_reduced_payload(
    payload: Mapping[str, object],
    *,
    index: ManifestIndex,
    manifests: Sequence[WordManifest],
    extensive_path: Path = DEFAULT_EXTENSIVE,
    dual_e7_path: Path = DEFAULT_DUAL_E7,
) -> None:
    if payload.get("schema_version") != 1:
        raise ValueError("reduced E9 schema version mismatch")
    if payload.get("kind") != "issue128_dual_e9_pairing":
        raise ValueError("unexpected reduced dual E9 artifact kind")
    if (
        payload.get("degree") != 9
        or payload.get("length") != E9_CONTRACTION_LENGTH
        or payload.get("shard_count") != index.shard_count
        or index.degree != 9
    ):
        raise ValueError("reduced E9 configuration mismatch")
    if payload.get("manifest_index_sha256") != _index_digest(index):
        raise ValueError("reduced manifest index digest mismatch")
    if payload.get("word_set_digest") != index.word_set_digest:
        raise ValueError("reduced word-set digest mismatch")
    if {manifest.shard_index for manifest in manifests} != set(
        range(index.shard_count)
    ):
        raise ValueError("reduced verifier manifest set is incomplete")
    for manifest in manifests:
        record = _index_record(index, manifest.shard_index)
        if (
            manifest.degree != index.degree
            or manifest.shard_count != index.shard_count
            or manifest.total_groups != index.total_groups
            or manifest.total_words != index.total_words
            or manifest.word_set_digest != index.word_set_digest
            or manifest.formula != index.formula
            or manifest.implementation_sources != index.implementation_sources
            or len(manifest.groups) != record.group_count
            or sum(len(group.records) for group in manifest.groups)
            != record.word_count
        ):
            raise ValueError("reduced verifier manifest configuration mismatch")
    parents = payload.get("parents")
    if not isinstance(parents, list) or len(parents) != index.shard_count:
        raise ValueError("reduced worker parent list is incomplete")
    for shard_index, parent in enumerate(parents):
        if not isinstance(parent, Mapping) or parent.get("shard_index") != shard_index:
            raise ValueError("reduced worker parent indices are malformed")
        _valid_digest(parent.get("sha256"), "worker parent digest")
        if parent.get("manifest_sha256") != index.shards[shard_index].sha256:
            raise ValueError("reduced parent manifest digest mismatch")
    inputs = payload.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("reduced E9 inputs are missing")
    if inputs.get("extensive_witness_sha256") != _sha256_file(extensive_path):
        raise ValueError("reduced extensive witness digest mismatch")
    if inputs.get("dual_e7_sha256") != _sha256_file(dual_e7_path):
        raise ValueError("reduced dual E7 digest mismatch")
    extensive = _load_json(extensive_path, "extensive witness artifact")
    dual_e7 = _load_json(dual_e7_path, "dual E7 artifact")
    verify_extensive_payload(extensive)
    verify_dual_e7_payload(dual_e7)
    q5 = _parse_cubic(
        extensive["stable_relations"]["tau_w_e5_per_cell"], "q5"
    )
    q7 = _parse_cubic(dual_e7["pairings"]["tau_w_per_cell"], "q7")
    if _parse_cubic(inputs.get("q5_per_cell"), "submitted q5") != q5:
        raise ValueError("reduced q5 input mismatch")
    if _parse_cubic(inputs.get("q7_per_cell"), "submitted q7") != q7:
        raise ValueError("reduced q7 input mismatch")
    coverage = payload.get("coverage")
    if not isinstance(coverage, Mapping):
        raise ValueError("reduced E9 coverage is missing")
    if (
        coverage.get("group_count") != index.total_groups
        or coverage.get("word_count") != index.total_words
        or coverage.get("reduction_order_check") != "forward_equals_reverse"
    ):
        raise ValueError("reduced E9 coverage proof is incomplete")
    for field in ("nonzero_word_count", "retained_term_count"):
        value = coverage.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"reduced E9 {field} is invalid")
    pairings = payload.get("pairings")
    if not isinstance(pairings, Mapping):
        raise ValueError("reduced E9 pairings are missing")
    tau_h = _parse_cubic(pairings.get("tau_h"), "reduced tau_h")
    tau_h2 = _parse_cubic(pairings.get("tau_h2"), "reduced tau_h2")
    tau_w = _parse_cubic(pairings.get("tau_w"), "reduced tau_w")
    if tau_w != tau_h2 + tau_h / 2:
        raise ValueError("reduced witness pairing mismatch")
    q9 = tau_w / E9_CONTRACTION_CELLS
    if _parse_cubic(pairings.get("tau_w_per_cell"), "reduced q9") != q9:
        raise ValueError("reduced per-cell E9 pairing mismatch")
    if _parse_cubic(pairings.get("q9_over_97_pow_8"), "scaled q9") != q9 / 97**8:
        raise ValueError("reduced scaled E9 pairing mismatch")
    combined = q5 / 97**4 + q7 / 97**6 + q9 / 97**8
    if (
        _parse_cubic(
            pairings.get("exact_e5_e7_e9_per_cell_at_r97"),
            "combined E5 E7 E9 pairing",
        )
        != combined
    ):
        raise ValueError("reduced combined E5 E7 E9 pairing mismatch")
    if pairings.get("tau_w_nonzero") != (tau_w != Cubic.zero()):
        raise ValueError("reduced E9 nonzero status mismatch")
    raw_sources = payload.get("implementation_sources")
    if not isinstance(raw_sources, Mapping):
        raise ValueError("reduced implementation sources are missing")
    sources = _source_tuple(raw_sources)
    if payload.get("implementation_sources_digest") != _digest(dict(sources)):
        raise ValueError("reduced implementation source digest mismatch")
    claim = payload.get("claim")
    if not isinstance(claim, Mapping):
        raise ValueError("reduced E9 claim is missing")
    if (
        claim.get("full_e9_operator") != "not_computed"
        or claim.get("dual_e9_pairing") != "exact"
    ):
        raise ValueError("reduced dual E9 claim mismatch")
    if claim.get("finite_step_status") != "inconclusive":
        raise ValueError("reduced finite-step claim exceeds evidence")
    if claim.get("missing") != "E11-and-higher dual tail":
        raise ValueError("reduced missing-evidence claim mismatch")
    if payload.get("mathematical_payload_sha256") != _digest(
        _reduced_mathematical_view(payload)
    ):
        raise ValueError("reduced mathematical payload digest mismatch")


def _peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _scheduler_metadata() -> dict[str, str] | None:
    mapping = {
        "job_id": "SLURM_JOB_ID",
        "array_task_id": "SLURM_ARRAY_TASK_ID",
        "cluster": "SLURM_CLUSTER_NAME",
        "host": "HOSTNAME",
    }
    result = {
        field: os.environ[variable]
        for field, variable in mapping.items()
        if os.environ.get(variable)
    }
    return result or None


def _load_all_manifests(index_path: Path, index: ManifestIndex) -> tuple[WordManifest, ...]:
    verify_manifest_index(index, index_path.parent)
    return tuple(
        load_word_manifest(index_path.parent / record.path, index)
        for record in index.shards
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--reduce", type=Path, nargs="+")
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    modes = sum((args.manifest is not None, args.reduce is not None, args.verify is not None))
    if modes != 1:
        raise SystemExit("choose exactly one of worker, reduce, or verify mode")
    index = load_manifest_index(args.index)

    if args.verify is not None:
        payload = _load_json(args.verify, "dual E9 artifact")
        if payload.get("kind") == "issue128_dual_e9_pairing_shard":
            shard_index = payload.get("shard_index")
            if not isinstance(shard_index, int):
                raise ValueError("worker shard index is invalid")
            record = _index_record(index, shard_index)
            manifest = load_word_manifest(args.index.parent / record.path, index)
            verify_worker_payload(payload, index, manifest)
        else:
            manifests = _load_all_manifests(args.index, index)
            verify_reduced_payload(payload, index=index, manifests=manifests)
        print(f"dual E9 artifact valid: {args.verify}", flush=True)
        return

    if args.output is None:
        raise SystemExit("--output is required in worker and reduce modes")
    if args.manifest is not None:
        manifest = load_word_manifest(args.manifest, index)
        start = time.perf_counter()

        def progress(
            completed: int,
            total: int,
            partial: ManifestPairingPartial,
        ) -> None:
            if completed % 100 == 0 or completed == total:
                print(
                    f"E9 shard {partial.shard_index}/{partial.shard_count}: "
                    f"groups={completed}/{total} words={partial.word_count} "
                    f"retained={partial.retained_term_count}",
                    flush=True,
                )

        partial = contract_word_manifest(
            manifest,
            length=E9_CONTRACTION_LENGTH,
            progress=progress,
        )
        payload = build_worker_payload(
            partial,
            index,
            manifest,
            _sha256_file(args.manifest),
            implementation_sources=_source_hashes(),
            runtime={
                "wall_seconds": f"{time.perf_counter() - start:.9f}",
                "peak_rss_bytes": _peak_rss_bytes(),
                "scheduler": _scheduler_metadata(),
            },
        )
    else:
        manifests = _load_all_manifests(args.index, index)
        worker_payloads = tuple(
            _load_json(path, "dual E9 worker artifact") for path in args.reduce
        )
        payload = reduce_worker_payloads(
            worker_payloads,
            index=index,
            manifests=manifests,
            parent_sha256=tuple(_sha256_file(path) for path in args.reduce),
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(_canonical_bytes(payload))
    print(f"wrote dual E9 artifact: {args.output}", flush=True)


if __name__ == "__main__":
    main()
