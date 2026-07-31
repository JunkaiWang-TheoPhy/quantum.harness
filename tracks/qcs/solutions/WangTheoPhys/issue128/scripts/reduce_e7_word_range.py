#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

from trottercert.cubic_field import Cubic
from trottercert.hpc_artifacts import (
    CoordinateCubicTerms,
    coordinate_decode_terms,
    coordinate_terms_to_json,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)


def _merge(target: CoordinateCubicTerms, raw_terms: object) -> None:
    for key, coefficient in coordinate_decode_terms(raw_terms).items():
        updated = target.get(key, Cubic.zero()) + coefficient
        if updated == Cubic.zero():
            target.pop(key, None)
        else:
            target[key] = updated


def reduce_word_range(
    *,
    input_root: Path,
    shard_count: int,
    shard_start: int,
    shard_stop: int,
    output: Path,
    manifest: Path,
) -> dict[str, object]:
    if not 0 <= shard_start < shard_stop <= shard_count:
        raise ValueError("invalid source shard range")
    started = time.perf_counter()
    running: dict[str, object] = {
        "schema_version": 1,
        "kind": "issue128_exact_e7_partial_manifest",
        "status": "running",
        "shard_count": shard_count,
        "shard_start": shard_start,
        "shard_stop": shard_stop,
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    write_manifest_atomic(manifest, running)
    try:
        terms: CoordinateCubicTerms = {}
        sources: list[dict[str, object]] = []
        metadata: tuple[str, str, int, int, int] | None = None
        expected_word_first: int | None = None
        word_first: int | None = None
        word_total: int | None = None
        for index in range(shard_start, shard_stop):
            cell = input_root / f"word-{index:02d}"
            payload_path = cell / "shard.json.gz"
            manifest_path = cell / "manifest.json"
            source_manifest = json.loads(manifest_path.read_text())
            if (
                source_manifest.get("status") != "complete"
                or source_manifest.get("output_sha256") != sha256_file(payload_path)
            ):
                raise ValueError(f"word shard {index} is incomplete or corrupt")
            payload = read_shard_gzip(payload_path)
            if (
                payload.get("kind") != "issue128_exact_e7_word_shard"
                or payload.get("shard_index") != index
                or payload.get("shard_count") != shard_count
            ):
                raise ValueError(f"word shard {index} metadata mismatch")
            current_metadata = (
                str(payload.get("formula_id")),
                str(payload.get("source_commit")),
                int(payload.get("stage_count", -1)),
                int(payload.get("degree", -1)),
                int(payload.get("word_total", -1)),
            )
            if metadata is None:
                metadata = current_metadata
                word_first = int(payload["word_first"])
                expected_word_first = word_first
                word_total = current_metadata[4]
            elif current_metadata != metadata:
                raise ValueError("word-shard source metadata mismatch")
            if int(payload["word_first"]) != expected_word_first:
                raise ValueError("word-shard coverage has a gap or overlap")
            expected_word_first = int(payload["word_last"])
            _merge(terms, payload["terms"])
            sources.append(
                {
                    "shard_index": index,
                    "payload_sha256": sha256_file(payload_path),
                    "manifest_sha256": sha256_file(manifest_path),
                }
            )
        assert metadata is not None
        assert word_first is not None
        assert expected_word_first is not None
        assert word_total is not None
        payload = {
            "schema_version": 1,
            "kind": "issue128_exact_e7_partial",
            "formula_id": metadata[0],
            "source_commit": metadata[1],
            "stage_count": metadata[2],
            "degree": metadata[3],
            "shard_count": shard_count,
            "shard_start": shard_start,
            "shard_stop": shard_stop,
            "word_first": word_first,
            "word_last": expected_word_first,
            "word_total": word_total,
            "sources": sources,
            "term_count": len(terms),
            "terms": coordinate_terms_to_json(terms),
        }
        write_shard_gzip(output, payload)
        complete = {
            **running,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - started,
            "source_commit": metadata[1],
            "word_first": word_first,
            "word_last": expected_word_first,
            "term_count": len(terms),
            "output_sha256": sha256_file(output),
        }
        write_manifest_atomic(manifest, complete)
        return complete
    except BaseException as error:
        write_manifest_atomic(
            manifest,
            {
                **running,
                "status": "failed",
                "error_type": type(error).__name__,
                "error": str(error),
            },
        )
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--shard-start", type=int, required=True)
    parser.add_argument("--shard-stop", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    result = reduce_word_range(
        input_root=args.input_root,
        shard_count=args.shard_count,
        shard_start=args.shard_start,
        shard_stop=args.shard_stop,
        output=args.output,
        manifest=args.manifest,
    )
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
