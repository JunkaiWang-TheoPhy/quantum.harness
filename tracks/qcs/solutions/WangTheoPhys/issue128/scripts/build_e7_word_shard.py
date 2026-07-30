#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
import os
from pathlib import Path
import resource
import sys
import time

from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import CubicTerms, _add_term, cubic_formula_log_series
from trottercert.hpc_artifacts import coordinate_encode_terms, sha256_file, write_manifest_atomic, write_shard_gzip
from trottercert.local_commutators import SymplecticDyadicLocalDensityEvaluator
from trottercert.refined_error import canonicalize_symplectic_unit_cell


def _max_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("word-shard index is out of range")
    started = time.perf_counter()
    source_commit = os.environ.get("ISSUE128_SOURCE_COMMIT", "unrecorded")
    running = {
        "schema_version": 1,
        "kind": "issue128_exact_e7_word_shard_manifest",
        "status": "running",
        "source_commit": source_commit,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    write_manifest_atomic(args.manifest, running)
    try:
        stages = fourth_order_suzuki_cubic_stages()
        logarithm = cubic_formula_log_series(stages, 7)
        words = list(logarithm[7].items())
        first = len(words) * args.shard_index // args.shard_count
        last = len(words) * (args.shard_index + 1) // args.shard_count
        evaluator = SymplecticDyadicLocalDensityEvaluator(shared_coordinates=True)
        registry = evaluator.registries[0]
        denominator = 7 * (1 << evaluator.denominator_exponent((0,) * 7))
        result: CubicTerms = {}
        for offset, (word, word_coefficient) in enumerate(words[first:last], start=1):
            for raw_pauli, numerator in evaluator.evaluate(word).items():
                pauli = canonicalize_symplectic_unit_cell(registry, raw_pauli)
                _add_term(result, pauli, word_coefficient * Fraction(numerator, denominator))
            evaluator.cache.pop(word, None)
            if offset % 128 == 0 or first + offset == last:
                print(f"shard={args.shard_index} words={offset}/{last-first} terms={len(result)} wall={time.perf_counter()-started:.3f}s", flush=True)
        payload = {
            "schema_version": 1,
            "kind": "issue128_exact_e7_word_shard",
            "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
            "source_commit": source_commit,
            "stage_count": len(stages),
            "degree": 7,
            "shard_index": args.shard_index,
            "shard_count": args.shard_count,
            "word_first": first,
            "word_last": last,
            "word_total": len(words),
            "terms": coordinate_encode_terms(registry, result),
        }
        write_shard_gzip(args.output, payload)
        complete = {
            **running,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": _max_rss_bytes(),
            "word_first": first,
            "word_last": last,
            "term_count": len(result),
            "output_sha256": sha256_file(args.output),
        }
        write_manifest_atomic(args.manifest, complete)
        print(json.dumps(complete, sort_keys=True), flush=True)
    except BaseException as error:
        write_manifest_atomic(args.manifest, {**running, "status": "failed", "error_type": type(error).__name__, "error": str(error)})
        raise


if __name__ == "__main__":
    main()
