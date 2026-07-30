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

from trottercert.hpc_artifacts import (
    coordinate_decode_terms,
    coordinate_terms_to_json,
    read_shard_gzip,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)
from trottercert.intervals import cube_root_four_interval


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _max_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--decimal-digits", type=int, default=30)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    arguments = parser.parse_args()
    started = time.perf_counter()
    extractor_commit = os.environ.get("ISSUE128_EXTRACTOR_COMMIT", "unrecorded")
    running = {
        "schema_version": 1,
        "kind": "issue128_exact_degree_extract_manifest",
        "status": "running",
        "degree": arguments.degree,
        "extractor_commit": extractor_commit,
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    write_manifest_atomic(arguments.manifest, running)
    try:
        parent = read_shard_gzip(arguments.input)
        raw_series = parent.get("series")
        if (
            parent.get("kind")
            != "issue128_exact_right_generator_monolithic"
            or not isinstance(raw_series, list)
            or arguments.degree < 0
            or arguments.degree >= len(raw_series)
        ):
            raise ValueError("parent artifact does not cover the requested degree")
        raw_terms = raw_series[arguments.degree]
        terms = coordinate_decode_terms(raw_terms)
        if coordinate_terms_to_json(terms) != raw_terms:
            raise ValueError("parent degree terms are not canonically sorted")
        root = cube_root_four_interval(arguments.decimal_digits)
        cell_l1 = sum(
            (coefficient.enclose(root).abs_upper() for coefficient in terms.values()),
            Fraction(),
        )
        payload = {
            "schema_version": 1,
            "kind": "issue128_exact_right_generator_degree",
            "formula_id": parent.get("formula_id"),
            "source_commit": parent.get("source_commit"),
            "stage_count": parent.get("stage_count"),
            "degree": arguments.degree,
            "coefficient_interval_decimal_digits": arguments.decimal_digits,
            "term_count": len(terms),
            "terms": raw_terms,
            "cell_pauli_l1_upper": _pair(cell_l1),
            "site_pauli_l1_upper": _pair(cell_l1 / 4),
            "parent_sha256": sha256_file(arguments.input),
        }
        write_shard_gzip(arguments.output, payload)
        complete = {
            **running,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": _max_rss_bytes(),
            "source_commit": parent.get("source_commit"),
            "term_count": len(terms),
            "cell_pauli_l1_upper": _pair(cell_l1),
            "site_pauli_l1_upper": _pair(cell_l1 / 4),
            "parent_sha256": payload["parent_sha256"],
            "output_sha256": sha256_file(arguments.output),
        }
        write_manifest_atomic(arguments.manifest, complete)
        print(json.dumps(complete, sort_keys=True), flush=True)
    except BaseException as error:
        write_manifest_atomic(
            arguments.manifest,
            {
                **running,
                "status": "failed",
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "wall_seconds": time.perf_counter() - started,
                "peak_rss_bytes": _max_rss_bytes(),
                "error_type": type(error).__name__,
                "error": str(error),
            },
        )
        raise


if __name__ == "__main__":
    main()
