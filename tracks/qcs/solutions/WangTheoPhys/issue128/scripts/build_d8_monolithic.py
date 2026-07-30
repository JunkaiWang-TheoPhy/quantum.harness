#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import os
from pathlib import Path
import resource
import sys
import time

from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import exact_right_generator_local_series
from trottercert.hpc_artifacts import (
    coordinate_encode_series,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)
from trottercert.intervals import cube_root_four_interval


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _decimal(value: Fraction) -> str:
    with localcontext() as context:
        context.prec = 18
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def _max_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    started_at = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    commit = os.environ.get("ISSUE128_SOURCE_COMMIT", "unrecorded")
    running = {
        "schema_version": 1,
        "kind": "issue128_d8_monolithic_manifest",
        "status": "running",
        "order": args.order,
        "source_commit": commit,
        "started_at": started_at,
    }
    write_manifest_atomic(args.manifest, running)

    def progress(index: int, series: object) -> None:
        counts = [len(degree) for degree in series]  # type: ignore[arg-type]
        print(f"stage={index}/30 terms={counts} wall={time.perf_counter()-started:.3f}s", flush=True)

    try:
        stages = fourth_order_suzuki_cubic_stages()
        registry, series = exact_right_generator_local_series(
            stages, args.order, progress=progress
        )
        root = cube_root_four_interval(30)
        cell_l1 = sum(
            (coefficient.enclose(root).abs_upper() for coefficient in series[args.order].values()),
            Fraction(),
        )
        payload = {
            "schema_version": 1,
            "kind": "issue128_exact_right_generator_monolithic",
            "formula_id": "five_copy_suzuki_fourth_order_exact_cubic",
            "source_commit": commit,
            "stage_count": len(stages),
            "order": args.order,
            "series": coordinate_encode_series(registry, series),
            "cell_pauli_l1_upper": _pair(cell_l1),
            "site_pauli_l1_upper": _pair(cell_l1 / 4),
        }
        write_shard_gzip(args.output, payload)
        complete = {
            **running,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": _max_rss_bytes(),
            "degree_term_counts": [len(degree) for degree in series],
            "cell_pauli_l1_upper": _pair(cell_l1),
            "cell_pauli_l1_upper_decimal": _decimal(cell_l1),
            "site_pauli_l1_upper": _pair(cell_l1 / 4),
            "site_pauli_l1_upper_decimal": _decimal(cell_l1 / 4),
            "output": args.output.name,
            "output_sha256": sha256_file(args.output),
        }
        write_manifest_atomic(args.manifest, complete)
        print(json.dumps(complete, sort_keys=True), flush=True)
    except BaseException as error:
        write_manifest_atomic(
            args.manifest,
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
