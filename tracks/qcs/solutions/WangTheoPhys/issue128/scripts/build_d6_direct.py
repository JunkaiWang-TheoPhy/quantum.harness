#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import (
    exact_d6_density,
    exact_log_e5_e7_densities,
)
from trottercert.hpc_artifacts import (
    coordinate_encode_terms,
    sha256_file,
    write_manifest_atomic,
    write_shard_gzip,
)
from trottercert.intervals import cube_root_four_interval


FORMULA_ID = "five_copy_suzuki_fourth_order_exact_cubic"


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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    arguments = parser.parse_args()
    started = time.perf_counter()
    source_commit = os.environ.get("ISSUE128_SOURCE_COMMIT", "unrecorded")
    running = {
        "schema_version": 1,
        "kind": "issue128_direct_d6_manifest",
        "status": "running",
        "source_commit": source_commit,
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    write_manifest_atomic(arguments.manifest, running)

    def progress(degree: int, index: int, total: int) -> None:
        if index == total or index % 256 == 0:
            print(
                f"log_degree={degree} words={index}/{total} "
                f"wall={time.perf_counter()-started:.3f}s",
                flush=True,
            )

    try:
        stages = fourth_order_suzuki_cubic_stages()
        registry, e5, e7 = exact_log_e5_e7_densities(stages, progress=progress)
        d6 = exact_d6_density(registry, e5, e7)
        root = cube_root_four_interval(30)
        cell_l1 = sum(
            (coefficient.enclose(root).abs_upper() for coefficient in d6.values()),
            Fraction(),
        )
        parent_payload = {
            "schema_version": 1,
            "kind": "issue128_direct_d6_parent",
            "formula_id": FORMULA_ID,
            "source_commit": source_commit,
            "stage_count": len(stages),
            "algorithm_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "e5_term_count": len(e5),
            "e7_term_count": len(e7),
            "identity": "D6=7*E7+(2/3)*ad_A^2(E5)",
        }
        write_manifest_atomic(arguments.parent, parent_payload)
        payload = {
            "schema_version": 1,
            "kind": "issue128_exact_right_generator_degree",
            "formula_id": FORMULA_ID,
            "source_commit": source_commit,
            "stage_count": len(stages),
            "degree": 6,
            "coefficient_interval_decimal_digits": 30,
            "term_count": len(d6),
            "terms": coordinate_encode_terms(registry, d6),
            "cell_pauli_l1_upper": _pair(cell_l1),
            "site_pauli_l1_upper": _pair(cell_l1 / 4),
            "parent_sha256": sha256_file(arguments.parent),
        }
        write_shard_gzip(arguments.output, payload)
        complete = {
            **running,
            "status": "complete",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": _max_rss_bytes(),
            "e5_term_count": len(e5),
            "e7_term_count": len(e7),
            "d6_term_count": len(d6),
            "cell_pauli_l1_upper": _pair(cell_l1),
            "cell_pauli_l1_upper_decimal": _decimal(cell_l1),
            "site_pauli_l1_upper": _pair(cell_l1 / 4),
            "site_pauli_l1_upper_decimal": _decimal(cell_l1 / 4),
            "parent": arguments.parent.name,
            "parent_sha256": sha256_file(arguments.parent),
            "output": arguments.output.name,
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
