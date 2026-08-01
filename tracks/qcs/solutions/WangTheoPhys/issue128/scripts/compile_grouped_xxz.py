"""Validate or profile the frozen grouped-XXZ setup.

Tasks 1--3 expose only a bounded profiling skeleton.  Scientific artifact
builds intentionally fail closed until the ledger and independent verifier
are completed in Tasks 4--8.
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction

from trottercert.grouped_xxz import (
    XXZCompileSpec,
    canonical_json_bytes,
    fraction_pair,
    verify_xxz_suzuki_schedule,
    xxz_suzuki_schedule,
)


def _exact_delta(raw: str) -> Fraction:
    parts = raw.split("/")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("delta must use canonical NUM/DEN syntax")
    try:
        numerator, denominator = (int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "delta must use canonical NUM/DEN syntax"
        ) from exc
    if denominator <= 0:
        raise argparse.ArgumentTypeError("delta denominator must be positive")
    value = Fraction(numerator, denominator)
    if raw != f"{value.numerator}/{value.denominator}":
        raise argparse.ArgumentTypeError("delta fraction must be reduced and canonical")
    return value


def _positive_integer(raw: str) -> int:
    try:
        value = int(raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be a positive integer") from exc
    if value < 1 or raw != str(value):
        raise argparse.ArgumentTypeError("value must be a canonical positive integer")
    return value


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--profile-only", action="store_true")
    mode.add_argument("--build", action="store_true")
    parser.add_argument("--delta", required=True, type=_exact_delta)
    parser.add_argument(
        "--pipeline",
        required=True,
        choices=("direct-theorem", "local-log"),
    )
    parser.add_argument("--max-records", type=_positive_integer)
    parser.add_argument("--summary")
    parser.add_argument("--witness")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    arguments = parser.parse_args(argv)
    if arguments.build:
        if arguments.max_records is not None:
            parser.error("--max-records is valid only with --profile-only")
        if arguments.summary is None or arguments.witness is None:
            parser.error("--build requires explicit --summary and --witness paths")
        parser.error("scientific artifact builds remain fail-closed until Task 8")
    if arguments.summary is not None or arguments.witness is not None:
        parser.error("profile mode does not accept artifact output paths")
    if arguments.pipeline == "local-log":
        parser.error("the local-log profiling pipeline remains fail-closed until Task 9")

    spec = XXZCompileSpec.pilot(arguments.delta)
    schedule = xxz_suzuki_schedule()
    verify_xxz_suzuki_schedule(schedule)
    payload = {
        "delta": fraction_pair(spec.delta),
        "max_records": arguments.max_records,
        "mode": "profile-only",
        "pipeline": arguments.pipeline,
        "setup": "validated",
        "stage_count": len(schedule),
    }
    sys.stdout.buffer.write(canonical_json_bytes(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
