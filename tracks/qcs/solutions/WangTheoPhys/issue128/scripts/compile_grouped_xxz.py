"""Profile or build one deterministic compressed grouped-XXZ artifact."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from trottercert.grouped_xxz import (
    XXZCompileSpec,
    canonical_json_bytes,
    fraction_pair,
)
from trottercert.grouped_xxz_artifact import (
    build_merged_pilot_artifact,
    build_per_delta_artifact,
    source_closure,
    verify_merged_pilot_artifact,
    verify_per_delta_artifact,
    write_merged_pilot_artifact,
    write_per_delta_artifact,
)
from trottercert.grouped_xxz_compressed import (
    build_compressed_xxz_ledger,
    compile_compressed_grouped_xxz,
    verify_compressed_xxz_ledger,
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
    mode.add_argument(
        "--merge-pilot",
        nargs=4,
        metavar=(
            "LEFT_SUMMARY",
            "LEFT_WITNESS",
            "RIGHT_SUMMARY",
            "RIGHT_WITNESS",
        ),
    )
    parser.add_argument("--delta", type=_exact_delta)
    parser.add_argument(
        "--pipeline",
        choices=("direct-theorem", "local-log"),
    )
    parser.add_argument("--max-records", type=_positive_integer)
    parser.add_argument("--summary")
    parser.add_argument("--witness")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    arguments = parser.parse_args(argv)
    if arguments.merge_pilot is not None:
        if arguments.delta is not None or arguments.pipeline is not None:
            parser.error("--merge-pilot does not accept --delta or --pipeline")
        if arguments.max_records is not None:
            parser.error("--merge-pilot does not accept --max-records")
        if arguments.summary is None or arguments.witness is None:
            parser.error(
                "--merge-pilot requires explicit --summary and --witness paths"
            )
        summary_path = Path(arguments.summary)
        witness_path = Path(arguments.witness)
        if summary_path.exists() or witness_path.exists():
            parser.error("artifact output already exists; overwrite is forbidden")
        input_paths = tuple(Path(path) for path in arguments.merge_pilot)
        if any(not path.is_file() for path in input_paths):
            parser.error("--merge-pilot input artifact is missing")
        closure = source_closure(ROOT)
        artifact = build_merged_pilot_artifact(
            input_paths[0].read_bytes(),
            input_paths[1].read_bytes(),
            input_paths[2].read_bytes(),
            input_paths[3].read_bytes(),
            closure,
        )
        verify_merged_pilot_artifact(
            artifact.summary_bytes,
            artifact.witness_gzip_bytes,
            expected_source_closure=closure,
        )
        write_merged_pilot_artifact(
            artifact,
            summary_path,
            witness_path,
        )
        sys.stdout.buffer.write(artifact.summary_bytes)
        return 0

    if arguments.delta is None or arguments.pipeline is None:
        parser.error("--profile-only and --build require --delta and --pipeline")
    if arguments.build:
        if arguments.max_records is not None:
            parser.error("--max-records is valid only with --profile-only")
        if arguments.summary is None or arguments.witness is None:
            parser.error("--build requires explicit --summary and --witness paths")
    if arguments.pipeline == "local-log":
        parser.error(
            "the local-log profiling pipeline remains fail-closed until Task 9"
        )

    spec = XXZCompileSpec.pilot(arguments.delta)
    if arguments.build:
        summary_path = Path(arguments.summary)
        witness_path = Path(arguments.witness)
        if summary_path.exists() or witness_path.exists():
            parser.error("artifact output already exists; overwrite is forbidden")

        def progress(completed: int, projected: int) -> None:
            print(
                f"raw_records={completed}/{projected}",
                file=sys.stderr,
                flush=True,
            )

        ledger = build_compressed_xxz_ledger(spec, progress=progress)
        certificate = compile_compressed_grouped_xxz(spec, ledger=ledger)
        closure = source_closure(ROOT)
        artifact = build_per_delta_artifact(ledger, certificate, closure)
        verify_per_delta_artifact(
            artifact.summary_bytes,
            artifact.witness_gzip_bytes,
            expected_source_closure=closure,
        )
        write_per_delta_artifact(
            artifact,
            summary_path,
            witness_path,
        )
        sys.stdout.buffer.write(artifact.summary_bytes)
        return 0

    if arguments.summary is not None or arguments.witness is not None:
        parser.error("profile mode does not accept artifact output paths")
    limit = 200 if arguments.max_records is None else arguments.max_records
    ledger = build_compressed_xxz_ledger(spec, max_records=limit)
    verify_compressed_xxz_ledger(ledger)
    pairs = sum(
        len(group.terms) == 2
        for block in ledger.representative_blocks
        for group in block.groups
    )
    singletons = sum(
        len(group.terms) == 1
        for block in ledger.representative_blocks
        for group in block.groups
    )
    payload = {
        "active_words": len(ledger.word_weights),
        "complete": ledger.complete,
        "delta": fraction_pair(spec.delta),
        "grouped_constant": fraction_pair(ledger.grouped_constant),
        "grouped_triangle_ratio": fraction_pair(
            ledger.grouped_constant / ledger.triangle_constant
        ),
        "max_records": limit,
        "mode": "profile-only",
        "pair_groups": pairs,
        "pipeline": arguments.pipeline,
        "projected_raw_records": ledger.projected_record_count,
        "raw_records": ledger.raw_record_count,
        "representative_blocks": len(ledger.representative_blocks),
        "representative_terms": sum(
            len(block.terms) for block in ledger.representative_blocks
        ),
        "singleton_groups": singletons,
        "triangle_constant": fraction_pair(ledger.triangle_constant),
    }
    sys.stdout.buffer.write(canonical_json_bytes(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
