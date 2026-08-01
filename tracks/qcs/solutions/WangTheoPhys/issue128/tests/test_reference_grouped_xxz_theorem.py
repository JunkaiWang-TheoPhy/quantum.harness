from __future__ import annotations

import ast
import copy
import hashlib
import json
import os
import subprocess
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import islice
from pathlib import Path

import pytest

import trottercert.grouped_xxz_compressed as compressed_module
from scripts.reference_grouped_xxz_theorem import (
    COMPRESSED_STREAM_DOMAIN,
    PROJECTED_RECORD_COUNT,
    CompressedRawTheoremSummary,
    Cubic,
    RawTheoremSummary,
    canonical_compressed_record_bytes,
    canonical_json_bytes,
    compressed_record_payload,
    decode_compressed_record_payload,
    decode_record_payload,
    iter_raw_theorem_records,
    load_canonical_json,
    record_payload,
    strict_fraction_pair,
    summarize_compressed_raw_theorem,
    summarize_raw_theorem,
    theorem_stages,
)
from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.grouped_xxz import (
    XXZCompileSpec,
    _iter_raw_theorem_records,
    _theorem_interval_stages,
)
from trottercert.grouped_xxz_compressed import (
    RepresentativeNormBlock,
    build_compressed_xxz_ledger,
    canonical_raw_record_bytes,
)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "scripts/reference_grouped_xxz_theorem.py"
PREFIX_DIGEST = "e57f4e2146c192658d89ae295d0bad5963e30ddf51c69f392dcfb794a3df378c"
FULL_DIGEST = "c89a7708ced8cc2465a9f1dd61747488ff43850d3e7a07c7787c1e5d501a8e24"
COMPRESSED_PREFIX_DIGEST = (
    "1ef82b794b735ba805c70d4bdf48b78d82b7784b3e9096800824f9758f84b912"
)
COMPRESSED_FULL_DIGEST = (
    "98a5a3c7bf404bce6dfa888f1ec3edeb715f4f134aa649487ccc5df8d9982d50"
)


def _primary_payload(record: object) -> dict[str, object]:
    weight = record.weight_interval.upper
    return {
        "record_id": record.record_id,
        "side": record.side,
        "partial_sum_index": record.partial_sum_index,
        "adjoint_stage_indices": list(record.adjoint_stage_indices),
        "composition": list(record.composition),
        "base_stage_index": record.base_stage_index,
        "fragment_word": list(record.fragment_word),
        "weight_upper": [weight.numerator, weight.denominator],
    }


def _primary_compressed_payload(record: object) -> dict[str, object]:
    interval = record.weight_interval
    return {
        "adjoint_stage_indices": list(record.adjoint_stage_indices),
        "base_stage_index": record.base_stage_index,
        "block_key": list(record.block_key),
        "composition": list(record.composition),
        "fragment_word": list(record.fragment_word),
        "partial_sum_index": record.partial_sum_index,
        "record_id": record.record_id,
        "side": record.side,
        "stream_domain": COMPRESSED_STREAM_DOMAIN,
        "weight_interval": [
            [interval.lower.numerator, interval.lower.denominator],
            [interval.upper.numerator, interval.upper.denominator],
        ],
    }


@lru_cache(maxsize=2)
def _primary_summaries(
    max_records: int | None,
) -> tuple[RawTheoremSummary, CompressedRawTheoremSummary]:
    legacy_digest = hashlib.sha256()
    compressed_digest = hashlib.sha256()
    weights: dict[tuple[int, ...], Fraction] = {}
    count = 0
    records = _iter_raw_theorem_records()
    selected = records if max_records is None else islice(records, max_records)
    for record in selected:
        legacy_digest.update(canonical_json_bytes(_primary_payload(record)))
        compressed_digest.update(canonical_raw_record_bytes(record))
        weights[record.fragment_word] = (
            weights.get(record.fragment_word, Fraction())
            + record.weight_interval.upper
        )
        count += 1
    legacy = RawTheoremSummary(
        record_count=count,
        projected_record_count=PROJECTED_RECORD_COUNT,
        complete=count == PROJECTED_RECORD_COUNT,
        stream_sha256=legacy_digest.hexdigest(),
        word_weights=tuple(sorted(weights.items())),
    )
    compressed = CompressedRawTheoremSummary(
        stream_domain=COMPRESSED_STREAM_DOMAIN,
        record_count=count,
        projected_record_count=PROJECTED_RECORD_COUNT,
        complete=count == PROJECTED_RECORD_COUNT,
        stream_sha256=compressed_digest.hexdigest(),
        word_weights=tuple(sorted(weights.items())),
    )
    return legacy, compressed


def _primary_summary(max_records: int | None) -> RawTheoremSummary:
    return _primary_summaries(max_records)[0]


def _primary_compressed_summary(
    max_records: int | None,
) -> CompressedRawTheoremSummary:
    return _primary_summaries(max_records)[1]


def test_reference_oracle_imports_only_python_standard_library() -> None:
    tree = ast.parse(REFERENCE.read_text(encoding="utf-8"))
    imported_roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_roots.add(node.module.split(".", 1)[0])

    assert imported_roots <= {
        "__future__",
        "argparse",
        "collections",
        "dataclasses",
        "fractions",
        "hashlib",
        "json",
        "math",
        "pathlib",
    }


def test_q_alpha_and_31_stage_schedule_match_primary_exactly() -> None:
    alpha = Cubic(0, 1, 0)
    assert alpha * alpha * alpha == Cubic(4, 0, 0)

    reference = theorem_stages()
    primary_cubic = tuple(reversed(fourth_order_suzuki_cubic_stages()))
    primary_interval = _theorem_interval_stages()

    assert len(reference) == 31
    assert [stage.fragment_index for stage in reference] == [
        stage.fragment_index for stage in primary_cubic
    ]
    assert [
        (
            stage.algebraic_coefficient.a0,
            stage.algebraic_coefficient.a1,
            stage.algebraic_coefficient.a2,
        )
        for stage in reference
    ] == [
        (stage.coefficient.a0, stage.coefficient.a1, stage.coefficient.a2)
        for stage in primary_cubic
    ]
    assert [
        (stage.coefficient_interval.lower, stage.coefficient_interval.upper)
        for stage in reference
    ] == [
        (stage.coefficient.lower, stage.coefficient.upper)
        for stage in primary_interval
    ]


def test_200_record_prefix_matches_primary_record_by_record() -> None:
    reference = tuple(iter_raw_theorem_records(200))
    primary = tuple(islice(_iter_raw_theorem_records(), 200))

    assert len(reference) == len(primary) == 200
    assert [record_payload(record) for record in reference] == [
        _primary_payload(record) for record in primary
    ]
    assert summarize_raw_theorem(200) == _primary_summary(200)
    assert summarize_raw_theorem(200).stream_sha256 == PREFIX_DIGEST


def test_compressed_200_prefix_matches_primary_schema_and_built_ledger() -> None:
    reference = tuple(iter_raw_theorem_records(200))
    primary = tuple(islice(_iter_raw_theorem_records(), 200))

    assert [compressed_record_payload(record) for record in reference] == [
        _primary_compressed_payload(record) for record in primary
    ]
    assert [canonical_compressed_record_bytes(record) for record in reference] == [
        canonical_raw_record_bytes(record) for record in primary
    ]
    summary = summarize_compressed_raw_theorem(200)
    assert summary == _primary_compressed_summary(200)
    assert summary.stream_sha256 == COMPRESSED_PREFIX_DIGEST

    ledger = build_compressed_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=200,
    )
    assert ledger.raw_record_count == summary.record_count
    assert ledger.raw_stream_sha256 == summary.stream_sha256
    assert {
        row.actual_word: row.raw_weight_upper for row in ledger.word_weights
    } == dict(summary.word_weights)


@pytest.mark.slow
def test_full_raw_stream_count_word_weights_and_digest_match_primary() -> None:
    reference = summarize_raw_theorem()
    primary = _primary_summary(None)

    assert reference.record_count == primary.record_count == 61_677
    assert reference.projected_record_count == PROJECTED_RECORD_COUNT
    assert reference.complete and primary.complete
    assert len(reference.word_weights) == len(primary.word_weights) == 4**5
    assert reference.word_weights == primary.word_weights
    assert reference.stream_sha256 == primary.stream_sha256 == FULL_DIGEST


@pytest.mark.slow
def test_full_compressed_stream_matches_primary_and_built_ledger(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference = summarize_compressed_raw_theorem()
    primary = _primary_compressed_summary(None)

    assert reference == primary
    assert reference.record_count == 61_677
    assert len(reference.word_weights) == 4**5
    assert reference.stream_sha256 == COMPRESSED_FULL_DIGEST

    def empty_block(representative: tuple[int, ...], _evaluator: object) -> object:
        return RepresentativeNormBlock(
            representative=representative,
            terms=(),
            groups=(),
            grouped_norm=Fraction(),
            triangle_norm=Fraction(),
        )

    # This test targets the raw stream boundary.  Avoid repeating the unrelated
    # finite-Pauli representative calculation in the full-stream slow lane.
    monkeypatch.setattr(compressed_module, "_representative_block", empty_block)
    ledger = build_compressed_xxz_ledger(XXZCompileSpec.pilot(Fraction(1, 2)))
    assert ledger.complete
    assert ledger.raw_record_count == reference.record_count
    assert ledger.raw_stream_sha256 == reference.stream_sha256
    assert {
        row.actual_word: row.raw_weight_upper for row in ledger.word_weights
    } == dict(reference.word_weights)


def test_record_decoder_replays_word_and_exact_positive_upper_weight() -> None:
    record = next(iter_raw_theorem_records())
    payload = record_payload(record)
    encoded = canonical_json_bytes(payload)

    assert load_canonical_json(encoded) == payload
    assert decode_record_payload(payload) == record
    assert record.weight_upper > 0


def test_compressed_decoder_replays_exact_interval_block_key_and_domain() -> None:
    record = next(iter_raw_theorem_records())
    payload = compressed_record_payload(record)

    assert decode_compressed_record_payload(payload) == record
    assert payload["stream_domain"] == COMPRESSED_STREAM_DOMAIN
    assert payload["block_key"] == payload["fragment_word"]
    assert payload["weight_interval"] == [
        [record.weight_lower.numerator, record.weight_lower.denominator],
        [record.weight_upper.numerator, record.weight_upper.denominator],
    ]


@pytest.mark.parametrize(
    "field,mutation",
    (
        ("stream_domain", "wrong-domain"),
        ("block_key", [3, 3, 3, 3, 3]),
        ("weight_interval", [[True, 1], [1, 1]]),
        ("weight_interval", [[1.0, 2], [1, 1]]),
        ("weight_interval", [[2, 2], [1, 1]]),
    ),
)
def test_compressed_decoder_rejects_schema_and_numeric_alias_mutations(
    field: str,
    mutation: object,
) -> None:
    payload = copy.deepcopy(
        compressed_record_payload(next(iter_raw_theorem_records()))
    )
    payload[field] = mutation

    with pytest.raises((TypeError, ValueError)):
        decode_compressed_record_payload(payload)


@pytest.mark.parametrize(
    "mutation",
    (
        lambda payload: payload.__setitem__("partial_sum_index", True),
        lambda payload: payload.__setitem__("partial_sum_index", 2.0),
        lambda payload: payload.__setitem__("unknown", 1),
        lambda payload: payload.__setitem__("weight_upper", [2, 2]),
    ),
)
def test_record_decoder_rejects_aliases_unknowns_and_noncanonical_fractions(
    mutation: object,
) -> None:
    payload = copy.deepcopy(record_payload(next(iter_raw_theorem_records())))
    mutation(payload)

    with pytest.raises((TypeError, ValueError)):
        decode_record_payload(payload)


def test_canonical_loader_rejects_bool_float_duplicate_and_pretty_json() -> None:
    payload = record_payload(next(iter_raw_theorem_records()))
    canonical = canonical_json_bytes(payload)
    duplicate = b'{"record_id":"shadow",' + canonical[1:]
    pretty = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()

    with pytest.raises(ValueError, match="duplicate JSON key"):
        load_canonical_json(duplicate)
    with pytest.raises(ValueError, match="canonical byte form"):
        load_canonical_json(pretty)
    for alias in (True, 2.0):
        forged = copy.deepcopy(payload)
        forged["partial_sum_index"] = alias
        with pytest.raises(TypeError, match="numeric aliases"):
            load_canonical_json(canonical_json_bytes(forged))


def test_strict_fraction_and_record_limit_reject_bool_aliases() -> None:
    with pytest.raises(TypeError):
        strict_fraction_pair([True, 1], "weight")
    with pytest.raises(ValueError):
        tuple(iter_raw_theorem_records(True))


def test_cli_runs_without_pythonpath_and_emits_canonical_prefix() -> None:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, str(REFERENCE), "--max-records", "200"],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr.decode()
    payload = load_canonical_json(completed.stdout)
    assert payload["record_count"] == 200
    assert payload["coverage_status"] == "bounded_prefix"
    assert payload["stream_sha256"] == PREFIX_DIGEST


def test_compressed_cli_runs_clean_and_emits_domain_separated_prefix() -> None:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [
            sys.executable,
            str(REFERENCE),
            "--stream-schema",
            "compressed-v1",
            "--max-records",
            "200",
        ],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr.decode()
    payload = load_canonical_json(completed.stdout)
    assert payload["stream_domain"] == COMPRESSED_STREAM_DOMAIN
    assert payload["summary_domain"] == (
        "grouped_xxz_compressed_raw_theorem_summary_v1"
    )
    assert payload["record_count"] == 200
    assert payload["stream_sha256"] == COMPRESSED_PREFIX_DIGEST
